# -*- coding: utf-8 -*-
# Archived study preprocessing source recovered from the pre-submission backup.
# Only BASE_PATH was made portable; Cinema 4D execution was not tested in this audit.
import c4d
from c4d import documents, utils
import os, math, time, traceback, gc
from datetime import datetime
from collections import deque, defaultdict

# ========= HIER ANPASSEN =========
BASE_PATH            = os.environ.get("WEEVIL_RAW_MESH_DIR", r"<RAW_MESH_DIRECTORY>")
OBJ_OUT_DIRNAME      = "Trochanter_clean_Obj"
PLY_OUT_DIRNAME      = "Trochanter_clean_Ply"
LOG_FILENAME         = "processing_log.txt"

# Pipeline-Parameter
REMESH_TARGET        = 10000     # Ziel-Polygone für ZRemesher
SMOOTH_ITERS         = 3         # Laplacian-Glättung: Iterationen
SMOOTH_ALPHA         = 0.40      # Laplacian-Glättung: Stärke

# „Holes“-Schritt: per Default AUS (robust für große Batches)
ENABLE_HOLES         = False     # True = versuchen, Loops zu schließen (mit Limits), False = Schritt überspringen

# Robustheits-Grenzen / Timeouts
HOLE_MAXLEN          = 4096      # Max. Loop-Länge fürs Capping (nur wenn ENABLE_HOLES=True)
MAX_BOUNDARY_EDGES   = 20000     # Wenn offene Randkanten > Grenze -> Holes skippen
MAX_LOOPS_TO_CAP     = 200       # Max. Anzahl Loops, die gecappt werden
MAX_SECONDS_HOLES    = 30.0      # Max. Zeit für Holes-Schritt (Sekunden)
MAX_SECONDS_PER_FILE = 600.0     # Harte Grenze pro Datei (Sekunden)
# =================================

# ---------- Logging ----------
class Logger:
    def __init__(self, base_path, filename):
        self.log_path = os.path.join(base_path, filename)
        os.makedirs(base_path, exist_ok=True)
        with open(self.log_path, "a", encoding="utf-8") as f:
            f.write("\n" + "="*84 + "\n")
            f.write(f"Run start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("="*84 + "\n")
    def write(self, msg):
        print(msg)
        with open(self.log_path, "a", encoding="utf-8") as f:
            f.write(msg + ("" if msg.endswith("\n") else "\n"))

log = Logger(BASE_PATH, LOG_FILENAME).write

# ---------- Utilities ----------
def ensure_dirs(base_path, obj_dirname, ply_dirname):
    obj_dir = os.path.join(base_path, obj_dirname)
    ply_dir = os.path.join(base_path, ply_dirname)
    os.makedirs(obj_dir, exist_ok=True)
    os.makedirs(ply_dir, exist_ok=True)
    return obj_dir, ply_dir

def iter_hierarchy(op):
    while op:
        yield op
        ch = op.GetDown()
        if ch:
            for c in iter_hierarchy(ch):
                yield c
        op = op.GetNext()

def list_all_objects(doc):
    root = doc.GetFirstObject()
    if not root: return []
    return list(iter_hierarchy(root))

def world_bbox(obj):
    if not isinstance(obj, c4d.PolygonObject): return None
    mg = obj.GetMg()
    pts = obj.GetAllPoints()
    if not pts: return None
    mn = c4d.Vector( 1e10,  1e10,  1e10)
    mx = c4d.Vector(-1e10, -1e10, -1e10)
    for p in pts:
        wp = mg * p
        mn.x = min(mn.x, wp.x); mn.y = min(mn.y, wp.y); mn.z = min(mn.z, wp.z)
        mx.x = max(mx.x, wp.x); mx.y = max(mx.y, wp.y); mx.z = max(mx.z, wp.z)
    return mn, mx

def bbox_center(obj):
    bb = world_bbox(obj)
    if not bb: return None
    mn, mx = bb
    return (mn + mx) * 0.5

def center_poly_at_origin(obj):
    c = bbox_center(obj)
    if c is None: return 0.0
    mg = obj.GetMg(); mg.off = mg.off - c
    obj.SetMg(mg)
    return (c.x*c.x + c.y*c.y + c.z*c.z) ** 0.5

def poly_counts(obj):
    return obj.GetPointCount(), obj.GetPolygonCount()

def count_tris_quads(obj):
    tris = others = 0
    for i in range(obj.GetPolygonCount()):
        p = obj.GetPolygon(i)
        if p.c == p.d: tris += 1
        else: others += 1
    return tris, others

def join_polys(doc, polys):
    if not polys: return None
    if len(polys) == 1: return polys[0]
    res = utils.SendModelingCommand(
        command=c4d.MCOMMAND_JOIN, list=polys,
        mode=c4d.MODELINGCOMMANDMODE_ALL, doc=doc,
        flags=0  # kein Undo => weniger RAM
    )
    if not res: return None
    return res[0] if isinstance(res, list) else res

def bake_to_global_clone(poly):
    clone = poly.GetClone()
    if clone is None: raise RuntimeError("Clone fehlgeschlagen.")
    mg = poly.GetMg()
    pts = clone.GetAllPoints()
    clone.SetAllPoints([mg * p for p in pts])
    clone.SetMg(c4d.Matrix())
    clone.Message(c4d.MSG_UPDATE)
    return clone

def check_timeout(start_time, limit_s, label):
    if (time.time() - start_time) > limit_s:
        raise TimeoutError(f"{label}: Zeitlimit {limit_s}s überschritten")

# ---------- Processing Steps ----------
def create_remesh(doc, poly, target_count=REMESH_TARGET):
    remesh = c4d.BaseObject(c4d.Oremesh)
    if remesh is None: raise RuntimeError("Remesh-Objekt konnte nicht erstellt werden.")
    remesh[c4d.ID_REMESHALGORITHM]                 = c4d.ID_REMESHALGORITHM_QUADREMESH
    remesh[c4d.ID_REMESH_POLYGONTARGET_MODE]       = c4d.ID_REMESH_POLYGONTARGET_MODE_ABSOLUTE
    remesh[c4d.ID_REMESH_POLYGONTARGET]            = int(target_count)
    remesh[c4d.ID_REMESH_QUADREMESH_HARDEDGEMODE]  = c4d.ID_REMESH_QUADREMESH_HARDEDGEMODE_AUTO
    remesh[c4d.ID_REMESH_QUADREMESH_ADAPTIVE]      = False
    remesh[c4d.ID_REMESH_SMOOTHITERATIONS]         = 0

    doc.InsertObject(remesh); poly.InsertUnder(remesh); c4d.EventAdd()
    result = utils.SendModelingCommand(
        command=c4d.MCOMMAND_MAKEEDITABLE, list=[remesh],
        mode=c4d.MODELINGCOMMANDMODE_ALL, doc=doc,
        flags=0
    )
    if not result: raise RuntimeError("Remesh fehlgeschlagen.")
    editable_root = result[0] if isinstance(result, list) else result
    polys = [editable_root] if isinstance(editable_root, c4d.PolygonObject) else \
            [o for o in iter_hierarchy(editable_root) if isinstance(o, c4d.PolygonObject)]
    final_poly = join_polys(doc, polys)
    if final_poly is None: raise RuntimeError("Kein Polygonobjekt nach Remesh gefunden.")

    for o in (editable_root, remesh, poly):
        try:
            if o and o.GetDocument() is not None and o != final_poly:
                o.Remove()
        except: pass
    final_poly.Remove(); doc.InsertObject(final_poly, pred=None, checknames=False)
    return final_poly

def triangulate_poly(doc, poly):
    return bool(utils.SendModelingCommand(
        command=c4d.MCOMMAND_TRIANGULATE, list=[poly],
        mode=c4d.MODELINGCOMMANDMODE_ALL, doc=doc, flags=0
    ))

def optimize_poly(doc, poly, tol):
    before = poly.GetPointCount()
    bc = c4d.BaseContainer()
    bc[c4d.MDATA_OPTIMIZE_TOLERANCE] = float(tol)
    bc[c4d.MDATA_OPTIMIZE_POINTS]    = True
    bc[c4d.MDATA_OPTIMIZE_POLYGONS]  = True
    utils.SendModelingCommand(
        command=c4d.MCOMMAND_OPTIMIZE, list=[poly],
        mode=c4d.MODELINGCOMMANDMODE_ALL, bc=bc, doc=doc, flags=0
    )
    after = poly.GetPointCount()
    return before, after

def build_point_adjacency(poly):
    cnt_p = poly.GetPointCount()
    adj = [set() for _ in range(cnt_p)]
    for pg in poly.GetAllPolygons():
        a,b,c,d = pg.a, pg.b, pg.c, pg.d
        adj[a].add(b); adj[b].add(a)
        adj[b].add(c); adj[c].add(b)
        adj[c].add(a); adj[a].add(c)
        if d != c:
            adj[c].add(d); adj[d].add(c)
            adj[d].add(a); adj[a].add(d)
    return adj

def laplacian_smooth_inplace(poly, iterations=SMOOTH_ITERS, strength=SMOOTH_ALPHA, logger=None):
    if iterations <= 0 or strength <= 0.0:
        if logger: logger("[Smooth] Übersprungen (iterations<=0 oder strength<=0).")
        return
    pts = list(poly.GetAllPoints()); n = len(pts)
    if n == 0:
        if logger: logger("[Smooth] Keine Punkte vorhanden.")
        return
    adj = build_point_adjacency(poly); alpha=float(strength)
    total_rms = 0.0
    for it in range(1, iterations+1):
        new_pts = [None]*n; sq=0.0; moved=0
        for i in range(n):
            nbrs = adj[i]
            if not nbrs: new_pts[i]=pts[i]; continue
            acc = c4d.Vector(0,0,0)
            for j in nbrs: acc += pts[j]
            avg = acc*(1.0/len(nbrs))
            p_old = pts[i]; p_new = p_old*(1.0-alpha) + avg*alpha
            new_pts[i]=p_new; d=p_new-p_old
            sq += d.x*d.x + d.y*d.y + d.z*d.z; moved += 1
        pts = new_pts
        rms = math.sqrt(sq/max(1,moved)); total_rms += rms
        if logger: logger(f"[Smooth] Iter {it}/{iterations}: RMS={rms:.6f}")
    poly.SetAllPoints(pts); poly.Message(c4d.MSG_UPDATE); c4d.EventAdd()
    if logger: logger(f"[Smooth] Gesamt-RMS={total_rms:.6f}")

def keep_largest_island(doc, poly, logger=None):
    pc = poly.GetPolygonCount()
    if pc == 0:
        if logger: logger("[Islands] Keine Polygone.")
        return 0,0,0
    polys = poly.GetAllPolygons()
    p2p = defaultdict(list)
    for i,pg in enumerate(polys):
        for idx in (pg.a,pg.b,pg.c,pg.d): p2p[idx].append(i)
    vis=[False]*pc; comps=[]
    for s in range(pc):
        if vis[s]: continue
        comp=[]; dq=deque([s]); vis[s]=True
        while dq:
            i=dq.popleft(); comp.append(i); pg=polys[i]
            for idx in (pg.a,pg.b,pg.c,pg.d):
                for j in p2p.get(idx,[]):
                    if not vis[j]: vis[j]=True; dq.append(j)
        comps.append(comp)
    n=len(comps)
    if n<=1:
        if logger: logger("[Islands] Nur 1 Insel – nichts zu löschen.")
        return 1, pc, 0
    comps.sort(key=len, reverse=True)
    keep=set(comps[0]); keep_n=len(keep)
    bs=poly.GetPolygonS(); bs.DeselectAll(); deleted=0
    for i in range(pc):
        if i not in keep: bs.Select(i); deleted+=1
    utils.SendModelingCommand(command=c4d.MCOMMAND_DELETE, list=[poly],
                              mode=c4d.MODELINGCOMMANDMODE_POLYGONSELECTION,
                              doc=doc, flags=0)
    optimize_poly(doc, poly, tol=0.0)
    if logger: logger(f"[Islands] Inseln gesamt: {n} | behalten: {keep_n} | gelöscht: {deleted}")
    return n, keep_n, deleted

def delete_degenerate_and_duplicate_tris(doc, poly, eps_rel=1e-12, logger=None):
    pts = poly.GetAllPoints(); polys = poly.GetAllPolygons(); n = poly.GetPolygonCount()
    if n == 0:
        if logger: logger("[Clean] Keine Polygone.")
        return 0,0
    bb = world_bbox(poly); diag2=1.0
    if bb:
        mn,mx = bb; diag2 = (mx-mn).GetSquaredLength() or 1.0
    area_eps = eps_rel*diag2
    dup=set(); to_del=[]
    for i,pg in enumerate(polys):
        a,b,c,d = pg.a,pg.b,pg.c,pg.d
        if d != c: continue
        pa,pb,pc = pts[a],pts[b],pts[c]
        area2 = (pb-pa).Cross(pc-pa).GetSquaredLength()*0.25
        if area2 <= area_eps: to_del.append(i); continue
        key=frozenset((a,b,c))
        if key in dup: to_del.append(i)
        else: dup.add(key)
    if not to_del:
        if logger: logger("[Clean] Degenerates/Duplicates: nichts zu löschen.")
        return 0,0
    bs=poly.GetPolygonS(); bs.DeselectAll()
    for i in to_del: bs.Select(i)
    utils.SendModelingCommand(command=c4d.MCOMMAND_DELETE, list=[poly],
                              mode=c4d.MODELINGCOMMANDMODE_POLYGONSELECTION,
                              doc=doc, flags=0)
    optimize_poly(doc, poly, tol=0.0)
    if logger: logger(f"[Clean] Gelöscht: {len(to_del)} (degenerate+duplicate).")
    return len(to_del),0

def close_holes_with_limits(doc, poly, start_time, logger=None):
    """Nur verwendet, wenn ENABLE_HOLES=True. Mit Grenzen & Zeitouts; wirft keine Exceptions."""
    try:
        t0 = time.time()
        polys=poly.GetAllPolygons(); pts=list(poly.GetAllPoints()); nF=len(polys)
        if nF==0:
            if logger: logger("[Holes] Keine Polygone.")
            return (0,0,0,"empty")

        # Kanten zählen
        edge_count=defaultdict(int)
        for pg in polys:
            tri=(pg.a,pg.b,pg.c) if pg.d==pg.c else (pg.a,pg.b,pg.c)
            for u,v in ((tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])):
                key=(min(u,v),max(u,v)); edge_count[key]+=1
        boundary_edges = sum(1 for _,cnt in edge_count.items() if cnt==1)
        if logger: logger(f"[Holes] Boundary-Edges={boundary_edges}")
        if boundary_edges==0:
            if logger: logger("[Holes] Keine offenen Ränder.")
            return (0,0,0,"no_boundary")
        if boundary_edges > MAX_BOUNDARY_EDGES:
            if logger: logger(f"[Holes] Übersprungen (Boundary>{MAX_BOUNDARY_EDGES}).")
            return (0,0,0,"too_many_edges")

        # Boundary-Graph
        boundary=defaultdict(list)
        for (u,v),cnt in edge_count.items():
            if cnt==1: boundary[u].append(v); boundary[v].append(u)

        visited=set(); loops=[]; loops_checked=0
        for start in list(boundary.keys()):
            if start in visited or not boundary[start]: continue
            check_timeout(start_time, MAX_SECONDS_PER_FILE, "Datei")
            if (time.time()-t0) > MAX_SECONDS_HOLES:
                if logger: logger(f"[Holes] Zeitlimit {MAX_SECONDS_HOLES}s erreicht → Abbruch.")
                break

            stack=[start]; comp=set()
            while stack:
                x=stack.pop()
                if x in comp: continue
                comp.add(x)
                for y in boundary.get(x,[]):
                    if y not in comp: stack.append(y)
            s=next((v for v in comp if len(boundary[v])==1), start)

            loop=[]; prev=None; cur=s; inner_t=time.time()
            while True:
                visited.add(cur); loop.append(cur)
                nbrs=[x for x in boundary[cur] if x!=prev]
                if not nbrs: loop=[]; break
                nxt=nbrs[0]
                if nxt==s:
                    loop.append(nxt); break
                prev,cur=cur,nxt
                if (time.time()-inner_t) > 2.0: loop=[]; break

            if not loop: continue
            uniq=[]; seen=set()
            for idx in loop:
                if idx not in seen: uniq.append(idx); seen.add(idx)
            if len(uniq)>=3 and uniq[0]==loop[-1]:
                loops.append(uniq[:-1])
                loops_checked += 1
                if loops_checked >= MAX_LOOPS_TO_CAP:
                    if logger: logger(f"[Holes] Max Loops {MAX_LOOPS_TO_CAP} erreicht → Rest skip.")
                    break

        found=len(loops)
        if found==0:
            if logger: logger("[Holes] Keine geschlossenen Rand-Loops gefunden.")
            return (0,0,0,"no_loops")

        if logger: logger(f"[Holes] Geschlossene Loops (begrenzt): {found}")

        added=0; capped=0
        for k,L in enumerate(loops,1):
            check_timeout(start_time, MAX_SECONDS_PER_FILE, "Datei")
            if (time.time()-t0) > MAX_SECONDS_HOLES:
                if logger: logger(f"[Holes] Zeitlimit erreicht beim Capping → Abbruch.")
                break
            m=len(L)
            if m<3:
                if logger: logger(f"[Holes]  - Loop {k}: m={m} → skip.")
                continue
            if m>HOLE_MAXLEN:
                if logger: logger(f"[Holes]  - Loop {k}: m={m}>{HOLE_MAXLEN} → skip.")
                continue
            center=c4d.Vector(0,0,0)
            for i in L: center+=pts[i]
            center*=1.0/m; pts.append(center)
            for i in range(1,m-1):
                a=L[0]; b=L[i]; c=L[i+1]
                poly.InsertPolygon(c4d.CPolygon(a,b,c)); added+=1
            capped+=1
            if logger: logger(f"[Holes]  - Loop {k}: m={m} → gecappt (+{m-2} Tris)")
        if added>0:
            poly.SetAllPoints(pts); poly.Message(c4d.MSG_UPDATE); c4d.EventAdd()
        if logger: logger(f"[Holes] Summary: Loops={found} | gecappt={capped} | neue Tris={added}")
        return (found,capped,added,"ok")
    except Exception as e:
        if logger: logger(f"[Holes] Fehler ignoriert: {e}")
        return (0,0,0,"error")

def signed_volume_of_trimesh(poly):
    pts=poly.GetAllPoints(); vol6=0.0
    for i in range(poly.GetPolygonCount()):
        pg=poly.GetPolygon(i)
        if pg.d!=pg.c: continue
        pa,pb,pc=pts[pg.a],pts[pg.b],pts[pg.c]
        vol6 += pa.Dot(pb.Cross(pc))
    return vol6/6.0

def ensure_normals_outward(poly, logger=None):
    vol = signed_volume_of_trimesh(poly)
    if vol==0.0:
        if logger: logger("[Normals] Volumen=0 → Richtung unklar (Mesh evtl. offen).")
        return False, vol
    if vol<0.0:
        polys=poly.GetAllPolygons()
        for i,pg in enumerate(polys):
            polys[i] = c4d.CPolygon(pg.a, pg.c, pg.b) if pg.d==pg.c \
                       else c4d.CPolygon(pg.a, pg.d, pg.c, pg.b)
        poly.SetAllPolygons(polys); poly.Message(c4d.MSG_UPDATE); c4d.EventAdd()
        if logger: logger("[Normals]  Volumen < 0 → Normalen invertiert.")
        return True, vol
    else:
        if logger: logger("[Normals] ✅ Volumen > 0 → Normalen plausibel nach außen.")
        return False, vol

# ---------- Export ----------
def export_obj_scene(doc, out_path):
    fmt = getattr(c4d, "FORMAT_OBJEXPORT", None)
    if fmt is None: return False
    try:
        ok = documents.SaveDocument(doc, out_path,
                                    c4d.SAVEDOCUMENTFLAGS_DONTADDTORECENTLIST, fmt)
        return bool(ok)
    except Exception:
        return False

def export_obj_fallback(poly_global, out_path):
    try:
        pts = poly_global.GetAllPoints()
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("# OBJ export (fallback)\n")
            for p in pts:
                f.write(f"v {p.x:.6f} {p.y:.6f} {p.z:.6f}\n")
            for i in range(poly_global.GetPolygonCount()):
                pg = poly_global.GetPolygon(i)
                if pg.d!=pg.c: continue
                a,b,c = pg.a+1, pg.b+1, pg.c+1
                f.write(f"f {a} {b} {c}\n")
        return True
    except Exception:
        return False

def export_ply_ascii(poly_global, out_path):
    try:
        pts = poly_global.GetAllPoints()
        faces=[]
        for i in range(poly_global.GetPolygonCount()):
            pg=poly_global.GetPolygon(i)
            if pg.d!=pg.c: continue
            faces.append((pg.a,pg.b,pg.c))
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("ply\nformat ascii 1.0\n")
            f.write(f"element vertex {len(pts)}\n")
            f.write("property float x\nproperty float y\nproperty float z\n")
            f.write(f"element face {len(faces)}\n")
            f.write("property list uchar int vertex_indices\nend_header\n")
            for p in pts: f.write(f"{p.x:.6f} {p.y:.6f} {p.z:.6f}\n")
            for (a,b,c) in faces: f.write(f"3 {a} {b} {c}\n")
        with open(out_path, "r", encoding="utf-8") as f:
            first=f.readline().strip(); second=f.readline().strip()
        return first=="ply" and second.startswith("format ascii")
    except Exception:
        return False

# ---------- Per-file pipeline in TEMP DOCUMENT ----------
def process_one_file(base_path, obj_dir, ply_dir, obj_filename):
    steps = []  # (stepname, ok)
    file_t0 = time.time()
    obj_path = os.path.join(base_path, obj_filename)
    base = os.path.splitext(obj_filename)[0]
    log(f"\n=== Start: {base} ===")
    success = True

    tdoc = documents.BaseDocument()
    try:
        check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")

        # Import
        try:
            ok = documents.MergeDocument(tdoc, obj_path, c4d.SCENEFILTER_OBJECTS)
            if not ok: raise RuntimeError("Import fehlgeschlagen.")
            c4d.EventAdd()
            steps.append(("Import", True))
        except Exception as e:
            steps.append(("Import", False)); raise

        # Polygonobjekt finden
        new_objs = list_all_objects(tdoc)
        poly = next((o for o in new_objs if isinstance(o, c4d.PolygonObject)), None)
        if poly is None:
            steps.append(("FindPoly", False)); raise RuntimeError("Kein Polygonobjekt gefunden.")
        poly.SetName(base); steps.append(("FindPoly", True))

        # Center
        try:
            mv = center_poly_at_origin(poly)
            log(f"[Center] Δ={mv:.4f}")
            steps.append(("Center", True))
        except Exception:
            steps.append(("Center", False))

        # Stats
        try:
            p0,f0 = poly_counts(poly)
            log(f"[Stats] Vor Remesh: Punkte={p0}, Polys={f0}")
        except: pass

        # Remesh
        try:
            check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
            rem = create_remesh(tdoc, poly, target_count=REMESH_TARGET); rem.SetName(base)
            mv2 = center_poly_at_origin(rem)
            p1,f1 = poly_counts(rem)
            log(f"[Remesh] Ziel={REMESH_TARGET} → Punkte={p1}, Polys={f1} (Δ={mv2:.4f})")
            steps.append(("Remesh", True))
        except Exception:
            steps.append(("Remesh", False)); raise

        # Triangulate
        try:
            check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
            tri_ok = triangulate_poly(tdoc, rem)
            tris,others = count_tris_quads(rem)
            log(f"[Triangulate] Dreiecke={tris}, Nicht-Dreiecke={others}")
            steps.append(("Triangulate", tri_ok and others==0))
        except Exception:
            steps.append(("Triangulate", False))

        # Optimize(0.01)
        try:
            check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
            b,a = optimize_poly(tdoc, rem, tol=0.01)
            log(f"[Optimize 0.01] Punkte {b} → {a}")
            steps.append(("Optimize_0.01", True))
        except Exception:
            steps.append(("Optimize_0.01", False))

        # Smooth
        try:
            check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
            p_before, f_before = poly_counts(rem)
            laplacian_smooth_inplace(rem, iterations=SMOOTH_ITERS, strength=SMOOTH_ALPHA, logger=log)
            p_after, f_after = poly_counts(rem)
            log(f"[Smooth] Punkte={p_after} (vorher {p_before}), Polys={f_after} (vorher {f_before})")
            steps.append(("Smooth", True))
        except Exception:
            steps.append(("Smooth", False))

        # Largest Island
        try:
            check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
            keep_largest_island(tdoc, rem, logger=log)
            steps.append(("LargestIsland", True))
        except Exception:
            steps.append(("LargestIsland", False))

        # Clean degenerates/duplicates
        try:
            check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
            delete_degenerate_and_duplicate_tris(tdoc, rem, logger=log)
            steps.append(("CleanDegDup", True))
        except Exception:
            steps.append(("CleanDegDup", False))

        # Holes (optional & robust; kann nie den Lauf sprengen)
        check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
        if ENABLE_HOLES:
            log("[Step] Holes …")
            _ = close_holes_with_limits(tdoc, rem, start_time=file_t0, logger=log)
            steps.append(("CloseHoles", True))
        else:
            log("[Step] Holes … (übersprungen)")
            steps.append(("CloseHoles", True))

        # Optimize(0.0)
        try:
            check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
            b2,a2 = optimize_poly(tdoc, rem, tol=0.0)
            log(f"[Optimize 0.0] Punkte {b2} → {a2}")
            steps.append(("Optimize_0.0", True))
        except Exception:
            steps.append(("Optimize_0.0", False))

        # Normals
        try:
            check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
            ensure_normals_outward(rem, logger=log)
            steps.append(("Normals", True))
        except Exception:
            steps.append(("Normals", False))

        # Export (baked global)
        check_timeout(file_t0, MAX_SECONDS_PER_FILE, "Datei")
        baked = bake_to_global_clone(rem)
        obj_out_path = os.path.join(obj_dir, f"{base}.obj")
        ply_out_path = os.path.join(ply_dir, f"{base}.ply")

        ok_obj = export_obj_scene(tdoc, obj_out_path)
        if not ok_obj:
            ok_obj = export_obj_fallback(baked.GetClone(), obj_out_path)
        log(f"[Export OBJ] {'OK' if ok_obj else 'FAIL'} → {obj_out_path}")
        steps.append(("Export_OBJ", ok_obj))

        ok_ply = export_ply_ascii(baked, ply_out_path)
        log(f"[Export PLY] {'OK' if ok_ply else 'FAIL'} → {ply_out_path}")
        steps.append(("Export_PLY", ok_ply))

    except TimeoutError as te:
        success = False
        steps.append(("TIMEOUT", False))
        log(f"[TIMEOUT] {te}")
    except Exception as e:
        success = False
        steps.append(("EXCEPTION", False))
        log(f"[ERROR] {type(e).__name__}: {e}")
        log(traceback.format_exc())
    finally:
        try:
            documents.KillDocument(tdoc)
        except: pass
        c4d.GeSyncMessage(c4d.EVMSG_DOCUMENTRECALCULATED)
        gc.collect()
        time.sleep(0.1)

    # Erfolgskriterium: Export OBJ & PLY müssen OK sein
    ok_all = success and all(k for n,k in steps if n in ("Export_OBJ","Export_PLY"))
    summary = "OK" if ok_all else "FAIL"
    detail  = ", ".join([f"{n}:{'✓' if k else '✗'}" for n,k in steps])
    dt = time.time()-file_t0
    log(f"[SUMMARY] {base} → {summary} ({dt:.2f}s) | {detail}")
    log(f"=== Ende: {base} ===")
    return ok_all

# ---------- Batch Main ----------
def main():
    obj_dir, ply_dir = ensure_dirs(BASE_PATH, OBJ_OUT_DIRNAME, PLY_OUT_DIRNAME)

    # Alle .obj im Ordner
    all_objs = sorted([f for f in os.listdir(BASE_PATH) if f.lower().endswith(".obj")])

    # Bereits fertig: nur überspringen, wenn OBJ **und** PLY existieren
    done = {os.path.splitext(f)[0] for f in os.listdir(obj_dir) if f.lower().endswith(".obj")}
    done &= {os.path.splitext(f)[0] for f in os.listdir(ply_dir) if f.lower().endswith(".ply")}

    todo = [f for f in all_objs if os.path.splitext(f)[0] not in done]
    log(f"\nGefundene OBJ: {len(all_objs)} | Bereits fertig: {len(done)} | ToDo: {len(todo)}")
    if not todo:
        log("Nichts zu tun – alle Dateien bereits exportiert.")
        return

    for i, fname in enumerate(todo, 1):
        log(f"\n--- [{i}/{len(todo)}] {fname} ---")
        process_one_file(BASE_PATH, obj_dir, ply_dir, fname)

    log("\nBatch abgeschlossen.")

if __name__ == "__main__":
    main()
