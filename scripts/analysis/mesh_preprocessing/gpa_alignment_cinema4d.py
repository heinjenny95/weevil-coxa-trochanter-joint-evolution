"""Cinema 4D GPA alignment script using three homologous landmarks.

This script loads paired OBJ meshes and landmark files, performs GPA
alignment with optional centroid-size normalization, aligns specimens to a
canonical three-landmark frame, and exports aligned meshes plus an optional
summary report.
"""

import os, math, c4d
from c4d import documents as docs, storage
# ===== Settings =====
USE_SCALE        = True     # Normalize by centroid size.
MAX_ITERS        = 50
TOL              = 1e-6
DRAW_MEAN_LM     = True
MEAN_SPHERE_DIV  = 300.0
LM_SPHERE_DIV    = 300.0
EXPORT_OBJS      = True
WRITE_REPORT     = True
PRINT_PROGRESS   = True
# ====================

# ---------- Math/Utils ----------
def dot(a,b): return a.x*b.x + a.y*b.y + a.z*b.z
def cross(a,b): return c4d.Vector(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
def norm(v):
    n = math.sqrt(dot(v,v))
    return v*(1.0/n) if n>0 else c4d.Vector(0)

def centroid(vs):
    s=c4d.Vector(0)
    for v in vs: s = s + v
    return s*(1.0/len(vs))

def center(vs):
    c=centroid(vs); return [v-c for v in vs], c

def centroid_size(vs):
    c=centroid(vs); acc=0.0
    for v in vs: d=v-c; acc += d.Dot(d)
    return math.sqrt(acc)

def normalize_centroid_size(vs):
    Xc,cx = center(vs); cs = centroid_size(vs)
    s = 1.0 if cs==0 else 1.0/cs
    return [p*s for p in Xc], cx, cs

def mat_mul_vec(M,v):
    return c4d.Vector(
        M[0][0]*v.x + M[0][1]*v.y + M[0][2]*v.z,
        M[1][0]*v.x + M[1][1]*v.y + M[1][2]*v.z,
        M[2][0]*v.x + M[2][1]*v.y + M[2][2]*v.z,
    )

def mat_mul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
           -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
           +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

def invert3(M):
    d = det3(M)
    if abs(d) < 1e-12: return [[1,0,0],[0,1,0],[0,0,1]]
    id = 1.0/d
    return [
        [(M[1][1]*M[2][2]-M[1][2]*M[2][1])*id, (M[0][2]*M[2][1]-M[0][1]*M[2][2])*id, (M[0][1]*M[1][2]-M[0][2]*M[1][1])*id],
        [(M[1][2]*M[2][0]-M[1][0]*M[2][2])*id, (M[0][0]*M[2][2]-M[0][2]*M[2][0])*id, (M[0][2]*M[1][0]-M[0][0]*M[1][2])*id],
        [(M[1][0]*M[2][1]-M[1][1]*M[2][0])*id, (M[0][1]*M[2][0]-M[0][0]*M[2][1])*id, (M[0][0]*M[1][1]-M[0][1]*M[1][0])*id],
    ]

def RMS_diff(A,B):
    acc=0.0
    for a,b in zip(A,B): d=a-b; acc += d.Dot(d)
    return math.sqrt(acc/len(A))

# ---------- IO ----------
def read_morphologika(path):
    pts=[]; section=None
    with open(path,"r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            s=line.strip()
            if not s: continue
            if s.startswith("["): section=s.lower(); continue
            if section=="[rawpoints]":
                if s.startswith("#"): continue
                p=s.split()
                if len(p)>=3: pts.append(c4d.Vector(float(p[0]),float(p[1]),float(p[2])))
    if len(pts)<3: raise ValueError(f"Too few landmarks in {os.path.basename(path)}")
    return pts

def read_obj_vertices_faces(path):
    verts=[]; faces=[]
    with open(path,"r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            if not line or line[0] in "#\r\n": continue
            parts=line.strip().split()
            if not parts: continue
            if parts[0]=='v' and len(parts)>=4:
                verts.append(c4d.Vector(float(parts[1]),float(parts[2]),float(parts[3])))
            elif parts[0]=='f' and len(parts)>=4:
                idx=[]
                for tok in parts[1:]:
                    vi = tok.split('/')[0]
                    if vi: idx.append(int(vi))
                if len(idx)>=3: faces.append(idx)
    if not verts or not faces: raise ValueError(f"OBJ missing vertices/faces: {os.path.basename(path)}")
    return verts, faces

def triangulate(faces):
    tris=[]; quads=[]
    for idx in faces:
        if len(idx)==3: tris.append(idx)
        elif len(idx)==4: quads.append(idx)
        else:
            v0=idx[0]
            for k in range(1,len(idx)-1): tris.append([v0,idx[k],idx[k+1]])
    return tris, quads

def build_poly(verts, faces):
    tris,quads = triangulate(faces)
    po=c4d.PolygonObject(len(verts), len(tris)+len(quads))
    for i,v in enumerate(verts): po.SetPoint(i,v)
    pi=0
    for a,b,c in tris:  po.SetPolygon(pi, c4d.CPolygon(a-1,b-1,c-1,c-1)); pi+=1
    for a,b,c,d in quads: po.SetPolygon(pi, c4d.CPolygon(a-1,b-1,c-1,d-1)); pi+=1
    po.Message(c4d.MSG_UPDATE); return po

def write_obj(path, points, faces):
    with open(path,"w",encoding="utf-8") as f:
        f.write("# GPA aligned export\n")
        for p in points: f.write(f"v {p.x:.6f} {p.y:.6f} {p.z:.6f}\n")
        for idx in faces:
            if len(idx)==3: f.write(f"f {idx[0]} {idx[1]} {idx[2]}\n")
            elif len(idx)==4: f.write(f"f {idx[0]} {idx[1]} {idx[2]} {idx[3]}\n")
            else:
                v0=idx[0]
                for k in range(1,len(idx)-1): f.write(f"f {v0} {idx[k]} {idx[k+1]}\n")

# ---------- Kabsch (Davenport) ----------
def _rotation_from_davenport(H):
    Sxx,Sxy,Sxz = H[0]; Syx,Syy,Syz=H[1]; Szx,Szy,Szz=H[2]
    tr=Sxx+Syy+Szz
    K=[[ tr,       Syz-Szy, Szx-Sxz, Sxy-Syx],
       [Syz-Szy,  Sxx-Syy-Szz, Sxy+Syx, Szx+Sxz],
       [Szx-Sxz,  Sxy+Syx, -Sxx+Syy-Szz, Syz+Szy],
       [Sxy-Syx,  Szx+Sxz, Syz+Szy, -Sxx-Syy+Szz]]
    q=[1.0,0.0,0.0,0.0]
    for _ in range(60):
        nq=[sum(K[i][j]*q[j] for j in range(4)) for i in range(4)]
        n=math.sqrt(sum(v*v for v in nq))
        if n<1e-12: break
        q=[v/n for v in nq]
    q0,qx,qy,qz=q
    R=[[0.0]*3 for _ in range(3)]
    R[0][0]=q0*q0+qx*qx-qy*qy-qz*qz; R[0][1]=2*(qx*qy-q0*qz);   R[0][2]=2*(qx*qz+q0*qy)
    R[1][0]=2*(qx*qy+q0*qz);          R[1][1]=q0*q0-qx*qx+qy*qy-qz*qz; R[1][2]=2*(qy*qz-q0*qx)
    R[2][0]=2*(qx*qz-q0*qy);          R[2][1]=2*(qy*qz+q0*qx);         R[2][2]=q0*q0-qx*qx-qy*qy+qz*qz
    return R

def kabsch(centered_X, centered_Y):
    H=[[0.0]*3 for _ in range(3)]
    for x,y in zip(centered_X, centered_Y):
        H[0][0]+=x.x*y.x; H[0][1]+=x.x*y.y; H[0][2]+=x.x*y.z
        H[1][0]+=x.y*y.x; H[1][1]+=x.y*y.y; H[1][2]+=x.y*y.z
        H[2][0]+=x.z*y.x; H[2][1]+=x.z*y.y; H[2][2]+=x.z*y.z
    return _rotation_from_davenport(H)

def apply_R(points,R):
    out=[]
    for p in points:
        out.append(c4d.Vector(
            R[0][0]*p.x + R[0][1]*p.y + R[0][2]*p.z,
            R[1][0]*p.x + R[1][1]*p.y + R[1][2]*p.z,
            R[2][0]*p.x + R[2][1]*p.y + R[2][2]*p.z
        ))
    return out

# ---------- 3-LM-Frame ----------
def frame_from_3lm(lm3):
    # Landmark order is fixed as [L1, L2, L3].
    p1,p2,p3 = lm3[0], lm3[1], lm3[2]
    X = norm(p2 - p1)
    v = p3 - p1
    v = v - X * dot(X, v)
    Y = norm(v) if v.GetLength()>0 else c4d.Vector(0,1,0)
    Z = norm(cross(X, Y))
    Y = norm(cross(Z, X))
    return X, Y, Z

def frame_matrix(X,Y,Z):
    return [[X.x,Y.x,Z.x],[X.y,Y.y,Z.y],[X.z,Y.z,Z.z]]

# ---------- Visual ----------
def make_group(name, kill_old=True):
    doc=docs.GetActiveDocument()
    if kill_old:
        old=doc.SearchObject(name)
        if old and old.CheckType(c4d.Onull): old.Remove()
    g=c4d.BaseObject(c4d.Onull); g.SetName(name); doc.InsertObject(g); return g

def spawn_spheres(parent, pts, radius, color):
    for v in pts:
        sp=c4d.BaseObject(c4d.Osphere)
        sp[c4d.PRIM_SPHERE_RAD]=radius
        sp.SetRelPos(v)
        sp[c4d.ID_BASEOBJECT_USECOLOR]=2
        sp[c4d.ID_BASEOBJECT_COLOR]=color
        sp.InsertUnder(parent)

def bbox_span(point_lists):
    xs=[]; ys=[]; zs=[]
    for V in point_lists:
        for p in V: xs.append(p.x); ys.append(p.y); zs.append(p.z)
    if not xs: return 1.0
    return max(max(xs)-min(xs), max(ys)-min(ys), max(zs)-min(zs))

# ---------- Pairing ----------
def scan_pairs_by_dirs(obj_dir, lm_dir):
    obj_files=[f for f in os.listdir(obj_dir) if f.lower().endswith(".obj")]
    lm_files =[f for f in os.listdir(lm_dir)  if f.lower().endswith(".txt")]
    obj_map={os.path.splitext(f)[0].lower(): os.path.join(obj_dir,f) for f in obj_files}
    lm_map ={os.path.splitext(f)[0].lower(): os.path.join(lm_dir,f)  for f in lm_files}
    keys=sorted(set(obj_map.keys()) & set(lm_map.keys()))
    pairs=[(k,obj_map[k],lm_map[k]) for k in keys]
    return pairs

# ---------- GPA ----------
def gpa_align_landmarks(all_lm, use_scale=True, max_iters=50, tol=1e-6):
    n=len(all_lm); m=len(all_lm[0])
    centered=[]; pre_cx=[]; pre_cs=[]
    for lm in all_lm:
        if use_scale: Xi,cx,cs = normalize_centroid_size(lm)
        else:         Xi,cx    = center(lm); cs=1.0
        centered.append(Xi); pre_cx.append(cx); pre_cs.append(cs)
    mean_shape=[c4d.Vector(p) for p in centered[0]]
    rot_list=[None]*n
    for it in range(max_iters):
        aligned=[]
        for i in range(n):
            R=kabsch(centered[i], mean_shape)
            rot_list[i]=R
            aligned.append(apply_R(centered[i], R))
        new_mean=[c4d.Vector(0) for _ in range(m)]
        for j in range(m):
            s=c4d.Vector(0)
            for i in range(n): s = s + aligned[i][j]
            new_mean[j]=s*(1.0/n)
        new_mean,_=center(new_mean)
        if use_scale:
            cs=centroid_size(new_mean)
            if cs>0: new_mean=[p*(1.0/cs) for p in new_mean]
        diff=RMS_diff(mean_shape,new_mean)
        mean_shape=new_mean
        print(f"[GPA] Iter {it+1}: RMS change {diff:.8e}")
        if diff<tol: break
    aligned_lm=[apply_R(centered[i], rot_list[i]) for i in range(n)]
    return aligned_lm, mean_shape, rot_list, pre_cx, pre_cs

# ---------- Mean-Mesh (nearest-vertex) ----------
def mean_mesh_from_nearest_vertex(template_pts, all_pts_lists):
    n=len(all_pts_lists); nv=len(template_pts); mean_pts=[]
    for vi in range(nv):
        refp=template_pts[vi]
        acc=c4d.Vector(0)
        for pts in all_pts_lists:
            best=None; bestd=None
            for p in pts:
                d=(p-refp).GetSquaredLength()
                if bestd is None or d<bestd: bestd=d; best=p
            acc = acc + best
        mean_pts.append(acc*(1.0/n))
        if PRINT_PROGRESS and (vi % max(1,nv//20) == 0): print(f"[MeanMesh] {vi+1}/{nv}")
    return mean_pts

# ---------- Pipeline ----------
def main():
    obj_dir = storage.LoadDialog(title="Select OBJ directory", flags=c4d.FILESELECT_DIRECTORY)
    if not obj_dir: return
    lm_dir  = storage.LoadDialog(title="Select landmark directory (.txt)", flags=c4d.FILESELECT_DIRECTORY)
    if not lm_dir: return
    out_dir = storage.LoadDialog(title="Select output directory", flags=c4d.FILESELECT_DIRECTORY)
    if not out_dir: return

    pairs = scan_pairs_by_dirs(obj_dir, lm_dir)
    if not pairs: raise RuntimeError("No matching OBJ/TXT pairs were found.")

    names=[]; meshes=[]; faces_list=[]; lm_list=[]
    for k,objp,lmp in pairs:
        V,F = read_obj_vertices_faces(objp)
        LM  = read_morphologika(lmp)
        if len(LM)<3: raise ValueError(f"Fewer than three landmarks: {k}")
        LM = LM[:3]
        names.append(k); meshes.append(V); faces_list.append(F); lm_list.append(LM)

    m=len(lm_list[0])
    for LM in lm_list:
        if len(LM)!=m: raise ValueError("All landmark sets must contain the same number of landmarks (3).")

    print(f"=== GPA START (3 LMs) === N={len(lm_list)} | SCALE={'ON' if USE_SCALE else 'OFF'}")
    aligned_lm, mean_lm, R_list, pre_cx, pre_cs = gpa_align_landmarks(
        lm_list, use_scale=USE_SCALE, max_iters=MAX_ITERS, tol=TOL
    )

    # Align all specimens to the reference three-landmark frame.
    refX,refY,refZ = frame_from_3lm(aligned_lm[0])
    Mr = frame_matrix(refX,refY,refZ)
    for i in range(len(aligned_lm)):
        X,Y,Z = frame_from_3lm(aligned_lm[i])

        if dot(Y, refY) < 0: Y = -Y

        Z = norm(cross(X, Y))
        if dot(Z, refZ) < 0: Z = -Z; Y = -Y
        Ms = frame_matrix(X,Y,Z)
        Rfix = mat_mul(Mr, invert3(Ms))
        aligned_lm[i] = [mat_mul_vec(Rfix, p) for p in aligned_lm[i]]
        R_list[i]     = mat_mul(Rfix, R_list[i])

    # Final scale and translation parameters.
    final_s=[]; final_t=[]
    for i in range(len(lm_list)):
        s = 1.0/pre_cs[i] if USE_SCALE else 1.0
        cx=pre_cx[i]; R=R_list[i]
        rcx = mat_mul_vec(R, cx)
        t   = c4d.Vector(-s*rcx.x, -s*rcx.y, -s*rcx.z)
        final_s.append(s); final_t.append(t)

    grp_mean=make_group("GPA_MEAN",True)
    grp_aln =make_group("GPA_ALIGNED",True)

    if DRAW_MEAN_LM:
        span=bbox_span(meshes)
        rad=max(span/MEAN_SPHERE_DIV, 0.02)
        spawn_spheres(grp_mean, mean_lm, rad, c4d.Vector(1.0,0.85,0.0))

    aligned_dir=os.path.join(out_dir,"Aligned_GPA")
    if EXPORT_OBJS and not os.path.isdir(aligned_dir): os.makedirs(aligned_dir,exist_ok=True)
    rep_path=os.path.join(out_dir,"GPA_report_3LM.csv") if WRITE_REPORT else None
    rep=None
    if WRITE_REPORT:
        rep=open(rep_path,"w",encoding="utf-8")
        rep.write("name,scale,tx,ty,tz,rms3_to_ref\n")

    transformed_points=[]
    for i,(key,V,F,LM) in enumerate(zip(names,meshes,faces_list,lm_list)):
        R=R_list[i]; s=final_s[i]; t=final_t[i]
        out_pts=[]
        for p in V:
            q=mat_mul_vec(R,p)
            out_pts.append(c4d.Vector(s*q.x + t.x, s*q.y + t.y, s*q.z + t.z))
        transformed_points.append(out_pts)

        po=build_poly(out_pts, F)
        po.SetName(f"GPA_{key}")
        po[c4d.ID_BASEOBJECT_USECOLOR]=2
        po[c4d.ID_BASEOBJECT_COLOR]=c4d.Vector(0.65,0.9,0.95)
        po.InsertUnder(grp_aln)

        # Landmark check spheres.
        LMc=[p - pre_cx[i] for p in LM]
        if USE_SCALE and pre_cs[i]!=0: LMc=[p*(1.0/pre_cs[i]) for p in LMc]
        LMf=[mat_mul_vec(R,p) for p in LMc]
        LMf=[p + t for p in LMf]
        rad_lm=max(bbox_span([V])/LM_SPHERE_DIV, 0.02)
        spawn_spheres(po, LMf, rad_lm, c4d.Vector(0.0,1.0,0.3))

        rms=RMS_diff(aligned_lm[i], aligned_lm[0])
        if WRITE_REPORT: rep.write(f"{key},{s:.10f},{t.x:.10f},{t.y:.10f},{t.z:.10f},{rms:.10f}\n")
        if EXPORT_OBJS: write_obj(os.path.join(aligned_dir,f"{key}_aligned.obj"), out_pts, F)
        print(f"[OK] {key} | RMS3={rms:.6f}")

    # Mean mesh using nearest-vertex correspondence.
    if transformed_points:
        template_pts = transformed_points[0]
        mean_pts = mean_mesh_from_nearest_vertex(template_pts, transformed_points)
        mean_po = build_poly(mean_pts, faces_list[0])
        mean_po.SetName("GPA_MEAN_MESH")
        mean_po[c4d.ID_BASEOBJECT_USECOLOR]=2
        mean_po[c4d.ID_BASEOBJECT_COLOR]=c4d.Vector(1.0,0.93,0.2)
        mean_po.InsertUnder(grp_mean)
        if EXPORT_OBJS:
            write_obj(os.path.join(aligned_dir,"GPA_MEAN_MESH.obj"), mean_pts, faces_list[0])

    if rep: rep.close()
    c4d.EventAdd(); c4d.CallCommand(12148)
    print("=== GPA DONE (3 LMs) ===")

if __name__ == "__main__":
    try: main()
    except Exception:
        import traceback; traceback.print_exc()
