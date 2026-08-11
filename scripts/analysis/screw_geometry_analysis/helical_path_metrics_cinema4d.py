"""
Estimate helical-path metrics from ordered 3D semilandmark CSV files.

The script searches an input directory for ordered semilandmark CSV files,
fits a cylindrical axis by minimizing radial circle-fit error in projected
cross-section, and exports winding metrics to a semicolon-delimited CSV table.
It does not fit a full three-dimensional helix or regress axial position on
angular position; axial pitch is derived downstream from endpoint axial span
and absolute winding angle.
"""

import c4d
from c4d import Vector
import os, glob, math

# =========================
# SETTINGS
# =========================
FOLDER_PATH = os.environ.get("WINDING_INPUT_DIR", ".")
CSV_GLOB = "*.csv"

SCALE = 1.0
OUTPUT_CSV = "winding_metrics.csv"
REVERSE_ORDER = False

COARSE_RANGE_DEG = 40.0
COARSE_STEP_DEG = 8.0
FINE_RANGE_DEG = 10.0
FINE_STEP_DEG = 2.0
ULTRAFINE_RANGE_DEG = 2.0
ULTRAFINE_STEP_DEG = 0.5
# =========================


# -------------------------------------------------
# BASIC HELPERS
# -------------------------------------------------
def clean(s):
    return s.replace("\ufeff", "").strip()

def clean_lower(s):
    return clean(s).lower()

def normalize(v):
    l = v.GetLength()
    if l <= 1e-12:
        return Vector(1, 0, 0)
    return v / l

def dot(a, b):
    return a.Dot(b)

def cross(a, b):
    return a.Cross(b)

def clamp(x, a, b):
    return max(a, min(b, x))


# -------------------------------------------------
# CSV PARSING
# -------------------------------------------------
def split_line(line, delim):
    return [p.strip() for p in line.split(delim)]

def detect_delimiter_for_line(line):
    sc = line.count(";")
    cc = line.count(",")
    if sc > cc:
        return ";"
    elif cc > sc:
        return ","
    else:
        # fallback
        return ";"

def canonical_colname(s):
    s = clean_lower(s)
    s = s.replace('"', "")
    s = s.replace("'", "")
    s = s.replace(" ", "")
    s = s.replace("_", "")
    return s

def match_xyz_columns(header_parts):
    """
    Accept several common header variants for X/Y/Z coordinates.
    """
    canon = [canonical_colname(h) for h in header_parts]

    x_aliases = set(["x", "pos.x", "posx", "px", "coordx"])
    y_aliases = set(["y", "pos.y", "posy", "py", "coordy"])
    z_aliases = set(["z", "pos.z", "posz", "pz", "coordz"])

    ix = iy = iz = None

    for i, h in enumerate(canon):
        if h in x_aliases and ix is None:
            ix = i
        elif h in y_aliases and iy is None:
            iy = i
        elif h in z_aliases and iz is None:
            iz = i

    return ix, iy, iz

def looks_like_numeric_xyz_row(parts):
    """
    Diagnostic helper used while identifying plausible coordinate rows.
    Used only as a diagnostic fallback.
    """
    nums = 0
    for p in parts:
        try:
            float(clean(p).replace(",", "."))
            nums += 1
        except:
            pass
    return nums >= 3

def find_header_and_delimiter(lines):
    """
    Locate a header row containing X/Y/Z columns.
    Unterstützt sep=; am Anfang.
    """
    start_idx = 0
    forced_delim = None

    if lines and clean_lower(lines[0]).startswith("sep="):
        forced_delim = clean(lines[0])[4:]
        start_idx = 1

    # suche Header innerhalb der ersten ~10 Zeilen
    max_check = min(len(lines), start_idx + 12)

    for i in range(start_idx, max_check):
        line = lines[i]
        candidates = [forced_delim] if forced_delim else [";", ","]

        for delim in candidates:
            parts = split_line(line, delim)
            ix, iy, iz = match_xyz_columns(parts)
            if ix is not None and iy is not None and iz is not None:
                return i, delim, ix, iy, iz

    # fallback: delimiter anhand erster sinnvoller Zeile erraten und Header nochmal suchen
    for i in range(start_idx, max_check):
        delim = detect_delimiter_for_line(lines[i])
        parts = split_line(lines[i], delim)
        ix, iy, iz = match_xyz_columns(parts)
        if ix is not None and iy is not None and iz is not None:
            return i, delim, ix, iy, iz

    # Debug-Ausgabe vorbereiten
    preview = "\n".join(lines[:5])
    raise RuntimeError("Could not identify X/Y/Z columns. First lines:\n" + preview)

def read_points_csv(path):
    with open(path, "r", encoding="utf-8-sig") as f:
        lines = [l.strip() for l in f if l.strip()]

    if len(lines) < 2:
        raise RuntimeError("CSV file is too short")

    header_idx, delim, ix, iy, iz = find_header_and_delimiter(lines)

    pts = []
    for l in lines[header_idx + 1:]:
        parts = split_line(l, delim)
        if max(ix, iy, iz) >= len(parts):
            continue

        try:
            x = float(clean(parts[ix]).replace(",", ".")) * SCALE
            y = float(clean(parts[iy]).replace(",", ".")) * SCALE
            z = float(clean(parts[iz]).replace(",", ".")) * SCALE
            pts.append(Vector(x, y, z))
        except:
            continue

    if len(pts) < 3:
        raise RuntimeError("Too few valid points were found")

    if REVERSE_ORDER:
        pts.reverse()

    return pts


# -------------------------------------------------
# GEOMETRY
# -------------------------------------------------
def mean_point(pts):
    s = Vector(0, 0, 0)
    for p in pts:
        s += p
    return s / float(len(pts))

def covariance(pts, c):
    cxx = cyy = czz = cxy = cxz = cyz = 0.0
    for p in pts:
        d = p - c
        cxx += d.x*d.x
        cyy += d.y*d.y
        czz += d.z*d.z
        cxy += d.x*d.y
        cxz += d.x*d.z
        cyz += d.y*d.z
    n = float(len(pts))
    return [
        [cxx/n, cxy/n, cxz/n],
        [cxy/n, cyy/n, cyz/n],
        [cxz/n, cyz/n, czz/n],
    ]

def orthonormal_basis(axis):
    a = normalize(axis)
    tmp = Vector(1, 0, 0)
    if abs(dot(a, tmp)) > 0.9:
        tmp = Vector(0, 1, 0)
    u = normalize(tmp - a * dot(a, tmp))
    v = normalize(cross(a, u))
    return u, v

def rotate_vector(v, axis, angle_rad):
    k = normalize(axis)
    ca = math.cos(angle_rad)
    sa = math.sin(angle_rad)
    return v * ca + cross(k, v) * sa + k * dot(k, v) * (1.0 - ca)

def project_point_to_plane_coords(p, origin, axis, u, v):
    d = p - origin
    return dot(d, u), dot(d, v)

def point_on_axis_from_2d_center(origin, u, v, cx, cy):
    return origin + u * cx + v * cy


# -------------------------------------------------
# 3x3 EIGENDECOMPOSITION
# -------------------------------------------------
def identity3():
    return [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ]

def jacobi_eigen_decomposition_3x3(a, max_iter=50):
    A = [
        [float(a[0][0]), float(a[0][1]), float(a[0][2])],
        [float(a[1][0]), float(a[1][1]), float(a[1][2])],
        [float(a[2][0]), float(a[2][1]), float(a[2][2])]
    ]
    V = identity3()

    def max_offdiag(M):
        pairs = [(0,1), (0,2), (1,2)]
        best = pairs[0]
        best_val = abs(M[0][1])
        for i, j in pairs[1:]:
            val = abs(M[i][j])
            if val > best_val:
                best_val = val
                best = (i, j)
        return best, best_val

    for _ in range(max_iter):
        (p, q), off = max_offdiag(A)
        if off < 1e-12:
            break

        app = A[p][p]
        aqq = A[q][q]
        apq = A[p][q]

        if abs(apq) < 1e-20:
            continue

        tau = (aqq - app) / (2.0 * apq)
        t = 1.0 / (abs(tau) + math.sqrt(1.0 + tau*tau))
        if tau < 0.0:
            t = -t
        c = 1.0 / math.sqrt(1.0 + t*t)
        s = t * c

        for k in range(3):
            if k != p and k != q:
                Akp = A[k][p]
                Akq = A[k][q]
                A[k][p] = c * Akp - s * Akq
                A[p][k] = A[k][p]
                A[k][q] = s * Akp + c * Akq
                A[q][k] = A[k][q]

        A[p][p] = c*c*app - 2.0*s*c*apq + s*s*aqq
        A[q][q] = s*s*app + 2.0*s*c*apq + c*c*aqq
        A[p][q] = 0.0
        A[q][p] = 0.0

        for k in range(3):
            Vkp = V[k][p]
            Vkq = V[k][q]
            V[k][p] = c * Vkp - s * Vkq
            V[k][q] = s * Vkp + c * Vkq

    eigvals = [A[0][0], A[1][1], A[2][2]]
    eigvecs = [
        Vector(V[0][0], V[1][0], V[2][0]),
        Vector(V[0][1], V[1][1], V[2][1]),
        Vector(V[0][2], V[1][2], V[2][2]),
    ]

    idx = sorted(range(3), key=lambda i: eigvals[i], reverse=True)
    eigvals = [eigvals[i] for i in idx]
    eigvecs = [normalize(eigvecs[i]) for i in idx]
    return eigvals, eigvecs


# -------------------------------------------------
# 2D CIRCLE FIT
# -------------------------------------------------
def solve_3x3(A, b):
    M = [
        [float(A[0][0]), float(A[0][1]), float(A[0][2]), float(b[0])],
        [float(A[1][0]), float(A[1][1]), float(A[1][2]), float(b[1])],
        [float(A[2][0]), float(A[2][1]), float(A[2][2]), float(b[2])]
    ]

    for col in range(3):
        pivot = col
        max_abs = abs(M[col][col])
        for r in range(col + 1, 3):
            if abs(M[r][col]) > max_abs:
                max_abs = abs(M[r][col])
                pivot = r

        if max_abs < 1e-14:
            raise RuntimeError("Singuläres 3x3-System")

        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]

        piv = M[col][col]
        for j in range(col, 4):
            M[col][j] /= piv

        for r in range(3):
            if r == col:
                continue
            fac = M[r][col]
            for j in range(col, 4):
                M[r][j] -= fac * M[col][j]

    return [M[0][3], M[1][3], M[2][3]]

def fit_circle_2d(coords):
    n = float(len(coords))
    if n < 3:
        raise RuntimeError("Too few points for circle fitting")

    sx = sy = sxx = syy = sxy = 0.0
    sxz = syz = sz = 0.0

    for x, y in coords:
        z = x*x + y*y
        sx += x
        sy += y
        sxx += x*x
        syy += y*y
        sxy += x*y
        sxz += x*z
        syz += y*z
        sz += z

    A = [
        [sxx, sxy, sx],
        [sxy, syy, sy],
        [sx,  sy,  n ]
    ]
    b = [-sxz, -syz, -sz]

    sol = solve_3x3(A, b)
    Acoef, Bcoef, Ccoef = sol

    cx = -Acoef / 2.0
    cy = -Bcoef / 2.0
    rad2 = cx*cx + cy*cy - Ccoef
    if rad2 <= 0:
        raise RuntimeError("Invalid fitted circle radius")
    r = math.sqrt(rad2)

    dists = []
    for x, y in coords:
        d = math.sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy))
        dists.append(d)

    mean_r = sum(dists) / len(dists)
    rms = math.sqrt(sum((d - mean_r)*(d - mean_r) for d in dists) / len(dists))

    return cx, cy, r, rms


# -------------------------------------------------
# ANGLES
# -------------------------------------------------
def unwrap_angles(angle_list):
    if len(angle_list) < 2:
        return angle_list[:]

    unwrapped = [angle_list[0]]
    for a in angle_list[1:]:
        prev = unwrapped[-1]
        while a - prev > math.pi:
            a -= 2.0 * math.pi
        while a - prev < -math.pi:
            a += 2.0 * math.pi
        unwrapped.append(a)
    return unwrapped

def signed_winding_angle_deg_from_center(coords, cx, cy):
    angles = []
    for x, y in coords:
        ang = math.atan2(y - cy, x - cx)
        angles.append(ang)

    if len(angles) < 2:
        return 0.0

    angles_unwrapped = unwrap_angles(angles)
    total = angles_unwrapped[-1] - angles_unwrapped[0]
    return math.degrees(total)


# -------------------------------------------------
# AXIS FITTING
# -------------------------------------------------
def evaluate_axis(pts, axis):
    axis = normalize(axis)
    origin = mean_point(pts)
    u, v = orthonormal_basis(axis)

    coords = [project_point_to_plane_coords(p, origin, axis, u, v) for p in pts]
    cx, cy, r, rms = fit_circle_2d(coords)

    axis_point = point_on_axis_from_2d_center(origin, u, v, cx, cy)
    coords_centered = [project_point_to_plane_coords(p, axis_point, axis, u, v) for p in pts]

    dists = [math.sqrt(x*x + y*y) for x, y in coords_centered]
    mean_r = sum(dists) / len(dists)
    rms2 = math.sqrt(sum((d - mean_r)*(d - mean_r) for d in dists) / len(dists))

    return {
        "axis": axis,
        "axis_point": axis_point,
        "coords": coords_centered,
        "radius": mean_r,
        "rms": rms2
    }

def refine_axis_grid(pts, base_axis, range_deg, step_deg):
    base_axis = normalize(base_axis)
    b1, b2 = orthonormal_basis(base_axis)

    best = None
    n = int(round((2.0 * range_deg) / step_deg)) + 1
    vals = [(-range_deg + i * step_deg) for i in range(n)]

    for ax in vals:
        for ay in vals:
            cand = rotate_vector(base_axis, b1, math.radians(ax))
            cand = rotate_vector(cand, b2, math.radians(ay))
            cand = normalize(cand)

            try:
                ev = evaluate_axis(pts, cand)
            except:
                continue

            if (best is None) or (ev["rms"] < best["rms"]):
                best = ev

    if best is None:
        raise RuntimeError("Axis search failed")
    return best

def fit_best_axis_and_center(pts):
    c = mean_point(pts)
    cov = covariance(pts, c)
    eigvals, eigvecs = jacobi_eigen_decomposition_3x3(cov)

    candidates = []
    for seed in eigvecs:
        for sign in [1.0, -1.0]:
            candidates.append(normalize(seed * sign))

    best = None
    for seed in candidates:
        try:
            ev = refine_axis_grid(pts, seed, COARSE_RANGE_DEG, COARSE_STEP_DEG)
            if (best is None) or (ev["rms"] < best["rms"]):
                best = ev
        except:
            pass

    if best is None:
        raise RuntimeError("No usable initial axis was found")

    best = refine_axis_grid(pts, best["axis"], FINE_RANGE_DEG, FINE_STEP_DEG)
    best = refine_axis_grid(pts, best["axis"], ULTRAFINE_RANGE_DEG, ULTRAFINE_STEP_DEG)

    return best["axis_point"], best["axis"], best["coords"], best["radius"], best["rms"]


# -------------------------------------------------
# METRICS
# -------------------------------------------------
def number_of_turns_from_signed_deg(angle_deg):
    return angle_deg / 360.0

def start_end_distance(pts):
    return (pts[-1] - pts[0]).GetLength()

def axial_span(pts, axis):
    return abs(dot(pts[-1] - pts[0], normalize(axis)))


# -------------------------------------------------
# FORMAT
# -------------------------------------------------
def fmt_de_float(x):
    return ("%.6f" % x).replace(".", ",")


# -------------------------------------------------
# MAIN
# -------------------------------------------------
def main():
    files = sorted(glob.glob(os.path.join(FOLDER_PATH, CSV_GLOB)))
    if not files:
        raise RuntimeError("No CSV files were found")

    # Output-Datei ausschließen
    output_name_lower = OUTPUT_CSV.lower()
    files = [p for p in files if os.path.basename(p).lower() != output_name_lower]

    # Auch alte Metrics-Dateien sicherheitshalber ausschließen
    files = [p for p in files if "winding_metrics" not in os.path.basename(p).lower()]

    out_path = os.path.join(FOLDER_PATH, OUTPUT_CSV)

    with open(out_path, "w", encoding="utf-8", newline="") as f:
        f.write("sep=;\n")
        f.write("specimen_id;signed_winding_angle_deg;abs_winding_angle_deg;n_turns_signed;n_turns_abs;start_end_dist;axial_span;fit_radius;fit_rms\n")

        for path in files:
            specimen = os.path.splitext(os.path.basename(path))[0]

            try:
                pts = read_points_csv(path)

                axis_point, axis, coords2d, radius, rms = fit_best_axis_and_center(pts)

                signed_angle = signed_winding_angle_deg_from_center(coords2d, 0.0, 0.0)
                abs_angle = abs(signed_angle)

                turns_signed = number_of_turns_from_signed_deg(signed_angle)
                turns_abs = abs(turns_signed)
                dist = start_end_distance(pts)
                axsp = axial_span(pts, axis)

                f.write(
                    specimen + ";" +
                    fmt_de_float(signed_angle) + ";" +
                    fmt_de_float(abs_angle) + ";" +
                    fmt_de_float(turns_signed) + ";" +
                    fmt_de_float(turns_abs) + ";" +
                    fmt_de_float(dist) + ";" +
                    fmt_de_float(axsp) + ";" +
                    fmt_de_float(radius) + ";" +
                    fmt_de_float(rms) + "\n"
                )

                print("{}: signed={:.2f}°, abs={:.2f}°, rms={:.4f}".format(
                    specimen, signed_angle, abs_angle, rms
                ))

            except Exception as e:
                print("{}: SKIP ({})".format(specimen, e))
                f.write(specimen + ";;;;;;;;\n")

    print("\nSaved:", out_path)


if __name__ == "__main__":
    main()
