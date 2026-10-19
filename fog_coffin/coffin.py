"""
Fog-machine coffin -- parametric model (CadQuery).

Coffin is a 6-sided "toe-pincher". Head at y=0, foot at +y, z=0 is the ground.
Run:  python3 coffin.py [--no-texture] [--out DIR]
Everything is driven by the PARAMS dict below.
"""
import math, os, sys, json
import numpy as np
import cadquery as cq
from shapely.geometry import Polygon, LineString, Point, MultiPolygon, box
from shapely.ops import unary_union
from shapely import affinity

PARAMS = dict(
    # ---- overall outer size (mm, measured at outside of walls) ----
    head_w=170.0,       # width of head end
    shoulder_w=260.0,   # width at the shoulders (widest point)
    foot_w=150.0,       # width of foot end
    shoulder_y=120.0,   # distance head -> shoulders
    length=430.0,       # head -> foot
    height=180.0,       # wall height (lid sits on top)
    wall=3.2,           # wall thickness (4 lines @ 0.8 nozzle)
    # ---- lid ----
    lid_t=5.0,          # lid plate thickness
    lid_overhang=6.0,   # lid lip beyond walls (not on hinge side)
    skirt_depth=10.0,   # locating skirt under lid
    skirt_t=2.4,
    skirt_clear=1.0,
    lid_split_y=215.0,  # lid printed in two halves split here
    # ---- hinge (lift-off barrel hinge on the right long wall) ----
    knuckle_r=5.5,
    hinge_clear=0.7,    # radial clearance barrel <-> wall groove
    hinge_gap=0.4,      # axial gaps (tested round 2)
    cradle_len=8.0,
    cradle_s=(22.0, 150.0, 280.0),   # cradle centres along hinge axis from shoulder corner
    pin_r=2.5,
    slot_w=5.0,         # cradle slot; 5.0 won fit test round 2 (5.4-5.8 sloppy)
    slot_angle=30.0,    # cradle slot opens inboard and this many degrees upward
    # ---- corner joints ----
    flange_t=5.0,
    flange_depth=18.0,
    flange_bolt_z=(20.0, 85.0, 150.0),
    # ---- M3 hardware ----
    m3_clear=3.5,
    nut_af=5.9,         # across flats incl. print allowance (M3 nut = 5.5)
    nut_depth=2.8,
    # ---- lid stay (ratchet) ----
    stay_s=40.0,        # axial position of lid lug
    lid_pin_uv=(20.0, -14.0),   # lid pivot, relative to hinge axis (u inboard, v up)
    wall_pin_uv=(16.0, -70.0),  # wall pin
    stay_notch_deg=(30.0, 55.0, 80.0),
    stay_max_deg=85.0,
    stay_w=17.0,
    stay_nut_depth=3.6,  # deeper nut pockets on stay lugs so M3x12 fully engages
    stay_washer_t=1.2,
    stay_boss_d=6.0,   # printed bushing the stay rides on (screw clamps the bushing, not the stay)
    stay_t=5.0,
    # ---- cord notch & plugs ----
    notch_w=30.0,
    notch_h=24.0,
    plug_flange=6.0,
    plug_flange_t=2.0,
    cord_slot_w=9.0,
    # ---- stakes ----
    stake_hole=6.5,
    # ---- texture ----
    grain_depth=0.6,
    grain_w=1.2,
    seam_depth=1.0,
    seam_w=1.8,
    lid_cross=True,
)

# ------------------------------------------------------------------ geometry helpers
def unit(v):
    v = np.asarray(v, float); return v / np.linalg.norm(v)

def V(p):
    p = list(p); p += [0.0] * (3 - len(p)); return cq.Vector(*map(float, p))

def outline(p):
    return np.array([[-p['head_w']/2, 0], [p['head_w']/2, 0],
                     [p['shoulder_w']/2, p['shoulder_y']], [p['foot_w']/2, p['length']],
                     [-p['foot_w']/2, p['length']], [-p['shoulder_w']/2, p['shoulder_y']]], float)

def edge_frames(P):
    n = len(P); E = []; M = []
    for i in range(n):
        e = unit(P[(i+1) % n] - P[i]); E.append(e); M.append(np.array([-e[1], e[0]]))  # inward (CCW)
    return np.array(E), np.array(M)

def line_isect(p1, d1, p2, d2):
    A = np.array([d1, -d2]).T; t = np.linalg.solve(A, p2 - p1); return p1 + t[0]*d1

def offset_poly(P, d):
    n = len(P); E, M = edge_frames(P)
    d = [d]*n if np.isscalar(d) else d
    out = []
    for i in range(n):
        j = (i-1) % n
        out.append(line_isect(P[j] + d[j]*M[j], E[j], P[i] + d[i]*M[i], E[i]))
    return np.array(out)

def prism(pts, z0, z1):
    return (cq.Workplane("XY").workplane(offset=z0)
            .polyline([tuple(map(float, q)) for q in pts]).close().extrude(z1 - z0).val())

def poly_solid(plane, poly, depth):
    """Extrude a shapely polygon (with holes) drawn in `plane` by depth (negative = against normal)."""
    wp = cq.Workplane(plane)
    polys = [poly] if isinstance(poly, Polygon) else list(poly.geoms)
    solids = []
    for pg in polys:
        if pg.is_empty or pg.area < 0.05: continue
        pg = pg.simplify(0.03)
        w = cq.Workplane(plane).polyline(list(pg.exterior.coords)[:-1]).close()
        for hole in pg.interiors:
            w = w.polyline(list(hole.coords)[:-1]).close()
        try:
            solids.append(w.extrude(depth).val())
        except Exception as ex:
            print("  (skipped a texture polygon:", ex, ")")
    return solids

def cut_many(shape, tools, batch=40):
    for k in range(0, len(tools), batch):
        shape = shape.cut(*tools[k:k+batch])
    return shape

def place(shape, origin, x3, y3, z3):
    """Map local coords (x,y,z) -> origin + x*x3 + y*y3 + z*z3 (must be a proper rotation)."""
    from OCP.gp import gp_Trsf
    t = gp_Trsf()
    o = np.asarray(origin, float)
    t.SetValues(x3[0], y3[0], z3[0], o[0],
                x3[1], y3[1], z3[1], o[1],
                x3[2], y3[2], z3[2], o[2])
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    return cq.Shape.cast(BRepBuilderAPI_Transform(shape.wrapped, t, True).Shape())

def to_local(shape, origin, x3, y3, z3):
    """Inverse of place(): world -> local where local axes are x3,y3,z3 at origin."""
    from OCP.gp import gp_Trsf
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    R = np.array([x3, y3, z3], float)          # rows
    o = -R @ np.asarray(origin, float)
    t = gp_Trsf()
    t.SetValues(R[0,0], R[0,1], R[0,2], o[0], R[1,0], R[1,1], R[1,2], o[1], R[2,0], R[2,1], R[2,2], o[2])
    return cq.Shape.cast(BRepBuilderAPI_Transform(shape.wrapped, t, True).Shape())

def hexpocket(origin, normal, depth, af):
    n = unit(normal); a = unit(np.cross(n, [0, 0, 1]) if abs(n[2]) < 0.9 else np.cross(n, [1, 0, 0]))
    pl = cq.Plane(origin=V(origin), xDir=V(a), normal=V(n))
    return cq.Workplane(pl).polygon(6, af / math.cos(math.pi/6)).extrude(depth).val()

def cyl(r, h, origin, d):
    return cq.Solid.makeCylinder(float(r), float(h), V(origin), V(unit(d)))

# ------------------------------------------------------------------ wood grain (2D)
def grain_pattern(x0, x1, y0, y1, spacing, seed, knots=3, seams=(), clip=None, avoid=None, p=PARAMS):
    rng = np.random.default_rng(seed)
    gw, sw = p['grain_w'], p['seam_w']
    kn = []
    for _ in range(knots):
        kx = rng.uniform(x0 + 0.15*(x1-x0), x1 - 0.15*(x1-x0)); ky = rng.uniform(y0 + 12, y1 - 12)
        kn.append((kx, ky, rng.uniform(4.0, 7.0)))
    xs = np.linspace(x0 - 5, x1 + 5, int((x1 - x0) / 2.5) + 4)
    lines = []
    yb = y0 + spacing * rng.uniform(0.3, 0.8)
    while yb < y1:
        a1, l1, f1 = rng.uniform(0.4, 1.4), rng.uniform(90, 220), rng.uniform(0, 6.3)
        a2, l2, f2 = rng.uniform(0.2, 0.7), rng.uniform(25, 60), rng.uniform(0, 6.3)
        y = yb + a1*np.sin(2*np.pi*xs/l1 + f1) + a2*np.sin(2*np.pi*xs/l2 + f2)
        for kx, ky, kr in kn:
            dy = yb - ky
            y = y + np.sign(dy if dy != 0 else 1) * kr * 1.9 * np.exp(-((xs-kx)/(kr*3.2))**2) * np.exp(-(dy/(kr*3.0))**2)
        lines.append(LineString(np.c_[xs, y]).buffer(gw/2, cap_style=2))
        yb += spacing * rng.uniform(0.75, 1.25)
    for kx, ky, kr in kn:
        for s in (1.0, 0.55):
            e = affinity.scale(Point(kx, ky).buffer(kr*s, 48), 1.7, 0.75)
            lines.append(e.exterior.buffer(gw/2))
    g = unary_union(lines)
    region = box(x0, y0, x1, y1) if clip is None else clip
    if avoid is not None: region = region.difference(avoid)
    g = g.intersection(region)
    deep = []
    for s in seams:   # plank seams (horizontal lines y = s)
        deep.append(LineString([(x0-5, s), (x1+5, s)]).buffer(sw/2, cap_style=2))
    deep = unary_union(deep).intersection(box(x0, y0, x1, y1) if clip is None else clip) if deep else None
    return g, deep

# ------------------------------------------------------------------ the model
def build(p=PARAMS, texture=True):
    H, T = p['height'], p['wall']
    P = outline(p); E, M = edge_frames(P); n = 6
    Q = offset_poly(P, T)
    O3 = [np.r_[-M[i], 0] for i in range(n)]       # outward normals (3D)
    E3 = [np.r_[E[i], 0] for i in range(n)]
    M3 = [np.r_[M[i], 0] for i in range(n)]
    Z3 = np.array([0, 0, 1.0])
    L_edge = [np.linalg.norm(P[(i+1) % n] - P[i]) for i in range(n)]

    # hinge frame on edge 2 (right long wall)
    rk = p['knuckle_r']; az = H + p['lid_t'] - rk
    he, hm = E3[2], M3[2]
    A0 = np.r_[P[2] + rk*M[2], az]
    def hpt(s, u, v): return A0 + s*he + u*hm + v*Z3
    def hplane(s): return cq.Plane(origin=V(hpt(s, 0, 0)), xDir=V(hm), normal=V(he))
    def hcyl(R, s0, s1): return cq.Workplane(hplane(s0)).circle(R).extrude(s1 - s0).val()
    def hpoly(pts, s0, s1): return cq.Workplane(hplane(s0)).polyline(pts).close().extrude(s1 - s0).val()
    def to_su(xy):   # world xy -> (s,u) in hinge frame
        d = np.asarray(xy) - A0[:2]; return np.array([d @ E[2], d @ M[2]])

    LP = offset_poly(P, [-p['lid_overhang']]*2 + [0] + [-p['lid_overhang']]*3)
    LPg = offset_poly(P, [-p['lid_overhang'] - p['hinge_clear']]*2 + [-p['hinge_clear']] + [-p['lid_overhang'] - p['hinge_clear']]*3)
    S_MIN, S_MAX = -60.0, L_edge[2] + 60.0
    cr_len, gap, hc = p['cradle_len'], p['hinge_gap'], p['hinge_clear']
    cradles = [(c - cr_len/2, c + cr_len/2) for c in p['cradle_s']]

    def intervals_excluding(a, b, holes):
        out = []; cur = a
        for h0, h1 in sorted(holes):
            if h0 > cur: out.append((cur, h0))
            cur = max(cur, h1)
        if cur < b: out.append((cur, b))
        return out

    # groove (cut into walls) everywhere along axis except cradle spans, clipped to lid footprint
    groove = [hcyl(rk + hc, a, b) for a, b in intervals_excluding(S_MIN, S_MAX, cradles)]
    lpg_prism = prism(LPg, az - 20, az + 20)
    groove = [g.intersect(lpg_prism) for g in groove]

    walls = []
    for i in range(n):
        j = (i + 1) % n
        w = prism([P[i], P[j], Q[j], Q[i]], 0, H)
        # --- corner flanges (start corner i: nut side; end corner j: head side)
        for corner, side_dir, nut in ((i, E[i], True), (j, -E[i], False)):
            Pc, Qc = P[corner], Q[corner]
            b = unit(Qc - Pc)
            nn = unit(side_dir - (side_dir @ b) * b)
            half = math.acos(np.clip(b @ side_dir, -1, 1))       # angle between bisector and wall
            t0 = p['flange_t'] / math.tan(half)
            mm = M[i]
            D = p['flange_depth']
            pts = [Qc - 0.6*mm, Qc + D*b, Qc + D*b + p['flange_t']*nn, Qc + t0*b + p['flange_t']*nn - 0.6*mm]
            fl = prism(pts, 0, H - 12.0)
            w = w.fuse(fl)
            for zb in p['flange_bolt_z']:
                c = np.r_[Qc + 10.0*b, zb]
                w = w.cut(cyl(p['m3_clear']/2, p['flange_t'] + 2, c - np.r_[nn, 0]*1.0, np.r_[nn, 0]))
                if nut:
                    w = w.cut(hexpocket(c + np.r_[nn, 0]*(p['flange_t'] - p['nut_depth']), np.r_[nn, 0], p['nut_depth'] + 1, p['nut_af']))
        def wplane(depth):
            return cq.Plane(origin=V(np.r_[P[i], 0] - depth*O3[i]), xDir=V(E3[i]), normal=V(O3[i]))
        # --- stake tabs on long walls
        if i in (2, 4):
            for f in (0.22, 0.78):
                sc = f * L_edge[i]
                tab = cq.Workplane(wplane(T - 0.5)).center(sc, 2.0).rect(22, 4).extrude(-18.5).val()
                tab = tab.cut(cq.Workplane(cq.Plane(origin=V(np.r_[P[i], 0] + sc*E3[i] + (T + 11)*M3[i]), xDir=V(E3[i]), normal=V(Z3)))
                              .circle(p['stake_hole']/2).extrude(10, both=True).val())
                w = w.fuse(tab)
        # --- cord notch + plug screw boss on head/foot walls
        if i in (0, 3):
            sc = L_edge[i] / 2; nw, nh = p['notch_w'], p['notch_h']
            notch = (cq.Workplane(wplane(-1)).center(sc, 0).moveTo(-nw/2, -1).lineTo(nw/2, -1)
                     .lineTo(nw/2, nh - nw/2).threePointArc((0, nh), (-nw/2, nh - nw/2)).close().extrude(-(T + 2)).val())
            w = w.cut(notch)
            zs = nh + 8
            boss = cq.Workplane(wplane(T - 0.5)).center(sc, zs).circle(5.0).extrude(-4.5).val()
            w = w.fuse(boss)
            w = w.cut(cq.Workplane(wplane(-1)).center(sc, zs).circle(p['m3_clear']/2).extrude(-(T + 6)).val())
            w = w.cut(hexpocket(np.r_[P[i], 0] + sc*E3[i] + zs*Z3 + (T + 4 - p['nut_depth'])*M3[i], M3[i], 5, p['nut_af']))
        walls.append(w)

    # --- hinge: groove all walls, add cradles + stay lug on wall 2
    walls = [cut_many(w, groove) for w in walls]
    w2 = walls[2]
    sa = math.radians(p['slot_angle']); dvec = np.array([math.cos(sa), math.sin(sa)])
    perp = np.array([-dvec[1], dvec[0]])
    r_slot = p['slot_w']/2
    arc = [(r_slot*math.cos(a), r_slot*math.sin(a)) for a in np.linspace(0, 2*math.pi, 40)]
    def slot_poly():
        back = [tuple(r_slot*np.array([math.cos(a), math.sin(a)])) for a in
                np.linspace(math.atan2(perp[1], perp[0]), math.atan2(perp[1], perp[0]) + math.pi, 24)]
        far = [tuple(dvec*14 - perp*r_slot), tuple(dvec*14 + perp*r_slot)]
        return back + far
    for a, b in cradles:
        k = hcyl(rk, a, b)
        # web down into the wall so the knuckle is solid with it
        ui = -(rk - T)   # wall inner face in hinge coords
        web = hpoly([(ui - 0.6, 0.0), (3.0, -3.0), (ui - 0.6, -rk - 10)], a, b)
        w2 = w2.fuse(k).fuse(web)
        w2 = w2.cut(hpoly(slot_poly(), a - 1, b + 1))
    # stay wall lug
    s_bar0 = p['stay_s'] + p['stay_t'] + 0.6
    s_wl = s_bar0 + p['stay_t'] + 0.6
    wu, wv = p['wall_pin_uv']
    lug = hpoly([(-2.3 - 1.0, wv + 12), (wu + 7, wv + 12), (wu + 7, wv - 8), (wu, wv - 12), (-2.3 - 1.0, wv - 26)], s_wl, s_wl + p['stay_t'])
    w2 = w2.fuse(lug)
    bl = 0.6 + p['stay_t'] + 0.3
    w2 = w2.fuse(cyl(p['stay_boss_d']/2, bl + 0.5, hpt(s_wl - bl, wu, wv), he))
    w2 = w2.cut(cyl(p['m3_clear']/2, 30, hpt(s_wl - 15, wu, wv), he))
    w2 = w2.cut(hexpocket(hpt(s_wl + p['stay_t'] - p['stay_nut_depth'], wu, wv), he, 5, p['nut_af']))
    walls[2] = w2

    # ------------------------------------------------------------- lid
    Ht = H + p['lid_t']
    lid = prism(LP, H, Ht)
    strip = prism([A0[:2] + S_MIN*E[2] - 60*M[2], A0[:2] + S_MAX*E[2] - 60*M[2],
                   A0[:2] + S_MAX*E[2], A0[:2] + S_MIN*E[2]], H - 30, Ht + 5)
    lid = lid.cut(strip)
    lp_prism = prism(LP, az - 20, az + 20)
    barrel = hcyl(rk, S_MIN, S_MAX).intersect(lp_prism)
    lid = lid.fuse(barrel)
    # skirt
    sk_out = offset_poly(P, T + p['skirt_clear']); sk_in = offset_poly(P, T + p['skirt_clear'] + p['skirt_t'])
    skirt = prism(sk_out, H - p['skirt_depth'], H + 0.5).cut(prism(sk_in, H - p['skirt_depth'] - 1, H + 1))
    keep = prism([A0[:2] + S_MIN*E[2] + 12*M[2], A0[:2] + S_MAX*E[2] + 12*M[2],
                  A0[:2] + S_MAX*E[2] + 400*M[2], A0[:2] + S_MIN*E[2] + 400*M[2]], H - 40, Ht)
    lid = lid.fuse(skirt.intersect(keep))
    # cradle clearances + pins
    for a, b in cradles:
        lid = lid.cut(hcyl(rk + hc, a - gap, b + gap))
        lid = lid.fuse(hcyl(p['pin_r'], a - gap - 0.5, a + 3.5)).fuse(hcyl(p['pin_r'], b - 3.5, b + gap + 0.5))
    # stay lid lug
    lu, lv = p['lid_pin_uv']; s_ll = p['stay_s']
    llug = hpoly([(lu - 8, 1.5), (lu + 8, 1.5), (lu + 8, lv - 2), (lu + 3, lv - 7), (lu - 3, lv - 7), (lu - 8, lv - 2)], s_ll, s_ll + p['stay_t'])
    lid = lid.fuse(llug)
    bl = 0.6 + p['stay_t'] + 0.3
    lid = lid.fuse(cyl(p['stay_boss_d']/2, bl + 0.5, hpt(s_ll + p['stay_t'] - 0.5, lu, lv), he))
    lid = lid.cut(cyl(p['m3_clear']/2, 30, hpt(s_ll - 10, lu, lv), he))
    lid = lid.cut(hexpocket(hpt(s_ll + p['stay_nut_depth'], lu, lv), -he, 5, p['nut_af']))

    # ------------------------------------------------------------- texture
    if texture:
        print("  texturing walls ...")
        seeds = [11, 23, 37, 41, 53, 67]
        for i in range(n):
            pl = cq.Plane(origin=V(np.r_[P[i], 0]), xDir=V(E3[i]), normal=V(O3[i]))
            top = H - (9 if i == 2 else 3)
            g, deep = grain_pattern(2.5, L_edge[i] - 2.5, 3, top, 6.5, seeds[i], knots=2 if L_edge[i] > 200 else 1,
                                    seams=(H/3, 2*H/3))
            tools = poly_solid(pl, g, -p['grain_depth']) + (poly_solid(pl, deep, -p['seam_depth']) if deep is not None else [])
            walls[i] = cut_many(walls[i], tools)
        print("  texturing lid ...")
        # lid local 2D frame: x' = y (along coffin), y' = x  (lines run head->foot)
        lid_poly = Polygon([(q[1], q[0]) for q in LP]).buffer(-3.0, join_style=2)
        # exclude the barrel band near the hinge
        su = lambda X, Y: to_su((Y, X))
        band = Polygon([(pt[1], pt[0]) for pt in [A0[:2] + S_MIN*E[2] + 7.5*M[2], A0[:2] + S_MAX*E[2] + 7.5*M[2],
                                                  A0[:2] + S_MAX*E[2] - 50*M[2], A0[:2] + S_MIN*E[2] - 50*M[2]]])
        region = lid_poly.difference(band)
        cross = None
        if p['lid_cross']:
            cx, cy0, cy1, arm_y, arm_w, bw = 0.0, 38.0, 196.0, 82.0, 96.0, 20.0
            cross = unary_union([box(cy0, cx - bw/2, cy1, cx + bw/2), box(arm_y - bw/2, cx - arm_w/2, arm_y + bw/2, cx + arm_w/2)])
        xs = [q[1] for q in LP]; ys = [q[0] for q in LP]
        avoid = cross.buffer(3.0, join_style=2) if cross is not None else None
        g, deep = grain_pattern(min(xs), max(xs), min(ys), max(ys), 7.0, 99, knots=4,
                                seams=(-p['shoulder_w']/6, p['shoulder_w']/6), clip=region, avoid=avoid)
        # lid top plane: x'=y, y'=x -> xDir=(0,1,0); normal must be such that yDir=(1,0,0): normal = x' cross y' = (0,0,-1)
        # use normal -z then extrude POSITIVE depth goes down (into lid)
        plid = cq.Plane(origin=V((0, 0, Ht)), xDir=V((0, 1, 0)), normal=V((0, 0, -1)))
        tools = poly_solid(plid, g, p['grain_depth'])
        if deep is not None:
            deep = deep.difference(avoid) if avoid is not None else deep
            tools += poly_solid(plid, deep, p['seam_depth'])
        if cross is not None:
            outline_g = cross.exterior.buffer(0.9, join_style=2).intersection(region)
            inner = cross.buffer(-4.0, join_style=2).exterior.buffer(0.6, join_style=2)
            tools += poly_solid(plid, outline_g, 1.2) + poly_solid(plid, inner, 0.6)
        lid = cut_many(lid, tools)

    # ------------------------------------------------------------- split lid + splice ribs
    Ys = p['lid_split_y']
    big = 2000
    head_box = prism([(-big, -big), (big, -big), (big, Ys), (-big, Ys)], 0, 400)
    foot_box = prism([(-big, Ys), (big, Ys), (big, big), (-big, big)], 0, 400)
    lid_head = lid.intersect(head_box); lid_foot = lid.cut(head_box)
    rib_region = Polygon(offset_poly(P, T + p['skirt_clear'] + 0.5)).intersection(
        Polygon([tuple(A0[:2] + S_MIN*E[2] + 16*M[2]), tuple(A0[:2] + S_MAX*E[2] + 16*M[2]),
                 tuple(A0[:2] + S_MAX*E[2] + 400*M[2]), tuple(A0[:2] + S_MIN*E[2] + 400*M[2])]))
    ribpoly_h = rib_region.intersection(box(-big, Ys - p['flange_t'], big, Ys))
    ribpoly_f = rib_region.intersection(box(-big, Ys, big, Ys + p['flange_t']))
    rib_h = prism(list(ribpoly_h.exterior.coords)[:-1], H - 12, H + 0.5)
    rib_f = prism(list(ribpoly_f.exterior.coords)[:-1], H - 12, H + 0.5)
    lid_head = lid_head.fuse(rib_h); lid_foot = lid_foot.fuse(rib_f)
    xb0, _, xb1, _ = ribpoly_h.bounds
    for xh in np.linspace(xb0 + 14, xb1 - 10, 3):
        c = np.array([xh, Ys - p['flange_t'] - 1, H - 5.5])
        hole = cyl(p['m3_clear']/2, 2*p['flange_t'] + 2, c, (0, 1, 0))
        lid_head = lid_head.cut(hole); lid_foot = lid_foot.cut(hole)
        lid_head = lid_head.cut(hexpocket(c, (0, 1, 0), p['nut_depth'] + 1.0, p['nut_af']))

    # ------------------------------------------------------------- small parts (built in print orientation)
    # stay bar
    def L_of(phi):
        c, s_ = math.cos(math.radians(phi)), math.sin(math.radians(phi))
        return np.array([lu*c - lv*s_, lu*s_ + lv*c])
    W = np.array(p['wall_pin_uv'])
    dist = lambda phi: np.linalg.norm(W - L_of(phi))
    d0, dmax = dist(0), dist(p['stay_max_deg'])
    rs = p['stay_boss_d']/2 + 0.3; bw = p['stay_w']; bt = p['stay_t']
    bar2d = LineString([(0, 0), (dmax, 0)]).buffer(bw/2, 48)
    slot2d = LineString([(d0 - 2, 0), (dmax, 0)]).buffer(rs, 32)
    notches = []
    for ang in p['stay_notch_deg']:
        dk = dist(ang)
        # pocket on -y side: steep face toward pivot, ramp away from it
        notches.append(Point(dk, -2.4).buffer(rs, 32))
        notches.append(Polygon([(dk - rs, -2.4), (dk - rs, 0), (dk + rs + 6, 0), (dk + rs, -2.4)]))
    bar2d = bar2d.difference(slot2d).difference(unary_union(notches)).difference(Point(0, 0).buffer(rs, 32))
    stay = poly_solid(cq.Plane.XY(), bar2d, bt)[0]
    washer = poly_solid(cq.Plane.XY(), Point(0, 0).buffer(5.5, 48).difference(Point(0, 0).buffer(p['m3_clear']/2, 32)), p['stay_washer_t'])[0]
    # plugs (print orientation: flange on bed, X along wall, Y = height, Z out of bed into wall)
    def make_plug(slotted):
        nw, nh, cl, fl, ft = p['notch_w'], p['notch_h'], 0.35, p['plug_flange'], p['plug_flange_t']
        tongue = (Polygon([(-nw/2 + cl, 0), (nw/2 - cl, 0), (nw/2 - cl, nh - nw/2)])
                  .union(box(-nw/2 + cl, 0, nw/2 - cl, nh - nw/2)).union(Point(0, nh - nw/2).buffer(nw/2 - cl, 48)))
        tongue = tongue.intersection(box(-100, 0, 100, 100))
        flange = (box(-nw/2 - fl, 0, nw/2 + fl, nh - nw/2).union(Point(0, nh - nw/2).buffer(nw/2 + fl, 48))
                  .union(Point(0, nh + 8).buffer(6.5, 48)).convex_hull).intersection(box(-100, 0, 100, 100))
        hole = Point(0, nh + 8).buffer(p['m3_clear']/2, 32)
        flange = flange.difference(hole)
        if slotted:
            sl = box(-p['cord_slot_w']/2, -1, p['cord_slot_w']/2, 9.0).union(Point(0, 9.0).buffer(p['cord_slot_w']/2, 32))
            tongue = tongue.difference(sl); flange = flange.difference(sl)
        s1 = poly_solid(cq.Plane.XY(), flange, ft)[0]
        s2 = poly_solid(cq.Plane(origin=(0, 0, ft), xDir=(1, 0, 0), normal=(0, 0, 1)), tongue, T - 0.3)[0]
        return s1.fuse(s2)
    plug_blank, plug_cord = make_plug(False), make_plug(True)

    info = dict(P=P, Q=Q, E=E, M=M, A0=A0, he=he, hm=hm, az=az, L_edge=L_edge, d0=d0, dmax=dmax,
                s_bar0=s_bar0, LP=LP, H=H, Ht=Ht, L_of=L_of, W=W)
    parts = dict(wall_0_head=walls[0], wall_1_right_shoulder=walls[1], wall_2_right_long_hinge=walls[2],
                 wall_3_foot=walls[3], wall_4_left_long=walls[4], wall_5_left_shoulder=walls[5],
                 lid_head_half=lid_head, lid_foot_half=lid_foot)
    small = dict(lid_stay=stay, stay_washer_a=washer, stay_washer_b=washer.translate(cq.Vector(14, 0, 0)), cord_plug_blank=plug_blank, cord_plug_slotted=plug_cord)
    return parts, small, info

# ------------------------------------------------------------------ orientation helpers
def wall_print(shape, info, i):
    P, E, M = info['P'], info['E'], info['M']
    e3 = np.r_[E[i], 0]; m3 = np.r_[M[i], 0]; z = np.array([0, 0, 1.0])
    return to_local(shape, np.r_[P[i], 0], -e3, z, m3)

def lid_print(shape, info):
    # rotate 180 deg about X: top face onto the bed
    return to_local(shape, np.array([0, 0, info['Ht']]), np.array([1, 0, 0.]), np.array([0, -1, 0.]), np.array([0, 0, -1.]))

def stay_world(shape, info, phi, p=PARAMS):
    L = info['L_of'](phi); W = info['W']
    x = unit(W - L); y = np.array([-x[1], x[0]])
    hm, he, A0 = info['hm'], info['he'], info['A0']
    org = A0 + info['s_bar0']*he + L[0]*hm + L[1]*np.array([0, 0, 1.])
    x3 = x[0]*hm + x[1]*np.array([0, 0, 1.]); y3 = y[0]*hm + y[1]*np.array([0, 0, 1.])
    return place(shape, org, x3, y3, he)

def plug_world(shape, info, wall_i, p=PARAMS):
    P, E, M = info['P'], info['E'], info['M']
    e3 = np.r_[E[wall_i], 0]; m3 = np.r_[M[wall_i], 0]
    org = np.r_[P[wall_i], 0] + info['L_edge'][wall_i]/2*e3 - p['plug_flange_t']*m3
    return place(shape, org, -e3, np.array([0, 0, 1.]), m3)

def lid_rotate(shape, info, phi):
    A0, he = info['A0'], info['he']
    return shape.rotate(V(A0), V(A0 + he), phi)

if __name__ == "__main__":
    texture = "--no-texture" not in sys.argv
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else "out"
    os.makedirs(out, exist_ok=True)
    parts, small, info = build(texture=texture)
    print("built")

# ------------------------------------------------------------------ anti-curl corner tabs
TAB = dict(len=24.0, out=10.0, overlap=2.0, t=0.5, r=4.0)   # t = one 0.5 mm first layer

def add_corner_tabs(shape, tab=TAB):
    """Print-oriented wall (outer face on z=0, length along X): add one-layer tabs at the 4 footprint
    corners. Tabs stay inside the part's X range (so long walls still fit 320 mm) and stick out in Y."""
    bb = shape.BoundingBox(); L = bb.xmax - bb.xmin
    tabs = []
    for x0 in (bb.xmin, bb.xmax - tab['len']):
        slab = cq.Solid.makeBox(tab['len'], bb.ylen + 40, 0.3, cq.Vector(x0, bb.ymin - 20, 0))
        fp = shape.intersect(slab).BoundingBox()          # first-layer footprint in this window
        for y_edge, sgn in ((fp.ymin, -1), (fp.ymax, +1)):
            y_in, y_out = y_edge - sgn*tab['overlap'], y_edge + sgn*tab['out']
            y0, y1 = min(y_in, y_out), max(y_in, y_out)
            t = (cq.Workplane("XY").center(x0 + tab['len']/2, (y0 + y1)/2)
                 .rect(tab['len'], y1 - y0).extrude(tab['t']).edges("|Z").fillet(tab['r']).val())
            tabs.append(t)
    return shape.fuse(*tabs).clean()
