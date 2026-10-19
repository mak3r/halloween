"""Hinge-cradle fit test: slices the real model around the middle cradle.
Makes one lid coupon (pins) and wall coupons in several slot widths, labelled on top."""
import os, copy, numpy as np, cadquery as cq, coffin

VARIANTS = (5.0, 5.1, 5.2, 5.3)     # round 2: all of 5.4/5.6/5.8 were sloppy (pin is 5.0 nominal)
GAP = 0.4                           # axial gap pin-barrel to cradle (was 0.6)
C = coffin.PARAMS['cradle_s'][1]    # middle cradle (on the foot lid half)
OUT = "out/hinge_test_r2"; os.makedirs(OUT, exist_ok=True)

def hbox(info, s0, s1, u0, u1, v0, v1):
    A0, he, hm = info['A0'], info['he'], info['hm']; z = np.array([0, 0, 1.0])
    pl = cq.Plane(origin=coffin.V(A0 + s0*he), xDir=coffin.V(hm), normal=coffin.V(he))
    return cq.Workplane(pl).polyline([(u0, v0), (u1, v0), (u1, v1), (u0, v1)]).close().extrude(s1 - s0).val()

def label(shape, text, z_top):
    bb = shape.BoundingBox()
    t = (cq.Workplane("XY").workplane(offset=z_top).center((bb.xmin + bb.xmax)/2, bb.ymin + 7)
         .text(text, 6, -0.6, kind="bold", halign="center", valign="center").val())
    return shape.cut(t)

def zero(s):
    bb = s.BoundingBox(); return s.translate(cq.Vector(-bb.xmin, -bb.ymin, -bb.zmin))

out = {}
for w in VARIANTS:
    p = copy.deepcopy(coffin.PARAMS); p['slot_w'] = w; p['hinge_gap'] = GAP
    parts, small, info = coffin.build(p, texture=False)
    wall = parts['wall_2_right_long_hinge'].intersect(hbox(info, C - 20, C + 20, -12, 30, -32, 12))
    wall = zero(coffin.wall_print(wall, info, 2))
    wall = label(wall, f"{w:.1f}", p['wall'])
    out[f"cradle_wall_slot_{w:.1f}"] = wall
    if w == VARIANTS[0]:
        lid = parts['lid_foot_half'].intersect(hbox(info, C - 20, C + 20, -12, 40, -12, 12))
        lid = zero(coffin.lid_print(lid, info))
        out["cradle_lid_pins"] = label(lid, "PIN", p['lid_t'])
        # sanity: nominal pair, lid swung 0..85 deg, no overlap
        wl = parts['wall_2_right_long_hinge'].intersect(hbox(info, C - 20, C + 20, -12, 30, -32, 12))
        ld = parts['lid_foot_half'].intersect(hbox(info, C - 20, C + 20, -12, 40, -12, 12))
        for phi in (0, 30, 60, 85):
            v = coffin.lid_rotate(ld, info, phi).intersect(wl).Volume()
            print(f"  swing {phi:2d} deg overlap {v:.3f} mm^3")

# plate: all coupons side by side
x = 0; asm = []
for k, s in out.items():
    assert s.isValid(), k
    cq.exporters.export(s, f"{OUT}/{k}.stl", tolerance=0.02, angularTolerance=0.1)
    bb = s.BoundingBox(); print(f"{k:24s} {bb.xlen:5.1f} x {bb.ylen:5.1f} x {bb.zlen:5.1f}  vol {s.Volume()/1000:.1f} cm3")
    asm.append(s.translate(cq.Vector(x, 0, 0))); x += bb.xlen + 8
plate = cq.Compound.makeCompound(asm)
cq.exporters.export(plate, f"{OUT}/hinge_test_r2_plate.stl", tolerance=0.02, angularTolerance=0.1)
