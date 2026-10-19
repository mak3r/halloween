import os, coffin, cadquery as cq, numpy as np, trimesh, json
parts, small, info = coffin.build(texture=True)
os.makedirs("out/stl", exist_ok=True); os.makedirs("out/plates", exist_ok=True)
printed = {}
for i,k in enumerate([k for k in parts if k.startswith('wall')]):
    printed[k] = coffin.add_corner_tabs(coffin.wall_print(parts[k], info, i))
for k in ('lid_head_half','lid_foot_half'):
    printed[k] = coffin.lid_print(parts[k], info)
printed.update(small)
for k,s in printed.items():
    assert s.isValid(), k
    bb = s.BoundingBox(); s = s.translate(cq.Vector(-bb.xmin, -bb.ymin, -bb.zmin)); printed[k] = s
    cq.exporters.export(s, f"out/stl/{k}.stl", tolerance=0.04, angularTolerance=0.15)
# assembly-position STEP (closed, untextured) for Fusion
up, us, ui = coffin.build(texture=False)
os.makedirs("out/step_untextured", exist_ok=True)
world = dict(up)
world['lid_stay'] = coffin.stay_world(us['lid_stay'], ui, 0)
world['cord_plug_slotted'] = coffin.plug_world(us['cord_plug_slotted'], ui, 0)
world['cord_plug_blank'] = coffin.plug_world(us['cord_plug_blank'], ui, 3)
asm = cq.Assembly(name="fog_coffin")
for k,s_ in world.items():
    cq.exporters.export(s_, f"out/step_untextured/{k}.step"); asm.add(s_, name=k)
asm.export("out/step_untextured/fog_coffin_assembly.step")
# plates (H2C left nozzle 325 x 320) -> 3MF scenes
plates = {
 "plate1_right_long_hinge_wall": ["wall_2_right_long_hinge"],
 "plate2_left_long_wall": ["wall_4_left_long"],
 "plate3_head_wall_and_plugs": ["wall_0_head", "cord_plug_blank", "cord_plug_slotted"],
 "plate4_foot_wall_and_stay": ["wall_3_foot", "lid_stay", "stay_washer_a", "stay_washer_b"],
 "plate5_shoulder_walls": ["wall_1_right_shoulder", "wall_5_left_shoulder"],
 "plate6_lid_head_half": ["lid_head_half"],
 "plate7_lid_foot_half": ["lid_foot_half"],
}
report = {}
for pname, names in plates.items():
    sc = trimesh.Scene(); x = 0.0; y = 0.0; rowh = 0; placed = []
    meshes = [trimesh.load(f"out/stl/{n}.stl") for n in names]
    # keep every plate inside 300 mm in X (area reachable by either H2C nozzle): turn long parts to run along Y
    for m in meshes:
        if m.extents[0] > 300 and m.extents[1] <= 320:
            m.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2, [0, 0, 1])); m.apply_translation(-m.bounds[0])
    # simple shelf packing: big part first, small ones to its right stacked
    big = meshes[0]; bx, by = big.extents[0], big.extents[1]
    offs = [(0, 0)]; cx = bx + 10; cy = 0
    for m in meshes[1:]:
        if cy + m.extents[1] > 320: cx += max(mm.extents[0] for mm in meshes[1:]) + 10; cy = 0
        offs.append((cx, cy)); cy += m.extents[1] + 10
    W = max(o[0] + m.extents[0] for o, m in zip(offs, meshes)); Hh = max(o[1] + m.extents[1] for o, m in zip(offs, meshes))
    ox, oy = (300 - W) / 2, (320 - Hh) / 2
    assert W <= 300 and Hh <= 320, (pname, W, Hh)
    for n, m, (dx, dy) in zip(names, meshes, offs):
        m = m.copy(); m.apply_translation([ox + dx, oy + dy, 0]); sc.add_geometry(m, node_name=n, geom_name=n)
    sc.export(f"out/plates/{pname}.3mf")
    report[pname] = dict(footprint=[round(W,1), round(Hh,1)], fits=bool(W <= 300 and Hh <= 320), parts=names,
                         height=round(max(m.extents[2] for m in meshes),1))
vol = {k: round(v.Volume()/1000,1) for k,v in printed.items()}
json.dump(dict(plates=report, volume_cm3=vol, total_cm3=round(sum(vol.values()),0)), open("out/report.json","w"), indent=1)
print(json.dumps(report, indent=1)); print(vol, sum(vol.values()))
