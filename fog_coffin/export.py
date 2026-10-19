import os, coffin, cadquery as cq, numpy as np, trimesh, json
os.chdir(os.path.dirname(os.path.abspath(__file__)))  # outputs go next to this script
parts, small, info = coffin.build(texture=True)
os.makedirs("stl", exist_ok=True); os.makedirs("plates", exist_ok=True)
printed = {}
for i,k in enumerate([k for k in parts if k.startswith('wall')]):
    printed[k] = coffin.add_corner_tabs(coffin.wall_print(parts[k], info, i))
for k in ('lid_head_half','lid_foot_half'):
    printed[k] = coffin.lid_print(parts[k], info)
printed.update(small)
for k,s in printed.items():
    assert s.isValid(), k
    bb = s.BoundingBox(); s = s.translate(cq.Vector(-bb.xmin, -bb.ymin, -bb.zmin)); printed[k] = s
    cq.exporters.export(s, f"stl/{k}.stl", tolerance=0.04, angularTolerance=0.15)
# assembly-position STEP (closed, untextured) for Fusion
up, us, ui = coffin.build(texture=False)
os.makedirs("step_untextured", exist_ok=True)
world = dict(up)
world['lid_stay'] = coffin.stay_world(us['lid_stay'], ui, 0)
world['cord_plug_slotted'] = coffin.plug_world(us['cord_plug_slotted'], ui, 0)
world['cord_plug_blank'] = coffin.plug_world(us['cord_plug_blank'], ui, 3)
asm = cq.Assembly(name="fog_coffin")
for k,s_ in world.items():
    cq.exporters.export(s_, f"step_untextured/{k}.step"); asm.add(s_, name=k)
asm.export("step_untextured/fog_coffin_assembly.step")
# plates: Bambu Studio projects with the tuned slicer settings from bambu/ (see bambu3mf.py)
import bambu3mf
report = bambu3mf.build_all()
vol = {k: round(v.Volume()/1000,1) for k,v in printed.items()}
json.dump(dict(plates=report, volume_cm3=vol, total_cm3=round(sum(vol.values()),0)), open("report.json","w"), indent=1)
print(json.dumps(report, indent=1)); print(vol, sum(vol.values()))
