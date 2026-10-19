import coffin, numpy as np, itertools
parts, small, info = coffin.build(texture=False)
walls = [parts[k] for k in parts if k.startswith('wall')]
lids = [parts['lid_head_half'], parts['lid_foot_half']]
def iv(a,b):
    try: return a.intersect(b).Volume()
    except Exception as e: return -1
print("== print-orientation bounding boxes (H2C left nozzle 325x320x325)")
for i,k in enumerate([k for k in parts if k.startswith('wall')]):
    s = coffin.wall_print(parts[k], info, i); bb=s.BoundingBox()
    print(f"  {k:28s} {bb.xlen:6.1f} x {bb.ylen:6.1f} x {bb.zlen:5.1f}  zmin={bb.zmin:.2f}")
for k in ('lid_head_half','lid_foot_half'):
    s = coffin.lid_print(parts[k], info); bb=s.BoundingBox()
    print(f"  {k:28s} {bb.xlen:6.1f} x {bb.ylen:6.1f} x {bb.zlen:5.1f}  zmin={bb.zmin:.2f}")
print("== wall/wall interference")
for a,b in itertools.combinations(range(6),2):
    v=iv(walls[a],walls[b])
    if v>0.01: print("  walls",a,b,round(v,2))
print("  lid halves:", round(iv(*lids),3))
print("== lid swing vs walls (mm^3 overlap)")
for phi in [0,5,10,20,30,45,55,65,80,85]:
    rl=[coffin.lid_rotate(l,info,phi) for l in lids]
    hits=[(wi,li,round(iv(rl[li],walls[wi]),2)) for wi in range(6) for li in range(2)]
    hits=[h for h in hits if h[2]>0.01]
    print(f"  {phi:3d} deg:", hits if hits else "clear")
print("== stay bar vs walls / lid")
for phi in [0,30,55,80,85]:
    bar=coffin.stay_world(small['lid_stay'],info,phi)
    rl=[coffin.lid_rotate(l,info,phi) for l in lids]
    hw=[(wi,round(iv(bar,walls[wi]),2)) for wi in range(6)]
    hl=[(li,round(iv(bar,rl[li]),2)) for li in range(2)]
    print(f"  {phi:3d} deg: walls {[h for h in hw if h[1]>0.01]} lids {[h for h in hl if h[1]>0.01]}")
print("== plugs vs walls")
for wi in (0,3):
    pl=coffin.plug_world(small['cord_plug_blank'],info,wi)
    print("  plug on wall",wi, round(iv(pl,walls[wi]),2))
