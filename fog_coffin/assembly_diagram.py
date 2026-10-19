"""Top-view assembly map for the fog coffin: walls, plates, corner joints (nut side vs screw-head side),
hinge cradles, lid stay, cord notches, stake tabs and the lid seam. Writes renders/assembly_map.png."""
import os, sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MPoly, RegularPolygon, Circle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coffin

p = coffin.PARAMS
P = coffin.outline(p); Q = coffin.offset_poly(P, p['wall']); E, M = coffin.edge_frames(P); n = 6
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "renders", "assembly_map.png")

INK, MUTED, WALL, FLANGE, NUT, HEAD, ACCENT = "#1f1f1f", "#6b6b6b", "#c9b9a6", "#a8927a", "#b3261e", "#1f5fa8", "#2e7d32"
walls = ["Head wall\nplate 3", "Right shoulder\nplate 5", "Hinge wall (right long)\nplate 1",
         "Foot wall\nplate 4", "Left long wall\nplate 2", "Left shoulder\nplate 5"]
letters = "ABCDEF"

fig, ax = plt.subplots(figsize=(8.2, 11))
ax.set_aspect("equal"); ax.axis("off")

for i in range(n):  # wall bands
    j = (i + 1) % n
    ax.add_patch(MPoly([P[i], P[j], Q[j], Q[i]], closed=True, fc=WALL, ec=INK, lw=1.2, zorder=2))

for i in range(n):  # corner flanges + screws
    j = (i + 1) % n
    for corner, side_dir, nut in ((i, E[i], True), (j, -E[i], False)):
        Pc, Qc = P[corner], Q[corner]
        b = coffin.unit(Qc - Pc); nn = coffin.unit(side_dir - (side_dir @ b) * b)
        half = np.arccos(np.clip(b @ side_dir, -1, 1)); t0 = p['flange_t'] / np.tan(half)
        D = p['flange_depth']
        pts = [Qc, Qc + D * b, Qc + D * b + p['flange_t'] * nn, Qc + t0 * b + p['flange_t'] * nn]
        ax.add_patch(MPoly(pts, closed=True, fc=FLANGE, ec=INK, lw=0.8, zorder=3))
        c = Qc + 10.0 * b + (p['flange_t'] + 1.6) * nn
        if nut:
            ax.add_patch(RegularPolygon(c, 6, radius=3.4, fc=NUT, ec="white", lw=0.6, zorder=5))
        else:
            ax.add_patch(Circle(c, 3.0, fc=HEAD, ec="white", lw=0.6, zorder=5))
    # wall label
    mid = (P[i] + P[j]) / 2; out = -M[i]
    rot = np.degrees(np.arctan2(E[i][1], E[i][0]))
    if rot > 90: rot -= 180
    if rot < -90: rot += 180
    off = 34 if i in (0, 3) else 46
    ax.text(*(mid + out * off), walls[i], ha="center", va="center", rotation=rot, fontsize=10, color=INK,
            fontweight="bold" if i == 2 else "normal", zorder=6)

for k in range(n):  # corner letters
    v = P[k] - np.mean(P, axis=0); v = v / np.linalg.norm(v)
    ax.add_patch(Circle(P[k] + v * 16, 8.5, fc="white", ec=INK, lw=1.2, zorder=6))
    ax.text(*(P[k] + v * 16), letters[k], ha="center", va="center", fontsize=11, fontweight="bold", zorder=7)

# hinge cradles + lid stay on the hinge wall
A0 = P[2]; e2, m2 = E[2], M[2]
for s in p['cradle_s']:
    c = A0 + s * e2 + 1.5 * m2
    ax.add_patch(MPoly([c - 4 * e2 - 3 * m2, c + 4 * e2 - 3 * m2, c + 4 * e2 + 6 * m2, c - 4 * e2 + 6 * m2],
                       closed=True, fc=ACCENT, ec=INK, lw=0.8, zorder=4))
c_st = A0 + (p['stay_s'] + 8) * e2 + 18 * m2
ax.add_patch(MPoly([c_st - 9 * e2 - 4 * m2, c_st + 9 * e2 - 4 * m2, c_st + 9 * e2 + 4 * m2, c_st - 9 * e2 + 4 * m2],
                   closed=True, fc="#e0a800", ec=INK, lw=0.8, zorder=4))
ax.annotate("lid stay", c_st, xytext=(c_st + m2 * 62 - e2 * 6), fontsize=9, color=INK,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8), ha="center", zorder=6)
ax.annotate("3 hinge cradles", A0 + p['cradle_s'][1] * e2 + 4 * m2, xytext=(A0 + p['cradle_s'][1] * e2 + 70 * m2),
            fontsize=9, color=ACCENT, ha="center", arrowprops=dict(arrowstyle="-", color=ACCENT, lw=0.8), zorder=6)

# stake tabs on both long walls
for wi in (2, 4):
    L = np.linalg.norm(P[(wi + 1) % n] - P[wi])
    for f in (0.22, 0.78):
        c = P[wi] + f * L * E[wi] + (p['wall'] + 9) * M[wi]
        ax.add_patch(MPoly([c - 11 * E[wi] - 9 * M[wi], c + 11 * E[wi] - 9 * M[wi], c + 11 * E[wi] + 9 * M[wi],
                            c - 11 * E[wi] + 9 * M[wi]], closed=True, fc="white", ec=MUTED, lw=0.8, ls="--", zorder=3))
        ax.add_patch(Circle(c, p['stake_hole'] / 2, fc="white", ec=MUTED, lw=0.8, zorder=4))

# cord notches (head + foot)
for wi, lab in ((0, "cord notch + plug"), (3, "cord notch + plug")):
    mid = (P[wi] + P[(wi + 1) % n]) / 2
    ax.plot(*np.array([mid - 15 * E[wi], mid + 15 * E[wi]]).T, color=INK, lw=4.5, solid_capstyle="butt", zorder=5)
    ax.plot(*np.array([mid - 15 * E[wi], mid + 15 * E[wi]]).T, color="white", lw=2.5, solid_capstyle="butt", zorder=6)
    ax.text(*(mid + M[wi] * 18), lab, ha="center", va="center", fontsize=8.5, color=MUTED, zorder=6)

# lid seam (lid shown as outline only)
LP = p and coffin.offset_poly(P, [-p['lid_overhang']] * 2 + [0] + [-p['lid_overhang']] * 3)
ax.add_patch(MPoly(LP, closed=True, fill=False, ec=MUTED, lw=0.8, ls=(0, (2, 3)), zorder=1))
y = p['lid_split_y']
ax.plot([-150, 150], [y, y], color=MUTED, lw=1.0, ls=(0, (6, 3)), zorder=1)
ax.text(0, y + 7, "lid seam: head half | foot half (3 screws, from underneath)", ha="center", fontsize=8.5,
        color=MUTED, zorder=6)
ax.text(0, 58, "INSIDE\n(no bottom)", ha="center", va="center", fontsize=9, color=MUTED)

# legend
lx, ly = -205, -78
ax.add_patch(RegularPolygon((lx, ly), 6, radius=4, fc=NUT, ec="white")); ax.text(lx + 9, ly, "nut goes in this flange's hex pocket", va="center", fontsize=9)
ax.add_patch(Circle((lx, ly - 15), 3.4, fc=HEAD, ec="white")); ax.text(lx + 9, ly - 15, "screw head on this flange (plain hole)", va="center", fontsize=9)
ax.text(lx - 6, ly - 33, "Each corner A–F: 3 × M3×12 at 20, 85 and 150 mm above the ground.", fontsize=9, color=INK)
ax.text(lx - 6, ly - 47, "Dashed boxes: stake tabs (inside, at ground level).", fontsize=9, color=MUTED)

ax.set_title("Fog coffin — assembly map (top view, looking down; head end at the bottom)", fontsize=12, pad=12)
ax.set_xlim(-215, 215); ax.set_ylim(-140, 470)
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=150, facecolor="white")
print("wrote", OUT)
