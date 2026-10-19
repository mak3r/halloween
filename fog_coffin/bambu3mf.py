"""Build Bambu Studio project files (.3mf) for every plate, using a hand-tuned project as the template.

The slicer settings (printer, process, filament, plate type, temperatures, fan, brim, elephant foot...)
are copied unchanged from TEMPLATE; only the parts on the plate are replaced. The template itself is
never modified.

    python3 bambu3mf.py            # rebuild plates/*.3mf from stl/*.stl

No third-party packages needed. export.py calls build_all() after regenerating the STLs.
"""
import os, re, struct, uuid, zipfile
from xml.sax.saxutils import escape

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "bambu", "plate1_hinge_wall_PETG_tuned.3mf")
STL_DIR = os.path.join(HERE, "stl")
OUT_DIR = os.path.join(HERE, "plates")

# usable area that every H2C nozzle can reach, and where it sits on the plate (taken from the template:
# Bambu Studio placed the 201 x 315 mm hinge wall centred at x=165, y=160)
AREA_W, AREA_D = 300.0, 320.0
CENTER = (165.0, 160.0)
GAP = 10.0

PLATES = {
    "plate1_right_long_hinge_wall": ["wall_2_right_long_hinge"],
    "plate2_left_long_wall": ["wall_4_left_long"],
    "plate3_head_wall_and_plugs": ["wall_0_head", "cord_plug_blank", "cord_plug_slotted"],
    "plate4_foot_wall_and_stay": ["wall_3_foot", "lid_stay", "stay_washer_a", "stay_washer_b"],
    "plate5_shoulder_walls": ["wall_1_right_shoulder", "wall_5_left_shoulder"],
    "plate6_lid_head_half": ["lid_head_half"],
    "plate7_lid_foot_half": ["lid_foot_half"],
}

# entries rebuilt per plate (everything else in the template, incl. project_settings.config, is copied as-is)
REBUILT = re.compile(r"^(3D/3dmodel\.model|3D/_rels/3dmodel\.model\.rels|3D/Objects/.*|_rels/\.rels|"
                     r"Metadata/model_settings\.config|Metadata/plate_\d+\.json|Metadata/.*\.png|"
                     r"Metadata/cut_information\.xml)$")

MODEL_NS = ('xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" '
            'xmlns:BambuStudio="http://schemas.bambulab.com/package/2021" '
            'xmlns:p="http://schemas.microsoft.com/3dmanufacturing/production/2015/06" requiredextensions="p"')


def read_stl(path):
    """Binary or ASCII STL -> (vertices [(x,y,z)], triangles [(i,j,k)]) with shared vertices merged."""
    data = open(path, "rb").read()
    tris_xyz = []
    n = struct.unpack_from("<I", data, 80)[0] if len(data) >= 84 else 0
    if len(data) == 84 + 50 * n:
        for t in range(n):
            f = struct.unpack_from("<12f", data, 84 + 50 * t)
            tris_xyz.append((f[3:6], f[6:9], f[9:12]))
    else:  # ASCII
        pts = [tuple(map(float, m)) for m in re.findall(rb"vertex\s+(\S+)\s+(\S+)\s+(\S+)", data)]
        tris_xyz = [pts[i:i + 3] for i in range(0, len(pts), 3)]
    index, verts, tris = {}, [], []
    for tri in tris_xyz:
        ids = []
        for v in tri:
            key = (round(v[0], 5), round(v[1], 5), round(v[2], 5))
            if key not in index:
                index[key] = len(verts); verts.append(key)
            ids.append(index[key])
        if len(set(ids)) == 3:
            tris.append(tuple(ids))
    return verts, tris


def bbox(verts):
    xs, ys, zs = zip(*verts)
    return min(xs), min(ys), min(zs), max(xs), max(ys), max(zs)


def load_part(name):
    verts, tris = read_stl(os.path.join(STL_DIR, name + ".stl"))
    x0, y0, z0, x1, y1, z1 = bbox(verts)
    if x1 - x0 > AREA_W and y1 - y0 <= AREA_D:          # long walls run the 320 mm direction
        verts = [(-y, x, z) for x, y, z in verts]
        x0, y0, z0, x1, y1, z1 = bbox(verts)
    cx, cy, cz = (x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2
    verts = [(x - cx, y - cy, z - cz) for x, y, z in verts]  # Bambu stores meshes about their own centre
    return dict(name=name, verts=verts, tris=tris, w=x1 - x0, d=y1 - y0, h=z1 - z0)


def layout(parts):
    """First part on the left, the rest stacked in columns to its right; whole group centred."""
    big = parts[0]; spots = [(0.0, 0.0)]
    cx, cy, colw = big["w"] + GAP, 0.0, 0.0
    for p in parts[1:]:
        if cy and cy + p["d"] > AREA_D:
            cx += colw + GAP; cy = 0.0; colw = 0.0
        spots.append((cx, cy)); cy += p["d"] + GAP; colw = max(colw, p["w"])
    W = max(sx + p["w"] for (sx, _), p in zip(spots, parts))
    D = max(sy + p["d"] for (_, sy), p in zip(spots, parts))
    if W > AREA_W or D > AREA_D:
        raise ValueError(f"plate is {W:.1f} x {D:.1f} mm, over {AREA_W} x {AREA_D}")
    ox, oy = CENTER[0] - W / 2, CENTER[1] - D / 2
    return [(ox + sx + p["w"] / 2, oy + sy + p["d"] / 2, p["h"] / 2) for (sx, sy), p in zip(spots, parts)]


def object_model(inner_id, part):
    out = [f'<?xml version="1.0" encoding="UTF-8"?>\n<model unit="millimeter" xml:lang="en-US" {MODEL_NS}>',
           ' <metadata name="BambuStudio:3mfVersion">1</metadata>', ' <resources>',
           f'  <object id="{inner_id}" p:UUID="{uuid.uuid4()}" type="model">', '   <mesh>', '    <vertices>']
    out += [f'     <vertex x="{x:.6g}" y="{y:.6g}" z="{z:.6g}"/>' for x, y, z in part["verts"]]
    out += ['    </vertices>', '    <triangles>']
    out += [f'     <triangle v1="{a}" v2="{b}" v3="{c}"/>' for a, b, c in part["tris"]]
    out += ['    </triangles>', '   </mesh>', '  </object>', ' </resources>', ' <build/>', '</model>', '']
    return "\n".join(out)


def build_plate(tpl, plate_name, part_names, extruder, app_meta):
    parts = [load_part(n) for n in part_names]
    places = layout(parts)
    files = {}
    res, items, rels, objs, insts, asm = [], [], [], [], [], []
    for k, (part, (px, py, pz)) in enumerate(zip(parts, places)):
        oid, inner = 2 * (k + 1), 2 * (k + 1) - 1
        path = f"3D/Objects/object_{oid}.model"
        files[path] = object_model(inner, part)
        tr = f"1 0 0 0 1 0 0 0 1 {px:.5f} {py:.5f} {pz:.5f}"
        res.append(f'  <object id="{oid}" p:UUID="{uuid.uuid4()}" type="model">\n   <components>\n'
                   f'    <component p:path="/{path}" objectid="{inner}" p:UUID="{uuid.uuid4()}" '
                   f'transform="1 0 0 0 1 0 0 0 1 0 0 0"/>\n   </components>\n  </object>')
        items.append(f'  <item objectid="{oid}" p:UUID="{uuid.uuid4()}" transform="{tr}" printable="1"/>')
        rels.append(f' <Relationship Target="/{path}" Id="rel-{k + 1}" '
                    'Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>')
        nm = escape(part["name"]); fc = len(part["tris"])
        objs.append(f'  <object id="{oid}">\n    <metadata key="name" value="{nm}"/>\n'
                    f'    <metadata key="extruder" value="{extruder}"/>\n    <metadata face_count="{fc}"/>\n'
                    f'    <part id="1" subtype="normal_part" uuid="{uuid.uuid4()}">\n'
                    f'      <metadata key="name" value="{nm}"/>\n'
                    '      <metadata key="matrix" value="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"/>\n'
                    f'      <mesh_stat face_count="{fc}" edges_fixed="0" degenerate_facets="0" facets_removed="0" '
                    'facets_reversed="0" backwards_edges="0"/>\n    </part>\n  </object>')
        insts.append(f'    <model_instance>\n      <metadata key="object_id" value="{oid}"/>\n'
                     f'      <metadata key="instance_id" value="0"/>\n'
                     f'      <metadata key="identify_id" value="{100 + k}"/>\n    </model_instance>')
        asm.append(f'   <assemble_item object_id="{oid}" instance_id="0" transform="{tr}" offset="0 0 0" />')

    files["3D/3dmodel.model"] = ("\n".join(
        [f'<?xml version="1.0" encoding="UTF-8"?>\n<model unit="millimeter" xml:lang="en-US" {MODEL_NS}>']
        + app_meta + [' <resources>'] + res + [' </resources>', f' <build p:UUID="{uuid.uuid4()}">']
        + items + [' </build>', '</model>', '']))
    files["3D/_rels/3dmodel.model.rels"] = (
        '<?xml version="1.0" encoding="UTF-8"?>\n<Relationships '
        'xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n' + "\n".join(rels)
        + '\n</Relationships>\n')
    files["_rels/.rels"] = (
        '<?xml version="1.0" encoding="UTF-8"?>\n<Relationships '
        'xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
        ' <Relationship Target="/3D/3dmodel.model" Id="rel-1" '
        'Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>\n</Relationships>\n')
    plate_meta = re.search(r"<plate>(.*?)<model_instance>", tpl["model_settings"], re.S).group(1)
    plate_meta = "\n".join(l for l in plate_meta.splitlines()
                           if l.strip() and not re.search(r'key="(thumbnail_file|thumbnail_no_light_file|top_file|pick_file)"', l))
    files["Metadata/model_settings.config"] = (
        '<?xml version="1.0" encoding="UTF-8"?>\n<config>\n' + "\n".join(objs) + "\n  <plate>\n"
        + plate_meta + "\n" + "\n".join(insts) + "\n  </plate>\n  <assemble>\n" + "\n".join(asm)
        + "\n  </assemble>\n</config>\n")

    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, plate_name + ".3mf")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in tpl["copied"]:
            z.writestr(name, data)
        for name, text in files.items():
            z.writestr(name, text)
    return out, parts, places


def load_template(path=TEMPLATE):
    z = zipfile.ZipFile(path)
    model = z.read("3D/3dmodel.model").decode()
    ms = z.read("Metadata/model_settings.config").decode()
    # Bambu project metadata (application, version, dates...) minus the thumbnail references
    app_meta = [l for l in model.splitlines() if l.strip().startswith("<metadata")
                and "Thumbnail" not in l]
    extruder = re.search(r'<metadata key="extruder" value="(\d+)"', ms).group(1)
    copied = [(i.filename, z.read(i.filename)) for i in z.infolist() if not REBUILT.match(i.filename)]
    return dict(model_settings=ms, app_meta=app_meta, extruder=extruder, copied=copied)


def build_all(template=TEMPLATE):
    tpl = load_template(template)
    report = {}
    for pname, names in PLATES.items():
        out, parts, places = build_plate(tpl, pname, names, tpl["extruder"], tpl["app_meta"])
        report[pname] = [(p["name"], round(p["w"], 1), round(p["d"], 1), (round(x, 1), round(y, 1)))
                         for p, (x, y, _) in zip(parts, places)]
        print(f"{os.path.relpath(out, HERE)}: " + ", ".join(
            f"{n} {w}x{d} @ {c}" for n, w, d, c in report[pname]))
    return report


if __name__ == "__main__":
    build_all()
