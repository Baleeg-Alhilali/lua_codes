import sys, math, argparse
sys.path.insert(0, "/opt/homebrew/opt/gmsh/lib")
import gmsh

parser = argparse.ArgumentParser(description="Generate a conforming three-domain nephron UGX mesh")
parser.add_argument("--swc", required=True)
parser.add_argument("--msh", required=True)
parser.add_argument("--ugx", required=True)
parser.add_argument("--lumen-radius", type=float, default=1.0)
parser.add_argument("--membrane-thickness", type=float, default=0.2)
parser.add_argument("--padding", type=float, default=10.0)
parser.add_argument("--near-size", type=float, default=0.5)
parser.add_argument("--far-size", type=float, default=12.0)
args = parser.parse_args()
swc = args.swc
points = []
for line in open(swc):
    fields = line.split("#", 1)[0].split()
    if len(fields) == 7:
        points.append(tuple(map(float, fields[2:5])))

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("nephron_v3")
occ = gmsh.model.occ

point_tags = [occ.addPoint(*p) for p in points]
curve = occ.addSpline(point_tags)
wire = occ.addWire([curve])

def swept_solid(radius):
    p0, p1 = points[0], points[1]
    tangent = [p1[i] - p0[i] for i in range(3)]
    length = math.sqrt(sum(x*x for x in tangent))
    tangent = [x / length for x in tangent]
    disk = occ.addDisk(*p0, radius, radius)
    z = (0.0, 0.0, 1.0)
    axis = (-tangent[1], tangent[0], 0.0)
    axis_len = math.sqrt(sum(x*x for x in axis))
    dot = max(-1.0, min(1.0, tangent[2]))
    angle = math.acos(dot)
    if axis_len > 1e-12:
        occ.rotate([(2, disk)], *p0, *(x / axis_len for x in axis), angle)
    elif dot < 0.0:
        occ.rotate([(2, disk)], *p0, 1, 0, 0, math.pi)
    result = occ.addPipe([(2, disk)], wire, "DiscreteTrihedron")
    volumes = [tag for dim, tag in result if dim == 3]
    if len(volumes) != 1:
        raise RuntimeError(f"Pipe produced {result}")
    return volumes[0]

outer_radius = args.lumen_radius + args.membrane_thickness
inner = swept_solid(args.lumen_radius)
outer = swept_solid(outer_radius)
margin = outer_radius + args.padding
mins = [min(p[i] for p in points) - margin for i in range(3)]
maxs = [max(p[i] for p in points) + margin for i in range(3)]
box = occ.addBox(*mins, *(maxs[i] - mins[i] for i in range(3)))
fragments, mapping = occ.fragment([(3, box)], [(3, outer), (3, inner)])
occ.synchronize()
print("FRAGMENTS", fragments)
print("MAPPING", mapping)

lumen, membrane, inter = 1, 3, 2
gmsh.model.addPhysicalGroup(3, [lumen], 1, "Lumen")
gmsh.model.addPhysicalGroup(3, [membrane], 2, "Membrane")
gmsh.model.addPhysicalGroup(3, [inter], 3, "Inter")

def boundary(volume):
    return {tag for dim, tag in gmsh.model.getBoundary([(3, volume)], False, False) if dim == 2}

b_lumen, b_membrane, b_inter = boundary(lumen), boundary(membrane), boundary(inter)
apical = sorted(b_lumen & b_membrane)
basolateral = sorted(b_membrane & b_inter)
outer_candidates = sorted(b_inter - b_membrane)
tol = 1e-6
outer_wall, lumen_ends = [], []
for surface in outer_candidates:
    bb = gmsh.model.getBoundingBox(2, surface)
    on_box = any(abs(bb[d] - mins[d]) < tol and abs(bb[d + 3] - mins[d]) < tol
                 or abs(bb[d] - maxs[d]) < tol and abs(bb[d + 3] - maxs[d]) < tol
                 for d in range(3))
    (outer_wall if on_box else lumen_ends).append(surface)
def center_of_mass(surface):
    return gmsh.model.occ.getCenterOfMass(2, surface)

def squared_distance(a, b):
    return sum((a[d] - b[d]) ** 2 for d in range(3))

inlet = [min(lumen_ends, key=lambda surface: squared_distance(center_of_mass(surface), points[0]))]
outlet = [surface for surface in lumen_ends if surface not in inlet]
print("SURFACES", "Apical", apical, "Basolateral", basolateral,
      "OuterWall", outer_wall, "Inlet", inlet, "Outlet", outlet)
if not apical or not basolateral or len(outer_wall) != 6 or len(lumen_ends) != 2:
    raise RuntimeError("Unexpected interface topology")
gmsh.model.addPhysicalGroup(2, apical, 4, "Apical")
gmsh.model.addPhysicalGroup(2, basolateral, 5, "Basolateral")
gmsh.model.addPhysicalGroup(2, outer_wall, 6, "OuterWall")
gmsh.model.addPhysicalGroup(2, inlet, 7, "Inlet")
gmsh.model.addPhysicalGroup(2, outlet, 8, "Outlet")

distance = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance, "SurfacesList", basolateral)
threshold = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
gmsh.model.mesh.field.setNumber(threshold, "SizeMin", args.near_size)
gmsh.model.mesh.field.setNumber(threshold, "SizeMax", args.far_size)
gmsh.model.mesh.field.setNumber(threshold, "DistMin", 1.0)
gmsh.model.mesh.field.setNumber(threshold, "DistMax", 15.0)
gmsh.model.mesh.field.setAsBackgroundMesh(threshold)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Binary", 0)
gmsh.model.mesh.generate(3)
gmsh.write(args.msh)
gmsh.finalize()

# Convert Gmsh MSH 2.2 to UGX and reproduce ProMesh's ordered subset closure:
# close/assign the three volumes first, then close/assign boundary surfaces.
physical, nodes = {}, {}
triangles, tetrahedra = [], []
tri_subsets, tet_subsets = {}, {}
with open(args.msh) as source:
    iterator = iter(source)
    for line in iterator:
        key = line.strip()
        if key == "$PhysicalNames":
            for _ in range(int(next(iterator))):
                dim, tag, name = next(iterator).split(maxsplit=2)
                physical[(int(dim), int(tag))] = name.strip().strip('"')
            assert next(iterator).strip() == "$EndPhysicalNames"
        elif key == "$Nodes":
            for _ in range(int(next(iterator))):
                fields = next(iterator).split()
                nodes[int(fields[0])] = tuple(map(float, fields[1:4]))
            assert next(iterator).strip() == "$EndNodes"
        elif key == "$Elements":
            for _ in range(int(next(iterator))):
                fields = list(map(int, next(iterator).split()))
                element_type, num_tags = fields[1], fields[2]
                physical_tag = fields[3] if num_tags else 0
                connectivity = fields[3 + num_tags:]
                if element_type == 2:
                    triangles.append(connectivity)
                    name = physical.get((2, physical_tag), "surface_%d" % physical_tag)
                    tri_subsets.setdefault(name, []).append(connectivity)
                elif element_type == 4:
                    index = len(tetrahedra)
                    tetrahedra.append(connectivity)
                    name = physical.get((3, physical_tag), "volume_%d" % physical_tag)
                    tet_subsets.setdefault(name, []).append(index)
            assert next(iterator).strip() == "$EndElements"

used_nodes = {node for element in triangles + tetrahedra for node in element}
node_ids = sorted(used_nodes)
node_index = {node_id: index for index, node_id in enumerate(node_ids)}
subset_names = ["Lumen", "Membrane", "Inter", "Apical", "Basolateral", "OuterWall", "Inlet", "Outlet"]
colors = ["0.2 0.5 1 1", "1 0.4 0.2 1", "0.7 0.7 0.7 1",
          "1 0 1 1", "0 0.8 0.2 1", "0.2 0.2 0.2 1", "0 1 1 1", "1 1 0 1"]

# Build every tetrahedron edge and face once. These are the elements that
# ProMesh's CloseSelection would add to a selected volume.
edge_keys, face_keys = [], []
edge_index, face_index = {}, {}
tet_edge_indices, tet_face_indices = [], []
for tet in tetrahedra:
    local_edges, local_faces = [], []
    for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        key = tuple(sorted((tet[a], tet[b])))
        if key not in edge_index:
            edge_index[key] = len(edge_keys)
            edge_keys.append(key)
        local_edges.append(edge_index[key])
    for a, b, c in ((0, 1, 2), (0, 3, 1), (1, 3, 2), (2, 3, 0)):
        key = tuple(sorted((tet[a], tet[b], tet[c])))
        if key not in face_index:
            face_index[key] = len(face_keys)
            face_keys.append((tet[a], tet[b], tet[c]))
        local_faces.append(face_index[key])
    tet_edge_indices.append(local_edges)
    tet_face_indices.append(local_faces)

# Assignment order matters: volumes first, named boundary faces afterward.
vertex_owner = [None] * len(node_ids)
edge_owner = [None] * len(edge_keys)
face_owner = [None] * len(face_keys)
volume_owner = [None] * len(tetrahedra)
for subset_id, name in enumerate(subset_names[:3]):
    for tet_id in tet_subsets.get(name, []):
        volume_owner[tet_id] = subset_id
        for node in tetrahedra[tet_id]:
            vertex_owner[node_index[node]] = subset_id
        for element_id in tet_edge_indices[tet_id]:
            edge_owner[element_id] = subset_id
        for element_id in tet_face_indices[tet_id]:
            face_owner[element_id] = subset_id

for subset_id, name in enumerate(subset_names[3:], start=3):
    for triangle in tri_subsets.get(name, []):
        face_owner[face_index[tuple(sorted(triangle))]] = subset_id
        for node in triangle:
            vertex_owner[node_index[node]] = subset_id
        for a, b in ((0, 1), (1, 2), (2, 0)):
            edge_owner[edge_index[tuple(sorted((triangle[a], triangle[b])))]] = subset_id

with open(args.ugx, "w") as output:
    output.write('<?xml version="1.0" encoding="utf-8"?>\n<grid name="defGrid">\n')
    output.write('\t<vertices coords="3">' + ' '.join(
        '%.17g' % value for node_id in node_ids for value in nodes[node_id]) + '</vertices>\n')
    output.write('\t<edges>' + ' '.join(
        str(node_index[node]) for element in edge_keys for node in element) + '</edges>\n')
    output.write('\t<triangles>' + ' '.join(
        str(node_index[node]) for element in face_keys for node in element) + '</triangles>\n')
    output.write('\t<tetrahedrons>' + ' '.join(
        str(node_index[node]) for element in tetrahedra for node in element) + '</tetrahedrons>\n')
    output.write('\t<subset_handler name="defSH">\n')
    for name, color in zip(subset_names, colors):
        subset_id = subset_names.index(name)
        output.write('\t\t<subset name="%s" color="%s" state="0">\n' % (name, color))
        for tag, owners in (("vertices", vertex_owner), ("edges", edge_owner),
                            ("faces", face_owner), ("volumes", volume_owner)):
            indices = [str(index) for index, owner in enumerate(owners) if owner == subset_id]
            if indices:
                output.write('\t\t\t<%s>%s</%s>\n' % (tag, ' '.join(indices), tag))
        output.write('\t\t</subset>\n')
    output.write('\t</subset_handler>\n</grid>\n')

print("UGX written:", args.ugx)
print("Mesh counts:", len(node_ids), "vertices,", len(edge_keys), "edges,",
      len(face_keys), "triangles,", len(tetrahedra), "tetrahedra")
print("Subset order: volumes Lumen/Membrane/Inter, then boundaries",
      "Apical/Basolateral/OuterWall/Inlet/Outlet")
