#!/usr/bin/env python3
"""Validate UGX volume orientation and face/subset adjacency.

Usage: python3 mesh_generation/validate_ugx_mesh.py mesh.ugx
Returns a non-zero status if an inverted/degenerate tetrahedron or a topology
error is found.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import math
import sys
import xml.etree.ElementTree as ET


FACE_LAYOUTS = {
    "tetrahedrons": ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)),
    "hexahedrons": (
        (0, 1, 2, 3), (4, 5, 6, 7),
        (0, 1, 5, 4), (1, 2, 6, 5),
        (2, 3, 7, 6), (3, 0, 4, 7),
    ),
    "prisms": (
        (0, 1, 2), (3, 4, 5),
        (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5),
    ),
    "pyramids": (
        (0, 1, 2, 3),
        (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4),
    ),
    "octahedrons": (
        (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4),
        (1, 0, 5), (2, 1, 5), (3, 2, 5), (0, 3, 5),
    ),
}

ELEMENT_ARITY = {
    "tetrahedrons": 4,
    "hexahedrons": 8,
    "prisms": 6,
    "pyramids": 5,
    "octahedrons": 6,
}

# This is also the order in which GridWriterUGX assigns face IDs.  The order
# matters because subset_handler face indices refer to this concatenated list.
FACE_TAGS = (
    ("triangles", 3, False),
    ("constraining_triangles", 3, False),
    ("constrained_triangles", 3, True),
    ("quadrilaterals", 4, False),
    ("constraining_quadrilaterals", 4, False),
    ("constrained_quadrilaterals", 4, True),
)


def records(root: ET.Element, tag: str, stride: int, keep: int | None = None):
    elem = root.find(tag)
    values = [] if elem is None else [int(x) for x in (elem.text or "").split()]
    if len(values) % stride:
        raise ValueError(f"<{tag}> has {len(values)} integers, not divisible by {stride}")
    width = stride if keep is None else keep
    return [tuple(values[i:i + width]) for i in range(0, len(values), stride)]


def constrained_face_records(root: ET.Element, tag: str, arity: int):
    """Read UGX constrained faces, whose records have variable length.

    A record contains ``arity`` vertex IDs, a constraining-object type, and an
    object ID unless the type is -1.  Treating these as fixed-width records
    shifts every later face ID and creates false topology/interface failures.
    """
    elem = root.find(tag)
    values = [] if elem is None else [int(x) for x in (elem.text or "").split()]
    result = []
    pos = 0
    while pos < len(values):
        if pos + arity >= len(values):
            raise ValueError(f"<{tag}> ends inside a constrained-face record")
        conn = tuple(values[pos:pos + arity])
        pos += arity
        object_type = values[pos]
        pos += 1
        if object_type != -1:
            if pos >= len(values):
                raise ValueError(f"<{tag}> is missing a constraining-object ID")
            pos += 1
        result.append(conn)
    return result


def subset_maps(root: ET.Element):
    handler = next(
        (h for h in root.findall("subset_handler") if h.get("name") == "defSH"),
        None,
    )
    if handler is None:
        raise ValueError("No subset_handler named defSH")
    by_dim = {"vertices": {}, "edges": {}, "faces": {}, "volumes": {}}
    counts = {}
    duplicates = []
    for subset in handler.findall("subset"):
        name = subset.get("name", "<unnamed>")
        counts[name] = {}
        for dim in by_dim:
            child = subset.find(dim)
            ids = [] if child is None else [int(x) for x in (child.text or "").split()]
            counts[name][dim] = len(ids)
            for idx in ids:
                if idx in by_dim[dim]:
                    duplicates.append((dim, idx, by_dim[dim][idx], name))
                by_dim[dim][idx] = name
    return by_dim, counts, duplicates


def signed_tet_volume(vertices, tet):
    a, b, c, d = (vertices[i] for i in tet)
    ba = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    ca = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    da = (d[0] - a[0], d[1] - a[1], d[2] - a[2])
    cross = (
        ca[1] * da[2] - ca[2] * da[1],
        ca[2] * da[0] - ca[0] * da[2],
        ca[0] * da[1] - ca[1] * da[0],
    )
    return (ba[0] * cross[0] + ba[1] * cross[1] + ba[2] * cross[2]) / 6.0


def center(vertices, conn):
    n = len(conn)
    return tuple(sum(vertices[i][axis] for i in conn) / n for axis in range(3))


def edge_ratio(vertices, conn):
    lengths = []
    for i in range(len(conn)):
        for j in range(i + 1, len(conn)):
            a, b = vertices[conn[i]], vertices[conn[j]]
            lengths.append(math.dist(a, b))
    return min(lengths) / max(lengths)


def normalized_tet_volume(vertices, tet, signed_volume):
    """Return |6 V| / L_max^3, a dimensionless tetrahedron flatness metric."""
    lengths = []
    for i in range(4):
        for j in range(i + 1, 4):
            lengths.append(math.dist(vertices[tet[i]], vertices[tet[j]]))
    longest = max(lengths)
    if longest == 0.0:
        return 0.0
    return abs(6.0 * signed_volume) / (longest ** 3)


def prism_signed_volume(vertices, prism):
    """Signed volume of a UG4 wedge as three consistently oriented tets."""
    a, b, c, d, e, f = prism
    return sum(
        signed_tet_volume(vertices, tet)
        for tet in ((a, b, c, d), (b, c, e, d), (c, e, f, d))
    )


def scaled_corner_jacobian(vertices, origin, first, second, third):
    edges = []
    for target in (first, second, third):
        edges.append(tuple(
            vertices[target][axis] - vertices[origin][axis]
            for axis in range(3)
        ))
    lengths = [math.sqrt(sum(value * value for value in edge)) for edge in edges]
    if min(lengths) == 0.0:
        return 0.0
    cross = (
        edges[1][1] * edges[2][2] - edges[1][2] * edges[2][1],
        edges[1][2] * edges[2][0] - edges[1][0] * edges[2][2],
        edges[1][0] * edges[2][1] - edges[1][1] * edges[2][0],
    )
    determinant = sum(edges[0][axis] * cross[axis] for axis in range(3))
    return determinant / math.prod(lengths)


def prism_min_scaled_jacobian(vertices, prism):
    a, b, c, d, e, f = prism
    corners = (
        (a, b, c, d), (b, c, a, e), (c, a, b, f),
        (d, f, e, a), (e, d, f, b), (f, e, d, c),
    )
    return min(scaled_corner_jacobian(vertices, *corner) for corner in corners)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh")
    parser.add_argument("--details", action="store_true")
    parser.add_argument(
        "--interface-contract", choices=("production", "none"),
        default="production",
        help=("validate the nephron production subset contract (default), or "
              "only geometry/topology for diagnostic refinement fixtures"),
    )
    parser.add_argument(
        "--min-edge-ratio", type=float, default=1.0e-4,
        help="fail below this dimensionless shortest/longest vertex distance",
    )
    parser.add_argument(
        "--min-normalized-tet-volume", type=float, default=1.0e-10,
        help="fail below |6V|/Lmax^3 (dimensionless tetrahedron flatness)",
    )
    parser.add_argument(
        "--min-prism-scaled-jacobian", type=float, default=1.0e-8,
        help="fail below the minimum signed corner Jacobian of any prism",
    )
    parser.add_argument(
        "--straight-lumen", nargs=3, type=float,
        metavar=("CENTER_X", "CENTER_Y", "RADIUS"),
        help=("also require every Apical vertex to lie on the specified "
              "z-aligned cylindrical lumen"),
    )
    parser.add_argument(
        "--lumen-radius-tolerance", type=float, default=1.0e-4,
        help="absolute radius tolerance used with --straight-lumen",
    )
    args = parser.parse_args()

    root = ET.parse(args.mesh).getroot()
    if root.tag != "grid":
        raise ValueError(f"Expected <grid> root, found <{root.tag}>")

    coords = [float(x) for x in (root.findtext("vertices") or "").split()]
    if len(coords) % 3:
        raise ValueError("Vertex coordinate count is not divisible by 3")
    vertices = [tuple(coords[i:i + 3]) for i in range(0, len(coords), 3)]

    volumes = []
    volume_kinds = []
    volume_counts = {}
    for kind, arity in ELEMENT_ARITY.items():
        elems = records(root, kind, arity)
        volume_counts[kind] = len(elems)
        volumes.extend(elems)
        volume_kinds.extend([kind] * len(elems))

    serialized_faces = []
    face_counts = {}
    for tag, arity, constrained in FACE_TAGS:
        faces = (constrained_face_records(root, tag, arity) if constrained
                 else records(root, tag, arity))
        face_counts[tag] = len(faces)
        serialized_faces.extend(faces)

    subsets, subset_counts, duplicate_assignments = subset_maps(root)
    volume_subset = subsets["volumes"]
    face_subset = subsets["faces"]

    invalid_indices = []
    adjacency = defaultdict(list)
    min_ratio = (float("inf"), -1, "")
    for volume_id, (kind, conn) in enumerate(zip(volume_kinds, volumes)):
        if any(i < 0 or i >= len(vertices) for i in conn):
            invalid_indices.append((volume_id, kind, conn))
            continue
        ratio = edge_ratio(vertices, conn)
        if ratio < min_ratio[0]:
            min_ratio = (ratio, volume_id, kind)
        for layout in FACE_LAYOUTS[kind]:
            key = tuple(sorted(conn[i] for i in layout))
            adjacency[key].append(volume_id)

    missing_face_geometry = []
    adjacency_hist = Counter()
    subset_pairs = defaultdict(Counter)
    face_keys_seen = Counter()
    for face_id, conn in enumerate(serialized_faces):
        key = tuple(sorted(conn))
        face_keys_seen[key] += 1
        adjacent = adjacency.get(key, [])
        adjacency_hist[len(adjacent)] += 1
        if not adjacent:
            missing_face_geometry.append((face_id, conn, face_subset.get(face_id)))
        name = face_subset.get(face_id, "<unassigned>")
        pair = tuple(sorted(volume_subset.get(v, "<unassigned>") for v in adjacent))
        subset_pairs[name][pair] += 1

    duplicate_face_geometry = [key for key, count in face_keys_seen.items() if count != 1]
    nonmanifold = [(key, ids) for key, ids in adjacency.items() if len(ids) > 2]

    expected_pairs = {
        "Apical": {("Lumen", "Membrane")},
        "Inlet": {("Lumen",)},
        "Outlet": {("Lumen",)},
        "OuterWall": {("Inter",)},
        "MembraneWall": {("Membrane",)},
    }
    interface_errors = {}
    if args.interface_contract == "production":
        for name, expected in expected_pairs.items():
            unexpected = {pair: count for pair, count in subset_pairs[name].items() if pair not in expected}
            if unexpected:
                interface_errors[name] = unexpected
        basolateral_observed = set(subset_pairs["Basolateral"])
        historical_basolateral = {("Inter",), ("Membrane",)}
        conforming_basolateral = {("Inter", "Membrane")}
        if not basolateral_observed or not (
            basolateral_observed <= historical_basolateral
            or basolateral_observed == conforming_basolateral
        ):
            interface_errors["Basolateral"] = dict(
                subset_pairs["Basolateral"]
            )

    tets = records(root, "tetrahedrons", 4)
    negative = []
    zero = []
    min_abs = (float("inf"), -1, 0.0)
    min_normalized = (float("inf"), -1)
    for tet_id, tet in enumerate(tets):
        value = signed_tet_volume(vertices, tet)
        if abs(value) < min_abs[0]:
            min_abs = (abs(value), tet_id, value)
        if value < 0.0:
            negative.append((tet_id, value, volume_subset.get(tet_id), center(vertices, tet)))
        elif value == 0.0:
            zero.append(tet_id)
        normalized = normalized_tet_volume(vertices, tet, value)
        if normalized < min_normalized[0]:
            min_normalized = (normalized, tet_id)

    prisms = records(root, "prisms", 6)
    negative_prisms = []
    zero_prisms = []
    min_prism_volume = (float("inf"), -1, 0.0)
    min_prism_jacobian = (float("inf"), -1)
    prism_volume_offset = volume_counts["tetrahedrons"]
    for prism_id, prism in enumerate(prisms):
        value = prism_signed_volume(vertices, prism)
        if abs(value) < min_prism_volume[0]:
            min_prism_volume = (abs(value), prism_id, value)
        if value < 0.0:
            negative_prisms.append((
                prism_id, value,
                volume_subset.get(prism_volume_offset + prism_id),
                center(vertices, prism),
            ))
        elif value == 0.0:
            zero_prisms.append(prism_id)
        jacobian = prism_min_scaled_jacobian(vertices, prism)
        if jacobian < min_prism_jacobian[0]:
            min_prism_jacobian = (jacobian, prism_id)

    print(f"mesh: {args.mesh}")
    print(f"vertices: {len(vertices)}")
    print("volumes:", ", ".join(f"{k}={v}" for k, v in volume_counts.items()))
    print("faces:", ", ".join(f"{k}={v}" for k, v in face_counts.items()))
    print(f"minimum all-pairs edge ratio: {min_ratio[0]:.8g} (volume {min_ratio[1]}, {min_ratio[2]})")
    print(f"tetrahedra: inverted={len(negative)}, exactly_zero={len(zero)}, min_abs_signed_volume={min_abs[0]:.12g}")
    print(f"minimum normalized tetrahedron volume: {min_normalized[0]:.12g} (tetrahedron {min_normalized[1]})")
    print(
        "prisms: "
        f"inverted={len(negative_prisms)}, exactly_zero={len(zero_prisms)}, "
        f"min_abs_signed_volume={min_prism_volume[0]:.12g}, "
        f"min_scaled_jacobian={min_prism_jacobian[0]:.12g} "
        f"(prism {min_prism_jacobian[1]})"
    )
    print("serialized face adjacency:", dict(sorted(adjacency_hist.items())))
    print(f"unmatched serialized faces: {len(missing_face_geometry)}")
    print(f"duplicate serialized face geometries: {len(duplicate_face_geometry)}")
    print(f"non-manifold volume faces (>2 neighbors): {len(nonmanifold)}")
    print(f"duplicate subset assignments: {len(duplicate_assignments)}")
    print(f"interface contract: {args.interface_contract}")
    print("interface adjacency:")
    for name in ("Apical", "Basolateral", "Inlet", "Outlet", "OuterWall", "MembraneWall"):
        pairs = ", ".join(f"{pair}={count}" for pair, count in subset_pairs[name].items())
        print(f"  {name}: {pairs or '<none>'}")

    lumen_geometry_error = False
    if args.straight_lumen:
        center_x, center_y, expected_radius = args.straight_lumen
        apical_vertex_ids = {
            vertex_id
            for face_id, conn in enumerate(serialized_faces)
            if face_subset.get(face_id) == "Apical"
            for vertex_id in conn
        }
        radial_errors = [
            abs(math.hypot(vertices[i][0] - center_x,
                           vertices[i][1] - center_y) - expected_radius)
            for i in apical_vertex_ids
        ]
        max_radial_error = max(radial_errors, default=float("inf"))
        lumen_geometry_error = (
            not apical_vertex_ids or
            max_radial_error > args.lumen_radius_tolerance
        )
        print(
            "straight-lumen conformance: "
            f"Apical_vertices={len(apical_vertex_ids)}, "
            f"expected_radius={expected_radius:.12g}, "
            f"max_abs_radius_error={max_radial_error:.12g}, "
            f"tolerance={args.lumen_radius_tolerance:.12g}"
        )

    if args.details and negative:
        print("inverted tetrahedra (global_volume_id, signed_volume, subset, center):")
        for item in negative:
            print(f"  {item[0]} {item[1]:.12e} {item[2]} {item[3]}")

    failures = (
        len(invalid_indices) + len(negative) + len(zero) +
        len(negative_prisms) + len(zero_prisms) +
        len(missing_face_geometry) + len(duplicate_face_geometry) +
        len(nonmanifold) + len(duplicate_assignments) + len(interface_errors) +
        int(lumen_geometry_error)
    )
    if min_ratio[0] < args.min_edge_ratio:
        failures += 1
        print(f"edge-ratio quality gate failed: {min_ratio[0]:.12g} < {args.min_edge_ratio:.12g}")
    if min_normalized[0] < args.min_normalized_tet_volume:
        failures += 1
        print("normalized-tetrahedron-volume quality gate failed: "
              f"{min_normalized[0]:.12g} < {args.min_normalized_tet_volume:.12g}")
    if min_prism_jacobian[0] < args.min_prism_scaled_jacobian:
        failures += 1
        print("prism scaled-Jacobian quality gate failed: "
              f"{min_prism_jacobian[0]:.12g} < "
              f"{args.min_prism_scaled_jacobian:.12g}")
    if failures:
        print("VALIDATION: FAIL")
        if interface_errors:
            print("unexpected interface pairs:", interface_errors)
        return 1
    print("VALIDATION: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
