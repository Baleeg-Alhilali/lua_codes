#!/usr/bin/env python3
"""Measure pressure roughness at inherited and newly refined VTK vertices."""

import argparse
import base64
import collections
import struct
import xml.etree.ElementTree as ET

import numpy as np


VTK_TYPES = {"Float32": "f", "Float64": "d", "Int32": "i", "Int8": "b"}


def read_array(node):
    raw = base64.b64decode("".join(node.itertext()).strip())
    size = struct.unpack("<I", raw[:4])[0]
    code = VTK_TYPES[node.attrib["type"]]
    item_size = struct.calcsize(code)
    values = struct.unpack("<" + code * (size // item_size), raw[4:4 + size])
    components = int(node.attrib.get("NumberOfComponents", "1"))
    result = np.asarray(values)
    return result.reshape((-1, components)) if components > 1 else result


def read_vtu(path, pressure_name):
    root = ET.parse(path).getroot()
    arrays = list(root.iter("DataArray"))
    points = read_array(next(a for a in arrays if "Name" not in a.attrib))
    named = {a.attrib["Name"]: read_array(a) for a in arrays if "Name" in a.attrib}
    connectivity = named["connectivity"].astype(np.int64)
    offsets = named["offsets"].astype(np.int64)
    types = named["types"].astype(np.int64)
    pressure = named[pressure_name]
    cells = []
    begin = 0
    for end, cell_type in zip(offsets, types):
        if cell_type in (10, 12, 13):
            cells.append((int(cell_type), tuple(connectivity[begin:end])))
        begin = end
    return points, pressure, cells


def boundary_graph(cells):
    faces = collections.Counter()
    for cell_type, vertices in cells:
        if cell_type == 10:  # tetrahedron
            a, b, c, d = vertices
            cell_faces = ((a, b, c), (a, b, d), (a, c, d), (b, c, d))
        elif cell_type == 12:  # hexahedron
            a, b, c, d, e, f, g, h = vertices
            cell_faces = (
                (a, b, c, d), (e, f, g, h), (a, b, f, e),
                (b, c, g, f), (c, d, h, g), (d, a, e, h),
            )
        else:  # VTK wedge
            a, b, c, d, e, f = vertices
            cell_faces = (
                (a, b, c), (d, e, f), (a, b, e, d),
                (b, c, f, e), (c, a, d, f),
            )
        for face in cell_faces:
            faces[tuple(sorted(face))] += 1
    boundary_faces = [face for face, count in faces.items() if count == 1]
    graph = collections.defaultdict(set)
    for face in boundary_faces:
        for i, vertex in enumerate(face):
            next_vertex = face[(i + 1) % len(face)]
            graph[vertex].add(next_vertex)
            graph[next_vertex].add(vertex)
    return graph


def coordinate_keys(points, tolerance):
    return {tuple(np.rint(point / tolerance).astype(np.int64)) for point in points}


def local_affine_residuals(points, pressure, graph):
    residual = {}
    for vertex, neighbors in graph.items():
        stencil = [vertex] + sorted(neighbors)
        xyz = points[stencil]
        design = np.column_stack((np.ones(len(stencil)), xyz))
        fit, *_ = np.linalg.lstsq(design, pressure[stencil], rcond=None)
        residual[vertex] = pressure[vertex] - design[0].dot(fit)
    return residual


def report(label, values):
    values = np.asarray(values, dtype=float)
    if not len(values):
        print(f"{label}: no vertices")
        return
    print(
        f"{label}: n={len(values)}, rms={np.sqrt(np.mean(values**2)):.6e}, "
        f"mean_abs={np.mean(np.abs(values)):.6e}, max_abs={np.max(np.abs(values)):.6e}"
    )


def ring_statistics(coarse_points, points, pressure):
    """Report circumferential pressure variation on repeated axial/radial rings."""
    scale = max(points[:, 2].max() - points[:, 2].min(), 1.0e-12)
    z_tolerance = scale * 1.0e-7
    radial_tolerance = max(np.ptp(points[:, 0]), np.ptp(points[:, 1])) * 2.0e-4
    coarse_z = {int(round(z / z_tolerance)) for z in coarse_points[:, 2]}
    z_groups = collections.defaultdict(list)
    for index, point in enumerate(points):
        z_groups[int(round(point[2] / z_tolerance))].append(index)

    inherited_ring_std, new_ring_std = [], []
    z_min, z_max = points[:, 2].min(), points[:, 2].max()
    middle_min = z_min + 0.2 * (z_max - z_min)
    middle_max = z_max - 0.2 * (z_max - z_min)
    for z_key, indices in z_groups.items():
        if len(indices) < 8:
            continue
        z_value = np.mean(points[indices, 2])
        if not middle_min <= z_value <= middle_max:
            continue
        center = points[indices, :2].mean(axis=0)
        radii = np.linalg.norm(points[indices, :2] - center, axis=1)
        shells = collections.defaultdict(list)
        for index, radius in zip(indices, radii):
            shells[int(round(radius / radial_tolerance))].append(index)
        for shell in shells.values():
            if len(shell) < 8:
                continue
            value = float(np.std(pressure[shell]))
            (inherited_ring_std if z_key in coarse_z else new_ring_std).append(value)

    report("inherited axial-ring circumferential std", inherited_ring_std)
    report("new axial-ring circumferential std", new_ring_std)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coarse_vtu")
    parser.add_argument("refined_vtu")
    parser.add_argument("pressure_name")
    parser.add_argument("--tolerance", type=float, default=2.0e-7)
    args = parser.parse_args()

    coarse_points, _, _ = read_vtu(args.coarse_vtu, args.pressure_name)
    points, pressure, cells = read_vtu(args.refined_vtu, args.pressure_name)
    inherited_keys = coordinate_keys(coarse_points, args.tolerance)
    graph = boundary_graph(cells)
    residuals = local_affine_residuals(points, pressure, graph)

    inherited, new = [], []
    inherited_middle, new_middle = [], []
    z_min, z_max = points[:, 2].min(), points[:, 2].max()
    middle_min = z_min + 0.2 * (z_max - z_min)
    middle_max = z_max - 0.2 * (z_max - z_min)
    for vertex, value in residuals.items():
        key = tuple(np.rint(points[vertex] / args.tolerance).astype(np.int64))
        is_inherited = key in inherited_keys
        (inherited if is_inherited else new).append(value)
        if middle_min <= points[vertex, 2] <= middle_max:
            (inherited_middle if is_inherited else new_middle).append(value)

    print(f"pressure range: [{pressure.min():.9e}, {pressure.max():.9e}]")
    print(f"boundary vertices: {len(graph)}")
    report("inherited boundary vertices", inherited)
    report("new boundary vertices", new)
    print(f"middle 60% z interval: [{middle_min:.9e}, {middle_max:.9e}]")
    report("inherited middle boundary vertices", inherited_middle)
    report("new middle boundary vertices", new_middle)
    ring_statistics(coarse_points, points, pressure)


if __name__ == "__main__":
    main()
