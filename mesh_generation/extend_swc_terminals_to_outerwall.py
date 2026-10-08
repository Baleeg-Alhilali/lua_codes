#!/usr/bin/env python3
"""Extend unbranched SWC terminal paths to a fixed padded axis-aligned box.

The input is expected to use the mesh coordinate unit already (millimetres in
the nephron pipeline).  The output carries UG4 metadata comments containing the
fixed box bounds.  Each terminal transition is a low-bending cubic Bezier curve
followed by a straight wall-normal tail.  The final SWC segments are therefore
normal to the selected box face.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable, Sequence


Vec = tuple[float, float, float]


@dataclass
class Node:
    old_id: int
    kind: int
    p: Vec
    radius: float
    parent: int


def add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def mul(a: Vec, value: float) -> Vec:
    return (a[0] * value, a[1] * value, a[2] * value)


def dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(a: Vec) -> float:
    return math.sqrt(dot(a, a))


def normalize(a: Vec) -> Vec:
    length = norm(a)
    if length <= 1e-14:
        raise ValueError("Cannot normalize a zero-length terminal direction")
    return mul(a, 1.0 / length)


def distance(a: Vec, b: Vec) -> float:
    return norm(sub(a, b))


def cross(a: Vec, b: Vec) -> Vec:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def bezier(p0: Vec, c1: Vec, c2: Vec, p3: Vec, t: float) -> Vec:
    s = 1.0 - t
    return (
        s**3 * p0[0] + 3*s*s*t*c1[0] + 3*s*t*t*c2[0] + t**3*p3[0],
        s**3 * p0[1] + 3*s*s*t*c1[1] + 3*s*t*t*c2[1] + t**3*p3[1],
        s**3 * p0[2] + 3*s*s*t*c1[2] + 3*s*t*t*c2[2] + t**3*p3[2],
    )


def bezier_derivatives(p0: Vec, c1: Vec, c2: Vec, p3: Vec, t: float) -> tuple[Vec, Vec]:
    s = 1.0 - t
    first = add(
        add(mul(sub(c1, p0), 3*s*s), mul(sub(c2, c1), 6*s*t)),
        mul(sub(p3, c2), 3*t*t),
    )
    second = add(
        mul(add(sub(c2, mul(c1, 2.0)), p0), 6*s),
        mul(add(sub(p3, mul(c2, 2.0)), c1), 6*t),
    )
    return first, second


def read_swc(path: Path) -> tuple[list[Node], list[str]]:
    nodes: list[Node] = []
    comments: list[str] = []
    for raw in path.read_text().splitlines():
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if not stripped.startswith("# UG4_OUTER_BOX") and not stripped.startswith("# UG4_TERMINALS_ON_OUTER_WALL"):
                comments.append(raw)
            continue
        payload = raw.split("#", 1)[0].split()
        if len(payload) != 7:
            raise ValueError(f"SWC row must contain seven columns: {raw}")
        nodes.append(Node(
            int(payload[0]), int(payload[1]),
            (float(payload[2]), float(payload[3]), float(payload[4])),
            float(payload[5]), int(payload[6]),
        ))
    if not nodes:
        raise ValueError("SWC contains no nodes")
    return nodes, comments


def ordered_chains(nodes: Sequence[Node]) -> list[list[Node]]:
    by_id = {node.old_id: node for node in nodes}
    if len(by_id) != len(nodes):
        raise ValueError("SWC node IDs must be unique")
    children: dict[int, list[int]] = {node.old_id: [] for node in nodes}
    roots: list[Node] = []
    for node in nodes:
        if node.parent == -1:
            roots.append(node)
        else:
            if node.parent not in by_id:
                raise ValueError(f"Missing parent {node.parent} for node {node.old_id}")
            children[node.parent].append(node.old_id)
    chains: list[list[Node]] = []
    visited: set[int] = set()
    for root in roots:
        chain = [root]
        visited.add(root.old_id)
        current = root
        while children[current.old_id]:
            if len(children[current.old_id]) != 1:
                raise ValueError("Terminal-to-wall extension currently requires unbranched SWC paths")
            current = by_id[children[current.old_id][0]]
            if current.old_id in visited:
                raise ValueError("SWC contains a cycle")
            visited.add(current.old_id)
            chain.append(current)
        if len(chain) < 2:
            raise ValueError("Each SWC path needs at least two nodes")
        chains.append(chain)
    if len(visited) != len(nodes):
        raise ValueError("SWC contains nodes not reachable from a root")
    return chains


def fixed_box(nodes: Sequence[Node], padding: Vec) -> tuple[Vec, Vec]:
    minimum = [math.inf, math.inf, math.inf]
    maximum = [-math.inf, -math.inf, -math.inf]
    for node in nodes:
        if node.radius <= 0.0:
            raise ValueError(f"Node {node.old_id} has nonpositive radius")
        for axis in range(3):
            minimum[axis] = min(minimum[axis], node.p[axis] - node.radius)
            maximum[axis] = max(maximum[axis], node.p[axis] + node.radius)
    return (
        tuple(minimum[axis] - padding[axis] for axis in range(3)),
        tuple(maximum[axis] + padding[axis] for axis in range(3)),
    )


def choose_face(p: Vec, tangent: Vec, box_min: Vec, box_max: Vec) -> tuple[int, int, Vec, float, float]:
    candidates = []
    for axis in range(3):
        for sign in (-1, 1):
            normal = tuple(float(sign if i == axis else 0) for i in range(3))
            plane = box_max[axis] if sign > 0 else box_min[axis]
            wall_distance = (plane - p[axis]) * sign
            alignment = dot(tangent, normal)
            if wall_distance > 1e-12 and alignment > 0.05:
                candidates.append((wall_distance, -alignment, axis, sign, normal, alignment))
    if not candidates:
        raise ValueError(
            "No outward box face is compatible with the terminal tangent; "
            "the trace endpoint may point back into the nephron"
        )
    min_distance = min(item[0] for item in candidates)
    near = [item for item in candidates if item[0] <= min_distance * (1.0 + 1e-8) + 1e-12]
    chosen = min(near, key=lambda item: item[1])
    wall_distance, _, axis, sign, normal, alignment = chosen
    return axis, sign, normal, wall_distance, alignment


def bending_objective(p0: Vec, c1: Vec, c2: Vec, p3: Vec, wall_normal: Vec, outward: Vec) -> float:
    energy = 0.0
    previous = p0
    for index in range(1, 81):
        t = index / 80.0
        point = bezier(p0, c1, c2, p3, t)
        first, second = bezier_derivatives(p0, c1, c2, p3, t)
        speed = norm(first)
        if speed <= 1e-12:
            return math.inf
        curvature = norm(cross(first, second)) / (speed**3)
        ds = distance(previous, point)
        energy += curvature * curvature * ds
        if dot(first, wall_normal) <= 1e-8:
            energy += 1e8
        if dot(sub(point, p0), outward) < -1e-9:
            energy += 1e8
        previous = point
    return energy + 1e-4 * distance(p0, p3)


def wall_cell_center(p: Vec, axis: int, sign: int, box_min: Vec, box_max: Vec,
                     box_refs: int) -> tuple[Vec, tuple[int, int]]:
    """Return the center of the structured OuterWall cell containing p's projection."""
    cells = 1 << box_refs
    target = list(p)
    target[axis] = box_max[axis] if sign > 0 else box_min[axis]
    tangential_axes = [dimension for dimension in range(3) if dimension != axis]
    indices = []
    for dimension in tangential_axes:
        extent = box_max[dimension] - box_min[dimension]
        scaled = (p[dimension] - box_min[dimension]) * cells / extent
        cell = max(0, min(cells - 1, int(math.floor(scaled))))
        target[dimension] = box_min[dimension] + (cell + 0.5) * extent / cells
        indices.append(cell)
    return (target[0], target[1], target[2]), (indices[0], indices[1])


def optimized_transition(p: Vec, tangent: Vec, normal: Vec, wall: Vec,
                         wall_distance: float, spacing: float,
                         radius: float) -> tuple[list[Vec], float]:
    tail = min(0.35 * wall_distance, max(3.0 * spacing, 2.0 * radius))
    tail = max(min(tail, 0.6 * wall_distance), 0.1 * wall_distance)
    transition_end = sub(wall, mul(normal, tail))
    chord = max(distance(p, transition_end), spacing)

    best = None
    fractions = [0.12 + 0.06 * i for i in range(19)]
    for f0 in fractions:
        for f1 in fractions:
            c1 = add(p, mul(tangent, f0 * chord))
            c2 = sub(transition_end, mul(normal, f1 * chord))
            objective = bending_objective(p, c1, c2, transition_end, normal, tangent)
            if best is None or objective < best[0]:
                best = (objective, c1, c2)
    if best is None or not math.isfinite(best[0]):
        raise ValueError("Could not construct a monotone low-curvature terminal transition")
    _, c1, c2 = best

    dense = [bezier(p, c1, c2, transition_end, i / 400.0) for i in range(401)]
    tail_steps = max(3, int(math.ceil(tail / max(spacing * 0.25, 1e-12))))
    dense.extend(add(transition_end, mul(normal, tail * i / tail_steps)) for i in range(1, tail_steps + 1))

    sampled = [dense[0]]
    carry = 0.0
    for a, b in zip(dense[:-1], dense[1:]):
        segment = distance(a, b)
        if segment <= 1e-14:
            continue
        while carry + segment >= spacing:
            fraction = (spacing - carry) / segment
            a = add(a, mul(sub(b, a), fraction))
            sampled.append(a)
            segment = distance(a, b)
            carry = 0.0
        carry += segment
    if distance(sampled[-1], wall) > 1e-10:
        sampled.append(wall)

    max_curvature = 0.0
    for i in range(1, 400):
        first, second = bezier_derivatives(p, c1, c2, transition_end, i / 400.0)
        speed = norm(first)
        if speed > 1e-12:
            max_curvature = max(max_curvature, norm(cross(first, second)) / speed**3)
    return sampled, max_curvature


def extend_chain(chain: Sequence[Node], box_min: Vec, box_max: Vec, spacing: float,
                 box_refs: int, label: str) -> tuple[list[tuple[Vec, float, int]], list[str]]:
    root_tangent = normalize(sub(chain[0].p, chain[1].p))
    tip_tangent = normalize(sub(chain[-1].p, chain[-2].p))
    diagnostics = []
    extensions = []
    for end_name, endpoint, tangent in (
        ("root", chain[0], root_tangent),
        ("tip", chain[-1], tip_tangent),
    ):
        axis, sign, normal, wall_distance, alignment = choose_face(endpoint.p, tangent, box_min, box_max)
        wall, cell = wall_cell_center(endpoint.p, axis, sign, box_min, box_max, box_refs)
        samples, max_curvature = optimized_transition(
            endpoint.p, tangent, normal, wall, wall_distance, spacing, endpoint.radius
        )
        face_name = "xyz"[axis] + ("max" if sign > 0 else "min")
        diagnostics.append(
            f"{label} {end_name}: {face_name}, distance={wall_distance:.9g}, "
            f"alignment={alignment:.6f}, cell={cell}, target={wall}, "
            f"added={len(samples)-1}, maxCurvature={max_curvature:.9g}"
        )
        extensions.append(samples)

    root_samples, tip_samples = extensions
    output: list[tuple[Vec, float, int]] = []
    for point in reversed(root_samples[1:]):
        output.append((point, chain[0].radius, chain[0].kind))
    output.extend((node.p, node.radius, node.kind) for node in chain)
    output.extend((point, chain[-1].radius, chain[-1].kind) for point in tip_samples[1:])
    return output, diagnostics


def write_swc(path: Path, chains: Sequence[Sequence[tuple[Vec, float, int]]], comments: Sequence[str], box_min: Vec, box_max: Vec) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("# UG4_TERMINALS_ON_OUTER_WALL 1\n")
        out.write(
            "# UG4_OUTER_BOX "
            + " ".join(f"{value:.17g}" for value in (
                box_min[0], box_max[0], box_min[1], box_max[1], box_min[2], box_max[2]
            ))
            + "\n"
        )
        out.write("# Terminal extensions use low-bending cubic transitions and straight wall-normal tails.\n")
        for comment in comments:
            out.write(comment + "\n")
        next_id = 1
        for chain in chains:
            parent = -1
            for index, (point, radius, kind) in enumerate(chain):
                # Neuro Collection interprets type 1 as a root/soma marker.
                # Only the first node of each chain may retain it; otherwise a
                # root extension made from the original root's kind would be
                # skipped as a run of soma points.
                node_kind = 1 if index == 0 else (3 if kind in (0, 1) else kind)
                out.write(
                    f"{next_id} {node_kind} {point[0]:.17g} {point[1]:.17g} "
                    f"{point[2]:.17g} {radius:.17g} {parent}\n"
                )
                parent = next_id
                next_id += 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("--padding", required=True, type=float)
    parser.add_argument("--padding-x", type=float, default=None,
                        help="X clearance; defaults to --padding")
    parser.add_argument("--padding-y", type=float, default=None,
                        help="Y clearance; defaults to --padding")
    parser.add_argument("--padding-z", type=float, default=None,
                        help="Z clearance; defaults to --padding")
    parser.add_argument("--spacing", required=True, type=float)
    parser.add_argument("--box-refs", type=int, default=4,
                        help="Outer-box refinement level used to center each terminal in its wall cell")
    parser.add_argument("--min-padding-factor", type=float, default=25.0,
                        help="minimum box padding as a multiple of the largest outer radius")
    args = parser.parse_args()
    requested_padding = (
        args.padding if args.padding_x is None else args.padding_x,
        args.padding if args.padding_y is None else args.padding_y,
        args.padding if args.padding_z is None else args.padding_z,
    )
    if (args.padding <= 0.0 or any(value <= 0.0 for value in requested_padding)
            or args.spacing <= 0.0 or args.min_padding_factor < 0.0):
        parser.error("padding/spacing must be positive and min-padding-factor nonnegative")
    if args.box_refs < 0 or args.box_refs > 10:
        parser.error("box-refs must be between 0 and 10")

    nodes, comments = read_swc(args.input)
    chains = ordered_chains(nodes)
    largest_radius = max(node.radius for node in nodes)
    safety_padding = args.min_padding_factor * largest_radius
    effective_padding = tuple(max(value, safety_padding) for value in requested_padding)
    if effective_padding != requested_padding:
        print(
            f"Raised terminal box padding from {requested_padding} to "
            f"{effective_padding} mesh units to limit terminal curvature."
        )
    box_min, box_max = fixed_box(nodes, effective_padding)
    extended = []
    diagnostics = []
    for index, chain in enumerate(chains, start=1):
        result, messages = extend_chain(
            chain, box_min, box_max, args.spacing, args.box_refs, f"nephron {index}"
        )
        extended.append(result)
        diagnostics.extend(messages)
    write_swc(args.output, extended, comments, box_min, box_max)
    print(f"Fixed OuterWall bounds: min={box_min}, max={box_max}")
    for message in diagnostics:
        print(message)
    print(f"Extended SWC: {args.output}")


if __name__ == "__main__":
    main()
