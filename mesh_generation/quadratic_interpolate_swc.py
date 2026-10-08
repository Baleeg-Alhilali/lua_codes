#!/usr/bin/env python3
"""Smooth and densely resample one unbranched SWC centerline.

Pipeline role
-------------
The Lua mesh workflows call this program before invoking UG4. A raw traced
centerline can contain short segments, isolated outliers, and abrupt direction
changes. Those defects produce unstable tube frames, sharp coarse cells, or
TetGen surface-recovery failures. This program creates the ``*_mesh_ready.swc``
used by the C++ neurite importer.

Processing stages
-----------------
1. Parse standard seven-column SWC records and validate parent references.
2. Require exactly one root and no branching; the nephron workflows treat each
   file as one continuous root-to-terminal trace.
3. Optionally soften control points whose direction-change angle exceeds a
   threshold. Repeated relaxation moves those points toward their neighbors.
4. Interpolate the adjusted controls:
   * order 1: linear interpolation
   * order 2: quadratic Bezier interpolation
   * order 3: Catmull-Rom cubic interpolation (production default)
5. Optionally resample the interpolated curve to nearly uniform arc-length
   spacing. Uniform spacing makes the requested UG4 anisotropy predictable.
6. Write sequential SWC IDs and a simple parent chain while preserving comments.
7. Report sharp-turn and spacing metrics for quality control.

The script changes only the trace representation. Physical Lumen and Membrane
radii are imposed later by the Lua pipeline.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path


DEFAULT_INPUT = Path("/Users/alhilaba/Desktop/1-4-range-filtered.swc")


@dataclass(frozen=True)
class SwcNode:
    """One immutable record exactly as represented by the SWC file."""
    node_id: int
    node_type: int
    x: float
    y: float
    z: float
    radius: float
    parent_id: int


@dataclass(frozen=True)
class Point:
    """Geometry/control-point representation used during interpolation."""
    x: float
    y: float
    z: float
    radius: float
    node_type: int


def read_swc(path: Path) -> tuple[list[str], dict[int, SwcNode]]:
    """Read comments and validated seven-column nodes keyed by node ID."""
    if not path.exists():
        raise FileNotFoundError(f"SWC file not found: {path}")

    comments: list[str] = []
    nodes: dict[int, SwcNode] = {}

    with path.open("r", encoding="utf-8") as swc_file:
        for line_number, raw_line in enumerate(swc_file, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                comments.append(raw_line.rstrip("\n"))
                continue

            parts = line.split()
            if len(parts) < 7:
                raise ValueError(
                    f"Line {line_number} has {len(parts)} columns; expected 7"
                )

            try:
                node = SwcNode(
                    node_id=int(parts[0]),
                    node_type=int(parts[1]),
                    x=float(parts[2]),
                    y=float(parts[3]),
                    z=float(parts[4]),
                    radius=float(parts[5]),
                    parent_id=int(parts[6]),
                )
            except ValueError as exc:
                raise ValueError(f"Could not parse line {line_number}: {line}") from exc

            if node.node_id in nodes:
                raise ValueError(f"Duplicate node id {node.node_id} on line {line_number}")
            nodes[node.node_id] = node

    if not nodes:
        raise ValueError(f"No SWC nodes found in {path}")

    return comments, nodes


def build_children(nodes: dict[int, SwcNode]) -> dict[int, list[int]]:
    """Build parent-to-children adjacency and reject missing parent IDs."""
    children = {node_id: [] for node_id in nodes}
    for node in nodes.values():
        if node.parent_id == -1:
            continue
        if node.parent_id not in nodes:
            raise ValueError(
                f"Node {node.node_id} references missing parent {node.parent_id}"
            )
        children[node.parent_id].append(node.node_id)
    return children


def single_trace_order(nodes: dict[int, SwcNode]) -> list[int]:
    """Return node ids for a single unbranched root-to-end trace."""
    children = build_children(nodes)
    roots = [node.node_id for node in nodes.values() if node.parent_id == -1]
    if len(roots) != 1:
        raise ValueError(f"Expected one SWC root, found {len(roots)}")

    branch_nodes = [node_id for node_id, child_ids in children.items() if len(child_ids) > 1]
    if branch_nodes:
        raise ValueError(
            "This quadratic interpolation script expects one unbranched path. "
            f"Branching node ids found: {branch_nodes[:10]}"
        )

    order = [roots[0]]
    current = roots[0]
    while children[current]:
        current = children[current][0]
        order.append(current)

    return order


def as_point(node: SwcNode) -> Point:
    return Point(node.x, node.y, node.z, node.radius, node.node_type)


def lerp(a: Point, b: Point, weight: float) -> Point:
    return Point(
        x=a.x + (b.x - a.x) * weight,
        y=a.y + (b.y - a.y) * weight,
        z=a.z + (b.z - a.z) * weight,
        radius=a.radius + (b.radius - a.radius) * weight,
        node_type=a.node_type,
    )


def quadratic_bezier(start: Point, control: Point, end: Point, t: float) -> Point:
    """Evaluate a 2nd-order Bezier polynomial at t in [0, 1]."""
    one_minus_t = 1.0 - t
    start_weight = one_minus_t * one_minus_t
    control_weight = 2.0 * one_minus_t * t
    end_weight = t * t
    return Point(
        x=start_weight * start.x + control_weight * control.x + end_weight * end.x,
        y=start_weight * start.y + control_weight * control.y + end_weight * end.y,
        z=start_weight * start.z + control_weight * control.z + end_weight * end.z,
        radius=(
            start_weight * start.radius
            + control_weight * control.radius
            + end_weight * end.radius
        ),
        node_type=control.node_type,
    )


def catmull_rom(
    previous_point: Point,
    start_point: Point,
    end_point: Point,
    next_point: Point,
    t: float,
) -> Point:
    """Evaluate a 3rd-order Catmull-Rom polynomial at t in [0, 1]."""
    t2 = t * t
    t3 = t2 * t

    def cubic(previous_value: float, start_value: float, end_value: float, next_value: float) -> float:
        return 0.5 * (
            (2.0 * start_value)
            + (-previous_value + end_value) * t
            + (2.0 * previous_value - 5.0 * start_value + 4.0 * end_value - next_value) * t2
            + (-previous_value + 3.0 * start_value - 3.0 * end_value + next_value) * t3
        )

    return Point(
        x=cubic(previous_point.x, start_point.x, end_point.x, next_point.x),
        y=cubic(previous_point.y, start_point.y, end_point.y, next_point.y),
        z=cubic(previous_point.z, start_point.z, end_point.z, next_point.z),
        radius=max(
            0.0,
            cubic(
                previous_point.radius,
                start_point.radius,
                end_point.radius,
                next_point.radius,
            ),
        ),
        node_type=start_point.node_type,
    )


def distance(a: Point, b: Point) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def turn_angle(previous_point: Point, point: Point, next_point: Point) -> float:
    """Return direction-change angle in degrees, where 0 is straight."""
    incoming = (
        point.x - previous_point.x,
        point.y - previous_point.y,
        point.z - previous_point.z,
    )
    outgoing = (
        next_point.x - point.x,
        next_point.y - point.y,
        next_point.z - point.z,
    )
    incoming_norm = math.sqrt(sum(value * value for value in incoming))
    outgoing_norm = math.sqrt(sum(value * value for value in outgoing))
    if incoming_norm == 0 or outgoing_norm == 0:
        return 180.0

    dot = sum(incoming[index] * outgoing[index] for index in range(3))
    cosine = max(-1.0, min(1.0, dot / (incoming_norm * outgoing_norm)))
    return math.degrees(math.acos(cosine))


def count_sharp_turns(points: list[Point], angle_threshold: float) -> int:
    return sum(
        1
        for index in range(1, len(points) - 1)
        if turn_angle(points[index - 1], points[index], points[index + 1])
        > angle_threshold
    )


def worst_turn_angle(points: list[Point]) -> float:
    if len(points) < 3:
        return 0.0
    return max(
        turn_angle(points[index - 1], points[index], points[index + 1])
        for index in range(1, len(points) - 1)
    )


def segment_spacing_range(points: list[Point]) -> tuple[float, float]:
    spacings = [distance(points[index - 1], points[index]) for index in range(1, len(points))]
    if not spacings:
        return 0.0, 0.0
    return min(spacings), max(spacings)


def append_if_distinct(points: list[Point], point: Point, min_distance: float = 1e-6) -> None:
    if not points or distance(points[-1], point) > min_distance:
        points.append(point)


def soften_problem_points(
    points: list[Point],
    angle_threshold: float,
    strength: float,
    passes: int,
) -> tuple[list[Point], int]:
    """Relax sharp interior controls while keeping both endpoints fixed.

    A point is edited only when its direction-change angle exceeds the selected
    threshold. ``strength`` blends it toward the midpoint of its neighbors, and
    multiple passes distribute a local correction smoothly along the trace.
    """
    """Move sharp interior control points toward local midpoints before interpolation."""
    if passes < 1:
        return points, 0
    if not 0.0 < angle_threshold < 180.0:
        raise ValueError("pre_smooth_angle must be between 0 and 180 degrees")
    if not 0.0 < strength <= 1.0:
        raise ValueError("pre_smooth_strength must be greater than 0 and at most 1")

    softened = list(points)
    total_moved = 0
    for _ in range(passes):
        next_points = list(softened)
        moved_this_pass = 0
        for index in range(1, len(softened) - 1):
            previous_point = softened[index - 1]
            point = softened[index]
            next_point = softened[index + 1]
            angle = turn_angle(previous_point, point, next_point)
            if angle <= angle_threshold:
                continue

            sharpness = (angle - angle_threshold) / (180.0 - angle_threshold)
            weight = strength * sharpness
            midpoint = lerp(previous_point, next_point, 0.5)
            next_points[index] = Point(
                x=point.x + (midpoint.x - point.x) * weight,
                y=point.y + (midpoint.y - point.y) * weight,
                z=point.z + (midpoint.z - point.z) * weight,
                radius=point.radius + (midpoint.radius - point.radius) * weight,
                node_type=point.node_type,
            )
            moved_this_pass += 1

        softened = next_points
        total_moved += moved_this_pass
        if moved_this_pass == 0:
            break

    return softened, total_moved


def interpolate_trace(
    original_points: list[Point],
    corner_fraction: float,
    samples_per_corner: int,
    order: int,
) -> list[Point]:
    """Generate dense samples from the ordered control points.

    Duplicate samples at neighboring interpolation windows are suppressed so
    the output never contains zero-length centerline segments.
    """
    """Round the trace with either 2nd- or 3rd-order polynomial interpolation."""
    if len(original_points) < 3:
        return original_points
    if order not in {1, 2, 3}:
        raise ValueError("order must be 1, 2, or 3")
    if samples_per_corner < 2:
        raise ValueError("samples_per_corner must be 2 or greater")

    if order == 1:
        return original_points

    if order == 3:
        smoothed: list[Point] = [original_points[0]]
        for index in range(len(original_points) - 1):
            previous_point = original_points[index - 1] if index > 0 else original_points[index]
            start_point = original_points[index]
            end_point = original_points[index + 1]
            next_point = (
                original_points[index + 2]
                if index + 2 < len(original_points)
                else original_points[index + 1]
            )

            for sample_index in range(1, samples_per_corner + 1):
                t = sample_index / samples_per_corner
                append_if_distinct(
                    smoothed,
                    catmull_rom(previous_point, start_point, end_point, next_point, t),
                )
        return smoothed

    if not 0.0 < corner_fraction < 0.5:
        raise ValueError("corner_fraction must be greater than 0 and less than 0.5")

    smoothed: list[Point] = [original_points[0]]
    for index in range(1, len(original_points) - 1):
        previous_point = original_points[index - 1]
        control_point = original_points[index]
        next_point = original_points[index + 1]

        before_corner = lerp(control_point, previous_point, corner_fraction)
        after_corner = lerp(control_point, next_point, corner_fraction)
        append_if_distinct(smoothed, before_corner)

        for sample_index in range(1, samples_per_corner + 1):
            t = sample_index / samples_per_corner
            append_if_distinct(
                smoothed,
                quadratic_bezier(before_corner, control_point, after_corner, t),
            )

    append_if_distinct(smoothed, original_points[-1])
    return smoothed


def resample_uniformly(points: list[Point], spacing: float) -> list[Point]:
    """Resample by accumulated arc length while preserving both endpoints."""
    """Resample a polyline at constant arc-length intervals."""
    if spacing <= 0.0:
        raise ValueError("resample_spacing must be greater than 0")
    if len(points) < 2:
        return points

    cumulative = [0.0]
    for index in range(1, len(points)):
        cumulative.append(cumulative[-1] + distance(points[index - 1], points[index]))

    total_length = cumulative[-1]
    if total_length <= 1e-12:
        return [points[0]]

    resampled = [points[0]]
    segment_index = 1
    target_distance = spacing
    while target_distance < total_length:
        while cumulative[segment_index] < target_distance:
            segment_index += 1

        segment_start_distance = cumulative[segment_index - 1]
        segment_length = cumulative[segment_index] - segment_start_distance
        if segment_length > 1e-12:
            weight = (target_distance - segment_start_distance) / segment_length
            append_if_distinct(
                resampled,
                lerp(points[segment_index - 1], points[segment_index], weight),
            )
        target_distance += spacing

    append_if_distinct(resampled, points[-1])
    return resampled


def write_swc(
    output_path: Path,
    comments: list[str],
    points: list[Point],
    source_path: Path,
    corner_fraction: float,
    samples_per_corner: int,
    order: int,
    pre_smooth_angle: float,
    pre_smooth_strength: float,
    pre_smooth_passes: int,
    resample_spacing: float,
) -> None:
    """Write comments followed by a sequential, unbranched SWC parent chain."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as output_file:
        for comment in comments:
            output_file.write(f"{comment}\n")
        output_file.write(
            "# SMOOTHED_BY:\tquadratic_interpolate_swc.py "
            f"source={source_path} order={order} corner_fraction={corner_fraction} "
            f"samples_per_corner={samples_per_corner} "
            f"pre_smooth_angle={pre_smooth_angle} "
            f"pre_smooth_strength={pre_smooth_strength} "
            f"pre_smooth_passes={pre_smooth_passes} "
            f"resample_spacing={resample_spacing}\n"
        )

        for index, point in enumerate(points, start=1):
            parent_id = -1 if index == 1 else index - 1
            output_file.write(
                f"{index:d} {point.node_type:d} "
                f"{point.x:.6f} {point.y:.6f} {point.z:.6f} "
                f"{point.radius:.6f} {parent_id:d}\n"
            )


def parse_args() -> argparse.Namespace:
    """Define CLI controls used directly and by both Lua mesh pipelines."""
    parser = argparse.ArgumentParser(
        description=(
            "Replace an unbranched SWC trace with a denser path made from "
            "polynomial interpolation."
        )
    )
    parser.add_argument(
        "swc_file",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input .swc file. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output .swc path. Default: input name with _cubic before .swc",
    )
    parser.add_argument(
        "--order",
        type=int,
        choices=(1, 2, 3),
        default=3,
        help=(
            "Polynomial order. 3 uses cubic Catmull-Rom interpolation; 2 keeps "
            "quadratic corner rounding; 1 preserves the input polyline for "
            "resampling without another interpolation pass."
        ),
    )
    parser.add_argument(
        "--corner-fraction",
        type=float,
        default=0.45,
        help="For --order 2 only: how far from each vertex to start/end the rounded curve. Use < 0.5",
    )
    parser.add_argument(
        "--samples-per-corner",
        type=int,
        default=10,
        help="Number of interpolated samples for each rounded vertex",
    )
    parser.add_argument(
        "--pre-smooth-angle",
        type=float,
        default=0.0,
        help=(
            "Before interpolation, soften control points sharper than this angle. "
            "Use 0 to disable."
        ),
    )
    parser.add_argument(
        "--pre-smooth-strength",
        type=float,
        default=0.5,
        help="How strongly to move detected problem points toward neighbor midpoints, 0-1",
    )
    parser.add_argument(
        "--pre-smooth-passes",
        type=int,
        default=0,
        help="Number of problem-point correction passes before interpolation",
    )
    parser.add_argument(
        "--resample-spacing",
        type=float,
        default=0.0,
        help=(
            "Uniformly resample the final trace at this arc-length spacing. "
            "Use 0 to preserve the existing interpolation output."
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Execute validation, smoothing, interpolation, resampling, and reporting."""
    args = parse_args()
    output = args.output
    if output is None:
        suffix = "quadratic" if args.order == 2 else "cubic"
        output = args.swc_file.with_name(f"{args.swc_file.stem}_{suffix}.swc")

    comments, nodes = read_swc(args.swc_file)
    order = single_trace_order(nodes)
    original_points = [as_point(nodes[node_id]) for node_id in order]
    control_points = original_points
    pre_smooth_moved = 0
    sharp_before = 0
    sharp_after = 0
    if args.pre_smooth_passes > 0 and args.pre_smooth_angle > 0.0:
        sharp_before = count_sharp_turns(control_points, args.pre_smooth_angle)
        control_points, pre_smooth_moved = soften_problem_points(
            control_points,
            angle_threshold=args.pre_smooth_angle,
            strength=args.pre_smooth_strength,
            passes=args.pre_smooth_passes,
        )
        sharp_after = count_sharp_turns(control_points, args.pre_smooth_angle)

    smoothed_points = interpolate_trace(
        control_points,
        corner_fraction=args.corner_fraction,
        samples_per_corner=args.samples_per_corner,
        order=args.order,
    )
    if args.resample_spacing < 0.0:
        raise ValueError("resample_spacing must be 0 or greater")
    if args.resample_spacing > 0.0:
        smoothed_points = resample_uniformly(smoothed_points, args.resample_spacing)

    write_swc(
        output,
        comments,
        smoothed_points,
        source_path=args.swc_file,
        corner_fraction=args.corner_fraction,
        samples_per_corner=args.samples_per_corner,
        order=args.order,
        pre_smooth_angle=args.pre_smooth_angle,
        pre_smooth_strength=args.pre_smooth_strength,
        pre_smooth_passes=args.pre_smooth_passes,
        resample_spacing=args.resample_spacing,
    )

    print(f"Read original nodes: {len(original_points)}")
    if args.pre_smooth_passes > 0 and args.pre_smooth_angle > 0.0:
        print(f"Problem turns before pre-smoothing: {sharp_before}")
        print(f"Problem turns after pre-smoothing:  {sharp_after}")
        print(f"Control-point edits across passes: {pre_smooth_moved}")
    print(f"Wrote interpolated nodes: {len(smoothed_points)}")
    minimum_spacing, maximum_spacing = segment_spacing_range(smoothed_points)
    print(f"Output segment spacing: {minimum_spacing:.6g} to {maximum_spacing:.6g}")
    print(f"Worst output turn: {worst_turn_angle(smoothed_points):.6g} degrees")
    print(f"Saved order-{args.order} interpolated SWC to: {output}")


if __name__ == "__main__":
    main()
