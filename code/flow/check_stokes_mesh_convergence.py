#!/usr/bin/env python3
"""Compare nested Stokes VTU solutions and enforce a mesh-convergence rule.

The comparison is performed at vertices shared by adjacent refinement levels,
so it does not introduce an interpolation error.  Pressure differences are
made gauge independent by removing the best constant offset.  Exact assembled
fluxes and pressure drops are supplied from the Lua diagnostic log because a
nodal VTU alone cannot reconstruct the FV1 boundary assembly exactly.
"""

from __future__ import annotations

import argparse
import base64
import math
import struct
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


FORMATS = {
    "Float32": ("f", 4),
    "Float64": ("d", 8),
    "Int32": ("i", 4),
    "UInt32": ("I", 4),
    "UInt8": ("B", 1),
}


def decode_array(array: ET.Element) -> tuple[float, ...]:
    if array.get("format") != "binary":
        return tuple(float(value) for value in (array.text or "").split())
    raw = base64.b64decode("".join((array.text or "").split()))
    size = struct.unpack("<I", raw[:4])[0]
    payload = raw[4:4 + size]
    code, width = FORMATS[array.get("type")]
    if len(payload) % width:
        raise ValueError("Binary VTK payload has an invalid byte count")
    return struct.unpack(f"<{len(payload) // width}{code}", payload)


def load_vtu(path: Path, coordinate_digits: int):
    root = ET.parse(path).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError(f"No UnstructuredGrid Piece in {path}")
    xyz = decode_array(piece.find("./Points/DataArray"))
    arrays = {
        array.get("Name"): decode_array(array)
        for array in piece.findall("./PointData/DataArray")
    }
    pressure = arrays["stokes_pressure_kinematic"]
    velocity_flat = arrays["stokes_velocity_mm_s"]
    points = tuple(tuple(xyz[i:i + 3]) for i in range(0, len(xyz), 3))
    velocity = tuple(
        tuple(velocity_flat[i:i + 3])
        for i in range(0, len(velocity_flat), 3)
    )
    keys = tuple(
        tuple(round(component, coordinate_digits) for component in point)
        for point in points
    )
    if len(set(keys)) != len(keys):
        raise ValueError(
            f"Coordinate rounding produced duplicate vertices in {path}; "
            "increase --coordinate-digits"
        )
    return points, velocity, pressure, {key: i for i, key in enumerate(keys)}


def relative_change(old: float, new: float) -> float:
    return abs(new - old) / max(abs(new), 1.0e-300)


def relative_l2(differences, reference) -> float:
    numerator = math.sqrt(sum(value * value for value in differences))
    denominator = math.sqrt(sum(value * value for value in reference))
    return numerator / max(denominator, 1.0e-300)


def compare_levels(coarse, fine, midplane: float, plane_tolerance: float):
    c_points, c_velocity, c_pressure, c_index = coarse
    f_points, f_velocity, f_pressure, f_index = fine
    common = sorted(set(c_index).intersection(f_index))
    if len(common) != len(c_points):
        raise ValueError(
            f"Only {len(common)} of {len(c_points)} coarse vertices occur in "
            "the finer VTU; inputs are not a nested refinement sequence"
        )

    velocity_difference = []
    velocity_reference = []
    pressure_difference = []
    pressure_reference = []
    pressure_offset = sum(
        f_pressure[f_index[key]] - c_pressure[c_index[key]] for key in common
    ) / len(common)
    for key in common:
        ci, fi = c_index[key], f_index[key]
        velocity_difference.extend(
            f_velocity[fi][component] - c_velocity[ci][component]
            for component in range(3)
        )
        velocity_reference.extend(f_velocity[fi])
        pressure_difference.append(
            f_pressure[fi] - c_pressure[ci] - pressure_offset
        )
        pressure_reference.append(f_pressure[fi])

    pressure_mean = sum(pressure_reference) / len(pressure_reference)
    centered_pressure = [value - pressure_mean for value in pressure_reference]
    plane_keys = [
        key for key in common if abs(c_points[c_index[key]][2] - midplane)
        <= plane_tolerance
    ]
    if not plane_keys:
        raise ValueError(f"No common vertices found on z={midplane:g}")
    plane_difference = []
    plane_reference = []
    for key in plane_keys:
        ci, fi = c_index[key], f_index[key]
        plane_difference.extend(
            f_velocity[fi][component] - c_velocity[ci][component]
            for component in range(3)
        )
        plane_reference.extend(f_velocity[fi])
    center_key = min(
        plane_keys,
        key=lambda key: math.hypot(
            c_points[c_index[key]][0], c_points[c_index[key]][1]
        ),
    )
    ci, fi = c_index[center_key], f_index[center_key]
    return {
        "common_points": len(common),
        "velocity_l2": relative_l2(velocity_difference, velocity_reference),
        "pressure_l2_gauge_free": relative_l2(
            pressure_difference, centered_pressure
        ),
        "midplane_points": len(plane_keys),
        "midplane_velocity_l2": relative_l2(
            plane_difference, plane_reference
        ),
        "centerline_axial_old": c_velocity[ci][2],
        "centerline_axial_new": f_velocity[fi][2],
        "centerline_axial_change": relative_change(
            c_velocity[ci][2], f_velocity[fi][2]
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path,
                        help="VTUs ordered from coarsest to finest")
    parser.add_argument("--pressure-drops",
                        help="comma-separated assembled pressure drops")
    parser.add_argument("--flows",
                        help="comma-separated assembled inlet-flow magnitudes")
    parser.add_argument("--tolerance", type=float, default=0.01)
    parser.add_argument("--coordinate-digits", type=int, default=8)
    parser.add_argument("--midplane", type=float, default=0.5)
    parser.add_argument("--plane-tolerance", type=float, default=1.0e-7)
    parser.add_argument("--require-converged", action="store_true")
    args = parser.parse_args()
    args.pressure_drops = (
        [float(value) for value in args.pressure_drops.split(",")]
        if args.pressure_drops else None
    )
    args.flows = (
        [float(value) for value in args.flows.split(",")]
        if args.flows else None
    )
    if len(args.files) < 2:
        parser.error("at least two ordered VTUs are required")
    for name, values in (("pressure drops", args.pressure_drops),
                         ("flows", args.flows)):
        if values is not None and len(values) != len(args.files):
            parser.error(f"{name} must contain one value per VTU")

    levels = [load_vtu(path, args.coordinate_digits) for path in args.files]
    comparisons = [
        compare_levels(levels[i - 1], levels[i], args.midplane,
                       args.plane_tolerance)
        for i in range(1, len(levels))
    ]
    scalar_changes = []
    for i, comparison in enumerate(comparisons, 1):
        print(f"level {i - 1} -> {i}")
        print(f"  shared vertices:             {comparison['common_points']}")
        print(f"  common-vertex velocity L2:   {comparison['velocity_l2']:.6%}")
        print("  gauge-free pressure L2:      "
              f"{comparison['pressure_l2_gauge_free']:.6%}")
        print(f"  midplane shared vertices:    {comparison['midplane_points']}")
        print("  midplane velocity L2:        "
              f"{comparison['midplane_velocity_l2']:.6%}")
        print("  centerline axial velocity:   "
              f"{comparison['centerline_axial_old']:.12g} -> "
              f"{comparison['centerline_axial_new']:.12g} mm/s "
              f"({comparison['centerline_axial_change']:.6%})")
        if args.pressure_drops:
            change = relative_change(args.pressure_drops[i - 1],
                                     args.pressure_drops[i])
            scalar_changes.append(("pressure drop", change))
            print(f"  assembled pressure-drop change: {change:.6%}")
        if args.flows:
            change = relative_change(args.flows[i - 1], args.flows[i])
            scalar_changes.append(("assembled flow", change))
            print(f"  assembled-flow change:       {change:.6%}")

    latest = comparisons[-1]
    checks = [
        ("common-vertex velocity L2", latest["velocity_l2"]),
        ("midplane velocity L2", latest["midplane_velocity_l2"]),
        ("centerline axial velocity", latest["centerline_axial_change"]),
    ]
    if args.pressure_drops:
        checks.append(("pressure drop", relative_change(
            args.pressure_drops[-2], args.pressure_drops[-1]
        )))
    if args.flows:
        checks.append(("assembled flow", relative_change(
            args.flows[-2], args.flows[-1]
        )))
    failed = [(name, value) for name, value in checks
              if value > args.tolerance]
    enough_levels = len(args.files) >= 3
    converged = enough_levels and not failed
    print(f"mesh convergence tolerance:    {args.tolerance:.6%}")
    print("mesh convergence:              " +
          ("PASS" if converged else "NOT ESTABLISHED"))
    if not enough_levels:
        print("  reason: fewer than three refinement levels")
    for name, value in failed:
        print(f"  above tolerance: {name} = {value:.6%}")
    return 2 if args.require_converged and not converged else 0


if __name__ == "__main__":
    sys.exit(main())
