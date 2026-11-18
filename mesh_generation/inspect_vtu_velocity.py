#!/usr/bin/env python3
"""Report a point-velocity field versus the dominant coordinate axis."""

from __future__ import annotations

import argparse
import base64
import glob
import math
import struct
import xml.etree.ElementTree as ET


FORMATS = {
    "Float32": ("f", 4),
    "Float64": ("d", 8),
    "Int32": ("i", 4),
    "UInt32": ("I", 4),
    "UInt8": ("B", 1),
}


def decode_array(array):
    if array.get("format") != "binary":
        return [float(x) for x in (array.text or "").split()]
    encoded = "".join((array.text or "").split())
    raw = base64.b64decode(encoded)
    size = struct.unpack("<I", raw[:4])[0]
    payload = raw[4:4 + size]
    code, width = FORMATS[array.get("type")]
    if len(payload) % width:
        raise ValueError("Binary VTK payload has an invalid byte count")
    return struct.unpack(f"<{len(payload) // width}{code}", payload)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    parser.add_argument("--field", default="stokes_velocity_mm_s")
    parser.add_argument("--bins", type=int, default=20)
    parser.add_argument("--center", nargs=2, type=float, metavar=("X", "Y"))
    args = parser.parse_args()

    paths = []
    for pattern in args.files:
        paths.extend(glob.glob(pattern))
    paths = sorted(set(paths))
    if not paths:
        raise SystemExit("No matching VTU files")

    points = []
    velocities = []
    for path in paths:
        root = ET.parse(path).getroot()
        piece = root.find("./UnstructuredGrid/Piece")
        point_array = piece.find("./Points/DataArray")
        velocity_array = next(
            a for a in piece.findall("./PointData/DataArray")
            if a.get("Name") == args.field
        )
        xyz = decode_array(point_array)
        uvw = decode_array(velocity_array)
        points.extend(tuple(xyz[i:i + 3]) for i in range(0, len(xyz), 3))
        velocities.extend(tuple(uvw[i:i + 3]) for i in range(0, len(uvw), 3))

    spans = [max(p[a] for p in points) - min(p[a] for p in points) for a in range(3)]
    axis = max(range(3), key=spans.__getitem__)
    labels = "xyz"
    lo = min(p[axis] for p in points)
    hi = max(p[axis] for p in points)
    mags = [math.sqrt(sum(c * c for c in v)) for v in velocities]
    nonzero = sum(m > 1.0e-10 for m in mags)

    print(f"files: {len(paths)}")
    print(f"points including partition duplicates: {len(points)}")
    print(f"coordinate spans: x={spans[0]:.9g}, y={spans[1]:.9g}, z={spans[2]:.9g}")
    print(f"dominant axis: {labels[axis]}, range=[{lo:.9g}, {hi:.9g}]")
    print(f"speed: min={min(mags):.12g}, max={max(mags):.12g}, nonzero(>1e-10)={nonzero}/{len(mags)}")
    suffix = " max_radius" if args.center and axis == 2 else ""
    print("bin_center count nonzero min_speed mean_speed max_speed" + suffix)
    for b in range(args.bins):
        left = lo + (hi - lo) * b / args.bins
        right = lo + (hi - lo) * (b + 1) / args.bins
        ids = [
            i for i, p in enumerate(points)
            if (left <= p[axis] < right) or (b == args.bins - 1 and p[axis] == right)
        ]
        values = [mags[i] for i in ids]
        nz = sum(v > 1.0e-10 for v in values)
        center = 0.5 * (left + right)
        radius_text = ""
        if args.center and axis == 2:
            cx, cy = args.center
            max_radius = max(math.hypot(points[i][0] - cx, points[i][1] - cy) for i in ids)
            radius_text = f" {max_radius:.9g}"
        print(
            f"{center:.9g} {len(values)} {nz} "
            f"{min(values):.9g} {sum(values) / len(values):.9g} {max(values):.9g}"
            + radius_text
        )


if __name__ == "__main__":
    main()
