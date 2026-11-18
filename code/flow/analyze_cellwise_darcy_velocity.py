#!/usr/bin/env python3
"""Audit element-centered Darcy velocity directly from VTK CellData."""

from __future__ import annotations

import argparse
import base64
import math
import statistics
import struct
import xml.etree.ElementTree as ET


FORMATS = {
    "Float32": ("f", 4),
    "Float64": ("d", 8),
    "Int32": ("i", 4),
    "Int64": ("q", 8),
    "UInt8": ("B", 1),
    "UInt32": ("I", 4),
    "UInt64": ("Q", 8),
}


def decode(array):
    text = "".join((array.text or "").split())
    if array.get("format") != "binary":
        return tuple(float(x) for x in text.split())
    raw = base64.b64decode(text)
    # UG4's inline VTK writer uses a UInt32 byte-count prefix.
    size = struct.unpack("<I", raw[:4])[0]
    payload = raw[4:4 + size]
    code, width = FORMATS[array.get("type")]
    return struct.unpack(f"<{len(payload) // width}{code}", payload)


def percentile(values, fraction):
    ordered = sorted(values)
    x = fraction * (len(ordered) - 1)
    lo = int(math.floor(x))
    hi = int(math.ceil(x))
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - x) + ordered[hi] * (x - lo)


def rms(values):
    return math.sqrt(sum(x * x for x in values) / len(values))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    parser.add_argument("--field", default="darcy_velocity_mm_s")
    args = parser.parse_args()

    piece = ET.parse(args.file).getroot().find("./UnstructuredGrid/Piece")
    points_raw = decode(piece.find("./Points/DataArray"))
    points = [points_raw[i:i + 3] for i in range(0, len(points_raw), 3)]
    arrays = {a.get("Name"): a for a in piece.findall("./Cells/DataArray")}
    connectivity = tuple(map(int, decode(arrays["connectivity"])))
    offsets = tuple(map(int, decode(arrays["offsets"])))
    velocity_array = next(
        a for a in piece.findall("./CellData/DataArray")
        if a.get("Name") == args.field
    )
    raw_velocity = decode(velocity_array)
    velocities = [raw_velocity[i:i + 3] for i in range(0, len(raw_velocity), 3)]
    if len(velocities) != len(offsets):
        raise SystemExit("CellData length does not match the number of cells")

    centroids = []
    begin = 0
    for end in offsets:
        vertices = [points[i] for i in connectivity[begin:end]]
        centroids.append(tuple(sum(p[d] for p in vertices) / len(vertices) for d in range(3)))
        begin = end

    speed = []
    radial = []
    tangential = []
    axial = []
    for (x, y, _), (vx, vy, vz) in zip(centroids, velocities):
        radius = math.hypot(x, y)
        if radius == 0.0:
            vr = vt = 0.0
        else:
            vr = (x * vx + y * vy) / radius
            vt = (-y * vx + x * vy) / radius
        speed.append(math.sqrt(vx * vx + vy * vy + vz * vz))
        radial.append(vr)
        tangential.append(vt)
        axial.append(vz)

    nonfinite = sum(not math.isfinite(x) for v in velocities for x in v)
    print(f"file: {args.file}")
    print(f"cell-centered vectors: {len(velocities)}; nonfinite components: {nonfinite}")
    print(
        "speed [mm/s]: "
        f"min={min(speed):.9e}, p50={percentile(speed, .5):.9e}, "
        f"p95={percentile(speed, .95):.9e}, max={max(speed):.9e}"
    )
    print(
        "cylindrical components [mm/s]: "
        f"radial mean={statistics.fmean(radial):.9e} rms={rms(radial):.9e}; "
        f"tangential mean={statistics.fmean(tangential):.9e} rms={rms(tangential):.9e}; "
        f"axial mean={statistics.fmean(axial):.9e} rms={rms(axial):.9e}"
    )
    transverse = math.sqrt(rms(radial) ** 2 + rms(tangential) ** 2)
    print(
        "directional roughness indicators: "
        f"tangential/transverse_rms={rms(tangential) / max(transverse, 1e-300):.9e}, "
        f"tangential/speed_rms={rms(tangential) / max(rms(speed), 1e-300):.9e}"
    )


if __name__ == "__main__":
    main()
