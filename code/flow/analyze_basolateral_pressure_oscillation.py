#!/usr/bin/env python3
"""Measure nodal Darcy-pressure roughness on a straight Basolateral cylinder."""

from __future__ import annotations

import argparse
import collections
import math

import numpy as np

from analyze_refined_pressure_mode import read_vtu


def rms(values):
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values))) if len(values) else float("nan")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vtu")
    parser.add_argument("--pressure-name", default="darcy_pressure_kinematic")
    parser.add_argument("--radius", type=float, default=0.017)
    parser.add_argument("--radius-tolerance", type=float, default=2.0e-6)
    parser.add_argument("--pressure-to-pa", type=float, default=1.0e-3)
    args = parser.parse_args()

    points, pressure, _ = read_vtu(args.vtu, args.pressure_name)
    radii = np.hypot(points[:, 0], points[:, 1])
    selected = np.flatnonzero(np.abs(radii - args.radius) <= args.radius_tolerance)
    if not len(selected):
        raise SystemExit("no points found on requested Basolateral radius")

    z_span = float(np.ptp(points[selected, 2]))
    z_tolerance = max(z_span * 1.0e-8, 1.0e-12)
    rings = collections.defaultdict(list)
    for index in selected:
        rings[int(round(points[index, 2] / z_tolerance))].append(index)
    ordered_rings = sorted(rings.values(), key=lambda ids: np.mean(points[ids, 2]))
    ordered_rings = [ids for ids in ordered_rings if len(ids) >= 8]

    ring_z = np.asarray([np.mean(points[ids, 2]) for ids in ordered_rings])
    ring_mean = np.asarray([np.mean(pressure[ids]) for ids in ordered_rings])
    circum_residuals = np.concatenate(
        [pressure[ids] - np.mean(pressure[ids]) for ids in ordered_rings]
    )
    circum_std = np.asarray([np.std(pressure[ids]) for ids in ordered_rings])

    normalized_z = 2.0 * (ring_z - ring_z.min()) / np.ptp(ring_z) - 1.0
    axial_fit = np.polynomial.legendre.legfit(normalized_z, ring_mean, deg=3)
    axial_residual = ring_mean - np.polynomial.legendre.legval(
        normalized_z, axial_fit
    )
    middle = (normalized_z >= -0.6) & (normalized_z <= 0.6)
    middle_fit = np.polynomial.legendre.legfit(
        normalized_z[middle], ring_mean[middle], deg=3
    )
    middle_axial_residual = ring_mean[middle] - np.polynomial.legendre.legval(
        normalized_z[middle], middle_fit
    )
    axial_signs = np.asarray([(-1.0) ** i for i in range(middle.sum())])
    axial_alternating = abs(
        float(np.dot(middle_axial_residual, axial_signs) / middle.sum())
    )

    nyquist_amplitudes = []
    for ids in ordered_rings:
        ids = sorted(ids, key=lambda i: math.atan2(points[i, 1], points[i, 0]))
        values = pressure[ids] - np.mean(pressure[ids])
        signs = np.asarray([(-1.0) ** i for i in range(len(ids))])
        nyquist_amplitudes.append(abs(float(np.dot(values, signs) / len(ids))))
    nyquist_amplitudes = np.asarray(nyquist_amplitudes)

    scale = args.pressure_to_pa
    print(f"file: {args.vtu}")
    print(f"Basolateral points={len(selected)}, axial rings={len(ordered_rings)}")
    print(
        "pressure range on Basolateral [Pa]: "
        f"{pressure[selected].min() * scale:.9e} .. "
        f"{pressure[selected].max() * scale:.9e}"
    )
    print(
        "circumferential residual after ring mean [Pa]: "
        f"rms={rms(circum_residuals) * scale:.9e}, "
        f"max={np.max(np.abs(circum_residuals)) * scale:.9e}"
    )
    print(
        "ring circumferential std [Pa]: "
        f"rms={rms(circum_std) * scale:.9e}, "
        f"max={circum_std.max() * scale:.9e}"
    )
    print(
        "axial ring-mean residual after cubic trend [Pa]: "
        f"rms={rms(axial_residual) * scale:.9e}, "
        f"max={np.max(np.abs(axial_residual)) * scale:.9e}"
    )
    print(
        "middle-60% axial residual after cubic trend [Pa]: "
        f"rms={rms(middle_axial_residual) * scale:.9e}, "
        f"max={np.max(np.abs(middle_axial_residual)) * scale:.9e}, "
        f"alternating={axial_alternating * scale:.9e}"
    )
    print(
        "circumferential alternating-mode amplitude [Pa]: "
        f"rms={rms(nyquist_amplitudes) * scale:.9e}, "
        f"max={nyquist_amplitudes.max() * scale:.9e}"
    )


if __name__ == "__main__":
    main()
