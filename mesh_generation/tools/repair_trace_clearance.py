#!/usr/bin/env python3
"""Repair a manually traced SWC centreline so a physical tube fits along it.

Purpose
-------
A hand-clicked centreline at true voxel scale can pass closer to itself than
the tube diameter and can contain corners much sharper than the tube radius.
The mesh pipeline's geometry preflight then rejects the requested radius, and
the historical workaround was to enlarge the whole centreline
(``-centerlineScale 24``), which made the nephron 24 times too long.

This tool keeps the true scale and changes the path only where it violates
the tube geometry:

* non-local clearance: points farther apart than ``--exclusion-arclength``
  along the trace must be at least ``--clearance`` apart in space;
* curvature: the three-point circumradius at the output sampling must be at
  least ``--min-curvature-radius``.

Algorithm (position-based relaxation on a densely resampled polyline):

1. Convert SWC units to micrometres with ``--voxel-size`` and resample the
   polyline at ``--working-spacing``.
2. Repeat until both constraints hold on the output sampling:
   a. push every violating point pair apart along its connecting line;
      the displacement field is smoothed along the curve so neighbouring
      points move together instead of creating kinks;
   b. relax points whose local circumradius is too small toward the midpoint
      of their neighbours;
   c. apply a weak global bending (Laplacian) step and re-space the points
      uniformly;
   d. pull free points weakly back toward the original trace at the same
      normalized arc length, so corrections stay local.  The endpoints are fixed; the length is free and is
      reported, because pinned endpoints plus fixed length forced buckling.
3. Resample at ``--output-spacing`` and audit with the same formulas as the
   Lua pipeline preflight.

The output is in micrometres, so the mesh pipeline must be run with
``-centerlineScale 1`` (and ``-skipSmoothing 1`` to keep this geometry).
To include the OuterWall terminal extensions, repair the extender's output
with ``--frozen-head``/``--frozen-tail`` covering the straight wall-normal
tails; ``# UG4_`` box metadata is carried over and the pipeline is then run
with ``-extendTerminals 0``.
The header records the source, parameters, and displacement/length metrics.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# SWC input/output
# ---------------------------------------------------------------------------

def read_unbranched_swc(path: Path) -> np.ndarray:
    """Return root-to-terminal points of an unbranched SWC (in SWC units)."""
    nodes: dict[int, tuple[float, float, float, int]] = {}
    for line in path.read_text().splitlines():
        payload = line.split("#", 1)[0].split()
        if not payload:
            continue
        if len(payload) != 7:
            raise ValueError(f"SWC line must have 7 columns: {line!r}")
        nid, parent = int(payload[0]), int(payload[6])
        nodes[nid] = (float(payload[2]), float(payload[3]),
                      float(payload[4]), parent)
    roots = [n for n, v in nodes.items() if v[3] == -1]
    if len(roots) != 1:
        raise ValueError(f"expected exactly one root, found {len(roots)}")
    children: dict[int, list[int]] = {}
    for n, v in nodes.items():
        if v[3] != -1:
            children.setdefault(v[3], []).append(n)
    if any(len(c) > 1 for c in children.values()):
        raise ValueError("branched SWC traces are not supported")
    order = [roots[0]]
    while order[-1] in children:
        order.append(children[order[-1]][0])
    if len(order) != len(nodes):
        raise ValueError("SWC contains disconnected nodes")
    return np.array([nodes[n][:3] for n in order], dtype=float)


def read_ug4_metadata(path: Path) -> list[str]:
    """Return ``# UG4_...`` metadata lines (OuterWall box, terminal flag).

    The repair never moves the endpoints, so box metadata written by
    extend_swc_terminals_to_outerwall.py stays valid and is carried over.
    """
    return [line.strip() for line in path.read_text().splitlines()
            if line.strip().startswith("# UG4_")]


def write_swc(path: Path, points: np.ndarray, header: list[str],
              metadata: list[str] | None = None) -> None:
    lines = list(metadata or []) + [f"# {h}" for h in header]
    for i, p in enumerate(points):
        parent = -1 if i == 0 else i
        node_type = 1 if i == 0 else 3
        lines.append(f"{i + 1} {node_type} {p[0]:.9f} {p[1]:.9f} {p[2]:.9f} "
                     f"0.000000 {parent}")
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def arclength(points: np.ndarray) -> np.ndarray:
    seg = np.linalg.norm(np.diff(points, axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(seg)])


def resample(points: np.ndarray, spacing: float) -> np.ndarray:
    s = arclength(points)
    n = max(2, int(round(s[-1] / spacing)) + 1)
    t = np.linspace(0.0, s[-1], n)
    return np.column_stack([np.interp(t, s, points[:, k]) for k in range(3)])


def circumradii(points: np.ndarray) -> np.ndarray:
    """Three-point circumradius at every interior point (inf if collinear)."""
    a, b, c = points[:-2], points[1:-1], points[2:]
    ab = np.linalg.norm(b - a, axis=1)
    bc = np.linalg.norm(c - b, axis=1)
    ac = np.linalg.norm(c - a, axis=1)
    twice_area = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = ab * bc * ac / (2.0 * twice_area)
    r[twice_area <= 1e-12] = np.inf
    return r


def nonlocal_pairs(points: np.ndarray, exclusion: float):
    """Distance matrix and mask of pairs farther apart than ``exclusion``."""
    s = arclength(points)
    arc = np.abs(s[:, None] - s[None, :])
    dist = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    mask = np.triu(arc > exclusion, k=1)
    return dist, mask


def audit(points: np.ndarray, exclusion: float) -> dict[str, float]:
    """Same quantities as the Lua pipeline geometry preflight."""
    dist, mask = nonlocal_pairs(points, exclusion)
    return {
        "length_um": float(arclength(points)[-1]),
        "min_curvature_radius_um": float(np.min(circumradii(points))),
        "min_nonlocal_distance_um":
            float(np.min(dist[mask])) if mask.any() else math.inf,
    }


def smooth_along_curve(field: np.ndarray, sigma_points: float) -> np.ndarray:
    """Gaussian smoothing of a per-point vector field along the curve."""
    if sigma_points <= 0.0:
        return field
    half = int(math.ceil(3.0 * sigma_points))
    offsets = np.arange(-half, half + 1)
    kernel = np.exp(-0.5 * (offsets / sigma_points) ** 2)
    out = np.zeros_like(field)
    weight = np.zeros(len(field))
    for off, w in zip(offsets, kernel):
        lo, hi = max(0, -off), min(len(field), len(field) - off)
        out[lo:hi] += w * field[lo + off:hi + off]
        weight[lo:hi] += w
    return out / weight[:, None]


# ---------------------------------------------------------------------------
# Repair
# ---------------------------------------------------------------------------

def split_frozen(points: np.ndarray, head: float, tail: float):
    """Split a polyline into frozen head, free middle, and frozen tail."""
    s = arclength(points)
    total = s[-1]
    if head + tail >= total:
        raise ValueError("frozen end lengths exceed the trace length")

    def point_at(t):
        return np.array([np.interp(t, s, points[:, k]) for k in range(3)])

    def piece(t0, t1):
        inner = points[(s > t0) & (s < t1)]
        return np.vstack([point_at(t0), inner, point_at(t1)])

    return piece(0.0, head), piece(head, total - tail), piece(total - tail,
                                                              total)


def repair(points: np.ndarray, args: argparse.Namespace):
    # Frozen ends (for example straight wall-normal tails produced by the
    # OuterWall terminal extender) are kept exactly; only the middle moves.
    head, middle, tail = split_frozen(points, args.frozen_head,
                                      args.frozen_tail)
    head = resample(head, args.working_spacing) if args.frozen_head > 0 \
        else head[:1]
    tail = resample(tail, args.working_spacing) if args.frozen_tail > 0 \
        else tail[-1:]

    def assemble(mid):
        return np.vstack([head[:-1], mid, tail[1:]])

    middle = resample(middle, args.working_spacing)
    reference_middle = middle.copy()
    reference_s = arclength(reference_middle)
    fidelity_sigma = args.min_curvature_radius / args.working_spacing
    work = assemble(middle)
    original = work.copy()
    n = len(work)
    n_head = len(head) - 1

    def free_mask(count):
        mask = np.zeros(count, dtype=bool)
        mask[n_head + 1:count - len(tail)] = True
        return mask
    sigma = args.displacement_smoothing / args.working_spacing
    curvature_stride = max(1, int(round(args.output_spacing
                                        / args.working_spacing)))
    clearance = args.clearance * (1.0 + args.margin)
    min_radius = args.min_curvature_radius * (1.0 + args.margin)

    history = []
    for iteration in range(1, args.max_iterations + 1):
        # (a) Non-local clearance.
        dist, mask = nonlocal_pairs(work, args.exclusion_arclength)
        violating = mask & (dist < clearance)
        push = np.zeros_like(work)
        count = np.zeros(n)
        for i, j in zip(*np.nonzero(violating)):
            d = max(dist[i, j], 1e-9)
            u = (work[i] - work[j]) / d
            deficit = 0.5 * (clearance - d)
            push[i] += deficit * u
            push[j] -= deficit * u
            count[i] += 1
            count[j] += 1
        active = count > 0
        push[active] /= count[active][:, None]
        push = smooth_along_curve(push, sigma)
        push[~free_mask(n)] = 0.0
        work += args.relaxation * push

        # (b) Curvature at the output sampling stride.
        k = curvature_stride
        for start in range(k):
            idx = np.arange(start + k, n - k, k)
            if len(idx) == 0:
                continue
            a, b, c = work[idx - k], work[idx], work[idx + k]
            ab = np.linalg.norm(b - a, axis=1)
            bc = np.linalg.norm(c - b, axis=1)
            ac = np.linalg.norm(c - a, axis=1)
            area2 = np.linalg.norm(np.cross(b - a, c - a), axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(area2 > 1e-12, ab * bc * ac / (2 * area2), np.inf)
            sharp = r < min_radius
            if sharp.any():
                move = np.zeros_like(work)
                mid = 0.5 * (a + c) - b
                move[idx[sharp]] = mid[sharp]
                move = smooth_along_curve(move, 0.5 * k)
                move[~free_mask(n)] = 0.0
                work += args.relaxation * move

        # (c) Weak bending stiffness, then uniform re-parametrization.  The
        # length is free: an expanding coil may lengthen, a rounded corner
        # shortens.  Holding length fixed with pinned endpoints buckled the
        # curve into sub-sample zigzags.
        if args.bending > 0.0:
            lap = np.zeros_like(work)
            lap[1:-1] = 0.5 * (work[:-2] + work[2:]) - work[1:-1]
            lap[~free_mask(n)] = 0.0
            work += args.bending * lap
        middle = resample(work[n_head:n - len(tail) + 1], args.working_spacing)

        # (d) Fidelity: pull every free point toward the original trace point
        # at the same normalized arc length.  Without it, curvature smoothing
        # acts like curve-shortening flow and the whole coil contracts.
        if args.fidelity > 0.0:
            s_mid = arclength(middle)
            t = s_mid / s_mid[-1] * reference_s[-1]
            target = np.column_stack([np.interp(t, reference_s,
                                                reference_middle[:, k])
                                      for k in range(3)])
            # Only the low-frequency part of the pull is applied: the raw
            # manual trace has clicked corners far sharper than the tube, and
            # pulling toward them pointwise would undo the curvature repair.
            pull = smooth_along_curve(target - middle, fidelity_sigma)
            middle[1:-1] += args.fidelity * pull[1:-1]
            middle = resample(middle, args.working_spacing)
        work = assemble(middle)
        n = len(work)

        out = resample(work, args.output_spacing)
        metrics = audit(out, args.exclusion_arclength)
        history.append(metrics)
        if args.verbose and (iteration % 25 == 0 or iteration == 1):
            print(f"iter {iteration:4d}: "
                  f"minRc={metrics['min_curvature_radius_um']:.2f} "
                  f"clearance={metrics['min_nonlocal_distance_um']:.2f} "
                  f"length={metrics['length_um']:.1f}")
        if (metrics["min_curvature_radius_um"] >= args.min_curvature_radius
                and metrics["min_nonlocal_distance_um"] >= args.clearance):
            return out, work, original, iteration, history
    raise RuntimeError(
        f"repair did not converge in {args.max_iterations} iterations: "
        f"{history[-1]}")


def displacement_metrics(work: np.ndarray, original: np.ndarray) -> dict:
    # The two curves are re-parametrized independently, so measure the
    # distance from every original sample to the nearest repaired sample.
    disp = np.min(np.linalg.norm(original[:, None, :] - work[None, :, :],
                                 axis=2), axis=1)
    return {
        "max_displacement_um": float(disp.max()),
        "mean_displacement_um": float(disp.mean()),
        "fraction_moved_over_2um": float(np.mean(disp > 2.0)),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("input", type=Path, help="raw SWC trace")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="repaired SWC in micrometres")
    parser.add_argument("--voxel-size", type=float, default=1.0,
                        help="micrometres per SWC coordinate unit")
    parser.add_argument("--outer-radius", type=float, default=20.0,
                        help="tube outer (Basolateral) radius in um")
    parser.add_argument("--inter-shell", type=float, default=3.0,
                        help="structured near-Inter shell thickness in um")
    parser.add_argument("--clearance", type=float, default=None,
                        help="minimum non-local centreline distance in um "
                             "(default 2*(outer+shell)/0.9)")
    parser.add_argument("--min-curvature-radius", type=float, default=None,
                        help="minimum circumradius in um (default "
                             "(outer + inter shell)/0.82)")
    parser.add_argument("--exclusion-arclength", type=float, default=None,
                        help="arc separation below which pairs are local "
                             "(default max(4*outer, 4*output spacing))")
    parser.add_argument("--working-spacing", type=float, default=1.0)
    parser.add_argument("--output-spacing", type=float, default=5.0)
    parser.add_argument("--displacement-smoothing", type=float, default=6.0,
                        help="Gaussian length (um) for smoothing pushes")
    parser.add_argument("--relaxation", type=float, default=0.5)
    parser.add_argument("--margin", type=float, default=0.03,
                        help="relative over-target used while iterating")
    parser.add_argument("--fidelity", type=float, default=0.0,
                        help="per-iteration low-pass pull toward the original "
                             "trace; off by default because it opposes the "
                             "loop expansion needed for clearance and did not "
                             "converge on 3rd_trace_243_third")
    parser.add_argument("--bending", type=float, default=0.05,
                        help="weight of the global Laplacian smoothing step")
    parser.add_argument("--max-iterations", type=int, default=2000)
    parser.add_argument("--frozen-head", type=float, default=0.0,
                        help="arc length (um) at the root kept unchanged")
    parser.add_argument("--frozen-tail", type=float, default=0.0,
                        help="arc length (um) at the terminal kept unchanged")
    parser.add_argument("--report", type=Path, default=None,
                        help="optional JSON report path")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    if args.clearance is None:
        args.clearance = 2.0 * (args.outer_radius + args.inter_shell) / 0.9
    if args.min_curvature_radius is None:
        # Calibrated on true-scale meshes with the all-hex shell and outward
        # caps (2026-11-23): a 20 um tube + 3 um shell passed r0-r2 at a
        # 28 um minimum bend radius (ratio 0.82) and failed at 26 um (0.88)
        # and 23 um (1.0), where the tube map folds on the inner side.
        args.min_curvature_radius = (args.outer_radius + args.inter_shell) / 0.82
    if args.exclusion_arclength is None:
        args.exclusion_arclength = max(4.0 * args.outer_radius,
                                       4.0 * args.output_spacing)

    raw = read_unbranched_swc(args.input) * args.voxel_size
    raw_out = resample(raw, args.output_spacing)
    before = audit(raw_out, args.exclusion_arclength)
    before["raw_polyline_length_um"] = float(arclength(raw)[-1])

    out, work, original, iterations, _ = repair(raw, args)
    after = audit(out, args.exclusion_arclength)
    disp = displacement_metrics(work, original)

    report = {
        "input": str(args.input),
        "voxel_size_um": args.voxel_size,
        "targets": {
            "clearance_um": args.clearance,
            "min_curvature_radius_um": args.min_curvature_radius,
            "exclusion_arclength_um": args.exclusion_arclength,
            "outer_radius_um": args.outer_radius,
            "inter_shell_um": args.inter_shell,
        },
        "iterations": iterations,
        "before": before,
        "after": after,
        "displacement": disp,
        "length_change_percent":
            100.0 * (after["length_um"] / before["length_um"] - 1.0),
    }

    header = [
        "Repaired centreline from repair_trace_clearance.py",
        f"SOURCE: {args.input}",
        "UNITS: micrometres (run the mesh pipeline with -centerlineScale 1)",
        f"VOXEL_SIZE_UM: {args.voxel_size}",
        f"TARGETS: clearance={args.clearance:.4f} "
        f"minCurvatureRadius={args.min_curvature_radius:.4f} "
        f"exclusionArclength={args.exclusion_arclength:.4f}",
        f"BEFORE: {json.dumps(before)}",
        f"AFTER: {json.dumps(after)}",
        f"DISPLACEMENT: {json.dumps(disp)}",
    ]
    write_swc(args.output, out, header, read_ug4_metadata(args.input))
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
