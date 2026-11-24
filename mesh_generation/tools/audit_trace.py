#!/usr/bin/env python3
"""Audit an SWC centreline against the tube geometry before meshing.

Reports, at the given voxel size, the same quantities as the Lua mesh
pipeline's geometry preflight (computed on a resampled polyline):

* length (raw polyline and resampled),
* minimum three-point curvature radius and how many points are sharper than
  the limit,
* minimum non-local self-distance (pairs farther apart than the exclusion arc
  length) and the closest pairs with their positions,
* bounding box.

Example (1 um voxels, 10 um lumen + 10 um membrane, 3 um Inter shell):

    python3 mesh_generation/tools/audit_trace.py "Nephron traces/4th_trace_167.swc" \
        --voxel-size 1 --outer-radius 20 --inter-shell 3

Limits printed for comparison: preflight curvature outer/0.9 and clearance
2*outer/0.9, the shell clearance 2*(outer+shell)/0.9, and the calibrated mesh
limit (outer+shell)/0.82 for the minimum bend radius (2026-11-23).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from repair_trace_clearance import (arclength, circumradii,  # noqa: E402
                                    read_unbranched_swc, resample)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("swc", type=pathlib.Path)
    ap.add_argument("--voxel-size", type=float, default=1.0,
                    help="micrometres per SWC unit (use 1000 for mm SWCs)")
    ap.add_argument("--outer-radius", type=float, default=20.0)
    ap.add_argument("--inter-shell", type=float, default=3.0)
    ap.add_argument("--spacing", type=float, default=5.0,
                    help="resampling spacing in um (pipeline uses 5)")
    ap.add_argument("--exclusion-arclength", type=float, default=None,
                    help="default max(4*outer, 4*spacing)")
    ap.add_argument("--closest", type=int, default=5)
    ap.add_argument("--json", type=pathlib.Path, default=None)
    a = ap.parse_args()

    raw = read_unbranched_swc(a.swc) * a.voxel_size
    pts = resample(raw, a.spacing)
    s = arclength(pts)
    excl = a.exclusion_arclength or max(4 * a.outer_radius, 4 * a.spacing)
    rc = circumradii(pts)
    arc = np.abs(s[:, None] - s[None, :])
    dist = np.linalg.norm(pts[:, None] - pts[None], axis=2)
    mask = np.triu(arc > excl, 1)
    i, j = np.nonzero(mask)
    order = np.argsort(dist[i, j])
    limits = {
        "preflight_min_curvature_um": a.outer_radius / 0.9,
        "preflight_clearance_um": 2 * a.outer_radius / 0.9,
        "shell_clearance_um": 2 * (a.outer_radius + a.inter_shell) / 0.9,
        "mesh_min_bend_um": (a.outer_radius + a.inter_shell) / 0.82,
    }
    out = {
        "swc": str(a.swc), "voxel_size_um": a.voxel_size,
        "raw_length_um": float(arclength(raw)[-1]),
        "resampled_length_um": float(s[-1]),
        "min_curvature_radius_um": float(rc.min()),
        "points_sharper_than_mesh_limit": int((rc < limits["mesh_min_bend_um"]).sum()),
        "min_nonlocal_distance_um": float(dist[mask].min()) if mask.any() else None,
        "pairs_closer_than_shell_clearance": int((dist[mask] < limits["shell_clearance_um"]).sum()),
        "bbox_min_um": raw.min(0).tolist(), "bbox_max_um": raw.max(0).tolist(),
        "limits": limits,
        "closest_pairs": [
            {"distance_um": float(dist[i[k], j[k]]),
             "arc_separation_um": float(arc[i[k], j[k]]),
             "s_a_um": float(s[i[k]]), "s_b_um": float(s[j[k]]),
             "a_um": pts[i[k]].round(2).tolist(), "b_um": pts[j[k]].round(2).tolist()}
            for k in order[:a.closest]],
    }
    print(json.dumps(out, indent=2))
    if a.json:
        a.json.write_text(json.dumps(out, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
