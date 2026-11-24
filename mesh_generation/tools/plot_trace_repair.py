#!/usr/bin/env python3
"""Plot an original SWC trace against its repaired centreline.

Three orthogonal projections; repaired points coloured by their distance from
the original trace; shaded discs show the tube outer radius.

    python3 mesh_generation/tools/plot_trace_repair.py "Nephron traces/4th_trace_167.swc" \
        repaired.swc -o comparison.png --voxel-size 1 --repaired-voxel-size 1 --outer-radius 20
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from repair_trace_clearance import arclength, read_unbranched_swc, resample  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("original", type=pathlib.Path)
    ap.add_argument("repaired", type=pathlib.Path)
    ap.add_argument("-o", "--output", type=pathlib.Path, required=True)
    ap.add_argument("--voxel-size", type=float, default=1.0,
                    help="um per unit of the original SWC")
    ap.add_argument("--repaired-voxel-size", type=float, default=1.0,
                    help="um per unit of the repaired SWC (1000 for mm)")
    ap.add_argument("--outer-radius", type=float, default=20.0)
    ap.add_argument("--title", default=None)
    a = ap.parse_args()

    raw = read_unbranched_swc(a.original) * a.voxel_size
    rep = read_unbranched_swc(a.repaired) * a.repaired_voxel_size
    dense = resample(raw, 1.0)
    disp = np.min(np.linalg.norm(rep[:, None] - dense[None], axis=2), axis=1)
    fig, axs = plt.subplots(1, 3, figsize=(18, 7))
    for ax, (i, j, lab) in zip(axs, [(0, 1, "x-y"), (0, 2, "x-z"), (1, 2, "y-z")]):
        for p in rep[::max(1, len(rep) // 150)]:
            ax.add_patch(plt.Circle((p[i], p[j]), a.outer_radius,
                                    color="C3", alpha=0.03, lw=0))
        ax.plot(raw[:, i], raw[:, j], "-o", color="0.5", ms=2, lw=1,
                label=f"original ({arclength(raw)[-1]:.0f} um)")
        sc = ax.scatter(rep[:, i], rep[:, j], c=disp, cmap="viridis", s=8,
                        zorder=3, label=f"repaired ({arclength(rep)[-1]:.0f} um)")
        ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_title(lab + " (um)")
    axs[0].legend(fontsize=8, loc="lower left")
    cb = fig.colorbar(sc, ax=axs, shrink=.6)
    cb.set_label("distance from original trace (um)")
    fig.suptitle(a.title or f"{a.original.name} vs {a.repaired.name}")
    fig.savefig(a.output, dpi=110, bbox_inches="tight")
    print(f"max move {disp.max():.1f} um, mean {disp.mean():.1f} um -> {a.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
