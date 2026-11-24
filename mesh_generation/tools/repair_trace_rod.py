#!/usr/bin/env python3
"""Length-preserving repair of an SWC centreline (inextensible-rod model).

Purpose
-------
``repair_trace_clearance.py`` enforces the mesh limits (minimum bend radius,
non-local clearance) by local relaxation with a free length.  Curvature
relaxation then acts like curve-shortening flow: the third-trace segment
(``3rd_trace_243_third.swc``) contracted from 349 um to 149 um at the
calibrated limits.  This tool treats the centreline as an inextensible rod:
the trace is discretized into ``N`` segments of equal, *fixed* length and the
rod is deformed in 3-D until it satisfies the limits.  Excess length has to
go somewhere, so the coil expands outward and out of its plane instead of
contracting.

Model
-----
Unknowns are the interior points ``x_1 .. x_{N-1}`` (the ends and optional
frozen head/tail stay fixed, or the ends are softly anchored with
``--free-ends``).  The energy minimized with L-BFGS (PyTorch, double) is

    E = w_len  * sum_i ((|x_{i+1}-x_i| - h)/h)^2           (inextensibility)
      + w_geo  * sum_i relu(kappa_i - 1/R*)^2 * R*^2       (bend radius)
      + w_geo  * sum_ij relu(C* - |x_i-x_j|)^2 / C*^2      (clearance, arc > excl)
      + w_fid  * sum_i |x_i - y_i|^2 / ell^2               (shape fidelity)
      + w_fair * sum_i |x_{i-1}-2x_i+x_{i+1}|^2 / h^2      (fairing)

``kappa_i`` is the three-point circumcurvature evaluated at the working
spacing and at the doubled spacing (the pipeline audits at 5 um), ``R*`` and
``C*`` are the targets with a small margin, and ``y_i`` is the denoised
input at the same arc length (arc length is a material coordinate of an
inextensible rod).  The constraint weight ``w_geo`` is increased by
continuation until the output passes the audit; the targets themselves are
ramped up from the input's values so the coil unfolds gradually.

Target length
-------------
The raw polyline length includes voxel-click noise.  By default the rod
length is the length of the input after a light Gaussian denoise
(``--denoise-sigma``, um), reported next to the raw length; ``--target-length``
overrides it.

Output is in micrometres (run the pipeline with ``-centerlineScale 1
-skipSmoothing 1``); ``# UG4_`` metadata is carried over when the endpoints
are fixed, so the frozen-tail pass on an OuterWall-extended trace keeps its
box metadata valid.  Requires numpy and torch.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from repair_trace_clearance import (arclength, audit,  # noqa: E402
                                    displacement_metrics, read_ug4_metadata,
                                    read_unbranched_swc, resample,
                                    smooth_along_curve, write_swc)

torch.set_default_dtype(torch.float64)


def resample_n(points: np.ndarray, n_seg: int) -> np.ndarray:
    s = arclength(points)
    t = np.linspace(0.0, s[-1], n_seg + 1)
    return np.column_stack([np.interp(t, s, points[:, k]) for k in range(3)])


def denoise(points: np.ndarray, sigma_um: float, spacing: float) -> np.ndarray:
    """Gaussian smoothing along the curve with fixed endpoints."""
    if sigma_um <= 0:
        return points
    p = resample(points, spacing)
    q = smooth_along_curve(p, sigma_um / spacing)
    q[0], q[-1] = p[0], p[-1]
    return q


def circumcurvature(a, b, c):
    ab = torch.linalg.norm(b - a, dim=1)
    bc = torch.linalg.norm(c - b, dim=1)
    ac = torch.linalg.norm(c - a, dim=1)
    cr = torch.linalg.cross(b - a, c - a, dim=1)
    area2 = torch.sqrt((cr * cr).sum(1) + 1e-18)
    return 2.0 * area2 / (ab * bc * ac + 1e-18)


def read_swc_chains(path: Path) -> list[np.ndarray]:
    """Root-to-terminal point chains of an SWC with one or more unbranched
    trees (a merged multi-nephron SWC), in file order of the roots."""
    nodes, order = {}, []
    for line in path.read_text().splitlines():
        payload = line.split("#", 1)[0].split()
        if not payload:
            continue
        if len(payload) != 7:
            raise ValueError(f"SWC line must have 7 columns: {line!r}")
        nid, parent = int(payload[0]), int(payload[6])
        nodes[nid] = (float(payload[2]), float(payload[3]), float(payload[4]),
                      parent)
        order.append(nid)
    children: dict[int, list[int]] = {}
    for n, v in nodes.items():
        if v[3] != -1:
            children.setdefault(v[3], []).append(n)
    if any(len(c) > 1 for c in children.values()):
        raise ValueError("branched SWC traces are not supported")
    chains = []
    for root in [n for n in order if nodes[n][3] == -1]:
        chain = [root]
        while chain[-1] in children:
            chain.append(children[chain[-1]][0])
        chains.append(np.array([nodes[n][:3] for n in chain], dtype=float))
    if sum(len(c) for c in chains) != len(nodes):
        raise ValueError("SWC contains disconnected nodes")
    return chains


def write_swc_chains(path: Path, chains, header, metadata=None) -> None:
    """One multi-root SWC: every chain is its own root-to-terminal tree."""
    lines = list(metadata or []) + [f"# {h}" for h in header]
    nid = 0
    for chain in chains:
        for i, p in enumerate(chain):
            nid += 1
            parent = -1 if i == 0 else nid - 1
            lines.append(f"{nid} {1 if i == 0 else 3} {p[0]:.9f} {p[1]:.9f} "
                         f"{p[2]:.9f} 0.000000 {parent}")
    path.write_text("\n".join(lines) + "\n")


class Rod:
    """One or several inextensible rods solved together.

    Points of all curves are stacked; ``offsets`` delimit the curves.  Length,
    bending, fidelity and fairing act per curve; clearance acts on every pair
    farther apart than the exclusion arc length on the same curve and on
    every pair of points on different curves (no exclusion)."""

    def __init__(self, inits, h_list, fixed_list, excl_steps_list, args):
        self.args = args
        self.offsets = np.cumsum([0] + [len(c) for c in inits])
        init = np.vstack(inits)
        self.x0 = torch.tensor(init)
        self.ref = torch.tensor(init.copy())
        self.h = h_list
        fixed = np.concatenate(fixed_list)
        self.free_idx = torch.nonzero(torch.tensor(~fixed)).flatten()
        pi, pj = [], []
        for c, ex in enumerate(excl_steps_list):
            o, n = self.offsets[c], len(inits[c])
            i, j = np.triu_indices(n, k=ex + 1)
            pi.append(i + o); pj.append(j + o)
            for d in range(c + 1, len(inits)):
                od, nd = self.offsets[d], len(inits[d])
                ii, jj = np.meshgrid(np.arange(n) + o, np.arange(nd) + od,
                                     indexing="ij")
                pi.append(ii.ravel()); pj.append(jj.ravel())
        self.pi = torch.tensor(np.concatenate(pi))
        self.pj = torch.tensor(np.concatenate(pj))
        ends = [self.offsets[c] for c in range(len(inits))] + \
               [self.offsets[c + 1] - 1 for c in range(len(inits))]
        self.ends = torch.tensor(ends)
        self.anchor = torch.tensor(init[ends])

    def curves(self, x):
        return [x[self.offsets[c]:self.offsets[c + 1]]
                for c in range(len(self.offsets) - 1)]

    def full(self, free: torch.Tensor) -> torch.Tensor:
        x = self.x0.clone()
        x[self.free_idx] = free
        return x

    def terms(self, x: torch.Tensor, rstar: float, cstar: float):
        e_len = e_bend = e_fair = 0.0
        kmax = 1.0 / rstar
        for xc, h in zip(self.curves(x), self.h):
            seg = torch.linalg.norm(xc[1:] - xc[:-1], dim=1)
            e_len = e_len + (((seg - h) / h) ** 2).sum()
            k1 = circumcurvature(xc[:-2], xc[1:-1], xc[2:])
            e_bend = e_bend + ((torch.relu(k1 - kmax) * rstar) ** 2).sum()
            if len(xc) > 4:
                k2 = circumcurvature(xc[:-4], xc[2:-2], xc[4:])
                e_bend = e_bend + ((torch.relu(k2 - kmax) * rstar) ** 2).sum()
            lap = xc[:-2] - 2 * xc[1:-1] + xc[2:]
            e_fair = e_fair + ((lap / h) ** 2).sum()
        d = torch.linalg.norm(x[self.pi] - x[self.pj], dim=1)
        e_clr = ((torch.relu(cstar - d) / cstar) ** 2).sum()
        e_fid = (((x - self.ref) / self.args.fidelity_length) ** 2).sum()
        e_end = (((x[self.ends] - self.anchor)
                  / self.args.fidelity_length) ** 2).sum()
        return e_len, e_bend, e_clr, e_fid, e_fair, e_end

    def solve(self, x: torch.Tensor, rstar, cstar, w_geo, iters):
        a = self.args
        free = x[self.free_idx].clone().requires_grad_(True)
        opt = torch.optim.LBFGS([free], lr=1.0, max_iter=iters,
                                history_size=50, tolerance_grad=1e-10,
                                tolerance_change=1e-14,
                                line_search_fn="strong_wolfe")

        def closure():
            opt.zero_grad()
            e_len, e_bend, e_clr, e_fid, e_fair, e_end = self.terms(
                self.full(free), rstar, cstar)
            e = (a.w_length * e_len + w_geo * (e_bend + e_clr)
                 + a.w_fidelity * e_fid + a.w_fair * e_fair
                 + a.w_end * e_end)
            e.backward()
            return e

        opt.step(closure)
        return self.full(free.detach())


def cross_distance(curves) -> float:
    """Minimum distance between points of different curves (inf for one)."""
    best = math.inf
    for c in range(len(curves)):
        for d in range(c + 1, len(curves)):
            dist = np.linalg.norm(curves[c][:, None] - curves[d][None], axis=2)
            best = min(best, float(dist.min()))
    return best


def prepare_curve(raw: np.ndarray, args):
    reference_dense = denoise(raw, args.denoise_sigma, 0.5)
    raw_len = float(arclength(raw)[-1])
    den_len = float(arclength(reference_dense)[-1])
    target = args.target_length or den_len * args.length_factor
    n_seg = max(8, int(round(target / args.working_spacing)))
    h = target / n_seg
    # Initial rod: the denoised trace re-parametrized to N equal pieces.
    init = resample_n(reference_dense, n_seg)
    fixed = np.zeros(n_seg + 1, dtype=bool)
    if not args.free_ends:
        fixed[0] = fixed[-1] = True
    head_n = int(round(args.frozen_head / h))
    tail_n = int(round(args.frozen_tail / h))
    if head_n or tail_n:
        if args.free_ends:
            raise ValueError("frozen ends require fixed endpoints")
        # Keep the frozen pieces exactly (rigid tails of the input).
        s_raw = arclength(raw)
        for idx in range(head_n + 1):
            t = idx * h
            init[idx] = [np.interp(t, s_raw, raw[:, k]) for k in range(3)]
        for idx in range(tail_n + 1):
            t = s_raw[-1] - idx * h
            init[n_seg - idx] = [np.interp(t, s_raw, raw[:, k])
                                 for k in range(3)]
        fixed[:head_n + 1] = True
        fixed[n_seg - tail_n:] = True
    info = {"raw_length_um": raw_len, "denoised_length_um": den_len,
            "target_length_um": target, "segments": n_seg,
            "segment_length_um": h}
    return init, h, fixed, int(math.ceil(args.exclusion_arclength / h)), info


def repair(raws, args):
    """Repair one or several curves together; returns per-curve outputs."""
    prepared = [prepare_curve(r, args) for r in raws]
    inits = [p[0] for p in prepared]
    rod = Rod(inits, [p[1] for p in prepared], [p[2] for p in prepared],
              [p[3] for p in prepared], args)

    rmin = args.min_curvature_radius * (1 + args.margin)
    cmin = args.clearance * (1 + args.margin)
    starts = [audit(resample(c, args.output_spacing), args.exclusion_arclength)
              for c in inits]
    r0 = max(1.0, min([rmin] + [s["min_curvature_radius_um"] for s in starts]))
    c0 = max(1.0, min([cmin, cross_distance(inits)]
                      + [s["min_nonlocal_distance_um"] for s in starts]))

    x = torch.tensor(np.vstack(inits))
    history = []
    t0 = time.time()
    # Ramp the targets (homotopy) so the coil unfolds gradually.
    for step in range(1, args.ramp_steps + 1):
        f = step / args.ramp_steps
        x = rod.solve(x, r0 + f * (rmin - r0), c0 + f * (cmin - c0),
                      args.w_geometry, args.lbfgs_iters)
    # Continuation on the constraint weight at the final targets.
    w = args.w_geometry
    for stage in range(args.max_stages):
        x = rod.solve(x, rmin, cmin, w, args.lbfgs_iters)
        rods = [c.numpy() for c in rod.curves(x)]
        outs = [resample(c, args.output_spacing) for c in rods]
        ms = [audit(o, args.exclusion_arclength) for o in outs]
        m = {"min_curvature_radius_um":
                 min(q["min_curvature_radius_um"] for q in ms),
             "min_nonlocal_distance_um":
                 min(q["min_nonlocal_distance_um"] for q in ms),
             "min_cross_curve_distance_um": cross_distance(outs),
             "length_um": [q["length_um"] for q in ms],
             "segment_length_rel_dev_max": max(
                 float(np.max(np.abs(np.linalg.norm(np.diff(c, axis=0), axis=1)
                                     / p[1] - 1)))
                 for c, p in zip(rods, prepared)),
             "weight": w}
        history.append(m)
        if args.verbose:
            print(f"stage {stage}: w={w:.3g} minR={m['min_curvature_radius_um']:.2f} "
                  f"clr={m['min_nonlocal_distance_um']:.2f} "
                  f"cross={m['min_cross_curve_distance_um']:.2f} "
                  f"L={[round(v, 1) for v in m['length_um']]} "
                  f"segdev={m['segment_length_rel_dev_max']:.2e} "
                  f"({time.time() - t0:.0f}s)", flush=True)
        if (m["min_curvature_radius_um"] >= args.min_curvature_radius
                and m["min_nonlocal_distance_um"] >= args.clearance
                and m["min_cross_curve_distance_um"] >= args.clearance):
            break
        w *= 4.0
    else:
        raise RuntimeError(f"rod repair did not meet the limits: {history[-1]}")
    infos = [p[4] for p in prepared]
    for i in infos:
        i["seconds"] = time.time() - t0
    return outs, rods, infos, history


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("input", type=Path, nargs="+",
                   help="one SWC, several traces repaired jointly (clearance "
                        "also enforced between traces), or one merged "
                        "multi-root SWC whose trees are repaired jointly")
    p.add_argument("-o", "--output", type=Path, nargs="+", required=True,
                   help="one output SWC per input")
    p.add_argument("--voxel-size", type=float, default=1.0,
                   help="micrometres per SWC unit (1000 for mm SWCs)")
    p.add_argument("--outer-radius", type=float, default=20.0)
    p.add_argument("--inter-shell", type=float, default=3.0)
    p.add_argument("--clearance", type=float, default=55.0,
                   help="minimum non-local distance (um), also between "
                        "traces; 55 is the value used for the passing "
                        "true-scale meshes")
    p.add_argument("--min-curvature-radius", type=float, default=None,
                   help="default (outer + shell)/0.72 = 32 um; the rod holds "
                        "long arcs at the limit, and 28 um (ratio 0.82, enough "
                        "for the old free-length repair) left one degenerate "
                        "cap at r2 on 3rd_trace_243_third (2026-11-24)")
    p.add_argument("--exclusion-arclength", type=float, default=None,
                   help="default max(4*outer, 4*output spacing)")
    p.add_argument("--working-spacing", type=float, default=2.5)
    p.add_argument("--output-spacing", type=float, default=5.0)
    p.add_argument("--denoise-sigma", type=float, default=2.0,
                   help="Gaussian sigma (um) defining the reference shape "
                        "and default target length")
    p.add_argument("--target-length", type=float, default=None,
                   help="single-trace runs only")
    p.add_argument("--length-factor", type=float, default=1.0,
                   help="target = factor * denoised length")
    p.add_argument("--free-ends", action="store_true",
                   help="softly anchor instead of fixing the endpoints")
    p.add_argument("--frozen-head", type=float, default=0.0)
    p.add_argument("--frozen-tail", type=float, default=0.0)
    p.add_argument("--margin", type=float, default=0.03)
    p.add_argument("--w-length", type=float, default=1e4)
    p.add_argument("--w-geometry", type=float, default=10.0)
    p.add_argument("--w-fidelity", type=float, default=1e-3)
    p.add_argument("--fidelity-length", type=float, default=10.0)
    p.add_argument("--w-fair", type=float, default=1e-3)
    p.add_argument("--w-end", type=float, default=10.0)
    p.add_argument("--ramp-steps", type=int, default=8)
    p.add_argument("--max-stages", type=int, default=12)
    p.add_argument("--lbfgs-iters", type=int, default=400)
    p.add_argument("--report", type=Path, default=None)
    p.add_argument("--verbose", action="store_true")
    a = p.parse_args(argv)
    if len(a.output) != len(a.input):
        p.error("give one output per input (one output for a merged SWC)")
    if a.target_length is not None and len(a.input) > 1:
        p.error("--target-length applies to single-trace runs")
    if a.min_curvature_radius is None:
        a.min_curvature_radius = (a.outer_radius + a.inter_shell) / 0.72
    if a.exclusion_arclength is None:
        a.exclusion_arclength = max(4 * a.outer_radius, 4 * a.output_spacing)

    merged = len(a.input) == 1 and len(read_swc_chains(a.input[0])) > 1
    if merged:
        # One merged multi-root SWC (e.g. the OuterWall-extended output of
        # swc_to_neurite_projector_multiple_nephrons.lua): repair its trees
        # jointly and write them back as one multi-root SWC.
        raws = [c * a.voxel_size for c in read_swc_chains(a.input[0])]
        sources = [Path(f"{a.input[0]}#tree{k + 1}") for k in range(len(raws))]
        outputs = [a.output[0]] * len(raws)
    else:
        raws = [read_unbranched_swc(f) * a.voxel_size for f in a.input]
        sources, outputs = a.input, a.output
    outs, rods, infos, history = repair(raws, a)
    reports = []
    headers = []
    for f, o, raw, out, rod, info in zip(sources, outputs, raws, outs, rods,
                                         infos):
        before = audit(resample(raw, a.output_spacing), a.exclusion_arclength)
        after = audit(out, a.exclusion_arclength)
        disp = displacement_metrics(rod, resample(raw, 1.0))
        report = {
            "input": str(f), "output": str(o),
            "method": "inextensible rod (repair_trace_rod.py)",
            "joint_with": [str(g) for g in sources if g != f],
            "targets": {"clearance_um": a.clearance,
                        "min_curvature_radius_um": a.min_curvature_radius,
                        "exclusion_arclength_um": a.exclusion_arclength},
            "before": before, "after": after, "displacement": disp,
            "length_retained_vs_raw": after["length_um"] / info["raw_length_um"],
            **info,
        }
        if len(sources) > 1:
            report["min_cross_trace_distance_um"] = \
                history[-1]["min_cross_curve_distance_um"]
            report["before_min_cross_trace_distance_um"] = cross_distance(
                [resample(r, a.output_spacing) for r in raws])
        reports.append(report)
        meta = [] if a.free_ends else read_ug4_metadata(a.input[0] if merged else f)
        header = [
            "Length-preserving repair from repair_trace_rod.py (inextensible rod)",
            f"SOURCE: {f}",
            *([f"JOINT_WITH: {' '.join(report['joint_with'])}"]
              if len(sources) > 1 else []),
            "UNITS: micrometres (run the mesh pipeline with -centerlineScale 1 -skipSmoothing 1)",
            f"VOXEL_SIZE_UM: {a.voxel_size}",
            f"TARGETS: clearance={a.clearance:.4f} minCurvatureRadius="
            f"{a.min_curvature_radius:.4f} exclusionArclength={a.exclusion_arclength:.4f}",
            f"LENGTH_UM: raw={info['raw_length_um']:.2f} denoised="
            f"{info['denoised_length_um']:.2f} repaired={after['length_um']:.2f}",
            f"BEFORE: {json.dumps(before)}",
            f"AFTER: {json.dumps(after)}",
            f"DISPLACEMENT: {json.dumps(disp)}",
        ]
        headers.append((header, meta))
        if not merged:
            write_swc(o, out, header, meta)
    if merged:
        header = ["Joint length-preserving repair of a merged multi-root SWC "
                  "(repair_trace_rod.py)",
                  "UNITS: micrometres"] + [
            f"TREE {k + 1}: {h[0][6]} {h[0][8]}" for k, h in enumerate(headers)]
        write_swc_chains(a.output[0], outs, header, headers[0][1])
    full = reports[0] if len(reports) == 1 else {"traces": reports,
                                                  "stages": history}
    if a.report:
        a.report.write_text(json.dumps(full, indent=2) + "\n")
    for r in reports:
        print(json.dumps({k: r.get(k) for k in (
            "input", "before", "after", "displacement", "length_retained_vs_raw",
            "min_cross_trace_distance_um", "seconds")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
