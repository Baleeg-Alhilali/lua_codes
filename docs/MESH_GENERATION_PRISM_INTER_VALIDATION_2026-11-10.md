# Straight-nephron prismatic Inter validation — 2026-11-10

## Result classification

The selectable straight-control prismatic Inter is a **geometry pass through
r2** and removes the refined hanging-interface pressure/velocity artifact. The
tested monolithic flow is an **algebraically converged diagnostic**, not a
fully physically validated or 1%-mesh-converged production model. ParaView was
therefore not auto-opened.

## Provenance and scope

- Serial run on one MPI rank with
  `/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell`.
- Active plugin binary SHA-256:
  `2dddb35b8855c2d17044f50f5d5b1ef6d9d8d18928eb8493d08861cc64920f53`.
- `neurites_from_swc.cpp` SHA-256:
  `cb743f10382868baf8c4ced21b860725e5df9ff8c451b8344e781625c67ca322`.
- `neurites_from_swc.h` SHA-256:
  `c6c54ace02f84a10f34294bb30a51e7f7d6fc9bb24ce6b09ec2abd8057c44e10`.
- `neuro_collection_plugin.cpp` SHA-256:
  `cbe0f243981c9d08885b82e69080b12fefb2e3dc94a504321745bd55703c31ce`.
- The active `neuro_collection` directory is not a Git worktree. Its complete
  modified source is preserved under
  `codex_history/checkpoints/2026-11-10_prism_inter_verified/neuro_collection`.
- Scope is deliberately restricted to one strictly straight nephron. Curved
  and multiple-nephron generation retain the default `-interMesh tet` path.

The original historical mesh was not overwritten:
`runs/straight_nephron_no_bend.ugx`, SHA-256
`5f79d67ad7dcfd936c0e7fb0ba12c6a231c7288b889b35e4c73dd5e66ac63894`.

The unchanged November source initially crashed while restoring structured
volumes after TetGen had erased their vertices. The repaired selectable `tet`
path now tetrahedralizes an isolated copy and retains the established
constrained Basolateral representation. Its new-file r0--r2 control passes the
same structural validator with zero inverted/degenerate volumes:
`runs/straight_nephron_current_tet_isolated{,_refined_1,_refined_2}.ugx`.
The respective SHA-256 values are `8d95590948fe9fec53fcf8ecc60593a99b65434500cdd981407ae37558543be7`,
`244721ac8b2f90ff632a0284e75fa457293cfa672c3621d116c3171e255128b5`,
and `886281fb2c940345a3c5939191a7bafedd278def7dd6ac3f09b07cbfd8dff0f7`.
This makes the default selector geometrically usable but does not cure its
refined flow oscillation; that is the reason for the prism option.

## Reproduction command

Run from `mesh_generation`:

```text
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex swc_to_neurite_projector_nephron_final_pipeline.lua -- \
  -swc "../Nephron traces/straight_nephron_outerwall_cell_centered.swc" \
  -radius 15 -membraneThickness 2 -centerlineScale 1 \
  -resampleSpacing 5 -BoxRefs 4 -useTraceEndpoints 1 \
  -padding 300 -terminalPaddingFactor 0 \
  -interMesh prism -numRefs 2 \
  -out ../runs/straight_nephron_prism_inter.ugx
```

The generator uses 24 axial stations, 10 Basolateral angular vertices, and a
graded square-minus-circle cross-section containing 144 triangles. Sweeping
the section creates 3312 Inter prisms. Basolateral consists of ordinary shared
Membrane/Inter quadrilaterals; there are no hanging constraints or prism/tet
transition cells.

## Geometry validation

| Level | Vertices | Hexes | Prisms | Basolateral quads | Minimum prism scaled Jacobian |
|---|---:|---:|---:|---:|---:|
| r0 | 2,760 | 690 | 3,542 | 230 | 0.213402 |
| r1 | 20,821 | 5,520 | 28,336 | 920 | 0.213402 |
| r2 | 161,913 | 44,160 | 226,688 | 3,680 | 0.213402 |

All three levels have zero tetrahedra, pyramids, inverted prisms, degenerate
prisms, unmatched faces, duplicate face geometries, non-manifold faces, and
constrained objects. All interface adjacency contracts pass. The r1/r2 radius
audit reports the same 330 tiny off-radius inherited points, with maximum
residual `1.23438e-9 mm`; no volume orientation failure results.

Mesh SHA-256 values:

- r0: `d9535ffd3243980cde5b3e8865223c2a9207cd4f84636450008af6933ea6a69e`
- r1: `a8ccaf02ccc9f2e3911f127a13cb20c9354c8bae4b2c13db9798957ab1184502`
- r2: `dc848d3dbc9f62878a72045e4ed4844b06ac50a9848d256aca7e107b9fb07775`

## Flow diagnostic

The historical and prism meshes used the same pressure-trace/FE-Darcy script,
with only the mesh stem changed. Parameters were `Kmem=1e-19 m2`,
`Kinter=1e-17 m2`, `Kbjs=1e-19 m2`, `bjsAlpha=0.2`, OuterWall pressure
`613.2812 Pa`, outlet pressure `2666.4474 Pa`, and Darcy algebraic pressure
scale `1e4`. The refined solver is BiCGStab with GMG/Uzawa and RAP; no global
SuperLU solve was used.

| Mesh | Level | Krylov steps by cycle | Lumen balance / inlet | Pressure drop (Pa) |
|---|---|---:|---:|---:|
| historical tet | r0 | 1 | about 16.8% | 29.116 |
| historical tet | r1 | 239 + 183 | 3.955% | 53.194 |
| historical tet | r2 | 307 + 149 | 0.796% | 53.774 |
| prism Inter | r0 | 1 | about 16.8% | 29.117 |
| prism Inter | r1 | 232 + 187 | 3.592% | 53.085 |
| prism Inter | r2 | 281 + 165 | 0.373% | 53.627 |

At prism r2 the final relative true defects were `2.12e-10` (u), `3.13e-10`
(v), `1.18e-9` (w), `9.19e-15` (Stokes pressure), and `1.13e-13` (scaled
Darcy pressure). The pressure-trace jump was `5.24e-10`. Inlet, outlet, and
Apical flows were `-1.051295e-3`, `1.054010e-3`, and `1.204776e-6 mm3/s`.

The r1-to-r2 prism changes are 1.01% in pressure drop, 2.59% in inlet flow,
1.22% in outlet flow, and 0.98% in mean radial membrane velocity. Thus the
three-level sequence shows convergence but does not meet the project's 1%
criterion for all required observables. The total Krylov work is also high
(419 steps at r1 and 446 at r2), so refinement-independent solver behavior has
not been established.

## Oscillation and cell-centered velocity audit

Nodal Basolateral pressure was analyzed separately from the element-centered
Darcy velocity. No cell-to-point conversion was used for the velocity audit.

| Mesh | Level | Circumferential pressure residual RMS (Pa) | Alternating pressure mode RMS (Pa) | Cell tangential/speed RMS |
|---|---|---:|---:|---:|
| historical tet | r0 | 3.854 | 0.769 | 0.0199% |
| historical tet | r1 | 1,462 | 1,216 | 43.10% |
| historical tet | r2 | 868.9 | 548.3 | 46.28% |
| prism Inter | r0 | 1.609 | 0.0710 | 0.00524% |
| prism Inter | r1 | 1.248 | 0.472 | 0.00861% |
| prism Inter | r2 | 0.807 | 0.118 | 0.00743% |

The historical constrained interface therefore creates a severe refined
circumferential mode. The conforming prism topology removes it in both nodal
pressure and direct cell-vector diagnostics.

## Outputs and remaining limits

- Meshes: `runs/straight_nephron_prism_inter{,_refined_1,_refined_2}.ugx`.
- Short r2 datasets:
  `Results_final/prism_inter_r2_{lumen,membrane,inter}_t0000.vtu`.
- Analysis tools:
  `code/flow/analyze_basolateral_pressure_oscillation.py` and
  `code/flow/analyze_cellwise_darcy_velocity.py`.
- No r3 run was attempted.
- This result validates the mesh fix, not the physical pressure-trace model.
  The underlying interface-law concerns recorded in `AGENTS.md` remain, the
  r0 conservation error is large, the r1-to-r2 flow change exceeds 1%, and
  GMG iteration behavior is not yet an acceptable scalable production result.
  These are the exact reasons ParaView auto-open remains gated.
