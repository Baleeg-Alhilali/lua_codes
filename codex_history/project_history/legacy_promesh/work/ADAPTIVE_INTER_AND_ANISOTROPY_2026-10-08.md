# Adaptive Inter refinement and curvature-safe anisotropy

Date: 2026-10-08

## Objective

The coarse mesh now uses two independent, scale-relative controls:

1. Inter tetrahedra are smaller near Basolateral and become progressively
   larger toward OuterWall.
2. Nephron axial cells use higher anisotropy on nearly straight spline
   intervals, but bends are subdivided toward anisotropy 1.

Only the Neuro Collection plugin and the mesh-generation Lua files were
changed. UG4 core and NavierStokes were not edited.

## Inter size field

After the closed Basolateral/OuterWall PLC is tetrahedralized, the lower
quartile of Basolateral boundary-edge lengths is used as `surfaceScale`.
Distances and target edges are multiples of this measured scale, so the rule
tracks changes in nephron radius and does not depend on a fixed number of mm.

Defaults:

- near distance: `1.5 * surfaceScale`
- far distance: `6.0 * surfaceScale`
- near target edge: `2.0 * surfaceScale`
- far target edge: `6.0 * surfaceScale`
- between the two distances, the target changes linearly
- beyond the far distance, no additional volume constraint is imposed

TetGen retetrahedralization preserves the recovered PLC boundary. Subsequent
UG global refinement refines Lumen, Membrane, Inter, and OuterWall together,
so the child grids remain conforming and form the normal multigrid hierarchy.

The first candidate used edge factors `0.75 -> 4` and produced 460,761 coarse
volumes. It was rejected as excessive. The accepted `2 -> 6` field produced
37,145 volumes for the large-radius validation case.

Measured grading for that accepted mesh:

| distance band | tetrahedra | mean longest edge (mm) |
|---|---:|---:|
| near | 8,966 | 0.0317102 |
| transition | 14,093 | 0.0507019 |
| far | 8,846 | 0.173954 |

## Curvature-safe anisotropy

`-anisotropy 4` is now the straight-path target. Each initial spline interval
is sampled at one-quarter, one-half, and three-quarter position. The largest
deviation from its straight chord is compared with
`bendChordErrorFactor * localRadius`.

If a bend is detected, that interval is recursively split until its curved
length is no larger than the axial target corresponding to
`-bendAnisotropy`. The default bend anisotropy is 1, which makes axial and
circumferential scales comparable around bends. Straight paths retain the
coarser anisotropy-4 station spacing.

For the validated third trace, anisotropy 4 initially requested 46 axial
segments. Curvature adaptation increased only the curved parts, producing 131
segments. A uniform anisotropy-1 grid would refine the entire path.

## Lua controls

```text
-anisotropy 4
-adaptiveAnisotropy 1
-bendAnisotropy 1
-bendChordErrorFactor 0.02
-gradedInter 1
-interNearEdgeFactor 2
-interNearDistanceFactor 1.5
-interFarDistanceFactor 6
-interFarEdgeFactor 6
```

Rollback/testing switches:

- old uniform safe tube spacing: `-anisotropy 1 -adaptiveAnisotropy 0`
- disable the new Inter field: `-gradedInter 0`

## Validation

Validation geometry used radius 15 um, membrane thickness 10 um,
centerline scale 6, and X/Y/Z padding 625/200/400 um.

- baseline, new controls disabled: 27,560 coarse volumes
- accepted adaptive coarse mesh: 37,145 volumes
- projected global refinement levels 1 and 2: zero indeterminate/degenerate
  volumes (level 2 corrected four orientations before serialization)
- Apical/Basolateral refinement radius audit: mean residual
  `6.31654e-11 mm`, maximum `1.07725e-09 mm`
- active Model pipeline smoke test reproduced 37,145 volumes and the same
  measured near/transition/far grading

Files:

- `runs/adaptive_refinement/baseline.ugx`
- `runs/adaptive_refinement/final_defaults.ugx`
- `runs/adaptive_refinement/final_defaults_refined_1.ugx`
- `runs/adaptive_refinement/final_defaults_ref2_refined_2.ugx`
- `runs/adaptive_refinement/final_defaults.log`

The single- and multiple-nephron Lua pipelines under
`/Users/alhilaba/Documents/Promesh_swc/Model/code` and `ProMeshFiles` contain
the tested defaults.
