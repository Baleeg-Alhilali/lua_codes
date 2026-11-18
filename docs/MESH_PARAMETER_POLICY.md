# Nephron mesh parameter and tolerance policy

The active mesh pipeline uses physical input lengths in micrometres and writes
UGX coordinates in millimetres. Geometric comparisons must therefore never use
an absolute coordinate epsilon that is meaningful for only one trace or tube
radius.

## Scale hierarchy

Each comparison now uses the smallest relevant local scale:

- SWC smoothing/deduplication: requested resampling distance and trace extent.
- Endpoint/OuterWall fitting: structured box-cell size and box diagonal.
- Terminal snapping: terminal radius on the coarse mesh and local terminal-edge
  length after refinement.
- Duplicate cleanup: minimum positive outer tube radius.
- Straight-tube detection: mean Basolateral radius.
- Angular comparisons: circumferential sector width `2*pi/ogrid`.
- Inter containment: the bounding scale of the individual recorded tube cell.

Every length tolerance is the maximum of a local relative tolerance and
`128 * machine_epsilon * global_scale`.

## Automatic Inter grading

The former fixed seven radial fractions were removed. The number of radial
layers is now computed from:

`layersPerSurfaceEdge * radialThickness / meanBasolateralEdgeLength`

and clamped to `interMinLayers..interMaxLayers`. Layer positions follow an
exponential distribution controlled by `interRadialGrowth`; zero growth gives
uniform spacing.

## Exposed dimensionless controls

| Option | Default | Safe range |
|---|---:|---:|
| `-coordinateRelTol` | `1e-10` | `(0, 1e-2]` |
| `-containmentRelTol` | `1e-6` | `(0, 1e-2]` |
| `-parameterTol` | `1e-7` | `(0, 1e-2]` |
| `-mergeRadiusFactor` | `1e-4` | `(0, 0.1]` |
| `-terminalSnapRadiusFactor` | `2` | `(0, 10]` |
| `-straightnessRadiusFactor` | `1e-7` | `(0, 0.1]` |
| `-angularSpacingToleranceFactor` | `1e-5` | `(0, 0.1]` |
| `-directionTolerance` | `1e-8` | positive |
| `-interMinLayers` | `4` | `2..64` |
| `-interMaxLayers` | `12` | `min..64` |
| `-interLayersPerSurfaceEdge` | `3` | `(0, 100]` |
| `-interRadialGrowth` | `2` | `[0, 20]` |

Normal mesh controls such as radius, membrane thickness, padding, O-grid count,
resampling distance, and Inter thickness remain physical inputs. Geometry
preflight rejects a requested radius that violates local curvature or
non-local centerline clearance.

## Regression coverage

The generalized policy was validated with:

1. Straight tube: radius 15 um, membrane 2 um, padding 300 um, O-grid 10,
   BoxRefs 4, through r2.
2. Curved nephron: radius 15 um, membrane 2 um, centerline scale 24, padding
   200 um, through r1.
3. Scaled straight tube: radius 8 um, membrane 1 um, centerline scale 1.5,
   padding 120 um, O-grid 12, BoxRefs 3, through r1.

All generated levels had zero inverted, degenerate, unmatched, duplicate, and
non-manifold elements/faces and a conforming Basolateral Membrane--Inter
interface.
