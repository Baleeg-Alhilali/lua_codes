#!/bin/sh
# Several traces in one shared Inter domain: joint length-preserving repair ->
# OuterWall extension -> joint frozen-tail repair -> multi-nephron mesh.
#
# usage: rod_repair_joint_workflow.sh <out_dir> <prefix> <raw1.swc> <raw2.swc> [...]
#   env: UG4 (default /Users/alhilaba/UG4_promesh_ogrid), PADDING (um, 100),
#        NUMREFS (generator levels, 2; never above 2 on the 36 GB workstation),
#        CHECK=1 to run check_refinement_quality.lua (mode 2) to NUMREFS,
#        ROD_OPTS options for repair_trace_rod.py (default: a 35 um bend
#        target; the joint 3rd+4th repair at 32 um inverted linear caps at r1
#        on the fourth trace, 35 um passed r1 with hex 0.971, 2026-11-24).
# The raw traces end close together at the image edge, so the coil repair uses
# --free-ends (soft end anchors) to separate the terminals before extension.
set -u
OUT=$1; P=$2; shift 2
UG4=${UG4:-/Users/alhilaba/UG4_promesh_ogrid}
PADDING=${PADDING:-100}
NUMREFS=${NUMREFS:-2}
ROD_OPTS=${ROD_OPTS:---min-curvature-radius 35}
HERE=$(cd "$(dirname "$0")" && pwd)
MG=$(cd "$HERE/.." && pwd)
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
MESH="-skipSmoothing 1 -centerlineScale 1 -radius 10 -membraneThickness 10 -ogrid 12 -BoxRefs 3 -interMesh hex -prismInterThickness 3"

COILS=""; SWCARGS=""; k=0
for raw in "$@"; do
    k=$((k + 1)); COILS="$COILS $OUT/${P}_coil_$k.swc"
done
# shellcheck disable=SC2086
python3 "$HERE/repair_trace_rod.py" "$@" -o $COILS --free-ends \
    --report "$OUT/${P}_coil.json" $ROD_OPTS > "$OUT/${P}_coil.log" 2>&1 \
    || { echo "joint coil repair failed"; tail -3 "$OUT/${P}_coil.log"; exit 1; }
k=0
for c in $COILS; do k=$((k + 1)); SWCARGS="$SWCARGS -swc$k $c"; done

# Pass 1 only writes <out>_outerwall_extended.swc (merged, multi-root, mm),
# before TetGen; a TetGen failure in this throwaway pass does not matter.
# shellcheck disable=SC2086
(cd "$MG" && "$UG4/bin/ugshell" -ex swc_to_neurite_projector_multiple_nephrons.lua -- \
    $SWCARGS $MESH -padding "$PADDING" -numRefs 0 -validate 0 \
    -out "$OUT/${P}_p1.ugx") > "$OUT/${P}_p1.log" 2>&1
[ -f "$OUT/${P}_p1_outerwall_extended.swc" ] || { echo "no extended trace"; grep -E "Assert|ERROR|failed" "$OUT/${P}_p1.log" | head; exit 1; }

# shellcheck disable=SC2086
python3 "$HERE/repair_trace_rod.py" "$OUT/${P}_p1_outerwall_extended.swc" --voxel-size 1000 \
    --frozen-head 40 --frozen-tail 40 -o "$OUT/${P}_ext.swc" \
    --report "$OUT/${P}_ext.json" $ROD_OPTS > "$OUT/${P}_ext.log" 2>&1 \
    || { echo "joint extended repair failed"; tail -3 "$OUT/${P}_ext.log"; exit 1; }

# shellcheck disable=SC2086
(cd "$MG" && "$UG4/bin/ugshell" -ex swc_to_neurite_projector_multiple_nephrons.lua -- \
    -mergedSwc "$OUT/${P}_ext.swc" $MESH -padding "$PADDING" -numRefs "$NUMREFS" \
    -out "$OUT/${P}.ugx") > "$OUT/${P}.log" 2>&1
GEN=$?
echo "$P generation exit=$GEN $(grep -E 'VALIDATION:' "$OUT/${P}.log" | tr '\n' ' ')"
grep -E "EXCEPTION|degenerate volume|Assert|  0: " "$OUT/${P}.log" | head -4
if [ "$GEN" = 0 ] && [ "${CHECK:-0}" = 1 ]; then
    (cd "$MG" && "$UG4/bin/ugshell" -ex check_refinement_quality.lua -grid "$OUT/${P}.ugx" \
        -refs "$NUMREFS" -mode 2) > "$OUT/${P}_check.log" 2>&1
    grep -E "Level 0|Local refinement|tets with|worst roid|EXCEPTION" "$OUT/${P}_check.log" | cut -c1-330
fi
exit $GEN
