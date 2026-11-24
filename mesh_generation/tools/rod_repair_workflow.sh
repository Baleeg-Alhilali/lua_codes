#!/bin/sh
# Length-preserving repair -> OuterWall extension -> frozen-tail repair -> mesh.
#
# usage: rod_repair_workflow.sh <raw.swc> <out_dir> <prefix> [rod options...]
#   env: UG4 (default /Users/alhilaba/UG4_promesh_ogrid), PADDING (um, 100),
#        NUMREFS (generator levels, 2; never above 2 on the 36 GB workstation),
#        CHECK=1 to run check_refinement_quality.lua (mode 2) to NUMREFS.
# Mesh settings: 10 um lumen + 10 um membrane, ogrid 12, BoxRefs 3, hex Inter,
# 3 um structured shell (the true-scale configuration of 2026-11-23/24).
set -u
RAW=$1; OUT=$2; P=$3; shift 3
UG4=${UG4:-/Users/alhilaba/UG4_promesh_ogrid}
PADDING=${PADDING:-100}
NUMREFS=${NUMREFS:-2}
HERE=$(cd "$(dirname "$0")" && pwd)
MG=$(cd "$HERE/.." && pwd)
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
MESH="-radius 10 -membraneThickness 10 -ogrid 12 -BoxRefs 3 -interMesh hex -prismInterThickness 3"

python3 "$HERE/repair_trace_rod.py" "$RAW" -o "$OUT/${P}_coil.swc" \
    --report "$OUT/${P}_coil.json" "$@" > "$OUT/${P}_coil.log" 2>&1 \
    || { echo "coil repair failed"; tail -3 "$OUT/${P}_coil.log"; exit 1; }

# Pass 1 only writes <name>_outerwall_extended.swc (written before TetGen;
# a TetGen failure in this throwaway pass does not matter).
(cd "$MG" && "$UG4/bin/ugshell" -ex swc_to_neurite_projector_nephron_final_pipeline.lua -- \
    -swc "$OUT/${P}_coil.swc" -skipSmoothing 1 -centerlineScale 1 -padding "$PADDING" \
    $MESH -numRefs 0 -validate 0 -out "$OUT/${P}_p1.ugx") > "$OUT/${P}_p1.log" 2>&1
[ -f "$OUT/${P}_p1_outerwall_extended.swc" ] || { echo "no extended trace"; exit 1; }

python3 "$HERE/repair_trace_rod.py" "$OUT/${P}_p1_outerwall_extended.swc" --voxel-size 1000 \
    --frozen-head 40 --frozen-tail 40 -o "$OUT/${P}_ext.swc" \
    --report "$OUT/${P}_ext.json" "$@" > "$OUT/${P}_ext.log" 2>&1 \
    || { echo "extended repair failed"; tail -3 "$OUT/${P}_ext.log"; exit 1; }

(cd "$MG" && "$UG4/bin/ugshell" -ex swc_to_neurite_projector_nephron_final_pipeline.lua -- \
    -swc "$OUT/${P}_ext.swc" -skipSmoothing 1 -centerlineScale 1 -extendTerminals 0 \
    $MESH -numRefs "$NUMREFS" -out "$OUT/${P}.ugx") > "$OUT/${P}.log" 2>&1
GEN=$?
echo "$P generation exit=$GEN $(grep -E 'VALIDATION:' "$OUT/${P}.log" | tr '\n' ' ')"
grep -E "EXCEPTION|degenerate volume|0: " "$OUT/${P}.log" | head -4
if [ "$GEN" = 0 ] && [ "${CHECK:-0}" = 1 ]; then
    (cd "$MG" && "$UG4/bin/ugshell" -ex check_refinement_quality.lua -grid "$OUT/${P}.ugx" \
        -refs "$NUMREFS" -mode 2) > "$OUT/${P}_check.log" 2>&1
    grep -E "Level 0|Local refinement|worst|EXCEPTION" "$OUT/${P}_check.log"
fi
exit $GEN
