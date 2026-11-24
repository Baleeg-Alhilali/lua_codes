#!/usr/bin/env bash
# Reproduce the UG4 build used for the 2026-11-24 nephron work on another Mac.
#
# Usage:  bash apply_overlay.sh /path/to/UG4_promesh_ogrid [--build]
#
# 1. Checks every UG4 git repository against the commit in MANIFEST.md.
#    A repository with local changes is first saved with `git stash` (recover
#    with `git stash list` / `git stash pop`), then checked out at the
#    manifest commit.
# 2. Applies patches/ugcore.patch and patches/plugin_NavierStokes.patch.
# 3. Replaces plugins/neuro_collection with the copy in this folder (the old
#    directory is moved to plugins/neuro_collection.backup_<timestamp>).
# 4. Restores the untracked EigenForUG4 files and the ug4 root CMakeLists.txt.
# 5. Verifies SHA-256 sums of every overlay-provided file.
# 6. With --build: configures build/ with the recorded CMake options and builds.
set -euo pipefail

UG4="${1:?usage: apply_overlay.sh /path/to/UG4_promesh_ogrid [--build]}"
DO_BUILD="${2:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"

# repository | commit | remote  (plain list: macOS ships bash 3.2)
REPOS="
ugcore 46b9b57a8c8398939bf1b6a4661396c863f6a723 https://github.com/UG4/ugcore.git
plugins/NavierStokes 8df0ce228312397f10638219ba39383dbf59b4ac https://github.com/UG4/plugin_NavierStokes.git
plugins/ConvectionDiffusion c55e90170fc1482ddf54f2989d3418b7dbb545a6 https://github.com/UG4/plugin_ConvectionDiffusion.git
plugins/Limex cb076d34a8f7fb68aa139486216468c2adf31533 https://github.com/UG4/plugin_Limex.git
plugins/LuaShell f4bfa0e6e3a61459330817399161ef2b94ac2238 https://github.com/UG4/plugin_LuaShell.git
plugins/ProMesh 9f3ee2ea2a4b27cab9be59e814bc57e354b1a8e3 https://github.com/UG4/plugin_ProMesh.git
plugins/SuperLU6 50b43fdac91aba42057afedafdfb0d32e6c099ce https://github.com/UG4/plugin_SuperLU6.git
externals/BoostForUG4 dbdd3c8cc7079c056f022d5f22545f18d5edf0f5 https://github.com/UG4/external_BoostForUG4.git
externals/EigenForUG4 0ebff70a3d12b664beefab1b5b9b73515330a288 https://github.com/UG4/external_EigenForUG4.git
"

echo "== UG4 root: $UG4"
echo "$REPOS" | while read -r repo want remote; do
  [ -z "$repo" ] && continue
  dir="$UG4/$repo"
  if [ ! -d "$dir/.git" ]; then
    echo "-- $repo missing: cloning $remote"
    git clone "$remote" "$dir"
  fi
  (
    cd "$dir"
    if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
      echo "-- $repo has local changes: saving them with git stash"
      git stash push -m "pre-overlay backup $STAMP"
    fi
    if [ "$(git rev-parse HEAD)" != "$want" ]; then
      echo "-- $repo: checking out $want"
      git fetch origin
      git checkout --quiet "$want"
    fi
    echo "-- $repo at $(git rev-parse --short HEAD)"
  )
done

echo "== Applying patches"
(cd "$UG4/ugcore" && git apply --check "$HERE/patches/ugcore.patch" && git apply "$HERE/patches/ugcore.patch")
(cd "$UG4/plugins/NavierStokes" && git apply --check "$HERE/patches/plugin_NavierStokes.patch" && git apply "$HERE/patches/plugin_NavierStokes.patch")

echo "== Installing neuro_collection"
if [ -d "$UG4/plugins/neuro_collection" ]; then
  mv "$UG4/plugins/neuro_collection" "$UG4/plugins/neuro_collection.backup_$STAMP"
  echo "-- previous plugin moved to plugins/neuro_collection.backup_$STAMP"
fi
mkdir -p "$UG4/plugins/neuro_collection"
rsync -a "$HERE/neuro_collection/" "$UG4/plugins/neuro_collection/"

echo "== Restoring untracked EigenForUG4 files and root CMakeLists.txt"
rsync -a "$HERE/externals_EigenForUG4_untracked/" "$UG4/externals/EigenForUG4/"
if [ -f "$UG4/CMakeLists.txt" ] && ! cmp -s "$UG4/CMakeLists.txt" "$HERE/ug4_root/CMakeLists.txt"; then
  cp "$UG4/CMakeLists.txt" "$UG4/CMakeLists.txt.backup_$STAMP"
fi
cp "$HERE/ug4_root/CMakeLists.txt" "$UG4/CMakeLists.txt"

echo "== Verifying file hashes"
(
  cd "$UG4"
  fail=0
  while read -r sum path; do
    actual="$(shasum -a 256 "$path" | cut -d' ' -f1)"
    if [ "$actual" != "$sum" ]; then echo "MISMATCH $path"; fail=1; fi
  done < "$HERE/SHA256SUMS_installed.txt"
  [ "$fail" = 0 ] && echo "-- all $(wc -l < "$HERE/SHA256SUMS_installed.txt" | tr -d ' ') files match"
  exit $fail
)

if [ "$DO_BUILD" = "--build" ]; then
  echo "== Configuring and building"
  mkdir -p "$UG4/build"
  cd "$UG4/build"
  cmake .. -DCMAKE_BUILD_TYPE=Release -DDIM=ALL -DCPU=ALL -DPRECISION=double \
    -DPARALLEL=ON -DLAPACK=ON -DBLAS=ON -DINTERNAL_BOOST=ON -DLINK_TETGEN=ON \
    -Dtetgen=ON -DCOMPILE_INFO=ON -DSTATIC_BUILD=OFF -DEMBEDDED_PLUGINS=OFF \
    -DENABLE_ALL_PLUGINS=OFF -DConvectionDiffusion=ON -DLimex=ON -DLuaShell=OFF \
    -DNavierStokes=ON -DProMesh=ON -DSuperLU6=ON -DUSE_INTERNAL_SUPERLU=ON \
    -Dneuro_collection=ON
  make -j"$(sysctl -n hw.ncpu)"
fi
echo "== Done"
