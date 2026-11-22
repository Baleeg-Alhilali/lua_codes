#!/usr/bin/env python3
"""Check refined diffusion VTK checkpoints for finite, bounded concentration."""

from __future__ import annotations

import argparse
import base64
import glob
import math
from pathlib import Path
import struct
import xml.etree.ElementTree as ET


def binary_float_array(element: ET.Element) -> tuple[float, ...]:
    payload = base64.b64decode("".join((element.text or "").split()))
    byte_count = struct.unpack("<I", payload[:4])[0]
    if byte_count % 4:
        raise ValueError("Float32 VTK payload length is not divisible by four")
    return struct.unpack(f"<{byte_count // 4}f", payload[4:4 + byte_count])


def concentration_values(files: list[Path]) -> list[float]:
    values: list[float] = []
    for path in files:
        root = ET.parse(path).getroot()
        array = next(
            (item for item in root.iter("DataArray")
             if item.get("Name") == "Urea"),
            None,
        )
        if array is not None:
            values.extend(binary_float_array(array))
    if not values:
        raise ValueError("No Urea point data found in the selected VTK files")
    return values


def compartment_series(prefix: Path, compartment: str) -> tuple[list[float], list[float], int]:
    """Read first/last serial or MPI pieces in a compartment PVD series."""
    pvd_path = Path(f"{prefix}_{compartment}.pvd")
    datasets = list(ET.parse(pvd_path).getroot().iter("DataSet"))
    if not datasets:
        raise ValueError(f"No time steps in {pvd_path}")

    def referenced_vtus(dataset: ET.Element) -> list[Path]:
        pvtu_path = pvd_path.parent / dataset.attrib["file"]
        pieces = list(ET.parse(pvtu_path).getroot().iter("Piece"))
        if not pieces:
            raise ValueError(f"No pieces referenced by {pvtu_path}")
        paths = [pvtu_path.parent / piece.attrib["Source"] for piece in pieces]
        missing = [path for path in paths if not path.is_file()]
        if missing:
            raise ValueError(
                f"{pvtu_path} references missing pieces: "
                + ", ".join(map(str, missing)))
        return paths

    initial = concentration_values(referenced_vtus(datasets[0]))
    final = concentration_values(referenced_vtus(datasets[-1]))
    return initial, final, len(datasets)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="output prefix passed to the Lua model")
    parser.add_argument("--bound-tolerance", type=float, default=1.0e-6)
    args = parser.parse_args()

    prefix = Path(args.prefix)
    compartment_pvds = [Path(f"{prefix}_{name}.pvd")
                        for name in ("lumen", "membrane", "inter")]
    if all(path.is_file() for path in compartment_pvds):
        expected_initial = {"lumen": 1.0, "membrane": 0.0, "inter": 0.0}
        step_counts: set[int] = set()
        final_bounds: dict[str, tuple[float, float]] = {}
        for compartment in ("lumen", "membrane", "inter"):
            initial, final, step_count = compartment_series(prefix, compartment)
            step_counts.add(step_count)
            values = initial + final
            if any(not math.isfinite(value) for value in values):
                raise SystemExit(f"{compartment}: non-finite concentration")
            target = expected_initial[compartment]
            if (min(initial) < target - args.bound_tolerance or
                    max(initial) > target + args.bound_tolerance):
                raise SystemExit(
                    f"{compartment}: initial field is not uniformly {target}")
            if (min(final) < -args.bound_tolerance or
                    max(final) > 1.0 + args.bound_tolerance):
                raise SystemExit(f"{compartment}: final field is outside [0, 1]")
            final_bounds[compartment] = (min(final), max(final))
            print(
                f"{compartment}: steps={step_count} "
                f"initial=[{min(initial):.9g}, {max(initial):.9g}] "
                f"final=[{min(final):.9g}, {max(final):.9g}]"
            )
        if len(step_counts) != 1:
            raise SystemExit("Compartment PVD files have different step counts")
        if final_bounds["membrane"][1] <= 0.0:
            raise SystemExit("No apical transfer into the membrane was detected")
        if final_bounds["inter"][1] <= 0.0:
            raise SystemExit("No basolateral transfer into the interstitium was detected")
        print("TRANSIENT_COMPARTMENT_VTK_BOUNDS: PASS")
        return 0

    initial_files = [
        Path(path) for pattern in (
            f"{prefix}_initial_s*_t0000.vtu",
            f"{prefix}_initial_p*_s*_t0000.vtu",
        ) for path in sorted(glob.glob(pattern))
    ]
    final_candidates = [
        Path(path) for pattern in (
            f"{prefix}_s*_t*.vtu",
            f"{prefix}_p*_s*_t*.vtu",
        ) for path in glob.glob(pattern)
    ]
    if not initial_files or not final_candidates:
        raise SystemExit("Missing initial or final VTK checkpoint files")

    final_step = max(int(path.stem.rsplit("_t", 1)[1]) for path in final_candidates)
    final_files = [
        path for path in final_candidates
        if int(path.stem.rsplit("_t", 1)[1]) == final_step
    ]
    initial = concentration_values(initial_files)
    final = concentration_values(final_files)

    nonfinite_initial = sum(not math.isfinite(value) for value in initial)
    nonfinite_final = sum(not math.isfinite(value) for value in final)
    initial_min, initial_max = min(initial), max(initial)
    final_min, final_max = min(final), max(final)
    lower_limit = initial_min - args.bound_tolerance
    upper_limit = initial_max + args.bound_tolerance

    print(f"initial_files={len(initial_files)} values={len(initial)}")
    print(f"final_step={final_step} files={len(final_files)} values={len(final)}")
    print(f"initial_bounds=[{initial_min:.17g}, {initial_max:.17g}]")
    print(f"final_bounds=[{final_min:.17g}, {final_max:.17g}]")
    print(f"nonfinite_initial={nonfinite_initial} nonfinite_final={nonfinite_final}")

    if nonfinite_initial or nonfinite_final:
        raise SystemExit("VTK concentration contains non-finite values")
    if final_min < lower_limit or final_max > upper_limit:
        raise SystemExit("Final concentration exceeds the initial bounds")

    initial_pvd = Path(f"{prefix}_initial_t0000.pvd")
    final_pvd = Path(f"{prefix}_t{final_step:04d}.pvd")
    if not initial_pvd.is_file() or not final_pvd.is_file():
        raise SystemExit("Missing initial or final PVD checkpoint index")

    print("TRANSPORT_DIFFUSION_VTK_BOUNDS: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
