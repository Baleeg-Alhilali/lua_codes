#!/usr/bin/env python3
"""Run the requested full-factorial nephron mesh sweep and delete mesh artifacts."""

from __future__ import annotations

import csv
import itertools
import mmap
import re
import shutil
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
RESULTS = MODEL / "sweep_results" / "mesh_parameter_sweep"
LOGS = RESULTS / "logs"
TEMP = MODEL / "runs" / "mesh_parameter_sweep_tmp"
CSV_PATH = RESULTS / "mesh_parameter_sweep.csv"
MD_PATH = RESULTS / "mesh_parameter_sweep.md"
UGSHELL = Path("/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell")
SCRIPT = "swc_to_neurite_projector_nephron_final_pipeline.lua"

OGRIDS = (6, 8, 10)
BOX_REFS = (1, 2, 3, 4)
ANISOTROPIES = (4, 6, 8, 10)
BEND_ANISOTROPIES = (1, 2, 4)
FIELDS = (
    "ogrid", "BoxRefs", "anisotropy", "bendAnisotropy",
    "r0_vertices", "r0_volumes", "r1_vertices", "r1_volumes",
    "r2_vertices", "r2_volumes", "status",
)
TOKEN_RE = re.compile(rb"\S+")
VOLUME_TAGS = {
    b"tetrahedrons": 4,
    b"pyramids": 5,
    b"prisms": 6,
    b"hexahedrons": 8,
    b"octahedrons": 6,
}


def count_tokens(mm: mmap.mmap, start: int, end: int) -> int:
    count = 0
    tail = b""
    pos = start
    chunk_size = 8 * 1024 * 1024
    while pos < end:
        stop = min(end, pos + chunk_size)
        data = tail + mm[pos:stop]
        pos = stop
        if pos < end:
            split = max(data.rfind(b" "), data.rfind(b"\n"), data.rfind(b"\t"))
            if split < 0:
                tail = data
                continue
            complete, tail = data[: split + 1], data[split + 1 :]
        else:
            complete, tail = data, b""
        count += sum(1 for _ in TOKEN_RE.finditer(complete))
    if tail:
        count += sum(1 for _ in TOKEN_RE.finditer(tail))
    return count


def tag_count(mm: mmap.mmap, tag: bytes, arity: int) -> int:
    opening = b"<" + tag
    start = mm.find(opening)
    if start < 0:
        return 0
    start = mm.find(b">", start) + 1
    end = mm.find(b"</" + tag + b">", start)
    if start <= 0 or end < 0:
        raise ValueError(f"Malformed UGX tag: {tag.decode()}")
    tokens = count_tokens(mm, start, end)
    if tokens % arity:
        raise ValueError(f"Non-integral {tag.decode()} count in UGX")
    return tokens // arity


def mesh_stats(path: Path) -> tuple[int, int]:
    with path.open("rb") as stream, mmap.mmap(
        stream.fileno(), 0, access=mmap.ACCESS_READ
    ) as mm:
        vertices = tag_count(mm, b"vertices", 3)
        volumes = sum(tag_count(mm, tag, arity) for tag, arity in VOLUME_TAGS.items())
    return vertices, volumes


def cleanup_temp() -> None:
    if not TEMP.exists():
        return
    for child in TEMP.iterdir():
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()


def read_completed() -> dict[tuple[int, int, int, int], dict[str, str]]:
    if not CSV_PATH.exists():
        return {}
    with CSV_PATH.open(newline="", encoding="utf-8") as stream:
        return {
            (
                int(row["ogrid"]), int(row["BoxRefs"]),
                int(row["anisotropy"]), int(row["bendAnisotropy"]),
            ): row
            for row in csv.DictReader(stream)
        }


def write_outputs(rows: list[dict[str, object]]) -> None:
    rows.sort(key=lambda row: (
        int(row["ogrid"]), int(row["BoxRefs"]),
        int(row["anisotropy"]), int(row["bendAnisotropy"]),
    ))
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    with MD_PATH.open("w", encoding="utf-8") as stream:
        stream.write("# Nephron mesh parameter sweep\n\n")
        stream.write("| ogrid | BoxRefs | anisotropy | bendAnisotropy | R0 vertices | R0 volumes | R1 vertices | R1 volumes | R2 vertices | R2 volumes | status |\n")
        stream.write("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|\n")
        for row in rows:
            stream.write(
                "| {ogrid} | {BoxRefs} | {anisotropy} | {bendAnisotropy} | "
                "{r0_vertices} | {r0_volumes} | {r1_vertices} | {r1_volumes} | "
                "{r2_vertices} | {r2_volumes} | {status} |\n".format(**row)
            )


def run_case(ogrid: int, box_refs: int, anisotropy: int, bend: int) -> dict[str, object]:
    case = f"o{ogrid}_b{box_refs}_a{anisotropy}_ba{bend}"
    output = TEMP / f"{case}.ugx"
    log_path = LOGS / f"{case}.log"
    command = [
        str(UGSHELL), "-ex", SCRIPT, "--",
        "-swc", "../Nephron traces/3rd_trace_243_third.swc",
        "-radius", "15", "-membraneThickness", "10",
        "-centerlineScale", "24", "-resampleSpacing", "5",
        "-extendTerminals", "1", "-useTraceEndpoints", "0",
        "-terminalSpacing", "20", "-terminalPaddingFactor", "0",
        "-padding", "200", "-BoxRefs", str(box_refs),
        "-anisotropy", str(anisotropy),
        "-bendAnisotropy", str(bend),
        "-bendChordErrorFactor", "0.02", "-ogrid", str(ogrid),
        "-interMesh", "tet", "-interNearEdgeFactor", "2",
        "-interFarEdgeFactor", "6", "-interNearDistanceFactor", "1.5",
        "-interFarDistanceFactor", "6", "-numRefs", "2",
        "-out", str(output),
    ]
    row: dict[str, object] = {
        "ogrid": ogrid, "BoxRefs": box_refs,
        "anisotropy": anisotropy, "bendAnisotropy": bend,
    }
    cleanup_temp()
    try:
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                command, cwd=HERE, stdout=log, stderr=subprocess.STDOUT,
                timeout=3600, check=False,
            )
        paths = [output, output.with_name(output.stem + "_refined_1.ugx"),
                 output.with_name(output.stem + "_refined_2.ugx")]
        if completed.returncode != 0 or not all(path.exists() for path in paths):
            raise RuntimeError(f"ugshell exit={completed.returncode}")
        stats = [mesh_stats(path) for path in paths]
        row.update({
            "r0_vertices": stats[0][0], "r0_volumes": stats[0][1],
            "r1_vertices": stats[1][0], "r1_volumes": stats[1][1],
            "r2_vertices": stats[2][0], "r2_volumes": stats[2][1],
            "status": "PASS",
        })
    except Exception as error:  # every failed configuration is a CRASH row
        row.update({
            "r0_vertices": "CRASH", "r0_volumes": "CRASH",
            "r1_vertices": "CRASH", "r1_volumes": "CRASH",
            "r2_vertices": "CRASH", "r2_volumes": "CRASH",
            "status": f"CRASH ({error})",
        })
    finally:
        cleanup_temp()
    return row


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    TEMP.mkdir(parents=True, exist_ok=True)
    completed = read_completed()
    rows: list[dict[str, object]] = list(completed.values())
    combinations = list(itertools.product(
        OGRIDS, BOX_REFS, ANISOTROPIES, BEND_ANISOTROPIES
    ))
    for index, values in enumerate(combinations, 1):
        if values in completed:
            continue
        print(f"[{index}/{len(combinations)}] ogrid={values[0]} BoxRefs={values[1]} "
              f"anisotropy={values[2]} bendAnisotropy={values[3]}", flush=True)
        row = run_case(*values)
        print(f"  {row['status']}", flush=True)
        rows.append(row)
        write_outputs(rows)
    cleanup_temp()
    TEMP.rmdir()
    write_outputs(rows)


if __name__ == "__main__":
    main()
