#!/usr/bin/env python3
"""Edge aspect ratio (longest / shortest edge) of every hexahedron in a UGX,
reported per subset of the first subset handler.

usage: hex_aspect_ratio.py mesh.ugx [--limit 10]
"""
import argparse
import re
import sys

import numpy as np

HEX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]
VOLUME_TAGS = ("tetrahedrons", "hexahedrons", "prisms", "pyramids", "octahedrons")
VOLUME_SIZES = {"tetrahedrons": 4, "hexahedrons": 8, "prisms": 6,
                "pyramids": 5, "octahedrons": 6}


def tag_body(text, tag, start=0):
    i = text.find("<" + tag, start)
    if i < 0:
        return None, -1
    j = text.find(">", i) + 1
    k = text.find("</" + tag + ">", j)
    return text[j:k], k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ugx")
    ap.add_argument("--limit", type=float, default=10.0)
    a = ap.parse_args()
    text = open(a.ugx).read()
    grid_end = text.find("<subset_handler")
    coords, _ = tag_body(text, "vertices")
    pos = np.array(coords.split(), float).reshape(-1, 3) * 1000.0  # mm -> um

    # Global volume numbering follows the order of the volume blocks in the file.
    blocks = []
    for tag in VOLUME_TAGS:
        i = text.find("<" + tag, 0, grid_end)
        if i >= 0:
            blocks.append((i, tag))
    blocks.sort()
    offset, hex_range, hexes = 0, None, None
    for _, tag in blocks:
        body, _ = tag_body(text, tag)
        n = len(body.split()) // VOLUME_SIZES[tag]
        if tag == "hexahedrons":
            hexes = np.array(body.split(), int).reshape(-1, 8)
            hex_range = (offset, offset + n)
        offset += n
    if hexes is None:
        print("no hexahedra"); return 0

    e = np.stack([np.linalg.norm(pos[hexes[:, p]] - pos[hexes[:, q]], axis=1)
                  for p, q in HEX_EDGES], axis=1)
    ratio = e.max(1) / e.min(1)

    sh, _ = tag_body(text, "subset_handler", grid_end - 1)
    bad_total = 0
    for m in re.finditer(r'<subset name="([^"]*)"[^>]*>(.*?)</subset>', sh, re.S):
        name, body = m.group(1), m.group(2)
        vol, _ = tag_body(body, "volumes")
        if not vol:
            continue
        ids = np.array(vol.split(), int)
        ids = ids[(ids >= hex_range[0]) & (ids < hex_range[1])] - hex_range[0]
        if ids.size == 0:
            continue
        r = ratio[ids]
        bad = int((r > a.limit).sum()); bad_total += bad
        k = ids[r.argmax()]
        print(f"{name:14s} hexes={ids.size:8d} aspect max={r.max():7.2f} "
              f"p99={np.percentile(r, 99):6.2f} median={np.median(r):5.2f} "
              f"above {a.limit:g}: {bad}  (worst edges um: min {e[k].min():.2f}, "
              f"max {e[k].max():.2f})")
    print(f"total hexes above {a.limit:g}: {bad_total}")
    return 1 if bad_total else 0


if __name__ == "__main__":
    sys.exit(main())
