#!/usr/bin/env python3
"""Make the refined boundary subsets of the transport UGX mesh P1-safe.

The original mesh represents the tetrahedral outer boundary with independent
quadrilaterals.  Their vertices and edges coincide with the volume mesh, but
the quadrilateral interiors do not.  Regular refinement therefore creates one
orphan face-centre vertex per quadrilateral.  P1 fields supported on OuterWall
then acquire zero-equation DoFs and produce a singular time-step Jacobian.

This utility removes only those detached face objects. It also moves the one
Inlet face and one Outlet face into lumen-only subsets. Their original subset
names remain on the shared vertices, but no longer generate membrane and
interstitial unknowns at the two lumen face centres. Face indices are remapped
in every subset handler.
"""

from __future__ import annotations

import argparse
import xml.etree.ElementTree as ET
from pathlib import Path


def integer_list(text: str | None) -> list[int]:
    return [] if not text else [int(value) for value in text.split()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()

    tree = ET.parse(args.source)
    root = tree.getroot()
    triangles = root.find("triangles")
    quadrilaterals = root.find("quadrilaterals")
    if triangles is None or quadrilaterals is None:
        raise RuntimeError("Expected triangles and quadrilaterals in UGX file")

    triangle_count = len(integer_list(triangles.text)) // 3
    quad_values = integer_list(quadrilaterals.text)
    quads = [quad_values[i : i + 4] for i in range(0, len(quad_values), 4)]

    default_handler = root.find("subset_handler")
    if default_handler is None:
        raise RuntimeError("UGX file has no subset handler")
    outer_wall = next(
        (subset for subset in default_handler if subset.get("name") == "OuterWall"),
        None,
    )
    if outer_wall is None:
        raise RuntimeError("Default subset handler has no OuterWall subset")

    removed_face_ids = set(integer_list(outer_wall.findtext("faces")))
    if not removed_face_ids:
        raise RuntimeError("OuterWall contains no faces")
    if min(removed_face_ids) < triangle_count:
        raise RuntimeError("OuterWall unexpectedly contains triangular faces")

    removed_quad_ids = {face_id - triangle_count for face_id in removed_face_ids}
    if max(removed_quad_ids) >= len(quads):
        raise RuntimeError("OuterWall face index is outside quadrilateral array")

    lumen_end_faces: dict[str, list[int]] = {}
    for old_name, new_name in (
        ("Inlet", "LumenInletFace"),
        ("Outlet", "LumenOutletFace"),
    ):
        old_subset = next(
            (subset for subset in default_handler if subset.get("name") == old_name),
            None,
        )
        if old_subset is None:
            raise RuntimeError(f"Default subset handler has no {old_name} subset")
        face_ids = integer_list(old_subset.findtext("faces"))
        if len(face_ids) != 1:
            raise RuntimeError(f"Expected exactly one {old_name} face")
        lumen_end_faces[new_name] = face_ids

    kept_quads = [quad for index, quad in enumerate(quads) if index not in removed_quad_ids]
    quadrilaterals.text = " ".join(str(vertex) for quad in kept_quads for vertex in quad)

    old_face_count = triangle_count + len(quads)
    face_id_map: dict[int, int] = {}
    next_id = 0
    for old_id in range(old_face_count):
        if old_id not in removed_face_ids:
            face_id_map[old_id] = next_id
            next_id += 1

    for handler in root.findall("subset_handler"):
        for subset in handler:
            faces = subset.find("faces")
            if faces is None:
                continue
            remapped = [face_id_map[i] for i in integer_list(faces.text) if i in face_id_map]
            if remapped:
                faces.text = " ".join(map(str, remapped))
            else:
                subset.remove(faces)

    # Inlet/Outlet vertices stay in their original subsets because several
    # compartments meet there. Only each lumen quadrilateral interior receives
    # a dedicated subset, preventing two orphan unknowns per species.
    for old_name in ("Inlet", "Outlet"):
        old_subset = next(
            subset for subset in default_handler if subset.get("name") == old_name
        )
        old_faces = old_subset.find("faces")
        if old_faces is not None:
            old_subset.remove(old_faces)

    for new_name, old_ids in lumen_end_faces.items():
        new_subset = ET.SubElement(
            default_handler,
            "subset",
            {"name": new_name, "color": "0.8 0.8 0.2 1", "state": "0"},
        )
        faces = ET.SubElement(new_subset, "faces")
        faces.text = " ".join(str(face_id_map[i]) for i in old_ids)

    args.destination.parent.mkdir(parents=True, exist_ok=True)
    tree.write(args.destination, encoding="UTF-8", xml_declaration=True)
    print(
        f"Removed {len(removed_face_ids)} detached OuterWall quadrilaterals and "
        "split the lumen Inlet/Outlet faces; "
        f"wrote {args.destination}"
    )


if __name__ == "__main__":
    main()
