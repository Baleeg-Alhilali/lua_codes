"""ParaView macro: plot apical urea leak flux versus distance from the inlet.

Load matching Lumen and Membrane ``.pvtu`` files, then run this script from
Tools > Python Shell > Run Script, or install it as a ParaView macro.  The
script finds the coincident Lumen/Membrane surface, computes

    J_apical = P_apical * (C_lumen - C_membrane),

along two radius-offset traces: the local inner-curvature side and the opposing
outer-curvature side. The x-axis is arc length on the UGX neurite-projector
centerline. Positive flux is directed from Lumen to Membrane. A CSV is written
beside the Lumen file and an XY chart is created in ParaView.

For command-line validation with pvpython:

    pvpython paraview_apical_flux_profile.py \
        --lumen path/to/lumen_t0100.pvtu \
        --membrane path/to/membrane_t0100.pvtu --no-chart
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import vtk
from paraview import servermanager
from paraview.simple import (
    AssignViewToLayout,
    CreateLayout,
    CreateView,
    GetSources,
    OpenDataFile,
    ProgrammableSource,
    Render,
    RenderAllViews,
    SetActiveSource,
    SetActiveView,
    Show,
)


DEFAULT_APICAL_PERMEABILITY_MM_S = 3.78e-3
DEFAULT_BINS = 160
DEFAULT_APICAL_RADIUS_MM = 1.5e-2
DEFAULT_GRID = (Path(__file__).resolve().parents[2]
                / "runs" / "third_trace_generalized3.ugx")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--lumen")
    parser.add_argument("--membrane")
    parser.add_argument("--permeability", type=float,
                        default=DEFAULT_APICAL_PERMEABILITY_MM_S)
    parser.add_argument("--bins", type=int, default=DEFAULT_BINS)
    parser.add_argument("--radius", type=float,
                        default=DEFAULT_APICAL_RADIUS_MM)
    parser.add_argument("--grid", default=str(DEFAULT_GRID))
    parser.add_argument("--csv")
    parser.add_argument("--no-chart", action="store_true")
    args, _ = parser.parse_known_args(sys.argv[1:])
    if args.permeability < 0.0:
        raise ValueError("Apical permeability must be non-negative")
    if args.bins < 2:
        raise ValueError("At least two distance bins are required")
    if args.radius <= 0.0:
        raise ValueError("Apical radius must be positive")
    return args


def source_filename(proxy) -> str | None:
    for property_name in ("FileName", "FileNames"):
        if not hasattr(proxy, property_name):
            continue
        value = getattr(proxy, property_name)
        if isinstance(value, (list, tuple)):
            return str(value[0]) if value else None
        if not isinstance(value, str):
            try:
                return str(value[0]) if len(value) else None
            except (TypeError, AttributeError, IndexError):
                pass
        if value:
            return str(value)
    return None


def find_loaded_source(token: str):
    candidates = []
    for (registration_name, _), proxy in GetSources().items():
        name = registration_name.lower()
        filename = (source_filename(proxy) or "").lower()
        if token in name or token in filename:
            candidates.append((registration_name, proxy))
    if not candidates:
        raise RuntimeError(
            f"No loaded {token} source found. Open matching .pvtu files first.")
    # The tNNNN suffix sorts chronologically, so prefer the latest loaded step.
    return sorted(candidates, key=lambda item: item[0])[-1][1]


def get_sources(args: argparse.Namespace):
    if bool(args.lumen) != bool(args.membrane):
        raise ValueError("Specify both --lumen and --membrane, or neither")
    if args.lumen:
        return OpenDataFile(args.lumen), OpenDataFile(args.membrane)
    return find_loaded_source("lumen"), find_loaded_source("membrane")


def iter_dataset_leaves(data_object):
    if isinstance(data_object, vtk.vtkDataSet):
        yield data_object
        return
    if isinstance(data_object, vtk.vtkCompositeDataSet):
        iterator = data_object.NewIterator()
        iterator.SkipEmptyNodesOn()
        iterator.InitTraversal()
        while not iterator.IsDoneWithTraversal():
            leaf = iterator.GetCurrentDataObject()
            if isinstance(leaf, vtk.vtkDataSet):
                yield leaf
            iterator.GoToNextItem()


def fetch_and_merge(proxy):
    fetched = servermanager.Fetch(proxy)
    leaves = list(iter_dataset_leaves(fetched))
    if not leaves:
        raise RuntimeError("The selected source contains no VTK data set")
    append = vtk.vtkAppendFilter()
    append.MergePointsOn()
    for leaf in leaves:
        append.AddInputData(leaf)
    append.Update()
    merged = append.GetOutput()

    # Parallel VTU pieces duplicate partition-boundary nodes.  Merge them before
    # extracting the exterior surface so rank interfaces are not mistaken for
    # physical membrane faces.
    if hasattr(vtk, "vtkStaticCleanUnstructuredGrid"):
        clean = vtk.vtkStaticCleanUnstructuredGrid()
        clean.SetInputData(merged)
        clean.Update()
        output = vtk.vtkUnstructuredGrid()
        output.ShallowCopy(clean.GetOutput())
        return output
    return merged


def distance_squared(a, b) -> float:
    return sum((a[i] - b[i]) ** 2 for i in range(3))


def normalized(vector):
    length = math.sqrt(sum(value * value for value in vector))
    if length <= 1.0e-15:
        return None
    return tuple(value / length for value in vector)


def vector_subtract(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def vector_add(a, b):
    return tuple(a[i] + b[i] for i in range(3))


def vector_scale(vector, scale):
    return tuple(scale * value for value in vector)


def load_projector_centerline(grid_path: Path):
    """Decode the cubic NeuriteProjector sections serialized in the UGX."""
    root = ET.parse(grid_path).getroot()
    handler = root.find("projection_handler")
    if handler is None:
        raise RuntimeError(f"No projection handler in {grid_path}")
    projector = next(
        (item for item in handler.findall("projector")
         if item.get("type") == "neurite"), None)
    if projector is None:
        raise RuntimeError(f"No neurite projector in {grid_path}")
    tokens = [float(value) for value in (projector.text or "").split()]
    if len(tokens) < 12:
        raise RuntimeError("Neurite projector serialization is incomplete")
    point_count = int(tokens[8])
    record_start = 11
    record_stride = 17
    required = record_start + point_count * record_stride
    if len(tokens) < required:
        raise RuntimeError("Unexpected neurite projector serialization layout")
    sections = []
    for index in range(point_count):
        start = record_start + index * record_stride
        end_parameter = tokens[start]
        coefficients = {axis: [] for axis in ("x", "y", "z", "r")}
        names = ("x", "y", "z", "r")
        for polynomial_order in range(4):
            for axis_index, axis in enumerate(names):
                coefficients[axis].append(
                    tokens[start + 1 + 4 * polynomial_order + axis_index])
        sections.append((end_parameter, coefficients))
    return sections


def sample_projector_centerline(sections, dense_count=5000):
    """Evaluate UG4's p(t)=a(te-t)^3+b(te-t)^2+c(te-t)+d."""
    points = []
    radii = []
    section_index = 0
    for index in range(dense_count):
        parameter = index / (dense_count - 1)
        while (section_index + 1 < len(sections)
               and sections[section_index][0] < parameter):
            section_index += 1
        end_parameter, coefficients = sections[section_index]
        local_parameter = end_parameter - parameter

        def evaluate(axis):
            a, b, c, d = coefficients[axis]
            return ((a * local_parameter + b) * local_parameter + c
                    ) * local_parameter + d

        points.append((evaluate("x"), evaluate("y"), evaluate("z")))
        radii.append(evaluate("r"))
    return points, radii


def resample_polyline(points, sample_count):
    cumulative = [0.0]
    for first, second in zip(points, points[1:]):
        cumulative.append(cumulative[-1] + math.dist(first, second))
    length = cumulative[-1]
    if length <= 0.0:
        raise RuntimeError("Centerline has zero length")
    result = []
    segment = 0
    for index in range(sample_count):
        distance = length * index / (sample_count - 1)
        while (segment + 1 < len(cumulative) - 1
               and cumulative[segment + 1] < distance):
            segment += 1
        span = cumulative[segment + 1] - cumulative[segment]
        fraction = ((distance - cumulative[segment]) / span
                    if span > 0.0 else 0.0)
        point = tuple(
            points[segment][axis]
            + fraction * (points[segment + 1][axis] - points[segment][axis])
            for axis in range(3))
        result.append((distance, point))
    return result, length


def extract_apical_profile(lumen, membrane, permeability: float,
                           sample_count: int, centerline_points,
                           apical_radius: float):
    lumen_urea = lumen.GetPointData().GetArray("Urea")
    membrane_urea = membrane.GetPointData().GetArray("Urea")
    if lumen_urea is None or membrane_urea is None:
        raise RuntimeError("Both inputs must contain point-data array 'Urea'")

    surface_filter = vtk.vtkDataSetSurfaceFilter()
    surface_filter.SetInputData(lumen)
    surface_filter.Update()
    surface = vtk.vtkPolyData()
    surface.ShallowCopy(surface_filter.GetOutput())
    surface_urea = surface.GetPointData().GetArray("Urea")
    if surface_urea is None:
        raise RuntimeError("Lumen surface did not retain the Urea point array")

    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(membrane)
    locator.BuildLocator()

    bounds = lumen.GetBounds()
    diagonal = math.sqrt(
        (bounds[1] - bounds[0]) ** 2
        + (bounds[3] - bounds[2]) ** 2
        + (bounds[5] - bounds[4]) ** 2)
    match_tolerance = max(diagonal * 2.0e-6, 1.0e-9)
    match_tolerance_squared = match_tolerance * match_tolerance

    point_count = surface.GetNumberOfPoints()
    matched = [False] * point_count
    membrane_values = [math.nan] * point_count
    for point_id in range(point_count):
        point = surface.GetPoint(point_id)
        membrane_id = locator.FindClosestPoint(point)
        if membrane_id < 0:
            continue
        membrane_point = membrane.GetPoint(membrane_id)
        if distance_squared(point, membrane_point) <= match_tolerance_squared:
            matched[point_id] = True
            membrane_values[point_id] = membrane_urea.GetTuple1(membrane_id)

    apical_cells = []
    used_points = set()
    for cell_id in range(surface.GetNumberOfCells()):
        cell = surface.GetCell(cell_id)
        ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
        if len(ids) >= 3 and all(matched[point_id] for point_id in ids):
            apical_cells.append(ids)
            used_points.update(ids)
    if not apical_cells:
        raise RuntimeError(
            "No coincident Lumen/Membrane surface was found; check that the "
            "matching compartment .pvtu files are loaded")

    # Identify the two open apical boundary rims so the projector centerline is
    # anchored to the physical inlet and outlet cap centers.
    edge_use_count = {}
    for ids in apical_cells:
        for index, point_id in enumerate(ids):
            neighbor_id = ids[(index + 1) % len(ids)]
            edge = tuple(sorted((point_id, neighbor_id)))
            edge_use_count[edge] = edge_use_count.get(edge, 0) + 1

    boundary_adjacency = {}
    for (first, second), count in edge_use_count.items():
        if count != 1:
            continue
        boundary_adjacency.setdefault(first, set()).add(second)
        boundary_adjacency.setdefault(second, set()).add(first)
    boundary_components = []
    unseen = set(boundary_adjacency)
    while unseen:
        seed = unseen.pop()
        component = {seed}
        stack = [seed]
        while stack:
            point_id = stack.pop()
            for neighbor_id in boundary_adjacency[point_id]:
                if neighbor_id in unseen:
                    unseen.remove(neighbor_id)
                    component.add(neighbor_id)
                    stack.append(neighbor_id)
        boundary_components.append(component)
    if len(boundary_components) < 2:
        raise RuntimeError(
            "Expected inlet and outlet rims on the apical surface, found "
            f"{len(boundary_components)} boundary component(s)")

    def component_center(component):
        count = len(component)
        return tuple(sum(surface.GetPoint(p)[axis] for p in component) / count
                     for axis in range(3))

    # The actual mesh inlet is the minimum-z cap. Of the remaining boundary
    # loops, the farthest centroid is the outlet cap (x=9.201 in this mesh).
    inlet_component = min(
        boundary_components, key=lambda component: component_center(component)[2])
    inlet_center = component_center(inlet_component)
    outlet_component = max(
        (component for component in boundary_components
         if component is not inlet_component),
        key=lambda component: distance_squared(
            component_center(component), inlet_center))
    outlet_center = component_center(outlet_component)

    centerline = list(centerline_points)
    if math.dist(inlet_center, centerline[0]) > 1.0e-8:
        centerline.insert(0, inlet_center)
    else:
        centerline[0] = inlet_center
    if math.dist(outlet_center, centerline[-1]) > 1.0e-8:
        centerline.append(outlet_center)
    else:
        centerline[-1] = outlet_center
    sampled_centerline, centerline_length = resample_polyline(
        centerline, sample_count)

    # Build a compact apical surface and interpolate on its triangles. Sampling
    # the continuous faces avoids the false longitudinal error produced by
    # snapping a target between axial rings to the nearest mesh vertex.
    apical_points = vtk.vtkPoints()
    apical_lumen_array = vtk.vtkDoubleArray()
    apical_lumen_array.SetName("LumenUrea")
    apical_membrane_array = vtk.vtkDoubleArray()
    apical_membrane_array.SetName("MembraneUrea")
    point_map = {}
    for new_id, point_id in enumerate(sorted(used_points)):
        point_map[point_id] = new_id
        apical_points.InsertNextPoint(surface.GetPoint(point_id))
        apical_lumen_array.InsertNextValue(surface_urea.GetTuple1(point_id))
        apical_membrane_array.InsertNextValue(membrane_values[point_id])
    apical_polys = vtk.vtkCellArray()
    for ids in apical_cells:
        apical_polys.InsertNextCell(len(ids))
        for point_id in ids:
            apical_polys.InsertCellPoint(point_map[point_id])
    apical_surface = vtk.vtkPolyData()
    apical_surface.SetPoints(apical_points)
    apical_surface.SetPolys(apical_polys)
    apical_surface.GetPointData().AddArray(apical_lumen_array)
    apical_surface.GetPointData().AddArray(apical_membrane_array)
    apical_locator = vtk.vtkStaticCellLocator()
    apical_locator.SetDataSet(apical_surface)
    apical_locator.BuildLocator()

    def sample_apical(target):
        closest = [0.0, 0.0, 0.0]
        cell_id = vtk.mutable(0)
        sub_id = vtk.mutable(0)
        squared_distance = vtk.mutable(0.0)
        apical_locator.FindClosestPoint(
            target, closest, cell_id, sub_id, squared_distance)
        cell = apical_surface.GetCell(int(cell_id))
        weights = [0.0] * cell.GetNumberOfPoints()
        evaluated = [0.0, 0.0, 0.0]
        parametric = [0.0, 0.0, 0.0]
        evaluate_sub_id = vtk.mutable(0)
        evaluate_distance = vtk.mutable(0.0)
        cell.EvaluatePosition(
            closest, evaluated, evaluate_sub_id, parametric,
            evaluate_distance, weights)
        lumen_value = 0.0
        membrane_value = 0.0
        for local_id, weight in enumerate(weights):
            point_id = cell.GetPointId(local_id)
            lumen_value += weight * apical_lumen_array.GetTuple1(point_id)
            membrane_value += weight * apical_membrane_array.GetTuple1(point_id)
        return (tuple(closest), lumen_value, membrane_value,
                math.sqrt(float(squared_distance)))

    centers = [point for _, point in sampled_centerline]
    tangents = []
    for index in range(sample_count):
        before = centers[max(0, index - 1)]
        after = centers[min(sample_count - 1, index + 1)]
        tangent = normalized(vector_subtract(after, before))
        if tangent is None:
            raise RuntimeError("Degenerate centerline tangent")
        tangents.append(tangent)

    curvature_normals = [None] * sample_count
    curvatures = [0.0] * sample_count
    previous_normal = None
    for index in range(1, sample_count - 1):
        change = vector_subtract(tangents[index + 1], tangents[index - 1])
        spacing = (sampled_centerline[index + 1][0]
                   - sampled_centerline[index - 1][0])
        curvature_vector = vector_scale(change, 1.0 / spacing)
        tangent = tangents[index]
        tangential_part = sum(curvature_vector[i] * tangent[i]
                              for i in range(3))
        normal_vector = vector_subtract(
            curvature_vector, vector_scale(tangent, tangential_part))
        normal = normalized(normal_vector)
        if normal is not None:
            previous_normal = normal
            curvatures[index] = math.sqrt(sum(
                value * value for value in normal_vector))
        elif previous_normal is not None:
            normal = previous_normal
        curvature_normals[index] = normal

    first_valid = next((normal for normal in curvature_normals
                        if normal is not None), None)
    if first_valid is None:
        raise RuntimeError("Centerline curvature is zero everywhere")
    for index in range(sample_count):
        if curvature_normals[index] is None:
            curvature_normals[index] = (
                curvature_normals[index - 1] if index > 0 else first_valid)
            curvatures[index] = (
                curvatures[index - 1] if index > 0 else curvatures[1])

    rows = []
    inner_points = []
    outer_points = []
    maximum_snap_distance = 0.0
    for index, (distance, center) in enumerate(sampled_centerline):
        normal = curvature_normals[index]
        _, _, _, center_to_apical_distance = sample_apical(center)
        inner_target = vector_add(center, vector_scale(normal, apical_radius))
        outer_target = vector_add(center, vector_scale(normal, -apical_radius))
        (inner_point, inner_lumen, inner_membrane,
         inner_snap_distance) = sample_apical(inner_target)
        (outer_point, outer_lumen, outer_membrane,
         outer_snap_distance) = sample_apical(outer_target)
        inner_points.append(inner_point)
        outer_points.append(outer_point)
        maximum_snap_distance = max(
            maximum_snap_distance,
            inner_snap_distance,
            outer_snap_distance)
        rows.append({
            "Centerline_distance_from_inlet_mm": distance,
            "Inner_apical_flux_conc_mm_s": permeability * (
                inner_lumen - inner_membrane),
            "Outer_apical_flux_conc_mm_s": permeability * (
                outer_lumen - outer_membrane),
            "Inner_lumen_concentration": inner_lumen,
            "Inner_membrane_concentration": inner_membrane,
            "Outer_lumen_concentration": outer_lumen,
            "Outer_membrane_concentration": outer_membrane,
            "Centerline_curvature_per_mm": curvatures[index],
            "Centerline_to_apical_distance_mm": center_to_apical_distance,
            "Inner_trace_snap_distance_mm": inner_snap_distance,
            "Outer_trace_snap_distance_mm": outer_snap_distance,
        })

    inner_length = sum(math.dist(a, b)
                       for a, b in zip(inner_points, inner_points[1:]))
    outer_length = sum(math.dist(a, b)
                       for a, b in zip(outer_points, outer_points[1:]))
    return (rows, len(apical_cells), len(used_points), centerline_length,
            inlet_center, outlet_center, len(boundary_components),
            inner_length, outer_length, maximum_snap_distance)


def write_csv(rows, csv_path: Path) -> None:
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def create_chart(rows) -> None:
    columns = list(rows[0])
    values = {column: [row[column] for row in rows] for column in columns}
    script_lines = [
        "import vtk",
        "output = self.GetOutput()",
    ]
    for column in columns:
        literal = repr(values[column])
        script_lines.extend([
            "array = vtk.vtkDoubleArray()",
            f"array.SetName({column!r})",
            f"for value in {literal}: array.InsertNextValue(value)",
            "output.AddColumn(array)",
        ])

    source = ProgrammableSource(registrationName="Apical Flux Profile")
    source.OutputDataSetType = "vtkTable"
    source.Script = "\n".join(script_lines)
    source.UpdatePipeline()

    chart = CreateView("XYChartView")
    layout = CreateLayout(name="Apical Flux Plot")
    AssignViewToLayout(view=chart, layout=layout, hint=0)
    chart.ViewSize = [1100, 650]
    display = Show(source, chart)
    display.UseIndexForXAxis = 0
    display.XArrayName = "Centerline_distance_from_inlet_mm"
    try:
        display.SeriesVisibility = [
            "Inner_apical_flux_conc_mm_s",
            "Outer_apical_flux_conc_mm_s",
        ]
    except Exception:
        pass
    try:
        chart.ChartTitle = "Inner versus outer apical flux along centerline"
        chart.LeftAxisTitle = "Apical flux (concentration mm/s)"
        chart.BottomAxisTitle = "Centerline distance from inlet (mm)"
    except Exception:
        pass
    SetActiveSource(source)
    SetActiveView(chart)
    Render(chart)
    RenderAllViews()
    print("Created and activated ParaView layout: Apical Flux Plot")


def main() -> None:
    args = parse_args()
    lumen_source, membrane_source = get_sources(args)
    lumen_filename = source_filename(lumen_source)
    lumen = fetch_and_merge(lumen_source)
    membrane = fetch_and_merge(membrane_source)
    projector_sections = load_projector_centerline(
        Path(args.grid).expanduser().resolve())
    centerline_points, projector_radii = sample_projector_centerline(
        projector_sections)
    (rows, cell_count, point_count, maximum_distance, inlet_center,
     outlet_center, boundary_count, inner_length, outer_length,
     maximum_snap_distance) = extract_apical_profile(
        lumen, membrane, args.permeability, args.bins,
        centerline_points, args.radius)

    if args.csv:
        csv_path = Path(args.csv).expanduser().resolve()
    elif lumen_filename:
        lumen_path = Path(lumen_filename).expanduser().resolve()
        csv_path = lumen_path.with_name(
            lumen_path.stem.replace("_lumen_", "_apical_flux_") + ".csv")
    else:
        csv_path = Path.cwd() / "apical_flux_profile.csv"
    write_csv(rows, csv_path)

    print(
        "Apical flux profile: "
        f"cells={cell_count}, points={point_count}, "
        f"distance=0..{maximum_distance:.9g} mm, bins={len(rows)}")
    print(f"Detected boundary rims: {boundary_count}")
    print(f"Inlet center: {inlet_center}")
    print(f"Outlet center: {outlet_center}")
    print(f"Projector cubic sections: {len(projector_sections)}")
    print(f"Projector radius range: {min(projector_radii):.9g}.."
          f"{max(projector_radii):.9g} mm")
    print(f"Apical trace offset radius: {args.radius:.9g} mm")
    print(f"Inner trace sampled length: {inner_length:.9g} mm")
    print(f"Outer trace sampled length: {outer_length:.9g} mm")
    print(f"Maximum trace-to-mesh snap distance: "
          f"{maximum_snap_distance:.9g} mm")
    print(f"Apical permeability: {args.permeability:.9g} mm/s")
    print(f"CSV: {csv_path}")
    if not args.no_chart:
        create_chart(rows)


# ParaView's Macro runner may execute a file with a non-standard ``__name__``.
# This file is an executable macro rather than an importable library, so invoke
# it unconditionally.
main()
