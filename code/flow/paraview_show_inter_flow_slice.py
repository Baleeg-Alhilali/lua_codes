"""ParaView macro: pressure contours and Darcy-flow arrows between two nephrons.

Run this after opening the calibrated Lumen, Membrane, and Inter PVTU files.
The macro finds the closest interior approach of the two lumens and builds a
shared cross-sectional slice there.  It does not depend on a hard-coded camera
position or slice origin.
"""

from paraview.simple import *
from paraview import servermanager
from vtkmodules.vtkCommonDataModel import vtkDataSet, vtkPolyData
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkFiltersCore import vtkAppendFilter, vtkConnectivityFilter
from vtkmodules.vtkCommonDataModel import vtkStaticPointLocator
import numpy as np


PRESSURE_ARRAY = "pressure_kinematic"
VELOCITY_ARRAY = "velocity_mm_s"
NUMBER_OF_CONTOURS = 18
NUMBER_OF_ARROWS = 450
BOUNDARY_EXCLUSION_MM = 0.05
TANGENT_NEIGHBORHOOD_MM = 0.12


def source_with_suffix(suffix):
    """Prefer the calibrated result, then accept any matching open source."""
    candidates = []
    for key, proxy in GetSources().items():
        name = key[0]
        file_name = ""
        if hasattr(proxy, "FileName"):
            value = proxy.FileName
            if isinstance(value, (list, tuple)):
                file_name = " ".join(str(item) for item in value)
            else:
                file_name = str(value)
        searchable_name = name + " " + file_name
        if suffix in searchable_name and "AUTO_" not in name:
            score = 2 if "multi_nephron_50_50" in searchable_name else 1
            candidates.append((score, name, proxy))
    if not candidates:
        raise RuntimeError(
            "Could not find an open source containing '{}'. Open the Lumen, "
            "Membrane, and Inter PVTU files first.".format(suffix)
        )
    candidates.sort(key=lambda item: (item[0], item[1]))
    return candidates[-1][2]


def fetched_dataset(proxy):
    """Fetch and flatten a distributed or composite ParaView dataset."""
    data = servermanager.Fetch(proxy)
    if isinstance(data, vtkDataSet):
        return data

    append = vtkAppendFilter()
    iterator = data.NewIterator()
    iterator.InitTraversal()
    while not iterator.IsDoneWithTraversal():
        block = iterator.GetCurrentDataObject()
        if isinstance(block, vtkDataSet) and block.GetNumberOfCells() > 0:
            append.AddInputData(block)
        iterator.GoToNextItem()
    append.Update()
    return append.GetOutput()


def connected_region_points(dataset):
    """Return point-coordinate arrays for the two largest lumen components."""
    connectivity = vtkConnectivityFilter()
    connectivity.SetInputData(dataset)
    connectivity.SetExtractionModeToAllRegions()
    connectivity.ColorRegionsOn()
    connectivity.Update()
    connected = connectivity.GetOutput()

    region_ids = connected.GetCellData().GetArray("RegionId")
    if region_ids is None:
        raise RuntimeError("Connectivity did not generate a cell RegionId array.")

    point_ids_by_region = {}
    for cell_id in range(connected.GetNumberOfCells()):
        region_id = int(region_ids.GetTuple1(cell_id))
        ids = point_ids_by_region.setdefault(region_id, set())
        cell = connected.GetCell(cell_id)
        for local_id in range(cell.GetNumberOfPoints()):
            ids.add(cell.GetPointId(local_id))

    largest = sorted(
        point_ids_by_region.items(), key=lambda item: len(item[1]), reverse=True
    )[:2]
    if len(largest) != 2:
        raise RuntimeError(
            "Expected two disconnected lumen regions, found {}.".format(
                len(point_ids_by_region)
            )
        )

    result = []
    for _, point_ids in largest:
        result.append(
            np.asarray([connected.GetPoint(point_id) for point_id in point_ids])
        )
    return result


def interior_points(points, bounds):
    """Remove points near OuterWall so terminal caps do not choose the slice."""
    lower = np.asarray([bounds[0], bounds[2], bounds[4]])
    upper = np.asarray([bounds[1], bounds[3], bounds[5]])
    mask = np.all(points > lower + BOUNDARY_EXCLUSION_MM, axis=1)
    mask &= np.all(points < upper - BOUNDARY_EXCLUSION_MM, axis=1)
    selected = points[mask]
    return selected if len(selected) >= 20 else points


def closest_pair(first, second):
    """Find the closest point pair with a VTK point locator."""
    vtk_points = vtkPoints()
    for point in second:
        vtk_points.InsertNextPoint(point)
    poly = vtkPolyData()
    poly.SetPoints(vtk_points)
    locator = vtkStaticPointLocator()
    locator.SetDataSet(poly)
    locator.BuildLocator()

    best_distance_squared = float("inf")
    best_first = None
    best_second = None
    for point in first:
        point_id = locator.FindClosestPoint(point)
        other = second[point_id]
        distance_squared = float(np.dot(point - other, point - other))
        if distance_squared < best_distance_squared:
            best_distance_squared = distance_squared
            best_first = point.copy()
            best_second = other.copy()
    return best_first, best_second


def local_tangent(points, center):
    """Estimate the local tube axis using PCA around a closest point."""
    distances = np.linalg.norm(points - center, axis=1)
    local = points[distances <= TANGENT_NEIGHBORHOOD_MM]
    if len(local) < 20:
        local = points[distances <= 2.0 * TANGENT_NEIGHBORHOOD_MM]
    if len(local) < 3:
        raise RuntimeError("Too few lumen points to estimate the local tangent.")
    centered = local - np.mean(local, axis=0)
    covariance = centered.T.dot(centered) / max(1, len(local) - 1)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    tangent = eigenvectors[:, int(np.argmax(eigenvalues))]
    return tangent / np.linalg.norm(tangent)


def slice_location(lumen_proxy, inter_proxy):
    lumen_data = fetched_dataset(lumen_proxy)
    inter_data = fetched_dataset(inter_proxy)
    regions = connected_region_points(lumen_data)
    first = interior_points(regions[0], inter_data.GetBounds())
    second = interior_points(regions[1], inter_data.GetBounds())
    first_point, second_point = closest_pair(first, second)

    first_tangent = local_tangent(regions[0], first_point)
    second_tangent = local_tangent(regions[1], second_point)
    if np.dot(first_tangent, second_tangent) < 0.0:
        second_tangent *= -1.0
    normal = first_tangent + second_tangent
    if np.linalg.norm(normal) < 1.0e-12:
        normal = first_tangent.copy()

    # Force the connecting line between the nephrons into the slice plane.
    connector = second_point - first_point
    connector /= np.linalg.norm(connector)
    normal -= np.dot(normal, connector) * connector
    normal /= np.linalg.norm(normal)
    origin = 0.5 * (first_point + second_point)
    separation = float(np.linalg.norm(second_point - first_point))
    return origin, normal, separation


def remove_previous_macro_objects():
    for key, proxy in list(GetSources().items()):
        if key[0].startswith("AUTO_"):
            Delete(proxy)


def plane_slice(input_proxy, origin, normal, name):
    result = Slice(registrationName=name, Input=input_proxy)
    result.SliceType = "Plane"
    result.SliceType.Origin = origin.tolist()
    result.SliceType.Normal = normal.tolist()
    return result


remove_previous_macro_objects()

inter_source = source_with_suffix("_inter_")
lumen_source = source_with_suffix("_lumen_")
membrane_source = source_with_suffix("_membrane_")

origin, normal, separation = slice_location(lumen_source, inter_source)

inter_slice = plane_slice(
    inter_source, origin, normal, "AUTO_Inter_closest_slice"
)
lumen_slice = plane_slice(
    lumen_source, origin, normal, "AUTO_Lumen_cross_sections"
)
membrane_slice = plane_slice(
    membrane_source, origin, normal, "AUTO_Membrane_cross_sections"
)
UpdatePipeline(proxy=inter_slice)

pressure_info = inter_slice.GetPointDataInformation().GetArray(PRESSURE_ARRAY)
if pressure_info is None:
    raise RuntimeError("Inter slice has no '{}' point array.".format(PRESSURE_ARRAY))
pressure_min, pressure_max = pressure_info.GetRange()
if pressure_max <= pressure_min:
    pressure_max = pressure_min + 1.0
contour_values = np.linspace(
    pressure_min, pressure_max, NUMBER_OF_CONTOURS + 2
)[1:-1].tolist()

pressure_contours = Contour(
    registrationName="AUTO_Pressure_contours", Input=inter_slice
)
pressure_contours.ContourBy = ["POINTS", PRESSURE_ARRAY]
pressure_contours.Isosurfaces = contour_values

cell_centers = CellCenters(
    registrationName="AUTO_Inter_cell_centers", Input=inter_slice
)
cell_centers.VertexCells = 1

# Project the 3-D Darcy velocity into the slice.  Without this step, vectors
# normal to the slice render as dots rather than useful in-plane arrows.
normal_expression = "({:.17g}*iHat + {:.17g}*jHat + {:.17g}*kHat)".format(
    normal[0], normal[1], normal[2]
)
normal_component = (
    "({0}_X*{1:.17g} + {0}_Y*{2:.17g} + {0}_Z*{3:.17g})"
).format(VELOCITY_ARRAY, normal[0], normal[1], normal[2])
velocity_projection = Calculator(
    registrationName="AUTO_In_plane_velocity", Input=cell_centers
)
velocity_projection.AttributeType = "Point Data"
velocity_projection.ResultArrayName = "velocity_in_slice"
velocity_projection.Function = (
    VELOCITY_ARRAY + " - " + normal_component + "*" + normal_expression
)
UpdatePipeline(proxy=velocity_projection)

velocity_info = velocity_projection.GetPointDataInformation().GetArray(
    "velocity_in_slice"
)
if velocity_info is None:
    raise RuntimeError(
        "Could not construct the in-plane Darcy velocity array."
    )
velocity_max = velocity_info.GetRange(-1)[1]
slice_bounds = inter_slice.GetDataInformation().GetBounds()
slice_diagonal = np.linalg.norm(
    [
        slice_bounds[1] - slice_bounds[0],
        slice_bounds[3] - slice_bounds[2],
        slice_bounds[5] - slice_bounds[4],
    ]
)
arrow_length = 0.022 * slice_diagonal

velocity_arrows = Glyph(
    registrationName="AUTO_Inter_velocity_arrows",
    Input=velocity_projection,
    GlyphType="Arrow",
)
velocity_arrows.OrientationArray = ["POINTS", "velocity_in_slice"]
velocity_arrows.ScaleArray = ["POINTS", "No scale array"]
velocity_arrows.ScaleFactor = arrow_length
velocity_arrows.GlyphMode = "Uniform Spatial Distribution (Bounds Based)"
velocity_arrows.MaximumNumberOfSamplePoints = NUMBER_OF_ARROWS
velocity_arrows.GlyphType.TipResolution = 12
velocity_arrows.GlyphType.ShaftResolution = 12
velocity_arrows.GlyphType.TipRadius = 0.18
velocity_arrows.GlyphType.ShaftRadius = 0.055

view = GetActiveViewOrCreate("RenderView")
for source in (inter_source, lumen_source, membrane_source):
    Hide(source, view)

slice_display = Show(inter_slice, view)
slice_display.Representation = "Surface"
slice_display.Opacity = 0.72
ColorBy(slice_display, ("POINTS", PRESSURE_ARRAY))
slice_display.RescaleTransferFunctionToDataRange(True, False)
pressure_lut = GetColorTransferFunction(PRESSURE_ARRAY)
pressure_lut.ApplyPreset("Cool to Warm", True)
pressure_lut.RescaleTransferFunction(pressure_min, pressure_max)
slice_display.SetScalarBarVisibility(view, True)

contour_display = Show(pressure_contours, view)
contour_display.ColorArrayName = [None, ""]
contour_display.AmbientColor = [0.05, 0.05, 0.05]
contour_display.DiffuseColor = [0.05, 0.05, 0.05]
contour_display.LineWidth = 2.5

arrow_display = Show(velocity_arrows, view)
arrow_display.ColorArrayName = [None, ""]
arrow_display.AmbientColor = [1.0, 0.82, 0.05]
arrow_display.DiffuseColor = [1.0, 0.82, 0.05]
arrow_display.Position = (0.002 * normal).tolist()

membrane_display = Show(membrane_slice, view)
membrane_display.ColorArrayName = [None, ""]
membrane_display.Representation = "Surface"
membrane_display.AmbientColor = [0.95, 0.45, 0.08]
membrane_display.DiffuseColor = [0.95, 0.45, 0.08]

lumen_display = Show(lumen_slice, view)
lumen_display.ColorArrayName = [None, ""]
lumen_display.Representation = "Surface"
lumen_display.AmbientColor = [0.12, 0.12, 0.12]
lumen_display.DiffuseColor = [0.12, 0.12, 0.12]

annotation = Text(registrationName="AUTO_Slice_information")
annotation.Text = (
    "Closest-nephron Inter slice\n"
    "separation = {:.4g} mm\n"
    "black: pressure contours | yellow: Darcy velocity"
).format(separation)
annotation_display = Show(annotation, view)
annotation_display.WindowLocation = "Upper Left Corner"
annotation_display.Color = [0.05, 0.05, 0.05]
annotation_display.FontSize = 12

view.CameraParallelProjection = 1
# Look directly along the slice normal; ResetCamera otherwise preserves the
# previous 3-D camera and can leave the new slice nearly edge-on.
reference_up = np.asarray([0.0, 0.0, 1.0])
if abs(float(np.dot(reference_up, normal))) > 0.9:
    reference_up = np.asarray([0.0, 1.0, 0.0])
view_up = reference_up - np.dot(reference_up, normal) * normal
view_up /= np.linalg.norm(view_up)
view.CameraFocalPoint = origin.tolist()
view.CameraPosition = (origin + normal).tolist()
view.CameraViewUp = view_up.tolist()
view.ResetCamera()
Render()

print("Created closest-nephron Inter flow slice")
print("  origin [mm]:", origin.tolist())
print("  normal:", normal.tolist())
print("  lumen surface separation [mm]:", separation)
print("  pressure range:", (pressure_min, pressure_max))
print("  maximum velocity [mm/s]:", velocity_max)
