/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-12-27
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "neurites_from_swc.h"
#include <functional>

#include "common/math/math_vector_matrix/math_vector_functions.h"  // VecScale
#include "common/util/file_util.h"  // FindFileInStandardPaths
#include "common/util/smart_pointer.h"  // SmartPtr
#include "common/util/string_util.h"  // TrimString etc.
#include "lib_algebra/small_algebra/small_algebra.h" // Invert
#include "lib_disc/domain_util.h"   // LoadDomain
#include "lib_disc/function_spaces/error_elem_marking_strategy.h" // GlobalMarking
#include "lib_disc/quadrature/gauss_legendre/gauss_legendre.h"
#include "lib_grid/algorithms/element_side_util.h" // GetOpposingSide
#include "lib_grid/algorithms/extrusion/extrusion.h" // Extrude
#include "lib_grid/algorithms/geom_obj_util/face_util.h" // CalculateNormal
#include "lib_grid/algorithms/geom_obj_util/edge_util.h" // AdjustEdgeOrientationToFaceOrientation
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h" // RemoveDoubles
#include "lib_grid/algorithms/grid_generation/icosahedron.h" // icosahedron
#include "lib_grid/algorithms/grid_generation/triangle_fill_sweep_line.h" // planar OuterWall ports
#include "lib_grid/algorithms/grid_generation/tetrahedralization.h" // Tetrahedralize
#include "lib_grid/algorithms/orientation_util.h" // FixFaceOrientation, FixOrientation
#include "lib_grid/algorithms/remove_duplicates_util.h" // RemoveDuplicates
#include "lib_grid/algorithms/subset_util.h" // SeparateSubsetsByLowerDimSubsets
#include "lib_grid/file_io/file_io_ugx.h"  // GridWriterUGX
#include "lib_grid/file_io/file_io.h"  // SaveGridHierarchyTransformed
#include "lib_grid/file_io/file_io_swc.h"  // FileReaderSWC
#include "lib_grid/global_attachments.h"
#include "lib_grid/grid/geometry.h" // MakeGeometry3d
#include "lib_grid/grid/neighborhood_util.h"  // for GetConnectedNeighbor
#include "lib_grid/refinement/global_multi_grid_refiner.h" // GlobalMultigridRefiner
#include "lib_grid/refinement/hanging_node_refiner_multi_grid.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"
#include "lib_grid/refinement/projectors/projection_handler.h"
#include "lib_grid/refinement/regular_refinement.h"  // Refine

#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <istream>
#include <map>
#include <limits>
#include <sstream>
#include <queue>
#include <set>
#include <stack>

namespace ug {
namespace neuro_collection {
namespace neurites_from_swc {

struct OuterWallTerminalMetadata
{
	bool terminalsOnOuterWall;
	bool hasFixedBox;
	vector3 boxMin;
	vector3 boxMax;
	OuterWallTerminalMetadata() : terminalsOnOuterWall(false), hasFixedBox(false) {}
};


struct AdaptiveNephronMeshingControls
{
	bool curvatureLimitedAnisotropy;
	number minimumBendAnisotropy;
	number bendChordErrorFactor;
	bool gradedInter;
	number interNearEdgeFactor;
	number interNearDistanceFactor;
	number interFarDistanceFactor;
	number interFarEdgeFactor;
	bool prismaticInter;
	AdaptiveNephronMeshingControls()
		: curvatureLimitedAnisotropy(false), minimumBendAnisotropy(1.0),
		  bendChordErrorFactor(0.02), gradedInter(false),
		  interNearEdgeFactor(2.0), interNearDistanceFactor(1.5),
		  interFarDistanceFactor(6.0), interFarEdgeFactor(6.0),
		  prismaticInter(false) {}
};

static AdaptiveNephronMeshingControls g_adaptiveNephronMeshing;


void configure_adaptive_nephron_meshing
(
	bool curvatureLimitedAnisotropy,
	number minimumBendAnisotropy,
	number bendChordErrorFactor,
	bool gradedInter,
	number interNearEdgeFactor,
	number interNearDistanceFactor,
	number interFarDistanceFactor,
	number interFarEdgeFactor,
	bool prismaticInter
)
{
	UG_COND_THROW(minimumBendAnisotropy < 1.0,
	              "Minimum bend anisotropy must be at least 1.");
	UG_COND_THROW(bendChordErrorFactor <= 0.0,
	              "Bend chord-error factor must be positive.");
	UG_COND_THROW(interNearEdgeFactor <= 0.0 || interFarEdgeFactor < interNearEdgeFactor,
	              "Inter edge factors must be positive and far >= near.");
	UG_COND_THROW(interNearDistanceFactor < 0.0
	              || interFarDistanceFactor <= interNearDistanceFactor,
	              "Inter distance factors require 0 <= near < far.");
	g_adaptiveNephronMeshing.curvatureLimitedAnisotropy = curvatureLimitedAnisotropy;
	g_adaptiveNephronMeshing.minimumBendAnisotropy = minimumBendAnisotropy;
	g_adaptiveNephronMeshing.bendChordErrorFactor = bendChordErrorFactor;
	g_adaptiveNephronMeshing.gradedInter = gradedInter;
	g_adaptiveNephronMeshing.interNearEdgeFactor = interNearEdgeFactor;
	g_adaptiveNephronMeshing.interNearDistanceFactor = interNearDistanceFactor;
	g_adaptiveNephronMeshing.interFarDistanceFactor = interFarDistanceFactor;
	g_adaptiveNephronMeshing.interFarEdgeFactor = interFarEdgeFactor;
	g_adaptiveNephronMeshing.prismaticInter = prismaticInter;
	UG_LOGN("Adaptive nephron meshing: curvatureLimitedAnisotropy="
	        << curvatureLimitedAnisotropy << ", minimumBendAnisotropy="
	        << minimumBendAnisotropy << ", bendChordErrorFactor="
	        << bendChordErrorFactor << ", gradedInter=" << gradedInter
	        << ", Inter edge factors=" << interNearEdgeFactor << " -> "
	        << interFarEdgeFactor << ", distance factors="
	        << interNearDistanceFactor << " -> " << interFarDistanceFactor
	        << ", Inter topology=" << (prismaticInter ? "prism" : "tet")
	        << ".");
}


static OuterWallTerminalMetadata read_outer_wall_terminal_metadata
(
	const std::string& fileName
)
{
	OuterWallTerminalMetadata metadata;
	std::ifstream input(fileName.c_str());
	if (!input) return metadata;
	std::string line;
	while (std::getline(input, line))
	{
		const std::string terminalMarker = "# UG4_TERMINALS_ON_OUTER_WALL";
		const std::string boxMarker = "# UG4_OUTER_BOX";
		if (line.find(terminalMarker) == 0)
		{
			std::istringstream values(line.substr(terminalMarker.size()));
			int enabled = 0;
			values >> enabled;
			metadata.terminalsOnOuterWall = enabled != 0;
		}
		else if (line.find(boxMarker) == 0)
		{
			std::istringstream values(line.substr(boxMarker.size()));
			values >> metadata.boxMin[0] >> metadata.boxMax[0]
			       >> metadata.boxMin[1] >> metadata.boxMax[1]
			       >> metadata.boxMin[2] >> metadata.boxMax[2];
			metadata.hasFixedBox = !values.fail();
		}
		else if (!line.empty() && line[0] != '#')
			break;
	}
	UG_COND_THROW(metadata.terminalsOnOuterWall && !metadata.hasFixedBox,
	              "Terminal-to-OuterWall SWC metadata is missing UG4_OUTER_BOX bounds.");
	if (metadata.hasFixedBox)
		for (size_t d = 0; d < 3; ++d)
			UG_COND_THROW(!(metadata.boxMin[d] < metadata.boxMax[d]),
			              "Invalid fixed OuterWall bounds in SWC metadata.");
	return metadata;
}


template <class TSubsetHandler>
static void snap_refined_terminal_vertices_to_outer_wall
(
	Grid& g,
	TSubsetHandler& sh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >&
		aaSurfParams,
	const OuterWallTerminalMetadata& metadata,
	size_t nephronCount,
	size_t refinementLevel
)
{
	if (!metadata.terminalsOnOuterWall) return;
	const int outerWallSubset = (int)(6 * nephronCount);
	const int tripleJunctionBase = outerWallSubset + 2;
	const int membraneWallBase = tripleJunctionBase + (int)nephronCount;
	number diagonal = 0.0;
	for (size_t d = 0; d < 3; ++d)
		diagonal += (metadata.boxMax[d] - metadata.boxMin[d])
			* (metadata.boxMax[d] - metadata.boxMin[d]);
	diagonal = std::sqrt(diagonal);
	const number numericalLengthTolerance = std::numeric_limits<number>::epsilon()
		* diagonal * 64.0;
	const number maximumSnapDistance = std::max(
		numericalLengthTolerance, diagonal * 2e-2);
	size_t snapped = 0;
	for (VertexIterator vit = g.begin<Vertex>(); vit != g.end<Vertex>(); ++vit)
	{
		Vertex* vertex = *vit;
		const int si = sh.get_subset_index(vertex);
		const bool capVertex = si >= 0 && si < outerWallSubset
			&& (si % 6 == 4 || si % 6 == 5);
		const bool membraneWallVertex = si >= membraneWallBase
			&& si < membraneWallBase + (int)nephronCount;
		const bool rimVertex = si >= tripleJunctionBase
			&& si < tripleJunctionBase + (int)nephronCount;
		// The terminal membrane annulus belongs to the Membrane subset rather
		// than Inlet/Outlet.  Its refined children are nevertheless on the same
		// OuterWall plane.  NeuriteProjector otherwise moves them back toward
		// the (slightly unsnapped) fitted spline endpoint, which can collapse
		// the first membrane cells for coarse angular O-grids.
		const bool nephronVertex = si >= 0 && si < outerWallSubset;
		const bool projectorVertex = nephronVertex || rimVertex || membraneWallVertex;
		const number axial = projectorVertex ? aaSurfParams[vertex].axial : 0.5;
		const number radial = projectorVertex ? aaSurfParams[vertex].radial : 0.0;
		const bool terminalParamVertex = projectorVertex
			&& radial > 1e-10
			&& (std::fabs(axial) < 1e-10 || std::fabs(axial - 1.0) < 1e-10);
		if (!capVertex && !rimVertex && !membraneWallVertex &&
		    !terminalParamVertex) continue;
		number nearestDistance = std::numeric_limits<number>::max();
		size_t nearestAxis = 0;
		size_t nearestSide = 0;
		for (size_t axis = 0; axis < 3; ++axis)
			for (size_t side = 0; side < 2; ++side)
			{
				const number plane = side ? metadata.boxMax[axis] : metadata.boxMin[axis];
				const number distance = std::fabs(aaPos[vertex][axis] - plane);
				if (distance < nearestDistance)
				{
					nearestDistance = distance;
					nearestAxis = axis;
					nearestSide = side;
				}
			}
		// A newly created interior vertex may temporarily carry the attachment's
		// default axial value (zero).  Treat it as terminal only when it is also
		// geometrically close to an OuterWall plane.  Named cap/rim subsets remain
		// strict because every one of those vertices must lie on that plane.
		if (nearestDistance > maximumSnapDistance)
		{
			UG_COND_THROW(capVertex || rimVertex || membraneWallVertex,
			              "Refined terminal vertex is too far from every OuterWall plane.");
			continue;
		}
		aaPos[vertex][nearestAxis] = nearestSide
			? metadata.boxMax[nearestAxis] : metadata.boxMin[nearestAxis];
		++snapped;
	}
	UG_LOGN("Refinement level " << refinementLevel << ": snapped " << snapped
	        << " Inlet/Outlet/TripleJunction/MembraneWall vertices to "
	        << "OuterWall planes.");
}


// A constraining Basolateral quadrilateral and its two constrained Inter
// triangles each create a level copy of their shared coarse corner.  UG4 passes
// npSurfParams to the copy made on the quadrilateral side, but the coincident
// copy made on the constrained-triangle side can retain the attachment default
// (radial == 0).  At the next refinement that false zero is averaged with the
// membrane radial coordinate and pulls a child vertex toward the centerline.
// Restore the parameters from the geometrically identical projector vertex
// before another level is marked/refined.
static size_t restore_coincident_neurite_surface_params
(
	MultiGrid& grid,
	ISubsetHandler& sh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >&
		aaSurfParams,
	size_t nephronCount,
	size_t refinementLevel
)
{
	typedef std::array<long long,3> PositionKey;
	vector3 gridMin(std::numeric_limits<number>::max());
	vector3 gridMax(-std::numeric_limits<number>::max());
	for (VertexIterator it = grid.begin<Vertex>(); it != grid.end<Vertex>(); ++it)
		for (size_t d = 0; d < 3; ++d)
		{
			gridMin[d] = std::min(gridMin[d], aaPos[*it][d]);
			gridMax[d] = std::max(gridMax[d], aaPos[*it][d]);
		}
	const number gridDiagonal = VecDistance(gridMin, gridMax);
	const number coincidenceTolerance = std::max(
		std::numeric_limits<number>::epsilon() * gridDiagonal * 64.0,
		gridDiagonal * 1e-11);
	const number keyScale = 1.0 / coincidenceTolerance;
	auto key = [&] (Vertex* vertex) -> PositionKey
	{
		PositionKey result = {{
			(long long)std::llround(aaPos[vertex][0] * keyScale),
			(long long)std::llround(aaPos[vertex][1] * keyScale),
			(long long)std::llround(aaPos[vertex][2] * keyScale)}};
		return result;
	};
	std::map<PositionKey,Vertex*> validAtPosition;
	for (VertexIterator it = grid.begin<Vertex>(); it != grid.end<Vertex>(); ++it)
		if (aaSurfParams[*it].radial > 1e-10)
			validAtPosition[key(*it)] = *it;

	size_t restored = 0;
	for (VertexIterator it = grid.begin<Vertex>(); it != grid.end<Vertex>(); ++it)
	{
		Vertex* vertex = *it;
		if (grid.get_level(vertex) != refinementLevel) continue;
		if (std::fabs((number)aaSurfParams[vertex].radial) > 1e-10) continue;
		std::map<PositionKey,Vertex*>::const_iterator found =
			validAtPosition.find(key(vertex));
		if (found == validAtPosition.end()) continue;
		vector3 difference;
		VecSubtract(difference, aaPos[vertex], aaPos[found->second]);
		if (VecLengthSq(difference) > coincidenceTolerance*coincidenceTolerance)
			continue;
		aaSurfParams[vertex] = aaSurfParams[found->second];
		++restored;
	}
	UG_LOGN("Refinement level " << refinementLevel << ": restored "
	        << restored << " coincident Basolateral projector parameters.");
	return restored;
}

// Build the same deterministic rotation-minimizing frame that is used by
// NeuriteProjector during refinement.  Keeping generation and projection on
// the identical frame field prevents newly projected midpoints from rotating
// away from the two rings that created their parent edge/face.
static void compute_nephron_parallel_transport_ONB
(
	vector3& tangentOut,
	vector3& secondOut,
	vector3& thirdOut,
	const NeuriteProjector::Neurite& neurite,
	number t
)
{
	const std::vector<NeuriteProjector::Section>& sections = neurite.vSec;
	UG_COND_THROW(sections.empty(), "Cannot transport a frame on an empty neurite.");

	auto velocityAt = [&] (number axial, vector3& velocity)
	{
		NeuriteProjector::Section cmp(axial);
		std::vector<NeuriteProjector::Section>::const_iterator secIt =
			std::lower_bound(sections.begin(), sections.end(), cmp,
			                 NeuriteProjector::CompareSections());
		if (secIt == sections.end()) secIt = sections.end() - 1;
		const number m = secIt->endParam - axial;
		const number* sx = &secIt->splineParamsX[0];
		const number* sy = &secIt->splineParamsY[0];
		const number* sz = &secIt->splineParamsZ[0];
		velocity[0] = (-3.0*sx[0]*m - 2.0*sx[1])*m - sx[2];
		velocity[1] = (-3.0*sy[0]*m - 2.0*sy[1])*m - sy[2];
		velocity[2] = (-3.0*sz[0]*m - 2.0*sz[1])*m - sz[2];
		VecNormalize(velocity, velocity);
	};

	vector3 tangent;
	velocityAt(0.0, tangent);
	number fac = VecProd(neurite.refDir, tangent);
	VecScaleAdd(secondOut, 1.0, neurite.refDir, -fac, tangent);
	if (VecNormSquared(secondOut) < 1e-20)
	{
		vector3 fallback(0.0);
		fallback[std::fabs(tangent[0]) < 0.8 ? 0 : 1] = 1.0;
		VecScaleAdd(secondOut, 1.0, fallback,
		            -VecProd(fallback, tangent), tangent);
	}
	VecNormalize(secondOut, secondOut);
	VecCross(thirdOut, tangent, secondOut);

	const number clampedT = std::max<number>(0.0, std::min<number>(1.0, t));
	const size_t fixedSteps = 128;
	const size_t completeSteps = (size_t)std::floor(fixedSteps * clampedT);
	for (size_t i = 1; i <= completeSteps; ++i)
	{
		velocityAt((number)i / (number)fixedSteps, tangent);
		vector3 transported;
		VecScaleAdd(transported, 1.0, secondOut,
		            -VecProd(secondOut, tangent), tangent);
		if (VecNormSquared(transported) < 1e-20)
		{
			VecCross(transported, thirdOut, tangent);
			if (VecNormSquared(transported) < 1e-20)
			{
				vector3 fallback(0.0);
				fallback[std::fabs(tangent[0]) < 0.8 ? 0 : 1] = 1.0;
				VecScaleAdd(transported, 1.0, fallback,
				            -VecProd(fallback, tangent), tangent);
			}
		}
		VecNormalize(secondOut, transported);
		VecCross(thirdOut, tangent, secondOut);
	}
	const number lastGridPoint = (number)completeSteps / (number)fixedSteps;
	if (clampedT > lastGridPoint + 1e-14)
	{
		velocityAt(clampedT, tangent);
		vector3 transported;
		VecScaleAdd(transported, 1.0, secondOut,
		            -VecProd(secondOut, tangent), tangent);
		if (VecNormSquared(transported) < 1e-20)
			VecCross(transported, thirdOut, tangent);
		VecNormalize(secondOut, transported);
		VecCross(thirdOut, tangent, secondOut);
	}
	tangentOut = tangent;
}


#if 0
static void smoothing(std::vector<swc_types::SWCPoint>& vPointsInOut, size_t n, number h, number gamma)
{
	// find neurite root vertices
	const size_t nP = vPointsInOut.size();
	std::vector<size_t> rootVrts;
	std::vector<bool> treated(nP, false);
	for (size_t i = 0; i < nP; ++i)
	{
		if (treated[i]) continue;
		treated[i] = true;
		if (vPointsInOut[i].type != swc_types::SWC_SOMA) continue;

		// here, we have a soma point;
		// find first non-soma point in all directions
		std::queue<size_t> q;
		q.push(i);

		while (!q.empty())
		{
			const size_t ind = q.front();
			const swc_types::SWCPoint& pt = vPointsInOut[ind];
			q.pop();

			if (pt.type == swc_types::SWC_SOMA)
			{
				const size_t nConn = pt.conns.size();
				for (size_t j = 0; j < nConn; ++j)
					if (!treated[pt.conns[j]])
						q.push(pt.conns[j]);
			}
			else
				rootVrts.push_back(ind);

			treated[ind] = true;
		}
	}

	// starting at root vertices, smooth the entire tree(s),
	// but leave out soma vertices as well as branching points
	std::vector<vector3> newPos(nP);
	for (size_t i = 0; i < n; ++i)
	{
		treated.clear();
		treated.resize(nP, false);

		std::stack<size_t> stack;
		for (size_t rv = 0; rv < rootVrts.size(); ++rv)
			stack.push(rootVrts[rv]);

		while (!stack.empty())
		{
			size_t ind = stack.top();
			stack.pop();
			const swc_types::SWCPoint& pt = vPointsInOut[ind];
			const vector3& x = pt.coords;
			vector3& x_new = newPos[ind];

			UG_COND_THROW(treated[ind], "Circle detected in supposedly tree-shaped neuron!\n"
				"Position: " << vPointsInOut[ind].coords);
			treated[ind] = true;

			// somata are not smoothed and not iterated over
			if (pt.type == swc_types::SWC_SOMA)
			{
				x_new = x;
				continue;
			}

			// branching points are not smoothed, but iterated over
			size_t connSz = pt.conns.size();
			for (size_t c = 0; c < connSz; ++c)
				if (!treated[pt.conns[c]])
					stack.push(pt.conns[c]);

			if (connSz != 2)
			{
				x_new = x;
				continue;
			}

			// here we have a non-branching, non-end, non-soma point: smooth
			const vector3& x1 = vPointsInOut[pt.conns[0]].coords;
			const vector3& x2 = vPointsInOut[pt.conns[1]].coords;

			number d1 = VecDistanceSq(x1, x);
			number d2 = VecDistanceSq(x2, x);
			number w1 = std::exp(-d1/(h*h));
			number w2 = std::exp(-d2/(h*h));

			// only really smooth if both adjacent edges are short
			number w = std::min(w1, w2);

			// correction
			vector3 corr;
			VecScaleAdd(corr, w, x1, -2*w, x, w, x2);
			VecScale(corr, corr, 1.0 / (1.0 + 2*w));

			// take only the part orthogonal to x1 - x2,
			// we do not want to shift x towards the nearer neighbor
			VecSubtract(x_new, x1, x2); // using x_new as intermediate variable
			number normSq = VecNormSquared(x_new);
			VecScaleAdd(corr, 1.0, corr, - VecProd(corr, x_new) / normSq, x_new);
			VecScaleAdd(x_new, 1.0, x, gamma, corr);
		}

		// assign new positions
		for (size_t p = 0; p < nP; ++p)
			if (treated[p]) // soma points may not have been treated
				vPointsInOut[p].coords = newPos[p];
	}

}



struct EdgeLengthCompare
{
	bool operator()(const std::pair<Edge*, number> e1, const std::pair<Edge*, number> e2)
	{return e1.second > e2.second;}
};


static void collapse_short_edges(Grid& g, SubsetHandler& sh)
{
	// get access to positions
	UG_COND_THROW(!g.has_vertex_attachment(aPosition), "Position attachment not attached to grid.")
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

	// get access to diameter attachment
	ANumber aDiam = GlobalAttachments::attachment<ANumber>("diameter");
	UG_COND_THROW(!g.has_vertex_attachment(aDiam), "No diameter attachment attached to grid.");
	Grid::AttachmentAccessor<Vertex, ANumber> aaDiam(g, aDiam);

	// a short edge is one that is shorter than its diameter
	// sort all short edges in a priority queue
	std::priority_queue<std::pair<Edge*, number>, std::vector<std::pair<Edge*, number> >, EdgeLengthCompare> pq;
	EdgeIterator eit = g.begin<Edge>();
	EdgeIterator edge_end = g.end<Edge>();
	for (; eit != edge_end; ++eit)
	{
		Edge* e = *eit;
		number length = EdgeLengthSq(e, aaPos);
		number diam = std::max(aaDiam[e->vertex(0)], aaDiam[e->vertex(1)]);
		diam = diam*diam;
		if (length < diam)
			pq.push(std::make_pair(e, length));
	}

	while (!pq.empty())
	{
		std::pair<Edge*, number> elp = pq.top();
		pq.pop();

		// edge length might not be up to date; if so: re-insert with correct length (if necessary)
		Edge* curEdge = elp.first;
		number curLen = EdgeLengthSq(curEdge, aaPos);
		if (curLen != elp.second)
		{
			number curDiam = std::max(aaDiam[curEdge->vertex(0)], aaDiam[curEdge->vertex(1)]);
			curDiam = curDiam*curDiam;
			if (curLen < curDiam)
			{
				elp.second = curLen;
				pq.push(elp);
			}
			continue;
		}

		// do not consider short edges connecting two branching points
		Vertex* v1 = curEdge->vertex(0);
		Vertex* v2 = curEdge->vertex(1);
		Grid::AssociatedEdgeIterator it = g.associated_edges_begin(v1);
		Grid::AssociatedEdgeIterator it_end = g.associated_edges_end(v1);
		size_t nAssV1 = 0;
		for (; it != it_end; ++it)
			++nAssV1;

		it = g.associated_edges_begin(v2);
		it_end = g.associated_edges_end(v2);
		size_t nAssV2 = 0;
		for (; it != it_end; ++it)
			++nAssV2;

		if (nAssV1 > 2 && nAssV2 > 2)
			continue;

		// otherwise, collapse edge

		// (a) calculate position (and radius) for new vertex --
		// if the old edge was parallel to one of the adjacent edges,
		// then add the complete edge length to that adjacent edge
		number newDiam;
		vector3 newPos;
		vector3 x1 = aaPos[v1];
		vector3 x2 = aaPos[v2];
		vector3 d0, d1, d2;
		VecSubtract(d0, x2, x1);
		VecNormalize(d0, d0);

		// never move branching points
		if (nAssV1 > 2)
		{
			newDiam = aaDiam[v1];
			newPos = x1;
		}
		else if (nAssV2 > 2)
		{
			newDiam = aaDiam[v2];
			newPos = x2;
		}

		// if the edge is a terminal edge, set new vertex at previous terminal vertex
		else if (nAssV1 == 1)
		{
			newDiam = aaDiam[v1];
			newPos = x1;
		}
		else if (nAssV2 == 1)
		{
			newDiam = aaDiam[v2];
			newPos = x2;
		}
		else
		{
			// calculate directions of adjacent edges
			it = g.associated_edges_begin(v1);
			if (*it != curEdge)
				VecSubtract(d1, x1, aaPos[GetOpposingSide(g, *it, v1)]);
			else
				VecSubtract(d1, x1, aaPos[GetOpposingSide(g, *(++it), v1)]);

			it = g.associated_edges_begin(v2);
			if (*it != curEdge)
				VecSubtract(d2, aaPos[GetOpposingSide(g, *it, v2)], x2);
			else
				VecSubtract(d2, aaPos[GetOpposingSide(g, *(++it), v2)], x2);

			VecNormalize(d1, d1);
			VecNormalize(d2, d2);
			number w1 = 1.0 - fabs(VecProd(d0, d1));
			number w2 = 1.0 - fabs(VecProd(d0, d2));

			// if all three directions are practically co-linear, choose middle
			if (w1 < 0.05 && w2 < 0.05)	// corresponds to a deviation of about 18 degrees
			{
				newDiam = 0.5 * (aaDiam[v1] + aaDiam[v2]);
				VecScaleAdd(newPos, 0.5, x1, 0.5, x2);
			}
			// otherwise, weighted sum
			else
			{
				newDiam = w1 * aaDiam[v1] + w2 * aaDiam[v2];
				newDiam /= (w1 + w2);
				VecScaleAdd(newPos, w1, x1, w2, x2);
				VecScale(newPos, newPos, 1.0 / (w1 + w2));
			}
		}

		// (b) actual collapse
		Vertex* newVrt = *g.create<RegularVertex>();
		sh.assign_subset(newVrt, sh.get_subset_index(curEdge));
		CollapseEdge(g, curEdge, newVrt);

		// (c) assign the new vertex its position and diameter
		aaPos[newVrt] = newPos;
		aaDiam[newVrt] = newDiam;
	}
}
#endif



static void convert_pointlist_to_neuritelist
(
    const std::vector<swc_types::SWCPoint>& vPoints,
    std::vector<std::vector<vector3> >& vPosOut,
    std::vector<std::vector<number> >& vRadOut,
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > >& vBPInfoOut,
    std::vector<size_t>& vRootNeuriteIndsOut
)
{
    // clear out vectors
    vPosOut.clear();
    vRadOut.clear();
    vBPInfoOut.clear();
    vRootNeuriteIndsOut.clear();

	size_t nPts = vPoints.size();
	std::vector<bool> ptProcessed(nPts, false);
	size_t nProcessed = 0;
	size_t curNeuriteInd = 0;

	while (nProcessed != nPts)
	{
		// find first soma's root point in geometry and save its index as i
		size_t i = 0;
		for (; i < nPts; ++i)
		{
			if (vPoints[i].type == swc_types::SWC_SOMA && !ptProcessed[i])
				break;
		}
		UG_COND_THROW(i == nPts, "No soma contained in (non-empty) list of unprocessed SWC points, \n"
				"i.e., there is at least one SWC point not connected to any soma.");

		// collect neurite root points
		std::vector<std::pair<size_t, size_t> > rootPts;
		std::queue<std::pair<size_t, size_t> > soma_queue;
		soma_queue.push(std::make_pair((size_t)-1, i));
		while (!soma_queue.empty())
		{
			size_t pind = soma_queue.front().first;
			size_t ind = soma_queue.front().second;
			soma_queue.pop();

			const swc_types::SWCPoint& pt = vPoints[ind];

			if (pt.type == swc_types::SWC_SOMA)
			{
				ptProcessed[ind] = true;
				++nProcessed;

				size_t nConn = pt.conns.size();
				for (size_t j = 0; j < nConn; ++j)
					if (pt.conns[j] != pind)
						soma_queue.push(std::make_pair(ind, pt.conns[j]));
			}
			else
				rootPts.push_back(std::make_pair(pind, ind));
		}

		vPosOut.resize(vPosOut.size() + rootPts.size());
		vRadOut.resize(vRadOut.size() + rootPts.size());
		vBPInfoOut.resize(vBPInfoOut.size() + rootPts.size());

		std::stack<std::pair<size_t, size_t> > processing_stack;
		for (size_t i = 0; i < rootPts.size(); ++i)
			processing_stack.push(rootPts[i]);

		vRootNeuriteIndsOut.push_back(curNeuriteInd);

		// helper map to be used to correctly save BPs:
		// maps branch root parent ID to neurite ID and BP ID
		std::map<size_t, std::pair<size_t, size_t> > helperMap;

		while (!processing_stack.empty())
		{
			size_t pind = processing_stack.top().first;
			size_t ind = processing_stack.top().second;
			processing_stack.pop();

			ptProcessed[ind] = true;
			++nProcessed;

			const swc_types::SWCPoint& pt = vPoints[ind];

			UG_COND_THROW(pt.type == swc_types::SWC_SOMA, "Detected neuron with more than one soma.");

			// push back coords and radius information to proper neurite
			vPosOut[curNeuriteInd].push_back(pt.coords);
			vRadOut[curNeuriteInd].push_back(pt.radius);

			size_t nConn = pt.conns.size();

			// branching point
			if (nConn > 2)
			{
				// branch with minimal angle will continue current branch
				vector3 parentDir;
				VecSubtract(parentDir, pt.coords, vPoints[pind].coords);
				VecNormalize(parentDir, parentDir);

				size_t parentToBeDiscarded = 0;
				size_t minAngleInd = 0;
				number minAngle = std::numeric_limits<number>::infinity();

				for (size_t i = 0; i < nConn; ++i)
				{
					if (pt.conns[i] == pind)
					{
						parentToBeDiscarded = i;
						continue;
					}

					vector3 dir;
					VecSubtract(dir, vPoints[pt.conns[i]].coords, pt.coords);
					VecNormalize(dir, dir);

					number angle = acos(VecProd(dir, parentDir));
					if (angle < minAngle)
					{
						minAngle = angle;
						minAngleInd = i;
					}
				}

				// branching point info
				std::pair<size_t, std::vector<size_t> > bp;
				bp.first = vPosOut[curNeuriteInd].size()-1; // BP is located at this index of the neurite

				// resize out vectors to accommodate new neurites starting here
				size_t newSize = vPosOut.size() + nConn - 2;
				vPosOut.resize(newSize);
				vRadOut.resize(newSize);
				vBPInfoOut.resize(newSize);

				for (size_t i = 0; i < nConn; ++i)
				{
					if (i == parentToBeDiscarded || i == minAngleInd)
						continue;

					// push new neurite starting point index to stack
					processing_stack.push(std::make_pair(ind, pt.conns[i]));

					// save current point ID with current neurite ID and BP ID
					// for later assignment of branching neurite ID to BP
					// NOT: bp.second.push_back(++futureNeuriteInd);
					helperMap[ind] = std::make_pair(curNeuriteInd, vBPInfoOut[curNeuriteInd].size());
				}

				// push next index of the current neurite to stack
				processing_stack.push(std::make_pair(ind, pt.conns[minAngleInd]));

				// add BP to BP vector
				vBPInfoOut[curNeuriteInd].push_back(bp);
			}

			// end point
			else if (nConn == 1)
			{
				// if the stack is not empty, the next ID on it will start a new neurite
				if (!processing_stack.empty())
				{
					++curNeuriteInd;

					size_t nextParentID = processing_stack.top().first;

					// is this a root neurite? (this is the case if helper map does not contain next parent)
					std::map<size_t, std::pair<size_t, size_t> >::const_iterator it = helperMap.find(nextParentID);
					if (it != helperMap.end())
					{
						// push back parent position and radius to new neurite
						vPosOut[curNeuriteInd].push_back(vPoints[nextParentID].coords);
						vRadOut[curNeuriteInd].push_back(vPoints[nextParentID].radius);

						// push back new neurite ID to BP at parent
						size_t nextParentNeuriteID = it->second.first;
						size_t nextParentBPID = it->second.second;
						vBPInfoOut[nextParentNeuriteID][nextParentBPID].second.push_back(curNeuriteInd);
					}
					// else: the next point is the root point of a root neurite
					else
					{
						vRootNeuriteIndsOut.push_back(curNeuriteInd);
					}
				}
			}

			// normal point
			else
			{
				for (size_t i = 0; i < nConn; ++i)
				{
					if (pt.conns[i] != pind)
						processing_stack.push(std::make_pair(ind, pt.conns[i]));
				}
			}
		}

		// last neurite of each neuron does not increase counter
		++curNeuriteInd;
	}

#if 0
    size_t numSomaPoints = vSomaPoints.size();
    UG_LOGN("Number of soma points: " << numSomaPoints);
    for (size_t i = 0; i < numSomaPoints; i++) {
    	UG_LOGN("Coords for soma: " << vSomaPoints[i].coords);
    }
#endif
}



static void create_spline_data_for_neurites
(
    std::vector<NeuriteProjector::Neurite>& vNeuritesOut,
    const std::vector<std::vector<vector3> >& vPos,
    const std::vector<std::vector<number> >& vR,
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > >* vBPInfo = NULL
)
{
    size_t nNeurites = vPos.size();
    vNeuritesOut.resize(nNeurites);

    std::vector<vector3> parentDirections(nNeurites);

    // first: reserve memory for branching region vectors
    // (we will point to their elements in BranchingPoints and do not want the vectors to reallocate!)
    if (vBPInfo)
        for (size_t n = 0; n < nNeurites; ++n)
            vNeuritesOut[n].vBR.reserve((*vBPInfo)[n].size()+1);

    for (size_t n = 0; n < nNeurites; ++n)
    {
        NeuriteProjector::Neurite& neuriteOut = vNeuritesOut[n];
        const std::vector<vector3>& pos = vPos[n];
        std::vector<number> r = vR[n];
        std::vector<std::pair<size_t, std::vector<size_t> > >* bpInfo = NULL;

        if (vBPInfo)
             bpInfo = &((*vBPInfo)[n]);

        // parameterize to achieve constant velocity on piece-wise linear geom
        size_t nVrt = pos.size();
        std::vector<number> tSuppPos(nVrt);
        std::vector<number> dt(nVrt);
        number totalLength = 0.0;
        for (size_t i = 0; i < nVrt-1; ++i)
        {
            tSuppPos[i] = totalLength;
            totalLength += VecDistance(pos[i], pos[i+1]);
        }
        for (size_t i = 0; i < nVrt-1; ++i)
            tSuppPos[i] /= totalLength;
        tSuppPos[nVrt-1] = 1.0;

        for (size_t i = 0; i < nVrt-1; ++i)
            dt[i+1] = tSuppPos[i+1] - tSuppPos[i];


        // now calculate cubic splines in all dimensions and for radius
        // using moments (as in script: Wittum, Numerik 0)
        DenseMatrix<VariableArray2<number> > mat;
        DenseVector<VariableArray1<number> > x0, x1, x2, xr, rhs;
        mat.resize(nVrt, nVrt);
        x1.resize(nVrt);
        rhs.resize(nVrt);

        for (size_t i = 0; i < nVrt; ++i)
            mat(i,i) = 2.0;
        for (size_t i = 1; i < nVrt-1; ++i)
        {
            number h2 = tSuppPos[i+1] - tSuppPos[i-1];
            mat(i,i+1) = dt[i+1] / h2;
            mat(i,i-1) = dt[i] / h2;
        }

        // boundary conditions (matrix)
        bool isRootNeurite = neuriteOut.vBR.size() == 0;
    	vector3 orthogStartDir(0.0);
    	if (!isRootNeurite)
    	{
    		// prepare boundary conditions
			vector3 startDir;
			VecSubtract(startDir, pos[1], pos[0]);
			const number compInParentDir = VecDot(parentDirections[n], startDir);
			VecScaleAdd(orthogStartDir, 1.0, startDir, -compInParentDir, parentDirections[n]);
			VecNormalize(orthogStartDir, orthogStartDir);
    	}

        if (!isRootNeurite)
        {
        	// add suitable flux boundary conditions to ensure orthogonal branching
        	mat(0,1) = 1.0;
        	mat(nVrt-1, nVrt-2) = 1.0;
        }

        // invert matrix
        UG_COND_THROW(!Invert(mat), "Failed to invert moment matrix for spline calculation.")

        // x coord
        for (size_t i = 1; i < nVrt-1; ++i)
            rhs[i] = 6.0 / (tSuppPos[i+1] - tSuppPos[i-1]) *
                     ((pos[i+1][0] - pos[i][0]) / dt[i+1]
                     - (pos[i][0] - pos[i-1][0]) / dt[i]);
        rhs[0] = 0.0;
        rhs[nVrt-1] = 0.0;
        if (!isRootNeurite)
        {
        	// add suitable flux boundary conditions to ensure orthogonal branching
        	number ddx0 = orthogStartDir[0] * totalLength; // deriv w.r.t. x in pt. 0
        	rhs[0] = 6.0 / dt[1] * ((pos[1][0] - pos[0][0]) / dt[1] - ddx0);
        }
        x0 = mat*rhs;

        // y coord
        for (size_t i = 1; i < nVrt-1; ++i)
            rhs[i] = 6.0 / (tSuppPos[i+1] - tSuppPos[i-1]) *
                     ((pos[i+1][1] - pos[i][1]) / dt[i+1]
                     - (pos[i][1] - pos[i-1][1]) / dt[i]);
        if (!isRootNeurite)
        {
        	// add suitable flux boundary conditions to ensure orthogonal branching
        	number ddy0 = orthogStartDir[1] * totalLength; // deriv w.r.t. y in pt. 0
        	rhs[0] = 6.0 / dt[1] * ((pos[1][1] - pos[0][1]) / dt[1] - ddy0);
        }
        x1 = mat*rhs;

        // z coord
        for (size_t i = 1; i < nVrt-1; ++i)
            rhs[i] = 6.0 / (tSuppPos[i+1] - tSuppPos[i-1]) *
                     ((pos[i+1][2] - pos[i][2]) / dt[i+1]
                     - (pos[i][2] - pos[i-1][2]) / dt[i]);
        if (!isRootNeurite)
		{
			// add suitable flux boundary conditions to ensure orthogonal branching
			number ddz0 = orthogStartDir[2] * totalLength; // deriv w.r.t. z in pt. 0
			rhs[0] = 6.0 / dt[1] * ((pos[1][2] - pos[0][2]) / dt[1] - ddz0);
		}
        x2 = mat*rhs;

        // radius
        for (size_t i = 1; i < nVrt-1; ++i)
            rhs[i] = 6.0 / (tSuppPos[i+1] - tSuppPos[i-1]) *
                     ((r[i+1] - r[i]) / dt[i+1]
                     - (r[i] - r[i-1]) / dt[i]);
        if (!isRootNeurite)
			rhs[0] = 0.0;

        xr = mat*rhs;

        // FIXME: find suitable permissible render vector
        vector3 neuriteDir;
        VecSubtract(neuriteDir, pos[nVrt-1], pos[0]);
        VecNormalize(neuriteDir, neuriteDir);
        if (fabs(neuriteDir[0]) < fabs(neuriteDir[1]))
        {
        	if (fabs(neuriteDir[0]) < fabs(neuriteDir[2]))
            	neuriteOut.refDir = vector3(1,0,0);
        	else
            	neuriteOut.refDir = vector3(0,0,1);
        }
        else
        {
			if (fabs(neuriteDir[1]) < fabs(neuriteDir[2]))
				neuriteOut.refDir = vector3(0,1,0);
			else
				neuriteOut.refDir = vector3(0,0,1);
		}

        //neuriteOut.refDir = vector3(0,1/sqrt(2),1/sqrt(2));
		neuriteOut.vSec.reserve(nVrt-1);

        // this will be 0 for root branches and 1 otherwise
        size_t brInd = neuriteOut.vBR.size();
        std::vector<std::pair<size_t, std::vector<size_t> > >::const_iterator brIt;
        std::vector<std::pair<size_t, std::vector<size_t> > >::const_iterator brIt_end;
        if (bpInfo)
        {
            // vBR may already contain an initial BR; now resize to accommodate all others
            neuriteOut.vBR.resize(brInd + bpInfo->size());
            brIt = bpInfo->begin();
            brIt_end = bpInfo->end();
        }

        for (size_t i = 0; i < nVrt-1; ++i)
        {
            NeuriteProjector::Section sec(tSuppPos[i+1]);
            number* param = &sec.splineParamsX[0];
            param[0] = (x0[i]-x0[i+1]) / (6.0 * dt[i+1]);
            param[1] = 0.5 * x0[i+1];
            param[2] = -(dt[i+1]/6.0 * (x0[i] + 2.0*x0[i+1]) + (pos[i+1][0] - pos[i][0]) / dt[i+1]);
            param[3] = pos[i+1][0];
            param = &sec.splineParamsY[0];
            param[0] = (x1[i]-x1[i+1]) / (6.0 * dt[i+1]);
            param[1] = 0.5 * x1[i+1];
            param[2] = -(dt[i+1]/6.0 * (x1[i] + 2.0*x1[i+1]) + (pos[i+1][1] - pos[i][1]) / dt[i+1]);
            param[3] = pos[i+1][1];
            param = &sec.splineParamsZ[0];
            param[0] = (x2[i]-x2[i+1]) / (6.0 * dt[i+1]);
            param[1] = 0.5 * x2[i+1];
            param[2] = -(dt[i+1]/6.0 * (x2[i] + 2.0*x2[i+1]) + (pos[i+1][2] - pos[i][2]) / dt[i+1]);
            param[3] = pos[i+1][2];
            param = &sec.splineParamsR[0];
            param[0] = (xr[i]-xr[i+1]) / (6.0 * dt[i+1]);
            param[1] = 0.5 * xr[i+1];
            param[2] = -(dt[i+1]/6.0 * (xr[i] + 2.0*xr[i+1]) + (r[i+1] - r[i]) / dt[i+1]);
            param[3] = r[i+1];

            // branching points?
            if (bpInfo && brIt != brIt_end && brIt->first == i+1)
            {
                NeuriteProjector::BranchingRegion& br = neuriteOut.vBR[brInd];
                std::vector<size_t>::const_iterator itBranch = brIt->second.begin();
                std::vector<size_t>::const_iterator itBranch_end = brIt->second.end();

                // create BP for parent neurite's BR
                br.bp = make_sp(new NeuriteProjector::BranchingPoint());

                // register parent BR at BP
                br.bp->vNid.push_back(n);
                br.bp->vRegions.push_back(&br);

                br.t = tSuppPos[i+1];

                // loop all child neurites starting at this BP
                for (; itBranch != itBranch_end; ++itBranch)
                {
                    size_t childID = *itBranch;

                    // save parent direction for child neurite
                    parentDirections[childID][0] = -sec.splineParamsX[2];
                    parentDirections[childID][1] = -sec.splineParamsY[2];
                    parentDirections[childID][2] = -sec.splineParamsZ[2];
                	VecNormalize(parentDirections[childID], parentDirections[childID]);

                    // save pointer to BP at child neurite's BR
                    NeuriteProjector::BranchingRegion newChildBR;
                    newChildBR.bp = br.bp;
                    vNeuritesOut[childID].vBR.push_back(newChildBR);
                    NeuriteProjector::BranchingRegion& childBR = vNeuritesOut[childID].vBR[0];

                    // register child BR at BP
                    br.bp->vNid.push_back(childID);
                    br.bp->vRegions.push_back(&childBR);

                    childBR.t = 0;
                }
                ++brInd;
                ++brIt;
            }

            neuriteOut.vSec.push_back(sec);
        }
    }

/*
    // debug
    // print out branching region data
    for (size_t n = 0; n < vNeuritesOut.size(); ++n)
    {
        UG_LOGN("Neurite " << n);
        const NeuriteProjector::Neurite& neurite = vNeuritesOut[n];
        const std::vector<NeuriteProjector::BranchingRegion> vBR = neurite.vBR;

        for (size_t b = 0; b < vBR.size(); ++b)
        {
            UG_LOGN("  BR " << b << ": " << vBR[b].tstart << ".." << vBR[b].tend);
            SmartPtr<NeuriteProjector::BranchingPoint> bp = vBR[b].bp;

            UG_LOGN("    associated BP data:")
            size_t bpSz = bp->vNid.size();
            if (bpSz != bp->vRegions.size())
            {
            	UG_LOGN(      "Size mismatch: vNid " << bpSz << ", vRegions " << bp->vRegions.size());
            }
            else
            {
				for (size_t i = 0; i < bpSz; ++i)
				{
					UG_LOGN("      " << bp->vNid[i] << " (" << bp->vRegions[i]->tstart
						<< ", " << bp->vRegions[i]->tend << ")");
				}
            }
        }
    }
*/
}


number calculate_length_over_radius
(
	number t_start,
	number t_end,
	const NeuriteProjector::Neurite& neurite,
	size_t startSec
)
{
	GaussLegendre gl(5);
	size_t nPts = gl.size();

	std::vector<NeuriteProjector::Section>::const_iterator sec_it = neurite.vSec.begin() + startSec;
	std::vector<NeuriteProjector::Section>::const_iterator sec_end = neurite.vSec.end();

	while (sec_it->endParam < t_start && sec_it->endParam < 1.0)
		++sec_it;

	// check that startSec was correct
	number sec_tstart = startSec > 0 ? (sec_it - 1)->endParam : 0.0;
	number sec_tend = sec_it->endParam;

	UG_COND_THROW(sec_tend < t_start || sec_tstart > t_start,
		"Wrong section iterator given to calc_length_over_radius().\n"
		"Section goes from " << (sec_it-1)->endParam << " to " << sec_it->endParam
		<< ", but t_start is " << t_start << ".");

	number integral = 0.0;
	while (sec_it != sec_end)
	{
		// integrate from t_start to min{t_end, sec_tend}
		const NeuriteProjector::Section& sec = *sec_it;
		sec_tstart = std::max(t_start, sec_it != neurite.vSec.begin() ? (sec_it - 1)->endParam : 0.0);
		sec_tend = std::min(t_end, sec.endParam);
		number dt = sec_tend - sec_tstart;
		number sec_integral = 0.0;
		for (size_t i = 0; i < nPts; ++i)
		{
			number t = sec.endParam - (sec_tstart + dt*gl.point(i)[0]);

			vector3 vel;
			const number* s = &sec.splineParamsX[0];
			number& v0 = vel[0];
			v0 = -3.0*s[0]*t - 2.0*s[1];
			v0 = v0*t - s[2];

			s = &sec.splineParamsY[0];
			number& v1 = vel[1];
			v1 = -3.0*s[0]*t - 2.0*s[1];
			v1 = v1*t - s[2];

			s = &sec.splineParamsZ[0];
			number& v2 = vel[2];
			v2 = -3.0*s[0]*t - 2.0*s[1];
			v2 = v2*t - s[2];

			s = &sec.splineParamsR[0];
			number r = s[0]*t + s[1];
			r = r*t + s[2];
			r = r*t + s[3];

			UG_COND_THROW(r*r <= VecNormSquared(vel)*1e-12, "r = " << r << " at t = " << t << "!");

			sec_integral += gl.weight(i) * sqrt(VecNormSquared(vel)) / r;
		}

		integral += dt * sec_integral;


		// update lower bound and iterator
		t_start = sec_tend;
		if (t_start >= t_end) break;
		++sec_it;
	}

	return integral;
}


void calculate_segment_axial_positions
(
	std::vector<number>& segAxPosOut,
	number t_start,
	number t_end,
	const NeuriteProjector::Neurite& neurite,
	size_t startSec,
	number segLength
)
{
	const size_t nSeg = segAxPosOut.size();

	GaussLegendre gl(5);
	size_t nPts = gl.size();

	std::vector<NeuriteProjector::Section>::const_iterator sec_it = neurite.vSec.begin() + startSec;
	std::vector<NeuriteProjector::Section>::const_iterator sec_end = neurite.vSec.end();

	while (sec_it->endParam < t_start && sec_it->endParam < 1.0)
		++sec_it;

	// check that startSec was correct
	number sec_tstart = startSec > 0 ? (sec_it - 1)->endParam : 0.0;
	number sec_tend = sec_it->endParam;
	UG_COND_THROW(sec_tend < t_start || sec_tstart > t_start,
		"Wrong section iterator given to calculate_segment_axial_positions().\n"
		"Section goes from " << (sec_it-1)->endParam << " to " << sec_it->endParam
		<< ", but t_start is " << t_start << ".");

	number integral = 0.0;
	size_t seg = 0;
	while (sec_it != sec_end)
	{
		// integrate from t_start to min{t_end, sec_tend}
		const NeuriteProjector::Section& sec = *sec_it;
		sec_tstart = std::max(t_start, sec_it != neurite.vSec.begin() ? (sec_it - 1)->endParam : 0.0);
		sec_tend = std::min(t_end, sec.endParam);
		number dt = sec_tend - sec_tstart;
		number sec_integral = 0.0;
		for (size_t i = 0; i < nPts; ++i)
		{
			number t = sec.endParam - (sec_tstart + dt*gl.point(i)[0]);

			vector3 vel;
			const number* s = &sec.splineParamsX[0];
			number& v0 = vel[0];
			v0 = -3.0*s[0]*t - 2.0*s[1];
			v0 = v0*t - s[2];

			s = &sec.splineParamsY[0];
			number& v1 = vel[1];
			v1 = -3.0*s[0]*t - 2.0*s[1];
			v1 = v1*t - s[2];

			s = &sec.splineParamsZ[0];
			number& v2 = vel[2];
			v2 = -3.0*s[0]*t - 2.0*s[1];
			v2 = v2*t - s[2];

			s = &sec.splineParamsR[0];
			number r = s[0]*t + s[1];
			r = r*t + s[2];
			r = r*t + s[3];

			UG_COND_THROW(r*r <= VecNormSquared(vel)*1e-12, "r = " << r << " at t = " << t << "!");

			sec_integral += gl.weight(i) * sqrt(VecNormSquared(vel)) / r;
		}
		integral += dt * sec_integral;

		// calculate exact position by linear interpolation, whenever integral has surpassed it
		while (integral >= (seg+1)*segLength)
		{
			number lastIntegral = integral - dt * sec_integral;
			segAxPosOut[seg] = t_start + ((seg+1)*segLength - lastIntegral) / sec_integral;
			++seg;
		}

		// update lower bound and iterator
		t_start = sec_tend;
		if (t_start >= t_end) break;
		++sec_it;
	}

	// rounding errors may make this necessary
	if (seg+1 == nSeg && (nSeg*segLength - integral)/integral < 1e-6)
	{
		segAxPosOut[nSeg-1] = t_end;
		++seg;
	}

	UG_ASSERT(seg == nSeg, "seg = " << seg << " != " << nSeg << " = nSeg");
}



static void create_neurite
(
	const std::vector<NeuriteProjector::Neurite>& vNeurites,
	const std::vector<std::vector<vector3> >& vPos,
	const std::vector<std::vector<number> >& vR,
	size_t nid,
	number anisotropy,
	Grid& g,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
	std::vector<Vertex*>* connectingVrts = NULL,
	std::vector<Edge*>* connectingEdges = NULL,
	std::vector<Face*>* connectingFaces = NULL,
	number initialOffset = 0.0
)
{
	const NeuriteProjector::Neurite& neurite = vNeurites[nid];
	const std::vector<vector3>& pos = vPos[nid];
	const std::vector<number>& r = vR[nid];

	number neurite_length = 0.0;
	for (size_t i = 1; i < pos.size(); ++i)
		neurite_length += VecDistance(pos[i], pos[i-1]);

	size_t nSec = neurite.vSec.size();

	const std::vector<NeuriteProjector::BranchingRegion>& vBR = neurite.vBR;
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit = vBR.begin();
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit_end = vBR.end();

	std::vector<Vertex*> vVrt;
	std::vector<Edge*> vEdge;
	std::vector<Face*> vFace;
	vVrt.resize(4);
	vEdge.resize(4);
	vFace.resize(1);

	vector3 vel;
	const NeuriteProjector::Section& sec = neurite.vSec[0];
	number h = sec.endParam;
	vel[0] = -3.0*sec.splineParamsX[0]*h*h - 2.0*sec.splineParamsX[1]*h - sec.splineParamsX[2];
	vel[1] = -3.0*sec.splineParamsY[0]*h*h - 2.0*sec.splineParamsY[1]*h - sec.splineParamsY[2];
	vel[2] = -3.0*sec.splineParamsZ[0]*h*h - 2.0*sec.splineParamsZ[1]*h - sec.splineParamsZ[2];

	vector3 projRefDir;
	VecNormalize(vel, vel);
	number fac = VecProd(neurite.refDir, vel);
	VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
	VecNormalize(projRefDir, projRefDir);
	vector3 thirdDir;
	VecCross(thirdDir, vel, projRefDir);
	vector3 previousVel = vel;

	number angleOffset = 0.0;

	number t_start = 0.0;
	number t_end = 0.0;

	if (connectingVrts && connectingEdges && connectingFaces)
	{
		vVrt = *connectingVrts;
		vEdge = *connectingEdges;
		vFace = *connectingFaces;

		// calculate start angle offset
		vector3 center(0.0);
		for (size_t i = 0; i < 4; ++i)
			VecAdd(center, center, aaPos[(*connectingVrts)[i]]);
		center /= 4;

		vector3 centerToFirst;
		VecSubtract(centerToFirst, aaPos[(*connectingVrts)[0]], center);

		vector2 relCoord;
		VecScaleAdd(centerToFirst, 1.0, centerToFirst, -VecProd(centerToFirst, vel), vel);
		relCoord[0] = VecProd(centerToFirst, projRefDir);
		VecScaleAdd(centerToFirst, 1.0, centerToFirst, -relCoord[0], projRefDir);
		relCoord[1] = VecProd(centerToFirst, thirdDir);
		VecNormalize(relCoord, relCoord);

		if (fabs(relCoord[0]) < 1e-8)
			angleOffset = relCoord[1] < 0 ? 1.5*PI : 0.5*PI;
		else
			angleOffset = relCoord[0] < 0 ? PI - atan(-relCoord[1]/relCoord[0]) : atan(relCoord[1]/relCoord[0]);
		if (angleOffset < 0) angleOffset += 2.0*PI;

		// ignore first branching region (the connecting region)
		++brit;

		// apply initial offset (to prevent first segment being shorter than the others)
		t_end = initialOffset / neurite_length;
	}
	else
	{
		// create first layer of vertices/edges //
		// surface vertices first
		for (size_t i = 0; i < 4; ++i)
		{
			Vertex* v = *g.create<RegularVertex>();
			vVrt[i] = v;
			number angle = 0.5*PI*i;
			VecScaleAdd(aaPos[v], 1.0, pos[0], r[0]*cos(angle), projRefDir, r[0]*sin(angle), thirdDir);

			aaSurfParams[v].neuriteID = nid;
			aaSurfParams[v].axial = 0.0;
			aaSurfParams[v].angular = angle;
			aaSurfParams[v].radial = 1.0;
		}

		// edges
		for (size_t i = 0; i < 4; ++i)
			vEdge[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[(i+1)%4]));

		// faces
		vFace[0] = *g.create<Quadrilateral>(QuadrilateralDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3]));
	}

	// Now create dendrite to the next branching point and iterate this process.
	// We want to create each of the segments with approx. the same aspect ratio.
	// To that end, we first calculate the length of the section to be created (in units of radius)
	// and then divide this number by 2^n where n is the number of anisotropic refinements to
	// be performed to make all segments (more or less) isotropic. The result is the number
	// of segments to be used for the section.
	vector3 lastPos = pos[0];
	size_t curSec = 0;

	while (true)
	{
		t_start = t_end;

		// if there is another BP, store its extensions here
		number bp_start = 1.0;
		number bp_end = 0.0;

		// initial branch offsets (to prevent the first segment being shorter than the following)
		std::vector<number> branchOffset;
		number surfBPoffset;

		// last section: create until tip
		if (brit == brit_end)
			t_end = 1.0;

		// otherwise: section goes to next branching point
		else
		{
			// calculate the exact position of the branching point,
			// i.e., the axial position of the intersection of the branching neurite's
			// spline with the surface of the current neurite
			// this is necessary esp. when the branching angle is very small
			const std::vector<uint32_t>& vBranchInd = brit->bp->vNid;
			size_t nBranches = vBranchInd.size();

			branchOffset.resize(brit->bp->vNid.size(), 0.0);
			for (size_t br = 1; br < nBranches; ++br)
			{
				uint32_t brInd = vBranchInd[br];

				// get position and radius of first point of branch
				const number brRadSeg1 = vR[brInd][0];

				// get position and radius of branching point
				const number bpTPos = brit->t;
				size_t brSec = curSec;
				for (; brSec < nSec; ++brSec)
				{
					const NeuriteProjector::Section& sec = neurite.vSec[brSec];
					if (bpTPos - sec.endParam < 1e-6*bpTPos)
						break;
				}
				UG_COND_THROW(brSec == nSec, "Could not find section containing branching point "
					"at t = " << bpTPos << ".");
				const number bpRad = vR[nid][brSec+1];

				// calculate branch and neurite directions
				vector3 branchDir;
				const NeuriteProjector::Section& childSec = vNeurites[brInd].vSec[0];
				number te = childSec.endParam;

				const number* s = &childSec.splineParamsX[0];
				number& v0 = branchDir[0];
				v0 = -3.0*s[0]*te - 2.0*s[1];
				v0 = v0*te - s[2];

				s = &childSec.splineParamsY[0];
				number& v1 = branchDir[1];
				v1 = -3.0*s[0]*te - 2.0*s[1];
				v1 = v1*te - s[2];

				s = &childSec.splineParamsZ[0];
				number& v2 = branchDir[2];
				v2 = -3.0*s[0]*te - 2.0*s[1];
				v2 = v2*te - s[2];

				VecNormalize(branchDir, branchDir);

				vector3 neuriteDir;
				const NeuriteProjector::Section& sec = neurite.vSec.at(brSec);
				vel[0] = -sec.splineParamsX[2];
				vel[1] = -sec.splineParamsY[2];
				vel[2] = -sec.splineParamsZ[2];
				number velNorm = sqrt(VecNormSquared(vel));
				VecScale(neuriteDir, vel, 1.0/velNorm);

				// calculate offset of true branch position, which is r1/sqrt(2)*cot(alpha)
				const number brScProd = VecProd(neuriteDir, branchDir);
				const number sinAlphaInv = 1.0 / sqrt(1.0 - brScProd*brScProd);
				surfBPoffset = 0.5*sqrt(2.0) * bpRad * brScProd * sinAlphaInv;

				// calculate offset of new branch, which is r1/sqrt(2)/sin(alpha)
				branchOffset[br] = 0.5*sqrt(2.0) * bpRad * sinAlphaInv;

				// calculate true half-length of BP, which is r2/sin(alpha)
				number surfBPhalfLength = brRadSeg1 * sinAlphaInv;

				// finally set bp start and end
				bp_start = std::min(bp_start, bpTPos + (surfBPoffset - surfBPhalfLength) / neurite_length);
				bp_end = std::max(bp_end, bpTPos + (surfBPoffset + surfBPhalfLength) / neurite_length);
			}

			t_end = bp_start;
		}

		// calculate total length in units of radius
		// = integral from t_start to t_end over: ||v(t)|| / r(t) dt
		number lengthOverRadius = calculate_length_over_radius(t_start, t_end, neurite, curSec);

		// to reach the desired anisotropy on the surface in the refinement limit,
		// it has to be multiplied by pi/2 h
		size_t nSeg = (size_t) round(lengthOverRadius / (anisotropy*0.5*PI));
		if (!nSeg)
			nSeg = 1;
		number segLength = lengthOverRadius / nSeg;	// segments are between 8 and 16 radii long
		std::vector<number> vSegAxPos(nSeg);
		calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neurite, curSec, segLength);

		// add the branching point to segment list (if present)
		if (brit != brit_end)
		{
			vSegAxPos.resize(nSeg+1);
			vSegAxPos[nSeg] = bp_end;
			++nSeg;
		}


		// in case we construct to a BP, find out the branching angle
		// in order to adjust this neurites offset bit by bit
		number addOffset = 0.0;
		size_t child_nid;
		size_t connFaceInd = 0;
		if (brit != brit_end)
		{
			// find branching child neurite
			SmartPtr<NeuriteProjector::BranchingPoint> bp = brit->bp;
			UG_COND_THROW(bp->vNid.size() > 2,
				"This implementation can only handle branching points with one branching child.");

			if (bp->vNid[0] != nid)
				child_nid = bp->vNid[0];
			else
				child_nid = bp->vNid[1];

			// find out branching child neurite initial direction
			vector3 childDir;
			const NeuriteProjector::Section& childSec = vNeurites[child_nid].vSec[0];
			number te = childSec.endParam;

			const number* sp = &childSec.splineParamsX[0];
			number& vc0 = childDir[0];
			vc0 = -3.0*sp[0]*te - 2.0*sp[1];
			vc0 = vc0*te - sp[2];

			sp = &childSec.splineParamsY[0];
			number& vc1 = childDir[1];
			vc1 = -3.0*sp[0]*te - 2.0*sp[1];
			vc1 = vc1*te - sp[2];

			sp = &childSec.splineParamsZ[0];
			number& vc2 = childDir[2];
			vc2 = -3.0*sp[0]*te - 2.0*sp[1];
			vc2 = vc2*te - sp[2];


			// find out neurite direction in next BP
			number bpAxPos = vSegAxPos[nSeg-1];
			size_t tmpSec = curSec;
			for (; tmpSec < nSec; ++tmpSec)
			{
				const NeuriteProjector::Section& sec = neurite.vSec[tmpSec];
				if (sec.endParam >= bpAxPos)
					break;
			}

			const NeuriteProjector::Section& sec = neurite.vSec[tmpSec];
			number monom = sec.endParam - bpAxPos;
			sp = &sec.splineParamsX[0];
			number& v0 = vel[0];
			v0 = -3.0*sp[0]*monom -2.0*sp[1];
			v0 = v0*monom - sp[2];

			sp = &sec.splineParamsY[0];
			number& v1 = vel[1];
			v1 = -3.0*sp[0]*monom -2.0*sp[1];
			v1 = v1*monom - sp[2];

			sp = &sec.splineParamsZ[0];
			number& v2 = vel[2];
			v2 = -3.0*sp[0]*monom -2.0*sp[1];
			v2 = v2*monom - sp[2];

			VecNormalize(vel, vel);

			// calculate offset
			number fac = VecProd(neurite.refDir, vel);
			VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
			VecNormalize(projRefDir, projRefDir);
			VecCross(thirdDir, vel, projRefDir);

			/*
			vector2 relCoord;
			VecScaleAppend(childDir, -VecProd(childDir, vel), vel);
			relCoord[0] = VecProd(childDir, projRefDir);
			VecScaleAppend(childDir, -relCoord[0], projRefDir);
			relCoord[1] = VecProd(childDir, thirdDir);
			VecNormalize(relCoord, relCoord);
			*/
			vector2 relCoord(VecProd(childDir, projRefDir), VecProd(childDir, thirdDir));
			VecNormalize(relCoord, relCoord);

			number branchOffset = 0.0;
			if (fabs(relCoord[0]) < 1e-8)
				branchOffset = relCoord[1] < 0 ? 1.5*PI : 0.5*PI;
			else
				branchOffset = relCoord[0] < 0 ? PI - atan(-relCoord[1]/relCoord[0]) : atan(relCoord[1]/relCoord[0]);

			addOffset = branchOffset - angleOffset;
			connFaceInd = floor(std::fmod(addOffset+4*PI, 2*PI) / (PI/2));
			addOffset = std::fmod(addOffset - (connFaceInd*PI/2 + PI/4) + 4*PI, 2*PI);
			if (addOffset > PI)
				addOffset -= 2*PI;
			addOffset /= nSeg - 1;
		}

		// create mesh for segments
		Selector sel(g);
		for (size_t s = 0; s < nSeg; ++s)
		{
			// get exact position, velocity and radius of segment end
			number segAxPos = vSegAxPos[s];
			for (; curSec < nSec; ++curSec)
			{
				const NeuriteProjector::Section& sec = neurite.vSec[curSec];
				if (sec.endParam >= segAxPos)
					break;
			}

			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			vector3 curPos;
			number monom = sec.endParam - segAxPos;
			const number* sp = &sec.splineParamsX[0];
			number& p0 = curPos[0];
			number& v0 = vel[0];
			p0 = sp[0]*monom + sp[1];
			p0 = p0*monom + sp[2];
			p0 = p0*monom + sp[3];
			v0 = -3.0*sp[0]*monom -2.0*sp[1];
			v0 = v0*monom - sp[2];

			sp = &sec.splineParamsY[0];
			number& p1 = curPos[1];
			number& v1 = vel[1];
			p1 = sp[0]*monom + sp[1];
			p1 = p1*monom + sp[2];
			p1 = p1*monom + sp[3];
			v1 = -3.0*sp[0]*monom -2.0*sp[1];
			v1 = v1*monom - sp[2];

			sp = &sec.splineParamsZ[0];
			number& p2 = curPos[2];
			number& v2 = vel[2];
			p2 = sp[0]*monom + sp[1];
			p2 = p2*monom + sp[2];
			p2 = p2*monom + sp[3];
			v2 = -3.0*sp[0]*monom -2.0*sp[1];
			v2 = v2*monom - sp[2];

			sp = &sec.splineParamsR[0];
			number radius;
			radius = sp[0]*monom + sp[1];
			radius = radius*monom + sp[2];
			radius = radius*monom + sp[3];

			// Orient the cross-section with the same center-to-center direction
			// that is used for the extrusion below. The analytic cubic-spline
			// derivative can become nearly zero and reverse abruptly at a cusp,
			// even while consecutive generated center points remain well behaved.
			vector3 chordVel;
			VecScaleAdd(chordVel, 1.0, curPos, -1.0, lastPos);
			number chordLength = VecLength(chordVel);
			if (chordLength > 1e-12)
				VecScale(vel, chordVel, 1.0 / chordLength);
			else
				VecNormalize(vel, vel);

			// Rotate the frame by the unique minimal rotation that maps the
			// previous tangent to the current tangent (discrete Bishop frame).
			vector3 rotationAxis;
			VecCross(rotationAxis, previousVel, vel);
			number sinAngle = VecLength(rotationAxis);
			number cosAngle = VecProd(previousVel, vel);
			cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
			if (cosAngle < 0.8660254037844386)
			{
				UG_LOGN(
					"WARNING: sharp generated centerline turn."
					<< " neurite=" << nid
					<< ", segment=" << s
					<< ", angleDeg=" << acos(cosAngle) * 180.0 / PI
					<< ", position=("
					<< curPos[0] << ", "
					<< curPos[1] << ", "
					<< curPos[2] << ")"
				);
			}

			vector3 transportedRefDir = projRefDir;
			if (sinAngle > 1e-12)
			{
				VecScale(rotationAxis, rotationAxis, 1.0 / sinAngle);
				vector3 axisCrossRef;
				VecCross(axisCrossRef, rotationAxis, projRefDir);
				number axisDotRef = VecProd(rotationAxis, projRefDir);
				VecScaleAdd(transportedRefDir,
				            cosAngle, projRefDir,
				            sinAngle, axisCrossRef,
				            (1.0-cosAngle)*axisDotRef, rotationAxis);
			}

			// Remove accumulated round-off in the normal direction.
			number normalComponent = VecProd(transportedRefDir, vel);
			VecScaleAdd(projRefDir, 1.0, transportedRefDir,
			            -normalComponent, vel);
			number frameLength = VecLength(projRefDir);
			if (frameLength <= 1e-12)
			{
				// This should only occur at an exact tangent reversal. Rebuild a
				// perpendicular direction from the neurite's fixed reference axis.
				number fallbackFac = VecProd(neurite.refDir, vel);
				VecScaleAdd(projRefDir, 1.0, neurite.refDir,
				            -fallbackFac, vel);
				frameLength = VecLength(projRefDir);
				UG_LOGN(
					"WARNING: rebuilt degenerate neurite frame."
					<< " neurite=" << nid
					<< ", segment=" << s
					<< ", position=("
					<< curPos[0] << ", "
					<< curPos[1] << ", "
					<< curPos[2] << ")"
				);
			}
			UG_COND_THROW(frameLength <= 1e-12,
			              "Could not construct neurite frame at segment " << s);
			VecNormalize(projRefDir, projRefDir);

			VecCross(thirdDir, vel, projRefDir);
			previousVel = vel;

			// usual segment: extrude
			if (s != nSeg - 1 || brit == brit_end)
			{
				// apply additional offset
				angleOffset = std::fmod(angleOffset + addOffset + 2*PI, 2*PI);

				// extrude from last pos to new pos
				vector3 extrudeDir;
				VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, lastPos);
				std::vector<Volume*> vVol;
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos, EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);

				// set new positions and param attachments
				for (size_t j = 0; j < 4; ++j)
				{
					number angle = 0.5*PI*j + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					Vertex* v = vVrt[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// ensure correct volume orientation
				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);
			}

			// BP segment: create BP with tetrahedra/pyramids and create whole branch
			else
			{
				std::vector<Vertex*> vNewVrt(4, NULL);

				// create all needed vertices (except branch mid point)
				for (size_t j = 0; j < 4; ++j)
				{
					Vertex* v = vNewVrt[j] = *g.create<RegularVertex>();

					number angle = 0.5*PI*j + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					vector3 radialVec;
					VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// correct offsets of non-connecting vertices (opposite direction than branching ones)
				VecScaleAppend(aaPos[vVrt[(connFaceInd+2)%4]], -2.0*surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd+3)%4]], -2.0*surfBPoffset, vel);
				VecScaleAppend(aaPos[vNewVrt[(connFaceInd+2)%4]], -2.0*surfBPoffset, vel);
				VecScaleAppend(aaPos[vNewVrt[(connFaceInd+3)%4]], -2.0*surfBPoffset, vel);

				aaSurfParams[vVrt[(connFaceInd+2)%4]].axial -= 2.0*surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd+3)%4]].axial -= 2.0*surfBPoffset / neurite_length;
				aaSurfParams[vNewVrt[(connFaceInd+2)%4]].axial -= 2.0*surfBPoffset / neurite_length;
				aaSurfParams[vNewVrt[(connFaceInd+3)%4]].axial -= 2.0*surfBPoffset / neurite_length;

				// prepare connectingVrts, -Edges, -Faces for branch
				std::vector<Vertex*> vBranchVrts(4);
				vBranchVrts[0] = vVrt[(connFaceInd+1)%4];
				vBranchVrts[1] = vNewVrt[(connFaceInd+1)%4];
				vBranchVrts[2] = vNewVrt[connFaceInd];
				vBranchVrts[3] = vVrt[connFaceInd];

				std::vector<Edge*> vBranchEdges(4);
				std::vector<Face*> vBranchFaces(1);
				for (size_t j = 0; j < 4; ++j)
				{
					if (j != 3)
						vBranchEdges[j] = *g.create<RegularEdge>(EdgeDescriptor(vBranchVrts[j], vBranchVrts[(j+1)%4]));
					else
						vBranchEdges[j] = vEdge[connFaceInd];
				}
				vBranchFaces[0] = *g.create<Quadrilateral>(QuadrilateralDescriptor(vBranchVrts[0], vBranchVrts[1], vBranchVrts[2], vBranchVrts[3]));

				// add branch neurite ID to its initial vertices
				for (size_t j = 0; j < 4; ++j)
				{
					Vertex* brVrt = vBranchVrts[j];
					aaSurfParams[brVrt].neuriteID += (brit - vBR.begin()) << 20;  // add branching region index
					aaSurfParams[brVrt].neuriteID += 1 << 28;  // add child ID (always 0, since there can only be one child here)
				}

				// recursively build branch
				create_neurite(vNeurites, vPos, vR, child_nid, anisotropy,
					g, aaPos, aaSurfParams, &vBranchVrts, &vBranchEdges, &vBranchFaces, branchOffset[1]);


				// prepare connectingVrts, -Edges, -Faces for remainder of own neurite
				for (size_t j = 0; j < 4; ++j)
				{
					if (j != connFaceInd)
						vEdge[j] = *g.create<RegularEdge>(EdgeDescriptor(vNewVrt[j], vNewVrt[(j+1)%4]));
					else
						vEdge[j] = vBranchEdges[1];
				}

				vFace[0] = *g.create<Quadrilateral>(QuadrilateralDescriptor(vNewVrt[0], vNewVrt[1], vNewVrt[2], vNewVrt[3]));


				// create all inner BP elements
				g.create<Hexahedron>(HexahedronDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3],
					vNewVrt[0], vNewVrt[1], vNewVrt[2], vNewVrt[3]));

				/*
				Vertex* innerVrt = *g.create<RegularVertex>();

				VecScaleAdd(aaPos[innerVrt], 0.5, curPos, 0.5, lastPos);
				aaSurfParams[innerVrt].neuriteID = nid;
				aaSurfParams[innerVrt].neuriteID += (brit - vBR.begin()) << 20;  // add branching region index
				aaSurfParams[innerVrt].neuriteID += 1 << 28;  // add child ID (always 0, since there can only be one child here)
				aaSurfParams[innerVrt].axial = 0.5*vSegAxPos[s] + 0.5*vSegAxPos[s-1];
				aaSurfParams[innerVrt].angular = 0.0;
				aaSurfParams[innerVrt].radial = 0.0;

				g.create<Pyramid>(PyramidDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3], innerVrt));
				g.create<Pyramid>(PyramidDescriptor(vNewVrt[3], vNewVrt[2], vNewVrt[1], vNewVrt[0], innerVrt));
				for (size_t j = 0; j < 4; ++j)
					g.create<Pyramid>(PyramidDescriptor(vVrt[j], vNewVrt[j], vNewVrt[(j+1)%4], vVrt[(j+1)%4], innerVrt));
				*/

				vVrt.swap(vNewVrt);
			}

			lastPos = curPos;
		}

		// update t_end and curSec
		if (brit != brit_end)
			t_end = bp_end;

		for (; curSec < nSec; ++curSec)
		{
			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			if (sec.endParam >= t_end)
				break;
		}

		// check whether tip has been reached
		if (brit == brit_end)
			break;
		else
			++brit;
	}

	/*
	// close the tip of the neurite
	const NeuriteProjector::Section& lastSec = neurite.vSec[nSec-1];
	vel = vector3(-lastSec.splineParamsX[2], -lastSec.splineParamsY[2], -lastSec.splineParamsZ[2]);
	number radius = lastSec.splineParamsR[3];
	VecScale(vel, vel, radius/sqrt(VecProd(vel, vel)));
	Vertex* tip = *g.create<RegularVertex>();
	for (size_t i = 0; i < 4; ++i)
		VecScaleAppend(aaPos[tip], 0.25, aaPos[vVrt[i]]);
	VecAdd(aaPos[tip], aaPos[tip], vel);
	g.create<Pyramid>(PyramidDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3], tip));

	aaSurfParams[tip].neuriteID = nid;
	aaSurfParams[tip].axial = 2.0;
	aaSurfParams[tip].angular = 0.0;
	aaSurfParams[tip].radial = 1.0;
	*/
}


static void create_neurite_with_er
(
	const std::vector<NeuriteProjector::Neurite>& vNeurites,
	const std::vector<std::vector<vector3> >& vPos,
	const std::vector<std::vector<number> >& vR,
	size_t nid,
	number erScaleFactor,
	number anisotropy,
	size_t nephronOGridVertices,
	bool coarseLumenCenterOnly,
	Grid& g,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
	SubsetHandler& sh,
	std::vector<Vertex*>* connectingVrts = NULL,
	std::vector<Edge*>* connectingEdges = NULL,
	std::vector<Face*>* connectingFaces = NULL,
	number initialOffset = 0.0
)
{
	const NeuriteProjector::Neurite& neurite = vNeurites[nid];
	const std::vector<vector3>& pos = vPos[nid];
	const std::vector<number>& r = vR[nid];

	number neurite_length = 0.0;
	for (size_t i = 1; i < pos.size(); ++i)
		neurite_length += VecDistance(pos[i], pos[i-1]);

	size_t nSec = neurite.vSec.size();

	const std::vector<NeuriteProjector::BranchingRegion>& vBR = neurite.vBR;
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit = vBR.begin();
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit_end = vBR.end();

	// A non-branching nephron uses matching lumen-midpoint, Apical,
	// membrane-midpoint, and Basolateral O-grid rings. The two midpoint rings
	// give both the coarse Lumen and Membrane an internal vertex layer.
	// Circumferential resolution remains user-controlled. The legacy neuronal
	// branch construction is untouched.
	const bool useMatchingNephronOGrid =
		vBR.empty() && !connectingVrts && !connectingEdges && !connectingFaces;
	const size_t nOGrid = useMatchingNephronOGrid ? nephronOGridVertices : 12;
	UG_COND_THROW(useMatchingNephronOGrid && nOGrid < 4,
		"The nephron O-grid requires at least four circumferential vertices.");
	const number lumenRadial = 0.5 * erScaleFactor;
	const number middleRadial = 0.5 * (erScaleFactor + 1.0);
	const size_t numRadialRings = coarseLumenCenterOnly ? 3 : 4;
	std::vector<number> ringRadial;
	std::vector<int> ringSubset;
	if (coarseLumenCenterOnly)
	{
		ringRadial.push_back(erScaleFactor);
		ringRadial.push_back(middleRadial);
		ringRadial.push_back(1.0);
		ringSubset.push_back(3); // Apical
		ringSubset.push_back(0); // Membrane midpoint
		ringSubset.push_back(2); // Basolateral
	}
	else
	{
		ringRadial.push_back(lumenRadial);
		ringRadial.push_back(erScaleFactor);
		ringRadial.push_back(middleRadial);
		ringRadial.push_back(1.0);
		ringSubset.push_back(1); // Lumen midpoint
		ringSubset.push_back(3); // Apical
		ringSubset.push_back(0); // Membrane midpoint
		ringSubset.push_back(2); // Basolateral
	}
	std::vector<size_t> ringBegin(numRadialRings);
	for (size_t layer = 0; layer < numRadialRings; ++layer)
		ringBegin[layer] = 1 + layer*nOGrid;

	std::vector<Vertex*> vVrt;
	std::vector<Edge*> vEdge;
	std::vector<Face*> vFace;
	vVrt.resize(useMatchingNephronOGrid ? 1 + numRadialRings*nOGrid : 16);
	vEdge.resize(useMatchingNephronOGrid ? 2*numRadialRings*nOGrid : 24);
	vFace.resize(useMatchingNephronOGrid ? numRadialRings*nOGrid : 9);

	vector3 vel;
	const NeuriteProjector::Section& sec = neurite.vSec[0];
	number h = sec.endParam;
	vel[0] = -3.0*sec.splineParamsX[0]*h*h - 2.0*sec.splineParamsX[1]*h - sec.splineParamsX[2];
	vel[1] = -3.0*sec.splineParamsY[0]*h*h - 2.0*sec.splineParamsY[1]*h - sec.splineParamsY[2];
	vel[2] = -3.0*sec.splineParamsZ[0]*h*h - 2.0*sec.splineParamsZ[1]*h - sec.splineParamsZ[2];

	vector3 projRefDir;
	vector3 thirdDir;
	if (useMatchingNephronOGrid)
		compute_nephron_parallel_transport_ONB(
			vel, projRefDir, thirdDir, neurite, 0.0);
	else
	{
		VecNormalize(vel, vel);
		number fac = VecProd(neurite.refDir, vel);
		VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
		VecNormalize(projRefDir, projRefDir);
		VecCross(thirdDir, vel, projRefDir);
	}

	number angleOffset = 0.0;

	number t_start = 0.0;
	number t_end = 0.0;

	if (connectingVrts && connectingEdges && connectingFaces)
	{
		vVrt = *connectingVrts;
		vEdge = *connectingEdges;
		vFace = *connectingFaces;

		// calculate start angle offset
		vector3 center(0.0);
		for (size_t i = 0; i < 4; ++i)
			VecAdd(center, center, aaPos[(*connectingVrts)[i]]);
		center /= 4;

		vector3 centerToFirst;
		VecSubtract(centerToFirst, aaPos[(*connectingVrts)[0]], center);

		vector2 relCoord;
		relCoord[0] = VecProd(centerToFirst, projRefDir);
		relCoord[1] = VecProd(centerToFirst, thirdDir);
		VecNormalize(relCoord, relCoord);

		if (fabs(relCoord[0]) < 1e-8)
			angleOffset = relCoord[1] < 0 ? 1.5*PI : 0.5*PI;
		else
			angleOffset = relCoord[0] < 0 ? PI - atan(-relCoord[1]/relCoord[0]) : atan(relCoord[1]/relCoord[0]);
		if (angleOffset < 0) angleOffset += 2.0*PI;

		// ignore first branching region (the connecting region)
		++brit;

		// apply initial offset (to prevent first segment being shorter than the others)
		t_end = initialOffset / neurite_length;
		for (size_t i = 0; i < 4; ++i)
		{
			aaSurfParams[(*connectingVrts)[i]].axial = t_end;
			aaSurfParams[(*connectingVrts)[i]].angular = 0.5*PI*i + angleOffset < 2*PI ? 0.5*PI*i + angleOffset : 0.5*PI*i + angleOffset - 2*PI;
			aaSurfParams[(*connectingVrts)[i]].radial = erScaleFactor;
		}
	}
	else
	{
		// create first layer of vertices/edges //
		if (useMatchingNephronOGrid)
		{
			// Vertex 0 is the lumen center. It is followed by matching rings
			// inside Lumen, on Apical, inside Membrane, and on Basolateral.
			Vertex* center = *g.create<RegularVertex>();
			vVrt[0] = center;
			aaPos[center] = pos[0];
			aaSurfParams[center].neuriteID = nid;
			aaSurfParams[center].axial = 0.0;
			aaSurfParams[center].angular = 0.0;
			aaSurfParams[center].radial = 0.0;
			sh.assign_subset(center, 1);

			for (size_t i = 0; i < nOGrid; ++i)
			{
				const number angle = 2.0*PI*((number)i/(number)nOGrid);
				for (size_t layer = 0; layer < numRadialRings; ++layer)
				{
					Vertex* v = *g.create<RegularVertex>();
					vVrt[ringBegin[layer]+i] = v;
					VecScaleAdd(aaPos[v], 1.0, pos[0],
						ringRadial[layer]*r[0]*cos(angle), projRefDir,
						ringRadial[layer]*r[0]*sin(angle), thirdDir);
					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = 0.0;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = ringRadial[layer];
					sh.assign_subset(v, ringSubset[layer]);
				}
			}

			for (size_t i = 0; i < nOGrid; ++i)
			{
				const size_t next = (i+1)%nOGrid;
				for (size_t layer = 0; layer < numRadialRings; ++layer)
				{
					Vertex* previous = layer == 0 ? center
						: vVrt[ringBegin[layer-1]+i];
					Vertex* current = vVrt[ringBegin[layer]+i];
					vEdge[(2*layer)*nOGrid+i] = *g.create<RegularEdge>(
						EdgeDescriptor(previous,current));
					vEdge[(2*layer+1)*nOGrid+i] = *g.create<RegularEdge>(
						EdgeDescriptor(current,vVrt[ringBegin[layer]+next]));
					sh.assign_subset(vEdge[(2*layer)*nOGrid+i],
					                 layer == 0 ? 1 : (layer == 1 && !coarseLumenCenterOnly ? 1 : 0));
					sh.assign_subset(vEdge[(2*layer+1)*nOGrid+i],ringSubset[layer]);
					if (layer == 0)
						vFace[i] = *g.create<Triangle>(TriangleDescriptor(
							center,current,vVrt[ringBegin[layer]+next]));
					else
						vFace[layer*nOGrid+i] = *g.create<Quadrilateral>(QuadrilateralDescriptor(
							vVrt[ringBegin[layer-1]+i],current,
							vVrt[ringBegin[layer]+next],vVrt[ringBegin[layer-1]+next]));
					sh.assign_subset(vFace[layer*nOGrid+i],
					                 layer <= (coarseLumenCenterOnly ? 0 : 1) ? 1 : 0);
				}
			}
		}
		else
		{

		// ER vertices
		for (size_t i = 0; i < 4; ++i)
		{
			Vertex* v = *g.create<RegularVertex>();
			vVrt[i] = v;
			number angle = 0.5*PI*i;
			VecScaleAdd(aaPos[v], 1.0, pos[0], erScaleFactor*r[0]*cos(angle), projRefDir, erScaleFactor*r[0]*sin(angle), thirdDir);

			aaSurfParams[v].neuriteID = nid;
			aaSurfParams[v].axial = 0.0;
			aaSurfParams[v].angular = angle;
			aaSurfParams[v].radial = erScaleFactor;
			sh.assign_subset(v, 3);
		}

		for (size_t i = 0; i < 12; ++i)
		{
			Vertex* v = *g.create<RegularVertex>();
			vVrt[i+4] = v;
			number angle = PI*((number)i/6);
			VecScaleAdd(aaPos[v], 1.0, pos[0], r[0]*cos(angle), projRefDir, r[0]*sin(angle), thirdDir);

			aaSurfParams[v].neuriteID = nid;
			aaSurfParams[v].axial = 0.0;
			aaSurfParams[v].angular = angle;
			aaSurfParams[v].radial = 1.0;
			sh.assign_subset(v, 2);
		}

		// edges
		for (size_t i = 0; i < 4; ++i)
		{
			vEdge[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[(i+1)%4]));
			vEdge[i+4] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[5+3*i]));
			vEdge[i+8] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[(i+1)%4], vVrt[6+3*i]));

			sh.assign_subset(vEdge[i], 3);
			sh.assign_subset(vEdge[i+4], 0);
		}
		for (size_t i = 0; i < 12; ++i)
		{
			vEdge[i+12] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i+4], vVrt[(i+1)%12+4]));
			sh.assign_subset(vEdge[i+12], 2);
		}

		// faces
		vFace[0] = *g.create<Quadrilateral>(QuadrilateralDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3]));
		sh.assign_subset(vFace[0], 1);
		for (size_t i = 0; i < 4; ++i)
		{
			vFace[i+1] = *g.create<Quadrilateral>(QuadrilateralDescriptor(vVrt[i], vVrt[(3*i+11)%12+4], vVrt[3*i+4], vVrt[3*i+5]));
			vFace[i+5] = *g.create<Quadrilateral>(QuadrilateralDescriptor(vVrt[i], vVrt[3*i+5], vVrt[3*i+6], vVrt[(i+1)%4]));
			sh.assign_subset(vFace[i+1], 0);
			sh.assign_subset(vFace[i+5], 0);
		}
		}
	}

	// Now create dendrite to the next branching point and iterate this process.
	// We want to create each of the segments with approx. the same aspect ratio.
	// To that end, we first calculate the length of the section to be created (in units of radius)
	// and then divide this number by 2^n where n is the number of anisotropic refinements to
	// be performed to make all segments (more or less) isotropic. The result is the number
	// of segments to be used for the section.
	vector3 prevPos = pos[0];
	size_t curSec = 0;

	while (true)
	{
		t_start = t_end;

		// if there is another BP, store its extensions here
		number bp_start = 1.0;
		number bp_end = 0.0;

		// initial branch offsets (to prevent the first segment being shorter than the following)
		std::vector<number> branchOffset;
		number surfBPoffset;

		// last section: create until tip
		if (brit == brit_end)
			t_end = 1.0;

		// otherwise: section goes to next branching point
		else
		{
			// calculate the exact position of the branching point,
			// i.e., the axial position of the intersection of the branching neurite's
			// spline with the surface of the current neurite
			// this is necessary esp. when the branching angle is very small
			const std::vector<uint32_t>& vBranchInd = brit->bp->vNid;
			size_t nBranches = vBranchInd.size();

			branchOffset.resize(brit->bp->vNid.size(), 0.0);
			for (size_t br = 1; br < nBranches; ++br)
			{
				uint32_t brInd = vBranchInd[br];

				// get position and radius of first point of branch
				const number brRadSeg1 = vR[brInd][0];

				// get position and radius of branching point
				const number bpTPos = brit->t;
				size_t brSec = curSec;
				for (; brSec < nSec; ++brSec)
				{
					const NeuriteProjector::Section& sec = neurite.vSec[brSec];
					if (bpTPos - sec.endParam < 1e-6*bpTPos)
						break;
				}
				UG_COND_THROW(brSec == nSec, "Could not find section containing branching point "
					"at t = " << bpTPos << ".");
				const number bpRad = vR[nid][brSec+1];

				// calculate branch and neurite directions
				vector3 branchDir;
				const NeuriteProjector::Section& childSec = vNeurites[brInd].vSec[0];
				number te = childSec.endParam;

				const number* s = &childSec.splineParamsX[0];
				number& v0 = branchDir[0];
				v0 = -3.0*s[0]*te - 2.0*s[1];
				v0 = v0*te - s[2];

				s = &childSec.splineParamsY[0];
				number& v1 = branchDir[1];
				v1 = -3.0*s[0]*te - 2.0*s[1];
				v1 = v1*te - s[2];

				s = &childSec.splineParamsZ[0];
				number& v2 = branchDir[2];
				v2 = -3.0*s[0]*te - 2.0*s[1];
				v2 = v2*te - s[2];

				VecNormalize(branchDir, branchDir);

				vector3 neuriteDir;
				const NeuriteProjector::Section& sec = neurite.vSec.at(brSec);
				vel[0] = -sec.splineParamsX[2];
				vel[1] = -sec.splineParamsY[2];
				vel[2] = -sec.splineParamsZ[2];
				number velNorm = sqrt(VecNormSquared(vel));
				VecScale(neuriteDir, vel, 1.0/velNorm);

				// calculate offset of true branch position, which is r1/sqrt(2)*cot(alpha)
				const number brScProd = VecProd(neuriteDir, branchDir);
				const number sinAlphaInv = 1.0 / sqrt(1.0 - brScProd*brScProd);
				surfBPoffset = 0.5*sqrt(2.0) * bpRad * brScProd * sinAlphaInv;

				// calculate offset of new branch, which is r1/sqrt(2)/sin(alpha)
				branchOffset[br] = 0.5*sqrt(2.0) * bpRad * sinAlphaInv;

				// calculate true half-length of BP, which is r2/sqrt(2)/sin(alpha)
				number surfBPhalfLength = 0.5*sqrt(2.0) * brRadSeg1 * sinAlphaInv;

				// finally set bp start and end
				bp_start = std::min(bp_start, bpTPos - surfBPhalfLength / neurite_length);
				bp_end = std::max(bp_end, bpTPos + surfBPhalfLength / neurite_length);
			}

			t_end = bp_start;
		}

		// calculate total length in units of radius
		// = integral from t_start to t_end over: ||v(t)|| / r(t) dt
		number lengthOverRadius = calculate_length_over_radius(t_start, t_end, neurite, curSec);

		// to reach the desired anisotropy on the surface in the refinement limit,
		// it has to be multiplied by pi/6 h
		size_t nSeg = (size_t) round(lengthOverRadius / (anisotropy*0.16666666*PI));
		if (!nSeg)
			nSeg = 1;
		number segLength = lengthOverRadius / nSeg;	// segments are between 8 and 16 radii long
		std::vector<number> vSegAxPos(nSeg);
		calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neurite, curSec, segLength);

		// Keep spline/chord midpoint error below 20% of local membrane thickness.
		// A denser 5% threshold produced nearly coplanar rings that destabilized
		// TetGen facet recovery; 20% remains small relative to the shell width.
		const bool enableCurvatureSubdivision =
			g_adaptiveNephronMeshing.curvatureLimitedAnisotropy
			&& anisotropy > g_adaptiveNephronMeshing.minimumBendAnisotropy;
		if (useMatchingNephronOGrid && enableCurvatureSubdivision)
		{
			auto centerAndRadiusAt = [&] (number axial, vector3& center, number& radius)
			{
				NeuriteProjector::Section cmp(axial);
				std::vector<NeuriteProjector::Section>::const_iterator secIt =
					std::lower_bound(neurite.vSec.begin(), neurite.vSec.end(), cmp,
					                 NeuriteProjector::CompareSections());
				if (secIt == neurite.vSec.end()) secIt = neurite.vSec.end() - 1;
				const number m = secIt->endParam - axial;
				const number* sx = &secIt->splineParamsX[0];
				const number* sy = &secIt->splineParamsY[0];
				const number* sz = &secIt->splineParamsZ[0];
				const number* sr = &secIt->splineParamsR[0];
				center[0] = ((sx[0]*m + sx[1])*m + sx[2])*m + sx[3];
				center[1] = ((sy[0]*m + sy[1])*m + sy[2])*m + sy[3];
				center[2] = ((sz[0]*m + sz[1])*m + sz[2])*m + sz[3];
				radius = ((sr[0]*m + sr[1])*m + sr[2])*m + sr[3];
			};
			std::vector<number> safePositions;
			std::function<void(number, number, size_t)> splitInterval;
			splitInterval = [&] (number a, number b, size_t depth)
			{
				const number mid = 0.5*(a+b);
				const number q1 = 0.5*(a+mid);
				const number q3 = 0.5*(mid+b);
				vector3 pa, pq1, pm, pq3, pb, chordQ1, chordMid, chordQ3;
				number ra, rq1, rm, rq3, rb;
				centerAndRadiusAt(a, pa, ra);
				centerAndRadiusAt(q1, pq1, rq1);
				centerAndRadiusAt(mid, pm, rm);
				centerAndRadiusAt(q3, pq3, rq3);
				centerAndRadiusAt(b, pb, rb);
				VecScaleAdd(chordQ1, 0.75, pa, 0.25, pb);
				VecScaleAdd(chordMid, 0.5, pa, 0.5, pb);
				VecScaleAdd(chordQ3, 0.25, pa, 0.75, pb);
				const number error = std::max(VecDistance(pq1, chordQ1),
					std::max(VecDistance(pm, chordMid), VecDistance(pq3, chordQ3)));
				const number curvedLength = VecDistance(pa, pq1)
					+ VecDistance(pq1, pm) + VecDistance(pm, pq3)
					+ VecDistance(pq3, pb);
				const number thickness =
					std::max<number>(1e-8, (1.0-erScaleFactor)*rm);
				const number bendTargetLength =
					g_adaptiveNephronMeshing.minimumBendAnisotropy
					* 2.0 * PI * rm / static_cast<number>(nOGrid);
				const bool bendDetected = error >
					g_adaptiveNephronMeshing.bendChordErrorFactor * rm;
				const bool bendTooLong = bendDetected
					&& curvedLength > bendTargetLength;
				const bool geometryTooCurved = error > 0.20*thickness;
				if ((bendTooLong || geometryTooCurved) && depth < 12)
				{
					splitInterval(a, mid, depth+1);
					splitInterval(mid, b, depth+1);
				}
				else safePositions.push_back(b);
			};
			number intervalStart = t_start;
			for (size_t i = 0; i < vSegAxPos.size(); ++i)
			{
				splitInterval(intervalStart, vSegAxPos[i], 0);
				intervalStart = vSegAxPos[i];
			}
			if (safePositions.size() != vSegAxPos.size())
				UG_LOGN("Curvature-aware O-grid subdivision: " << vSegAxPos.size()
				        << " -> " << safePositions.size() << " axial segments.");
			vSegAxPos.swap(safePositions);
			nSeg = vSegAxPos.size();
		}

		// add the branching point to segment list (if present)
		if (brit != brit_end)
		{
			vSegAxPos.resize(nSeg+1);
			vSegAxPos[nSeg] = bp_end;
			++nSeg;
		}


		// in case we construct to a BP, find out the branching angle
		// in order to adjust this neurites offset bit by bit
		number addOffset = 0.0;
		size_t child_nid;
		size_t connFaceInd = 0;
		if (brit != brit_end)
		{
			// find branching child neurite
			SmartPtr<NeuriteProjector::BranchingPoint> bp = brit->bp;
			UG_COND_THROW(bp->vNid.size() > 2,
				"This implementation can only handle branching points with one branching child.");

			if (bp->vNid[0] != nid)
				child_nid = bp->vNid[0];
			else
				child_nid = bp->vNid[1];

			// find out branching child neurite initial direction
			vector3 childDir;
			const NeuriteProjector::Section& childSec = vNeurites[child_nid].vSec[0];
			number te = childSec.endParam;

			const number* sp = &childSec.splineParamsX[0];
			number& vc0 = childDir[0];
			vc0 = -3.0*sp[0]*te - 2.0*sp[1];
			vc0 = vc0*te - sp[2];

			sp = &childSec.splineParamsY[0];
			number& vc1 = childDir[1];
			vc1 = -3.0*sp[0]*te - 2.0*sp[1];
			vc1 = vc1*te - sp[2];

			sp = &childSec.splineParamsZ[0];
			number& vc2 = childDir[2];
			vc2 = -3.0*sp[0]*te - 2.0*sp[1];
			vc2 = vc2*te - sp[2];


			// find out neurite direction in next BP
			number bpAxPos = vSegAxPos[nSeg-1];
			size_t tmpSec = curSec;
			for (; tmpSec < nSec; ++tmpSec)
			{
				const NeuriteProjector::Section& sec = neurite.vSec[tmpSec];
				if (sec.endParam >= bpAxPos)
					break;
			}

			const NeuriteProjector::Section& sec = neurite.vSec[tmpSec];
			number monom = sec.endParam - bpAxPos;
			sp = &sec.splineParamsX[0];
			number& v0 = vel[0];
			v0 = -3.0*sp[0]*monom -2.0*sp[1];
			v0 = v0*monom - sp[2];

			sp = &sec.splineParamsY[0];
			number& v1 = vel[1];
			v1 = -3.0*sp[0]*monom -2.0*sp[1];
			v1 = v1*monom - sp[2];

			sp = &sec.splineParamsZ[0];
			number& v2 = vel[2];
			v2 = -3.0*sp[0]*monom -2.0*sp[1];
			v2 = v2*monom - sp[2];

			VecNormalize(vel, vel);

			// calculate offset
			number fac = VecProd(neurite.refDir, vel);
			VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
			VecNormalize(projRefDir, projRefDir);
			VecCross(thirdDir, vel, projRefDir);

			vector2 relCoord;
			VecScaleAppend(childDir, -VecProd(childDir, vel), vel);
			relCoord[0] = VecProd(childDir, projRefDir);
			VecScaleAppend(childDir, -relCoord[0], projRefDir);
			relCoord[1] = VecProd(childDir, thirdDir);
			VecNormalize(relCoord, relCoord);

			number branchOffset = 0.0;
			if (fabs(relCoord[0]) < 1e-8)
				branchOffset = relCoord[1] < 0 ? 1.5*PI : 0.5*PI;
			else
				branchOffset = relCoord[0] < 0 ? PI - atan(-relCoord[1]/relCoord[0]) : atan(relCoord[1]/relCoord[0]);

			addOffset = branchOffset - angleOffset;
			connFaceInd = floor(std::fmod(addOffset+4*PI, 2*PI) / (PI/2));
			addOffset = std::fmod(addOffset - (connFaceInd*PI/2 + PI/4) + 4*PI, 2*PI);
			if (addOffset > PI)
				addOffset -= 2*PI;
			addOffset /= nSeg - 1;
		}

		// create mesh for segments
		Selector sel(g);
		for (size_t s = 0; s < nSeg; ++s)
		{
			vector3 curPos;
			number radius;

			auto computePosVelRad = [&] (number axPos)
			{
				for (; curSec < nSec; ++curSec)
				{
					const NeuriteProjector::Section& sec = neurite.vSec[curSec];
					if (sec.endParam >= axPos)
						break;
				}

				const NeuriteProjector::Section& sec = neurite.vSec[curSec];
				number monom = sec.endParam - axPos;
				const number* sp = &sec.splineParamsX[0];
				number& p0 = curPos[0];
				number& v0 = vel[0];
				p0 = sp[0]*monom + sp[1];
				p0 = p0*monom + sp[2];
				p0 = p0*monom + sp[3];
				v0 = -3.0*sp[0]*monom -2.0*sp[1];
				v0 = v0*monom - sp[2];

				sp = &sec.splineParamsY[0];
				number& p1 = curPos[1];
				number& v1 = vel[1];
				p1 = sp[0]*monom + sp[1];
				p1 = p1*monom + sp[2];
				p1 = p1*monom + sp[3];
				v1 = -3.0*sp[0]*monom -2.0*sp[1];
				v1 = v1*monom - sp[2];

				sp = &sec.splineParamsZ[0];
				number& p2 = curPos[2];
				number& v2 = vel[2];
				p2 = sp[0]*monom + sp[1];
				p2 = p2*monom + sp[2];
				p2 = p2*monom + sp[3];
				v2 = -3.0*sp[0]*monom -2.0*sp[1];
				v2 = v2*monom - sp[2];

				sp = &sec.splineParamsR[0];
				radius = sp[0]*monom + sp[1];
				radius = radius*monom + sp[2];
				radius = radius*monom + sp[3];
			};

			auto computeONB = [&] (number axial)
			{
				if (useMatchingNephronOGrid)
				{
					compute_nephron_parallel_transport_ONB(
						vel, projRefDir, thirdDir, neurite, axial);
					return;
				}
				VecNormalize(vel, vel);

				// calculate reference dir projected to normal plane of velocity
				number fac = VecProd(neurite.refDir, vel);
				VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
				VecNormalize(projRefDir, projRefDir);

				VecCross(thirdDir, vel, projRefDir);
			};


			// usual segment: extrude
			if (s != nSeg - 1 || brit == brit_end)
			{
				// get exact position, velocity, radius and ONB of segment end
				number segAxPos = vSegAxPos[s];
				computePosVelRad(segAxPos);
				computeONB(segAxPos);

				// apply additional offset
				angleOffset = std::fmod(angleOffset + addOffset + 2*PI, 2*PI);

				// extrude from last pos to new pos
				vector3 extrudeDir;
				VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, prevPos);
				std::vector<Volume*> vVol;
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos, EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);

				// set new positions and param attachments
				if (useMatchingNephronOGrid)
				{
					Vertex* center = vVrt[0];
					aaPos[center] = curPos;
					aaSurfParams[center].neuriteID = nid;
					aaSurfParams[center].axial = segAxPos;
					aaSurfParams[center].angular = 0.0;
					aaSurfParams[center].radial = 0.0;
					for (size_t j = 0; j < nOGrid; ++j)
					{
						number angle =
							2.0*PI*((number)j/(number)nOGrid) + angleOffset;
						if (angle >= 2*PI) angle -= 2*PI;
						for (size_t layer = 0; layer < numRadialRings; ++layer)
						{
							const number radial = ringRadial[layer];
							Vertex* v = vVrt[ringBegin[layer] + j];
							vector3 radialVec;
							VecScaleAdd(radialVec,
								radial*radius*cos(angle), projRefDir,
								radial*radius*sin(angle), thirdDir);
							VecAdd(aaPos[v], curPos, radialVec);
							aaSurfParams[v].neuriteID = nid;
							aaSurfParams[v].axial = segAxPos;
							aaSurfParams[v].angular = angle;
							aaSurfParams[v].radial = radial;
						}
					}
				}
				else
				{
					for (size_t j = 0; j < 4; ++j)
					{
						number angle = 0.5*PI*j + angleOffset;
						if (angle > 2*PI) angle -= 2*PI;
						Vertex* v = vVrt[j];
						vector3 radialVec;
						VecScaleAdd(radialVec, erScaleFactor*radius*cos(angle), projRefDir, erScaleFactor*radius*sin(angle), thirdDir);
						VecAdd(aaPos[v], curPos, radialVec);
						aaSurfParams[v].neuriteID = nid;
						aaSurfParams[v].axial = segAxPos;
						aaSurfParams[v].angular = angle;
						aaSurfParams[v].radial = erScaleFactor;
					}
					for (size_t j = 0; j < 12; ++j)
					{
						number angle = PI*((number)j/6) + angleOffset;
						if (angle > 2*PI) angle -= 2*PI;
						Vertex* v = vVrt[j+4];
						vector3 radialVec;
						VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
						VecAdd(aaPos[v], curPos, radialVec);
						aaSurfParams[v].neuriteID = nid;
						aaSurfParams[v].axial = segAxPos;
						aaSurfParams[v].angular = angle;
						aaSurfParams[v].radial = 1.0;
					}
				}

				// ensure correct volume orientation
				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);
			}

			// BP segment: create BP with tetrahedra/pyramids and create whole branch
			else
			{
				auto rotateToBranchInclination = [&] ()
				{
					VecScaleAppend(aaPos[vVrt[(connFaceInd)%4]], erScaleFactor*surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[(connFaceInd+1)%4]], erScaleFactor*surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[(connFaceInd+2)%4]], -erScaleFactor*surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[(connFaceInd+3)%4]], -erScaleFactor*surfBPoffset, vel);
					aaSurfParams[vVrt[(connFaceInd)%4]].axial += erScaleFactor*surfBPoffset / neurite_length;
					aaSurfParams[vVrt[(connFaceInd+1)%4]].axial += erScaleFactor*surfBPoffset / neurite_length;
					aaSurfParams[vVrt[(connFaceInd+2)%4]].axial -= erScaleFactor*surfBPoffset / neurite_length;
					aaSurfParams[vVrt[(connFaceInd+3)%4]].axial -= erScaleFactor*surfBPoffset / neurite_length;

					VecScaleAppend(aaPos[vVrt[4+3*((connFaceInd)%4)]], surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[4+3*((connFaceInd+1)%4)]], surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[4+3*((connFaceInd+2)%4)]], -surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[4+3*((connFaceInd+3)%4)]], -surfBPoffset, vel);
					aaSurfParams[vVrt[4+3*((connFaceInd)%4)]].axial += surfBPoffset / neurite_length;
					aaSurfParams[vVrt[4+3*((connFaceInd+1)%4)]].axial += surfBPoffset / neurite_length;
					aaSurfParams[vVrt[4+3*((connFaceInd+2)%4)]].axial -= surfBPoffset / neurite_length;
					aaSurfParams[vVrt[4+3*((connFaceInd+3)%4)]].axial -= surfBPoffset / neurite_length;

					VecScaleAppend(aaPos[vVrt[5+3*((connFaceInd)%4)]], 1.366*surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[6+3*((connFaceInd)%4)]], 1.366*surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[5+3*((connFaceInd+2)%4)]], -1.366*surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[6+3*((connFaceInd+2)%4)]], -1.366*surfBPoffset, vel);
					aaSurfParams[vVrt[5+3*((connFaceInd)%4)]].axial += 1.366*surfBPoffset / neurite_length;
					aaSurfParams[vVrt[6+3*((connFaceInd)%4)]].axial += 1.366*surfBPoffset / neurite_length;
					aaSurfParams[vVrt[5+3*((connFaceInd+2)%4)]].axial -= 1.366*surfBPoffset / neurite_length;
					aaSurfParams[vVrt[6+3*((connFaceInd+2)%4)]].axial -= 1.366*surfBPoffset / neurite_length;

					VecScaleAppend(aaPos[vVrt[5+3*((connFaceInd+1)%4)]], 0.366*surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[6+3*((connFaceInd+1)%4)]], -0.366*surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[5+3*((connFaceInd+3)%4)]], -0.366*surfBPoffset, vel);
					VecScaleAppend(aaPos[vVrt[6+3*((connFaceInd+3)%4)]], 0.366*surfBPoffset, vel);
					aaSurfParams[vVrt[5+3*((connFaceInd+1)%4)]].axial += 0.366*surfBPoffset / neurite_length;
					aaSurfParams[vVrt[6+3*((connFaceInd+1)%4)]].axial -= 0.366*surfBPoffset / neurite_length;
					aaSurfParams[vVrt[5+3*((connFaceInd+3)%4)]].axial -= 0.366*surfBPoffset / neurite_length;
					aaSurfParams[vVrt[6+3*((connFaceInd+3)%4)]].axial += 0.366*surfBPoffset / neurite_length;
				};


				std::vector<Volume*> vBPVols;
				vBPVols.reserve(27);

				// correct vertex offsets to reflect angle at which child branches
				rotateToBranchInclination();

				// add branch neurite ID to its initial vertices
				const uint32 brId = ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[4+3*((connFaceInd)%4)]].neuriteID |= brId;
				aaSurfParams[vVrt[4+3*((connFaceInd+1)%4)]].neuriteID |= brId;
				aaSurfParams[vVrt[5+3*((connFaceInd)%4)]].neuriteID |= brId;
				aaSurfParams[vVrt[6+3*((connFaceInd)%4)]].neuriteID |= brId;


				// prepare branch vertices
				std::vector<Vertex*> vBranchVrts(16);
				vBranchVrts[4] = vVrt[4+3*((connFaceInd+1)%4)];
				vBranchVrts[13] = vVrt[4+3*((connFaceInd)%4)];
				vBranchVrts[14] = vVrt[5+3*((connFaceInd)%4)];
				vBranchVrts[15] = vVrt[6+3*((connFaceInd)%4)];


				// extrude to first third of BP //
				number segAxPos = 0.5*(1.0 + erScaleFactor)*vSegAxPos[s-1] + 0.5*(1.0 - erScaleFactor)*vSegAxPos[s];
				computePosVelRad(segAxPos);
				computeONB(segAxPos);

				vector3 extrudeDir;
				VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, prevPos);
				std::vector<Volume*> vVol;
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos, EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);
				for (size_t j = 0; j < 9; ++j)
					vBPVols.push_back(vVol[j]);

				// set new positions and param attachments
				for (size_t j = 0; j < 4; ++j)
				{
					number angle = 0.5*PI*j + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					Vertex* v = vVrt[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, erScaleFactor*radius*cos(angle), projRefDir, erScaleFactor*radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = erScaleFactor;
				}
				for (size_t j = 0; j < 12; ++j)
				{
					number angle = PI*((number)j/6) + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					Vertex* v = vVrt[j+4];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// correct vertex offsets to reflect angle at which child branches
				rotateToBranchInclination();

				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);

				// add branch neurite ID to its initial vertices
				aaSurfParams[vVrt[(connFaceInd)%4]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[(connFaceInd+1)%4]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);

				aaSurfParams[vVrt[4+3*((connFaceInd)%4)]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[4+3*((connFaceInd+1)%4)]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[5+3*((connFaceInd)%4)]].neuriteID = child_nid;
				aaSurfParams[vVrt[6+3*((connFaceInd)%4)]].neuriteID = child_nid;

				// prepare branch vertices
				vBranchVrts[0] = vVrt[6+3*((connFaceInd)%4)];
				vBranchVrts[3] = vVrt[5+3*((connFaceInd)%4)];
				vBranchVrts[5] = vVrt[4+3*((connFaceInd+1)%4)];
				vBranchVrts[12] = vVrt[4+3*((connFaceInd)%4)];


				// extrude to second third of BP //
				segAxPos = 0.5*(1.0 - erScaleFactor)*vSegAxPos[s-1] + 0.5*(1.0 + erScaleFactor)*vSegAxPos[s];
				prevPos = curPos;
				computePosVelRad(segAxPos);
				computeONB(segAxPos);

				VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, prevPos);
				vVol.clear();
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos, EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);
				for (size_t j = 0; j < 9; ++j)
					vBPVols.push_back(vVol[j]);

				// set new positions and param attachments
				for (size_t j = 0; j < 4; ++j)
				{
					number angle = 0.5*PI*j + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					Vertex* v = vVrt[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, erScaleFactor*radius*cos(angle), projRefDir, erScaleFactor*radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = erScaleFactor;
				}
				for (size_t j = 0; j < 12; ++j)
				{
					number angle = PI*((number)j/6) + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					Vertex* v = vVrt[j+4];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// correct vertex offsets to reflect angle at which child branches
				rotateToBranchInclination();

				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);

				// add branch neurite ID to its initial vertices
				aaSurfParams[vVrt[(connFaceInd)%4]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[(connFaceInd+1)%4]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);

				aaSurfParams[vVrt[4+3*((connFaceInd)%4)]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[4+3*((connFaceInd+1)%4)]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[5+3*((connFaceInd)%4)]].neuriteID = child_nid;
				aaSurfParams[vVrt[6+3*((connFaceInd)%4)]].neuriteID = child_nid;

				// prepare branch vertices
				vBranchVrts[1] = vVrt[6+3*((connFaceInd)%4)];
				vBranchVrts[2] = vVrt[5+3*((connFaceInd)%4)];
				vBranchVrts[6] = vVrt[4+3*((connFaceInd+1)%4)];
				vBranchVrts[11] = vVrt[4+3*((connFaceInd)%4)];

				// extrude to end of BP //
				segAxPos = vSegAxPos[s];
				prevPos = curPos;
				computePosVelRad(segAxPos);
				computeONB(segAxPos);

				VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, prevPos);
				vVol.clear();
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos, EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);
				for (size_t j = 0; j < 9; ++j)
					vBPVols.push_back(vVol[j]);

				// set new positions and param attachments
				for (size_t j = 0; j < 4; ++j)
				{
					number angle = 0.5*PI*j + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					Vertex* v = vVrt[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, erScaleFactor*radius*cos(angle), projRefDir, erScaleFactor*radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = erScaleFactor;
				}
				for (size_t j = 0; j < 12; ++j)
				{
					number angle = PI*((number)j/6) + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					Vertex* v = vVrt[j+4];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// correct vertex offsets to reflect angle at which child branches
				rotateToBranchInclination();

				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);

				// add branch neurite ID to its initial vertices
				aaSurfParams[vVrt[4+3*((connFaceInd)%4)]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[4+3*((connFaceInd+1)%4)]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[5+3*((connFaceInd)%4)]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[6+3*((connFaceInd)%4)]].neuriteID += ((brit - vBR.begin()) << 20) + (1 << 28);

				// prepare branch vertices
				vBranchVrts[7] = vVrt[4+3*((connFaceInd+1)%4)];
				vBranchVrts[8] = vVrt[6+3*((connFaceInd)%4)];
				vBranchVrts[9] = vVrt[5+3*((connFaceInd)%4)];
				vBranchVrts[10] = vVrt[4+3*((connFaceInd)%4)];


				// prepare connectingEdges, -Faces for branch
				std::vector<EdgeDescriptor> vED(24);
				for (size_t i = 0; i < 4; ++i)
				{
					vED[i] = EdgeDescriptor(vBranchVrts[i], vBranchVrts[(i+1)%4]);
					vED[i+4] = EdgeDescriptor(vBranchVrts[i], vBranchVrts[5+3*i]);
					vED[i+8] = EdgeDescriptor(vBranchVrts[(i+1)%4], vBranchVrts[6+3*i]);
				}
				for (size_t i = 0; i < 12; ++i)
					vED[i+12] = EdgeDescriptor(vBranchVrts[i+4], vBranchVrts[(i+1)%12+4]);

				std::vector<FaceDescriptor> vFD(9);
				vFD[0] = FaceDescriptor(vBranchVrts[0], vBranchVrts[1], vBranchVrts[2], vBranchVrts[3]);
				for (size_t i = 0; i < 4; ++i)
				{
					vFD[i+1] = FaceDescriptor(vBranchVrts[i], vBranchVrts[(3*i+11)%12+4], vBranchVrts[3*i+4], vBranchVrts[3*i+5]);
					vFD[i+5] = FaceDescriptor(vBranchVrts[i], vBranchVrts[3*i+5], vBranchVrts[3*i+6], vBranchVrts[(i+1)%4]);
				}

				typedef Grid::traits<Face>::secure_container faceCont;
				std::vector<Face*> vBranchFaces(9);
				for (size_t j = 0; j < 9; ++j)
				{
					const FaceDescriptor& qDesc = vFD[j];

					faceCont fl;
					bool found = false;
					for (size_t k = 0; k < 27; ++k)
					{
						g.associated_elements(fl, vBPVols[k]);
						const size_t flSz = fl.size();
						for (size_t f = 0; f < flSz; ++f)
						{
							if (CompareVertices(fl[f], &qDesc))
							{
								vBranchFaces[j] = fl[f];
								found = true;
								break;
							}
						}
						if (found)
							break;
					}
					UG_COND_THROW(!found, "Connecting face " << j << " not found.")
				}

				typedef Grid::traits<Edge>::secure_container edgeCont;
				std::vector<Edge*> vBranchEdges(24);
				for (size_t j = 0; j < 24; ++j)
				{
					const EdgeDescriptor& eDesc = vED[j];

					edgeCont el;
					bool found = false;
					for (size_t k = 0; k < 9; ++k)
					{
						g.associated_elements(el, vBranchFaces[k]);
						const size_t elSz = el.size();
						for (size_t e = 0; e < elSz; ++e)
						{
							if (CompareVertices(el[e], &eDesc))
							{
								vBranchEdges[j] = el[e];
								found = true;
								break;
							}
						}
						if (found)
							break;
					}
					UG_COND_THROW(!found, "Connecting edge " << j << " not found.")
				}

				// correct connecting volume subset indices
				{
					Volume* connVol = vBPVols[connFaceInd+14];
					sh.assign_subset(connVol, 1);

					typedef Grid::traits<Face>::secure_container faceCont;
					faceCont fl;
					g.associated_elements(fl, connVol);
					const size_t flSz = fl.size();
					for (size_t f = 0; f < flSz; ++f)
					{
						Face* sideFace = fl[f];
						Volume* opp = GetConnectedNeighbor(g, sideFace, connVol);
						if (!opp || sh.get_subset_index(opp) == 1)
							sh.assign_subset(sideFace, 1);
						else
						{
							sh.assign_subset(sideFace, 3);

							typedef Grid::traits<Edge>::secure_container edgeCont;
							edgeCont el;
							g.associated_elements(el, sideFace);
							const size_t elSz = el.size();
							for (size_t e = 0; e < elSz; ++e)
							{
								Edge* sideEdge = el[e];
								sh.assign_subset(sideEdge, 3);
								sh.assign_subset(sideEdge->vertex(0), 3);
								sh.assign_subset(sideEdge->vertex(1), 3);
							}
						}
					}
				}
				for (size_t j = 1; j < 9; ++j)
				{
					Face* connFace = vBranchFaces[j];
					sh.assign_subset(connFace, 0);
				}
				for (size_t j = 4; j < 12; ++j)
				{
					Edge* connEdge = vBranchEdges[j];
					sh.assign_subset(connEdge, 0);
				}


				// recursively build branch
				create_neurite_with_er(vNeurites, vPos, vR, child_nid, erScaleFactor, anisotropy,
					nephronOGridVertices, coarseLumenCenterOnly,
					g, aaPos, aaSurfParams, sh, &vBranchVrts, &vBranchEdges, &vBranchFaces, branchOffset[1]);
			}

			prevPos = curPos;
		}

		// update t_end and curSec
		if (brit != brit_end)
			t_end = bp_end;

		for (; curSec < nSec; ++curSec)
		{
			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			if (sec.endParam >= t_end)
				break;
		}

		// check whether tip has been reached
		if (brit == brit_end)
			break;
		else
			++brit;
	}
}



static void create_neurite_1d
(
    const std::vector<NeuriteProjector::Neurite>& vNeurites,
    const std::vector<std::vector<vector3> >& vPos,
    const std::vector<std::vector<number> >& vR,
    size_t nid,
	number anisotropy,
    Grid& g,
    Grid::VertexAttachmentAccessor<APosition>& aaPos,
    Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
	Grid::VertexAttachmentAccessor<Attachment<number> >& aaDiam,
    Vertex* connectingVrt = NULL
)
{
    const NeuriteProjector::Neurite& neurite = vNeurites[nid];
    const std::vector<vector3>& pos = vPos[nid];
    const std::vector<number>& r = vR[nid];

    number neurite_length = 0.0;
    for (size_t i = 1; i < pos.size(); ++i)
    	neurite_length += VecDistance(pos[i], pos[i-1]);

    size_t nSec = neurite.vSec.size();

    const std::vector<NeuriteProjector::BranchingRegion>& vBR = neurite.vBR;
    std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit = vBR.begin();
    std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit_end = vBR.end();

    if (connectingVrt)
    {
        // ignore first branching region (the connecting region)
        ++brit;
    }
    else
    {
        // create first vertex
    	connectingVrt = *g.create<RegularVertex>();
		aaPos[connectingVrt] = pos[0];
		aaDiam[connectingVrt] = r[0];

		aaSurfParams[connectingVrt].neuriteID = nid;
		aaSurfParams[connectingVrt].axial = 0.0;
		aaSurfParams[connectingVrt].angular = 0.0;
		aaSurfParams[connectingVrt].radial = 0.0;
    }

    // Now create dendrite to the next branching point and iterate this process.
    // We want to create each of the segments with approx. the same aspect ratio.
    // To that end, we first calculate the length of the section to be created (in units of radius)
    // and then divide this number by 2^n where n is the number of anisotropic refinements to
    // be performed to make all segments (more or less) isotropic. The result is the number
    // of segments to be used for the section.
    number t_start = 0.0;
    number t_end = 0.0;
    size_t curSec = 0;

    while (true)
    {
    	t_start = t_end;

    	// last section: create until tip
    	if (brit == brit_end)
    		t_end = 1.0;

    	// otherwise: section goes to next branching point
    	else
    		t_end = brit->t;

    	// calculate total length in units of radius
    	// = integral from t_start to t_end over: ||v(t)|| / r(t) dt
    	number lengthOverRadius = calculate_length_over_radius(t_start, t_end, neurite, curSec);

    	// to reach the desired anisotropy on the surface in the refinement limit,
		// it has to be multiplied by pi/2 h
		size_t nSeg = (size_t) floor(lengthOverRadius / (anisotropy*0.5*PI));
    	if (nSeg < 1 || lengthOverRadius < 0)
    		nSeg = 1;
    	number segLength = lengthOverRadius / nSeg;	// segments are between 8 and 16 radii long
    	std::vector<number> vSegAxPos(nSeg);
    	calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neurite, curSec, segLength);

    	// create mesh for segments
    	Selector sel(g);
    	for (size_t s = 0; s < nSeg; ++s)
    	{
    		// get exact position and radius of segment end
    		number segAxPos = vSegAxPos[s];
    		for (; curSec < nSec; ++curSec)
    		{
				const NeuriteProjector::Section& sec = neurite.vSec[curSec];
				if (sec.endParam >= segAxPos)
					break;
    		}

    		const NeuriteProjector::Section& sec = neurite.vSec[curSec];
    		vector3 curPos;
    		number monom = sec.endParam - segAxPos;
    		const number* sp = &sec.splineParamsX[0];
			number& p0 = curPos[0];
			p0 = sp[0]*monom + sp[1];
			p0 = p0*monom + sp[2];
			p0 = p0*monom + sp[3];

    		sp = &sec.splineParamsY[0];
			number& p1 = curPos[1];
			p1 = sp[0]*monom + sp[1];
			p1 = p1*monom + sp[2];
			p1 = p1*monom + sp[3];

    		sp = &sec.splineParamsZ[0];
			number& p2 = curPos[2];
			p2 = sp[0]*monom + sp[1];
			p2 = p2*monom + sp[2];
			p2 = p2*monom + sp[3];

    		number curRad;
    		sp = &sec.splineParamsR[0];
    		curRad = sp[0]*monom + sp[1];
    		curRad = curRad*monom + sp[2];
    		curRad = curRad*monom + sp[3];


			// create new vertex and connect with edge
			Vertex* newVrt = *g.create<RegularVertex>();
			*g.create<RegularEdge>(EdgeDescriptor(connectingVrt, newVrt));

			// position and diameter
			aaPos[newVrt] = curPos;
			aaDiam[newVrt] = 2*curRad;

			// set new param attachments
			aaSurfParams[newVrt].neuriteID = nid;
			aaSurfParams[newVrt].axial = segAxPos;
			aaSurfParams[newVrt].angular = 0.0;
			aaSurfParams[newVrt].radial = 0.0;

			// update
			connectingVrt = newVrt;
    	}

    	// connect branching neurites if present
		if (brit != brit_end)
		{
			// find branching child neurite
			SmartPtr<NeuriteProjector::BranchingPoint> bp = brit->bp;
			UG_COND_THROW(bp->vNid.size() > 2,
				"This implementation can only handle branching points with one branching child.");

			size_t child_nid;
			if (bp->vNid[0] != nid)
				child_nid = bp->vNid[0];
			else
				child_nid = bp->vNid[1];


			// add branch neurite ID to its initial vertex
			aaSurfParams[connectingVrt].neuriteID += (brit - vBR.begin()) << 20;  // add branching region index
			aaSurfParams[connectingVrt].neuriteID += 1 << 28;  // add child ID (always 0, since there can only be one child here)

			create_neurite_1d(vNeurites, vPos, vR, child_nid, anisotropy,
				g, aaPos, aaSurfParams, aaDiam, connectingVrt);
		}

    	// update curSec
    	for (; curSec < nSec; ++curSec)
    	{
			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			if (sec.endParam >= t_end)
				break;
    	}

    	// check whether tip has been reached
    	if (brit == brit_end)
    		break;
    	else
			++brit;
    }
}



void export_to_swc
(
	Grid& g,
	ISubsetHandler& sh,
	const std::string& fileName
)
{
	// get access to positions
	UG_COND_THROW(!g.has_vertex_attachment(aPosition), "Position attachment not attached to grid.")
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

	// get access to diameter attachment
	ANumber aDiam = GlobalAttachments::attachment<ANumber>("diameter");
	UG_COND_THROW(!g.has_vertex_attachment(aDiam), "No diameter attachment attached to grid.");
	Grid::AttachmentAccessor<Vertex, ANumber> aaDiam(g, aDiam);

    // analyze subset names to find out corresponding swc-types
    size_t nss = sh.num_subsets();
    std::vector<size_t> vType(nss);
    bool soma_subset_present = false;
    for (size_t i = 0; i < nss; ++i)
    {
        std::string name(sh.get_subset_name(i));
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        if (name.find("SOMA") != std::string::npos)
        {
            soma_subset_present = true;
            vType[i] = 1;
        }
        else if (name.find("AXON") != std::string::npos)
            vType[i] = 2;
        else if (name.find("APIC") != std::string::npos)
            vType[i] = 4;
        else if (name.find("DEND") != std::string::npos)
            vType[i] = 3;
        else vType[i] = 0;
    }

    if (!soma_subset_present)
        UG_LOGN("Warning: No somatic subset could be identified.")

	if (g.begin<Vertex>() == g.end<Vertex>())
	{
		UG_LOGN("Warning: No vertices contained in grid.")
		return;
	}

    // find soma vertex (if identifiable)
    Vertex* start = *g.begin<Vertex>();
    if (soma_subset_present)
    {
        g.begin_marking();
        std::queue<Vertex*> q; // corresponds to breadth-first
        q.push(start);
        while (!q.empty())
        {
            Vertex* v = q.front();
            if (vType[sh.get_subset_index(v)] == 1) break;
            g.mark(v);
            q.pop();

            // push neighboring elems to queue
            Grid::traits<Edge>::secure_container edges;
            g.associated_elements(edges, v);

            size_t sz = edges.size();
            for (size_t e = 0; e < sz; ++e)
            {
                Vertex* otherEnd = GetOpposingSide(g, edges[e], v);
                if (!g.is_marked(otherEnd))
                    q.push(otherEnd);
            }
        }
        g.end_marking();

        if (q.empty())
            UG_LOGN("Warning: No soma vertex could be found in the requested neuron.")
        else
            start = q.front();
    }

    // write the neuron to file
    std::ofstream outFile(fileName.c_str(), std::ios::out);
    UG_COND_THROW(!outFile.is_open(), "Could not open output file '" << fileName << "'.");

    outFile << "# This file has been generated by UG4." << std::endl;

    std::stack<std::pair<Vertex*, int> > stack; // corresponds to depth-first
    stack.push(std::make_pair(start, -1));

    g.begin_marking();
    int ind = 0;   // by convention, swc starts with index 1
    bool all_types_identified = true;
    while (!stack.empty())
    {
        // get all infos regarding vertex
        std::pair<Vertex*, int>& info = stack.top();
        Vertex* v = info.first;
        int conn = info.second;
        stack.pop();

        // mark curr vrt
        g.mark(v);

        size_t type = vType[sh.get_subset_index(v)];
        if (!type) all_types_identified = false;

        const Domain3d::position_type& coord = aaPos[v];

        number radius = 0.5*aaDiam[v];

        // write line to file
        outFile << ++ind << " " << type << " "
            << coord[0] << " " << coord[1] << " " << coord[2] << " "
            << radius << " " << conn << std::endl;

        // push neighboring elems to queue
        Grid::traits<Edge>::secure_container edges;
        g.associated_elements(edges, v);

        size_t sz = edges.size();
        for (size_t e = 0; e < sz; ++e)
        {
            Vertex* otherEnd = GetOpposingSide(g, edges[e], v);
            if (!g.is_marked(otherEnd))
                stack.push(std::make_pair(otherEnd, ind));
        }
    }
    g.end_marking();

    if (!all_types_identified)
        UG_LOGN("WARNING: Some vertex type(s) - soma, dendrite, axon, etc. -\n"
            "could not be identified using the subset names.\n"
            << "To ensure correct types in the resulting swc file, the ugx subset names\n"
            "need to contain one of the strings \"SOMA\", \"AXON\", \"DEND\", \"APIC\",\n"
            "upper/lower case can be ignored.");

    outFile.close();
}


void swc_points_to_grid
(
	const std::vector<swc_types::SWCPoint>& vPts,
	Grid& g,
	SubsetHandler& sh,
	number scale_length = 1.0
)
{
	if (!g.has_vertex_attachment(aPosition))
		g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

	ANumber aDiam = GlobalAttachments::attachment<ANumber>("diameter");
	if (!g.has_vertex_attachment(aDiam))
		g.attach_to_vertices(aDiam);
	Grid::AttachmentAccessor<Vertex, ANumber> aaDiam(g, aDiam);

	// create grid
	const size_t nP = vPts.size();
	std::vector<Vertex*> vrts(nP, NULL);
	for (size_t i = 0; i < nP; ++i)
	{
		const swc_types::SWCPoint& pt = vPts[i];

		// create vertex and save
		Vertex* v = vrts[i] = *g.create<RegularVertex>();
		VecScale(aaPos[v], pt.coords, scale_length);
		sh.assign_subset(v, pt.type - 1);
		aaDiam[v] = 2 * pt.radius * scale_length;

		// create edge connections to already created vertices
		for (size_t j = 0; j < pt.conns.size(); ++j)
		{
			if (pt.conns[j] < i)
			{
				Edge* e = *g.create<RegularEdge>(EdgeDescriptor(vrts[pt.conns[j]], v));
				sh.assign_subset(e, vPts[pt.conns[j]].type - 1);
			}
		}
	}

	// final subset managment
	AssignSubsetColors(sh);
	sh.set_subset_name("soma", 0);
	sh.set_subset_name("axon", 1);
	sh.set_subset_name("dend", 2);
	sh.set_subset_name("apic", 3);
	EraseEmptySubsets(sh);
}



static void create_padded_box_surface
(
	Grid& g,
	SubsetHandler& sh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	number padding,
	int subsetIndex,
	size_t numSurfaceRefs = 0,
	const OuterWallTerminalMetadata* terminalMetadata = NULL,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >*
		aaSurfParams = NULL
)
{
	UG_COND_THROW(padding <= 0.0, "Box padding must be positive.");
	UG_COND_THROW(g.num<Vertex>() == 0, "Cannot create a box around an empty grid.");
	UG_COND_THROW(numSurfaceRefs > 10,
	              "Outer-box surface refinement is unreasonably large.");

	const bool terminalPorts = terminalMetadata
		&& terminalMetadata->terminalsOnOuterWall;
	vector3 boxMin;
	vector3 boxMax;
	if (terminalPorts)
	{
		UG_COND_THROW(!terminalMetadata->hasFixedBox,
		              "Terminal ports require fixed OuterWall bounds.");
		boxMin = terminalMetadata->boxMin;
		boxMax = terminalMetadata->boxMax;
	}
	else
	{
		VertexIterator vit = g.begin<Vertex>();
		boxMin = aaPos[*vit];
		boxMax = boxMin;
		for (; vit != g.end<Vertex>(); ++vit)
		{
			const vector3& p = aaPos[*vit];
			for (size_t d = 0; d < 3; ++d)
			{
				boxMin[d] = std::min(boxMin[d], p[d]);
				boxMax[d] = std::max(boxMax[d], p[d]);
			}
		}
		for (size_t d = 0; d < 3; ++d)
		{
			boxMin[d] -= padding;
			boxMax[d] += padding;
		}
	}

	/*
	 * Build the empty box as a conforming structured surface grid before
	 * TetGen sees it. numSurfaceRefs controls the regular subdivision of the
	 * original six-quadrilateral shell. No Inter volume exists yet, so only
	 * OuterWall is refined.
	 */
	const size_t n = static_cast<size_t>(1) << numSurfaceRefs;
	typedef std::array<size_t, 3> BoxKey;
	std::map<BoxKey, Vertex*> vertices;
	std::set<std::pair<Vertex*, Vertex*> > edges;
	std::vector<Edge*> portalEdges[6];
	std::vector<size_t> portalPlanes;
	if (terminalPorts)
	{
		UG_COND_THROW(!aaSurfParams,
		              "Terminal ports require npSurfParams for exact cap identification.");
		number diagonal = 0.0;
		for (size_t d = 0; d < 3; ++d)
			diagonal += (boxMax[d] - boxMin[d]) * (boxMax[d] - boxMin[d]);
		diagonal = std::sqrt(diagonal);
		// Identify terminal sections from the projector coordinates, not from a
		// geometric tolerance. This includes the center, lumen ring, Apical ring,
		// membrane-middle ring, and Basolateral rim of both ends.
		struct AxialRange
		{
			number minimum;
			number maximum;
			AxialRange() : minimum(std::numeric_limits<number>::max()),
			               maximum(-std::numeric_limits<number>::max()) {}
		};
		std::map<uint32,AxialRange> axialRanges;
		for (VertexIterator vit = g.begin<Vertex>(); vit != g.end<Vertex>(); ++vit)
		{
			const NeuriteProjector::SurfaceParams& sp = (*aaSurfParams)[*vit];
			const uint32 neuriteID = sp.neuriteID & ((1u << 20) - 1u);
			AxialRange& range = axialRanges[neuriteID];
			range.minimum = std::min(range.minimum, (number)sp.axial);
			range.maximum = std::max(range.maximum, (number)sp.axial);
		}
		typedef std::pair<uint32,int> TerminalKey;
		std::map<TerminalKey,std::vector<Vertex*> > terminalVertices;
		const number axialTolerance = 1e-7;
		for (VertexIterator vit = g.begin<Vertex>(); vit != g.end<Vertex>(); ++vit)
		{
			const NeuriteProjector::SurfaceParams& sp = (*aaSurfParams)[*vit];
			const uint32 neuriteID = sp.neuriteID & ((1u << 20) - 1u);
			const AxialRange& range = axialRanges[neuriteID];
			if (std::fabs((number)sp.axial - range.minimum) <= axialTolerance)
				terminalVertices[TerminalKey(neuriteID,0)].push_back(*vit);
			if (std::fabs((number)sp.axial - range.maximum) <= axialTolerance)
				terminalVertices[TerminalKey(neuriteID,1)].push_back(*vit);
		}
		size_t snappedVertices = 0;
		for (std::map<TerminalKey,std::vector<Vertex*> >::iterator groupIt =
			 terminalVertices.begin(); groupIt != terminalVertices.end(); ++groupIt)
		{
			std::vector<Vertex*>& group = groupIt->second;
			UG_COND_THROW(group.empty(), "An identified terminal section is empty.");
			Vertex* center = group[0];
			for (size_t i = 1; i < group.size(); ++i)
				if (std::fabs((number)(*aaSurfParams)[group[i]].radial)
					< std::fabs((number)(*aaSurfParams)[center].radial))
					center = group[i];
			number nearestDistance = std::numeric_limits<number>::max();
			size_t nearestAxis = 0;
			size_t nearestSide = 0;
			for (size_t axis = 0; axis < 3; ++axis)
				for (size_t side = 0; side < 2; ++side)
				{
					const number plane = side ? boxMax[axis] : boxMin[axis];
					const number distance = std::fabs(aaPos[center][axis] - plane);
					if (distance < nearestDistance)
					{
						nearestDistance = distance;
						nearestAxis = axis;
						nearestSide = side;
					}
				}
			// The coarse generator intentionally offsets its first cross-section
			// from the root SWC sample to avoid a short initial cell. A sub-cell
			// correction is therefore expected at the root; the long straight
			// wall-normal tail makes this snap shape-preserving.
			const number numericalLengthTolerance =
				std::numeric_limits<number>::epsilon() * diagonal * 64.0;
			const number maximumSnapDistance = std::max(
				numericalLengthTolerance, diagonal * 2e-2);
			UG_COND_THROW(nearestDistance > maximumSnapDistance,
			              "Terminal center is " << nearestDistance
			              << " from its nearest OuterWall plane; extension metadata and spline disagree.");
			const number plane = nearestSide ? boxMax[nearestAxis] : boxMin[nearestAxis];
			// The root cross-section may be generated a short distance away from the
			// first SWC sample.  Center the actual terminal section—not only the SWC
			// endpoint—in its structured OuterWall cell.  This keeps the circular rim
			// away from cell corners and gives the local constrained triangulation a
			// balanced transition on both the root and tip sides.
			vector3 tangentialShift;
			VecSet(tangentialShift, 0.0);
			for (size_t dimension = 0; dimension < 3; ++dimension)
			{
				if (dimension == nearestAxis) continue;
				const number extent = boxMax[dimension] - boxMin[dimension];
				const number scaled = (aaPos[center][dimension] - boxMin[dimension])
					* static_cast<number>(n) / extent;
				const size_t cell = std::min(n-1, static_cast<size_t>(std::max(
					(number)0.0, std::floor(scaled))));
				const number cellCenter = boxMin[dimension]
					+ (static_cast<number>(cell) + 0.5) * extent
					/ static_cast<number>(n);
				tangentialShift[dimension] = cellCenter - aaPos[center][dimension];
			}
			for (size_t i = 0; i < group.size(); ++i)
			{
				VecAdd(aaPos[group[i]], aaPos[group[i]], tangentialShift);
				aaPos[group[i]][nearestAxis] = plane;
				++snappedVertices;
			}
			UG_LOGN("Snapped terminal neurite=" << groupIt->first.first
			        << " end=" << groupIt->first.second << " (" << group.size()
			        << " vertices) to " << "xyz"[nearestAxis]
			        << (nearestSide ? "max" : "min")
			        << "; pre-snap center distance=" << nearestDistance
			        << ", tangential center correction="
			        << std::sqrt(VecLengthSq(tangentialShift)) << ".");
		}
		const number planeTolerance = std::max(
			std::numeric_limits<number>::epsilon() * diagonal * 64.0,
			diagonal * 1e-10);
		UG_LOGN("Snapped " << snappedVertices
		        << " terminal-section vertices exactly onto OuterWall planes.");
		for (EdgeIterator eit = g.begin<Edge>(); eit != g.end<Edge>(); ++eit)
		{
			Edge* edge = *eit;
			if (sh.get_subset_index(edge) != 2) continue;
			for (size_t axis = 0; axis < 3; ++axis)
				for (size_t side = 0; side < 2; ++side)
				{
					const number plane = side ? boxMax[axis] : boxMin[axis];
					if (std::fabs(aaPos[edge->vertex(0)][axis] - plane) <= planeTolerance
						&& std::fabs(aaPos[edge->vertex(1)][axis] - plane) <= planeTolerance)
						portalEdges[2*axis + side].push_back(edge);
				}
		}
		for (size_t plane = 0; plane < 6; ++plane)
			if (!portalEdges[plane].empty()) portalPlanes.push_back(plane);
		if (portalPlanes.empty())
		{
			number nearest = std::numeric_limits<number>::max();
			size_t subsetTwoEdges = 0;
			for (EdgeIterator eit = g.begin<Edge>(); eit != g.end<Edge>(); ++eit)
			{
				Edge* edge = *eit;
				if (sh.get_subset_index(edge) != 2) continue;
				++subsetTwoEdges;
				for (size_t axis = 0; axis < 3; ++axis)
					for (size_t side = 0; side < 2; ++side)
					{
						const number plane = side ? boxMax[axis] : boxMin[axis];
						const number edgeDistance = std::max(
							std::fabs(aaPos[edge->vertex(0)][axis] - plane),
							std::fabs(aaPos[edge->vertex(1)][axis] - plane));
						nearest = std::min(nearest, edgeDistance);
					}
			}
			UG_LOGN("Terminal rim diagnostic: subset-2 edges=" << subsetTwoEdges
			        << ", nearest common-plane distance=" << nearest
			        << ", detection tolerance=" << planeTolerance << ".");
		}
		UG_COND_THROW(portalPlanes.empty(),
		              "No terminal Basolateral rim lies on the fixed OuterWall.");
		size_t rimEdgeCount = 0;
		for (size_t i = 0; i < portalPlanes.size(); ++i)
			rimEdgeCount += portalEdges[portalPlanes[i]].size();
		UG_LOGN("Detected " << rimEdgeCount
		        << " terminal rim edges on " << portalPlanes.size()
		        << " OuterWall plane(s).");
	}

	auto vertexAt = [&] (size_t ix, size_t iy, size_t iz) -> Vertex*
	{
		const BoxKey key = {{ix, iy, iz}};
		std::map<BoxKey, Vertex*>::iterator found = vertices.find(key);
		if (found != vertices.end()) return found->second;

		Vertex* v = *g.create<RegularVertex>();
		const number fx = static_cast<number>(ix) / static_cast<number>(n);
		const number fy = static_cast<number>(iy) / static_cast<number>(n);
		const number fz = static_cast<number>(iz) / static_cast<number>(n);
		aaPos[v] = vector3(
			boxMin[0] + fx * (boxMax[0] - boxMin[0]),
			boxMin[1] + fy * (boxMax[1] - boxMin[1]),
			boxMin[2] + fz * (boxMax[2] - boxMin[2]));
		sh.assign_subset(v, subsetIndex);
		vertices[key] = v;
		return v;
	};

	auto createEdge = [&] (Vertex* a, Vertex* b)
	{
		if (std::less<Vertex*>()(b, a)) std::swap(a, b);
		const std::pair<Vertex*, Vertex*> key(a, b);
		if (!edges.insert(key).second) return;
		Edge* edge = g.get_edge(a, b);
		if (!edge) edge = *g.create<RegularEdge>(EdgeDescriptor(a, b));
		sh.assign_subset(edge, subsetIndex);
	};

	auto createQuad = [&] (Vertex* a, Vertex* b, Vertex* c, Vertex* d)
	{
		createEdge(a, b);
		createEdge(b, c);
		createEdge(c, d);
		createEdge(d, a);
		createEdge(a, c);
		Face* first = *g.create<Triangle>(TriangleDescriptor(a, b, c));
		Face* second = *g.create<Triangle>(TriangleDescriptor(a, c, d));
		sh.assign_subset(first, subsetIndex);
		sh.assign_subset(second, subsetIndex);
	};

	auto hasPortal = [&] (size_t plane) -> bool
	{
		return !portalEdges[plane].empty();
	};

	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
		{
			// z-min and z-max
			if (!hasPortal(4))
				createQuad(vertexAt(i,j,0), vertexAt(i,j+1,0),
				           vertexAt(i+1,j+1,0), vertexAt(i+1,j,0));
			if (!hasPortal(5))
				createQuad(vertexAt(i,j,n), vertexAt(i+1,j,n),
				           vertexAt(i+1,j+1,n), vertexAt(i,j+1,n));

			// y-min and y-max
			if (!hasPortal(2))
				createQuad(vertexAt(i,0,j), vertexAt(i+1,0,j),
				           vertexAt(i+1,0,j+1), vertexAt(i,0,j+1));
			if (!hasPortal(3))
				createQuad(vertexAt(i+1,n,j), vertexAt(i,n,j),
				           vertexAt(i,n,j+1), vertexAt(i+1,n,j+1));

			// x-min and x-max
			if (!hasPortal(0))
				createQuad(vertexAt(0,i+1,j), vertexAt(0,i,j),
				           vertexAt(0,i,j+1), vertexAt(0,i+1,j+1));
			if (!hasPortal(1))
				createQuad(vertexAt(n,i,j), vertexAt(n,i+1,j),
				           vertexAt(n,i+1,j+1), vertexAt(n,i,j+1));
		}

	// A portal-bearing box face keeps the same structured triangular grid as the
	// other box faces.  Only the grid-aligned patch pierced by a terminal is
	// retriangulated with the exact Basolateral rim as an interior constraint.
	// This is important: triangulating the complete face from its perimeter and
	// the small rim alone creates long triangles and leaves all BoxRefs vertices
	// in the face interior unused.
	for (size_t portalIndex = 0; portalIndex < portalPlanes.size(); ++portalIndex)
	{
		const size_t planeIndex = portalPlanes[portalIndex];
		const size_t axis = planeIndex / 2;
		const size_t side = planeIndex % 2;
		const size_t u = (axis + 1) % 3;
		const size_t v = (axis + 2) % 3;
		auto boxVertexOnPlane = [&] (size_t iu, size_t iv) -> Vertex*
		{
			size_t index[3] = {0,0,0};
			index[axis] = side ? n : 0;
			index[u] = iu;
			index[v] = iv;
			return vertexAt(index[0], index[1], index[2]);
		};

		struct PortalPatch
		{
			size_t uMin, uMax, vMin, vMax;
			std::vector<Edge*> rimEdges;
		};
		std::vector<PortalPatch> patches;
		std::set<Edge*> unusedRimEdges(portalEdges[planeIndex].begin(),
		                                    portalEdges[planeIndex].end());
		while (!unusedRimEdges.empty())
		{
			PortalPatch patch;
			std::set<Vertex*> componentVertices;
			std::vector<Edge*> pending;
			pending.push_back(*unusedRimEdges.begin());
			unusedRimEdges.erase(unusedRimEdges.begin());
			while (!pending.empty())
			{
				Edge* edge = pending.back();
				pending.pop_back();
				patch.rimEdges.push_back(edge);
				componentVertices.insert(edge->vertex(0));
				componentVertices.insert(edge->vertex(1));
				for (std::set<Edge*>::iterator it = unusedRimEdges.begin();
				     it != unusedRimEdges.end(); )
				{
					Edge* candidate = *it;
					if (componentVertices.count(candidate->vertex(0))
						|| componentVertices.count(candidate->vertex(1)))
					{
						pending.push_back(candidate);
						unusedRimEdges.erase(it++);
					}
					else ++it;
				}
			}

			number minimumU = std::numeric_limits<number>::max();
			number maximumU = -std::numeric_limits<number>::max();
			number minimumV = std::numeric_limits<number>::max();
			number maximumV = -std::numeric_limits<number>::max();
			for (size_t i = 0; i < patch.rimEdges.size(); ++i)
				for (size_t endpoint = 0; endpoint < 2; ++endpoint)
				{
					const vector3& p = aaPos[patch.rimEdges[i]->vertex(endpoint)];
					minimumU = std::min(minimumU, p[u]);
					maximumU = std::max(maximumU, p[u]);
					minimumV = std::min(minimumV, p[v]);
					maximumV = std::max(maximumV, p[v]);
				}
			const number scaleU = static_cast<number>(n) / (boxMax[u] - boxMin[u]);
			const number scaleV = static_cast<number>(n) / (boxMax[v] - boxMin[v]);
			patch.uMin = std::min(n-1, static_cast<size_t>(std::max((number)0.0,
				std::floor((minimumU - boxMin[u]) * scaleU))));
			patch.uMax = std::max(patch.uMin + 1, std::min(n, static_cast<size_t>(
				std::ceil((maximumU - boxMin[u]) * scaleU))));
			patch.vMin = std::min(n-1, static_cast<size_t>(std::max((number)0.0,
				std::floor((minimumV - boxMin[v]) * scaleV))));
			patch.vMax = std::max(patch.vMin + 1, std::min(n, static_cast<size_t>(
				std::ceil((maximumV - boxMin[v]) * scaleV))));
			patches.push_back(patch);
		}

		// Coalesce overlapping patches so no box cell is generated twice.
		for (size_t i = 0; i < patches.size(); ++i)
			for (size_t j = i + 1; j < patches.size(); )
			{
				const bool overlap = patches[i].uMin < patches[j].uMax
					&& patches[j].uMin < patches[i].uMax
					&& patches[i].vMin < patches[j].vMax
					&& patches[j].vMin < patches[i].vMax;
				if (!overlap) {++j; continue;}
				patches[i].uMin = std::min(patches[i].uMin, patches[j].uMin);
				patches[i].uMax = std::max(patches[i].uMax, patches[j].uMax);
				patches[i].vMin = std::min(patches[i].vMin, patches[j].vMin);
				patches[i].vMax = std::max(patches[i].vMax, patches[j].vMax);
				patches[i].rimEdges.insert(patches[i].rimEdges.end(),
				                           patches[j].rimEdges.begin(), patches[j].rimEdges.end());
				patches.erase(patches.begin() + j);
				j = i + 1;
			}

		// Build a conforming, graded 4 -> 2 -> 1 surface grid around every
		// terminal patch. Coordinates are represented on a four-times finer
		// integer lattice so all transition vertices are shared exactly.
		const size_t adaptiveScale = 4;
		std::map<std::pair<size_t,size_t>,Vertex*> adaptiveVertices;
		auto adaptiveVertex = [&] (size_t fineU, size_t fineV) -> Vertex*
		{
			const std::pair<size_t,size_t> key(fineU,fineV);
			std::map<std::pair<size_t,size_t>,Vertex*>::iterator found =
				adaptiveVertices.find(key);
			if (found != adaptiveVertices.end()) return found->second;
			if (fineU % adaptiveScale == 0 && fineV % adaptiveScale == 0)
			{
				Vertex* coarse = boxVertexOnPlane(fineU / adaptiveScale,
				                                  fineV / adaptiveScale);
				adaptiveVertices[key] = coarse;
				return coarse;
			}
			Vertex* vertex = *g.create<RegularVertex>();
			vector3 point;
			point[axis] = side ? boxMax[axis] : boxMin[axis];
			point[u] = boxMin[u] + (boxMax[u] - boxMin[u])
				* static_cast<number>(fineU)
				/ static_cast<number>(n * adaptiveScale);
			point[v] = boxMin[v] + (boxMax[v] - boxMin[v])
				* static_cast<number>(fineV)
				/ static_cast<number>(n * adaptiveScale);
			aaPos[vertex] = point;
			sh.assign_subset(vertex, subsetIndex);
			adaptiveVertices[key] = vertex;
			return vertex;
		};

		auto cellLevel = [&] (size_t iu, size_t iv) -> size_t
		{
			size_t level = 0;
			for (size_t p = 0; p < patches.size(); ++p)
			{
				const PortalPatch& patch = patches[p];
				const size_t distanceU = iu < patch.uMin ? patch.uMin - iu
					: (iu >= patch.uMax ? iu - patch.uMax + 1 : 0);
				const size_t distanceV = iv < patch.vMin ? patch.vMin - iv
					: (iv >= patch.vMax ? iv - patch.vMax + 1 : 0);
				const size_t distance = std::max(distanceU,distanceV);
				if (distance == 0) return 2; // port cell: four subdivisions
				if (distance == 1) level = std::max(level,(size_t)1); // neighbor: two
			}
			return level;
		};

		auto edgeDivisions = [&] (size_t iu, size_t iv, int du, int dv) -> size_t
		{
			const size_t ownLevel = cellLevel(iu,iv);
			const int neighborU = static_cast<int>(iu) + du;
			const int neighborV = static_cast<int>(iv) + dv;
			// A face-border edge is shared with a perpendicular box face. Keep it
			// at the original segmentation unless that face is refined as well;
			// otherwise the added points would be hanging vertices on the box seam.
			if (neighborU < 0 || neighborV < 0
				|| neighborU >= static_cast<int>(n) || neighborV >= static_cast<int>(n))
				return 1;
			const size_t neighborLevel = cellLevel(static_cast<size_t>(neighborU),
			                                             static_cast<size_t>(neighborV));
			return static_cast<size_t>(1) << std::max(ownLevel,neighborLevel);
		};

		for (size_t iu = 0; iu < n; ++iu)
			for (size_t iv = 0; iv < n; ++iv)
			{
				if (cellLevel(iu,iv) == 2) continue; // filled with the circular port below
				const size_t bottomDivisions = edgeDivisions(iu,iv,0,-1);
				const size_t rightDivisions = edgeDivisions(iu,iv,1,0);
				const size_t topDivisions = edgeDivisions(iu,iv,0,1);
				const size_t leftDivisions = edgeDivisions(iu,iv,-1,0);
				if (bottomDivisions == 1 && rightDivisions == 1
					&& topDivisions == 1 && leftDivisions == 1)
				{
					createQuad(boxVertexOnPlane(iu,iv), boxVertexOnPlane(iu+1,iv),
					           boxVertexOnPlane(iu+1,iv+1), boxVertexOnPlane(iu,iv+1));
					continue;
				}

				std::vector<Vertex*> boundary;
				for (size_t k = 0; k < bottomDivisions; ++k)
					boundary.push_back(adaptiveVertex(adaptiveScale*iu
						+ k*adaptiveScale/bottomDivisions, adaptiveScale*iv));
				for (size_t k = 0; k < rightDivisions; ++k)
					boundary.push_back(adaptiveVertex(adaptiveScale*(iu+1), adaptiveScale*iv
						+ k*adaptiveScale/rightDivisions));
				for (size_t k = 0; k < topDivisions; ++k)
					boundary.push_back(adaptiveVertex(adaptiveScale*(iu+1)
						- k*adaptiveScale/topDivisions, adaptiveScale*(iv+1)));
				for (size_t k = 0; k < leftDivisions; ++k)
					boundary.push_back(adaptiveVertex(adaptiveScale*iu, adaptiveScale*(iv+1)
						- k*adaptiveScale/leftDivisions));
				Vertex* center = adaptiveVertex(adaptiveScale*iu + adaptiveScale/2,
				                                adaptiveScale*iv + adaptiveScale/2);
				for (size_t k = 0; k < boundary.size(); ++k)
				{
					Vertex* a = boundary[k];
					Vertex* b = boundary[(k+1) % boundary.size()];
					createEdge(a,b);
					createEdge(a,center);
					createEdge(b,center);
					Triangle* face = *g.create<Triangle>(TriangleDescriptor(a,b,center));
					sh.assign_subset(face, subsetIndex);
				}
			}

		for (size_t patchIndex = 0; patchIndex < patches.size(); ++patchIndex)
		{
			PortalPatch& patch = patches[patchIndex];
			std::vector<std::pair<Vertex*,Vertex*> > constraints;
			const size_t portSegments = 4;
			for (size_t iu = patch.uMin; iu < patch.uMax; ++iu)
			{
				const size_t bottomSegments = patch.vMin == 0 ? 1 : portSegments;
				const size_t topSegments = patch.vMax == n ? 1 : portSegments;
				for (size_t segment = 0; segment < bottomSegments; ++segment)
				{
					constraints.push_back(std::make_pair(
						adaptiveVertex(adaptiveScale*iu
							+ segment*adaptiveScale/bottomSegments, adaptiveScale*patch.vMin),
						adaptiveVertex(adaptiveScale*iu
							+ (segment+1)*adaptiveScale/bottomSegments, adaptiveScale*patch.vMin)));
				}
				for (size_t segment = 0; segment < topSegments; ++segment)
				{
					constraints.push_back(std::make_pair(
						adaptiveVertex(adaptiveScale*iu
							+ segment*adaptiveScale/topSegments, adaptiveScale*patch.vMax),
						adaptiveVertex(adaptiveScale*iu
							+ (segment+1)*adaptiveScale/topSegments, adaptiveScale*patch.vMax)));
				}
			}
			for (size_t iv = patch.vMin; iv < patch.vMax; ++iv)
			{
				const size_t leftSegments = patch.uMin == 0 ? 1 : portSegments;
				const size_t rightSegments = patch.uMax == n ? 1 : portSegments;
				for (size_t segment = 0; segment < leftSegments; ++segment)
				{
					constraints.push_back(std::make_pair(
						adaptiveVertex(adaptiveScale*patch.uMin, adaptiveScale*iv
							+ segment*adaptiveScale/leftSegments),
						adaptiveVertex(adaptiveScale*patch.uMin, adaptiveScale*iv
							+ (segment+1)*adaptiveScale/leftSegments)));
				}
				for (size_t segment = 0; segment < rightSegments; ++segment)
				{
					constraints.push_back(std::make_pair(
						adaptiveVertex(adaptiveScale*patch.uMax, adaptiveScale*iv
							+ segment*adaptiveScale/rightSegments),
						adaptiveVertex(adaptiveScale*patch.uMax, adaptiveScale*iv
							+ (segment+1)*adaptiveScale/rightSegments)));
				}
			}

			// Guide the annular triangulation explicitly.  Each square-boundary
			// vertex is connected to the Basolateral rim vertex with the closest
			// polar direction about the port center.  The two convex cycles are
			// concentric, so these ordered spokes do not cross and avoid the poor
			// corner-biased connections selected by an unconstrained sweep.
			std::set<Vertex*> boundaryVertices;
			for (size_t i = 0; i < constraints.size(); ++i)
			{
				boundaryVertices.insert(constraints[i].first);
				boundaryVertices.insert(constraints[i].second);
			}
			std::set<Vertex*> rimVertexSet;
			vector3 rimCenter;
			VecSet(rimCenter, 0.0);
			for (size_t i = 0; i < patch.rimEdges.size(); ++i)
				for (size_t endpoint = 0; endpoint < 2; ++endpoint)
					rimVertexSet.insert(patch.rimEdges[i]->vertex(endpoint));
			for (std::set<Vertex*>::iterator it = rimVertexSet.begin();
			     it != rimVertexSet.end(); ++it)
				VecAdd(rimCenter, rimCenter, aaPos[*it]);
			VecScale(rimCenter, rimCenter,
			         1.0 / static_cast<number>(rimVertexSet.size()));
			for (std::set<Vertex*>::iterator outer = boundaryVertices.begin();
			     outer != boundaryVertices.end(); ++outer)
			{
				Vertex* best = NULL;
				number bestAlignment = -std::numeric_limits<number>::max();
				const number outerU = aaPos[*outer][u] - rimCenter[u];
				const number outerV = aaPos[*outer][v] - rimCenter[v];
				const number outerLength = std::sqrt(outerU*outerU + outerV*outerV);
				for (std::set<Vertex*>::iterator inner = rimVertexSet.begin();
				     inner != rimVertexSet.end(); ++inner)
				{
					const number innerU = aaPos[*inner][u] - rimCenter[u];
					const number innerV = aaPos[*inner][v] - rimCenter[v];
					const number innerLength = std::sqrt(innerU*innerU + innerV*innerV);
					if (outerLength <= 0.0 || innerLength <= 0.0) continue;
					const number alignment = (outerU*innerU + outerV*innerV)
						/ (outerLength*innerLength);
					if (alignment > bestAlignment)
					{
						bestAlignment = alignment;
						best = *inner;
					}
				}
				if (best) constraints.push_back(std::make_pair(*outer,best));
			}

			for (size_t i = 0; i < patch.rimEdges.size(); ++i)
				constraints.push_back(std::make_pair(patch.rimEdges[i]->vertex(0),
				                                     patch.rimEdges[i]->vertex(1)));

			std::map<Vertex*,int> vertexIndices;
			std::vector<Vertex*> indexedVertices;
			std::vector<vector3> positions;
			std::vector<int> constrainedIndices;
			auto indexOf = [&] (Vertex* vertex) -> int
			{
				std::map<Vertex*,int>::iterator found = vertexIndices.find(vertex);
				if (found != vertexIndices.end()) return found->second;
				const int index = (int)indexedVertices.size();
				vertexIndices[vertex] = index;
				indexedVertices.push_back(vertex);
				positions.push_back(aaPos[vertex]);
				return index;
			};
			for (size_t i = 0; i < constraints.size(); ++i)
			{
				constrainedIndices.push_back(indexOf(constraints[i].first));
				constrainedIndices.push_back(indexOf(constraints[i].second));
			}
			std::vector<int> triangles;
			UG_COND_THROW(!TriangleFill_SweepLine(triangles, positions, constrainedIndices),
			              "Could not triangulate a local OuterWall terminal patch.");

			auto insideTerminalRim = [&] (const vector3& point) -> bool
			{
				bool inside = false;
				for (size_t i = 0; i < patch.rimEdges.size(); ++i)
				{
					const vector3& a = aaPos[patch.rimEdges[i]->vertex(0)];
					const vector3& b = aaPos[patch.rimEdges[i]->vertex(1)];
					const bool crosses = ((a[v] > point[v]) != (b[v] > point[v]));
					if (!crosses) continue;
					const number crossingU = a[u] + (point[v] - a[v])
						* (b[u] - a[u]) / (b[v] - a[v]);
					if (point[u] < crossingU) inside = !inside;
				}
				return inside;
			};
			for (size_t i = 0; i < triangles.size(); i += 3)
			{
				Vertex* a = indexedVertices[triangles[i]];
				Vertex* b = indexedVertices[triangles[i+1]];
				Vertex* c = indexedVertices[triangles[i+2]];
				vector3 center;
				VecAdd(center, aaPos[a], aaPos[b]);
				VecAdd(center, center, aaPos[c]);
				VecScale(center, center, 1.0 / 3.0);
				if (insideTerminalRim(center)) continue;
				Triangle* face = *g.create<Triangle>(TriangleDescriptor(a,b,c));
				sh.assign_subset(face, subsetIndex);
			}
			for (size_t i = 0; i < constraints.size(); ++i)
			{
				Vertex* a = constraints[i].first;
				Vertex* b = constraints[i].second;
				Edge* edge = g.get_edge(a,b);
				if (!edge) edge = *g.create<RegularEdge>(EdgeDescriptor(a,b));
				if (std::find(patch.rimEdges.begin(), patch.rimEdges.end(), edge)
					== patch.rimEdges.end())
					sh.assign_subset(edge, subsetIndex);
			}
		}
		UG_LOGN("OuterWall plane " << planeIndex << ": retained structured grid and "
		        << "retriangulated " << patches.size()
		        << " local terminal patch(es) with an aligned constrained 4-to-2-to-1 transition.");
	}

}


static void tetrahedralize_lumen_and_inter
(
	Grid& g,
	SubsetHandler& sh
)
{
	// Retain the closed neurite surface, but remove its coarse hexahedra so
	// TetGen can create one non-overlapping conforming tetrahedralization.
	std::vector<Face*> interiorFaces;
	for (FaceIterator fit = g.begin<Face>(); fit != g.end<Face>(); ++fit)
	{
		const int numAssocVols = NumAssociatedVolumes(g, *fit);
		if (numAssocVols > 1)
			interiorFaces.push_back(*fit);
		else if (numAssocVols == 1 && sh.get_subset_index(*fit) == 0)
			sh.assign_subset(*fit, 2); // dedicated nonzero Apical TetGen marker
	}
	g.erase(g.begin<Volume>(), g.end<Volume>());
	g.erase(interiorFaces.begin(), interiorFaces.end());
	UG_COND_THROW(!Tetrahedralize(g, sh, 10.0, true, true, aPosition, 1),
	              "Failed to tetrahedralize the Lumen/Inter box geometry.");

	// TetGen assigns ordinary internal tetrahedron faces marker 0. Only the
	// explicitly marked OuterWall (1) and Apical (2) faces are separators.
	Selector boundarySelector(g);
	for (FaceIterator fit = sh.begin<Face>(1); fit != sh.end<Face>(1); ++fit)
		boundarySelector.select(*fit);
	for (FaceIterator fit = sh.begin<Face>(2); fit != sh.end<Face>(2); ++fit)
		boundarySelector.select(*fit);
	SeparateSubsetsByLowerDimSelection<Volume>(g, sh, boundarySelector, true);

	int interSubset = -1;
	std::vector<Volume*> adjacentVolumes;
	for (FaceIterator fit = sh.begin<Face>(1); fit != sh.end<Face>(1); ++fit)
	{
		CollectAssociated(adjacentVolumes, g, *fit);
		if (!adjacentVolumes.empty())
		{
			interSubset = sh.get_subset_index(adjacentVolumes[0]);
			break;
		}
	}
	UG_COND_THROW(interSubset < 0,
	              "Could not identify the Inter volume adjacent to OuterWall.");

	for (int si = 3; si < sh.num_subsets(); ++si)
	{
		if (sh.num<Volume>(si) == 0)
			continue;
		if (si == interSubset)
			sh.set_subset_name("Inter", si);
		else
			sh.set_subset_name("Lumen", si);
	}

	// First close the two volume subsets down to faces, edges, and vertices.
	// Then restore boundary precedence and move Apical back to subset 0.
	std::vector<Face*> outerFaces(sh.begin<Face>(1), sh.end<Face>(1));
	std::vector<Face*> apicalFaces(sh.begin<Face>(2), sh.end<Face>(2));
	CopySubsetIndicesToSides(sh, false);

	auto assignBoundaryClosure = [&] (const std::vector<Face*>& faces, int si)
	{
		std::vector<Edge*> edges;
		for (size_t i = 0; i < faces.size(); ++i)
		{
			Face* face = faces[i];
			sh.assign_subset(face, si);
			for (size_t j = 0; j < face->num_vertices(); ++j)
				sh.assign_subset(face->vertex(j), si);
			CollectAssociated(edges, g, face);
			for (size_t j = 0; j < edges.size(); ++j)
				sh.assign_subset(edges[j], si);
		}
	};
	assignBoundaryClosure(outerFaces, 1);
	assignBoundaryClosure(apicalFaces, 0);
	EraseEmptySubsets(sh); // removes the temporary Apical marker subset 2

	UG_LOGN("Created conforming tetrahedral regions; Inter subset=" << interSubset);
}


static void tetrahedralize_lumen_membrane_and_inter
(
	Grid& g,
	SubsetHandler& sh
)
{
	// create_neurite_with_er supplies the three closed interfaces as follows:
	// subset 2: outer tube (Basolateral), subset 3: inner tube (Apical),
	// subset 4: box (OuterWall). Remove the old hexahedral regions and their
	// cross-section faces, then let TetGen fill all three nested regions at once.
	std::vector<Face*> interiorFaces;
	for (FaceIterator fit = g.begin<Face>(); fit != g.end<Face>(); ++fit)
	{
		// The Apical interface (subset 3) is intentionally shared by the old
		// Lumen and Membrane volumes and must survive as a TetGen constraint.
		// Other two-volume faces are only coarse axial cross-sections.
		if (NumAssociatedVolumes(g, *fit) > 1 && sh.get_subset_index(*fit) != 3)
			interiorFaces.push_back(*fit);
	}
	g.erase(g.begin<Volume>(), g.end<Volume>());
	g.erase(interiorFaces.begin(), interiorFaces.end());

	UG_COND_THROW(!Tetrahedralize(g, sh, 10.0, true, true, aPosition, 1),
	              "Failed to tetrahedralize the Lumen/Membrane/Inter geometry.");

	// Marker-zero faces are ordinary faces created by TetGen. Only the three
	// explicit closed interfaces are allowed to separate volume components.
	// TetGen/UG4 may represent a constrained triangle and a tetrahedron side as
	// distinct, coincident Face objects. Select separators by their vertex keys
	// so the flood fill cannot cross an unmarked duplicate of a marked surface.
	typedef std::array<std::uintptr_t, 3> TriangleKey;
	auto vertexTriangleKey = [] (Vertex* a, Vertex* b, Vertex* c) -> TriangleKey
	{
		TriangleKey key = {{
			reinterpret_cast<std::uintptr_t>(a),
			reinterpret_cast<std::uintptr_t>(b),
			reinterpret_cast<std::uintptr_t>(c)
		}};
		std::sort(key.begin(), key.end());
		return key;
	};
	auto triangleKey = [] (Face* face) -> TriangleKey
	{
		UG_COND_THROW(face->num_vertices() != 3,
		              "Expected triangular faces after tetrahedralization.");
		TriangleKey key = {{
			reinterpret_cast<std::uintptr_t>(face->vertex(0)),
			reinterpret_cast<std::uintptr_t>(face->vertex(1)),
			reinterpret_cast<std::uintptr_t>(face->vertex(2))
		}};
		std::sort(key.begin(), key.end());
		return key;
	};

	std::set<TriangleKey> separatorKeys;
	for (int si = 2; si <= 4; ++si)
		for (FaceIterator fit = sh.begin<Face>(si); fit != sh.end<Face>(si); ++fit)
			separatorKeys.insert(triangleKey(*fit));

	// Build tetrahedron adjacency directly from vertex triples. This avoids the
	// duplicate Face-object ambiguity in Grid::get_side for TetGen constraints.
	std::map<TriangleKey, std::vector<Volume*> > faceVolumes;
	for (VolumeIterator vit = g.begin<Volume>(); vit != g.end<Volume>(); ++vit)
	{
		Volume* vol = *vit;
		UG_COND_THROW(vol->num_vertices() != 4,
		              "Expected tetrahedra after tetrahedralization.");
		for (size_t omit = 0; omit < 4; ++omit)
		{
			Vertex* fv[3];
			size_t k = 0;
			for (size_t j = 0; j < 4; ++j)
				if (j != omit) fv[k++] = vol->vertex(j);
			faceVolumes[vertexTriangleKey(fv[0], fv[1], fv[2])].push_back(vol);
		}
		sh.assign_subset(vol, -1);
	}

	int nextRegionSubset = sh.num_subsets();
	std::vector<Volume*> floodStack;
	for (VolumeIterator vit = g.begin<Volume>(); vit != g.end<Volume>(); ++vit)
	{
		Volume* seed = *vit;
		if (sh.get_subset_index(seed) != -1)
			continue;
		floodStack.push_back(seed);
		while (!floodStack.empty())
		{
			Volume* vol = floodStack.back();
			floodStack.pop_back();
			if (sh.get_subset_index(vol) != -1)
				continue;
			sh.assign_subset(vol, nextRegionSubset);
			for (size_t omit = 0; omit < 4; ++omit)
			{
				Vertex* fv[3];
				size_t k = 0;
				for (size_t j = 0; j < 4; ++j)
					if (j != omit) fv[k++] = vol->vertex(j);
				const TriangleKey key = vertexTriangleKey(fv[0], fv[1], fv[2]);
				if (separatorKeys.count(key) != 0)
					continue;
				const std::vector<Volume*>& neighbors = faceVolumes[key];
				for (size_t j = 0; j < neighbors.size(); ++j)
					if (sh.get_subset_index(neighbors[j]) == -1)
						floodStack.push_back(neighbors[j]);
			}
		}
		++nextRegionSubset;
	}

	for (int marker = 2; marker <= 4; ++marker)
	{
		std::set<int> adjacentRegionSubsets;
		for (FaceIterator fit = sh.begin<Face>(marker);
		     fit != sh.end<Face>(marker); ++fit)
		{
			const std::vector<Volume*>& markerVols = faceVolumes[triangleKey(*fit)];
			for (size_t i = 0; i < markerVols.size(); ++i)
				adjacentRegionSubsets.insert(sh.get_subset_index(markerVols[i]));
		}
		UG_LOG("Nested-region marker " << marker << ": faces="
		       << sh.num<Face>(marker) << ", adjacent volume subsets=");
		for (std::set<int>::const_iterator it = adjacentRegionSubsets.begin();
		     it != adjacentRegionSubsets.end(); ++it)
			UG_LOG(" " << *it);
		UG_LOGN("");
	}

	auto adjacentSubsetNot = [&] (int faceSubset, int excludedSubset) -> int
	{
		for (FaceIterator fit = sh.begin<Face>(faceSubset);
		     fit != sh.end<Face>(faceSubset); ++fit)
		{
			const std::vector<Volume*>& vols = faceVolumes[triangleKey(*fit)];
			for (size_t i = 0; i < vols.size(); ++i)
			{
				const int si = sh.get_subset_index(vols[i]);
				if (si != excludedSubset)
					return si;
			}
		}
		return -1;
	};

	const int interSubset = adjacentSubsetNot(4, -1);
	UG_COND_THROW(interSubset < 0,
	              "Could not identify Inter adjacent to OuterWall.");
	const int membraneSubset = adjacentSubsetNot(2, interSubset);
	UG_COND_THROW(membraneSubset < 0 || membraneSubset == interSubset,
	              "Could not identify Membrane inside Basolateral.");
	const int lumenSubset = adjacentSubsetNot(3, membraneSubset);
	UG_COND_THROW(lumenSubset < 0 || lumenSubset == membraneSubset ||
	              lumenSubset == interSubset,
	              "Could not identify Lumen inside Apical.");

	sh.set_subset_name("Inter", interSubset);
	sh.set_subset_name("Membrane", membraneSubset);
	sh.set_subset_name("Lumen", lumenSubset);

	std::vector<Face*> basolateralFaces(sh.begin<Face>(2), sh.end<Face>(2));
	std::vector<Face*> apicalFaces(sh.begin<Face>(3), sh.end<Face>(3));
	std::vector<Face*> outerFaces(sh.begin<Face>(4), sh.end<Face>(4));

	// Give each volume region a complete lower-dimensional closure first, then
	// restore the named physical boundaries with boundary precedence.
	CopySubsetIndicesToSides(sh, false);
	auto assignBoundaryClosure = [&] (const std::vector<Face*>& faces, int si)
	{
		std::vector<Edge*> edges;
		for (size_t i = 0; i < faces.size(); ++i)
		{
			Face* face = faces[i];
			sh.assign_subset(face, si);
			for (size_t j = 0; j < face->num_vertices(); ++j)
				sh.assign_subset(face->vertex(j), si);
			CollectAssociated(edges, g, face);
			for (size_t j = 0; j < edges.size(); ++j)
				sh.assign_subset(edges[j], si);
		}
	};
	assignBoundaryClosure(basolateralFaces, 2);
	assignBoundaryClosure(apicalFaces, 3);
	assignBoundaryClosure(outerFaces, 4);
	sh.set_subset_name("Basolateral", 2);
	sh.set_subset_name("Apical", 3);
	sh.set_subset_name("OuterWall", 4);
	EraseEmptySubsets(sh);

	UG_LOGN("Created conforming nested regions: Lumen=" << lumenSubset
	        << ", Membrane=" << membraneSubset << ", Inter=" << interSubset);
}


static void add_inter_to_version2_nephron
(
	Grid& g,
	SubsetHandler& sh,
	SubsetHandler& psh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	number tetQuality
)
{
	UG_COND_THROW(tetQuality <= 0.0,
	              "Inter tetrahedral quality must be positive.");

	typedef std::array<std::uintptr_t, 3> TriangleKey;
	auto triKey = [] (Vertex* a, Vertex* b, Vertex* c) -> TriangleKey
	{
		TriangleKey key = {{reinterpret_cast<std::uintptr_t>(a),
		                    reinterpret_cast<std::uintptr_t>(b),
		                    reinterpret_cast<std::uintptr_t>(c)}};
		std::sort(key.begin(), key.end());
		return key;
	};

	struct TubeVolumeRecord
	{
		Volume* original;
		std::vector<Vertex*> vrts;
		std::vector<std::vector<Vertex*> > faces;
		std::vector<int> faceSubsets;
		vector3 boxMin;
		vector3 boxMax;
		int subset;
	};
	std::vector<TubeVolumeRecord> oldTubeVolumes;
	for (VolumeIterator vit = g.begin<Volume>(); vit != g.end<Volume>(); ++vit)
	{
		Volume* vol = *vit;
		const int subset = sh.get_subset_index(vol);
		const bool isHex = vol->num_vertices() == 8 && vol->num_faces() == 6;
		const bool isPrism = vol->num_vertices() == 6 && vol->num_faces() == 5;
		UG_COND_THROW(!isHex && !(subset == 1 && isPrism),
		              "The conforming nephron path supports membrane hexahedra "
		              "and lumen hexahedra/prisms only.");
		TubeVolumeRecord rec;
		rec.original = vol;
		rec.subset = subset;
		rec.vrts.resize(vol->num_vertices());
		for (size_t i = 0; i < vol->num_vertices(); ++i)
		{
			rec.vrts[i] = vol->vertex(i);
			if (i == 0) rec.boxMin = rec.boxMax = aaPos[rec.vrts[i]];
			else for (size_t d = 0; d < 3; ++d)
			{
				rec.boxMin[d] = std::min(rec.boxMin[d], aaPos[rec.vrts[i]][d]);
				rec.boxMax[d] = std::max(rec.boxMax[d], aaPos[rec.vrts[i]][d]);
			}
		}
		rec.faces.resize(vol->num_faces());
		rec.faceSubsets.resize(vol->num_faces());
		for (size_t i = 0; i < vol->num_faces(); ++i)
		{
			Face* face = g.get_face(vol, i);
			UG_COND_THROW(!face ||
			              (face->num_vertices() != 3 && face->num_vertices() != 4),
			              "Could not record a triangular/quadrilateral tube face.");
			rec.faces[i].resize(face->num_vertices());
			for (size_t j = 0; j < face->num_vertices(); ++j)
				rec.faces[i][j] = face->vertex(j);
			rec.faceSubsets[i] = sh.get_subset_index(face);
		}
		oldTubeVolumes.push_back(rec);
	}

	// Keep only the two version-2 tube surfaces and the box surface. TetGen then
	// fills precisely the exterior region between Basolateral and OuterWall.
	g.erase(g.begin<Volume>(), g.end<Volume>());
	std::vector<Face*> nonBoundaryFaces;
	for (FaceIterator fit = g.begin<Face>(); fit != g.end<Face>(); ++fit)
	{
		const int si = sh.get_subset_index(*fit);
		if (si != 2 && si != 3 && si != 4)
			nonBoundaryFaces.push_back(*fit);
	}
	g.erase(nonBoundaryFaces.begin(), nonBoundaryFaces.end());
	UG_COND_THROW(!Tetrahedralize(g, sh, tetQuality, true, true, aPosition, 1),
	              "Failed to tetrahedralize the exterior Inter region.");

	// TetGen fills the complete box. Remove every newly created tetrahedron whose
	// center lies in one of the original version-2 tube volumes. A small spatial bin
	// map keeps this containment test local along the long nephron trace.
	number binSize = 0.0;
	for (size_t h = 0; h < oldTubeVolumes.size(); ++h)
		for (size_t d = 0; d < 3; ++d)
			binSize = std::max(binSize,
			                   oldTubeVolumes[h].boxMax[d] - oldTubeVolumes[h].boxMin[d]);
	UG_COND_THROW(binSize <= 0.0,
	              "Cannot build a spatial index for zero-size tube volumes.");
	binSize *= 1.01;
	typedef std::array<int, 3> BinKey;
	std::map<BinKey, std::vector<size_t> > tubeBins;
	auto binKey = [&] (const vector3& p) -> BinKey
	{
		return BinKey{{(int)floor(p[0]/binSize),
		               (int)floor(p[1]/binSize),
		               (int)floor(p[2]/binSize)}};
	};
	for (size_t h = 0; h < oldTubeVolumes.size(); ++h)
	{
		const BinKey lo = binKey(oldTubeVolumes[h].boxMin);
		const BinKey hi = binKey(oldTubeVolumes[h].boxMax);
		for (int i = lo[0]; i <= hi[0]; ++i)
			for (int j = lo[1]; j <= hi[1]; ++j)
				for (int k = lo[2]; k <= hi[2]; ++k)
					tubeBins[BinKey{{i,j,k}}].push_back(h);
	}

	auto pointInRecordedVolume = [&] (const TubeVolumeRecord& rec,
	                                  const vector3& p) -> bool
	{
		number recordScale = 0.0;
		for (size_t d = 0; d < 3; ++d)
			recordScale = std::max(recordScale, rec.boxMax[d] - rec.boxMin[d]);
		const number containmentTolerance = std::max(
			std::numeric_limits<number>::epsilon() * recordScale * 64.0,
			recordScale * 1e-6);
		vector3 volumeCenter(0.0);
		for (size_t i = 0; i < rec.vrts.size(); ++i)
			volumeCenter += aaPos[rec.vrts[i]];
		volumeCenter *= 1.0 / rec.vrts.size();
		for (size_t f = 0; f < rec.faces.size(); ++f)
		{
			const vector3& a = aaPos[rec.faces[f][0]];
			const vector3& b = aaPos[rec.faces[f][1]];
			const vector3& c = aaPos[rec.faces[f][2]];
			vector3 ab, ac, normal, faceCenter(0.0), toCenter, toPoint;
			VecSubtract(ab, b, a);
			VecSubtract(ac, c, a);
			VecCross(normal, ab, ac);
			for (size_t j = 0; j < rec.faces[f].size(); ++j)
				faceCenter += aaPos[rec.faces[f][j]];
			faceCenter *= 1.0 / rec.faces[f].size();
			VecSubtract(toCenter, volumeCenter, faceCenter);
			if (VecProd(normal, toCenter) > 0.0) normal *= -1.0;
			VecSubtract(toPoint, p, faceCenter);
			if (VecProd(normal, toPoint) > containmentTolerance
				* sqrt(VecNormSquared(normal)))
				return false;
		}
		return true;
	};
	std::vector<Volume*> overlappingTetgenVolumes;
	for (VolumeIterator vit = g.begin<Volume>(); vit != g.end<Volume>(); ++vit)
	{
		Volume* vol = *vit;
		vector3 center(0.0);
		for (size_t i = 0; i < vol->num_vertices(); ++i) center += aaPos[vol->vertex(i)];
		center *= 1.0 / vol->num_vertices();
		const std::vector<size_t>& candidates = tubeBins[binKey(center)];
		for (size_t i = 0; i < candidates.size(); ++i)
		{
			const TubeVolumeRecord& rec = oldTubeVolumes[candidates[i]];
			bool inBox = true;
			for (size_t d = 0; d < 3; ++d)
				inBox = inBox && center[d] >= rec.boxMin[d] && center[d] <= rec.boxMax[d];
			if (inBox && pointInRecordedVolume(rec, center))
			{
				overlappingTetgenVolumes.push_back(vol);
				break;
			}
		}
	}
	g.erase(overlappingTetgenVolumes.begin(), overlappingTetgenVolumes.end());
	// TetGen's bulk removal can leave both the vertex-to-face lookup and the
	// per-volume face cache with erased face pointers.  Volume auto-generation
	// uses the former while the closure pass below uses the latter.  Rebuild both
	// attachments from the live grid before creating any new volume.
	const bool storedVertexFaces =
		g.option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES);
	const bool storedVolumeFaces =
		g.option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES);
	if (storedVolumeFaces)
		g.disable_options(VOLOPT_STORE_ASSOCIATED_FACES);
	if (storedVertexFaces)
		g.disable_options(VRTOPT_STORE_ASSOCIATED_FACES);
	if (storedVertexFaces)
		g.enable_options(VRTOPT_STORE_ASSOCIATED_FACES);
	if (storedVolumeFaces)
		g.enable_options(VOLOPT_STORE_ASSOCIATED_FACES);
	for (VolumeIterator vit = g.begin<Volume>(); vit != g.end<Volume>(); ++vit)
		sh.assign_subset(*vit, 5);
	UG_LOGN("Removed " << overlappingTetgenVolumes.size()
	        << " overlapping TetGen volumes from inside Lumen/Membrane.");

	std::set<TriangleKey> tetgenTriangles;
	for (FaceIterator fit = g.begin<Face>(); fit != g.end<Face>(); ++fit)
	{
		Face* face = *fit;
		if (face->num_vertices() == 3)
			tetgenTriangles.insert(triKey(face->vertex(0), face->vertex(1), face->vertex(2)));
	}

	std::vector<Volume*> restoredVolumes;
	for (size_t h = 0; h < oldTubeVolumes.size(); ++h)
	{
		const TubeVolumeRecord& rec = oldTubeVolumes[h];
		if (rec.subset == 1)
		{
			Volume* lumenVolume;
			if (rec.vrts.size() == 8)
				lumenVolume = *g.create<Hexahedron>(HexahedronDescriptor(
					rec.vrts[0], rec.vrts[1], rec.vrts[2], rec.vrts[3],
					rec.vrts[4], rec.vrts[5], rec.vrts[6], rec.vrts[7]));
			else
				lumenVolume = *g.create<Prism>(PrismDescriptor(
					rec.vrts[0], rec.vrts[1], rec.vrts[2],
					rec.vrts[3], rec.vrts[4], rec.vrts[5]));
			sh.assign_subset(lumenVolume, 1);
			restoredVolumes.push_back(lumenVolume);
			continue;
		}
		UG_COND_THROW(rec.vrts.size() != 8 || rec.faces.size() != 6,
		              "Membrane reconstruction requires hexahedral cells.");

		RegularVertex* center = *g.create<RegularVertex>();
		aaPos[center] = vector3(0.0);
		for (size_t i = 0; i < 8; ++i) aaPos[center] += aaPos[rec.vrts[i]];
		aaPos[center] *= 0.125;
		sh.assign_subset(center, rec.subset);

		for (size_t f = 0; f < 6; ++f)
		{
			Vertex* a = rec.faces[f][0];
			Vertex* b = rec.faces[f][1];
			Vertex* c = rec.faces[f][2];
			Vertex* d = rec.faces[f][3];
			const bool tetgenDiagAC =
				tetgenTriangles.count(triKey(a,b,c)) &&
				tetgenTriangles.count(triKey(a,c,d));
			const bool tetgenDiagBD =
				tetgenTriangles.count(triKey(a,b,d)) &&
				tetgenTriangles.count(triKey(b,c,d));

			bool useAC;
			if (tetgenDiagAC || tetgenDiagBD)
				useAC = tetgenDiagAC;
			else
			{
				const std::array<std::uintptr_t,2> ac = {{
					reinterpret_cast<std::uintptr_t>(std::min(a,c)),
					reinterpret_cast<std::uintptr_t>(std::max(a,c))}};
				const std::array<std::uintptr_t,2> bd = {{
					reinterpret_cast<std::uintptr_t>(std::min(b,d)),
					reinterpret_cast<std::uintptr_t>(std::max(b,d))}};
				useAC = ac < bd;
			}

			// Only Basolateral must match TetGen's triangular interface. All
			// other membrane sides stay quadrilateral pyramid bases, preserving
			// conformity with the version-2 Lumen and neighboring membrane cells.
			if (rec.faceSubsets[f] != 2)
			{
				Pyramid* pyramid = *g.create<Pyramid>(PyramidDescriptor(a,b,c,d,center));
				sh.assign_subset(pyramid, rec.subset);
				restoredVolumes.push_back(pyramid);
				continue;
			}

			Tetrahedron* t0;
			Tetrahedron* t1;
			if (useAC)
			{
				t0 = *g.create<Tetrahedron>(TetrahedronDescriptor(center,a,b,c));
				t1 = *g.create<Tetrahedron>(TetrahedronDescriptor(center,a,c,d));
			}
			else
			{
				t0 = *g.create<Tetrahedron>(TetrahedronDescriptor(center,a,b,d));
				t1 = *g.create<Tetrahedron>(TetrahedronDescriptor(center,b,c,d));
			}
			sh.assign_subset(t0, rec.subset);
			sh.assign_subset(t1, rec.subset);
			restoredVolumes.push_back(t0);
			restoredVolumes.push_back(t1);
		}
	}
	FixOrientation(g, restoredVolumes.begin(), restoredVolumes.end(), aaPos);

	std::vector<Face*> basolateralFaces(sh.begin<Face>(2), sh.end<Face>(2));
	std::vector<Face*> apicalFaces(sh.begin<Face>(3), sh.end<Face>(3));
	std::vector<Face*> outerFaces(sh.begin<Face>(4), sh.end<Face>(4));
	// Do not use CopySubsetIndicesToSides here. TetGen has just removed and
	// recreated many volumes/sides, and its attached-element lists may still
	// contain stale entries until the grid's next full attachment rebuild. That
	// made O-grids with more than four circumferential vertices dereference an
	// erased side. Close each live volume explicitly through Grid associations.
	std::vector<Face*> closureFaces;
	std::vector<Edge*> closureEdges;
	for (VolumeIterator vit = g.begin<Volume>(); vit != g.end<Volume>(); ++vit)
	{
		Volume* vol = *vit;
		const int si = sh.get_subset_index(vol);
		CollectAssociated(closureFaces, g, vol);
		for (size_t f = 0; f < closureFaces.size(); ++f)
		{
			Face* face = closureFaces[f];
			sh.assign_subset(face, si);
			for (size_t j = 0; j < face->num_vertices(); ++j)
				sh.assign_subset(face->vertex(j), si);
			CollectAssociated(closureEdges, g, face);
			for (size_t e = 0; e < closureEdges.size(); ++e)
				sh.assign_subset(closureEdges[e], si);
		}
	}
	auto restoreBoundary = [&] (const std::vector<Face*>& faces, int si)
	{
		std::vector<Edge*> edges;
		for (size_t i = 0; i < faces.size(); ++i)
		{
			Face* face = faces[i];
			sh.assign_subset(face, si);
			psh.assign_subset(face, si == 4 ? 1 : 0);
			for (size_t j = 0; j < face->num_vertices(); ++j)
			{
				sh.assign_subset(face->vertex(j), si);
				psh.assign_subset(face->vertex(j), si == 4 ? 1 : 0);
			}
			CollectAssociated(edges, g, face);
			for (size_t j = 0; j < edges.size(); ++j)
			{
				sh.assign_subset(edges[j], si);
				psh.assign_subset(edges[j], si == 4 ? 1 : 0);
			}
		}
	};
	restoreBoundary(basolateralFaces, 2);
	restoreBoundary(apicalFaces, 3);
	restoreBoundary(outerFaces, 4);
	sh.set_subset_name("Membrane", 0);
	sh.set_subset_name("Lumen", 1);
	sh.set_subset_name("Basolateral", 2);
	sh.set_subset_name("Apical", 3);
	sh.set_subset_name("OuterWall", 4);
	sh.set_subset_name("Inter", 5);
	UG_LOGN("Added conforming exterior Inter to the version-2 nephron geometry.");
}


static void add_straight_prismatic_inter_to_version2_nephron
(
	Grid& g,
	SubsetHandler& sh,
	SubsetHandler& psh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >&
		aaSurfParams,
	const OuterWallTerminalMetadata& terminalMetadata
)
{
	UG_COND_THROW(!terminalMetadata.hasFixedBox ||
	              !terminalMetadata.terminalsOnOuterWall,
	              "Prismatic Inter requires fixed OuterWall bounds and terminal "
	              "ports on the two end planes.");

	// The straight prototype sweeps one fixed 2-D triangulation through every
	// Basolateral axial station. Group the existing projector-backed ring
	// vertices by their exact axial parameter, then sort each ring by angle.
	std::map<float, std::set<Vertex*> > ringVertexSets;
	std::vector<Face*> basolateralQuads;
	for (FaceIterator fit = sh.begin<Face>(2); fit != sh.end<Face>(2); ++fit)
	{
		Face* face = *fit;
		if (face->num_vertices() != 4) continue;
		basolateralQuads.push_back(face);
		for (size_t i = 0; i < 4; ++i)
		{
			Vertex* vertex = face->vertex(i);
			ringVertexSets[aaSurfParams[vertex].axial].insert(vertex);
		}
	}
	UG_COND_THROW(ringVertexSets.size() < 2 || basolateralQuads.empty(),
	              "Prismatic Inter found no swept Basolateral quadrilateral mesh.");

	std::vector<std::vector<Vertex*> > rings;
	std::vector<float> ringAxial;
	for (std::map<float,std::set<Vertex*> >::const_iterator it =
	     ringVertexSets.begin(); it != ringVertexSets.end(); ++it)
	{
		std::vector<Vertex*> ring(it->second.begin(), it->second.end());
		vector3 localCenter(0.0);
		for (size_t i = 0; i < ring.size(); ++i) localCenter += aaPos[ring[i]];
		localCenter *= 1.0 / ring.size();
		std::sort(ring.begin(), ring.end(), [&] (Vertex* a, Vertex* b)
		{
			number angleA = std::atan2(aaPos[a][1] - localCenter[1],
			                            aaPos[a][0] - localCenter[0]);
			number angleB = std::atan2(aaPos[b][1] - localCenter[1],
			                            aaPos[b][0] - localCenter[0]);
			if (angleA < 0.0) angleA += 2.0 * PI;
			if (angleB < 0.0) angleB += 2.0 * PI;
			return angleA < angleB;
		});
		rings.push_back(ring);
		ringAxial.push_back(it->first);
	}
	const size_t angularCount = rings[0].size();
	UG_COND_THROW(angularCount < 4,
	              "Prismatic Inter requires at least four angular vertices.");
	std::vector<vector3> ringCenters(rings.size(), vector3(0.0));
	for (size_t s = 0; s < rings.size(); ++s)
	{
		for (size_t i = 0; i < rings[s].size(); ++i)
			ringCenters[s] += aaPos[rings[s][i]];
		ringCenters[s] *= 1.0 / rings[s].size();
	}
	for (size_t s = 0; s < rings.size(); ++s)
	{
		UG_COND_THROW(rings[s].size() != angularCount,
		              "Prismatic Inter requires the same angular ring at every "
		              "axial station.");
		for (size_t i = 0; i < angularCount; ++i)
		{
			number angle = std::atan2(
				aaPos[rings[s][i]][1] - ringCenters[s][1],
				aaPos[rings[s][i]][0] - ringCenters[s][0]);
			number referenceAngle = std::atan2(
				aaPos[rings[0][i]][1] - ringCenters[0][1],
				aaPos[rings[0][i]][0] - ringCenters[0][0]);
			const number angularDifference = std::fabs(std::atan2(
				std::sin(angle - referenceAngle),
				std::cos(angle - referenceAngle)));
			UG_COND_THROW(angularDifference > 1e-5,
			              "Prismatic Inter is restricted to a non-twisting straight "
			              "tube in this prototype.");
		}
	}
	UG_COND_THROW(basolateralQuads.size() !=
	              (rings.size() - 1) * angularCount,
	              "Basolateral quadrilateral count does not match the swept rings.");

	const number boxSpan = std::max(
		terminalMetadata.boxMax[0] - terminalMetadata.boxMin[0],
		terminalMetadata.boxMax[1] - terminalMetadata.boxMin[1]);
	const number straightTolerance = std::max<number>(1e-12, boxSpan * 1e-8);
	for (size_t s = 1; s < ringCenters.size(); ++s)
	{
		UG_COND_THROW(std::fabs(ringCenters[s][0] - ringCenters[0][0]) >
		              straightTolerance ||
		              std::fabs(ringCenters[s][1] - ringCenters[0][1]) >
		              straightTolerance,
		              "Prismatic Inter prototype supports only a z-aligned straight "
		              "tube with a fixed cross-section center.");
	}
	UG_COND_THROW(std::fabs(ringCenters.front()[2] -
	              terminalMetadata.boxMin[2]) > straightTolerance ||
	              std::fabs(ringCenters.back()[2] -
	              terminalMetadata.boxMax[2]) > straightTolerance,
	              "Prismatic Inter axial stations do not span the two OuterWall "
	              "terminal planes.");

	// Remove the old box surface. The prism sweep below recreates all six box
	// planes as ordinary one-sided Inter faces, including the two port planes.
	std::vector<Face*> oldOuterFaces(sh.begin<Face>(4), sh.end<Face>(4));
	g.erase(oldOuterFaces.begin(), oldOuterFaces.end());

	// Seven graded radial layers retain a near-membrane scale comparable to the
	// Basolateral edge length while reaching the rectangular far boundary. The
	// final four corner triangles make the outer polygon exactly the original box,
	// rather than a chordal polygon inscribed in it.
	const number radialFractions[] =
		{0.04, 0.10, 0.20, 0.35, 0.55, 0.78, 1.0};
	const size_t addedLayerCount =
		sizeof(radialFractions) / sizeof(radialFractions[0]);
	const size_t layerCount = addedLayerCount + 1;
	const size_t outerLayer = layerCount - 1;
	const number twoPi = 2.0 * PI;

	std::vector<number> angles(angularCount);
	for (size_t i = 0; i < angularCount; ++i)
	{
		angles[i] = std::atan2(
			aaPos[rings[0][i]][1] - ringCenters[0][1],
			aaPos[rings[0][i]][0] - ringCenters[0][0]);
		while (angles[i] < 0.0) angles[i] += twoPi;
		while (angles[i] >= twoPi) angles[i] -= twoPi;
	}

	const vector3& referenceCenter = ringCenters[0];
	vector3 corners[4] = {
		vector3(terminalMetadata.boxMin[0], terminalMetadata.boxMin[1], 0.0),
		vector3(terminalMetadata.boxMax[0], terminalMetadata.boxMin[1], 0.0),
		vector3(terminalMetadata.boxMax[0], terminalMetadata.boxMax[1], 0.0),
		vector3(terminalMetadata.boxMin[0], terminalMetadata.boxMax[1], 0.0)
	};
	std::vector<number> cornerAngles(4);
	for (size_t c = 0; c < 4; ++c)
	{
		cornerAngles[c] = std::atan2(corners[c][1] - referenceCenter[1],
		                               corners[c][0] - referenceCenter[0]);
		if (cornerAngles[c] < 0.0) cornerAngles[c] += twoPi;
	}

	// Abstract 2-D node indices are identical at every axial station. The first
	// layer reuses Basolateral vertices; later nodes are newly created Inter
	// vertices. Corner nodes are appended unless an angular ray already hits the
	// corner exactly.
	std::vector<size_t> cornerNode(4, std::numeric_limits<size_t>::max());
	std::vector<int> cornerSector(4, -1);
	const size_t regularNodeCount = layerCount * angularCount;
	size_t abstractNodeCount = regularNodeCount;
	for (size_t c = 0; c < 4; ++c)
	{
		for (size_t i = 0; i < angularCount; ++i)
		{
			const number angularDistance = std::fabs(std::atan2(
				std::sin(cornerAngles[c] - angles[i]),
				std::cos(cornerAngles[c] - angles[i])));
			if (angularDistance < 1e-8)
			{
				cornerNode[c] = outerLayer * angularCount + i;
				break;
			}
		}
		if (cornerNode[c] != std::numeric_limits<size_t>::max()) continue;
		cornerNode[c] = abstractNodeCount++;
		for (size_t i = 0; i < angularCount; ++i)
		{
			const size_t j = (i + 1) % angularCount;
			number a = angles[i];
			number b = angles[j];
			if (j == 0 || b <= a) b += twoPi;
			number x = cornerAngles[c];
			if (x <= a) x += twoPi;
			if (x > a + 1e-8 && x < b - 1e-8)
			{
				cornerSector[c] = (int)i;
				break;
			}
		}
		UG_COND_THROW(cornerSector[c] < 0,
		              "Could not place a rectangular OuterWall corner between "
		              "two prismatic angular rays.");
	}

	std::vector<std::vector<Vertex*> > stationNodes(
		rings.size(), std::vector<Vertex*>(abstractNodeCount, NULL));
	for (size_t s = 0; s < rings.size(); ++s)
	{
		const vector3& center = ringCenters[s];
		for (size_t i = 0; i < angularCount; ++i)
		{
			stationNodes[s][i] = rings[s][i];
			vector3 direction;
			VecSubtract(direction, aaPos[rings[s][i]], center);
			direction[2] = 0.0;
			const number directionLength = VecLength(direction);
			UG_COND_THROW(directionLength <= 0.0,
			              "Zero-radius Basolateral vertex in prismatic Inter.");
			direction *= 1.0 / directionLength;
			number hitDistance = std::numeric_limits<number>::max();
			if (direction[0] > 1e-14)
				hitDistance = std::min(hitDistance,
					(terminalMetadata.boxMax[0] - center[0]) / direction[0]);
			else if (direction[0] < -1e-14)
				hitDistance = std::min(hitDistance,
					(terminalMetadata.boxMin[0] - center[0]) / direction[0]);
			if (direction[1] > 1e-14)
				hitDistance = std::min(hitDistance,
					(terminalMetadata.boxMax[1] - center[1]) / direction[1]);
			else if (direction[1] < -1e-14)
				hitDistance = std::min(hitDistance,
					(terminalMetadata.boxMin[1] - center[1]) / direction[1]);
			UG_COND_THROW(!std::isfinite(hitDistance) ||
			              hitDistance <= directionLength,
			              "OuterWall ray does not extend beyond Basolateral.");

			for (size_t layer = 1; layer < layerCount; ++layer)
			{
				RegularVertex* vertex = *g.create<RegularVertex>();
				const number distance = directionLength +
					radialFractions[layer - 1] *
					(hitDistance - directionLength);
				VecScaleAdd(aaPos[vertex], 1.0, center, distance, direction);
				aaSurfParams[vertex] = aaSurfParams[rings[s][i]];
				aaSurfParams[vertex].radial =
					(float)(distance / directionLength);
				sh.assign_subset(vertex, 5);
				psh.assign_subset(vertex, 1);
				stationNodes[s][layer * angularCount + i] = vertex;
			}
		}
		for (size_t c = 0; c < 4; ++c)
		{
			if (cornerNode[c] < regularNodeCount) continue;
			RegularVertex* vertex = *g.create<RegularVertex>();
			aaPos[vertex] = vector3(corners[c][0], corners[c][1], center[2]);
			aaSurfParams[vertex] = aaSurfParams[rings[s][0]];
			aaSurfParams[vertex].angular = (float)cornerAngles[c];
			aaSurfParams[vertex].radial = 2.0f;
			sh.assign_subset(vertex, 5);
			psh.assign_subset(vertex, 1);
			stationNodes[s][cornerNode[c]] = vertex;
		}
	}

	typedef std::array<size_t,3> PrismTriangle;
	std::vector<PrismTriangle> crossSectionTriangles;
	for (size_t layer = 0; layer + 1 < layerCount; ++layer)
	{
		const size_t inner = layer * angularCount;
		const size_t outer = (layer + 1) * angularCount;
		for (size_t i = 0; i < angularCount; ++i)
		{
			const size_t j = (i + 1) % angularCount;
			crossSectionTriangles.push_back(PrismTriangle{{
				inner + i, outer + i, inner + j}});
			crossSectionTriangles.push_back(PrismTriangle{{
				inner + j, outer + i, outer + j}});
		}
	}
	for (size_t c = 0; c < 4; ++c)
	{
		if (cornerSector[c] < 0) continue;
		const size_t i = (size_t)cornerSector[c];
		const size_t j = (i + 1) % angularCount;
		crossSectionTriangles.push_back(PrismTriangle{{
			outerLayer * angularCount + i, cornerNode[c],
			outerLayer * angularCount + j}});
	}

	// Normalize every 2-D triangle to counter-clockwise order. Sweeping the same
	// oriented connectivity through consecutive z stations gives one conforming
	// triangular prism per cross-section triangle and preserves each original
	// Basolateral quad as a rectangular prism face.
	for (size_t t = 0; t < crossSectionTriangles.size(); ++t)
	{
		PrismTriangle& tri = crossSectionTriangles[t];
		const vector3& a = aaPos[stationNodes[0][tri[0]]];
		const vector3& b = aaPos[stationNodes[0][tri[1]]];
		const vector3& c = aaPos[stationNodes[0][tri[2]]];
		const number area2 = (b[0]-a[0])*(c[1]-a[1])
			- (b[1]-a[1])*(c[0]-a[0]);
		UG_COND_THROW(std::fabs(area2) <= straightTolerance * straightTolerance,
		              "Degenerate 2-D triangle in prismatic Inter cross-section.");
		if (area2 < 0.0) std::swap(tri[1], tri[2]);
	}

	std::vector<Volume*> interPrisms;
	for (size_t s = 0; s + 1 < stationNodes.size(); ++s)
		for (size_t t = 0; t < crossSectionTriangles.size(); ++t)
		{
			const PrismTriangle& tri = crossSectionTriangles[t];
			Prism* prism = *g.create<Prism>(PrismDescriptor(
				stationNodes[s][tri[0]], stationNodes[s][tri[1]],
				stationNodes[s][tri[2]], stationNodes[s+1][tri[0]],
				stationNodes[s+1][tri[1]], stationNodes[s+1][tri[2]]));
			sh.assign_subset(prism, 5);
			psh.assign_subset(prism, 1);
			interPrisms.push_back(prism);
		}
	FixOrientation(g, interPrisms.begin(), interPrisms.end(), aaPos);

	for (size_t s = 0; s + 1 < rings.size(); ++s)
		for (size_t i = 0; i < angularCount; ++i)
		{
			const size_t j = (i + 1) % angularCount;
			FaceDescriptor fd;
			fd.set_num_vertices(4);
			fd.set_vertex(0, rings[s][i]);
			fd.set_vertex(1, rings[s][j]);
			fd.set_vertex(2, rings[s+1][j]);
			fd.set_vertex(3, rings[s+1][i]);
			Face* face = g.get_face(fd);
			UG_COND_THROW(!face || NumAssociatedVolumes(g, face) != 2,
			              "A prismatic Basolateral quadrilateral is not one ordinary "
			              "two-volume face.");
			sh.assign_subset(face, 2);
			psh.assign_subset(face, 0);
		}

	sh.set_subset_name("Membrane", 0);
	sh.set_subset_name("Lumen", 1);
	sh.set_subset_name("Basolateral", 2);
	sh.set_subset_name("Apical", 3);
	sh.set_subset_name("OuterWall", 4);
	sh.set_subset_name("Inter", 5);
	UG_LOGN("Added full swept-prism Inter for straight control: axial stations="
	        << rings.size() << ", angular vertices=" << angularCount
	        << ", cross-section triangles=" << crossSectionTriangles.size()
	        << ", prisms=" << interPrisms.size() << ".");
}


static void add_isolated_inter_to_version2_nephron
(
	Grid& g,
	SubsetHandler& sh,
	SubsetHandler& psh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
	number tetQuality,
	bool terminalPorts = false,
	bool conformingMembraneTransition = false
)
{
	/*
	 * EXTERIOR INTER TETRAHEDRALIZATION
	 * ---------------------------------
	 * The nested Lumen/Membrane tube created by create_neurite_with_er is the
	 * trusted, projector-compatible geometry. TetGen is deliberately NOT run on
	 * that grid because TetGen may reorder, insert, suppress, or move vertices.
	 * Earlier versions that tetrahedralized everything together corrupted the
	 * tube and generated overlapping Inter cells inside Lumen/Membrane.
	 *
	 * This routine therefore:
	 *   1. Copies only the Basolateral tube surface and OuterWall box into tg.
	 *   2. Runs TetGen on that isolated temporary grid.
	 *   3. Builds tetrahedral face adjacency and flood-fills from OuterWall,
	 *      treating every Basolateral triangle as an uncrossable separator.
	 *   4. Deletes all components not connected to OuterWall.
	 *   5. Copies the surviving exterior tetrahedra back as Inter.
	 *   6. Corrects implicitly created Inter faces and edges, which otherwise
	 *      inherit subset 0 (Membrane) from the default subset handler.
	 *
	 * Boundary faces/edges already marked Basolateral or OuterWall are preserved.
	 * The result is one conforming exterior volume without mutating the tube.
	 */
	UG_COND_THROW(tetQuality <= 0.0,
	              "Inter tetrahedral quality must be positive.");

	// The O-grid membrane is hexahedral, while TetGen necessarily triangulates
	// its exterior surface. Record the membrane cells now; after the exterior
	// Inter mesh has been copied back, each membrane hex is replaced by pyramids
	// and tetrahedra. Its Basolateral triangles can then be the very same face
	// objects used by Inter, so later refinement remains conforming.
	struct MembraneHexRecord
	{
		Volume* original;
		std::array<Vertex*, 8> vrts;
		std::array<std::array<Vertex*, 4>, 6> faces;
		std::array<int, 6> faceSubsets;
	};
	std::vector<MembraneHexRecord> membraneHexes;
	for (VolumeIterator vit = g.begin<Volume>(); vit != g.end<Volume>(); ++vit)
	{
		Volume* vol = *vit;
		if (sh.get_subset_index(vol) != 0) continue;
		UG_COND_THROW(vol->num_vertices() != 8 || vol->num_faces() != 6,
		              "The conforming O-grid path requires hexahedral membrane cells.");
		MembraneHexRecord rec;
		rec.original = vol;
		for (size_t i = 0; i < 8; ++i) rec.vrts[i] = vol->vertex(i);
		for (size_t i = 0; i < 6; ++i)
		{
			Face* face = g.get_face(vol, i);
			UG_COND_THROW(!face || face->num_vertices() != 4,
			              "Could not record a quadrilateral membrane face.");
			for (size_t j = 0; j < 4; ++j)
				rec.faces[i][j] = face->vertex(j);
			rec.faceSubsets[i] = sh.get_subset_index(face);
		}
		membraneHexes.push_back(rec);
	}

	// Record the exact pre-existing physical boundary objects. Pointer identity
	// is the reliable distinction here: TetGen later creates new triangles and
	// diagonals between boundary vertices, but those new objects are Inter sides,
	// not parts of the original Basolateral/OuterWall mesh.
	std::set<Face*> trustedBoundaryFaces;
	std::set<Edge*> trustedBoundaryEdges;
	std::set<Vertex*> trustedBoundaryVertices;
	for (FaceIterator fit = g.begin<Face>(); fit != g.end<Face>(); ++fit)
	{
		const int si = sh.get_subset_index(*fit);
		if (si != 2 && si != 4) continue;
		trustedBoundaryFaces.insert(*fit);
		Grid::traits<Edge>::secure_container edges;
		g.associated_elements(edges, *fit);
		for (size_t i = 0; i < edges.size(); ++i)
			trustedBoundaryEdges.insert(edges[i]);
		for (size_t i = 0; i < (*fit)->num_vertices(); ++i)
			trustedBoundaryVertices.insert((*fit)->vertex(i));
	}

	// TetGen operates only on this temporary grid. The projector-backed
	// version-2 Lumen/Membrane grid and all of its vertex positions stay intact.
	Grid tg;
	SubsetHandler tsh(tg);
	tsh.set_default_subset_index(0);
	tg.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> taaPos(tg, aPosition);
	std::map<Vertex*, Vertex*> origToTmp;
	std::map<Vertex*, Vertex*> tmpToOrig;

	auto tempVertex = [&] (Vertex* ov) -> Vertex*
	{
		std::map<Vertex*,Vertex*>::iterator it = origToTmp.find(ov);
		if (it != origToTmp.end()) return it->second;
		Vertex* tv = *tg.create<RegularVertex>();
		taaPos[tv] = aaPos[ov];
		origToTmp[ov] = tv;
		tmpToOrig[tv] = ov;
		return tv;
	};

	for (FaceIterator fit = g.begin<Face>(); fit != g.end<Face>(); ++fit)
	{
		Face* face = *fit;
		const int si = sh.get_subset_index(face);
		const bool isOuterWall = si == 4;
		const bool isNephronExterior = !isOuterWall
			&& NumAssociatedVolumes(g, face) == 1
			&& (!terminalPorts || si == 2);
		if (!isOuterWall && !isNephronExterior) continue;

		FaceDescriptor fd;
		fd.set_num_vertices((uint)face->num_vertices());
		for (size_t i = 0; i < face->num_vertices(); ++i)
			fd.set_vertex((uint)i, tempVertex(face->vertex(i)));
		const int boundarySubset = isOuterWall ? 4 : 2;
		if (face->num_vertices() == 3)
		{
			Face* tf = *tg.create<Triangle>(TriangleDescriptor(
				fd.vertex(0), fd.vertex(1), fd.vertex(2)));
			tsh.assign_subset(tf, boundarySubset);
		}
		else if (terminalPorts)
		{
			Face* first = *tg.create<Triangle>(TriangleDescriptor(
				fd.vertex(0), fd.vertex(1), fd.vertex(2)));
			Face* second = *tg.create<Triangle>(TriangleDescriptor(
				fd.vertex(0), fd.vertex(2), fd.vertex(3)));
			tsh.assign_subset(first, boundarySubset);
			tsh.assign_subset(second, boundarySubset);
		}
		else
		{
			Face* tf = *tg.create<Quadrilateral>(QuadrilateralDescriptor(
				fd.vertex(0), fd.vertex(1), fd.vertex(2), fd.vertex(3)));
			tsh.assign_subset(tf, boundarySubset);
		}
	}

	// The independently generated tube and wall patches need one globally
	// consistent winding across their shared portal rims before TetGen receives
	// the PLC.
	FixFaceOrientation(tg, tg.begin<Face>(), tg.end<Face>());

	// TetGen may segfault on an open/non-manifold PLC instead of returning an
	// error. Validate the temporary shell explicitly and report the defect at
	// the source boundary before invoking it.
	typedef std::pair<std::uintptr_t,std::uintptr_t> TmpEdgeKey;
	std::map<TmpEdgeKey,size_t> plcEdgeIncidence;
	for (FaceIterator fit = tg.begin<Face>(); fit != tg.end<Face>(); ++fit)
	{
		Face* face = *fit;
		for (size_t i = 0; i < face->num_vertices(); ++i)
		{
			std::uintptr_t a = reinterpret_cast<std::uintptr_t>(face->vertex(i));
			std::uintptr_t b = reinterpret_cast<std::uintptr_t>(
				face->vertex((i+1)%face->num_vertices()));
			if (b < a) std::swap(a,b);
			++plcEdgeIncidence[TmpEdgeKey(a,b)];
		}
	}
	size_t openPlcEdges = 0;
	size_t nonManifoldPlcEdges = 0;
	size_t reportedOpenPlcEdges = 0;
	for (std::map<TmpEdgeKey,size_t>::const_iterator it = plcEdgeIncidence.begin();
	     it != plcEdgeIncidence.end(); ++it)
	{
		if (it->second == 1)
		{
			++openPlcEdges;
			if (reportedOpenPlcEdges < 20)
			{
				Vertex* a = reinterpret_cast<Vertex*>(it->first.first);
				Vertex* b = reinterpret_cast<Vertex*>(it->first.second);
				UG_LOGN("Open Inter PLC edge " << reportedOpenPlcEdges
				        << ": " << taaPos[a] << " -> " << taaPos[b]);
				++reportedOpenPlcEdges;
			}
		}
		else if (it->second != 2) ++nonManifoldPlcEdges;
	}
	UG_COND_THROW(openPlcEdges || nonManifoldPlcEdges,
	              "Invalid Inter PLC before TetGen: " << openPlcEdges
	              << " open edges and " << nonManifoldPlcEdges
	              << " non-manifold edges.");
	UG_COND_THROW(!Tetrahedralize(tg, tsh, tetQuality, true, true, aPosition, 1),
	              "Failed to tetrahedralize the isolated exterior Inter grid.");

	// Optionally refine the tetrahedral field by a continuous distance-to-nephron
	// size rule.  The first TetGen pass supplies a valid conforming PLC mesh.  Its
	// Basolateral edge scale defines the near-field target, so this remains
	// independent of whether coordinates are stored in um or mm.  TetGen's
	// retetrahedralization consumes one maximum-volume value per existing tet and
	// inserts interior points without changing the trusted boundary vertices.
	if (g_adaptiveNephronMeshing.gradedInter)
	{
		std::set<Vertex*> tubeBoundaryVertices;
		std::vector<number> tubeBoundaryEdgeLengths;
		for (FaceIterator fit = tg.begin<Face>(); fit != tg.end<Face>(); ++fit)
		{
			Face* face = *fit;
			if (tsh.get_subset_index(face) != 2) continue;
			for (size_t i = 0; i < face->num_vertices(); ++i)
			{
				Vertex* a = face->vertex(i);
				Vertex* b = face->vertex((i + 1) % face->num_vertices());
				tubeBoundaryVertices.insert(a);
				tubeBoundaryVertices.insert(b);
				const number length = VecDistance(taaPos[a], taaPos[b]);
				if (length > 0.0) tubeBoundaryEdgeLengths.push_back(length);
			}
		}
		UG_COND_THROW(tubeBoundaryVertices.empty() || tubeBoundaryEdgeLengths.empty(),
		              "Graded Inter sizing found no Basolateral boundary geometry.");
		std::sort(tubeBoundaryEdgeLengths.begin(), tubeBoundaryEdgeLengths.end());
		const number surfaceScale = tubeBoundaryEdgeLengths[
			(tubeBoundaryEdgeLengths.size() - 1) / 4];
		const number nearDistance =
			g_adaptiveNephronMeshing.interNearDistanceFactor * surfaceScale;
		const number farDistance =
			g_adaptiveNephronMeshing.interFarDistanceFactor * surfaceScale;
		const number nearEdge =
			g_adaptiveNephronMeshing.interNearEdgeFactor * surfaceScale;
		const number farEdge =
			g_adaptiveNephronMeshing.interFarEdgeFactor * surfaceScale;

		ANumber aVolumeConstraint;
		tg.attach_to_volumes_dv(aVolumeConstraint, -1.0, true);
		Grid::VolumeAttachmentAccessor<ANumber> aaVolumeConstraint(
			tg, aVolumeConstraint);
		size_t constrainedTetrahedra = 0;
		for (VolumeIterator vit = tg.begin<Volume>(); vit != tg.end<Volume>(); ++vit)
		{
			Volume* volume = *vit;
			vector3 center(0.0);
			for (size_t i = 0; i < volume->num_vertices(); ++i)
				center += taaPos[volume->vertex(i)];
			center *= 1.0 / static_cast<number>(volume->num_vertices());
			number distanceToTube = std::numeric_limits<number>::max();
			for (std::set<Vertex*>::const_iterator it = tubeBoundaryVertices.begin();
			     it != tubeBoundaryVertices.end(); ++it)
				distanceToTube = std::min(distanceToTube,
				                              VecDistance(center, taaPos[*it]));
			if (distanceToTube >= farDistance)
			{
				aaVolumeConstraint[volume] = -1.0;
				continue;
			}
			const number blend = farDistance > nearDistance
				? std::max<number>(0.0, std::min<number>(1.0,
					(distanceToTube - nearDistance) / (farDistance - nearDistance)))
				: 1.0;
			const number targetEdge = nearEdge + blend * (farEdge - nearEdge);
			// A regular tetrahedron has volume h^3/(6*sqrt(2)).  Using 0.10
			// leaves a small safety margin below that ideal value.
			aaVolumeConstraint[volume] = 0.10 * targetEdge * targetEdge * targetEdge;
			++constrainedTetrahedra;
		}
		UG_LOGN("Graded Inter sizing: surface scale=" << surfaceScale
		        << ", distance=" << nearDistance << " -> " << farDistance
		        << ", target edge=" << nearEdge << " -> " << farEdge
		        << ", constrained initial tetrahedra=" << constrainedTetrahedra << ".");
		UG_COND_THROW(!Retetrahedralize(tg, tsh, aVolumeConstraint, tetQuality,
		                                    true, true, aPosition, true, 1),
		              "Failed graded Inter retetrahedralization.");

		// Report the achieved spatial grading, not only the requested size field.
		// This makes over-refinement visible in ordinary build logs and keeps the
		// defaults tunable across nephron radii and box sizes.
		size_t bandCount[3] = {0, 0, 0};
		number bandEdgeSum[3] = {0.0, 0.0, 0.0};
		for (VolumeIterator vit = tg.begin<Volume>(); vit != tg.end<Volume>(); ++vit)
		{
			Volume* volume = *vit;
			vector3 center(0.0);
			for (size_t i = 0; i < volume->num_vertices(); ++i)
				center += taaPos[volume->vertex(i)];
			center *= 1.0 / static_cast<number>(volume->num_vertices());
			number distanceToTube = std::numeric_limits<number>::max();
			for (std::set<Vertex*>::const_iterator it = tubeBoundaryVertices.begin();
			     it != tubeBoundaryVertices.end(); ++it)
				distanceToTube = std::min(distanceToTube,
				                              VecDistance(center, taaPos[*it]));
			const size_t band = distanceToTube < nearDistance ? 0
				: (distanceToTube < farDistance ? 1 : 2);
			number longestEdge = 0.0;
			for (size_t i = 0; i < volume->num_vertices(); ++i)
				for (size_t j = i + 1; j < volume->num_vertices(); ++j)
					longestEdge = std::max(longestEdge,
						VecDistance(taaPos[volume->vertex(i)], taaPos[volume->vertex(j)]));
			++bandCount[band];
			bandEdgeSum[band] += longestEdge;
		}
		UG_LOGN("Achieved Inter grading (near/transition/far): tetrahedra="
		        << bandCount[0] << "/" << bandCount[1] << "/" << bandCount[2]
		        << ", mean longest edge="
		        << (bandCount[0] ? bandEdgeSum[0]/bandCount[0] : 0.0) << "/"
		        << (bandCount[1] ? bandEdgeSum[1]/bandCount[1] : 0.0) << "/"
		        << (bandCount[2] ? bandEdgeSum[2]/bandCount[2] : 0.0) << ".");
		tg.detach_from_volumes(aVolumeConstraint);
	}

	// Do not delete tetrahedra using center/near-vertex point samples here.
	// TetGen has already recovered the Basolateral PLC facets.  Sampling close
	// to that curved interface is numerically ambiguous and used to delete valid
	// exterior tetrahedra, leaving thousands of artificial holes next to the
	// nephron.  Those one-sided Inter faces were subsequently mislabeled as
	// OuterWall.  The exact face-key flood fill below is the authoritative region
	// extraction: it cannot cross a recovered Basolateral facet and therefore
	// removes the sealed Lumen/Membrane-side component without damaging Inter.

	// Exact region extraction: retain only the tetrahedral component connected
	// to OuterWall, never crossing a Basolateral triangle. This removes the
	// entire inner TetGen component even when large tetrahedra defeat sampling.
	typedef std::array<std::uintptr_t,3> TmpTriKey;
	auto tmpTriKey = [] (Vertex* a, Vertex* b, Vertex* c) -> TmpTriKey
	{
		TmpTriKey key = {{reinterpret_cast<std::uintptr_t>(a),
		                  reinterpret_cast<std::uintptr_t>(b),
		                  reinterpret_cast<std::uintptr_t>(c)}};
		std::sort(key.begin(), key.end());
		return key;
	};
	std::set<TmpTriKey> basolateralKeys;
	std::set<TmpTriKey> outerWallKeys;
	for (FaceIterator fit = tg.begin<Face>(); fit != tg.end<Face>(); ++fit)
	{
		Face* face = *fit;
		if (face->num_vertices() != 3) continue;
		const int si = tsh.get_subset_index(face);
		const TmpTriKey key = tmpTriKey(face->vertex(0), face->vertex(1), face->vertex(2));
		if (si == 2) basolateralKeys.insert(key);
		else if (si == 4) outerWallKeys.insert(key);
	}
	std::map<TmpTriKey,std::vector<Volume*> > tmpFaceVolumes;
	for (VolumeIterator vit = tg.begin<Volume>(); vit != tg.end<Volume>(); ++vit)
	{
		Volume* vol = *vit;
		for (size_t omit=0; omit<4; ++omit)
		{
			Vertex* fv[3]; size_t k=0;
			for (size_t j=0; j<4; ++j) if (j != omit) fv[k++] = vol->vertex(j);
			tmpFaceVolumes[tmpTriKey(fv[0],fv[1],fv[2])].push_back(vol);
		}
	}
	Volume* outerSeed = NULL;
	for (std::set<TmpTriKey>::const_iterator it=outerWallKeys.begin();
	     it!=outerWallKeys.end() && !outerSeed; ++it)
	{
		const std::vector<Volume*>& vols = tmpFaceVolumes[*it];
		if (!vols.empty()) outerSeed = vols[0];
	}
	UG_COND_THROW(!outerSeed, "Could not find an Inter tetrahedron adjacent to OuterWall.");
	std::set<Volume*> exteriorComponent;
	std::vector<Volume*> componentStack(1, outerSeed);
	while (!componentStack.empty())
	{
		Volume* vol = componentStack.back(); componentStack.pop_back();
		if (!exteriorComponent.insert(vol).second) continue;
		for (size_t omit=0; omit<4; ++omit)
		{
			Vertex* fv[3]; size_t k=0;
			for (size_t j=0; j<4; ++j) if (j != omit) fv[k++] = vol->vertex(j);
			const TmpTriKey key = tmpTriKey(fv[0],fv[1],fv[2]);
			if (basolateralKeys.count(key)) continue;
			const std::vector<Volume*>& neighbors = tmpFaceVolumes[key];
			for (size_t j=0; j<neighbors.size(); ++j)
				if (!exteriorComponent.count(neighbors[j])) componentStack.push_back(neighbors[j]);
		}
	}
	std::vector<Volume*> nonExteriorVolumes;
	for (VolumeIterator vit=tg.begin<Volume>(); vit!=tg.end<Volume>(); ++vit)
		if (!exteriorComponent.count(*vit)) nonExteriorVolumes.push_back(*vit);
	tg.erase(nonExteriorVolumes.begin(), nonExteriorVolumes.end());

	std::map<Vertex*,Vertex*> tmpToDest = tmpToOrig;
	auto destVertex = [&] (Vertex* tv) -> Vertex*
	{
		std::map<Vertex*,Vertex*>::iterator it = tmpToDest.find(tv);
		if (it != tmpToDest.end()) return it->second;
		Vertex* nv = *g.create<RegularVertex>();
		aaPos[nv] = taaPos[tv];
		sh.assign_subset(nv, 5);
		psh.assign_subset(nv, 1);
		tmpToDest[tv] = nv;
		return nv;
	};

	std::vector<Volume*> copiedInter;
	for (VolumeIterator vit = tg.begin<Volume>(); vit != tg.end<Volume>(); ++vit)
	{
		Volume* tv = *vit;
		Tetrahedron* tet = *g.create<Tetrahedron>(TetrahedronDescriptor(
			destVertex(tv->vertex(0)), destVertex(tv->vertex(1)),
			destVertex(tv->vertex(2)), destVertex(tv->vertex(3))));
		sh.assign_subset(tet, 5);
		psh.assign_subset(tet, 1);
		copiedInter.push_back(tet);
	}
	FixOrientation(g, copiedInter.begin(), copiedInter.end(), aaPos);

	// Map TetGen's recovered Basolateral triangles to destination-grid vertex
	// keys. These keys select the same diagonal when the membrane side is split.
	typedef std::array<std::uintptr_t,3> DestTriKey;
	auto destTriKey = [] (Vertex* a, Vertex* b, Vertex* c) -> DestTriKey
	{
		DestTriKey key = {{reinterpret_cast<std::uintptr_t>(a),
		                   reinterpret_cast<std::uintptr_t>(b),
		                   reinterpret_cast<std::uintptr_t>(c)}};
		std::sort(key.begin(), key.end());
		return key;
	};
	std::set<DestTriKey> destinationBasolateralKeys;
	for (FaceIterator fit = tg.begin<Face>(); fit != tg.end<Face>(); ++fit)
	{
		Face* face = *fit;
		if (face->num_vertices() != 3 || tsh.get_subset_index(face) != 2)
			continue;
		destinationBasolateralKeys.insert(destTriKey(
			destVertex(face->vertex(0)), destVertex(face->vertex(1)),
			destVertex(face->vertex(2))));
	}

	if (conformingMembraneTransition)
	{
		// TetGen ran on the isolated copy, so every original tube vertex is
		// still alive. Replace only the Membrane hexahedra: non-Basolateral
		// quadrilateral sides become pyramid bases, while each Basolateral quad
		// is split along TetGen's recovered diagonal into two tetrahedra. The
		// resulting triangles are ordinary shared faces of Membrane and Inter.
		std::vector<Volume*> oldMembraneVolumes;
		for (size_t h = 0; h < membraneHexes.size(); ++h)
			oldMembraneVolumes.push_back(membraneHexes[h].original);
		g.erase(oldMembraneVolumes.begin(), oldMembraneVolumes.end());

		std::vector<Face*> oldBasolateralQuads;
		for (FaceIterator fit = sh.begin<Face>(2); fit != sh.end<Face>(2); ++fit)
			if ((*fit)->num_vertices() == 4)
				oldBasolateralQuads.push_back(*fit);
		g.erase(oldBasolateralQuads.begin(), oldBasolateralQuads.end());

		std::set<DestTriKey> actualBasolateralKeys;
		std::vector<Volume*> restoredMembraneVolumes;
		for (size_t h = 0; h < membraneHexes.size(); ++h)
		{
			const MembraneHexRecord& rec = membraneHexes[h];
			RegularVertex* center = *g.create<RegularVertex>();
			aaPos[center] = vector3(0.0);
			for (size_t i = 0; i < 8; ++i)
				aaPos[center] += aaPos[rec.vrts[i]];
			aaPos[center] *= 0.125;
			aaSurfParams[center] = aaSurfParams[rec.vrts[0]];
			aaSurfParams[center].axial = 0.0;
			aaSurfParams[center].radial = 0.0;
			for (size_t i = 0; i < 8; ++i)
			{
				aaSurfParams[center].axial +=
					0.125 * aaSurfParams[rec.vrts[i]].axial;
				aaSurfParams[center].radial +=
					0.125 * aaSurfParams[rec.vrts[i]].radial;
			}
			sh.assign_subset(center, 0);
			psh.assign_subset(center, 0);

			for (size_t f = 0; f < 6; ++f)
			{
				Vertex* a = rec.faces[f][0];
				Vertex* b = rec.faces[f][1];
				Vertex* c = rec.faces[f][2];
				Vertex* d = rec.faces[f][3];
				if (rec.faceSubsets[f] != 2)
				{
					Pyramid* pyramid = *g.create<Pyramid>(
						PyramidDescriptor(a,b,c,d,center));
					sh.assign_subset(pyramid, 0);
					psh.assign_subset(pyramid, 0);
					restoredMembraneVolumes.push_back(pyramid);
					continue;
				}

				const DestTriKey ac0 = destTriKey(a,b,c);
				const DestTriKey ac1 = destTriKey(a,c,d);
				const DestTriKey bd0 = destTriKey(a,b,d);
				const DestTriKey bd1 = destTriKey(b,c,d);
				const bool useAC = destinationBasolateralKeys.count(ac0)
					&& destinationBasolateralKeys.count(ac1);
				const bool useBD = destinationBasolateralKeys.count(bd0)
					&& destinationBasolateralKeys.count(bd1);
				UG_COND_THROW(useAC == useBD,
				              "Could not select one recovered diagonal for a "
				              "Basolateral quadrilateral.");

				Tetrahedron* t0;
				Tetrahedron* t1;
				if (useAC)
				{
					t0 = *g.create<Tetrahedron>(
						TetrahedronDescriptor(center,a,b,c));
					t1 = *g.create<Tetrahedron>(
						TetrahedronDescriptor(center,a,c,d));
					actualBasolateralKeys.insert(ac0);
					actualBasolateralKeys.insert(ac1);
				}
				else
				{
					t0 = *g.create<Tetrahedron>(
						TetrahedronDescriptor(center,a,b,d));
					t1 = *g.create<Tetrahedron>(
						TetrahedronDescriptor(center,b,c,d));
					actualBasolateralKeys.insert(bd0);
					actualBasolateralKeys.insert(bd1);
				}
				sh.assign_subset(t0, 0);
				sh.assign_subset(t1, 0);
				psh.assign_subset(t0, 0);
				psh.assign_subset(t1, 0);
				restoredMembraneVolumes.push_back(t0);
				restoredMembraneVolumes.push_back(t1);
			}
		}
		FixOrientation(g, restoredMembraneVolumes.begin(),
		               restoredMembraneVolumes.end(), aaPos);

		std::vector<Face*> basolateralFaces;
		for (FaceIterator fit = g.begin<Face>(); fit != g.end<Face>(); ++fit)
		{
			Face* face = *fit;
			if (face->num_vertices() != 3) continue;
			if (actualBasolateralKeys.count(destTriKey(
				face->vertex(0), face->vertex(1), face->vertex(2))))
				basolateralFaces.push_back(face);
		}
		UG_COND_THROW(basolateralFaces.size() != actualBasolateralKeys.size(),
		              "Recovered Basolateral triangle count does not match the "
		              "conforming Membrane transition.");
		std::vector<Face*> apicalFaces(sh.begin<Face>(3), sh.end<Face>(3));
		std::vector<Face*> outerFaces(sh.begin<Face>(4), sh.end<Face>(4));

		// Close each live volume explicitly, then restore physical-interface
		// precedence. TetGen never touched this destination grid, so these
		// associations contain no temporary-grid or erased-tube pointers.
		std::vector<Face*> closureFaces;
		std::vector<Edge*> closureEdges;
		for (VolumeIterator vit = g.begin<Volume>(); vit != g.end<Volume>(); ++vit)
		{
			Volume* vol = *vit;
			const int si = sh.get_subset_index(vol);
			CollectAssociated(closureFaces, g, vol);
			for (size_t f = 0; f < closureFaces.size(); ++f)
			{
				Face* face = closureFaces[f];
				sh.assign_subset(face, si);
				for (size_t j = 0; j < face->num_vertices(); ++j)
					sh.assign_subset(face->vertex(j), si);
				CollectAssociated(closureEdges, g, face);
				for (size_t e = 0; e < closureEdges.size(); ++e)
					sh.assign_subset(closureEdges[e], si);
			}
		}
		auto restoreBoundary = [&] (const std::vector<Face*>& faces, int si)
		{
			std::vector<Edge*> edges;
			for (size_t i = 0; i < faces.size(); ++i)
			{
				Face* face = faces[i];
				sh.assign_subset(face, si);
				psh.assign_subset(face, si == 4 ? 1 : 0);
				for (size_t j = 0; j < face->num_vertices(); ++j)
				{
					sh.assign_subset(face->vertex(j), si);
					psh.assign_subset(face->vertex(j), si == 4 ? 1 : 0);
				}
				CollectAssociated(edges, g, face);
				for (size_t j = 0; j < edges.size(); ++j)
				{
					sh.assign_subset(edges[j], si);
					psh.assign_subset(edges[j], si == 4 ? 1 : 0);
				}
			}
		};
		restoreBoundary(basolateralFaces, 2);
		restoreBoundary(apicalFaces, 3);
		restoreBoundary(outerFaces, 4);
		sh.set_subset_name("Membrane", 0);
		sh.set_subset_name("Lumen", 1);
		sh.set_subset_name("Basolateral", 2);
		sh.set_subset_name("Apical", 3);
		sh.set_subset_name("OuterWall", 4);
		sh.set_subset_name("Inter", 5);
		UG_LOGN("Added isolated, conforming exterior Inter with ordinary shared "
		        << "Basolateral triangles.");
		return;
	}

	// Keep the structured membrane hexahedra and their quadrilateral ring strips.
	// TetGen represents every Basolateral quad by two triangles.  Register those
	// triangles as constrained children of the original quad instead of replacing
	// the membrane hex by a center fan.  This distinction is essential for
	// cross-sectional refinement: refining a hex splits both neighboring rings
	// coherently and creates the desired midpoint-to-midpoint axial edge, whereas
	// refining two independent triangles connects each new midpoint diagonally to
	// an old coarse-ring corner.
	std::set<Face*> basolateralQuads;
	for (size_t h = 0; h < membraneHexes.size(); ++h)
	{
		for (size_t f = 0; f < 6; ++f)
		{
			if (membraneHexes[h].faceSubsets[f] != 2) continue;
			Face* face = g.get_face(membraneHexes[h].original, f);
			if (face) basolateralQuads.insert(face);
		}
	}

	std::map<DestTriKey, Face*> interBoundaryTriangles;
	for (FaceIterator fit = g.begin<Face>(); fit != g.end<Face>(); ++fit)
	{
		Face* face = *fit;
		if (face->num_vertices() != 3) continue;
		const DestTriKey key = destTriKey(
			face->vertex(0), face->vertex(1), face->vertex(2));
		if (!destinationBasolateralKeys.count(key)) continue;
		UG_COND_THROW(interBoundaryTriangles.count(key),
		              "Duplicate recovered Basolateral triangle in destination grid.");
		interBoundaryTriangles[key] = face;
	}

	size_t constrainedBasolateralPairs = 0;
	for (std::set<Face*>::iterator it = basolateralQuads.begin();
	     it != basolateralQuads.end(); ++it)
	{
		Face* oldQuad = *it;
		Vertex* a = oldQuad->vertex(0);
		Vertex* b = oldQuad->vertex(1);
		Vertex* c = oldQuad->vertex(2);
		Vertex* d = oldQuad->vertex(3);
		DestTriKey k0, k1;
		if (destinationBasolateralKeys.count(destTriKey(a,b,c)) &&
		    destinationBasolateralKeys.count(destTriKey(a,c,d)))
		{
			k0 = destTriKey(a,b,c); k1 = destTriKey(a,c,d);
		}
		else if (destinationBasolateralKeys.count(destTriKey(a,b,d)) &&
		         destinationBasolateralKeys.count(destTriKey(b,c,d)))
		{
			k0 = destTriKey(a,b,d); k1 = destTriKey(b,c,d);
		}
		else UG_THROW("Could not match a Basolateral quad to TetGen's two triangles.");

		UG_COND_THROW(!interBoundaryTriangles.count(k0) ||
		              !interBoundaryTriangles.count(k1),
		              "Recovered Basolateral triangles are missing from the Inter grid.");
		trustedBoundaryFaces.erase(oldQuad);
		ConstrainingQuadrilateral* quad =
			*g.create_and_replace<ConstrainingQuadrilateral>(oldQuad);
		sh.assign_subset(quad, 2);
		psh.assign_subset(quad, 0);
		trustedBoundaryFaces.insert(quad);

		const DestTriKey pairKeys[2] = {k0, k1};
		for (size_t t = 0; t < 2; ++t)
		{
			Face* oldTri = interBoundaryTriangles[pairKeys[t]];
			ConstrainedTriangle* tri =
				*g.create_and_replace<ConstrainedTriangle>(oldTri);
			tri->set_constraining_object(quad);
			quad->add_constrained_object(tri);
			sh.assign_subset(tri, 2);
			psh.assign_subset(tri, 1);
			trustedBoundaryFaces.insert(tri);
			interBoundaryTriangles[pairKeys[t]] = tri;
		}
		++constrainedBasolateralPairs;
	}
	UG_LOGN("Preserved " << membraneHexes.size()
	        << " structured Membrane hexahedra and linked "
	        << constrainedBasolateralPairs
	        << " Basolateral quads to their Inter triangle pairs.");

	// Faces created implicitly with the copied tetrahedra inherit the subset
	// handler's default subset (0 == Membrane).  Correct their classification:
	// every ordinary face belonging to an Inter tetrahedron is Inter.  Preserve
	// the explicitly constructed physical boundaries, which must remain named
	// Basolateral and OuterWall.
	std::set<Face*> interFaces;
	for (size_t i = 0; i < copiedInter.size(); ++i)
	{
		Grid::traits<Face>::secure_container faces;
		g.associated_elements(faces, copiedInter[i]);
		for (size_t j = 0; j < faces.size(); ++j)
			interFaces.insert(faces[j]);
	}
	for (std::set<Face*>::iterator it = interFaces.begin();
	     it != interFaces.end(); ++it)
	{
		Face* face = *it;
		if (!trustedBoundaryFaces.count(face))
		{
			sh.assign_subset(face, 5);
			psh.assign_subset(face, 1);
		}
	}

	// Apply the same correction to the edge closure of Inter. Implicitly
	// created tetrahedron edges otherwise retain subset 0 (Membrane). Do not
	// trust an edge's inherited subset here: an internal edge can inherit 2 or 4
	// from a temporary boundary construction. Preserve it as a physical boundary
	// only when it is actually incident to a Basolateral or OuterWall face.
	std::set<Edge*> interEdges;
	for (std::set<Face*>::iterator it = interFaces.begin();
	     it != interFaces.end(); ++it)
	{
		Grid::traits<Edge>::secure_container edges;
		g.associated_elements(edges, *it);
		for (size_t j = 0; j < edges.size(); ++j)
			interEdges.insert(edges[j]);
	}
	for (std::set<Edge*>::iterator it = interEdges.begin();
	     it != interEdges.end(); ++it)
	{
		if (trustedBoundaryEdges.count(*it)) continue;
		sh.assign_subset(*it, 5);
		psh.assign_subset(*it, 1);
	}

	// Vertices need an explicit closure pass as well. Without this pass, vertices
	// introduced by TetGen or inherited from implicitly created Inter sides can
	// remain in subset 0 (Membrane). On refinement those stale vertex subsets are
	// propagated and appear as isolated Membrane points and edges around bends.
	std::set<Vertex*> interVertices;
	for (size_t i = 0; i < copiedInter.size(); ++i)
		for (size_t j = 0; j < copiedInter[i]->num_vertices(); ++j)
			interVertices.insert(copiedInter[i]->vertex(j));
	for (std::set<Vertex*>::iterator it = interVertices.begin();
	     it != interVertices.end(); ++it)
	{
		if (trustedBoundaryVertices.count(*it)) continue;
		sh.assign_subset(*it, 5);
		psh.assign_subset(*it, 1);
	}

	UG_LOGN("Inter closure classified from exact pre-existing boundary topology: "
	        << interFaces.size() << " faces, " << interEdges.size()
	        << " edges, " << interVertices.size() << " vertices checked.");
	sh.set_subset_name("Membrane", 0);
	sh.set_subset_name("Lumen", 1);
	sh.set_subset_name("Basolateral", 2);
	sh.set_subset_name("Apical", 3);
	sh.set_subset_name("OuterWall", 4);
	sh.set_subset_name("Inter", 5);
}


static void cleanup_completed_nephron_mesh_before_final_subsets
(
	Grid& g,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	number mergeTolerance
)
{
	UG_COND_THROW(mergeTolerance <= 0.0,
	              "Mesh cleanup tolerance must be positive.");

	const size_t v0 = g.num<Vertex>();
	const size_t e0 = g.num<Edge>();
	const size_t f0 = g.num<Face>();
	const size_t c0 = g.num<Volume>();

	// Orient the completed volume mesh and its explicit sides.  The second pass
	// below is intentional: merging coincident vertices and deleting duplicate
	// sides can change adjacency, so orientation is validated again afterwards.
	const int fixedBefore =
		FixOrientation(g, g.begin<Volume>(), g.end<Volume>(), aaPos);
	FixFaceOrientation(g, g.begin<Face>(), g.end<Face>());
	AdjustEdgeOrientationToFaceOrientation(g, g.begin<Edge>(), g.end<Edge>());

	// Merge geometrically coincident vertices.  Do not call RemoveDuplicates on
	// faces or edges after volumes exist: Grid::erase(Face*) also erases volumes
	// incident to that face.  On the nephron/Inter handoff that old cleanup
	// removed valid tetrahedra and opened holes beside the membrane.  TetGen has
	// already removed duplicate PLC sides before tetrahedralization, while the
	// orphan-constraint sweep in assign_indexed_nephron_subsets safely removes
	// surface objects that have no adjacent volume.
	RemoveDoubles<3>(g, g.begin<Vertex>(), g.end<Vertex>(), aaPos,
	                 mergeTolerance);
	UG_COND_THROW(g.num<Volume>() != c0,
	              "Merging coincident vertices changed the volume count from "
	              << c0 << " to " << g.num<Volume>()
	              << "; refusing to continue with an opened volume mesh.");

	const int fixedAfter =
		FixOrientation(g, g.begin<Volume>(), g.end<Volume>(), aaPos);
	FixFaceOrientation(g, g.begin<Face>(), g.end<Face>());
	AdjustEdgeOrientationToFaceOrientation(g, g.begin<Edge>(), g.end<Edge>());

	UG_LOGN("Pre-subset mesh cleanup (tolerance=" << mergeTolerance << "): "
	        << "vertices " << v0 << " -> " << g.num<Vertex>() << ", edges "
	        << e0 << " -> " << g.num<Edge>() << ", faces " << f0 << " -> "
	        << g.num<Face>() << ", volumes " << c0 << " -> "
	        << g.num<Volume>() << ", volume orientations fixed "
	        << fixedBefore << "+" << fixedAfter << ".");
}


static void validate_and_repair_orientation_after_projected_refinement
(
	Grid& g,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	ISubsetHandler& sh,
	size_t refinementLevel,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >*
		aaSurfParams = NULL
)
{
	const int corrected =
		FixOrientation(g, g.begin<Volume>(), g.end<Volume>(), aaPos);

	size_t indeterminate = 0;
	std::map<int, size_t> indeterminateBySubset;
	for (VolumeIterator it = g.begin<Volume>(); it != g.end<Volume>(); ++it)
	{
		Volume* vol = *it;
		const bool beforeFlip = CheckOrientation(vol, aaPos);
		g.flip_orientation(vol);
		const bool afterFlip = CheckOrientation(vol, aaPos);
		g.flip_orientation(vol);
		if (beforeFlip == afterFlip)
		{
			++indeterminate;
			++indeterminateBySubset[sh.get_subset_index(vol)];
			if (aaSurfParams && indeterminate <= 30)
			{
				vector3 center(0.0, 0.0, 0.0);
				float minAxial = std::numeric_limits<float>::max();
				float maxAxial = -std::numeric_limits<float>::max();
				float minRadial = std::numeric_limits<float>::max();
				float maxRadial = -std::numeric_limits<float>::max();
				for (size_t v = 0; v < vol->num_vertices(); ++v)
				{
					VecAdd(center, center, aaPos[vol->vertex(v)]);
					const NeuriteProjector::SurfaceParams& sp =
						(*aaSurfParams)[vol->vertex(v)];
					minAxial = std::min(minAxial, sp.axial);
					maxAxial = std::max(maxAxial, sp.axial);
					minRadial = std::min(minRadial, sp.radial);
					maxRadial = std::max(maxRadial, sp.radial);
				}
				VecScale(center, center, 1.0 / vol->num_vertices());
				UG_LOGN("  degenerate detail: subset=" << sh.get_subset_index(vol)
				        << ", vertices=" << vol->num_vertices() << ", center="
				        << center << ", axial=[" << minAxial << "," << maxAxial
				        << "], radial=[" << minRadial << "," << maxRadial << "].");
			}
		}
	}
	for (std::map<int, size_t>::const_iterator it = indeterminateBySubset.begin();
	     it != indeterminateBySubset.end(); ++it)
		UG_LOGN("  degenerate volume subset " << it->first << " ("
		        << sh.get_subset_name(it->first) << "): " << it->second);

	UG_LOGN("Projected refinement " << refinementLevel
	        << " orientation check: corrected " << corrected
	        << " volumes; indeterminate/degenerate volumes: "
	        << indeterminate << ".");
	UG_COND_THROW(indeterminate != 0,
	              "Projected refinement " << refinementLevel << " created "
	              << indeterminate
	              << " volumes whose orientation cannot be determined.");
}


static size_t mark_nephron_cross_section_volumes
(
	HangingNodeRefiner_MultiGrid& refiner,
	MultiGrid& grid,
	ISubsetHandler& sh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >&
		aaSurfParams,
	size_t nephronCount,
	bool allowMembraneRefinement,
	std::set<Edge*>& transverseEdges,
	std::set<Edge*>& axialEdges
)
{
	const number axialSectionTolerance = 1e-6;
	const int outerWallSubset = (int)(6 * nephronCount);
	transverseEdges.clear();
	axialEdges.clear();
	size_t markedVolumes = 0;
	size_t skippedTransitionVolumes = 0;
	size_t markedMembraneTransitions = 0;
	size_t skippedWarpedMembranePyramids = 0;
	size_t skippedWithoutTransverseEdges = 0;
	for (size_t level = 0; level <= grid.top_level(); ++level)
	for (VolumeIterator it = grid.begin<Volume>(level);
	     it != grid.end<Volume>(level); ++it)
	{
		Volume* volume = *it;
		if (grid.has_children(volume)) continue;
		const int si = sh.get_subset_index(volume);
		if (si < 0 || si >= outerWallSubset) continue;
		const int localSI = si % 6;
		if (localSI != 0 && localSI != 1) continue; // Membrane or Lumen volume
		int localMark = 0;
		std::set<Edge*> localTransverseEdges;
		std::set<Edge*> localAxialEdges;
		EdgeDescriptor edge;
		for (size_t e = 0; e < volume->num_edges(); ++e)
		{
			volume->edge_desc(e, edge);
			const number a0 = aaSurfParams[edge.vertex(0)].axial;
			const number a1 = aaSurfParams[edge.vertex(1)].axial;
			Edge* gridEdge = grid.get_edge(volume, e);
			if (std::fabs(a0 - a1) < axialSectionTolerance)
			{
				localTransverseEdges.insert(gridEdge);
				localMark |= (1 << e);
			}
			else
				localAxialEdges.insert(gridEdge);
		}
		if (localMark == 0)
		{
			// A protected transition parent can leave active children whose
			// available edges are all axial. They are not candidates for this
			// cross-sectional pass and must simply remain at their current level.
			++skippedWithoutTransverseEdges;
			continue;
		}
		// Reconstructed Membrane cells are pyramids/tetrahedra. Their local
		// rules are recursive rather than advertised as regular, but with valid
		// neurite parameters on the reconstruction center their transverse edge
		// masks can be refined without introducing an axial section.
		bool unsafeWarpedPyramid = false;
		if (localSI == 0 && dynamic_cast<Pyramid*>(volume) != NULL)
		{
			vector3 ab, ac, ad, normal;
			VecSubtract(ab, aaPos[volume->vertex(1)], aaPos[volume->vertex(0)]);
			VecSubtract(ac, aaPos[volume->vertex(2)], aaPos[volume->vertex(0)]);
			VecSubtract(ad, aaPos[volume->vertex(3)], aaPos[volume->vertex(0)]);
			VecCross(normal, ab, ac);
			const number scale = std::max(VecLength(ab), VecLength(ac));
			const number baseWarp = VecLength(normal) > 0.0 && scale > 0.0
				? std::fabs(VecProd(ad, normal)) / (VecLength(normal) * scale)
				: std::numeric_limits<number>::max();
			number baseMinAxial = aaSurfParams[volume->vertex(0)].axial;
			number baseMaxAxial = baseMinAxial;
			for (size_t bv = 1; bv < 4; ++bv)
			{
				const number a = aaSurfParams[volume->vertex(bv)].axial;
				baseMinAxial = std::min(baseMinAxial, a);
				baseMaxAxial = std::max(baseMaxAxial, a);
			}
			const bool terminalBase = baseMaxAxial - baseMinAxial < axialSectionTolerance &&
				(std::fabs(baseMinAxial) < axialSectionTolerance ||
				 std::fabs(baseMinAxial - 1.0) < axialSectionTolerance);
			const bool coarseExactProjectionRisk = grid.get_level(volume) == 0 &&
				((localMark == 15 && terminalBase) ||
				 (localMark == 5 && baseWarp > 0.04));
			const bool refinedFullSplitRisk = grid.get_level(volume) > 0 &&
				localMark == 255 && baseWarp > 0.05;
			unsafeWarpedPyramid = coarseExactProjectionRisk || refinedFullSplitRisk;
		}
		if (unsafeWarpedPyramid)
		{
			// UG4's full split of a pyramid with a strongly warped base can
			// create a zero-volume child. Keep this isolated transition cell at
			// its current level while refining the surrounding Membrane cells.
			++skippedWarpedMembranePyramids;
			continue;
		}
		const bool firstCrossSectionLevel = allowMembraneRefinement;
		const bool membraneTransition = localSI == 0 && firstCrossSectionLevel &&
			(dynamic_cast<Tetrahedron*>(volume) != NULL ||
			 dynamic_cast<Pyramid*>(volume) != NULL);
		if (!membraneTransition && !volume->is_regular_ref_rule(localMark))
		{
			++skippedTransitionVolumes;
			continue;
		}
		transverseEdges.insert(localTransverseEdges.begin(),
		                       localTransverseEdges.end());
		axialEdges.insert(localAxialEdges.begin(), localAxialEdges.end());
		// Supply the matching local mark on every face explicitly.  Inferring
		// face marks from the Basolateral constraining-quad/triangle relation can
		// otherwise promote a structured hex to RM_FULL.  RM_FULL asks for a
		// center on every axial side face and recreates diagonal connections to
		// coarse corners.  Here cross-sectional faces receive a full edge mask,
		// while axial strip faces receive only their two opposite transverse
		// edges, producing one straight midpoint-to-midpoint child edge.
		for (size_t f = 0; f < volume->num_faces(); ++f)
		{
			Face* face = grid.get_face(volume, f);
			if (!face) continue;
			int faceMark = 0;
			Grid::traits<Edge>::secure_container faceEdges;
			grid.associated_elements_sorted(faceEdges, face);
			for (size_t fe = 0; fe < faceEdges.size(); ++fe)
			{
				const number a0 = aaSurfParams[faceEdges[fe]->vertex(0)].axial;
				const number a1 = aaSurfParams[faceEdges[fe]->vertex(1)].axial;
				if (std::fabs(a0 - a1) < axialSectionTolerance)
					faceMark |= (1 << fe);
			}
			if (faceMark) refiner.mark_local(face, faceMark);
		}
		refiner.mark_local(volume, localMark);
		++markedVolumes;
		if (membraneTransition) ++markedMembraneTransitions;
	}

	UG_COND_THROW(transverseEdges.empty(),
	              "Cross-sectional refinement found no active transverse edges.");
	UG_LOGN("Cross-sectional volume masks: marked " << markedVolumes
	        << " cells, including " << markedMembraneTransitions
	        << " reconstructed Membrane pyramid/tetrahedron cells; retained "
	        << skippedTransitionVolumes
	        << " transition cells at their current level ("
	        << skippedWarpedMembranePyramids
	        << " warped Membrane pyramids and "
	        << skippedWithoutTransverseEdges
	        << " cells without transverse edges excluded)." );
	return transverseEdges.size();
}


static std::map<uint32, std::vector<number> > collect_nephron_axial_sections
(
	MultiGrid& grid,
	ISubsetHandler& sh,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >&
		aaSurfParams,
	size_t nephronCount
)
{
	const number axialSectionTolerance = 1e-6;
	const int outerWallSubset = (int)(6 * nephronCount);
	std::map<uint32, std::vector<number> > sections;
	std::set<Vertex*> tubeVertices;
	for (size_t level = 0; level <= grid.top_level(); ++level)
	for (VolumeIterator it = grid.begin<Volume>(level);
	     it != grid.end<Volume>(level); ++it)
	{
		if (grid.has_children(*it)) continue;
		const int si = sh.get_subset_index(*it);
		if (si < 0 || si >= outerWallSubset) continue;
		const int localSI = si % 6;
		if (localSI != 0 && localSI != 1) continue;
		for (size_t v = 0; v < (*it)->num_vertices(); ++v)
			tubeVertices.insert((*it)->vertex(v));
	}
	for (std::set<Vertex*>::iterator it = tubeVertices.begin();
	     it != tubeVertices.end(); ++it)
	{
		const NeuriteProjector::SurfaceParams& sp = aaSurfParams[*it];
		const uint32 neuriteID = sp.neuriteID & ((1 << 20) - 1);
		std::vector<number>& values = sections[neuriteID];
		bool known = false;
		for (size_t i = 0; i < values.size(); ++i)
			if (std::fabs(values[i] - sp.axial) < axialSectionTolerance)
			{
				known = true;
				break;
			}
		if (!known) values.push_back(sp.axial);
	}
	return sections;
}


static size_t count_vertices_on_new_axial_sections
(
	MultiGrid& grid,
	ISubsetHandler& sh,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >&
		aaSurfParams,
	size_t nephronCount,
	const std::map<uint32, std::vector<number> >& allowedSections
)
{
	const number axialSectionTolerance = 1e-6;
	const int outerWallSubset = (int)(6 * nephronCount);
	std::set<Vertex*> unexpected;
	for (size_t level = 0; level <= grid.top_level(); ++level)
	for (VolumeIterator it = grid.begin<Volume>(level);
	     it != grid.end<Volume>(level); ++it)
	{
		if (grid.has_children(*it)) continue;
		const int si = sh.get_subset_index(*it);
		if (si < 0 || si >= outerWallSubset) continue;
		const int localSI = si % 6;
		if (localSI != 0 && localSI != 1) continue;
		for (size_t v = 0; v < (*it)->num_vertices(); ++v)
		{
			Vertex* vertex = (*it)->vertex(v);
			const NeuriteProjector::SurfaceParams& sp = aaSurfParams[vertex];
			const uint32 neuriteID = sp.neuriteID & ((1 << 20) - 1);
			std::map<uint32, std::vector<number> >::const_iterator found =
				allowedSections.find(neuriteID);
			bool known = false;
			number nearestDifference = std::numeric_limits<number>::max();
			if (found != allowedSections.end())
				for (size_t j = 0; j < found->second.size(); ++j)
				{
					nearestDifference = std::min(nearestDifference,
						std::fabs(found->second[j] - sp.axial));
					if (std::fabs(found->second[j] - sp.axial) < axialSectionTolerance)
					{
						known = true;
						break;
					}
				}
			if (!known && unexpected.insert(vertex).second && unexpected.size() <= 8)
				UG_LOGN("  new axial sample: neuriteID=" << neuriteID
				        << ", axial=" << sp.axial
				        << ", nearest prior difference=" << nearestDifference);
		}
	}
	return unexpected.size();
}


static size_t restore_active_basolateral_children
(
	MultiGrid& grid,
	ISubsetHandler& sh,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >&
		aaSurfParams,
	size_t nephronCount,
	size_t refinementLevel
)
{
	// Hanging-node refinement creates the outer children as sides of Membrane
	// volumes. Without this pass those faces inherit the Membrane volume subset,
	// even though their parent is the Basolateral constraint. Geometrically this
	// interface is unambiguous: all of its vertices have the outer neurite radial
	// parameter, whereas internal Membrane faces contain at least one smaller
	// radial value.
	const int outerWallSubset = (int)(6 * nephronCount);
	std::vector<std::set<Edge*> > boundaryEdges(nephronCount);
	std::vector<std::set<Vertex*> > boundaryVertices(nephronCount);
	size_t restoredFaces = 0;
	for (size_t level = 0; level <= grid.top_level(); ++level)
	for (FaceIterator it = grid.begin<Face>(level);
	     it != grid.end<Face>(level); ++it)
	{
		Face* face = *it;
		if (grid.has_children(face)) continue;
		bool outerRadial = true;
		for (size_t v = 0; v < face->num_vertices(); ++v)
			outerRadial = outerRadial &&
				aaSurfParams[face->vertex(v)].radial > 1.0 - 1e-6;
		if (!outerRadial) continue;

		int owner = -1;
		Grid::traits<Volume>::secure_container volumes;
		grid.associated_elements(volumes, face);
		for (size_t v = 0; v < volumes.size(); ++v)
		{
			const int volumeSubset = sh.get_subset_index(volumes[v]);
			if (volumeSubset < 0 || volumeSubset >= outerWallSubset ||
			    volumeSubset % 6 != 0)
				continue;
			const int candidate = volumeSubset / 6;
			UG_COND_THROW(owner >= 0 && owner != candidate,
			              "Outer Membrane face touches different nephrons.");
			owner = candidate;
		}
		if (owner < 0) continue;

		const int basolateralSubset = 6 * owner + 2;
		if (sh.get_subset_index(face) != basolateralSubset)
		{
			sh.assign_subset(face, basolateralSubset);
			++restoredFaces;
		}
		Grid::traits<Edge>::secure_container edges;
		grid.associated_elements(edges, face);
		for (size_t e = 0; e < edges.size(); ++e)
		{
			boundaryEdges[(size_t)owner].insert(edges[e]);
			boundaryVertices[(size_t)owner].insert(edges[e]->vertex(0));
			boundaryVertices[(size_t)owner].insert(edges[e]->vertex(1));
		}
	}

	for (size_t owner = 0; owner < nephronCount; ++owner)
	{
		const int subset = (int)(6 * owner + 2);
		for (std::set<Edge*>::iterator it = boundaryEdges[owner].begin();
		     it != boundaryEdges[owner].end(); ++it)
			sh.assign_subset(*it, subset);
		for (std::set<Vertex*>::iterator it = boundaryVertices[owner].begin();
		     it != boundaryVertices[owner].end(); ++it)
			sh.assign_subset(*it, subset);
	}
	UG_LOGN("Refinement level " << refinementLevel << ": restored "
	        << restoredFaces << " active Basolateral child faces and their "
	        << "boundary edges/vertices.");
	return restoredFaces;
}


static void report_active_terminal_face_geometry
(
	MultiGrid& grid,
	ISubsetHandler& sh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	size_t nephronCount,
	size_t refinementLevel
)
{
	for (size_t n = 0; n < nephronCount; ++n)
	for (int terminal = 0; terminal < 2; ++terminal)
	{
		const int subset = (int)(6 * n + 4 + terminal);
		std::set<Vertex*> vertices;
		size_t allFaces = 0;
		size_t activeFaces = 0;
		number maxEdgeLength = 0.0;
		for (size_t level = 0; level <= grid.top_level(); ++level)
		for (FaceIterator it = grid.begin<Face>(level);
		     it != grid.end<Face>(level); ++it)
		{
			Face* face = *it;
			if (sh.get_subset_index(face) != subset) continue;
			++allFaces;
			if (grid.has_children(face)) continue;
			++activeFaces;
			for (size_t v = 0; v < face->num_vertices(); ++v)
				vertices.insert(face->vertex(v));
			EdgeDescriptor edge;
			for (size_t e = 0; e < face->num_edges(); ++e)
			{
				face->edge_desc(e, edge);
				maxEdgeLength = std::max(maxEdgeLength,
					VecDistance(aaPos[edge.vertex(0)], aaPos[edge.vertex(1)]));
			}
		}
		number bboxDiagonal = 0.0;
		if (!vertices.empty())
		{
			vector3 minPos = aaPos[*vertices.begin()];
			vector3 maxPos = minPos;
			for (std::set<Vertex*>::iterator it = vertices.begin();
			     it != vertices.end(); ++it)
				for (size_t d = 0; d < 3; ++d)
				{
					minPos[d] = std::min(minPos[d], aaPos[*it][d]);
					maxPos[d] = std::max(maxPos[d], aaPos[*it][d]);
				}
			bboxDiagonal = VecDistance(minPos, maxPos);
		}
		UG_LOGN("Refinement level " << refinementLevel << " nephron " << n + 1
		        << (terminal == 0 ? " Inlet" : " Outlet")
		        << " active-cap check: faces=" << activeFaces << "/" << allFaces
		        << ", vertices=" << vertices.size()
		        << ", max edge=" << maxEdgeLength
		        << ", bounding diagonal=" << bboxDiagonal << ".");
	}
}


static void remove_orphan_objects_after_projected_refinement
(
	MultiGrid& g,
	size_t refinementLevel
)
{
	// The completed nephron mesh contains two coincident descriptions of some
	// terminal interfaces: the structured tube cap and TetGen's conforming
	// triangular boundary.  They share the coarse vertices, but refinement can
	// create child sides for the description that is not used by any volume.
	// If those detached children survive, the next refinement passes their
	// edges to NeuriteProjector and can collapse nearby tetrahedra at axial 0/1.
	//
	// Erase strictly in descending dimension.  A face with no adjacent volume
	// is never part of the physical 3-D mesh.  Once such faces are gone, edges
	// and vertices with no remaining incidence are equally safe to remove.
	std::vector<Face*> orphanFaces;
	for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
		if (NumAssociatedVolumes(g, *it) == 0) orphanFaces.push_back(*it);
	for (size_t i = 0; i < orphanFaces.size(); ++i) g.erase(orphanFaces[i]);

	std::vector<Edge*> orphanEdges;
	for (EdgeIterator it = g.begin<Edge>(); it != g.end<Edge>(); ++it)
		if (NumAssociatedFaces(g, *it) == 0) orphanEdges.push_back(*it);
	for (size_t i = 0; i < orphanEdges.size(); ++i) g.erase(orphanEdges[i]);

	std::vector<Vertex*> orphanVertices;
	for (VertexIterator it = g.begin<Vertex>(); it != g.end<Vertex>(); ++it)
		if (NumAssociatedEdges(g, *it) == 0) orphanVertices.push_back(*it);
	for (size_t i = 0; i < orphanVertices.size(); ++i)
		g.erase(orphanVertices[i]);

	UG_LOGN("Projected refinement " << refinementLevel
	        << " orphan cleanup: removed " << orphanFaces.size()
	        << " faces, " << orphanEdges.size() << " edges, and "
	        << orphanVertices.size() << " vertices.");

	// This is a hard topology invariant for every filled nephron level.  It is
	// intentionally checked here so a bad level cannot become the parent of the
	// next projected refinement.
	for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
		UG_COND_THROW(NumAssociatedVolumes(g, *it) == 0,
		              "Projected refinement " << refinementLevel
		              << " retained an orphan face after cleanup.");
}


static void report_nephron_surface_projection_residual
(
	Grid& g,
	ISubsetHandler& sh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >&
		aaSurfParams,
	NeuriteProjector& projector,
	size_t numNephrons,
	size_t refinementLevel
)
{
	const number endpointTolerance = 1e-7;
	const number relativeResidualTolerance = 1e-8;
	size_t checked = 0;
	size_t limited = 0;
	number maximumResidual = 0.0;
	number residualSum = 0.0;
	number maximumPositionDisplacement = 0.0;
	number positionDisplacementSum = 0.0;
	MultiGrid* multigrid = dynamic_cast<MultiGrid*>(&g);
	for (VertexIterator it = g.begin<Vertex>(); it != g.end<Vertex>(); ++it)
	{
		Vertex* vertex = *it;
		if (multigrid && multigrid->get_level(vertex) != refinementLevel)
			continue;
		const int si = sh.get_subset_index(vertex);
		bool surfaceVertex = false;
		for (size_t n = 0; n < numNephrons; ++n)
			if (si == (int)(6*n + 2) || si == (int)(6*n + 3))
			{
				surfaceVertex = true;
				break;
			}
		if (!surfaceVertex) continue;
		const uint32 rawNeuriteID = aaSurfParams[vertex].neuriteID;
		const uint32 plainNeuriteID = (rawNeuriteID << 12) >> 12;
		if (plainNeuriteID >= projector.neurites().size()) continue;
		const number axial = aaSurfParams[vertex].axial;
		if (axial <= endpointTolerance || axial >= 1.0 - endpointTolerance)
			continue;

		const vector3 positionBefore = aaPos[vertex];
		const NeuriteProjector::SurfaceParams paramsBefore = aaSurfParams[vertex];
		projector.project(vertex);
		const vector3 exactPosition = aaPos[vertex];
		const number positionDisplacement =
			VecDistance(positionBefore, exactPosition);
		positionDisplacementSum += positionDisplacement;
		maximumPositionDisplacement = std::max(
			maximumPositionDisplacement, positionDisplacement);
		aaPos[vertex] = positionBefore;
		aaSurfParams[vertex] = paramsBefore;
		aaSurfParams[vertex].radial = 0.0;
		projector.project(vertex);
		const vector3 exactCenter = aaPos[vertex];
		const number requestedRadius = VecDistance(exactPosition, exactCenter);
		const number actualRadius = VecDistance(positionBefore, exactCenter);
		const number residual = std::fabs(actualRadius - requestedRadius);
		aaPos[vertex] = positionBefore;
		aaSurfParams[vertex] = paramsBefore;
		++checked;
		residualSum += residual;
		maximumResidual = std::max(maximumResidual, residual);
		const number residualTolerance = std::max(
			std::numeric_limits<number>::epsilon() * requestedRadius * 64.0,
			requestedRadius * relativeResidualTolerance);
		if (residual > residualTolerance) ++limited;
	}
	UG_LOGN("Refinement level " << refinementLevel
	        << " Apical/Basolateral radius audit: checked=" << checked
	        << ", off-radius=" << limited
	        << ", mean residual=" << (checked ? residualSum/checked : 0.0)
	        << ", max residual=" << maximumResidual
	        << ", mean full-position displacement="
	        << (checked ? positionDisplacementSum/checked : 0.0)
	        << ", max full-position displacement="
	        << maximumPositionDisplacement << ".");
}


static void import_neurites_from_swc_impl
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number anisotropy,
	size_t numRefs,
	number boxPadding,
	bool fillInter
)
{
	// read in file to intermediate structure
    std::string inFileName = FindFileInStandardPaths(fileNameIn.c_str());
    UG_COND_THROW(inFileName == "", "File '" << fileNameIn
    	<< "' could not be located in standard paths.");

	FileReaderSWC swcFileReader;
	swcFileReader.load_file(inFileName.c_str());
	std::vector<swc_types::SWCPoint>& vPoints = swcFileReader.swc_points();

	// preconditioning
	//smoothing(vPoints, 5, 1.0, 1.0);

	// convert intermediate structure to neurite data
	std::vector<std::vector<vector3> > vPos;
	std::vector<std::vector<number> > vRad;
	std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
	std::vector<size_t> vRootNeuriteIndsOut;

	convert_pointlist_to_neuritelist(vPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

	// prepare grid and projector
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	Selector sel(g);


	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
			"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	if (!g.has_vertex_attachment(aSP))
		g.attach_to_vertices(aSP);

	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	aaSurfParams.access(g, aSP);


	ProjectionHandler projHandler(&sh);
	SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
	projHandler.set_geometry(geom3d);

	SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
	projHandler.set_projector(0, neuriteProj);

	// create spline data
	std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
	create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

	// create coarse grid
	for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i)
		create_neurite(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
			anisotropy, g, aaPos, aaSurfParams, NULL, NULL);

	// at branching points, we have not computed the correct positions yet,
	// so project the complete geometry using the projector
	VertexIterator vit = g.begin<Vertex>();
	VertexIterator vit_end = g.end<Vertex>();
	for (; vit != vit_end; ++vit)
		neuriteProj->project(*vit);

	if (boxPadding > 0.0)
		create_padded_box_surface(g, sh, aaPos, boxPadding, 1);
	if (fillInter)
		tetrahedralize_lumen_and_inter(g, sh);

	// assign subset
	AssignSubsetColors(sh);
	sh.set_subset_name(fillInter ? "Apical" : "neurites", 0);
	if (boxPadding > 0.0)
		sh.set_subset_name("OuterWall", 1);

	// output
	std::string outFileNameBase = FilenameAndPathWithoutExtension(fileNameOut);
	std::string outFileName = outFileNameBase + ".ugx";
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	ugxWriter.add_projection_handler(projHandler, "defPH", 0);
	if (!ugxWriter.write_to_file(outFileName.c_str()))
		UG_THROW("Grid could not be written to file '" << outFileName << "'.");

	if (numRefs == 0)
		return;

	// refinement
	Domain3d dom;
	dom.create_additional_subset_handler("projSH");
	try {LoadDomain(dom, outFileName.c_str());}
	UG_CATCH_THROW("Failed loading domain from '" << outFileName << "'.");
	number offset = 5.0;
	std::string curFileName = outFileName.substr(0, outFileName.size()-4) + "_refined_0.ugx";
	try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
	UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

	GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
	Grid::VertexAttachmentAccessor<APosition> refinedPos(*dom.grid(), aPosition);
	for (uint i = 0; i < numRefs; ++i)
	{
		ref.refine();
		if (fillInter)
			remove_orphan_objects_after_projected_refinement(
				*dom.grid(), i + 1);
		validate_and_repair_orientation_after_projected_refinement(
			*dom.grid(), refinedPos, *dom.subset_handler(), i + 1);

		std::ostringstream oss;
		oss << "_refined_" << i+1 << ".ugx";
		curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
		try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
		UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
	}
}


void import_neurites_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number anisotropy,
	size_t numRefs
)
{
	import_neurites_from_swc_impl(fileNameIn, fileNameOut, anisotropy, numRefs, -1.0, false);
}


void import_neurites_with_box_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number anisotropy,
	size_t numRefs,
	number padding
)
{
	import_neurites_from_swc_impl(fileNameIn, fileNameOut, anisotropy, numRefs, padding, false);
}


void import_neurites_with_inter_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number anisotropy,
	size_t numRefs,
	number padding
)
{
	import_neurites_from_swc_impl(fileNameIn, fileNameOut, anisotropy, numRefs, padding, true);
}



static size_t assign_indexed_nephron_subsets
(
	Grid& g,
	SubsetHandler& sh,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	Grid::VertexAttachmentAccessor<
		Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
	size_t expectedNephrons
)
{
	// Remove every orphan face before component indexing. A completed 3-D mesh
	// must associate each physical boundary face with one volume and every
	// internal interface face with two. Faces with no volume are legacy surface
	// constraints left behind by reconstruction; projecting/refining them creates
	// exactly the detached triangular shards seen near bends and terminal caps.
	std::vector<Face*> detachedBoundaryFaces;
	for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
	{
		Face* face = *it;
		// Any face with no adjacent volume is a detached reconstruction artifact,
		// irrespective of its old subset. In particular, preserving old subset-1
		// terminal fragments creates the spikes seen outside the nephron.
		if (NumAssociatedVolumes(g, face) != 0) continue;
		detachedBoundaryFaces.push_back(face);
	}
	for (size_t i = 0; i < detachedBoundaryFaces.size(); ++i)
		g.erase(detachedBoundaryFaces[i]);

	// Sweep the complete grid, not only the closure of the faces above. Earlier
	// reconstruction stages can leave standalone radial edges that were never
	// attached to a face; those appear as the star-shaped spikes at either cap.
	std::vector<Edge*> detachedEdges;
	for (EdgeIterator it = g.begin<Edge>(); it != g.end<Edge>(); ++it)
		if (NumAssociatedFaces(g, *it) == 0) detachedEdges.push_back(*it);
	for (size_t i = 0; i < detachedEdges.size(); ++i) g.erase(detachedEdges[i]);

	std::vector<Vertex*> detachedVertices;
	for (VertexIterator it = g.begin<Vertex>(); it != g.end<Vertex>(); ++it)
		if (NumAssociatedEdges(g, *it) == 0) detachedVertices.push_back(*it);
	for (size_t i = 0; i < detachedVertices.size(); ++i) g.erase(detachedVertices[i]);
	if (!detachedBoundaryFaces.empty())
		UG_LOGN("Removed " << detachedBoundaryFaces.size()
		        << " orphan faces, " << detachedEdges.size()
		        << " orphan edges, and " << detachedVertices.size()
		        << " orphan vertices before subset indexing.");

	/*
	 * MULTI-NEPHRON SUBSET INDEXING
	 * -----------------------------
	 * The merged SWC contains several disconnected root trees. The common tube
	 * generator initially places all trees into the same legacy subsets:
	 *   0 Membrane, 1 Lumen, 2 Basolateral, 3 Apical, 4 OuterWall, 5 Inter.
	 *
	 * This post-processing pass separates disconnected Membrane and Lumen volume
	 * components and assigns six stable subset slots per input nephron:
	 *   6*i+0 Membrane_i
	 *   6*i+1 Lumen_i
	 *   6*i+2 Basolateral_i
	 *   6*i+3 Apical_i
	 *   6*i+4 Inlet_i
	 *   6*i+5 Outlet_i
	 * followed by one shared OuterWall, one shared Inter, one
	 * TripleJunction_i subset per nephron, and one MembraneWall_i subset per
	 * nephron for the planar membrane annulus touching the OuterWall.
	 *
	 * Lumen boundary faces that are not Apical are terminal caps. Their mean
	 * npSurfParams.axial value distinguishes the root/low-axial Inlet from the
	 * terminal/high-axial Outlet. Cap faces and their interior vertices receive
	 * the cap subset. Edges and vertices on the shared outer circular rim receive
	 * a dedicated TripleJunction subset. The same NeuriteProjector is attached to
	 * that subset, keeping refinement conforming without classifying the
	 * Inter-touching rim as part of the Apical transport interface.
	 *
	 * The order of disconnected components follows the order in which root trees
	 * were generated from the merged SWC, which is the -swc1, -swc2, ... order.
	 */
	typedef std::map<Volume*, size_t> ComponentMap;

	auto findComponents = [&] (int sourceSubset, ComponentMap& componentOf)
	{
		size_t component = 0;
		for (VolumeIterator vit = g.begin<Volume>(); vit != g.end<Volume>(); ++vit)
		{
			Volume* seed = *vit;
			if (sh.get_subset_index(seed) != sourceSubset ||
			    componentOf.count(seed)) continue;

			std::vector<Volume*> stack(1, seed);
			while (!stack.empty())
			{
				Volume* vol = stack.back();
				stack.pop_back();
				if (componentOf.count(vol) ||
				    sh.get_subset_index(vol) != sourceSubset) continue;
				componentOf[vol] = component;

				Grid::traits<Face>::secure_container faces;
				g.associated_elements(faces, vol);
				for (size_t i = 0; i < faces.size(); ++i)
				{
					Volume* neighbor = GetConnectedNeighbor(g, faces[i], vol);
					if (neighbor && !componentOf.count(neighbor) &&
					    sh.get_subset_index(neighbor) == sourceSubset)
						stack.push_back(neighbor);
				}
			}
			++component;
		}
		return component;
	};

	ComponentMap membraneComponent;
	ComponentMap lumenComponent;
	const size_t numMembranes = findComponents(0, membraneComponent);
	const size_t numLumens = findComponents(1, lumenComponent);
	UG_COND_THROW(numMembranes != expectedNephrons ||
	              numLumens != expectedNephrons,
	              "Expected " << expectedNephrons
	              << " disconnected nephron components, but found "
	              << numMembranes << " Membrane and " << numLumens
	              << " Lumen components.");

	const int outerWallSubset = (int)(6 * expectedNephrons);
	const int interSubset = outerWallSubset + 1;
	const int tripleJunctionBase = interSubset + 1;
	const int membraneWallBase = tripleJunctionBase + (int)expectedNephrons;

	std::vector<Volume*> originalInterVolumes;
	for (VolumeIterator it = g.begin<Volume>(); it != g.end<Volume>(); ++it)
		if (sh.get_subset_index(*it) == 5)
			originalInterVolumes.push_back(*it);

	std::map<Face*, int> oldFaceSubset;
	std::map<Edge*, int> oldEdgeSubset;
	std::map<Vertex*, int> oldVertexSubset;
	for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
		oldFaceSubset[*it] = sh.get_subset_index(*it);
	for (EdgeIterator it = g.begin<Edge>(); it != g.end<Edge>(); ++it)
		oldEdgeSubset[*it] = sh.get_subset_index(*it);
	for (VertexIterator it = g.begin<Vertex>(); it != g.end<Vertex>(); ++it)
		oldVertexSubset[*it] = sh.get_subset_index(*it);

	// The padded box encloses every other object, so its six planes are the global
	// coordinate bounds of the completed grid.  Use those geometric planes as
	// the physical OuterWall identity.  The original subset-4 constraint faces
	// may have been replaced by TetGen triangles during duplicate cleanup, so
	// their old subset labels are not a reliable identifier at this point.
	// Topology alone is not sufficient either: a hole accidentally opened next
	// to the nephron is also a one-sided Inter boundary.
	vector3 boxLo(std::numeric_limits<number>::max());
	vector3 boxHi(-std::numeric_limits<number>::max());
	bool foundGridVertex = false;
	for (VertexIterator it = g.begin<Vertex>(); it != g.end<Vertex>(); ++it)
	{
		const vector3& p = aaPos[*it];
		for (size_t d = 0; d < 3; ++d)
		{
			boxLo[d] = std::min(boxLo[d], p[d]);
			boxHi[d] = std::max(boxHi[d], p[d]);
		}
		foundGridVertex = true;
	}
	UG_COND_THROW(!foundGridVertex, "Cannot classify an empty nephron grid.");
	number boxSpan = 0.0;
	for (size_t d = 0; d < 3; ++d)
		boxSpan = std::max(boxSpan, boxHi[d] - boxLo[d]);
	const number boxPlaneTolerance = std::max<number>(
		std::numeric_limits<number>::epsilon()
			* std::max<number>(boxSpan, 1.0) * 64.0,
		1e-8 * boxSpan);
	auto liesOnOuterBoxPlane = [&] (Face* face) -> bool
	{
		for (size_t d = 0; d < 3; ++d)
		{
			bool allAtLow = true;
			bool allAtHigh = true;
			for (size_t i = 0; i < face->num_vertices(); ++i)
			{
				const number x = aaPos[face->vertex(i)][d];
				allAtLow = allAtLow && fabs(x - boxLo[d]) <= boxPlaneTolerance;
				allAtHigh = allAtHigh && fabs(x - boxHi[d]) <= boxPlaneTolerance;
			}
			if (allAtLow || allAtHigh) return true;
		}
		return false;
	};

	for (ComponentMap::const_iterator it = membraneComponent.begin();
	     it != membraneComponent.end(); ++it)
		sh.assign_subset(it->first, (int)(6 * it->second));
	for (ComponentMap::const_iterator it = lumenComponent.begin();
	     it != lumenComponent.end(); ++it)
		sh.assign_subset(it->first, (int)(6 * it->second + 1));
	for (size_t i = 0; i < originalInterVolumes.size(); ++i)
		sh.assign_subset(originalInterVolumes[i], interSubset);

	auto adjacentMembrane = [&] (Face* face) -> int
	{
		Grid::traits<Volume>::secure_container vols;
		g.associated_elements(vols, face);
		for (size_t i = 0; i < vols.size(); ++i)
		{
			ComponentMap::const_iterator found = membraneComponent.find(vols[i]);
			if (found != membraneComponent.end()) return (int)found->second;
		}
		return -1;
	};

	auto adjacentLumen = [&] (Face* face) -> int
	{
		Grid::traits<Volume>::secure_container vols;
		g.associated_elements(vols, face);
		for (size_t i = 0; i < vols.size(); ++i)
		{
			ComponentMap::const_iterator found = lumenComponent.find(vols[i]);
			if (found != lumenComponent.end()) return (int)found->second;
		}
		return -1;
	};

	// TetGen triangulates each quadrilateral Basolateral constraint.  Those
	// Inter-side triangles share the membrane's boundary vertices but not its
	// quadrilateral face object, so they have only one face-adjacent Inter volume.
	// Recover their nephron owner from membrane volumes incident to the vertices.
	// -1 means no membrane owner; -2 means an invalid connection to two nephrons.
	auto vertexAdjacentMembrane = [&] (Face* face) -> int
	{
		int component = -1;
		for (size_t i = 0; i < face->num_vertices(); ++i)
		{
			Grid::traits<Volume>::secure_container vols;
			g.associated_elements(vols, face->vertex(i));
			for (size_t j = 0; j < vols.size(); ++j)
			{
				ComponentMap::const_iterator found = membraneComponent.find(vols[j]);
				if (found == membraneComponent.end()) continue;
				const int foundComponent = (int)found->second;
				if (component >= 0 && component != foundComponent) return -2;
				component = foundComponent;
			}
		}
		return component;
	};

	// Return a Lumen component only when every face vertex is incident to that
	// same Lumen.  This identifies TetGen's triangular counterpart of a terminal
	// Lumen-cap quadrilateral.  Using only one shared apical-ring vertex would be
	// too permissive and would incorrectly include the membrane terminal annulus.
	auto allVerticesAdjacentLumen = [&] (Face* face) -> int
	{
		int component = -1;
		for (size_t i = 0; i < face->num_vertices(); ++i)
		{
			int vertexComponent = -1;
			Grid::traits<Volume>::secure_container vols;
			g.associated_elements(vols, face->vertex(i));
			for (size_t j = 0; j < vols.size(); ++j)
			{
				ComponentMap::const_iterator found = lumenComponent.find(vols[j]);
				if (found == lumenComponent.end()) continue;
				const int foundComponent = (int)found->second;
				if (vertexComponent >= 0 && vertexComponent != foundComponent)
					return -2;
				vertexComponent = foundComponent;
			}
			if (vertexComponent < 0) return -1;
			if (component >= 0 && component != vertexComponent) return -2;
			component = vertexComponent;
		}
		return component;
	};

	// Return the owning Lumen component only for a true terminal interface:
	// exactly one adjacent Lumen, no adjacent Membrane, and any other adjacent
	// volume must be Inter. A cap may therefore have either one volume (open
	// exterior side) or two volumes (Lumen--Inter after the box is filled).
	auto terminalLumen = [&] (Face* face) -> int
	{
		Grid::traits<Volume>::secure_container vols;
		g.associated_elements(vols, face);
		int component = -1;
		size_t lumenCount = 0;
		for (size_t i = 0; i < vols.size(); ++i)
		{
			ComponentMap::const_iterator lumen = lumenComponent.find(vols[i]);
			if (lumen != lumenComponent.end())
			{
				component = (int)lumen->second;
				++lumenCount;
				continue;
			}
			if (membraneComponent.find(vols[i]) != membraneComponent.end())
				return -1;
			if (sh.get_subset_index(vols[i]) != interSubset)
				return -1;
		}
		return lumenCount == 1 ? component : -1;
	};

	std::vector<number> minCapAxial(expectedNephrons,
		std::numeric_limits<number>::max());
	std::vector<number> maxCapAxial(expectedNephrons,
		-std::numeric_limits<number>::max());
	std::map<Face*, number> capAxial;
	std::map<Face*, int> interSideCapOwner;
	// Automatically classify only the two true exterior terminal boundaries.
	// Internal Lumen faces are explicitly excluded below.
	const bool assignTerminalCapSubsets = true;
	size_t invalidOpenInterFaces = 0;
	for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
	{
		Face* face = *it;
		if (!assignTerminalCapSubsets) continue;
		// The extruder does not reliably preserve legacy subset 1 on the far-end
		// cap, so neighboring volume domains are authoritative. Internal Lumen
		// faces have two Lumen neighbors and are explicitly excluded.
		int component = terminalLumen(face);
		if (component < 0) continue;
		number axial = 0.0;
		for (size_t i = 0; i < face->num_vertices(); ++i)
			axial += aaSurfParams[face->vertex(i)].axial;
		axial /= face->num_vertices();
		capAxial[face] = axial;
		minCapAxial[component] = std::min(minCapAxial[component], axial);
		maxCapAxial[component] = std::max(maxCapAxial[component], axial);
	}
	for (size_t component = 0; component < expectedNephrons; ++component)
	{
		if (!assignTerminalCapSubsets) continue;
		UG_COND_THROW(minCapAxial[component] == std::numeric_limits<number>::max(),
		              "Could not identify lumen terminal faces for nephron "
		              << component + 1 << ".");
		UG_LOGN("Lumen terminal axial range for nephron " << component + 1
		        << ": " << minCapAxial[component] << " to "
		        << maxCapAxial[component] << ".");
	}

	// A quadrilateral Lumen-cap face and TetGen's two Inter-side triangles occupy
	// the same physical patch but are distinct grid faces.  Add the triangle
	// counterparts to the same Inlet/Outlet set.  Requiring all vertices to be
	// Lumen-adjacent excludes the surrounding membrane terminal annulus.
	for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
	{
		Face* face = *it;
		if (!assignTerminalCapSubsets || capAxial.count(face)) continue;
		Grid::traits<Volume>::secure_container vols;
		g.associated_elements(vols, face);
		if (vols.size() != 1 || sh.get_subset_index(vols[0]) != interSubset)
			continue;
		const int component = allVerticesAdjacentLumen(face);
		if (component < 0) continue;
		number axial = 0.0;
		for (size_t i = 0; i < face->num_vertices(); ++i)
			axial += aaSurfParams[face->vertex(i)].axial;
		axial /= face->num_vertices();
		const number axialTolerance = 1e-6 * std::max<number>(
			1.0, fabs(maxCapAxial[component] - minCapAxial[component]));
		if (fabs(axial - minCapAxial[component]) > axialTolerance &&
		    fabs(axial - maxCapAxial[component]) > axialTolerance)
			continue;
		capAxial[face] = axial;
		interSideCapOwner[face] = component;
	}

	auto capOwner = [&] (Face* face) -> int
	{
		const int directOwner = terminalLumen(face);
		if (directOwner >= 0) return directOwner;
		std::map<Face*, int>::const_iterator found = interSideCapOwner.find(face);
		return found == interSideCapOwner.end() ? -1 : found->second;
	};

	for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
	{
		Face* face = *it;
		const int oldSI = oldFaceSubset[face];
		int newSI = -1;
		std::map<Face*, number>::const_iterator cap = capAxial.find(face);
		if (cap != capAxial.end())
		{
			int component = capOwner(face);
			UG_COND_THROW(component < 0,
			              "Could not associate terminal cap with a Lumen volume.");
			const number split = 0.5 *
				(minCapAxial[component] + maxCapAxial[component]);
			newSI = 6 * component + (cap->second <= split ? 4 : 5);
		}
		else
		{
			const int membrane = adjacentMembrane(face);
			const int lumen = adjacentLumen(face);
			Grid::traits<Volume>::secure_container vols;
			g.associated_elements(vols, face);
			bool hasInter = false;
			for (size_t i = 0; i < vols.size(); ++i)
				hasInter = hasInter || sh.get_subset_index(vols[i]) == interSubset;
			if (membrane >= 0 && lumen < 0 && liesOnOuterBoxPlane(face))
				newSI = membraneWallBase + membrane; // terminal Membrane--OuterWall contact
			else if (membrane >= 0 && lumen >= 0)
				newSI = 6 * membrane + 3; // Lumen--Membrane: Apical
			else if (membrane >= 0 && hasInter)
				newSI = 6 * membrane + 2; // outer Membrane: Basolateral
			else if (oldSI == 2 && !liesOnOuterBoxPlane(face))
			{
				// A structured Basolateral constraining quad is adjacent only to
				// Membrane, while each constrained triangular counterpart is
				// adjacent only to Inter. Their constraint relation—not ordinary
				// face adjacency—makes them one physical interface.
				const int owner = membrane >= 0 ? membrane
					: vertexAdjacentMembrane(face);
				if (owner >= 0) newSI = 6 * owner + 2;
			}
			else if (hasInter && membrane < 0 && lumen < 0 && vols.size() == 1)
			{
				if (liesOnOuterBoxPlane(face))
					newSI = outerWallSubset;
				else
				{
					const int owner = vertexAdjacentMembrane(face);
					if (owner >= 0)
					{
						// The conforming cylindrical interface has already acquired
						// both adjacent volumes. A remaining one-sided Inter face
						// whose vertices touch Membrane is part of a terminal annulus,
						// not a Basolateral transport interface.
						newSI = interSubset;
					}
					else
					{
						// Keep the assignment harmless until the complete mesh has
						// been inspected, then fail with the aggregate count below.
						newSI = interSubset;
						++invalidOpenInterFaces;
					}
				}
			}
			else if (oldSI == 4 && liesOnOuterBoxPlane(face))
				newSI = outerWallSubset;
			else if (oldSI == 1)
			{
				// Retain unassigned terminal faces in the owning Lumen subset.
				int component = lumen;
				if (component < 0 && expectedNephrons == 1) component = 0;
				if (component >= 0) newSI = 6 * component + 1;
			}
			else if (!vols.empty())
			{
				newSI = sh.get_subset_index(vols[0]);
			}
		}
		if (newSI >= 0) sh.assign_subset(face, newSI);
	}
	UG_COND_THROW(invalidOpenInterFaces != 0,
	              "Detected " << invalidOpenInterFaces
	              << " open Inter boundary faces away from the six OuterWall "
	              << "box planes. The exterior tetrahedralization contains "
	              << "holes next to a nephron.");

	auto indexedBoundarySubset = [&] (int si) -> bool
	{
		if (si == outerWallSubset) return true;
		if (si < 0 || si >= outerWallSubset) return false;
		const int role = si % 6;
		return role >= 2 && role <= 5;
	};

	for (EdgeIterator it = g.begin<Edge>(); it != g.end<Edge>(); ++it)
	{
		Edge* edge = *it;
		int fallbackSI = -1;
		int boundarySI = -1;
		int capSI = -1;
		int apicalSI = -1;
		int basolateralSI = -1;
		int membraneWallSI = -1;
		Grid::traits<Face>::secure_container faces;
		g.associated_elements(faces, edge);
		for (size_t i = 0; i < faces.size(); ++i)
		{
			const int faceSI = sh.get_subset_index(faces[i]);
			const int role = faceSI >= 0 && faceSI < outerWallSubset
				? faceSI % 6 : -1;
			if (faceSI >= membraneWallBase &&
			    faceSI < membraneWallBase + (int)expectedNephrons)
			{
				membraneWallSI = faceSI;
				continue;
			}
			// A terminal cap shares its outer rim with Apical. Keep only that
			// shared rim edge Apical so NeuriteProjector propagates npSurfParams
			// to its refined midpoint. Every other cap edge must be Inlet/Outlet;
			// otherwise its midpoint inherits Lumen, Basolateral, or OuterWall
			// and a Lumen element contains a vertex without flow DoFs.
			if (role == 3)
			{
				apicalSI = faceSI;
				continue;
			}
			if (role == 2)
			{
				UG_COND_THROW(basolateralSI >= 0 && basolateralSI != faceSI,
				              "Edge belongs to Basolateral faces of different nephrons.");
				basolateralSI = faceSI;
				continue;
			}
			if (role == 4 || role == 5)
			{
				UG_COND_THROW(capSI >= 0 && capSI != faceSI,
				              "Edge belongs to different terminal cap subsets.");
				capSI = faceSI;
				continue;
			}
			if (indexedBoundarySubset(faceSI))
			{
				if (boundarySI < 0) boundarySI = faceSI;
				continue;
			}
			if (fallbackSI < 0) fallbackSI = faceSI;
		}
		// An edge shared by an Apical face and an Inlet/Outlet cap is the
		// three-domain endpoint rim. Keep it separate from the physical Apical
		// interface while retaining a deterministic projector-aware subset.
		int tripleJunctionSI = -1;
		if (apicalSI >= 0 && capSI >= 0)
		{
			UG_COND_THROW(apicalSI / 6 != capSI / 6,
			              "Terminal rim joins subsets from different nephrons.");
			tripleJunctionSI = tripleJunctionBase + apicalSI / 6;
		}

		// Do not let iteration order decide the edge subset. In particular,
		// Basolateral/OuterWall used to win before a later cap face was seen.
		const int newSI = tripleJunctionSI >= 0 ? tripleJunctionSI :
		                  apicalSI >= 0 ? apicalSI :
		                  capSI >= 0 ? capSI :
		                  membraneWallSI >= 0 ? membraneWallSI :
		                  basolateralSI >= 0 ? basolateralSI :
		                  boundarySI >= 0 ? boundarySI : fallbackSI;
		if (newSI >= 0) sh.assign_subset(edge, newSI);
	}

	for (VertexIterator it = g.begin<Vertex>(); it != g.end<Vertex>(); ++it)
	{
		Vertex* vertex = *it;
		int fallbackSI = -1;
		int boundarySI = -1;
		int capSI = -1;
		int apicalSI = -1;
		int basolateralSI = -1;
		int membraneWallSI = -1;
		int tripleJunctionSI = -1;
		Grid::traits<Edge>::secure_container edges;
		g.associated_elements(edges, vertex);
		for (size_t i = 0; i < edges.size(); ++i)
		{
			const int edgeSI = sh.get_subset_index(edges[i]);
			if (edgeSI >= membraneWallBase &&
			    edgeSI < membraneWallBase + (int)expectedNephrons)
			{
				membraneWallSI = edgeSI;
				continue;
			}
			if (edgeSI >= tripleJunctionBase &&
			    edgeSI < tripleJunctionBase + (int)expectedNephrons)
			{
				UG_COND_THROW(tripleJunctionSI >= 0 &&
				              tripleJunctionSI != edgeSI,
				              "Vertex belongs to TripleJunction rims of different nephrons.");
				tripleJunctionSI = edgeSI;
				continue;
			}
			const int role = edgeSI >= 0 && edgeSI < outerWallSubset
				? edgeSI % 6 : -1;
			// The outer circular boundary of a terminal cap is already identified
			// by its Apical edge.  Give that edge's end vertices the same Apical
			// ownership; otherwise a later cap edge encountered in iteration order
			// incorrectly turns the wall/cap rim into Inlet or Outlet vertices.
			if (role == 3)
			{
				UG_COND_THROW(apicalSI >= 0 && apicalSI != edgeSI,
				              "Vertex belongs to Apical rims of different nephrons.");
				apicalSI = edgeSI;
				continue;
			}
			if (role == 2)
			{
				UG_COND_THROW(basolateralSI >= 0 && basolateralSI != edgeSI,
				              "Vertex belongs to Basolateral rims of different nephrons.");
				basolateralSI = edgeSI;
				continue;
			}
			if (role == 4 || role == 5)
			{
				UG_COND_THROW(capSI >= 0 && capSI != edgeSI,
				              "Vertex belongs to different terminal cap subsets.");
				capSI = edgeSI;
				continue;
			}
			if (indexedBoundarySubset(edgeSI))
			{
				if (boundarySI < 0) boundarySI = edgeSI;
				continue;
			}
			if (fallbackSI < 0) fallbackSI = edgeSI;
		}
		const int newSI = tripleJunctionSI >= 0 ? tripleJunctionSI :
		                  apicalSI >= 0 ? apicalSI :
		                  capSI >= 0 ? capSI :
		                  membraneWallSI >= 0 ? membraneWallSI :
		                  basolateralSI >= 0 ? basolateralSI :
		                  boundarySI >= 0 ? boundarySI : fallbackSI;
		if (newSI >= 0) sh.assign_subset(vertex, newSI);
	}

	// The outer circular edge of the terminal membrane annulus is shared by a
	// structured Basolateral constraint and the triangulated OuterWall.  The
	// constrained-face relation can hide the Basolateral face from the ordinary
	// edge-to-face adjacency query above, causing OuterWall to win.  Identify this
	// rim from its exact neurite parameters and membrane-volume ownership.  Only
	// the rim edge/vertices become Basolateral; the surrounding planar faces stay
	// OuterWall.  This lets their new angular midpoints use NeuriteProjector.
	for (EdgeIterator it = g.begin<Edge>(); it != g.end<Edge>(); ++it)
	{
		Edge* edge = *it;
		const NeuriteProjector::SurfaceParams& p0 = aaSurfParams[edge->vertex(0)];
		const NeuriteProjector::SurfaceParams& p1 = aaSurfParams[edge->vertex(1)];
		const bool outerRadius = p0.radial > 1.0 - 1e-6
			&& p1.radial > 1.0 - 1e-6;
		const bool oneTerminal =
			(std::fabs((number)p0.axial) < 1e-7 &&
			 std::fabs((number)p1.axial) < 1e-7) ||
			(std::fabs((number)p0.axial - 1.0) < 1e-7 &&
			 std::fabs((number)p1.axial - 1.0) < 1e-7);
		if (!outerRadius || !oneTerminal) continue;
		int owner = -1;
		Grid::traits<Volume>::secure_container volumes;
		g.associated_elements(volumes, edge);
		for (size_t i = 0; i < volumes.size(); ++i)
		{
			ComponentMap::const_iterator found = membraneComponent.find(volumes[i]);
			if (found == membraneComponent.end()) continue;
			owner = (int)found->second;
			break;
		}
		if (owner < 0) continue;
		const int basolateralSubset = 6 * owner + 2;
		sh.assign_subset(edge, basolateralSubset);
		sh.assign_subset(edge->vertex(0), basolateralSubset);
		sh.assign_subset(edge->vertex(1), basolateralSubset);
	}

	// Validate the cap closure before serialization. A cap edge may remain in
	// TripleJunction only when it is the shared outer rim; every other cap edge
	// must use the same Inlet/Outlet subset as its cap face. A cap vertex is in
	// TripleJunction exactly when it touches one of those rim edges; every other
	// cap vertex is assigned to its Inlet/Outlet face. These rules also give all
	// refined midpoint vertices a deterministic boundary subset.
	std::vector<std::set<Edge*> > capClosureEdges(2 * expectedNephrons);
	std::vector<std::set<Edge*> > capTripleJunctionEdges(2 * expectedNephrons);
	std::vector<std::set<Vertex*> > capTripleJunctionVertices(2 * expectedNephrons);
	for (std::map<Face*, number>::const_iterator it = capAxial.begin();
	     it != capAxial.end(); ++it)
	{
		Face* face = it->first;
		const int capSubset = sh.get_subset_index(face);
		const int role = capSubset % 6;
		const int component = capSubset / 6;
		UG_COND_THROW(role != 4 && role != 5,
		              "Terminal face is not assigned to Inlet/Outlet.");
		const size_t capIndex = 2 * (size_t)component + (size_t)(role - 4);
		Grid::traits<Edge>::secure_container edges;
		g.associated_elements(edges, face);
		for (size_t i = 0; i < edges.size(); ++i)
		{
			const int edgeSI = sh.get_subset_index(edges[i]);
			const int tripleJunctionSI = tripleJunctionBase + component;
			UG_COND_THROW(edgeSI != capSubset && edgeSI != tripleJunctionSI,
			              "Terminal cap edge escaped Inlet/Outlet/TripleJunction closure.");
			if (edgeSI == capSubset) capClosureEdges[capIndex].insert(edges[i]);
			else capTripleJunctionEdges[capIndex].insert(edges[i]);
		}
		for (size_t i = 0; i < face->num_vertices(); ++i)
		{
			Vertex* vertex = face->vertex(i);
			bool onTripleJunction = false;
			Grid::traits<Edge>::secure_container vertexEdges;
			g.associated_elements(vertexEdges, vertex);
			for (size_t j = 0; j < vertexEdges.size(); ++j)
			{
				if (sh.get_subset_index(vertexEdges[j]) ==
				    tripleJunctionBase + component)
				{
					onTripleJunction = true;
					break;
				}
			}
			const int expectedVertexSubset =
				onTripleJunction ? tripleJunctionBase + component : capSubset;
			UG_COND_THROW(sh.get_subset_index(vertex) != expectedVertexSubset,
			              "Terminal cap vertex has the wrong Inlet/Outlet/TripleJunction subset.");
			if (onTripleJunction)
				capTripleJunctionVertices[capIndex].insert(vertex);
		}
	}
	for (size_t component = 0; component < expectedNephrons; ++component)
	{
		UG_LOGN("Nephron " << component + 1
		        << " cap edge closure: Inlet="
		        << capClosureEdges[2 * component].size()
		        << " plus " << capTripleJunctionEdges[2 * component].size()
		        << " TripleJunction edges / "
		        << capTripleJunctionVertices[2 * component].size()
		        << " TripleJunction vertices, Outlet="
		        << capClosureEdges[2 * component + 1].size()
		        << " plus " << capTripleJunctionEdges[2 * component + 1].size()
		        << " TripleJunction edges / "
		        << capTripleJunctionVertices[2 * component + 1].size()
		        << " TripleJunction vertices.");
	}

	// When manual cap assignment is requested, guarantee that the Inlet/Outlet
	// subsets are completely empty in every grid dimension. Older closure rules
	// may otherwise leave terminal rim vertices or edges in roles 4 and 5 even
	// though no cap faces were assigned.
	if (!assignTerminalCapSubsets)
	{
		for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
		{
			const int si = sh.get_subset_index(*it);
			if (si >= 0 && si < outerWallSubset && (si % 6 == 4 || si % 6 == 5))
				sh.assign_subset(*it, si - si % 6 + 1);
		}
		for (EdgeIterator it = g.begin<Edge>(); it != g.end<Edge>(); ++it)
		{
			const int si = sh.get_subset_index(*it);
			if (si >= 0 && si < outerWallSubset && (si % 6 == 4 || si % 6 == 5))
				sh.assign_subset(*it, si - si % 6 + 1);
		}
		for (VertexIterator it = g.begin<Vertex>(); it != g.end<Vertex>(); ++it)
		{
			const int si = sh.get_subset_index(*it);
			if (si >= 0 && si < outerWallSubset && (si % 6 == 4 || si % 6 == 5))
				sh.assign_subset(*it, si - si % 6 + 1);
		}
	}
	else
	{
		// Fail during generation if a cap face is not a true Lumen terminal
		// interface. This guards against the former bug where internal Lumen
		// faces inherited Inlet/Outlet merely because their legacy subset was 1.
		std::vector<size_t> inletFaces(expectedNephrons, 0);
		std::vector<size_t> outletFaces(expectedNephrons, 0);
		for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
		{
			Face* face = *it;
			const int si = sh.get_subset_index(face);
			if (si < 0 || si >= outerWallSubset) continue;
			const int role = si % 6;
			if (role != 4 && role != 5) continue;

			const size_t component = (size_t)(si / 6);
			UG_COND_THROW(capOwner(face) != (int)component,
			              "Inlet/Outlet face is not a terminal interface of its owning Lumen.");
			if (role == 4) ++inletFaces[component];
			else ++outletFaces[component];
		}
		for (size_t i = 0; i < expectedNephrons; ++i)
		{
			UG_COND_THROW(inletFaces[i] == 0 || outletFaces[i] == 0,
			              "Missing physical Inlet or Outlet faces for nephron "
			              << i + 1 << ".");
			UG_LOGN("Nephron " << i + 1 << " cap faces: Inlet="
			        << inletFaces[i] << ", Outlet=" << outletFaces[i] << ".");
		}
	}

	for (size_t i = 0; i < expectedNephrons; ++i)
	{
		std::ostringstream suffixStream;
		if (expectedNephrons > 1)
			suffixStream << "_" << i + 1;
		const std::string suffix = suffixStream.str();
		const std::string membraneName = "Membrane" + suffix;
		const std::string lumenName = "Lumen" + suffix;
		const std::string basolateralName = "Basolateral" + suffix;
		const std::string apicalName = "Apical" + suffix;
		const std::string inletName = "Inlet" + suffix;
		const std::string outletName = "Outlet" + suffix;
		sh.set_subset_name(membraneName.c_str(), (int)(6 * i));
		sh.set_subset_name(lumenName.c_str(), (int)(6 * i + 1));
		sh.set_subset_name(basolateralName.c_str(), (int)(6 * i + 2));
		sh.set_subset_name(apicalName.c_str(), (int)(6 * i + 3));
		sh.set_subset_name(inletName.c_str(), (int)(6 * i + 4));
		sh.set_subset_name(outletName.c_str(), (int)(6 * i + 5));
		const std::string tripleJunctionName = "TripleJunction" + suffix;
		sh.set_subset_name(tripleJunctionName.c_str(),
		                   tripleJunctionBase + (int)i);
		const std::string membraneWallName = "MembraneWall" + suffix;
		sh.set_subset_name(membraneWallName.c_str(),
		                   membraneWallBase + (int)i);
	}
	sh.set_subset_name("OuterWall", outerWallSubset);
	sh.set_subset_name("Inter", interSubset);
	return expectedNephrons;
}


static void import_er_neurites_from_swc_impl
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number erScaleFactor,
	number anisotropy,
	size_t numRefs,
	number boxPadding,
	bool fillInter,
	number interTetQuality,
	size_t indexedNephrons = 0,
	size_t nephronOGridVertices = 12,
	bool serializeProjector = true,
	size_t boxSurfaceRefs = 2,
	bool coarseLumenCenterOnly = false,
	bool crossSectionRefinement = false
)
{
	/*
	 * COMMON NESTED-TUBE IMPORT PIPELINE
	 * ----------------------------------
	 * Shared implementation used by the original ER/neurite generator, the
	 * single-nephron Lumen/Membrane/Inter importer, and the indexed multi-nephron
	 * importer.
	 *
	 * Major stages:
	 *   A. Parse SWC points and convert each root tree into spline-ready neurites.
	 *   B. Create a Grid, defSH/projSH handlers, npSurfParams attachment, and
	 *      NeuriteProjector containing the serialized centerline splines.
	 *   C. Generate nested anisotropic hexahedral tubes. The SWC radius is the
	 *      outer radius; erScaleFactor/lumenScaleFactor defines the inner radius.
	 *   D. Project initial tube vertices once to the analytical spline geometry.
	 *   E. When fillInter is enabled, create the padded box and call the isolated
	 *      TetGen routine above.
	 *   F. For multi-nephron input, split and rename disconnected components,
	 *      create Inlet_i/Outlet_i, and attach projectors only to
	 *      Basolateral_i/Apical_i.
	 *   G. Serialize defGrid, defSH, projSH, defPH, and npSurfParams into UGX.
	 *
	 * Projector scope is intentionally narrow. Projecting volume, cap, Inter, or
	 * OuterWall subsets onto a neurite surface causes collapsed/fanned refinement.
	 */
	// read in file to intermediate structure
	std::string inFileName = FindFileInStandardPaths(fileNameIn.c_str());
	UG_COND_THROW(inFileName == "", "File '" << fileNameIn
		<< "' could not be located in standard paths.");
	const OuterWallTerminalMetadata terminalMetadata =
		read_outer_wall_terminal_metadata(inFileName);
	UG_COND_THROW(terminalMetadata.terminalsOnOuterWall && !fillInter,
	              "Terminal-to-OuterWall SWC input requires the Inter/OuterWall importer.");

	FileReaderSWC swcFileReader;
	swcFileReader.load_file(inFileName.c_str());
	std::vector<swc_types::SWCPoint>& vPoints = swcFileReader.swc_points();

	// smoothing
	//smoothing(vPoints, 5, 1.0, 1.0);

	// convert intermediate structure to neurite data
	std::vector<std::vector<vector3> > vPos;
	std::vector<std::vector<number> > vRad;
	std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
	std::vector<size_t> vRootNeuriteIndsOut;

	convert_pointlist_to_neuritelist(vPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

	// prepare grid and projector
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);


	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
			"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	if (!g.has_vertex_attachment(aSP))
		g.attach_to_vertices(aSP);

	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	aaSurfParams.access(g, aSP);

	SubsetHandler psh(g);
	psh.set_default_subset_index(0);

	// Combined nephron/Inter grids must serialize the projector against defSH.
	// proMesh reloads defSH reliably, whereas a separate projSH can be lost as
	// the active projector subset handler during Mesh-based refinement.
	ProjectionHandler projHandler(fillInter ? &sh : &psh);
	SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
	projHandler.set_geometry(geom3d);

	SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
	if (fillInter && indexedNephrons == 0)
	{
		projHandler.set_projector(2, neuriteProj); // Basolateral
		projHandler.set_projector(3, neuriteProj); // Apical
	}
	else if (!fillInter)
		projHandler.set_projector(0, neuriteProj);

	// create spline data
	std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
	create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

	// create coarse grid
	for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i)
		create_neurite_with_er(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
			erScaleFactor, anisotropy, nephronOGridVertices, coarseLumenCenterOnly,
			g, aaPos, aaSurfParams, sh, NULL, NULL);

	// at branching points, we have not computed the correct positions yet,
	// so project the complete geometry using the projector
	VertexIterator vit = g.begin<Vertex>();
	VertexIterator vit_end = g.end<Vertex>();
	for (; vit != vit_end; ++vit)
		neuriteProj->project(*vit);

	// All original nested-tube elements belong to projector subset 0. Elements
	// created from this point onward (box and tetrahedral bulk) must not be
	// projected onto the neurite surface during later refinement.
	if (fillInter)
		psh.set_default_subset_index(1);
	if (boxPadding > 0.0)
		create_padded_box_surface(g, sh, aaPos, boxPadding, 4, boxSurfaceRefs,
		                          &terminalMetadata, &aaSurfParams);
	if (fillInter)
	{
		if (g_adaptiveNephronMeshing.prismaticInter)
			add_straight_prismatic_inter_to_version2_nephron(
				g, sh, psh, aaPos, aaSurfParams, terminalMetadata);
		else
			// Tetrahedralize Inter on an isolated copy so TetGen cannot delete or
			// move projector-backed tube vertices. Retain the established
			// constrained Basolateral representation for the selectable tet path;
			// the experimental pyramid/tetra transition is not refinement-safe.
			add_isolated_inter_to_version2_nephron(g, sh, psh, aaPos,
				aaSurfParams, interTetQuality, true, false);
	}
	if (fillInter)
	{
		// Keep the duplicate-vertex tolerance proportional to the input SWC
		// radius so it has the same physical meaning for um and mm meshes.
		number minOuterRadius = std::numeric_limits<number>::max();
		for (size_t i = 0; i < vRad.size(); ++i)
			for (size_t j = 0; j < vRad[i].size(); ++j)
				if (vRad[i][j] > 0.0)
					minOuterRadius = std::min(minOuterRadius, vRad[i][j]);
		UG_COND_THROW(minOuterRadius == std::numeric_limits<number>::max(),
		              "Nephron SWC has no positive radius for mesh cleanup.");
		cleanup_completed_nephron_mesh_before_final_subsets(g, aaPos,
			minOuterRadius * 1e-4);
	}

	if (indexedNephrons > 0)
	{
		assign_indexed_nephron_subsets(g, sh, aaPos, aaSurfParams, indexedNephrons);
		const int tripleJunctionBase = (int)(6 * indexedNephrons + 2);
		const int membraneWallBase = tripleJunctionBase + (int)indexedNephrons;
		for (size_t i = 0; i < indexedNephrons; ++i)
		{
			// Use the same exact neurite map for the internal radial layers as for
			// Apical and Basolateral. Moving only the two interfaces leaves their
			// supporting Lumen/Membrane vertices on straight chords.
			projHandler.set_projector((int)(6 * i + 0), neuriteProj);
			projHandler.set_projector((int)(6 * i + 1), neuriteProj);
			projHandler.set_projector((int)(6 * i + 2), neuriteProj);
			projHandler.set_projector((int)(6 * i + 3), neuriteProj);
			// Terminal cap refinement must also propagate npSurfParams.  Leaving
			// Inlet/Outlet on the default linear projector creates cap midpoint
			// vertices with zero axial/radial parameters.  A child edge shared with
			// Apical is then projected on the next level using the average of a
			// valid surface radius and zero, collapsing cells toward the lumen
			// center.  All vertices of each cap have one terminal axial coordinate,
			// so NeuriteProjector keeps the cap planar while interpolating its
			// radial and angular coordinates correctly.
			projHandler.set_projector((int)(6 * i + 4), neuriteProj);
			projHandler.set_projector((int)(6 * i + 5), neuriteProj);
			// The terminal rim is a three-domain junction, not an Apical
			// transport subset. It nevertheless needs the NeuriteProjector so
			// refined rim vertices inherit valid surface parameters.
			projHandler.set_projector(tripleJunctionBase + (int)i, neuriteProj);
			// The planar terminal membrane annulus uses the same parameter map as
			// the tube during refinement, then is snapped back to the OuterWall
			// plane. Keeping it separate permits a distinct wall-contact BC.
			projHandler.set_projector(membraneWallBase + (int)i, neuriteProj);
		}
	}

	// assign subset
	AssignSubsetColors(sh);
	if (!fillInter)
	{
		sh.set_subset_name("cyt", 0);
		sh.set_subset_name("er", 1);
		sh.set_subset_name("pm", 2);
		sh.set_subset_name("erm", 3);
	}

	// output
	std::string outFileNameBase = FilenameAndPathWithoutExtension(fileNameOut);
	std::string outFileName = outFileNameBase + ".ugx";
	// The no-projector pathway needs surface parameters only while positioning
	// the initial O-grid. Remove them before writing a clean linear-refinement
	// mesh that cannot invoke NeuriteProjector later.
	if (serializeProjector)
		report_nephron_surface_projection_residual(
			g, sh, aaPos, aaSurfParams, *neuriteProj,
			indexedNephrons, 0);
	if (!serializeProjector && g.has_vertex_attachment(aSP))
		g.detach_from_vertices(aSP);
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	if (serializeProjector)
	{
		ugxWriter.add_subset_handler(psh, "projSH", 0);
		ugxWriter.add_projection_handler(projHandler, "defPH", 0);
	}
	if (!ugxWriter.write_to_file(outFileName.c_str()))
		UG_THROW("Grid could not be written to file '" << outFileName << "'.");

	if (numRefs == 0)
		return;

	// refinement
	Domain3d dom;
	dom.create_additional_subset_handler("projSH");
	try {LoadDomain(dom, outFileName.c_str());}
	UG_CATCH_THROW("Failed loading domain from '" << outFileName << "'.");
	ProjectionHandler* loadedProjectionHandler =
		dynamic_cast<ProjectionHandler*>(dom.refinement_projector().get());
	UG_COND_THROW(!loadedProjectionHandler,
	              "Loaded nephron domain has no ProjectionHandler.");
	NeuriteProjector* loadedNeuriteProjector = dynamic_cast<NeuriteProjector*>(
		loadedProjectionHandler->projector(2).get());
	UG_COND_THROW(!loadedNeuriteProjector,
	              "Loaded Basolateral subset has no NeuriteProjector.");

	SmartPtr<IRefiner> ref;
	HangingNodeRefiner_MultiGrid* crossRef = NULL;
	if (crossSectionRefinement)
	{
		crossRef = new HangingNodeRefiner_MultiGrid(
			*dom.grid(), dom.refinement_projector());
		// Structured Lumen/Membrane hexahedra need face-center dependencies on
		// their two existing axial rings.  Those paired face centers are joined by
		// the new straight axial child edge; disabling order 1 forced the former
		// midpoint-to-coarse-corner triangular topology instead.
		crossRef->enable_node_dependency_order_1(true);
		ref = SmartPtr<IRefiner>(crossRef);
	}
	else
		ref = SmartPtr<IRefiner>(new GlobalMultiGridRefiner(
			*dom.grid(), dom.refinement_projector()));
	Grid::VertexAttachmentAccessor<APosition> refinedPos(*dom.grid(), aPosition);
	Grid::VertexAttachmentAccessor<Attachment<NPSP> > refinedSurfParams(
		*dom.grid(), aSP);
	for (uint i = 0; i < numRefs; ++i)
	{
		std::set<Edge*> transverseEdges;
		std::set<Edge*> axialEdges;
		std::map<uint32, std::vector<number> > allowedAxialSections;
		if (crossSectionRefinement)
		{
			crossRef->clear_marks();
			allowedAxialSections = collect_nephron_axial_sections(
				*dom.grid(), *dom.subset_handler(), refinedSurfParams,
				indexedNephrons > 0 ? indexedNephrons : 1);
			const size_t marked = mark_nephron_cross_section_volumes(
				*crossRef, *dom.grid(), *dom.subset_handler(), refinedPos, refinedSurfParams,
				indexedNephrons > 0 ? indexedNephrons : 1, true,
				transverseEdges, axialEdges);
			UG_LOGN("Cross-sectional refinement " << i + 1 << ": marked "
			        << marked << " unique transverse tube edges; axial edges "
			        << "remain unrefined.");
		}
		ref->refine();
		if (crossSectionRefinement)
			restore_active_basolateral_children(
				*dom.grid(), *dom.subset_handler(), refinedSurfParams,
				indexedNephrons > 0 ? indexedNephrons : 1, i + 1);
		if (crossSectionRefinement)
			restore_coincident_neurite_surface_params(
				*dom.grid(), *dom.subset_handler(), refinedPos, refinedSurfParams,
				indexedNephrons > 0 ? indexedNephrons : 1, i + 1);
		if (crossSectionRefinement)
		{
			const size_t newAxialVertices = count_vertices_on_new_axial_sections(
				*dom.grid(), *dom.subset_handler(), refinedSurfParams,
				indexedNephrons > 0 ? indexedNephrons : 1,
				allowedAxialSections);
			UG_LOGN("Cross-sectional axial-section check on level " << i + 1
			        << ": " << newAxialVertices
			        << " tube vertices on newly introduced axial parameters.");
			UG_COND_THROW(newAxialVertices != 0,
			              "Cross-sectional refinement introduced "
			              << newAxialVertices
			              << " tube vertices between the existing axial rings.");
		}
		if (fillInter && terminalMetadata.terminalsOnOuterWall)
			snap_refined_terminal_vertices_to_outer_wall(
				*dom.grid(), *dom.subset_handler(), refinedPos, refinedSurfParams,
				terminalMetadata,
				indexedNephrons > 0 ? indexedNephrons : 1, i + 1);
		if (fillInter && !crossSectionRefinement)
			remove_orphan_objects_after_projected_refinement(
				*dom.grid(), i + 1);
		if (fillInter)
			report_nephron_surface_projection_residual(
				*dom.grid(), *dom.subset_handler(), refinedPos, refinedSurfParams,
				*loadedNeuriteProjector,
				indexedNephrons > 0 ? indexedNephrons : 1, i + 1);
		validate_and_repair_orientation_after_projected_refinement(
			*dom.grid(), refinedPos, *dom.subset_handler(), i + 1,
			&refinedSurfParams);
		if (crossSectionRefinement)
			report_active_terminal_face_geometry(
				*dom.grid(), *dom.subset_handler(), refinedPos,
				indexedNephrons > 0 ? indexedNephrons : 1, i + 1);

		std::ostringstream oss;
		oss << "_refined_" << i+1 << ".ugx";
		const std::string curFileName = outFileNameBase + oss.str();
		// Match O-grid_working_FINAL: each generated refinement is a standalone,
		// flat mesh containing only that level.  The flow solver must load the
		// projector-backed coarse UGX and creates its geometric-multigrid hierarchy
		// in memory with RefineAndRebalanceDomain.  Serializing coarse parents here
		// makes ProMesh draw inactive parent faces on top of the refined mesh.
		const bool saved = SaveGridLevelToFile(
			*dom.grid(), *dom.subset_handler(), (int)i + 1,
			curFileName.c_str());
		UG_COND_THROW(!saved,
		              "Grid level writer failed for '" << curFileName << "'.");
		UG_LOGN("Wrote flat projected refinement level: " << curFileName);

	}
}


void import_er_neurites_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number erScaleFactor,
	number anisotropy,
	size_t numRefs
)
{
	import_er_neurites_from_swc_impl(fileNameIn, fileNameOut, erScaleFactor,
	                                 anisotropy, numRefs, -1.0, false, 2.0);
}


void import_nephron_with_membrane_and_inter_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number lumenScaleFactor,
	number anisotropy,
	size_t numRefs,
	number padding,
	number interTetQuality,
	size_t boxSurfaceRefs,
	size_t oGridVertices,
	bool coarseLumenCenterOnly,
	bool crossSectionRefinement
)
{
	UG_COND_THROW(lumenScaleFactor <= 0.0 || lumenScaleFactor >= 1.0,
	              "Lumen scale factor must be strictly between 0 and 1.");
	UG_COND_THROW(oGridVertices < 4,
	              "The nephron O-grid requires at least four circumferential vertices.");
	import_er_neurites_from_swc_impl(fileNameIn, fileNameOut, lumenScaleFactor,
	                                 anisotropy, numRefs, padding, true,
	                                 interTetQuality, 1, oGridVertices, true,
	                                 boxSurfaceRefs, coarseLumenCenterOnly,
	                                 crossSectionRefinement);
}


void repair_nephron_transport_boundary_subsets
(
	SmartPtr<Domain3d> dom
)
{
	UG_COND_THROW(!dom.valid(), "Cannot repair a null nephron domain.");
	MultiGrid& g = *dom->grid();
	ISubsetHandler& sh = *dom->subset_handler();
	const int membraneSI = sh.get_subset_index("Membrane");
	const int interSI = sh.get_subset_index("Inter");
	const int basolateralSI = sh.get_subset_index("Basolateral");
	UG_COND_THROW(membraneSI < 0 || interSI < 0 || basolateralSI < 0,
	              "The nephron mesh lacks Membrane, Inter, or Basolateral.");

	size_t repairedFaces = 0;
	for (FaceIterator it = g.begin<Face>(); it != g.end<Face>(); ++it)
	{
		Face* face = *it;
		if (sh.get_subset_index(face) != basolateralSI) continue;
		Grid::traits<Volume>::secure_container vols;
		g.associated_elements(vols, face);
		if (vols.size() != 1) continue;
		const int volumeSI = sh.get_subset_index(vols[0]);
		if (volumeSI != membraneSI && volumeSI != interSI) continue;
		sh.assign_subset(face, volumeSI);
		++repairedFaces;
	}

	auto incidentCompartmentMask = [&] (GridObject* obj) -> int
	{
		Grid::traits<Volume>::secure_container vols;
		g.associated_elements(vols, obj);
		int mask = 0;
		for (size_t i = 0; i < vols.size(); ++i)
		{
			const int si = sh.get_subset_index(vols[i]);
			if (si == membraneSI) mask |= 1;
			else if (si == interSI) mask |= 2;
		}
		return mask;
	};

	size_t repairedEdges = 0;
	for (EdgeIterator it = g.begin<Edge>(); it != g.end<Edge>(); ++it)
	{
		Edge* edge = *it;
		if (sh.get_subset_index(edge) != basolateralSI) continue;
		const int mask = incidentCompartmentMask(edge);
		if (mask == 1 || mask == 2)
		{
			sh.assign_subset(edge, mask == 1 ? membraneSI : interSI);
			++repairedEdges;
		}
	}

	size_t repairedVertices = 0;
	for (VertexIterator it = g.begin<Vertex>(); it != g.end<Vertex>(); ++it)
	{
		Vertex* vertex = *it;
		if (sh.get_subset_index(vertex) != basolateralSI) continue;
		const int mask = incidentCompartmentMask(vertex);
		if (mask == 1 || mask == 2)
		{
			sh.assign_subset(vertex, mask == 1 ? membraneSI : interSI);
			++repairedVertices;
		}
	}

	UG_LOGN("Nephron terminal-boundary repair: reclassified "
	        << repairedFaces << " one-sided Basolateral faces, "
	        << repairedEdges << " edges, and " << repairedVertices
	        << " vertices.");
}


void import_nephron_with_membrane_and_inter_no_projector_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number lumenScaleFactor,
	number anisotropy,
	size_t oGridVertices,
	number padding,
	number interTetQuality
)
{
	UG_COND_THROW(lumenScaleFactor <= 0.0 || lumenScaleFactor >= 1.0,
	              "Lumen scale factor must be strictly between 0 and 1.");
	UG_COND_THROW(oGridVertices < 4,
	              "The nephron O-grid requires at least four circumferential vertices.");

	// Position the initial curved nested tube once, then write no projection
	// handler and no neurite surface-parameter attachment. Later refinement is
	// therefore ordinary linear refinement of the generated mesh.
	import_er_neurites_from_swc_impl(fileNameIn, fileNameOut, lumenScaleFactor,
	                                 anisotropy, 0, padding, true,
	                                 interTetQuality, 1,
	                                 oGridVertices, false);
}


void import_multiple_nephrons_with_membrane_and_inter_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	size_t numNephrons,
	number lumenScaleFactor,
	number anisotropy,
	size_t numRefs,
	number padding,
	number interTetQuality,
	size_t boxSurfaceRefs,
	size_t oGridVertices,
	bool coarseLumenCenterOnly,
	bool crossSectionRefinement
)
{
	/*
	 * Lua-facing multi-nephron entry point. The Lua script has already smoothed
	 * each trace, assigned a constant outer radius, remapped node IDs, and merged
	 * the traces into one multi-root SWC. This function validates the requested
	 * number of disconnected trees and delegates to the common implementation.
	 */
	UG_COND_THROW(numNephrons < 2,
	              "The multi-nephron importer requires at least two SWC trees.");
	UG_COND_THROW(lumenScaleFactor <= 0.0 || lumenScaleFactor >= 1.0,
	              "Lumen scale factor must be strictly between 0 and 1.");
	UG_COND_THROW(oGridVertices < 4,
	              "The nephron O-grid requires at least four circumferential vertices.");
	import_er_neurites_from_swc_impl(fileNameIn, fileNameOut, lumenScaleFactor,
	                                 anisotropy, numRefs, padding, true,
	                                 interTetQuality, numNephrons, oGridVertices, true,
	                                 boxSurfaceRefs, coarseLumenCenterOnly,
	                                 crossSectionRefinement);
}



void import_1d_neurites_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number anisotropy,
	size_t numRefs,
	number scale
)
{
	// read in file to intermediate structure
	std::string inFileName = FindFileInStandardPaths(fileNameIn.c_str());
	UG_COND_THROW(inFileName == "", "File '" << fileNameIn
		<< "' could not be located in standard paths.");

	FileReaderSWC swcFileReader;
	swcFileReader.load_file(inFileName.c_str());
	std::vector<swc_types::SWCPoint>& vPoints = swcFileReader.swc_points();

	// scale
	const size_t sz = vPoints.size();
	for (size_t i = 0; i < sz; ++i)
	{
		VecScale(vPoints[i].coords, vPoints[i].coords, scale);
		vPoints[i].radius *= scale;
	}

	// smoothing
	//smoothing(vPoints, 5, 1.0, 1.0);

	// convert intermediate structure to neurite data
	std::vector<std::vector<vector3> > vPos;
	std::vector<std::vector<number> > vRad;
	std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
	std::vector<size_t> vRootNeuriteIndsOut;

	convert_pointlist_to_neuritelist(vPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

	// prepare grid and projector
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	Selector sel(g);


	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
		"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	if (!g.has_vertex_attachment(aSP))
		g.attach_to_vertices(aSP);

	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	aaSurfParams.access(g, aSP);

	UG_COND_THROW(!GlobalAttachments::is_declared("diameter"),
		"GlobalAttachment 'diameter' not declared.");
	Attachment<number> aDiam = GlobalAttachments::attachment<Attachment<number> >("diameter");
	if (!g.has_vertex_attachment(aDiam))
		g.attach_to_vertices(aDiam);

	Grid::VertexAttachmentAccessor<Attachment<number> > aaDiam;
	aaDiam.access(g, aDiam);
	ProjectionHandler projHandler(&sh);
	SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
	projHandler.set_geometry(geom3d);

	SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
	projHandler.set_projector(0, neuriteProj);

	// create spline data
	std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
	create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

	// create coarse grid
	for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i)
		create_neurite_1d(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
			anisotropy, g, aaPos, aaSurfParams, aaDiam, NULL);


	// subsets
	AssignSubsetColors(sh);
	sh.set_subset_name("neurites", 0);

	// export geometry
	std::string outFileNameBase = FilenameAndPathWithoutExtension(fileNameOut);
	std::string outFileName = outFileNameBase + ".ugx";
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	ugxWriter.add_projection_handler(projHandler, "defPH", 0);
	if (!ugxWriter.write_to_file(outFileName.c_str()))
		UG_THROW("Grid could not be written to file '" << outFileName << "'.");


	// refinements
	Domain3d dom;
	try {LoadDomain(dom, outFileName.c_str());}
	UG_CATCH_THROW("Failed loading domain from '" << outFileName << "'.");

	number offset = 2.0;
	std::string curFileName = outFileName.substr(0, outFileName.size()-4) + "_refined_0.ugx";
	try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
	UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

	GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
	for (size_t i = 0; i < numRefs; ++i)
	{
		ref.refine();
		std::ostringstream oss;
		oss << "_refined_" << i+1 << ".ugx";
		curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
		try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
		UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
	}
}




} // namespace neurites_from_swc
} // namespace neuro_collection
} // namespace ug
