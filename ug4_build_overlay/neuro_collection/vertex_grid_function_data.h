/*
 * VertexGridFunctionNumberData: value of a P1 grid-function component at a
 * mesh vertex, looked up by the vertex's exact position.
 *
 * Intended for data evaluated at vertex positions, e.g. the density function
 * of MembraneTransportFV1 (evaluated at the boundary-face corner = node).
 * GlobalGridFunctionNumberData locates the point by an element search; its
 * edge variant tests only the x coordinate and its volume variant fails at
 * curved hexahedral/prismatic corners, so they return wrong values or throw
 * there.  This class maps every local vertex (of the subsets where the
 * component is defined) by its quantized coordinates to its DoF, and returns
 * the current nodal value.  Points that are not a vertex (within the
 * tolerance) throw.
 *
 * Added 2026-11-27 for the nephron solvent-drag water flux.
 */

#ifndef __H__UG__PLUGINS__NEURO_COLLECTION__VERTEX_GRID_FUNCTION_DATA_H__
#define __H__UG__PLUGINS__NEURO_COLLECTION__VERTEX_GRID_FUNCTION_DATA_H__

#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "common/common.h"
#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug {
namespace neuro_collection {

template <typename TGridFunction>
class VertexGridFunctionNumberData
	: public StdGlobPosData<VertexGridFunctionNumberData<TGridFunction>,
	                        number, TGridFunction::dim>
{
	public:
		static const int dim = TGridFunction::dim;
		typedef std::array<long long, dim> key_type;

		VertexGridFunctionNumberData(SmartPtr<TGridFunction> spGridFct,
		                             const char* cmp)
			: m_spGridFct(spGridFct), m_tol(1.0e-9)
		{
			m_fct = spGridFct->fct_id_by_name(cmp);
			typename TGridFunction::domain_type::position_accessor_type& aaPos
				= spGridFct->domain()->position_accessor();

			std::vector<DoFIndex> ind;
			for (int si = 0; si < spGridFct->num_subsets(); ++si)
			{
				if (!spGridFct->is_def_in_subset(m_fct, si)) continue;
				typedef typename TGridFunction::template traits<Vertex>::const_iterator iter_t;
				iter_t iter = spGridFct->template begin<Vertex>(si);
				iter_t iterEnd = spGridFct->template end<Vertex>(si);
				for (; iter != iterEnd; ++iter)
				{
					Vertex* vrt = *iter;
					spGridFct->inner_dof_indices(vrt, m_fct, ind);
					if (ind.empty()) continue;
					Entry e;
					e.pos = aaPos[vrt];
					e.dof = ind[0];
					m_map[key(e.pos)].push_back(e);
				}
			}
		}

		inline void evaluate(number& value, const MathVector<dim>& x,
		                     number time, int si) const
		{
			const Entry* e = find(x);
			UG_COND_THROW(e == NULL, "VertexGridFunctionNumberData: no vertex at "
				<< x << " (tolerance " << m_tol << ").");
			value = DoFRef(*m_spGridFct, e->dof);
		}

		inline number evaluate(const MathVector<dim>& x) const
		{
			number value;
			evaluate(value, x, 0.0, -1);
			return value;
		}

		/// number of mapped vertices on this process
		size_t num_vertices() const
		{
			size_t n = 0;
			for (typename map_type::const_iterator it = m_map.begin();
			     it != m_map.end(); ++it)
				n += it->second.size();
			return n;
		}

		bool continuous() const {return true;}
		bool requires_grid_fct() const {return false;}

	private:
		struct Entry
		{
			MathVector<dim> pos;
			DoFIndex dof;
		};
		typedef std::map<key_type, std::vector<Entry> > map_type;

		key_type key(const MathVector<dim>& x) const
		{
			key_type k;
			for (int d = 0; d < dim; ++d)
				k[d] = (long long) std::floor(x[d] / m_tol);
			return k;
		}

		const Entry* find(const MathVector<dim>& x) const
		{
			const key_type k0 = key(x);
			// own cell and its neighbours (a point may sit on a cell border)
			const int nNb = (dim == 3) ? 27 : (dim == 2 ? 9 : 3);
			for (int n = 0; n < nNb; ++n)
			{
				key_type k = k0;
				int r = n;
				for (int d = 0; d < dim; ++d) {k[d] += (r % 3) - 1; r /= 3;}
				typename map_type::const_iterator it = m_map.find(k);
				if (it == m_map.end()) continue;
				for (size_t i = 0; i < it->second.size(); ++i)
					if (VecDistance(it->second[i].pos, x) <= m_tol)
						return &it->second[i];
			}
			return NULL;
		}

		SmartPtr<TGridFunction> m_spGridFct;
		size_t m_fct;
		number m_tol;
		map_type m_map;
};

} // namespace neuro_collection
} // namespace ug

#endif
