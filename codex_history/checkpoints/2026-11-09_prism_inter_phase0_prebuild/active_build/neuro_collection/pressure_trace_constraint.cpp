#include "pressure_trace_constraint.h"

#include "common/error.h"
#include "common/util/string_util.h"
#include "lib_disc/common/multi_index.h"
#include "lib_grid/tools/subset_group.h"

#include <algorithm>
#include <map>

#ifdef UG_PARALLEL
#include "pcl/pcl_process_communicator.h"
#endif

namespace ug {
namespace neuro_collection {

template <typename TDomain, typename TAlgebra>
PressureTraceConstraint<TDomain, TAlgebra>::PressureTraceConstraint(
    const char* stokesPressure, const char* darcyStoredPressure,
    const char* interfaceSubset)
    : m_stokesPressure(stokesPressure),
      m_darcyStoredPressure(darcyStoredPressure),
      m_interfaceSubset(interfaceSubset),
      m_pressureScale(1.0),
      m_captureNextInterfaceDefect(false),
      m_interfaceDefectSums(3, 0.0)
{}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::set_pressure_scale(number scale)
{
    UG_COND_THROW(scale <= 0.0,
        "PressureTraceConstraint: pressure scale must be positive.");
    m_pressureScale = scale;
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::
capture_next_interface_defect()
{
    m_captureNextInterfaceDefect = true;
    std::fill(m_interfaceDefectSums.begin(),
              m_interfaceDefectSums.end(), 0.0);
}

template <typename TDomain, typename TAlgebra>
std::vector<number> PressureTraceConstraint<TDomain, TAlgebra>::
interface_defect_sums() const
{
    return m_interfaceDefectSums;
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::collect_pairs(
    ConstSmartPtr<DoFDistribution> dd,
    std::vector<std::pair<size_t, size_t> >& pairs) const
{
    UG_COND_THROW(TAlgebra::blockSize != 1,
        "PressureTraceConstraint currently requires scalar algebra blocks.");

    const size_t fctS = dd->fct_id_by_name(m_stokesPressure.c_str());
    const size_t fctD = dd->fct_id_by_name(m_darcyStoredPressure.c_str());
    SubsetGroup subsets(
        dd->subset_handler(), TokenizeTrimString(m_interfaceSubset));
    UG_COND_THROW(subsets.size() == 0,
        "PressureTraceConstraint requires at least one interface subset.");

    pairs.clear();
    std::vector<DoFIndex> indS;
    std::vector<DoFIndex> indD;
    typedef typename DoFDistribution::traits<Vertex>::const_iterator iterator;
    for (size_t si = 0; si < subsets.size(); ++si)
    {
        iterator it = dd->begin<Vertex>(subsets[si]);
        iterator itEnd = dd->end<Vertex>(subsets[si]);
        for (; it != itEnd; ++it)
        {
            dd->dof_indices(*it, fctS, indS, true, true);
            dd->dof_indices(*it, fctD, indD, true, true);
            UG_COND_THROW(indS.size() != 1 || indD.size() != 1,
                "PressureTraceConstraint expects one P1 DoF per pressure at "
                "each interface vertex.");
            pairs.push_back(std::make_pair(indS[0][0], indD[0][0]));
        }
    }
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::constrain_matrix(
    matrix_type& A, ConstSmartPtr<DoFDistribution> dd) const
{
    typedef typename matrix_type::value_type block_type;
    typedef typename matrix_type::row_iterator row_iterator;
    const number c = 1.0 / m_pressureScale;
    std::vector<std::pair<size_t, size_t> > pairs;
    collect_pairs(dd, pairs);

    // Eliminate all interface slaves simultaneously.  Doing this pair by
    // pair is incorrect: eliminating one slave can insert an entry into a
    // different slave column which is not present in that slave's row.  A
    // later row-pattern-based elimination then misses the new entry.  The
    // resulting matrix depends on interface-vertex ordering and can imprint
    // coarse-grid patches on the pressure solution.
    std::map<size_t, size_t> slaveToMaster;
    std::vector<std::vector<std::pair<size_t, block_type> > > slaveRows(
        pairs.size());
    for (size_t p = 0; p < pairs.size(); ++p)
    {
        const size_t master = pairs[p].first;
        const size_t slave = pairs[p].second;
        UG_COND_THROW(master >= A.num_rows() || slave >= A.num_rows(),
            "PressureTraceConstraint: interface DoF is outside the matrix.");
        UG_COND_THROW(slaveToMaster.find(slave) != slaveToMaster.end(),
            "PressureTraceConstraint: duplicate Darcy pressure slave DoF.");
        slaveToMaster[slave] = master;
        for (row_iterator it = A.begin_row(slave); it != A.end_row(slave); ++it)
            slaveRows[p].push_back(std::make_pair(it.index(), it.value()));
    }

    // Retain both conservation laws in physical residual units by adding
    // every original Darcy row to its matching Stokes row.  All source rows
    // were snapshotted above, so no earlier constraint can contaminate a
    // later one.
    for (size_t p = 0; p < pairs.size(); ++p)
    {
        const size_t master = pairs[p].first;
        for (size_t k = 0; k < slaveRows[p].size(); ++k)
        {
            A(master, slaveRows[p][k].first) += slaveRows[p][k].second;
        }
    }

    // Substitute p_D,stored = p_S/scale in every unconstrained equation.
    // Scan every matrix row instead of assuming that a slave-row sparsity
    // pattern also contains every entry in the corresponding slave column.
    for (size_t i = 0; i < A.num_rows(); ++i)
    {
        if (slaveToMaster.find(i) != slaveToMaster.end()) continue;

        std::vector<std::pair<size_t, block_type> > replacements;
        std::vector<size_t> slaveColumns;
        for (row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
        {
            typename std::map<size_t, size_t>::const_iterator found =
                slaveToMaster.find(it.index());
            if (found == slaveToMaster.end()) continue;
            replacements.push_back(
                std::make_pair(found->second, c * it.value()));
            slaveColumns.push_back(it.index());
        }
        for (size_t k = 0; k < slaveColumns.size(); ++k)
            A(i, slaveColumns[k]) = 0.0;
        for (size_t k = 0; k < replacements.size(); ++k)
            A(i, replacements[k].first) += replacements[k].second;
    }

    // Finally replace every Darcy trace row by the exact continuity equation.
    for (size_t p = 0; p < pairs.size(); ++p)
    {
        const size_t master = pairs[p].first;
        const size_t slave = pairs[p].second;
        SetRow(A, slave, 0.0);
        A(slave, slave) = 1.0;
        A(slave, master) = -c;
    }
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::reduce_vector(
    vector_type& vec, ConstSmartPtr<DoFDistribution> dd) const
{
    std::vector<std::pair<size_t, size_t> > pairs;
    collect_pairs(dd, pairs);
    for (size_t p = 0; p < pairs.size(); ++p)
    {
        vec[pairs[p].first] += vec[pairs[p].second];
        vec[pairs[p].second] = 0.0;
    }
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::set_slave_values(
    vector_type& vec, ConstSmartPtr<DoFDistribution> dd) const
{
    const number c = 1.0 / m_pressureScale;
    std::vector<std::pair<size_t, size_t> > pairs;
    collect_pairs(dd, pairs);
    for (size_t p = 0; p < pairs.size(); ++p)
        vec[pairs[p].second] = c * vec[pairs[p].first];
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::adjust_linear(
    matrix_type& mat, vector_type& rhs, ConstSmartPtr<DoFDistribution> dd,
    int, number)
{
    constrain_matrix(mat, dd);
    reduce_vector(rhs, dd);
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::adjust_jacobian(
    matrix_type& J, const vector_type&, ConstSmartPtr<DoFDistribution> dd,
    int, number, ConstSmartPtr<VectorTimeSeries<vector_type> >, const number)
{
    constrain_matrix(J, dd);
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::adjust_defect(
    vector_type& d, const vector_type& u, ConstSmartPtr<DoFDistribution> dd,
    int, number, ConstSmartPtr<VectorTimeSeries<vector_type> >,
    const std::vector<number>*, const std::vector<number>*)
{
    if (m_captureNextInterfaceDefect)
    {
        vector_type uniqueDefect(d);
#ifdef UG_PARALLEL
        UG_COND_THROW(!uniqueDefect.change_storage_type(PST_UNIQUE),
            "PressureTraceConstraint: could not create a unique interface "
            "defect for diagnostics.");
#endif

        std::vector<std::pair<size_t, size_t> > diagnosticPairs;
        collect_pairs(dd, diagnosticPairs);
        number localSums[2] = {0.0, 0.0};
        for (size_t p = 0; p < diagnosticPairs.size(); ++p)
        {
            localSums[0] += BlockRef(
                uniqueDefect[diagnosticPairs[p].first], 0);
            localSums[1] += BlockRef(
                uniqueDefect[diagnosticPairs[p].second], 0);
        }

        number globalSums[2] = {localSums[0], localSums[1]};
#ifdef UG_PARALLEL
        pcl::ProcessCommunicator communicator;
        communicator.allreduce(
            localSums, globalSums, 2, PCL_DT_DOUBLE, PCL_RO_SUM);
#endif
        m_interfaceDefectSums[0] = globalSums[0];
        m_interfaceDefectSums[1] = globalSums[1];
        m_interfaceDefectSums[2] = globalSums[0] + globalSums[1];
        m_captureNextInterfaceDefect = false;
    }

    reduce_vector(d, dd);
    const number c = 1.0 / m_pressureScale;
    std::vector<std::pair<size_t, size_t> > pairs;
    collect_pairs(dd, pairs);
    for (size_t p = 0; p < pairs.size(); ++p)
        d[pairs[p].second] = c * u[pairs[p].first] - u[pairs[p].second];
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::adjust_rhs(
    vector_type& rhs, const vector_type&, ConstSmartPtr<DoFDistribution> dd,
    int, number)
{
    reduce_vector(rhs, dd);
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::adjust_solution(
    vector_type& u, ConstSmartPtr<DoFDistribution> dd, int, number)
{
    set_slave_values(u, dd);
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::adjust_correction(
    vector_type& c, ConstSmartPtr<DoFDistribution> dd, int, number)
{
    set_slave_values(c, dd);
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::adjust_prolongation(
    matrix_type& P, ConstSmartPtr<DoFDistribution> ddFine,
    ConstSmartPtr<DoFDistribution> ddCoarse, int, number)
{
    // A fine-grid correction must remain in the constrained space. Standard
    // component-wise P1 transfer otherwise interpolates the two pressure
    // traces independently at newly created Apical vertices.
    typedef typename matrix_type::value_type block_type;
    typedef typename matrix_type::row_iterator row_iterator;
    const number c = 1.0 / m_pressureScale;
    std::vector<std::pair<size_t, size_t> > finePairs;
    std::vector<std::pair<size_t, size_t> > coarsePairs;
    collect_pairs(ddFine, finePairs);
    collect_pairs(ddCoarse, coarsePairs);
    std::map<size_t, size_t> coarseSlaveToMaster;
    for (size_t p = 0; p < coarsePairs.size(); ++p)
        coarseSlaveToMaster[coarsePairs[p].second] = coarsePairs[p].first;

    for (size_t p = 0; p < finePairs.size(); ++p)
    {
        const size_t master = finePairs[p].first;
        const size_t slave = finePairs[p].second;
        std::vector<std::pair<size_t, block_type> > masterRow;
        std::vector<std::pair<size_t, block_type> > slaveRow;
        for (row_iterator it = P.begin_row(master); it != P.end_row(master); ++it)
            masterRow.push_back(std::make_pair(it.index(), it.value()));
        for (row_iterator it = P.begin_row(slave); it != P.end_row(slave); ++it)
            slaveRow.push_back(std::make_pair(it.index(), it.value()));

        // Keep the ordinary Darcy interpolation so the coarse slave remains
        // an algebraic unknown. Remove its value on the constrained subspace,
        // then add c times the Stokes interpolation. For p_D,c=c*p_S,c this
        // reduces exactly to p_D,f=c*p_S,f.
        for (size_t k = 0; k < slaveRow.size(); ++k)
        {
            typename std::map<size_t, size_t>::const_iterator found =
                coarseSlaveToMaster.find(slaveRow[k].first);
            if (found != coarseSlaveToMaster.end())
                P(slave, found->second) -= c * slaveRow[k].second;
        }
        for (size_t k = 0; k < masterRow.size(); ++k)
            P(slave, masterRow[k].first) += c * masterRow[k].second;
    }
}

template <typename TDomain, typename TAlgebra>
void PressureTraceConstraint<TDomain, TAlgebra>::adjust_restriction(
    matrix_type& R, ConstSmartPtr<DoFDistribution> ddCoarse,
    ConstSmartPtr<DoFDistribution>, int, number)
{
    // Retain standard component-wise restriction. The coarse slave row is the
    // exact pressure-continuity equation and must not be eliminated here.
}

#ifdef UG_CPU_1
#ifdef UG_DIM_1
template class PressureTraceConstraint<Domain1d, CPUAlgebra>;
#endif
#ifdef UG_DIM_2
template class PressureTraceConstraint<Domain2d, CPUAlgebra>;
#endif
#ifdef UG_DIM_3
template class PressureTraceConstraint<Domain3d, CPUAlgebra>;
#endif
#endif
#ifdef UG_CPU_2
#ifdef UG_DIM_1
template class PressureTraceConstraint<Domain1d, CPUBlockAlgebra<2> >;
#endif
#ifdef UG_DIM_2
template class PressureTraceConstraint<Domain2d, CPUBlockAlgebra<2> >;
#endif
#ifdef UG_DIM_3
template class PressureTraceConstraint<Domain3d, CPUBlockAlgebra<2> >;
#endif
#endif
#ifdef UG_CPU_3
#ifdef UG_DIM_1
template class PressureTraceConstraint<Domain1d, CPUBlockAlgebra<3> >;
#endif
#ifdef UG_DIM_2
template class PressureTraceConstraint<Domain2d, CPUBlockAlgebra<3> >;
#endif
#ifdef UG_DIM_3
template class PressureTraceConstraint<Domain3d, CPUBlockAlgebra<3> >;
#endif
#endif

} // namespace neuro_collection
} // namespace ug
