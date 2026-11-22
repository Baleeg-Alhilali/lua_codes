#ifndef UG__PLUGINS__NEURO_COLLECTION__PRESSURE_TRACE_CONSTRAINT_H
#define UG__PLUGINS__NEURO_COLLECTION__PRESSURE_TRACE_CONSTRAINT_H

#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/domain.h"

#include <string>
#include <utility>
#include <vector>

namespace ug {
namespace neuro_collection {

/// Eliminates the stored Darcy pressure trace by p_D,stored = p_S / scale.
///
/// The eliminated Darcy conservation row is added to the Stokes continuity
/// row without row scaling. Only the Darcy unknown/column carries 1/scale.
/// This retains both conservation laws
/// without introducing a zero-diagonal Lagrange multiplier.
template <typename TDomain, typename TAlgebra>
class PressureTraceConstraint
    : public IDomainConstraint<TDomain, TAlgebra>
{
public:
    typedef typename TAlgebra::matrix_type matrix_type;
    typedef typename TAlgebra::vector_type vector_type;

    PressureTraceConstraint(const char* stokesPressure,
                            const char* darcyStoredPressure,
                            const char* interfaceSubset);

    void set_pressure_scale(number scale);
    void capture_next_interface_defect();
    std::vector<number> interface_defect_sums() const;
    virtual int type() const { return CT_CONSTRAINTS; }

    virtual void adjust_jacobian(
        matrix_type& J, const vector_type& u,
        ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
        ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
        const number s_a0 = 1.0);

    virtual void adjust_defect(
        vector_type& d, const vector_type& u,
        ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
        ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
        const std::vector<number>* vScaleMass = NULL,
        const std::vector<number>* vScaleStiff = NULL);

    virtual void adjust_linear(
        matrix_type& mat, vector_type& rhs,
        ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0);

    virtual void adjust_rhs(
        vector_type& rhs, const vector_type& u,
        ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0);

    virtual void adjust_solution(
        vector_type& u, ConstSmartPtr<DoFDistribution> dd,
        int type, number time = 0.0);

    virtual void adjust_correction(
        vector_type& c, ConstSmartPtr<DoFDistribution> dd,
        int type, number time = 0.0);

    virtual void adjust_prolongation(
        matrix_type& P, ConstSmartPtr<DoFDistribution> ddFine,
        ConstSmartPtr<DoFDistribution> ddCoarse, int type,
        number time = 0.0);

    virtual void adjust_restriction(
        matrix_type& R, ConstSmartPtr<DoFDistribution> ddCoarse,
        ConstSmartPtr<DoFDistribution> ddFine, int type,
        number time = 0.0);

private:
    void collect_pairs(
        ConstSmartPtr<DoFDistribution> dd,
        std::vector<std::pair<size_t, size_t> >& pairs) const;
    void constrain_matrix(matrix_type& mat,
                          ConstSmartPtr<DoFDistribution> dd) const;
    void reduce_vector(vector_type& vec,
                       ConstSmartPtr<DoFDistribution> dd) const;
    void set_slave_values(vector_type& vec,
                          ConstSmartPtr<DoFDistribution> dd) const;

private:
    std::string m_stokesPressure;
    std::string m_darcyStoredPressure;
    std::string m_interfaceSubset;
    number m_pressureScale;
    bool m_captureNextInterfaceDefect;
    std::vector<number> m_interfaceDefectSums;
};

} // namespace neuro_collection
} // namespace ug

#endif
