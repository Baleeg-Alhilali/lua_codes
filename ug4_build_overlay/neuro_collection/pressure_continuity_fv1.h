#ifndef UG__PLUGINS__NEURO_COLLECTION__PRESSURE_CONTINUITY_FV1_H
#define UG__PLUGINS__NEURO_COLLECTION__PRESSURE_CONTINUITY_FV1_H

#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"

namespace ug {
namespace neuro_collection {

/// FV1 saddle-point constraint for p_stokes = scale * p_darcy on an interface.
/**
 * Function order is: Stokes pressure, stored Darcy pressure, multiplier.
 * The multiplier preserves the separate Stokes and Darcy conservation rows;
 * unlike sharing one pressure DoF, no interface equation is discarded.
 */
template <typename TDomain>
class PressureContinuityFV1
    : public FV1InnerBoundaryElemDisc<PressureContinuityFV1<TDomain>, TDomain>
{
public:
    typedef TDomain domain_type;
    static const int dim = domain_type::dim;

private:
    typedef FV1InnerBoundaryElemDisc<PressureContinuityFV1<TDomain>, TDomain> base_type;
    typedef typename base_type::FluxCond FluxCond;
    typedef typename base_type::FluxDerivCond FluxDerivCond;
    number m_pressureScale;
    number m_multiplierRegularization;
    number m_augmentedLagrangian;

public:
    PressureContinuityFV1(const char* functions, const char* subsets);
    void set_pressure_scale(number scale);
    void set_multiplier_regularization(number epsilon);
    void set_augmented_lagrangian(number gamma);

    bool fluxDensityFct(
        const std::vector<LocalVector::value_type>& u,
        GridObject*, const MathVector<dim>&, int, FluxCond& fc);
    bool fluxDensityFct(
        const std::vector<LocalVector::value_type>& u,
        const std::vector<LocalVector::value_type>&,
        GridObject* e, const MathVector<dim>& x, int si, FluxCond& fc);
    bool fluxDensityDerivFct(
        const std::vector<LocalVector::value_type>& u,
        GridObject*, const MathVector<dim>&, int, FluxDerivCond& fdc);
    bool fluxDensityDerivFct(
        const std::vector<LocalVector::value_type>& u,
        const std::vector<LocalVector::value_type>&,
        GridObject* e, const MathVector<dim>& x, int si, FluxDerivCond& fdc);
};

} // namespace neuro_collection
} // namespace ug

#endif
