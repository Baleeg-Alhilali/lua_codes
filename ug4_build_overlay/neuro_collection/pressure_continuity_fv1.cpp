#include "pressure_continuity_fv1.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary_impl.h"

namespace ug {
namespace neuro_collection {

template <typename TDomain>
PressureContinuityFV1<TDomain>::PressureContinuityFV1(
    const char* functions, const char* subsets)
    : base_type(functions, subsets), m_pressureScale(1.0),
      m_multiplierRegularization(0.0), m_augmentedLagrangian(0.0)
{}

template <typename TDomain>
void PressureContinuityFV1<TDomain>::set_pressure_scale(number scale)
{
    if (scale <= 0.0)
        UG_THROW("PressureContinuityFV1: pressure scale must be positive.");
    m_pressureScale = scale;
}

template <typename TDomain>
void PressureContinuityFV1<TDomain>::set_multiplier_regularization(
    number epsilon)
{
    if (epsilon < 0.0)
        UG_THROW("PressureContinuityFV1: multiplier regularization must be "
                 "non-negative.");
    m_multiplierRegularization = epsilon;
}

template <typename TDomain>
void PressureContinuityFV1<TDomain>::set_augmented_lagrangian(number gamma)
{
    if (gamma < 0.0)
        UG_THROW("PressureContinuityFV1: augmented-Lagrangian coefficient "
                 "must be non-negative.");
    m_augmentedLagrangian = gamma;
}

template <typename TDomain>
bool PressureContinuityFV1<TDomain>::fluxDensityFct(
    const std::vector<LocalVector::value_type>& u,
    GridObject*, const MathVector<dim>&, int, FluxCond& fc)
{
    if (u.size() != 3)
        UG_THROW("PressureContinuityFV1 needs exactly three functions: "
                 "p_stokes, p_darcy_stored, lambda.");

    fc.flux.resize(3);
    fc.from.resize(3);
    fc.to.assign(3, (size_t)-1);

    // The constraint row uses physical p_D = scale*p_Dstored. The Darcy
    // conservation equation itself is not row-scaled when its unknown is
    // changed, so its multiplier coupling remains -lambda. This makes the
    // physical solution invariant under algebraic Darcy-pressure scaling.
    const LocalVector::value_type constraint =
        u[0] - m_pressureScale * u[1];
    // Adding +/- gamma*constraint to the two conservation rows is an
    // augmented-Lagrangian stabilization: it changes the linear algebra but
    // not the converged solution because the third row enforces constraint=0.
    fc.flux[0] = u[2] + m_augmentedLagrangian * constraint;
    fc.from[0] = 0;
    fc.flux[1] = -u[2] - m_augmentedLagrangian * constraint;
    fc.from[1] = 1;
    // epsilon = 0 gives the exact Lagrange-multiplier constraint.  A positive
    // epsilon gives the standard stabilized saddle-point form
    // p_S - p_D - epsilon*lambda = 0 and supplies a multiplier diagonal.
    fc.flux[2] = u[0] - m_pressureScale * u[1]
                 - m_multiplierRegularization * u[2];
    fc.from[2] = 2;
    return true;
}

template <typename TDomain>
bool PressureContinuityFV1<TDomain>::fluxDensityFct(
    const std::vector<LocalVector::value_type>& u,
    const std::vector<LocalVector::value_type>&,
    GridObject* e, const MathVector<dim>& x, int si, FluxCond& fc)
{
    return fluxDensityFct(u, e, x, si, fc);
}

template <typename TDomain>
bool PressureContinuityFV1<TDomain>::fluxDensityDerivFct(
    const std::vector<LocalVector::value_type>& u,
    GridObject*, const MathVector<dim>&, int, FluxDerivCond& fdc)
{
    if (u.size() != 3)
        UG_THROW("PressureContinuityFV1 needs exactly three functions.");

    fdc.fluxDeriv.resize(3);
    fdc.from.resize(3);
    fdc.to.assign(3, (size_t)-1);

    fdc.from[0] = 0;
    fdc.fluxDeriv[0].push_back(std::make_pair((size_t)2, 1.0));
    if (m_augmentedLagrangian > 0.0) {
        fdc.fluxDeriv[0].push_back(
            std::make_pair((size_t)0, m_augmentedLagrangian));
        fdc.fluxDeriv[0].push_back(
            std::make_pair((size_t)1,
                           -m_augmentedLagrangian * m_pressureScale));
    }
    fdc.from[1] = 1;
    fdc.fluxDeriv[1].push_back(std::make_pair((size_t)2, -1.0));
    if (m_augmentedLagrangian > 0.0) {
        fdc.fluxDeriv[1].push_back(
            std::make_pair((size_t)0, -m_augmentedLagrangian));
        fdc.fluxDeriv[1].push_back(
            std::make_pair((size_t)1,
                           m_augmentedLagrangian * m_pressureScale));
    }
    fdc.from[2] = 2;
    fdc.fluxDeriv[2].push_back(std::make_pair((size_t)0, 1.0));
    fdc.fluxDeriv[2].push_back(
        std::make_pair((size_t)1, -m_pressureScale));
    if (m_multiplierRegularization > 0.0)
        fdc.fluxDeriv[2].push_back(
            std::make_pair((size_t)2, -m_multiplierRegularization));
    return true;
}

template <typename TDomain>
bool PressureContinuityFV1<TDomain>::fluxDensityDerivFct(
    const std::vector<LocalVector::value_type>& u,
    const std::vector<LocalVector::value_type>&,
    GridObject* e, const MathVector<dim>& x, int si, FluxDerivCond& fdc)
{
    return fluxDensityDerivFct(u, e, x, si, fdc);
}

#ifdef UG_DIM_1
template class PressureContinuityFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class PressureContinuityFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class PressureContinuityFV1<Domain3d>;
#endif

} // namespace neuro_collection
} // namespace ug
