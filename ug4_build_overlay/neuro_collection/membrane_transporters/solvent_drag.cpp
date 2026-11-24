#include "solvent_drag.h"

namespace ug { namespace neuro_collection {

SolventDrag::SolventDrag(const std::vector<std::string>& fcts)
    : IMembraneTransporter(fcts), m_sigma(1.0) {}

SolventDrag::SolventDrag(const char* fcts)
    : IMembraneTransporter(fcts), m_sigma(1.0) {}

void SolventDrag::set_reflection_coefficient(number sigma)
{
    UG_COND_THROW(sigma < 0.0 || sigma > 1.0,
        "SolventDrag reflection coefficient must be in [0,1].");
    m_sigma = sigma;
}

void SolventDrag::calc_flux(const std::vector<number>& u, GridObject*,
                            std::vector<number>& flux) const
{
    flux[0] = (1.0 - m_sigma) * 0.5 *
        (u[_SOURCE_] + u[_TARGET_]);
}

void SolventDrag::calc_flux_deriv(const std::vector<number>&, GridObject*,
    std::vector<std::vector<std::pair<size_t, number> > >& deriv) const
{
    const number d = 0.5 * (1.0 - m_sigma);
    size_t j = 0;
    if (!has_constant_value(_SOURCE_))
        deriv[0][j++] = std::make_pair(local_fct_index(_SOURCE_), d);
    if (!has_constant_value(_TARGET_))
        deriv[0][j++] = std::make_pair(local_fct_index(_TARGET_), d);
}

size_t SolventDrag::n_dependencies() const
{
    return (has_constant_value(_SOURCE_) ? 0 : 1) +
           (has_constant_value(_TARGET_) ? 0 : 1);
}

const std::pair<size_t,size_t> SolventDrag::flux_from_to(size_t) const
{
    return std::make_pair(local_fct_index(_SOURCE_),
                          local_fct_index(_TARGET_));
}

void SolventDrag::check_supplied_functions() const
{
    if (m_vFct.size() != 2)
        UG_THROW("SolventDrag requires {source concentration, target concentration}.");
}

void SolventDrag::print_units() const
{
    UG_LOGN("SolventDrag: concentration times signed normal water velocity.");
}

} }
