#include "nak_pump.h"
#include <cmath>

namespace ug { namespace neuro_collection {

NaKPump::NaKPump(const std::vector<std::string>& fcts)
    : IMembraneTransporter(fcts), m_maxFlux(1.0e-22), m_kNa(10.0), m_kK(1.5) {}

NaKPump::NaKPump(const char* fcts)
    : IMembraneTransporter(fcts), m_maxFlux(1.0e-22), m_kNa(10.0), m_kK(1.5) {}

number NaKPump::cycle_rate(const std::vector<number>& u) const
{
    const number naSat = u[_NAI_] / (u[_NAI_] + m_kNa);
    const number kSat = u[_KO_] / (u[_KO_] + m_kK);
    return m_maxFlux * naSat * naSat * naSat * kSat * kSat;
}

void NaKPump::calc_flux(const std::vector<number>& u, GridObject*,
                        std::vector<number>& flux) const
{
    const number rate = cycle_rate(u);
    flux[0] = 3.0 * rate; // Na_i -> Na_o
    flux[1] = 2.0 * rate; // K_o -> K_i
}

void NaKPump::calc_flux_deriv(const std::vector<number>& u, GridObject*,
    std::vector<std::vector<std::pair<size_t, number> > >& deriv) const
{
    const number nai = u[_NAI_], ko = u[_KO_];
    const number naDenom = nai + m_kNa;
    const number kDenom = ko + m_kK;
    const number naSat = nai / naDenom;
    const number kSat = ko / kDenom;
    const number dRateNa = m_maxFlux * 3.0 * naSat * naSat *
        m_kNa / (naDenom * naDenom) * kSat * kSat;
    const number dRateK = m_maxFlux * naSat * naSat * naSat * 2.0 * kSat *
        m_kK / (kDenom * kDenom);

    size_t j = 0;
    if (!has_constant_value(_NAI_)) deriv[0][j++] = std::make_pair(local_fct_index(_NAI_), 3.0*dRateNa);
    if (!has_constant_value(_KO_))  deriv[0][j++] = std::make_pair(local_fct_index(_KO_), 3.0*dRateK);
    j = 0;
    if (!has_constant_value(_NAI_)) deriv[1][j++] = std::make_pair(local_fct_index(_NAI_), 2.0*dRateNa);
    if (!has_constant_value(_KO_))  deriv[1][j++] = std::make_pair(local_fct_index(_KO_), 2.0*dRateK);
}

size_t NaKPump::n_dependencies() const
{
    return (has_constant_value(_NAI_) ? 0 : 1) + (has_constant_value(_KO_) ? 0 : 1);
}

const std::pair<size_t,size_t> NaKPump::flux_from_to(size_t i) const
{
    if (i == 0) return std::make_pair(local_fct_index(_NAI_), local_fct_index(_NAO_));
    return std::make_pair(local_fct_index(_KO_), local_fct_index(_KI_));
}

void NaKPump::check_supplied_functions() const
{
    if (m_vFct.size() != 4) UG_THROW("NaKPump requires {Na_i, Na_o, K_i, K_o}.");
}

void NaKPump::print_units() const
{
    UG_LOGN("NaKPump: concentrations and K values in mM; max flux is the cycle flux density before 3:2 stoichiometric scaling.");
}

} }
