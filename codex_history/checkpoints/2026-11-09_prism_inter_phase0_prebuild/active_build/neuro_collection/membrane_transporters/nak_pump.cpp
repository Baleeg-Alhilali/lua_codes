#include "nak_pump.h"
#include <cmath>

namespace ug { namespace neuro_collection {

NaKPump::NaKPump(const std::vector<std::string>& fcts)
    : IMembraneTransporter(fcts), m_maxFlux(1.0e-22), m_kNa(10.0), m_kK(1.5) {}

NaKPump::NaKPump(const char* fcts)
    : IMembraneTransporter(fcts), m_maxFlux(1.0e-22), m_kNa(10.0), m_kK(1.5) {}

number NaKPump::cycle_rate(const std::vector<number>& u) const
{
    const number nai3 = u[_NAI_] * u[_NAI_] * u[_NAI_];
    const number ko2 = u[_KO_] * u[_KO_];
    const number kna3 = m_kNa * m_kNa * m_kNa;
    const number kk2 = m_kK * m_kK;
    return m_maxFlux * nai3 / (kna3 + nai3) * ko2 / (kk2 + ko2);
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
    const number nai2 = nai*nai, nai3 = nai2*nai;
    const number ko2 = ko*ko;
    const number kna3 = m_kNa*m_kNa*m_kNa, kk2 = m_kK*m_kK;
    const number dRateNa = m_maxFlux * 3.0*nai2*kna3 /
        ((kna3+nai3)*(kna3+nai3)) * ko2/(kk2+ko2);
    const number dRateK = m_maxFlux * nai3/(kna3+nai3) * 2.0*ko*kk2 /
        ((kk2+ko2)*(kk2+ko2));

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
    UG_LOGN("NaKPump: concentrations and K values in mM; max flux in mol/s per pump.");
}

} }
