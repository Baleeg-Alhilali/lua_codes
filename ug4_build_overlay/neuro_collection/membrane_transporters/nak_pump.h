#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__NAK_PUMP_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__NAK_PUMP_H

#include "membrane_transporter_interface.h"

namespace ug { namespace neuro_collection {

/// Saturating Na+/K+-ATPase: 3 Na_i -> Na_o and 2 K_o -> K_i per cycle.
class NaKPump : public IMembraneTransporter
{
public:
    enum {_NAI_ = 0, _NAO_, _KI_, _KO_};

    NaKPump(const std::vector<std::string>& fcts);
    NaKPump(const char* fcts);
    virtual ~NaKPump() {}

    void set_max_flux(number value) {m_maxFlux = value;}
    void set_k_na(number value) {m_kNa = value;}
    void set_k_k(number value) {m_kK = value;}

    virtual void calc_flux(const std::vector<number>& u, GridObject* e,
                           std::vector<number>& flux) const;
    virtual void calc_flux_deriv(const std::vector<number>& u, GridObject* e,
        std::vector<std::vector<std::pair<size_t, number> > >& deriv) const;
    virtual size_t n_dependencies() const;
    virtual size_t n_fluxes() const {return 2;}
    virtual const std::pair<size_t,size_t> flux_from_to(size_t i) const;
    virtual const std::string name() const {return "NaKPump";}
    virtual void check_supplied_functions() const;
    virtual void print_units() const;

private:
    number cycle_rate(const std::vector<number>& u) const;
    number m_maxFlux;
    number m_kNa;
    number m_kK;
};

} }
#endif
