#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__SOLVENT_DRAG_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__SOLVENT_DRAG_H

#include "membrane_transporter_interface.h"

namespace ug { namespace neuro_collection {

/// Solvent-drag concentration factor (1-sigma)*(c_source+c_target)/2.
/** MembraneTransportFV1 supplies the signed normal water velocity through
 *  set_flux_scale_normal_velocity(). */
class SolventDrag : public IMembraneTransporter
{
public:
    enum {_SOURCE_ = 0, _TARGET_};

    SolventDrag(const std::vector<std::string>& fcts);
    SolventDrag(const char* fcts);
    virtual ~SolventDrag() {}

    void set_reflection_coefficient(number sigma);

    virtual void calc_flux(const std::vector<number>& u, GridObject* e,
                           std::vector<number>& flux) const;
    virtual void calc_flux_deriv(const std::vector<number>& u, GridObject* e,
        std::vector<std::vector<std::pair<size_t, number> > >& deriv) const;
    virtual size_t n_dependencies() const;
    virtual size_t n_fluxes() const {return 1;}
    virtual const std::pair<size_t,size_t> flux_from_to(size_t i) const;
    virtual const std::string name() const {return "SolventDrag";}
    virtual void check_supplied_functions() const;
    virtual void print_units() const;

private:
    number m_sigma;
};

} }
#endif
