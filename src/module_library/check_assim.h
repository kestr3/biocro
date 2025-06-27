#ifndef CHECK_ASSIM_H
#define CHECK_ASSIM_H

#include "../framework/module.h"
#include "../framework/state_map.h"
#include "c3_temperature_response.h"  // for c3_temperature_response_parameters
#include "c3photo.h"

namespace standardBML
{

class c3_assim_check : public direct_module
{
   public:
    c3_assim_check(
        state_map const& input_quantities,
        state_map* output_quantities)
        : direct_module{},

          // Get pointers to input quantities
          atmospheric_pressure{get_input(input_quantities, "atmospheric_pressure")},
          b0{get_input(input_quantities, "b0")},
          b1{get_input(input_quantities, "b1")},
          beta_PSII{get_input(input_quantities, "beta_PSII")},
          Catm{get_input(input_quantities, "Catm")},
          Ci{get_input(input_quantities, "Ci")},
          electrons_per_carboxylation{get_input(input_quantities, "electrons_per_carboxylation")},
          electrons_per_oxygenation{get_input(input_quantities, "electrons_per_oxygenation")},
          gbw{get_input(input_quantities, "gbw")},
          Gs_min{get_input(input_quantities, "Gs_min")},
          Gstar_c{get_input(input_quantities, "Gstar_c")},
          Gstar_Ea{get_input(input_quantities, "Gstar_Ea")},
          Jmax_at_25{get_input(input_quantities, "Jmax_at_25")},
          Jmax_c{get_input(input_quantities, "Jmax_c")},
          Jmax_Ea{get_input(input_quantities, "Jmax_Ea")},
          Kc_c{get_input(input_quantities, "Kc_c")},
          Kc_Ea{get_input(input_quantities, "Kc_Ea")},
          Ko_c{get_input(input_quantities, "Ko_c")},
          Ko_Ea{get_input(input_quantities, "Ko_Ea")},
          O2{get_input(input_quantities, "O2")},
          phi_PSII_0{get_input(input_quantities, "phi_PSII_0")},
          phi_PSII_1{get_input(input_quantities, "phi_PSII_1")},
          phi_PSII_2{get_input(input_quantities, "phi_PSII_2")},
          Qabs{get_input(input_quantities, "Qabs")},
          rh{get_input(input_quantities, "rh")},
          RL_at_25{get_input(input_quantities, "RL_at_25")},
          RL_c{get_input(input_quantities, "RL_c")},
          RL_Ea{get_input(input_quantities, "RL_Ea")},
          StomataWS{get_input(input_quantities, "StomataWS")},
          Tambient{get_input(input_quantities, "temp")},
          theta_0{get_input(input_quantities, "theta_0")},
          theta_1{get_input(input_quantities, "theta_1")},
          theta_2{get_input(input_quantities, "theta_2")},
          Tleaf{get_input(input_quantities, "Tleaf")},
          Tp_at_25{get_input(input_quantities, "Tp_at_25")},
          Tp_c{get_input(input_quantities, "Tp_c")},
          Tp_Ha{get_input(input_quantities, "Tp_Ha")},
          Tp_Hd{get_input(input_quantities, "Tp_Hd")},
          Tp_S{get_input(input_quantities, "Tp_S")},
          Vcmax_at_25{get_input(input_quantities, "Vcmax_at_25")},
          Vcmax_c{get_input(input_quantities, "Vcmax_c")},
          Vcmax_Ea{get_input(input_quantities, "Vcmax_Ea")},

          // Get pointers to output quantities
          Assim_check_op{get_op(output_quantities, "Assim_check")}
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "c3_assim_check"; }

   private:
    // References to input quantities
    double const& atmospheric_pressure;
    double const& b0;
    double const& b1;
    double const& beta_PSII;
    double const& Catm;
    double const& Ci;
    double const& electrons_per_carboxylation;
    double const& electrons_per_oxygenation;
    double const& gbw;
    double const& Gs_min;
    double const& Gstar_c;
    double const& Gstar_Ea;
    double const& Jmax_at_25;
    double const& Jmax_c;
    double const& Jmax_Ea;
    double const& Kc_c;
    double const& Kc_Ea;
    double const& Ko_c;
    double const& Ko_Ea;
    double const& O2;
    double const& phi_PSII_0;
    double const& phi_PSII_1;
    double const& phi_PSII_2;
    double const& Qabs;
    double const& rh;
    double const& RL_at_25;
    double const& RL_c;
    double const& RL_Ea;
    double const& StomataWS;
    double const& Tambient;
    double const& theta_0;
    double const& theta_1;
    double const& theta_2;
    double const& Tleaf;
    double const& Tp_at_25;
    double const& Tp_c;
    double const& Tp_Ha;
    double const& Tp_Hd;
    double const& Tp_S;
    double const& Vcmax_at_25;
    double const& Vcmax_c;
    double const& Vcmax_Ea;

    // Pointers to output quantities
    double* Assim_op;
    double* Assim_check_op;
    double* Assim_conductance_op;
    double* Ci_op;
    double* Cs_op;
    double* GrossAssim_op;
    double* Gs_op;
    double* RHs_op;
    double* RL_op;
    double* Rp_op;
    double* iterations_op;

    // Main operation
    void do_operation() const;
};

string_vector c3_assim_check::get_inputs()
{
    return {
        "atmospheric_pressure",  // Pa
        "b0",                    // mol / m^2 / s
        "b1",                    // dimensionless
        "beta_PSII",             // dimensionless (fraction of absorbed light that reaches photosystem II)
        "Catm",                  // micromol / mol
        "Ci",
        "electrons_per_carboxylation",  // self-explanatory units
        "electrons_per_oxygenation",    // self-explanatory units
        "gbw",                          // mol / m^2 / s
        "Gs_min",                       // mol / m^2 / s
        "Gstar_c",                      // dimensionless
        "Gstar_Ea",                     // J / mol
        "Jmax_at_25",                   // micromol / m^2 / s
        "Jmax_c",                       // dimensionless
        "Jmax_Ea",                      // J / mol
        "Kc_c",                         // dimensionless
        "Kc_Ea",                        // J / mol
        "Ko_c",                         // dimensionless
        "Ko_Ea",                        // J / mol
        "O2",                           // millimol / mol
        "phi_PSII_0",                   // dimensionless
        "phi_PSII_1",                   // (degrees C)^(-1)
        "phi_PSII_2",                   // (degrees C)^(-2)
        "Qabs",                         // micromol / m^2 / s
        "rh",                           // dimensionless
        "RL_at_25",                     // micromol / m^2 / s
        "RL_c",                         // dimensionless
        "RL_Ea",                        // J / mol
        "StomataWS",                    // dimensionless
        "temp",                         // degrees C
        "theta_0",                      // dimensionless
        "theta_1",                      // (degrees C)^(-1)
        "theta_2",                      // (degrees C)^(-2)
        "Tleaf",                        // degrees C
        "Tp_at_25",                     // micromol / m^2 / s
        "Tp_c",                         // dimensionless
        "Tp_Ha",                        // J / mol
        "Tp_Hd",                        // J / mol
        "Tp_S",                         // J / K / mol
        "Vcmax_at_25",                  // micromol / m^2 / s
        "Vcmax_c",                      // dimensionless
        "Vcmax_Ea",                     // J / mol
    };
}

string_vector c3_assim_check::get_outputs()
{
    return {
        "Assim_check"};
}

void c3_assim_check::do_operation() const
{
    // Combine temperature response parameters
    c3_temperature_response_parameters const tr_param{
        Gstar_c,
        Gstar_Ea,
        Jmax_c,
        Jmax_Ea,
        Kc_c,
        Kc_Ea,
        Ko_c,
        Ko_Ea,
        phi_PSII_0,
        phi_PSII_1,
        phi_PSII_2,
        RL_c,
        RL_Ea,
        theta_0,
        theta_1,
        theta_2,
        Tp_c,
        Tp_Ha,
        Tp_Hd,
        Tp_S,
        Vcmax_c,
        Vcmax_Ea};

    double check = check_assim(
        Ci,
        tr_param,
        Qabs,
        Tleaf,
        Tambient,
        rh,
        Vcmax_at_25,
        Jmax_at_25,
        Tp_at_25,
        RL_at_25,
        b0,
        b1,
        Gs_min,
        Catm,
        atmospheric_pressure,
        O2,
        StomataWS,
        electrons_per_carboxylation,
        electrons_per_oxygenation,
        beta_PSII,
        gbw);

    // Update the output quantity list
    update(Assim_check_op, check);
}

}  // namespace standardBML
#endif
