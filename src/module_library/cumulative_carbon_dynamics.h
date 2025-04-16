#ifndef STANDARDBML_CUMULATIVE_CARBON_DYNAMICS_H
#define STANDARDBML_CUMULATIVE_CARBON_DYNAMICS_H

#include "../framework/module.h"
#include "../framework/state_map.h"

namespace standardBML
{
/**
 * @class cumulative_carbon_dynamics
 *
 * @brief Enables calculations of cumulative water dynamics.
 *
 * Cumulative water dynamics will be included in the simulation output as
 * differential quantities called:
 * - ``'canopy_assimilation'``: The cumulative net CO2 assimilation
 * - ``'canopy_gross_assimilation'``: The cumulative gross CO2 assimilation
 * - ``'canopy_non_photorespiratory_CO2_release'``: The cumulative CO2 lost to
 *   non-photorespiratory CO2 release
 * - ``'canopy_photorespiration'``: The cumulative CO2 lost to photorespiration
 */
class cumulative_carbon_dynamics : public differential_module
{
   public:
    cumulative_carbon_dynamics(
        state_map const& input_quantities,
        state_map* output_quantities)
        : differential_module{},

          // Get references to input quantities
          canopy_assimilation_rate{get_input(input_quantities, "canopy_assimilation_rate")},
          canopy_gross_assimilation_rate{get_input(input_quantities, "canopy_gross_assimilation_rate")},
          canopy_photorespiration_rate{get_input(input_quantities, "canopy_photorespiration_rate")},
          canopy_RL_rate{get_input(input_quantities, "canopy_non_photorespiratory_CO2_release_rate")},

          // Get pointers to output quantities
          canopy_assimilation_op{get_op(output_quantities, "canopy_assimilation")},
          canopy_gross_assimilation_op{get_op(output_quantities, "canopy_gross_assimilation")},
          canopy_photorespiration_op{get_op(output_quantities, "canopy_photorespiration")},
          canopy_RL_op{get_op(output_quantities, "canopy_non_photorespiratory_CO2_release")}
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "cumulative_carbon_dynamics"; }

   private:
    // References to input quantities
    double const& canopy_assimilation_rate;
    double const& canopy_gross_assimilation_rate;
    double const& canopy_photorespiration_rate;
    double const& canopy_RL_rate;

    // Pointers to output quantities
    double* canopy_assimilation_op;
    double* canopy_gross_assimilation_op;
    double* canopy_photorespiration_op;
    double* canopy_RL_op;

    // Main operation
    void do_operation() const;
};

string_vector cumulative_carbon_dynamics::get_inputs()
{
    return {
        "canopy_assimilation_rate",                      // Mg / ha / hr
        "canopy_gross_assimilation_rate",                // Mg / ha / hr
        "canopy_non_photorespiratory_CO2_release_rate",  // Mg / ha / hr
        "canopy_photorespiration_rate"                   // Mg / ha / hr
    };
}

string_vector cumulative_carbon_dynamics::get_outputs()
{
    return {
        "canopy_assimilation",                      // Mg / ha / hr
        "canopy_gross_assimilation",                // Mg / ha / hr
        "canopy_non_photorespiratory_CO2_release",  // Mg / ha / hr
        "canopy_photorespiration"                   // Mg / ha / hr
    };
}

void cumulative_carbon_dynamics::do_operation() const
{
    update(canopy_assimilation_op, canopy_assimilation_rate);
    update(canopy_gross_assimilation_op, canopy_gross_assimilation_rate);
    update(canopy_photorespiration_op, canopy_photorespiration_rate);
    update(canopy_RL_op, canopy_RL_rate);
}

}  // namespace standardBML
#endif
