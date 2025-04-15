#ifndef STANDARDBML_CUMULATIVE_CARBON_DYNAMICS_H
#define STANDARDBML_CUMULATIVE_CARBON_DYNAMICS_H

#include "../framework/module.h"
#include "../framework/state_map.h"

namespace standardBML
{
/**
 * @class cumulative_carbon_dynamics
 *
 * @brief Enables calculations of cumulative carbon dynamics.
 *
 * Cumulative canopy gas exchange will be included in the simulation output as
 * differential quantities called:
 * - ``'canopy_assimilation'``: The cumulative net CO2 assimilation
 * - ``'canopy_gross_assimilation'``: The cumulative gross CO2 assimilation
 * - ``'canopy_non_photorespiratory_CO2_release'``: The cumulative CO2 lost to
 *   non-photorespiratory CO2 release
 * - ``'canopy_photorespiration'``: The cumulative CO2 lost to photorespiration
 *
 * Cumulative maintenance respiration will be included in the simulation output
 * as differential quantities called ``'Leaf_mr'``, ``'Stem_mr'``, etc, where
 * "mr" standard for "maintenance respiration."
 *
 * Cumulative growth respiration will be included in the simulation output
 * as differential quantities called ``'Stem_gr'``, ``'Root_gr'``, etc, where
 * "gr" standard for "growth respiration."
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
          Grain_mrr{get_input(input_quantities, "Grain_mrr")},
          Leaf_grr{get_input(input_quantities, "Leaf_grr")},
          Leaf_mrr{get_input(input_quantities, "Leaf_mrr")},
          Rhizome_grr{get_input(input_quantities, "Rhizome_grr")},
          Rhizome_mrr{get_input(input_quantities, "Rhizome_mrr")},
          Root_grr{get_input(input_quantities, "Root_grr")},
          Root_mrr{get_input(input_quantities, "Root_mrr")},
          Shell_mrr{get_input(input_quantities, "Shell_mrr")},
          Stem_grr{get_input(input_quantities, "Stem_grr")},
          Stem_mrr{get_input(input_quantities, "Stem_mrr")},

          // Get pointers to output quantities
          canopy_assimilation_op{get_op(output_quantities, "canopy_assimilation")},
          canopy_gross_assimilation_op{get_op(output_quantities, "canopy_gross_assimilation")},
          canopy_photorespiration_op{get_op(output_quantities, "canopy_photorespiration")},
          canopy_RL_op{get_op(output_quantities, "canopy_non_photorespiratory_CO2_release")},
          Grain_mr_op{get_op(output_quantities, "Grain_mr")},
          Leaf_gr_op{get_op(output_quantities, "Leaf_gr")},
          Leaf_mr_op{get_op(output_quantities, "Leaf_mr")},
          Rhizome_gr_op{get_op(output_quantities, "Rhizome_gr")},
          Rhizome_mr_op{get_op(output_quantities, "Rhizome_mr")},
          Root_gr_op{get_op(output_quantities, "Root_gr")},
          Root_mr_op{get_op(output_quantities, "Root_mr")},
          Shell_mr_op{get_op(output_quantities, "Shell_mr")},
          Stem_gr_op{get_op(output_quantities, "Stem_gr")},
          Stem_mr_op{get_op(output_quantities, "Stem_mr")}
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
    double const& Grain_mrr;
    double const& Leaf_grr;
    double const& Leaf_mrr;
    double const& Rhizome_grr;
    double const& Rhizome_mrr;
    double const& Root_grr;
    double const& Root_mrr;
    double const& Shell_mrr;
    double const& Stem_grr;
    double const& Stem_mrr;

    // Pointers to output quantities
    double* canopy_assimilation_op;
    double* canopy_gross_assimilation_op;
    double* canopy_photorespiration_op;
    double* canopy_RL_op;
    double* Grain_mr_op;
    double* Leaf_gr_op;
    double* Leaf_mr_op;
    double* Rhizome_gr_op;
    double* Rhizome_mr_op;
    double* Root_gr_op;
    double* Root_mr_op;
    double* Shell_mr_op;
    double* Stem_gr_op;
    double* Stem_mr_op;

    // Main operation
    void do_operation() const;
};

string_vector cumulative_carbon_dynamics::get_inputs()
{
    return {
        "canopy_assimilation_rate",                      // Mg / ha / hr
        "canopy_gross_assimilation_rate",                // Mg / ha / hr
        "canopy_non_photorespiratory_CO2_release_rate",  // Mg / ha / hr
        "canopy_photorespiration_rate",                  // Mg / ha / hr
        "Grain_mrr",                                     // Mg / ha / hr
        "Leaf_grr",                                      // Mg / ha / hr
        "Leaf_mrr",                                      // Mg / ha / hr
        "Rhizome_grr",                                   // Mg / ha / hr
        "Rhizome_mrr",                                   // Mg / ha / hr
        "Root_grr",                                      // Mg / ha / hr
        "Root_mrr",                                      // Mg / ha / hr
        "Shell_mrr",                                     // Mg / ha / hr
        "Stem_grr",                                      // Mg / ha / hr
        "Stem_mrr"                                       // Mg / ha / hr
    };
}

string_vector cumulative_carbon_dynamics::get_outputs()
{
    return {
        "canopy_assimilation",                      // Mg / ha / hr
        "canopy_gross_assimilation",                // Mg / ha / hr
        "canopy_non_photorespiratory_CO2_release",  // Mg / ha / hr
        "canopy_photorespiration",                  // Mg / ha / hr
        "Grain_mr",                                 // Mg / ha / hr
        "Leaf_gr",                                  // Mg / ha / hr
        "Leaf_mr",                                  // Mg / ha / hr
        "Rhizome_gr",                               // Mg / ha / hr
        "Rhizome_mr",                               // Mg / ha / hr
        "Root_gr",                                  // Mg / ha / hr
        "Root_mr",                                  // Mg / ha / hr
        "Shell_mr",                                 // Mg / ha / hr
        "Stem_gr",                                  // Mg / ha / hr
        "Stem_mr"                                   // Mg / ha / hr
    };
}

void cumulative_carbon_dynamics::do_operation() const
{
    update(canopy_assimilation_op, canopy_assimilation_rate);
    update(canopy_gross_assimilation_op, canopy_gross_assimilation_rate);
    update(canopy_photorespiration_op, canopy_photorespiration_rate);
    update(canopy_RL_op, canopy_RL_rate);
    update(Grain_mr_op, Grain_mrr);
    update(Leaf_gr_op, Leaf_grr);
    update(Leaf_mr_op, Leaf_mrr);
    update(Rhizome_gr_op, Rhizome_grr);
    update(Rhizome_mr_op, Rhizome_mrr);
    update(Root_gr_op, Root_grr);
    update(Root_mr_op, Root_mrr);
    update(Shell_mr_op, Shell_mrr);
    update(Stem_gr_op, Stem_grr);
    update(Stem_mr_op, Stem_mrr);
}

}  // namespace standardBML
#endif
