#ifndef MAINTENANCE_RESPIRATION_CALCULATOR_H
#define MAINTENANCE_RESPIRATION_CALCULATOR_H

#include "../framework/module.h"
#include "../framework/state_map.h"
#include "temperature_response_functions.h"  // for Q10_temprature_response

namespace standardBML
{
/**
 * @class maintenance_respiration_calculator
 *
 * @brief Calculates the rate of change in plant organ biomasses due to
 * maintenance respiration; these are referred to as "maintenance respiration
 * rates," and the associated quantity names end with `_mr_rate`.
 *
 * This module is intended to be used along with the `maintenance_respiration`
 * module.
 *
 * The amount that each plant component respires is determined as a percentage
 * of its current biomass. This ideas is from this paper:
 * https://doi.org/10.1016/j.fcr.2010.07.007
 */
class maintenance_respiration_calculator : public direct_module
{
   public:
    maintenance_respiration_calculator(
        state_map const& input_quantities,
        state_map* output_quantities)
        : direct_module{},

          // Get references to input quantities
          Grain{get_input(input_quantities, "Grain")},
          Leaf{get_input(input_quantities, "Leaf")},
          mrc_grain{get_input(input_quantities, "mrc_grain")},
          mrc_leaf{get_input(input_quantities, "mrc_leaf")},
          mrc_rhizome{get_input(input_quantities, "mrc_rhizome")},
          mrc_root{get_input(input_quantities, "mrc_root")},
          mrc_shell{get_input(input_quantities, "mrc_shell")},
          mrc_stem{get_input(input_quantities, "mrc_stem")},
          Rhizome{get_input(input_quantities, "Rhizome")},
          Root{get_input(input_quantities, "Root")},
          Shell{get_input(input_quantities, "Shell")},
          Stem{get_input(input_quantities, "Stem")},
          temp{get_input(input_quantities, "temp")},

          // Get pointers to output quantities
          Grain_mr_rate_op{get_op(output_quantities, "Grain_mr_rate")},
          Leaf_mr_rate_op{get_op(output_quantities, "Leaf_mr_rate")},
          Rhizome_mr_rate_op{get_op(output_quantities, "Rhizome_mr_rate")},
          Root_mr_rate_op{get_op(output_quantities, "Root_mr_rate")},
          Shell_mr_rate_op{get_op(output_quantities, "Shell_mr_rate")},
          Stem_mr_rate_op{get_op(output_quantities, "Stem_mr_rate")}
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "maintenance_respiration_calculator"; }

   private:
    // References to input quantities
    const double& Grain;
    const double& Leaf;
    const double& mrc_grain;
    const double& mrc_leaf;
    const double& mrc_rhizome;
    const double& mrc_root;
    const double& mrc_shell;
    const double& mrc_stem;
    const double& Rhizome;
    const double& Root;
    const double& Shell;
    const double& Stem;
    const double& temp;

    // Pointers to output quantities
    double* Grain_mr_rate_op;
    double* Leaf_mr_rate_op;
    double* Rhizome_mr_rate_op;
    double* Root_mr_rate_op;
    double* Shell_mr_rate_op;
    double* Stem_mr_rate_op;

    // Implement the pure virtual function do_operation():
    void do_operation() const override final;
};

string_vector maintenance_respiration_calculator::get_inputs()
{
    return {
        "Grain",        // Mg / ha
        "Leaf",         // Mg / ha
        "mrc_grain",    // kg / kg / hr
        "mrc_leaf",     // kg / kg / hr
        "mrc_rhizome",  // kg / kg / hr
        "mrc_root",     // kg / kg / hr
        "mrc_shell",    // kg / kg / hr
        "mrc_stem",     // kg / kg / hr
        "Rhizome",      // Mg / ha
        "Root",         // Mg / ha
        "Shell",        // Mg / ha
        "Stem",         // Mg / ha
        "temp"          // degree C
    };
}

string_vector maintenance_respiration_calculator::get_outputs()
{
    return {
        "Grain_mr_rate",    // Mg / ha / hr
        "Leaf_mr_rate",     // Mg / ha / hr
        "Rhizome_mr_rate",  // Mg / ha / hr
        "Root_mr_rate",     // Mg / ha / hr
        "Shell_mr_rate",    // Mg / ha / hr
        "Stem_mr_rate"      // Mg / ha / hr
    };
}

void maintenance_respiration_calculator::do_operation() const
{
    // Define the reference temperature for the Q10 function
    double constexpr Tref = 25.0;

    // Calculate the multiplier for the Q10 response
    double Q10 = Q10_temperature_response(temp, Tref);  // dimensionless

    // Calculate each maintenance respiration rate (mrr)
    double const Grain_mr_rate = Grain * mrc_grain * Q10;        // Mg / ha / hr
    double const Leaf_mr_rate = Leaf * mrc_leaf * Q10;           // Mg / ha / hr
    double const Rhizome_mr_rate = Rhizome * mrc_rhizome * Q10;  // Mg / ha / hr
    double const Root_mr_rate = Root * mrc_root * Q10;           // Mg / ha / hr
    double const Shell_mr_rate = Shell * mrc_shell * Q10;        // Mg / ha / hr
    double const Stem_mr_rate = Stem * mrc_stem * Q10;           // Mg / ha / hr

    update(Grain_mr_rate_op, Grain_mr_rate);      // Mg / ha / hr
    update(Leaf_mr_rate_op, Leaf_mr_rate);        // Mg / ha / hr
    update(Rhizome_mr_rate_op, Rhizome_mr_rate);  // Mg / ha / hr
    update(Root_mr_rate_op, Root_mr_rate);        // Mg / ha / hr
    update(Shell_mr_rate_op, Shell_mr_rate);      // Mg / ha / hr
    update(Stem_mr_rate_op, Stem_mr_rate);        // Mg / ha / hr
}

}  // namespace standardBML
#endif
