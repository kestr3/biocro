#ifndef MAINTENANCE_RESPIRATION_H
#define MAINTENANCE_RESPIRATION_H

#include "../framework/module.h"
#include "../framework/state_map.h"
#include "temperature_response_functions.h"  // for Q10_temprature_response

namespace standardBML
{
/**
 * @class maintenance_respiration
 *
 * @brief Calculates the change in plant organ biomasses due to maintenance respiration.
 *
 * The amount that each plant component respires is determined as a
 * percentage of its current biomass.
 * This ideas is from this paper:
 * (https://doi.org/10.1016/j.fcr.2010.07.007)
 */
class maintenance_respiration : public differential_module
{
   public:
    maintenance_respiration(
        state_map const& input_quantities,
        state_map* output_quantities)
        : differential_module{},

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
          Grain_op{get_op(output_quantities, "Grain")},
          Leaf_op{get_op(output_quantities, "Leaf")},
          Rhizome_op{get_op(output_quantities, "Rhizome")},
          Root_op{get_op(output_quantities, "Root")},
          Shell_op{get_op(output_quantities, "Shell")},
          Stem_op{get_op(output_quantities, "Stem")}
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "maintenance_respiration"; }

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
    double* Grain_op;
    double* Leaf_op;
    double* Rhizome_op;
    double* Root_op;
    double* Shell_op;
    double* Stem_op;

    // Implement the pure virtual function do_operation():
    void do_operation() const override final;
};

string_vector maintenance_respiration::get_inputs()
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

string_vector maintenance_respiration::get_outputs()
{
    return {
        "Grain",    // Mg / ha
        "Leaf",     // Mg / ha
        "Rhizome",  // Mg / ha
        "Root",     // Mg / ha
        "Shell",    // Mg / ha
        "Stem"      // Mg / ha
    };
}

void maintenance_respiration::do_operation() const
{
    double constexpr Tref = 25.0;  // reference temperature for the Q10 function

    double const dLeaf = -Leaf * mrc_leaf * Q10_temperature_response(temp, Tref);  // Mg / ha

    double const dStem = -Stem * mrc_stem * Q10_temperature_response(temp, Tref);  // Mg / ha

    double const dRoot = -Root * mrc_root * Q10_temperature_response(temp, Tref);  // Mg / ha

    double const dRhizome = -Rhizome * mrc_rhizome * Q10_temperature_response(temp, Tref);  // Mg / ha

    double const dGrain = -Grain * mrc_grain * Q10_temperature_response(temp, Tref);  // Mg / ha

    double const dShell = -Shell * mrc_shell * Q10_temperature_response(temp, Tref);  // Mg / ha

    update(Grain_op, dGrain);      // Mg / ha
    update(Leaf_op, dLeaf);        // Mg / ha
    update(Rhizome_op, dRhizome);  // Mg / ha
    update(Root_op, dRoot);        // Mg / ha
    update(Shell_op, dShell);      // Mg / ha
    update(Stem_op, dStem);        // Mg / ha
}

}  // namespace standardBML
#endif
