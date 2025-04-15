#ifndef MAINTENANCE_RESPIRATION_H
#define MAINTENANCE_RESPIRATION_H

#include "../framework/module.h"
#include "../framework/state_map.h"

namespace standardBML
{
/**
 * @class maintenance_respiration
 *
 * @brief Includes maintenance respiration losses in the rate of change of each
 * organ biomass.
 *
 * This module is intended to be used along with the
 * `maintenance_respiration_calculator` module; see that module for more
 * information.
 */
class maintenance_respiration : public differential_module
{
   public:
    maintenance_respiration(
        state_map const& input_quantities,
        state_map* output_quantities)
        : differential_module{},

          // Get references to input quantities
          Grain_mrr{get_input(input_quantities, "Grain_mrr")},
          Leaf_mrr{get_input(input_quantities, "Leaf_mrr")},
          Rhizome_mrr{get_input(input_quantities, "Rhizome_mrr")},
          Root_mrr{get_input(input_quantities, "Root_mrr")},
          Shell_mrr{get_input(input_quantities, "Shell_mrr")},
          Stem_mrr{get_input(input_quantities, "Stem_mrr")},

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
    const double& Grain_mrr;
    const double& Leaf_mrr;
    const double& Rhizome_mrr;
    const double& Root_mrr;
    const double& Shell_mrr;
    const double& Stem_mrr;

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
        "Grain_mrr",    // Mg / ha / hr
        "Leaf_mrr",     // Mg / ha / hr
        "Rhizome_mrr",  // Mg / ha / hr
        "Root_mrr",     // Mg / ha / hr
        "Shell_mrr",    // Mg / ha / hr
        "Stem_mrr"      // Mg / ha / hr
    };
}

string_vector maintenance_respiration::get_outputs()
{
    return {
        "Grain",    // Mg / ha / hr
        "Leaf",     // Mg / ha / hr
        "Rhizome",  // Mg / ha / hr
        "Root",     // Mg / ha / hr
        "Shell",    // Mg / ha / hr
        "Stem"      // Mg / ha / hr
    };
}

void maintenance_respiration::do_operation() const
{
    update(Grain_op, -Grain_mrr);      // Mg / ha / hr
    update(Leaf_op, -Leaf_mrr);        // Mg / ha / hr
    update(Rhizome_op, -Rhizome_mrr);  // Mg / ha / hr
    update(Root_op, -Root_mrr);        // Mg / ha / hr
    update(Shell_op, -Shell_mrr);      // Mg / ha / hr
    update(Stem_op, -Stem_mrr);        // Mg / ha / hr
}

}  // namespace standardBML
#endif
