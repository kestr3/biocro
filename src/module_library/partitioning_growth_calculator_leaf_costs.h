#ifndef PARTITIONING_GROWTH_CALCULATOR_LEAF_COSTS_H
#define PARTITIONING_GROWTH_CALCULATOR_LEAF_COSTS_H

#include "../framework/module.h"
#include "../framework/state_map.h"
#include "respiration.h"  // for growth_resp_Q10

namespace standardBML
{
/**
 *  @class partitioning_growth_calculator_leaf_costs
 *
 *  @brief Uses a set of partitioning coefficients to determine net assimilation
 *  rates due to photosynthesis and respiration for several plant organs.
 *
 *  ### Partitioning overview
 *
 *  BioCro includes several partitioning growth calculators that determine these
 *  rates using slightly different methods. The different modules can be
 *  distinguished by the sets of tissues they use, the ways they apply
 *  respiration and water stress, and their responses to negative canopy
 *  assimilation rates. (A negative canopy assimilation rate indicates that the
 *  leaves are respiring more than they are photosynthesizing.)
 *
 *  In all partitioning growth calculators, the base growth rate for an organ is
 *  determined from the net canopy assimilation rate and a coefficient that
 *  determines the fraction of the net assimilation that is "partitioned" to
 *  that organ. Then, further modifications may take place to account for water
 *  stress, maintenance respiration, or other processes that affect the amount
 *  of carbon available to the organ for growth. Note that losses due to
 *  senescence and gains due to remobilized carbon from other organs are handled
 *  elsewhere and are not included here.
 *
 *  Respiration is included via the `growth_resp_Q10()` function, which
 *  implements an empirical rule for determining the fraction of energy spent on
 *  respiration at a particular temperature. See the following paper for a
 *  general discussion of the importance of respiration in understanding plant
 *  growth: [Amthor, J. S. "The role of maintenance respiration in plant growth"
 *  Plant, Cell & Environment 7, 561–569 (1984)]
 *  (https://doi.org/10.1111/1365-3040.ep11591833).
 *
 *  ### Specifics of this module
 *
 *  When the canopy assimilation rate is negative, this module ensures that all
 *  losses come from the leaf, and the growth rates for the other organs are set
 *  to zero since no carbon is flowing out of the leaf.
 *
 *  This module does not attempt to explicitly include any effect due to water
 *  stress.
 *
 *  This module includes five organs:
 *  - `Leaf`: The leaf growth rate is *not* modified by respiration because the
 *     net canopy assimilation rate already includes it.
 *  - `Stem`: The stem growth rate is modified by respiration.
 *  - `Root`: The root growth rate is modified by respiration.
 *  - `Rhizome`: The rhizome growth rate is modified by respiration.
 *  - `Grain`: The grain growth rate is *not* modified by respiration.
 *  - `Shell`: The shell growth rate is *not* modified by respiration.
 *
 *  Along with the growth rate of each tissue, this module also calculates
 *  growth respiration rates; the associated quantity names end with `_gr_rate`.
 */
class partitioning_growth_calculator_leaf_costs : public direct_module
{
   public:
    partitioning_growth_calculator_leaf_costs(
        state_map const& input_quantities,
        state_map* output_quantities)
        : direct_module{},

          // Get references to input quantities
          canopy_assim{get_input(input_quantities, "canopy_assimilation_rate")},
          grc_rhizome{get_input(input_quantities, "grc_rhizome")},
          grc_root{get_input(input_quantities, "grc_root")},
          grc_stem{get_input(input_quantities, "grc_stem")},
          kGrain{get_input(input_quantities, "kGrain")},
          kLeaf{get_input(input_quantities, "kLeaf")},
          kRhizome{get_input(input_quantities, "kRhizome")},
          kRoot{get_input(input_quantities, "kRoot")},
          kShell{get_input(input_quantities, "kShell")},
          kStem{get_input(input_quantities, "kStem")},
          temp{get_input(input_quantities, "temp")},

          // Get pointers to output quantities
          Grain_gr_rate_op{get_op(output_quantities, "Grain_gr_rate")},
          Leaf_gr_rate_op{get_op(output_quantities, "Leaf_gr_rate")},
          Leaf_WS_loss_rate_op{get_op(output_quantities, "Leaf_WS_loss_rate")},
          net_assimilation_rate_grain_op{get_op(output_quantities, "net_assimilation_rate_grain")},
          net_assimilation_rate_leaf_op{get_op(output_quantities, "net_assimilation_rate_leaf")},
          net_assimilation_rate_rhizome_op{get_op(output_quantities, "net_assimilation_rate_rhizome")},
          net_assimilation_rate_root_op{get_op(output_quantities, "net_assimilation_rate_root")},
          net_assimilation_rate_shell_op{get_op(output_quantities, "net_assimilation_rate_shell")},
          net_assimilation_rate_stem_op{get_op(output_quantities, "net_assimilation_rate_stem")},
          Rhizome_gr_rate_op{get_op(output_quantities, "Rhizome_gr_rate")},
          Root_gr_rate_op{get_op(output_quantities, "Root_gr_rate")},
          Shell_gr_rate_op{get_op(output_quantities, "Shell_gr_rate")},
          Stem_gr_rate_op{get_op(output_quantities, "Stem_gr_rate")}
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "partitioning_growth_calculator_leaf_costs"; }

   private:
    // Pointers to input quantities
    // References to input quantities
    const double& canopy_assim;
    const double& grc_rhizome;
    const double& grc_root;
    const double& grc_stem;
    const double& kGrain;
    const double& kLeaf;
    const double& kRhizome;
    const double& kRoot;
    const double& kShell;
    const double& kStem;
    const double& temp;

    // Pointers to output quantities
    double* Grain_gr_rate_op;
    double* Leaf_gr_rate_op;
    double* Leaf_WS_loss_rate_op;
    double* net_assimilation_rate_grain_op;
    double* net_assimilation_rate_leaf_op;
    double* net_assimilation_rate_rhizome_op;
    double* net_assimilation_rate_root_op;
    double* net_assimilation_rate_shell_op;
    double* net_assimilation_rate_stem_op;
    double* Rhizome_gr_rate_op;
    double* Root_gr_rate_op;
    double* Shell_gr_rate_op;
    double* Stem_gr_rate_op;

    // Main operation
    void do_operation() const;
};

string_vector partitioning_growth_calculator_leaf_costs::get_inputs()
{
    return {
        "canopy_assimilation_rate",  // Mg / ha / hour
        "grc_rhizome",               // dimensionless
        "grc_root",                  // dimensionless
        "grc_stem",                  // dimensionless
        "kGrain",                    // dimensionless
        "kLeaf",                     // dimensionless
        "kRhizome",                  // dimensionless
        "kRoot",                     // dimensionless
        "kShell",                    // dimensionless
        "kStem",                     // dimensionless
        "temp"                       // degrees C
    };
}

string_vector partitioning_growth_calculator_leaf_costs::get_outputs()
{
    return {
        "Grain_gr_rate",                  // Mg / ha / hour
        "Leaf_gr_rate",                   // Mg / ha / hour
        "Leaf_WS_loss_rate",              // Mg / ha / hour
        "net_assimilation_rate_grain",    // Mg / ha / hour
        "net_assimilation_rate_leaf",     // Mg / ha / hour
        "net_assimilation_rate_rhizome",  // Mg / ha / hour
        "net_assimilation_rate_root",     // Mg / ha / hour
        "net_assimilation_rate_shell",    // Mg / ha / hour
        "net_assimilation_rate_stem",     // Mg / ha / hour
        "Rhizome_gr_rate",                // Mg / ha / hour
        "Root_gr_rate",                   // Mg / ha / hour
        "Shell_gr_rate",                  // Mg / ha / hour
        "Stem_gr_rate"                    // Mg / ha / hour
    };
}

void partitioning_growth_calculator_leaf_costs::do_operation() const
{
    // Specify the base temperature for Q10 growth respiration
    double constexpr Tref = 0.0;  // degrees C

    // Calculate the rate of new leaf production, including _all_ losses due to
    // photorespiration and RL, but not accounting for water stress or any
    // additional growth respiration (Mg / ha / hr)
    double const base_rate_leaf{
        kLeaf <= 0         ? 0
        : canopy_assim < 0 ? canopy_assim
                           : canopy_assim * kLeaf};
    double constexpr Leaf_WS_loss_rate{0.0};
    double constexpr Leaf_gr_rate{0.0};

    // Determine the carbon flux to use for the non-leaf organs (Mg / ha / hr)
    double const nonleaf_carbon_flux{canopy_assim > 0 ? canopy_assim : 0};

    // Calculate the base rate of new stem production and the associated
    // respiratory costs (Mg / ha / hr)
    double const base_rate_stem{kStem > 0 ? nonleaf_carbon_flux * kStem : 0};
    double const Stem_gr_rate{growth_resp_Q10(base_rate_stem, grc_stem, temp, Tref)};

    // Calculate the base rate of new root production and the associated
    // respiratory costs (Mg / ha / hr)
    double const base_rate_root{kRoot > 0 ? nonleaf_carbon_flux * kRoot : 0};
    double const Root_gr_rate{growth_resp_Q10(base_rate_root, grc_root, temp, Tref)};

    // Calculate the base rate of new rhizome production and the associated
    // respiratory costs (Mg / ha / hr)
    double const base_rate_rhizome{kRhizome > 0 ? nonleaf_carbon_flux * kRhizome : 0};
    double const Rhizome_gr_rate{growth_resp_Q10(base_rate_rhizome, grc_rhizome, temp, Tref)};

    // Calculate the base rate of new grain production (Mg / ha / hr)
    double const base_rate_grain{kGrain > 0 ? nonleaf_carbon_flux * kGrain : 0};
    double const Grain_gr_rate{0.0};

    // Calculate the base rate of new shell production (Mg / ha / hr)
    double const base_rate_shell{kShell > 0 ? nonleaf_carbon_flux * kShell : 0};
    double const Shell_gr_rate{0.0};

    // Update the output quantity list
    update(Grain_gr_rate_op, Grain_gr_rate);
    update(Leaf_gr_rate_op, Leaf_gr_rate);
    update(Leaf_WS_loss_rate_op, Leaf_WS_loss_rate);
    update(net_assimilation_rate_grain_op, base_rate_grain - Grain_gr_rate);
    update(net_assimilation_rate_leaf_op, base_rate_leaf - Leaf_gr_rate);
    update(net_assimilation_rate_rhizome_op, base_rate_rhizome - Rhizome_gr_rate);
    update(net_assimilation_rate_root_op, base_rate_root - Root_gr_rate);
    update(net_assimilation_rate_shell_op, base_rate_shell - Shell_gr_rate);
    update(net_assimilation_rate_stem_op, base_rate_stem - Stem_gr_rate);
    update(Rhizome_gr_rate_op, Rhizome_gr_rate);
    update(Root_gr_rate_op, Root_gr_rate);
    update(Shell_gr_rate_op, Shell_gr_rate);
    update(Stem_gr_rate_op, Stem_gr_rate);
}

}  // namespace standardBML
#endif
