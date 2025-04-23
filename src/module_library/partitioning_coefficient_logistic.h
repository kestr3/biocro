#ifndef PARTITIONING_COEFFICIENT_LOGISTIC_H
#define PARTITIONING_COEFFICIENT_LOGISTIC_H

#include "../framework/module.h"
#include "../framework/state_map.h"
#include <cmath>  // for exp

namespace standardBML
{
double strength_term(double const alpha, double const beta, double const DVI);

/**
 * @class partitioning_coefficient_logistic
 *
 * @brief Calculates carbon partitioning coefficients based on logistic-based
 * functions and development index using the logistic-based functions from
 * Osborne et al 2015.
 *
 * Intended to be used with any of the following modules:
 * - `partitioning_growth_calculator_leaf_costs`
 * - `partitioning_growth_calculator`
 *
 * Using the following function, calculates the percentage of carbon allocated
 * to the root, stem, leaf, shell, and grain at a given development index.
 *
 * \f[ k_i = \frac{\exp{(\alpha_i+\beta_i x)}}  {\exp{(\alpha_R+\beta_R x)} +
 * \exp{(\alpha_S+\beta_S x)} + \exp{(\alpha_L+\beta_L x)} +
 * \exp{(\alpha_{Sh}+\beta_{Sh} x)} + 1}, \f]
 *
 * where \f$ i = {R, S, L, Sh} \f$ for root, stem, leaf, and shell respectively,
 * and \f$ x \f$ is the development index. For the grain,
 *
 * \f[ k_G = \frac{1}{\exp{(\alpha_R+\beta_R x)} + \exp{(\alpha_S+\beta_S x)} +
 * \exp{(\alpha_L+\beta_L x)} + \exp{(\alpha_{Sh}+\beta_{Sh})} + 1}. \f]
 *
 * See Matthews et al. for more description of how this module was used in
 * Soybean-BioCro and for details on the parameter fitting to identify the
 * \f$ \alpha \text{ and } \beta \f$ parameters. Note that the original model
 * did not include a shell component.
 *
 * Although it is not used in the soybean model, this module also includes an
 * option for a rhizome to contribute carbon to other organs during emergence.
 * See comments in the code for more details.
 *
 * ### References:
 *
 * [Matthews, M. L.et al. 2021. "Soybean-BioCro: a semi-mechanistic model of
 * soybean growth." in silico Plants 4, diab032.]
 * (https://doi.org/10.1093/insilicoplants/diab032)
 *
 * [Osborne, T. et al. 2015. "JULES-Crop: A Parametrisation of Crops in the Joint
 * UK Land Environment Simulator." Geoscientific Model Development 8(4): 1139–55.]
 * (https://doi.org/10.5194/gmd-8-1139-2015)
 */
class partitioning_coefficient_logistic : public direct_module
{
   public:
    partitioning_coefficient_logistic(
        state_map const& input_quantities,
        state_map* output_quantities)
        : direct_module{},

          // Get references to input quantities
          alphaLeaf{get_input(input_quantities, "alphaLeaf")},
          alphaRoot{get_input(input_quantities, "alphaRoot")},
          alphaShell{get_input(input_quantities, "alphaShell")},
          alphaStem{get_input(input_quantities, "alphaStem")},
          betaLeaf{get_input(input_quantities, "betaLeaf")},
          betaRoot{get_input(input_quantities, "betaRoot")},
          betaShell{get_input(input_quantities, "betaShell")},
          betaStem{get_input(input_quantities, "betaStem")},
          DVI{get_input(input_quantities, "DVI")},
          kRhizome_emr{get_input(input_quantities, "kRhizome_emr")},

          // Get pointers to output quantities
          kGrain_op{get_op(output_quantities, "kGrain")},
          kLeaf_op{get_op(output_quantities, "kLeaf")},
          kRhizome_op{get_op(output_quantities, "kRhizome")},
          kRoot_op{get_op(output_quantities, "kRoot")},
          kShell_op{get_op(output_quantities, "kShell")},
          kStem_op{get_op(output_quantities, "kStem")}
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "partitioning_coefficient_logistic"; }

   private:
    // Pointers to input quantities
    const double& alphaLeaf;
    const double& alphaRoot;
    const double& alphaShell;
    const double& alphaStem;
    const double& betaLeaf;
    const double& betaRoot;
    const double& betaShell;
    const double& betaStem;
    const double& DVI;
    const double& kRhizome_emr;

    // Pointers to output quantities
    double* kGrain_op;
    double* kLeaf_op;
    double* kRhizome_op;
    double* kRoot_op;
    double* kShell_op;
    double* kStem_op;

    // Implement the pure virtual function do_operation():
    void do_operation() const override final;
};

string_vector partitioning_coefficient_logistic::get_inputs()
{
    return {
        "alphaLeaf",    // dimensionless
        "alphaRoot",    // dimensionless
        "alphaShell",   // dimensionless
        "alphaStem",    // dimensionless
        "betaLeaf",     // dimensionless
        "betaRoot",     // dimensionless
        "betaShell",    // dimensionless
        "betaStem",     // dimensionless
        "DVI",          // dimensionless
        "kRhizome_emr"  // dimensionless
    };
}

string_vector partitioning_coefficient_logistic::get_outputs()
{
    return {
        "kGrain",    // dimensionless
        "kLeaf",     // dimesnionless
        "kRhizome",  // dimensionless
        "kRoot",     // dimensionless
        "kShell",    // dimensionless
        "kStem"      // dimensionless
    };
}

void partitioning_coefficient_logistic::do_operation() const
{
    // Determine partitioning coefficients using multinomial logistic equations
    // from Osborne et al., 2015 JULES-crop https://doi.org/10.5194/gmd-8-1139-2015

    // Calculate the sink strength of each tissue (relative to grain)
    double const leaf_strength = strength_term(alphaLeaf, betaLeaf, DVI);
    double const root_strength = strength_term(alphaRoot, betaRoot, DVI);
    double const shell_strength = strength_term(alphaShell, betaShell, DVI);
    double const stem_strength = strength_term(alphaStem, betaStem, DVI);
    double constexpr grain_strength = 1.0;
    double constexpr rhizome_strength = 0.0;

    // Calculate the total sink strength
    double const total_strength =
        leaf_strength + rhizome_strength + root_strength + shell_strength +
        stem_strength + grain_strength;

    // The k values are the fraction of total demand from each tissue
    double const kGrain{grain_strength / total_strength};                               // dimensionless
    double const kLeaf{leaf_strength / total_strength};                                 // dimensionless
    double const kRhizome{DVI < 0 ? kRhizome_emr : rhizome_strength / total_strength};  // dimensionless
    double const kRoot{root_strength / total_strength};                                 // dimensionless
    double const kShell{shell_strength / total_strength};                               // dimensionless
    double const kStem{stem_strength / total_strength};                                 // dimensionless

    // Update the output quantities
    update(kGrain_op, kGrain);      // dimensionless
    update(kLeaf_op, kLeaf);        // dimensionless
    update(kRhizome_op, kRhizome);  // dimensionless
    update(kRoot_op, kRoot);        // dimensionless
    update(kShell_op, kShell);      // dimensionless
    update(kStem_op, kStem);        // dimensionless
}

double strength_term(double const alpha, double const beta, double const DVI)
{
    return exp(alpha + beta * DVI);  // dimensionless
}

}  // namespace standardBML
#endif
