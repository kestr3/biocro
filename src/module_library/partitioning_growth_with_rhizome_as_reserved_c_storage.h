#ifndef PARTITIONING_GROWTH_WITH_RHIZOME_AS_RESERVED_C_STORAGE_H
#define PARTITIONING_GROWTH_WITH_RHIZOME_AS_RESERVED_C_STORAGE_H

#include "../framework/module.h"
#include "../framework/state_map.h"
#include "BioCro.h" // include for resp()

// This is the same as partitioning_growth_module/ no leaf_resp_partitioning_growth_calculator except that 
// 1. CanopyA already includes losses from leaf respiration, so it should not be removed again from leaf.
// 2. When canopy assimilation rate is below zero, R

namespace standardBML 
{

class partitioning_growth_with_rhizome_as_reserved_c_storage : public direct_module
{    
    public:
        partitioning_growth_with_rhizome_as_reserved_c_storage(
            //const std::unordered_map<std::string, double>* input_parameters, std::unordered_map<std::string, double>* output_parameters) :
            state_map const& input_quantities,
            state_map* output_quantities)
            : direct_module{},

            // Get pointers to input parameters
            kLeaf{get_input(input_quantities, "kLeaf")},
            kStem{get_input(input_quantities, "kStem")},
            kRoot{get_input(input_quantities, "kRoot")},
            kRhizome{get_input(input_quantities, "kRhizome")},
            canopy_assimilation_rate{get_input(input_quantities, "canopy_assimilation_rate")},
            mrc1{get_input(input_quantities, "mrc1")},
            mrc2{get_input(input_quantities, "mrc2")},
            temp{get_input(input_quantities, "temp")},

            // Get pointers to output parameters
            net_assimilation_rate_leaf_op{get_op(output_quantities, "net_assimilation_rate_leaf")},
            net_assimilation_rate_stem_op{get_op(output_quantities, "net_assimilation_rate_stem")},
            net_assimilation_rate_root_op{get_op(output_quantities, "net_assimilation_rate_root")},
            net_assimilation_rate_rhizome_op{get_op(output_quantities, "net_assimilation_rate_rhizome")}
        {
        }
        static string_vector get_inputs();
        static string_vector get_outputs();
        static std::string get_name() {return "partitioning_coefficient_logistic"; }

    private:
        // Pointers to input parameters
        const double& kLeaf;
        const double& kStem;
        const double& kRoot;
        const double& kRhizome;
        const double& canopy_assimilation_rate;
        const double& mrc1;
        const double& mrc2;
        const double& temp;
        // Pointers to output parameters
        double* net_assimilation_rate_leaf_op;
        double* net_assimilation_rate_stem_op;
        double* net_assimilation_rate_root_op;
        double* net_assimilation_rate_rhizome_op;
        // Main operation
        void do_operation() const override final;
};

string_vector partitioning_growth_with_rhizome_as_reserved_c_storage::get_inputs()
{   
    return {
        "kLeaf",
        "kStem",
        "kRoot",
        "kRhizome",
        "canopy_assimilation_rate",
        "mrc1",
        "mrc2",
        "temp"
    };
}

string_vector partitioning_growth_with_rhizome_as_reserved_c_storage::get_outputs()
{
    return {
        "net_assimilation_rate_leaf",     // Mg / ha / hour
        "net_assimilation_rate_stem",     // Mg / ha / hour
        "net_assimilation_rate_root",     // Mg / ha / hour
        "net_assimilation_rate_rhizome",  // Mg / ha / hour
    };
}

void partitioning_growth_with_rhizome_as_reserved_c_storage::do_operation() const 
{
    // Collect inputs and make calculations
    
    // double kLeaf = kLeaf_ip;
    // double kStem = kStem_ip;
    // double kRoot = kRoot_ip;
    // double kRhizome = kRhizome_ip;
    // double canopy_assimilation_rate = canopy_assimilation_rate_ip;
    // double mrc1 = mrc1_ip;
    // double mrc2 = mrc2_ip;
    // double temp = temp_ip;
    

    // this is where im leaving off 10/2/25
    double net_assimilation_rate_leaf;
    double net_assimilation_rate_stem; 
    double net_assimilation_rate_root;
    double net_assimilation_rate_rhizome;
    //double net_assimilation_rate_grain;
    
    double nonrhizome_carbon_flux; // nonleaf_carbon_flux;
    if(canopy_assimilation_rate < 0) nonrhizome_carbon_flux = 0.0;
    else nonrhizome_carbon_flux = canopy_assimilation_rate;
    
    // Calculate the amount of new leaf produced
    if(kLeaf > 0) {
        net_assimilation_rate_leaf = nonrhizome_carbon_flux * kLeaf;
 //     newLeafcol = resp(newLeafcol, mrc1, temp);// leaf respiration is also included in canopy assimlation, so no need to add it here again.
    }

    
    // Calculate the amount of new stem produced
    if(kStem >= 0) {
        net_assimilation_rate_stem = nonrhizome_carbon_flux * kStem;
        net_assimilation_rate_stem = resp(net_assimilation_rate_stem, mrc1, temp);
    }
//    else throw std::range_error("Thrown in partitioning_growth_with_rhizome_as_reserved_c_storage: kStem should be positive"); MLM removed this error message 04/22/2020; can cause issues with integration
    
    // Calculate the amount of new root produced
    if(kRoot > 0) {
        net_assimilation_rate_root = nonrhizome_carbon_flux * kRoot;
        net_assimilation_rate_root = resp(net_assimilation_rate_root, mrc2, temp);
    }
    else net_assimilation_rate_root = 0.0;
    
    // Calculate the amount of new rhizome produced
    if(kRhizome > 0) {
        if(canopy_assimilation_rate < 0) 
         {
          net_assimilation_rate_rhizome = canopy_assimilation_rate;// Negative assimilation rate means reduction in rhizome biomass
          net_assimilation_rate_rhizome = resp(0.0, mrc2, temp); // when canopy assimilation is below zero, there is no growth in rhizome mass so respiration is calculated assuming newrhizomecol = 0
         }
         else 
        {
          net_assimilation_rate_rhizome = canopy_assimilation_rate * kRhizome;
          net_assimilation_rate_rhizome = resp(net_assimilation_rate_rhizome, mrc2, temp);
         }
    }
    else net_assimilation_rate_rhizome = 0.0;
    
    // Grain has no respiration or senescence at the moment, so we don't need to calculate
    //  the amount of new grain here
    
    // Update the output parameter list
    update(net_assimilation_rate_leaf_op, net_assimilation_rate_leaf);
    update(net_assimilation_rate_stem_op, net_assimilation_rate_stem);
    update(net_assimilation_rate_root_op, net_assimilation_rate_root);
    update(net_assimilation_rate_rhizome_op, net_assimilation_rate_rhizome);
}

// double resp(double base_rate, double grc, double temp)
// {
//     double ans = base_rate * (1 - (grc * Q10_temperature_response(temp, 0.0)));

//     if (ans < 0) ans = 0;

//     return ans;
// }

#endif
}