#ifndef BROYDEN_TEST_H
#define BROYDEN_TEST_H

#include "../framework/module.h"
#include "../framework/state_map.h"
#include "broyden_method.h"

#include <array>

namespace standardBML
{

std::array<double, 2> test_function(std::array<double, 2> x)
{
    std::array<double, 2> y;
    y[0] = x[0] / (x[0] * x[0] + 1) - x[1] + 0.5;
    y[1] = 2 - x[0] - x[1];
    return y;
}

class broyden_test : public direct_module
{
    struct result {
        double* root_op;
        double* residual_op;
        double* iteration_op;
        double* flag_op;

        result(
            state_map* output_quantities,
            std::string&& name)
            : root_op{get_op(output_quantities, name + "_root")},
              residual_op{get_op(output_quantities, name + "_residual")},
              iteration_op{get_op(output_quantities, name + "_iteration")},
              flag_op{get_op(output_quantities, name + "_flag")}
        {
        }
    };

   public:
    broyden_test(
        state_map const& input_quantities, state_map* output_quantities)
        : direct_module{},

          // Get pointers to input quantities
          max_iterations{get_input(input_quantities, "max_iterations")},
          abs_tol{get_input(input_quantities, "abs_tol")},
          rel_tol{get_input(input_quantities, "rel_tol")},

    // Get pointers to output quantities

    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "broyden_test"; }

   private:
    // Pointers to input quantities
    const double& max_iterations;
    const double& abs_tol;
    const double& rel_tol;

    // Pointers to output quantities

    // Main operation
    void do_operation() const;

    static string_vector make_qname(std::string& name)
    {
        return {
            name + "_root",
            name + "_residual",
            name + "_iteration",
            name + "_flag"};
    }

    void inline update_result(
        const result& r, root_algorithm::result_t& result) const
    {
        update(r.root_op, result.root);
        update(r.residual_op, result.residual);
        update(r.iteration_op, result.iteration);
        update(r.flag_op, static_cast<int>(result.flag));
    }
};

string_vector broyden_test::get_inputs()
{
    return {
        "ecc",
        "answer",
        "max_iterations",
        "abs_tol",
        "rel_tol",
        "lower_bracket",
        "upper_bracket",
        "single_guess"};
}

string_vector broyden_test::get_outputs()
{
    string_vector out;
    const string_vector methods = {
        "secant", "fixed_point", "newton", "halley", "steffensen",
        "bisection", "regula_falsi", "ridder",
        "illinois", "pegasus", "anderson_bjorck"};
    for (auto name : methods) {
        string_vector sv = make_qname(name);
        out.insert(out.end(), sv.begin(), sv.end());
    }

    return out;
}

void broyden_test::do_operation() const
{
}

}  // namespace standardBML
#endif
