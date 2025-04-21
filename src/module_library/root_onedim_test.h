#ifndef ROOT_ONEDIM_TEST_H
#define ROOT_ONEDIM_TEST_H

#include "../framework/module.h"
#include "../framework/state_map.h"
#include "root_onedim.h"

#include <functional>  // for std::function
#include <map>
namespace standardBML
{

namespace root_test_functions
{

// root at 1 for all tests.

struct hard_test {
    double operator()(double x)
    {
        return x < 0 ? -1 : std::pow(x, 10) - 1;
    }

    double derivative(double x)
    {
        return x < 0 ? 0.0 : 10 * std::pow(x, 9);
    }

    double second_derivative(double x)
    {
        return x < 0 ? 0.0 : 90 * std::pow(x, 8);
    }
};

struct linear {
    double operator()(double x)
    {
        return (1 - x);
    }

    double derivative(double x)
    {
        return -1;
    }

    double second_derivative(double x)
    {
        return 0;
    }
};

struct quadratic {
    double operator()(double x)
    {
        return (1 - x) * (x + 1);
    }

    double derivative(double x)
    {
        return -2 * x;
    }

    double second_derivative(double x)
    {
        return -2;
    }
};

struct double_root {
    double operator()(double x)
    {
        return std::pow(x - 1, 2);
    }

    double derivative(double x)
    {
        return 2 * (x - 1);
    }

    double second_derivative(double x)
    {
        return -2;
    }
};

struct triple_root {
    double operator()(double x)
    {
        return std::pow(x - 1, 3);
    }

    double derivative(double x)
    {
        return 3 * std::pow(x - 1, 2);
    }

    double second_derivative(double x)
    {
        return 6 * (x - 1);
    }
};

struct step_function {
    double operator()(double x)
    {
        return x >= 1 ? 1 : -1;
    }

    double derivative(double x)
    {
        return 0;
    }

    double second_derivative(double x)
    {
        return 0;
    }
};

struct singularity {
    double operator()(double x)
    {
        return 1 / (x - 1);
    }

    double derivative(double x)
    {
        return -1 / std::pow(x - 1, 2);
    }

    double second_derivative(double x)
    {
        return 2 / std::pow(x - 1, 3);
    }
};
const std::map<std::string, std::function<double(double)>>
    tests = {
        {"hard_test", hard_test{}},
        {"linear", linear{}},
        {"quadratic", quadratic{}},
        {"double_root", double_root{}},
        {"triple_root", triple_root{}},
        {"step_function", step_function{}},
        {"singularity", singularity{}}};
};  // namespace root_test_functions

class root_onedim_test : public direct_module
{
    struct result {
        double* root_op;
        double* residual_op;
        double* iteration_op;
        double* flag_op;
    };

   public:
    root_onedim_test(state_map const& input_quantities, state_map* output_quantities)
        : direct_module{},

          // Get pointers to input quantities

          // Get pointers to output quantities
          results{get_result_op(output_quantities)}
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "root_onedim_test"; }

   private:
    // Pointers to input quantities

    // Pointers to output quantities
    std::map<std::string, result> results;

    // Main operation
    void do_operation() const;

    std::map<std::string, result> get_result_op(state_map* output_quantities)
    {
        std::map<std::string, result> out;
        for (auto& test : root_test_functions::tests) {
            auto name = test.first;
            out[name] = result{
                get_op(output_quantities, name + "_root"),
                get_op(output_quantities, name + "_residual"),
                get_op(output_quantities, name + "_iteration"),
                get_op(output_quantities, name + "_flag")};
        }
        return out;
    }
};

string_vector root_onedim_test::get_inputs()
{
    return {
        // None
    };
}

string_vector root_onedim_test::get_outputs()
{
    string_vector out;
    for (auto& test : root_test_functions::tests) {
        auto name = test.first;
        out.push_back(name + "_root");
        out.push_back(name + "_residual");
        out.push_back(name + "_iteration");
        out.push_back(name + "_flag");
    }
    return out;
}

void root_onedim_test::do_operation() const
{
    // Collect inputs and make calculations
    root_algorithm::root_finder<root_algorithm::bisection> solver(100, 1e-14, 1e-12);
    root_algorithm::result_t result;

    for (auto& test : root_test_functions::tests) {
        auto& name = test.first;
        auto& func = test.second;
        result = solver.solve(func, 0, 3);

        auto& r = results.at(name);

        update(r.root_op, result.root);
        update(r.residual_op, result.residual);
        update(r.iteration_op, result.iteration);
        update(r.flag_op, static_cast<int>(result.flag));
    }

    // Update the output quantity list
}

}  // namespace standardBML
#endif
