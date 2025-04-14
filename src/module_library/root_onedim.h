#ifndef ROOTS_ONEDIM_H
#define ROOTS_ONEDIM_H

#include <cmath>
#include <limits>
#include <string>

namespace root_algorithm
{

// Helper function declarations
inline bool is_close(double x, double y, double tol, double rtol);
inline bool is_zero(double x, double tol);
inline bool same_signs(double x, double y);
inline bool opposite_signs(double x, double y);
inline bool smaller(double x, double y);


struct graph_t {
    double x;
    double y;
};

enum class Flag {
    valid,  // don't terminate
    residual_zero,
    delta_root_zero,
    bracket_width_zero,
    max_iterations,
    invalid_bracket,
    division_by_zero
};

bool successful_termination(Flag flag);
std::string flag_message(Flag flag);


template <typename T>
struct with_flag {
    T state;
    Flag flag;

    with_flag() {}
    with_flag(T x) : state(x), flag(Flag::valid) {}
    with_flag(T x, Flag f) : state(x), flag(f) {}
    with_flag(Flag f) : flag(f) {}


};

template <typename T>
inline bool is_valid(with_flag<T> x){ return x.flag == Flag::valid; }


/**
 * @class result_t
 *
 * @brief Result from root finding algorithm.
 *
 * @param root The identified root.
 *
 * @param [in] residual The output value (root); should equal zero if
 * successful.
 *
 * @param [in] iteration The number of iterations performed.
 *
 * @return result_t
 *
 * @details
 *
 *
 */
struct result_t {
    double root;
    double residual;
    size_t iteration;
    Flag flag;
};

/**
 * @class root_finder
 *
 * @brief Function object for finding the root or zero of a function.
 *
 * @param [in] max_iter The total number of iterations allowed.
 *
 * @param [in] abs_tol The absolute tolerance, the error threshold for
 * floating-point zero. E.g., abs_tol = 1e-12 means |x| < 1e-12, or x equals
 * zero to 12 digits. Used for convergence testing.
 *
 * @param [in] rel_tol The relative tolerance, the error threshold for
 * comparing floating point number equality. E.g., rel_tol = 1e-12 means
 * x equals y if the first 12 digits match. Used for convergernce testing.
 *
 * @return result_t
 *
 * @details This class creates a function object whose call method (operator())
 * is an interface for using any of the root finding methods. A class template
 * so an instance of this class must be declared with a method as the template
 * argument. Example usage:
 *
 * @code{.cpp}
 * // Zero is the sqrt of 3
 * double f(double x) { return x*x - 3; };
 *
 * //Declare the solver
 * root_algorithm::root_finder<root_algorithm::secant> solver;
 *
 * // Call the solve method to find a root.
 * root_algorithm::result_t result = solver.solve(f, 1., 2.);
 * @endcode
 *
 * The function can be a c++ function, lambda function, or function object.
 * However, for algorithms which explicitly evaluate derivatives, you must
 * pass a function object with the appropriate derivative method.
 *
 */
template <typename Method>
struct root_finder : public Method {

    size_t max_iterations = 100;
    double _abs_tol = std::numeric_limits<double>::epsilon();
    double _rel_tol = 2 * std::numeric_limits<double>::epsilon();

    root_finder() {}
    root_finder(size_t max_iter, double abs_tol, double rel_tol) : max_iterations{max_iter}, _abs_tol{abs_tol}, _rel_tol{rel_tol} {}

    using flagged_state = with_flag<typename Method::state>;

    template <typename F, typename... Args>
    result_t solve(F&& func, Args&&... args)
    {
        flagged_state s = Method::initialize(
            std::forward<F>(func), std::forward<Args>(args)...);

        for (size_t i = 0; i < (max_iterations + 1); ++i) {
            if (is_valid(s)) {
                s.flag = Method::check_convergence(s.state, _abs_tol, _rel_tol);
            }

            if (s.flag != Flag::valid) {
                return make_result(s, i);
            }

            s = Method::iterate(std::forward<F>(func), s.state);
        }

        s.flag = Flag::max_iterations;
        return make_result(s, max_iterations);
    }

    private:

    inline result_t make_result(const flagged_state& s, size_t iteration)
    {
        return result_t{Method::root(s.state), Method::residual(s.state), iteration, s.flag};
    }

};

// Methods
// Multi-step methods

struct two_step_method {

    struct state {
        graph_t last;
        graph_t current;
    };

    template <typename F>
    with_flag<state> initialize(F&& fun, double x0, double x1)
    {
        graph_t first{x0, fun(x0)};
        graph_t second{x1, fun(x1)};
        return state{first, second};
    }


    Flag check_convergence(const state& s,  double abs_tol, double rel_tol){


        if (is_zero(s.current.y, abs_tol)) {
            return Flag::residual_zero;
        }

        if (is_close(s.last.x, s.current.x, abs_tol, rel_tol)) {
            return Flag::delta_root_zero;
        }

        return Flag::valid;
    }

    inline double root(const state& s)
    {
        return s.current.x;
    }

    inline double residual(const state& s)
    {
        return s.current.y;
    }

};

struct secant : public two_step_method {


    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        double secant_slope = (s.current.y - s.last.y) / (s.current.x - s.last.x);
        if (secant_slope == 0.0) {
            return {s, Flag::division_by_zero};
        }
        double x2 = s.current.x - s.current.y / secant_slope;
        graph_t new_guess{x2, fun(x2)};
        return state{s.current, new_guess};
    }

};

// Householder methods

struct newton : public two_step_method {


    template <typename F>
    with_flag<state> initialize(F&& fun, double x0)
    {
        state s;
        s.current = {x0, fun(x0)};
        return iterate(std::forward<F>(fun), s);
    }

    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        double slope = fun.derivative(s.current.x);
        if (slope == 0.0) { return {s, Flag::division_by_zero}; }
        double x_new = s.current.x - s.current.y /slope;
        graph_t new_guess{x_new, fun(x_new)};
        return state{s.current, new_guess};
    }

};

struct steffensen : public two_step_method {


    template <typename F>
    with_flag<state> initialize(F&& fun, double x0)
    {
        state s;
        s.current = {x0, fun(x0)};
        return iterate(std::forward<F>(fun), s);
    }

    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        double gx = fun(s.current.x + s.current.y) / s.current.y - 1;
        if (gx == 0.0) { return {s, Flag::division_by_zero}; }
        double x_new = s.current.x - s.current.y / gx;
        graph_t new_guess{x_new, fun(x_new)};
        return state{s.current, new_guess};
    }

};

// Bracketing methods
struct bracket_method {

    struct state {
        graph_t left;
        graph_t right;
    };

    template <typename F>
    with_flag<state> initialize(F&& fun, double a, double b)
    {
        graph_t left{a, fun(a)};
        graph_t right{b, fun(b)};
        if (same_signs(left.y, right.y)){
            return {state{left, right}, Flag::invalid_bracket};
        }
        return state{left, right};
    }

    Flag check_convergence(const state& s,  double abs_tol, double rel_tol){
        bool a = is_zero(s.left.y, abs_tol);
        bool b = is_zero(s.right.y, abs_tol);
        if (a || b) return Flag::residual_zero;

        bool c = is_close(s.left.x, s.right.x, abs_tol, rel_tol);
        if(c) return Flag::bracket_width_zero;

        return Flag::valid;
    }

};

struct bisection : public bracket_method {


    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        graph_t midpoint;
        midpoint.x = (s.left.x + s.right.x) / 2;
        midpoint.y = fun(midpoint.x);

        if (same_signs(s.left.y, midpoint.y)) {
            return state{midpoint, s.right};
        } else {
            return state{s.left, midpoint};
        }
    }

    inline double root(const state& s)
    {
        return (s.left.x + s.right.x) / 2;
    }

    inline double residual(const state& s)
    {
        return (s.left.y + s.right.y) / 2;
    }

};

struct regula_falsi : public bracket_method {

    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        graph_t new_point;
        double slope = (s.right.y - s.left.y) / (s.right.x - s.left.x);
        if (slope == 0.0) { return {s , Flag::division_by_zero}; }
        new_point.x = s.right.x - s.right.y / slope;
        new_point.y = fun(new_point.x);

        if (same_signs(s.left.y, new_point.y)) {
            return state{new_point, s.right};
        }
        return state{s.left, new_point};
    }

    inline double root(const state& s)
    {
        return smaller(s.left.y, s.right.y) ? s.left.x : s.right.x;
    }

    inline double residual(const state& s)
    {
        return smaller(s.left.y, s.right.y) ? s.left.y : s.right.y;
    }

};

struct ridder : public bracket_method {

    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        graph_t midpoint;
        midpoint.x = (s.left.x + s.right.x) / 2;
        midpoint.y = fun(midpoint.x);

        double d = midpoint.x - s.left.x;
        double a = midpoint.y / s.left.y;
        double b = s.right.y / s.left.y;

        graph_t new_point;
        new_point.x = midpoint.x + d * a / std::sqrt(a * a - b);
        if (!std::isfinite(new_point.x )) {
            return {s, Flag::division_by_zero};
        }
        new_point.y = fun(new_point.x);

        if (same_signs(s.left.y, new_point.y)) {
            return state{new_point, s.right};
        }
        return state{s.left, new_point};
    }

    inline double root(const state& s)
    {
        return smaller(s.left.y, s.right.y) ? s.left.x : s.right.x;
    }

    inline double residual(const state& s)
    {
        return smaller(s.left.y, s.right.y) ? s.left.y : s.right.y;
    }

};


// Helper function definitions.
inline bool is_close(double x, double y, double tol, double rtol)
{
    return std::abs(x - y) <= tol + rtol * std::abs(y);
}

inline bool is_zero(double x, double tol)
{
    return std::abs(x) <= tol;
}

inline bool same_signs(double x, double y)
{
    return x * y > 0;
}

inline bool opposite_signs(double x, double y)
{
    return x * y < 0;
}

inline bool smaller(double x, double y)
{
    return std::abs(x) < std::abs(y);
}

bool successful_termination(Flag flag)
{
    switch (flag) {
        case Flag::residual_zero:
            return true;
        case Flag::delta_root_zero:
            return true;
        case Flag::bracket_width_zero:
            return true;
        default:
            return false;
    }
}

std::string flag_message(Flag flag){
    switch (flag) {
        case Flag::residual_zero:
            return "Residual is zero.";
        case Flag::delta_root_zero:
            return "Change in guess is zero.";
        case Flag::bracket_width_zero:
            return "Bracket width is zero.";
        case Flag::invalid_bracket:
            return "Bracket is invalid; Function has same signs at both endpoints.";
        case Flag::max_iterations:
            return "Reached the maximum number of iterations.";
        case Flag::division_by_zero:
            return "Division by zero occurred.";
        case Flag::valid:
            return "Valid state, but convergence not reached.";
        default:
            return "Flag not recognized.";
    }
}


// Test functions
struct hard_test {

    double answer() { return 1.; }

    double operator()(double x){
        if (x < 0) return -0.25;
        return (std::pow(x, 4) - 1) / 4;
    }

    double derivative(double x){
        if (x < 0) return 0.;
        return std::pow(x, 3);
    }

};

struct easy_test {
    double operator()(double x) {
        return (1 - x) * (x + 1);
    }

    double derivative(double x){
        return - 2 * x ;
    }

};

struct double_root {

    double operator()(double x) {
        return std::pow(x - 1, 2);
    }

    double derivative(double x){
        return 2 * (x - 1) ;
    }
};


struct triple_root {

    double operator()(double x) {
        return std::pow(x - 1, 3);
    }

    double derivative(double x){
        return 3 * std::pow(x - 1, 2);
    }

};

}  // namespace root_algorithm
#endif
