#ifndef ROOTS_ONEDIM_H
#define ROOTS_ONEDIM_H

#include <cmath>
#include <limits>
#include <string>

namespace root_algorithm
{

// Helper function declarations
inline bool is_close(double x, double y, double tol, double rtol);
inline bool is_close(double x, double y);
inline bool is_zero(double x, double tol);
inline bool is_zero(double x);
inline bool same_signs(double x, double y);
inline bool opposite_signs(double x, double y);
inline bool smaller(double x, double y);
inline bool is_between(double x, double a, double b);
inline bool all_distinct(double x, double y, double z);

// For error handling.
enum class Flag {
    valid,  // don't terminate
    residual_zero,
    delta_root_zero,
    bracket_width_zero,
    max_iterations,
    invalid_bracket,
    division_by_zero,
    halley_no_cross
};
bool successful_termination(Flag flag);
std::string flag_message(Flag flag);

/*
For error handling, this struct adds flags to a state used in the
iterations. This is designed to be a monad, meaning the following rules
should be respected: For any function/function object with the signature:
`T f(T x)` there is a function `with_flag<T> f(with_flag<T> x)` such that
if (is_valid(x)) return with_flag<T>{f(x.state)}; else return x;

*/
template <typename T>
struct with_flag {
    T state;
    Flag flag;

    with_flag() = default;
    with_flag(T x) : state(x), flag(Flag::valid) {}
    with_flag(T x, Flag f) : state(x), flag(f) {}
    with_flag(Flag f) : state{}, flag(f) {}
};

template <typename T>
inline bool is_valid(with_flag<T> x)
{
    return x.flag == Flag::valid;
}

/**
 * @class result_t
 *
 * @brief Result from root finding algorithm.
 *
 * @param root The identified root.
 *
 * @param residual The output value (root); should equal zero if
 * successful.
 *
 * @param iteration The number of iterations performed.
 *
 * @param flag Indicates the reason for termination.
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
 * @details This class creates a function object whose solve or call method
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
 * We will say function object for all of these options.
 * For algorithms which explicitly evaluate derivatives, you must
 * pass a function object with a `double derivative(double x)` method.
 *
 * See documentation for each method for details.
 *
 * List of methods:
 *
 * + newton
 * + halley
 * + steffensen
 * + secant
 * + bisection
 * + regula_falsi
 * + ridder
 *
 * The robustness and efficiency of each method is problem-dependent.
 * Root-bracketing methods are typically more robust, but methods that use
 * derivatives almost always take fewer iterations. Computational efficiency
 * depends on how expensive evaluating the function (or its derivatives)
 * is. The function calls per iteration can be counted.
 *
 * A new algorithm can be added by creating a class or struct with the
 * following methods:
 *
 * @code{.cpp}
 * //holds whatever state saved between iterations.
 * struct state;
 * //takes the initial problem info (e.g., a bracket, an
 * // initial guess) and instantiates a value of type `state`. If the
 * provided info is invalid, then return a flag indicating the error.
 * with_flag<state> initialize(F&& f, Args... args);
 * //maps state to state. The actual formula of most methods is
 * // implemented here. If the iteration encounters an error, return the
 * // flag indicating the error.
 * // Some methods handle errors that other methods do not.
 * with_flag<state> iterate(F&& f, const state& s);
 * // checks the state to see if the algorithm has found a root to within
 * // tolerance; return the appropriate flag, else return `Flag::valid`
 * Flag check_convergence(const state& s, _abs_tol, _rel_tol);
 * // Extract a single value for the root.
 * double root(const state& s);
 * // Extract the residual from the state.
 * double residual(const state& s);
 * @endcode
 *
 * Several methods have the same convergence checks and state; so a struct
 * provides a template for those methods. New methods can (but do not have
 * to) inherit from these templates if they are useful.
 *
 * `initialize` takes the initial problem info (like the two guesses in
 * the secant method) and creates an object that holds the state that is
 * kept from one iteration to the next. `iterate` contains the meat of each
 *  method; it refines the guess, etc. `check_convergence` checks if state has
 *  converged; this means the `state` object has to have whatever info is
 * used to make this decision. `root` and `residual` extract the last best
 *  guess and the function's residual at that root from the `state`.
 *
 *
 */
template <typename Method>
struct root_finder : public Method {
    size_t max_iterations = 100;
    double _abs_tol = std::numeric_limits<double>::epsilon();
    double _rel_tol = 2 * std::numeric_limits<double>::epsilon();

    root_finder() {}
    root_finder(size_t max_iter, double abs_tol, double rel_tol)
        : max_iterations{max_iter},
          _abs_tol{abs_tol},
          _rel_tol{rel_tol} {}

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

            if (!is_valid(s)) {
                return make_result(s, i);
            }

            s = Method::iterate(std::forward<F>(func), s.state);
        }

        s.flag = Flag::max_iterations;
        return make_result(s, max_iterations);
    }

    template <typename F, typename... Args>
    result_t operator()(F&& func, Args&&... args)
    {
        return solve(std::forward<F>(func), std::forward<Args>(args)...);
    }

   private:
    inline result_t make_result(const flagged_state& s, size_t iteration)
    {
        return result_t{Method::root(s.state), Method::residual(s.state), iteration, s.flag};
    }
};

struct graph_t {
    double x;
    double y;
};

// Householder and multistep methods

/**
 * @brief Secant Method. Provide two initial guesses for root.
 *
 * @details
 */
struct secant {
    struct state {
        graph_t last;
        graph_t best;
    };

    template <typename F>
    with_flag<state> initialize(F&& fun, double x0, double x1)
    {
        graph_t first = {x0, fun(x0)};
        graph_t second = {x1, fun(x1)};
        return state{first, second};
    }

    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        double r = s.best.y / s.last.y;
        double p = (s.best.x - s.last.x) * r;
        double q = 1 - r;
        if (is_zero(q)) {
            return {s, Flag::division_by_zero};
        }
        double x2 = s.best.x + p / q;
        graph_t new_guess = {x2, fun(x2)};
        return state{s.best, new_guess};
    }

    Flag check_convergence(const state& s, double abs_tol, double rel_tol)
    {
        if (is_zero(s.best.y, abs_tol)) {
            return Flag::residual_zero;
        }

        if (is_close(s.last.x, s.best.x, abs_tol, rel_tol)) {
            return Flag::delta_root_zero;
        }

        return Flag::valid;
    }

    inline double root(const state& s)
    {
        return s.best.x;
    }

    inline double residual(const state& s)
    {
        return s.best.y;
    }
};

// Householder methods

// methods common to newton, halley, etc.
struct one_step_method {
    using state = graph_t;

    template <typename F>
    with_flag<state> initialize(F&& fun, double x0)
    {
        return state{x0, fun(x0)};
    }

    Flag check_convergence(const state& s, double abs_tol, double rel_tol)
    {
        return is_zero(s.y, abs_tol) ? Flag::residual_zero : Flag::valid;
    }

    inline double root(const state& s)
    {
        return s.x;
    }

    inline double residual(const state& s)
    {
        return s.y;
    }
};

/**
 * @brief Newton's Method. Provide a function object with `double derivative(double x)`
 * method implementing the derivative. Provide an initual guesss.
 *
 * @details
 */
struct newton : public one_step_method {
    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        double slope = fun.derivative(s.x);
        if (is_zero(slope)) {
            return {s, Flag::division_by_zero};
        }
        double x_new = s.x - s.y / slope;
        return graph_t{x_new, fun(x_new)};
    }
};

/**
 * @brief Halley's Method. Provide a function object with methods
 * `double derivative(double x)` and `double second_derivative(double x)`
 * implementing the first and second derivatives. Provide an initual guesss.
 *
 * @details
 */
struct halley : public one_step_method {
    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        double df = fun.derivative(s.x);
        if (is_zero(df)) {
            return {s, Flag::division_by_zero};
        }
        double df2 = fun.second_derivative(s.x);
        double a = s.y / df;
        double b = df2 / (2 * df);
        double ab = a * b;
        if (is_close(ab, 1)) {
            return {s, Flag::halley_no_cross};
        }
        double x_new = s.x - a / (1. - a * b);  // is division by zero possible here? Only if no solution
        return graph_t{x_new, fun(x_new)};
    }
};

/**
 * @brief Steffensen's Method. Provide a function object and
 * an initial guesss.
 *
 * @details
 */
struct steffensen : public one_step_method {
    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        double g = fun(s.x + s.y) / s.y - 1;
        if (is_zero(g)) {
            return {s, Flag::division_by_zero};
        }
        double x_new = s.x - s.y / g;
        return graph_t{x_new, fun(x_new)};
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
        if (same_signs(left.y, right.y)) {
            return {state{left, right}, Flag::invalid_bracket};
        }
        return state{left, right};
    }

    Flag check_convergence(const state& s, double abs_tol, double rel_tol)
    {
        bool a = is_zero(s.left.y, abs_tol);
        bool b = is_zero(s.right.y, abs_tol);
        if (a || b) return Flag::residual_zero;

        bool c = is_close(s.left.x, s.right.x, abs_tol, rel_tol);
        if (c) return Flag::bracket_width_zero;

        return Flag::valid;
    }

    inline double root(const state& s)
    {
        return smaller(s.left.y, s.right.y) ? s.left.x : s.right.x;
    }

    inline double residual(const state& s)
    {
        return smaller(s.left.y, s.right.y) ? s.left.y : s.right.y;
    }

    // update bracket so that sign difference is maintained
    inline state update_bracket(const state& current, const graph_t& new_point)
    {
        return same_signs(current.left.y, new_point.y) ? state{new_point, current.right} : state{current.left, new_point};
    }
};

/**
 * @brief The Bisection Method. Provide a function object and an initial
 * valid bracket. Bracket is valid if the sign of the function differs at
 * the end points.
 *
 * @details
 */
struct bisection : public bracket_method {
    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        graph_t midpoint = {
            (s.left.x + s.right.x) / 2,
            fun(midpoint.x)};

        return update_bracket(s, midpoint);
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

/**
 * @brief Regula falsi, or the false position method.
 * Provide a function object and an initial valid bracket.
 * Bracket is valid if the sign of the function differs at the end points.
 *
 * @details
 */
struct regula_falsi : public bracket_method {
    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        double slope = (s.right.y - s.left.y) / (s.right.x - s.left.x);

        graph_t new_point;
        new_point.x = s.right.x - s.right.y / slope;
        new_point.y = fun(new_point.x);

        return update_bracket(s, new_point);
    }
};

/**
 * @brief The Bisection Method. Provide a function object and an initial
 * valid bracket. Bracket is valid if the sign of the function differs at
 * the end points.
 *
 * @details
 */
struct ridder : public bracket_method {
    template <typename F>
    with_flag<state> iterate(F&& fun, const state& s)
    {
        graph_t midpoint = {
            (s.left.x + s.right.x) / 2,
            fun(midpoint.x)};

        double d = midpoint.x - s.left.x;
        double a = midpoint.y / s.left.y;
        double b = s.right.y / s.left.y;

        graph_t new_point;
        new_point.x = midpoint.x + d * a / std::sqrt(a * a - b);
        if (!std::isfinite(new_point.x)) {
            return {s, Flag::division_by_zero};  // is this state possible ?
        }
        new_point.y = fun(new_point.x);
        return update_bracket(s, new_point);
    }
};

// Helper function definitions.
inline bool is_close(double x, double y, double tol, double rtol)
{
    return std::abs(x - y) <= tol + rtol * std::abs(y);
}

inline bool is_close(double x, double y)
{
    constexpr double tol = 2 * std::numeric_limits<double>::epsilon();
    return is_close(x, y, tol, tol);
}

inline bool is_zero(double x, double tol)
{
    return std::abs(x) <= tol;
}

inline bool is_zero(double x)
{
    constexpr double tol = 2 * std::numeric_limits<double>::epsilon();
    return is_zero(x, tol);
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

inline bool is_between(double x, double a, double b)
{
    return ((x >= a) && (x <= b)) || ((x <= a) && (x >= b));
}

inline bool all_distinct(double x, double y, double z)
{
    constexpr double tol = std::numeric_limits<double>::epsilon();
    bool a = is_close(x, y, tol, tol);
    a |= is_close(y, z, tol, tol);
    a |= is_close(x, y, tol, tol);
    return !a;
}

bool successful_termination(Flag flag)
{
    switch (flag) {
        case Flag::residual_zero:
            return true;
        case Flag::bracket_width_zero:
            return true;
        default:
            return false;
    }
}

std::string flag_message(Flag flag)
{
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
        case Flag::halley_no_cross:
            return "Halley update failed; local quadratic does not cross zero.";
        case Flag::valid:
            return "Valid state, but convergence not reached.";
        default:
            return "Flag not recognized.";
    }
}

// Test functions
struct hard_test {
    double answer() { return 1.; }

    double operator()(double x)
    {
        return x < 0 ? -0.25 : (std::pow(x, 4) - 1) / 4;
    }

    double derivative(double x)
    {
        return x < 0 ? 0.0 : std::pow(x, 3);
    }

    double second_derivative(double x)
    {
        return x < 0 ? 0.0 : 3 * std::pow(x, 2);
    }
};

struct easy_test {
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
}  // namespace root_algorithm
#endif
