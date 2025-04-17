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

template <typename T>
inline bool is_valid(T x)
{
    return x.flag == Flag::valid;
}
inline bool successful_termination(Flag flag);
inline std::string flag_message(Flag flag);

/*
For error handling, this struct adds flags to a state used in the
iterations. This is designed to be a monad, meaning the following rules
should be respected: For any function/function object with the signature:
`T f(T x)` there is a function `with_flag<T> f(with_flag<T> x)` such that
if (is_valid(x)) return with_flag<T>{f(x.state)}; else return x;

*/
// template <typename T>
// struct with_flag {
//     T state;
//     Flag flag;

//     with_flag() = default;
//     with_flag(T x) : state(x), flag(Flag::valid) {}
//     with_flag(T x, Flag f) : state(x), flag(f) {}
//     with_flag(Flag f) : state{}, flag(f) {}
// };


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
 * state initialize(F&& f, Args... args);
 * //maps state to state. The actual formula of most methods is
 * // implemented here. If the iteration encounters an error, return the
 * // flag indicating the error.
 * // Some methods handle errors that other methods do not.
 * state& iterate(F&& f, state& s);
 * // iterate is implemented like a compound assignment operator,
 * // updating in-place rather than using an immutable maping.
 * // It should be equivalent to state iterate(f, const state& s)
 * // Where it maps a state to state, and adds a new flag if an error is
 * // encountered.
 *
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
    static constexpr double eps = std::numeric_limits<double>::epsilon();
    size_t max_iterations = 100;
    double _abs_tol = 2 * eps;
    double _rel_tol = 2 * eps;

    root_finder() {}
    root_finder(size_t max_iter) : max_iterations{max_iter},
          _abs_tol{2 * eps},
          _rel_tol{2 * eps} {}
    root_finder(size_t max_iter, double abs_tol, double rel_tol)
        : max_iterations{max_iter},
          _abs_tol{abs_tol},
          _rel_tol{rel_tol} {}

    using state = typename Method::state;

    template <typename F, typename... Args>
    result_t solve(F&& func, Args&&... args)
    {
        state s = Method::initialize(
            std::forward<F>(func), std::forward<Args>(args)...,
            _abs_tol, _rel_tol
        );

        for (size_t i = 0; i < (max_iterations + 1); ++i) {

            if (!is_valid(s)) {
                return make_result(s, i);
            }

            s = Method::iterate(std::forward<F>(func), s, _abs_tol, _rel_tol);

            s = Method::check_convergence(s, _abs_tol, _rel_tol);

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
    inline result_t make_result(const state& s, size_t iteration)
    {
        return result_t{Method::root(s), Method::residual(s), iteration, s.flag};
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
        Flag flag;
        graph_t last;
        graph_t best;
    };

    template <typename F>
    inline state initialize(F&& fun, double x0, double x1,
        double abs_tol, double rel_tol)
    {
        graph_t first = {x0, fun(x0)};
        graph_t second = {x1, fun(x1)};
        if (is_zero(first.y, abs_tol)) {
            return state{Flag::residual_zero, second, first};
        }
        if (is_zero(second.y, abs_tol)) {
            return state{Flag::residual_zero, first, second};
        }
        return state{Flag::valid, first, second};
    }

    template <typename F>
    inline state& iterate(F&& fun, state& s, double abs_tol, double rel_tol)
    {

        auto& best = s.best;
        auto& last = s.last;
        double r = best.y / last.y;
        double p = (best.x - last.x) * r;
        double q = 1 - r;
        if (is_zero(q, abs_tol)) {
            s.flag = Flag::division_by_zero;
            return s;
        }
        last = best;
        best.x += p / q;
        best.y = fun(best.x);
        return s;
    }

    inline state& check_convergence(state& s, double abs_tol, double rel_tol)
    {
        if (is_zero(s.best.y, abs_tol)) {
            s.flag = Flag::residual_zero;
            return s;
        }

        if (is_close(s.last.x, s.best.x, abs_tol, rel_tol)) {
            s.flag = Flag::delta_root_zero;
            return s;
        }

        return s;
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
    struct state {
        Flag flag;
        double x;
        double y;
    };

    template <typename F>
    state initialize(F&& fun, double x0, double abs_tol, double rel_tol)
    {
        double y0 = fun(x0);
        if(is_zero(y0, abs_tol)) {
            return state{Flag::residual_zero, x0, fun(x0)};
        }
        return state{Flag::valid, x0, fun(x0)};
    }

    inline state& check_convergence(state& s, double abs_tol, double rel_tol)
    {
        if (is_zero(s.y, abs_tol)){
            s.flag = Flag::residual_zero;
        }
        return s;
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
    inline state& iterate(F&& fun, state& s, double abs_tol, double rel_tol)
    {
        double& x = s.x;
        double& y = s.y;
        double slope = fun.derivative(x);
        if (is_zero(slope, abs_tol)) {
            s.flag = Flag::division_by_zero;
            return s;
        }
        x -= y / slope;
        y = fun(x);
        return s;
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
    inline state& iterate(F&& fun, state& s, double abs_tol, double rel_tol)
    {
        double& x = s.x;
        double& y = s.y;
        double df = fun.derivative(x);
        if (is_zero(df, abs_tol)) {
            s = Flag::division_by_zero;
            return s;
        }
        double df2 = fun.second_derivative(x);
        double a = y / df;
        double b = df2 / (2 * df);
        b *= a;
        if (is_close(b, 1, abs_tol, rel_tol)) {
            s = Flag::halley_no_cross;
            return s;
        }
        x -= a / (1. - b);  // division by zero Only if no solution
        y = fun(x);
        return s;
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
    inline state& iterate(F&& fun, state& s, double abs_tol, double rel_tol)
    {
        double& x = s.x;
        double& y = s.y;
        double g = fun(x + y) / y - 1;
        if (is_zero(g, abs_tol)) {
            s.flag = Flag::division_by_zero;
            return s;
        }
        x -= y / g;
        y = fun(x);
        return s;
    }
};

// Bracketing methods
template<typename UpdateMethod>
struct bracket_method {
    struct state {
        Flag flag;
        double lower;
        double width;
        double proposal;
        double fun_proposal;
    };

    template <typename F>
    state initialize(F&& fun, double a, double b, double abs_tol, double rel_tol)
    {
        double ya = fun(a);
        if (is_zero(ya, abs_tol)) {
            return state{Flag::residual_zero, a, 0, a, ya};
        }

        double yb = fun(b);
        if (is_zero(right.y, abs_tol)) {
            return state{Flag::residual_zero, b, 0, b, yb};
        }

        if (same_signs(ya, yb)) {
            return state{Flag::invalid_bracket, a, b - a, b, yb};
        }

        if (ya < 0) {
            return state{Flag::valid, a, b - a, b, yb};
        }

        return state{Flag::valid, b, a - b, a, yb};
    }

    template <typename F>
    inline state& iterate(F&& fun, state& s, double abs_tol, double rel_tol)
    {
        s.fun_proposal = fun(s.proposal = UpdateMethod::propose(state& s));
        if (s.fun_proposal <= 0.0){
            s.lower = s.proposal;
        }
        return s;
    }

    inline state& check_convergence(state& s, double abs_tol, double rel_tol)
    {

        if (is_zero(s.fun_proposal, abs_tol)) {
            s.flag = Flag::residual_zero;
            return s;
        }

        if (is_zero(s.width, abs_tol))
            s.flag = Flag::bracket_width_zero;

        return s;
    }

    inline double root(const state& s)
    {
        return s.proposal;
    }

    inline double residual(const state& s)
    {
        return s.fun_proposal;
    }

    // update bracket so that sign difference is maintained

};

/**
 * @brief The Bisection Method. Provide a function object and an initial
 * valid bracket. Bracket is valid if the sign of the function differs at
 * the end points.
 *
 * @details
 */
struct bisection : public bracket_method<bisection> {
    template <typename F>
    inline double propose(state& s)
    {
        return lower + (s.width /= 2);
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
    inline state& iterate(F&& fun, state& s, double abs_tol, double rel_tol)
    {
        graph_t& left = s.left;
        graph_t& right = s.right;
        double r = right.y / left.y; // equivalent to a secant update.
        double p = (right.x - left.x) * r;
        double q = 1 - r;
        if (is_zero(q, abs_tol)) {
            s.flag = Flag::division_by_zero;
            return s;
        }
        double c = right.x + p / q;
        double fc = fun(c);
        s = update_bracket(s, c, fc, abs_tol);
        return s;
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
    inline state& iterate(F&& fun, state& s, double abs_tol, double rel_tol)
    {
        graph_t& left = s.left;
        graph_t& right = s.right;
        double c = (left.x + right.x) / 2;
        double yc = fun(c);

        double d = c - left.x;
        double a = yc / left.y;
        double b = right.y / left.y;
        double denom = a * a - b;
        if (is_zero(denom, abs_tol)) {
            s.flag = Flag::division_by_zero;  // is this state possible ?
            return s;
        }

        c  += d * a / std::sqrt(denom);
        yc = fun(c);
        s = update_bracket(s, c, yc, abs_tol);
        return s;
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
    return x * y >= 0;
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
