#ifndef ROOTS_ONEDIM_H
#define ROOTS_ONEDIM_H

#include <cmath>
#include <limits>

namespace root_algorithm
{

// Helper function declarations
inline bool is_close(double x, double y, double tol, double rtol);
inline bool is_zero(double x, double tol);
inline bool same_signs(double x, double y);
inline bool opposite_signs(double x, double y);

struct graph_t {
    double x;
    double y;
};

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
struct root_finder {
    size_t max_iterations = 100;
    double _abs_tol = std::numeric_limits<double>::epsilon();
    double _rel_tol = 2 * std::numeric_limits<double>::epsilon();

    root_finder(size_t max_iter, double abs_tol, double rel_tol) : max_iterations{max_iter}, _abs_tol{abs_tol}, _rel_tol{rel_tol} {}

    template <typename F, typename... Args>
    result_t solve(const F& func, Args&&... args)
    {
        typename Method::state s = Method::initialize(func, std::forward<Args>(args)...);

        for (size_t i = 0; i < (max_iterations + 1); ++i) {
            if (Method::has_converged(s, _abs_tol, _rel_tol)) {
                return result_t{Method::root(s), Method::residual(s), i};
            };

            s = Method::iterate(func, s);
        }
        return result_t{Method::root(s), Method::residual(s), max_iterations};
    }
};

// Methods
// Multi-step methods

struct secant {
    struct state {
        graph_t last;
        graph_t current;
    };

    template <typename F>
    static state initialize(const F& fun, double x0, double x1)
    {
        graph_t first{x0, fun(x0)};
        graph_t second{x1, fun(x1)};
        return state{first, second};
    }

    template <typename F>
    static state iterate(const F& fun, const state& s)
    {
        double secant_slope = (s.current.y - s.last.y) / (s.current.x - s.last.x);
        double x2 = s.current.x - s.current.y / secant_slope;
        graph_t new_guess{x2, fun(x2)};
        return state{s.current, new_guess};
    }

    static inline double root(const state& s)
    {
        return s.current.x;
    }

    static inline double residual(const state& s)
    {
        return s.current.y;
    }

    static inline bool has_converged(const state& s, double abs_tol, double rel_tol)
    {
        bool a = is_zero(s.current.y, abs_tol);
        bool b = is_close(s.last.x, s.current.x, abs_tol, rel_tol);
        return a || b;
    }
};

struct steffensen {
    struct state {
        graph_t last;
        graph_t current;
    };

    template <typename F>
    static state initialize(F&& fun, double x0)
    {
        state s;
        s.current = {x0, fun(x0)};
        return iterate(fun, s);
    }

    template <typename F>
    static state iterate(F&& fun, const state& s)
    {
        double gx = fun(s.current.x + s.current.y) / s.current.y - 1;
        double x_new = s.current.x - s.current.y / gx;
        graph_t new_guess{x_new, fun(x_new)};
        return state{s.current, new_guess};
    }

    static inline double root(const state& s)
    {
        return s.current.x;
    }

    static inline double residual(const state& s)
    {
        return s.current.y;
    }

    static inline bool has_converged(const state& s, double abs_tol, double rel_tol)
    {
        bool a = is_zero(s.current.y, abs_tol);
        bool b = is_close(s.last.x, s.current.x, abs_tol, rel_tol);
        return a || b;
    }
};

// Householder methods

struct newton {
    struct state {
        graph_t last;
        graph_t current;
    };

    template <typename F>
    static state initialize(F&& fun, double x0)
    {
        state s;
        s.current = {x0, fun(x0)};
        return iterate(fun, s);
    }

    template <typename F>
    static state iterate(F&& fun, const state& s)
    {
        double x_new = s.current.x - s.current.y / fun.derivative(s.current.x);
        graph_t new_guess{x_new, fun(x_new)};
        return state{s.current, new_guess};
    }

    static inline double root(const state& s)
    {
        return s.current.x;
    }

    static inline double residual(const state& s)
    {
        return s.current.y;
    }

    static inline bool has_converged(const state& s, double abs_tol, double rel_tol)
    {
        bool a = is_zero(s.current.y, abs_tol);
        bool b = is_close(s.last.x, s.current.x, abs_tol, rel_tol);
        return a || b;
    }
};

// Bracketing methods

struct bisection {
    struct state {
        graph_t left;
        graph_t right;
    };

    template <typename F>
    static state initialize(F&& fun, double a, double b)
    {
        graph_t left{a, fun(a)};
        graph_t right{b, fun(b)};
        return state{left, right};
    }

    template <typename F>
    static state iterate(F&& fun, const state& s)
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

    static inline double root(const state& s)
    {
        return (s.left.x + s.right.x) / 2;
    }

    static inline double residual(const state& s)
    {
        return (s.left.y + s.right.y) / 2;
    }

    static inline bool has_converged(const state& s, double abs_tol, double rel_tol)
    {
        bool a = is_zero(s.left.y, abs_tol);
        bool b = is_zero(s.right.y, abs_tol);
        bool c = is_close(s.left.x, s.right.x, abs_tol, rel_tol);
        return a || b || c;
    }
};

// Helper function definitions
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

// Test functions
double testfun(double x)
{
    if (x < 0) return -0.25;
    return (std::pow(x, 4) - 1) / 4;
}

double testfunjac(double x)
{
    if (x < 0) return 0.;
    return std::pow(x, 3);
}

}  // namespace root_algorithm
#endif
