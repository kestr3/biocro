#ifndef ROOT_NEWTON_H
#define ROOT_NEWTON_H

#include "roots.h"

namespace root_finding
{
/**
 * @brief Newton's Method. Provide a function object with `double derivative(double x)`
 * method implementing the derivative. Provide an initual guesss.
 *
 * @details Newton's method uses the zero of a first-order Taylor series
 * approximation of `f`to refine guesses for the root.
 *
 * Newton's method is quite fast, especially if the derivative and
 * function can be evaluated simultaneously. Convergence under optimal
 * conditions is quadratic.
 *
 * Regions where the first derivative is zero can cause divergent behavior.
 * Chaotic behavior can occur at the boundary between basins of attraction.
 * Non-simple roots tend to slow convergence.
 *
 * Use Newton's method preferentially if the derivatives can be evaluated. If not,
 * use the secant method, rather than approximate the first derivative.
 *
 */
struct newton : public root_finding_method<newton> {
    using root_finding_method::root_finding_method;
    graph_t current;

    template <typename F>
    bool initialize(F&& fun, double x0)
    {
        current.x = x0;
        current.y = fun(x0);
        if (is_zero(current.y)) {
            flag = Flag::residual_zero;

            return false;
        }
        return true;
    }

    template <typename F>
    bool iterate(F&& fun)
    {
        double slope = fun.derivative(current.x);
        if (is_zero(slope)) {
            flag = Flag::division_by_zero;
            return false;
        }
        current.x -= current.y / slope;
        current.y = fun(current.x);
        return true;
    }

    bool has_converged()
    {
        if (is_zero(current.y)) {
            flag = Flag::residual_zero;
            return true;
        }
        return false;
    }

    double root() const
    {
        return current.x;
    }

    double residual() const
    {
        return current.y;
    }
};

/**
 * @brief Halley's Method. Provide a function object with methods
 * `double derivative(double x)` and `double second_derivative(double x)`
 * implementing the first and second derivatives. Provide an initual guesss.
 *
 * @details Halley's rational method uses a first and second order Taylor series
 * approximation to refine a guess for the root. Essentially Newton's method
 * is used to pick a direction, and the second derivative is used to pick a
 * good step size (reducing the step size when high curvature is present).
 *
 * Halley's method is even faster than Newton's method, but requires evaluating
 * the second derivative. Not usually feasible for problems where `f`
 * involves complicated expressions, but fast for polynomials or simple
 * transcendental equations.
 *
 * Only use if derivatives can be evaluated directly, prefer the `secant`
 * method if derivatives cannot be evaluated.
 *
 * This method is Halley's rational method. Halley's irrational method solves the
 * quadratic of the Taylor series, using the Newton step to select a root of the
 * quadratic.
 */
struct halley : public newton {
    using newton::newton;

    template <typename F>
    bool iterate(F&& fun)
    {
        double df = fun.derivative(current.x);
        if (is_zero(df)) {
            flag = Flag::division_by_zero;
            return false;
        }
        double df2 = fun.second_derivative(current.x);
        double a = current.y / df;
        double b = df2 / (2 * df);
        b *= a;
        if (is_close(b, 1)) {
            flag = Flag::halley_no_cross;
            return false;
        }
        current.x -= a / (1. - b);  // division by zero Only if no solution
        current.y = fun(current.x);
        return true;
    }
};

/**
 * @brief Steffensen's Method. Provide a function object and
 * an initial guesss.
 *
 * @details Steffensen's method approximates the function `f` with
 * a first order Taylor series, using the formula:
 *
 * \f[ f'(x) \approx \frac{f(x + h) - h}{h} \qquad h = f(x) \f]
 *
 * to approximate the first derivative. This expression requires two function
 * evaluations per iteration. This experession achieves quadratic convergence
 * rates under optimal conditions. I.e., it can be as fast as Newton's method
 * however, it is far less stable than either Newton's method or the secant
 * method if the initial guesss is not close to a root. Moreover, the rate
 * of convergence per function call is higher for the secant method.
 *
 */
struct steffensen : public newton {
    using newton::newton;

    template <typename F>
    bool iterate(F&& fun)
    {
        double g = fun(current.x + current.y) / current.y - 1;
        if (is_zero(g)) {
            flag = Flag::division_by_zero;
            return false;
        }
        current.x -= current.y / g;
        current.y = fun(current.x);
        return false;
    }
};

}  // namespace root_finding

#endif
