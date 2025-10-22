#ifndef ROOT_BRACKET_H
#define ROOT_BRACKET_H

#include "roots.h"

namespace root_finding
{
/**
 * @brief The Bisection Method. Provide a function object and an initial
 * valid bracket. Bracket is valid if the sign of the function differs at
 * the end points.
 *
 * @details The bisection method is a bracketing method: if `(a, b)` is
 * either side of a simple root then `f(a)` and `f(b)` will have different
 * signs. The bisection method finds roots by finding where `f` switches
 * sign between positive and negative. At each iteration, this method
 * divides the interval in half `c = (a + b) / 2` and then replace `a` with
 * `c` if the sign of `f(c)` is the same as `f(a)`, otherwise `c` replaces
 * `b`.
 *
 * For a continuous function, this method is guaranteed to converge to
 * a root, provided a valid bracket. However, its convergence is quite slow
 * as it only halves the interval at each step. Regardless of problem,
 * 50-150 iterations is common depending on the desired accuracy.
 *
 * The bisection method is not guaranteed to converge to a double root.
 * It does work for any number of roots provided that there is an odd number
 * inside the bracket. It can also converge to discontinuities and singularities
 * if the function switches sign at those points. For instance, the bisection
 * method will correctly identify `x = 0` as the sign switch point for
 * `f(x) = 1 / x` or for `f(x) = signum(x)`.
 *
 * Note `abs_tol` is used to determine if a residual is zero or if the
 * bracket width is effectively zero. `rel_tol` sets the tolerance for
 * deciding if the function is continuous at the estimated root.
 *
 * References:
 *
 * - Press et al. (2007). Numerical recipes, 3rd edition.
 *   Cambridge University Press. https://numerical.recipes/book.html
 *
 * Common to all bracketing methods defined here  In general, a bracket method maintains
 * an interval, a pair of points below and above the root.
 *
 * The `state` struct holds the `proposal` so that each new point can be checked
 * before it's used to update the bracket.
 *
 * These methods use `abs_tol` to test whether the residual is zero and whether
 * the bracket width is zero. For the latter, `abs_tol` is both the absolute and
 * relative tolerance. It is possible that the bracket zeroes in on a root without
 * reaching the absolute tolerance. The `rel_tol` is only used to determine
 * whether or not the bracket of zero-width contains a zero, by testing continuity.
 */

struct bisection : public root_finding_method<bisection> {
    using root_finding_method::root_finding_method;

    graph_t left;
    graph_t right;
    graph_t proposal;

    template <typename F>
    bool initialize(F&& fun, double a, double b)
    {
        left.x = a;
        right.x = b;
        left.y = fun(a);
        right.y = fun(b);
        proposal = left;

        if (is_zero(left.y)) {
            flag = Flag::residual_zero;
            return false;
        }

        if (is_zero(right.y)) {
            flag = Flag::residual_zero;
            proposal = right;
            return false;
        }

        if (same_signs(left.y, right.y)) {
            flag = Flag::invalid_bracket;
            return false;
        }

        return false;
    }

    template <typename F>
    bool iterate(F&& fun)
    {
        proposal.x = get_midpoint(left, right);
        proposal.y = fun(proposal.x);
        update_bracket();
        return true;
    }

    bool has_converged()
    {
        if (is_zero(proposal.y)) {
            flag = Flag::residual_zero;
            return true;
        }

        if (is_close(left.x, right.x)) {
            flag = Flag::bracket_width_zero;
            return true;
        }

        return false;
    }

    double root() const
    {
        return proposal.x;
    }

    double residual() const
    {
        return proposal.y;
    }

    void update_bracket()
    {
        if (same_signs(left.y, proposal.y)) {
            left = proposal;
        } else {
            right = proposal;
        }
    }
};

/**
 * @brief Regula falsi, or the false position method.
 * Provide a function object and an initial valid bracket.
 * Bracket is valid if the sign of the function differs at the end points.
 *
 * @details Regula falsi, or the false position method, is a bracketing
 * method: if `(a, b)` is either side of a simple root then `f(a)` and
 * `f(b)` will have different signs. The method finds roots by finding
 * where `f` switches sign between positive and negative.
 * This method computes a new bracket using the secant of the end points,
 * replacing the end points based on the signs.
 *
 * Regula falsi has the same safety as the bisection method but with better
 * speed. However, it usually is not as fast as the secant method. The secant
 * method uses the best guesses found so far to achieve a better convergence
 * rate, while regula falsi sacrifices some speed for robustness.
 *
 * Regula falsi can converge much slower than the bisection method for a
 * function with large curvature within the bracket. Regula falsi doesn't
 * necessarily produce a sequence of brackets that shrinks to zero.
 *
 * Thus, try regula falsi and use bisection if regula falsi fails to converge.
 *
 * References:
 *
 * - Press et al. (2007). Numerical recipes, 3rd edition.
 *   Cambridge University Press. https://numerical.recipes/book.html
 *
 */
struct regula_falsi : public bisection {
    using bisection::bisection;

    template <typename F>
    bool iterate(F&& fun)
    {
        proposal.x = get_secant_update(left, right);
        proposal.y = fun(proposal.x);
        update_bracket();
        return true;
    }
};

/**
 * @brief Ridder's Method. Provide a function object and an initial
 * valid bracket. Bracket is valid if the sign of the function differs at
 * the end points.
 *
 * @details Ridder's method is a bracketing method: if `(a, b)` is either
 * side of a simple root then `f(a)` and `f(b)` will have different signs.
 * The method finds roots by finding where `f` switches sign between
 * positive and negative.
 *
 * In essence, Ridder's method fits an exponential to the end points and
 * the mid point to estimate `m` in:
 *
 * \f[ h(x) = f(x) \exp(m x) \f]
 *
 * The exponential fit corrects high curvature, especially for exponential
 * like functions. A secant `h(x)` is then used to generate a new bracket
 * endpoint.
 *
 * Ridder's method can perform well on problems that frustrate regula falsi,
 * but it is slower on quadratics.
 *
 * References:
 * - Ridders, C. (1979). "A new algorithm for computing a single root of
 *   a real continuous function". IEEE Transactions on Circuits and Systems.
 *   26 (11): 979–980. doi:10.1109/TCS.1979.1084580
 *
 * - Press et al. (2007). Numerical recipes, 3rd edition.
 *   Cambridge University Press. https://numerical.recipes/book.html
 *
 */
struct ridder : public bisection {
    using bisection::bisection;

    template <typename F>
    bool iterate(F&& fun)
    {
        proposal.x = get_midpoint(left, right);
        proposal.y = fun(proposal.x);

        double d = proposal.x - left.x;
        double a = proposal.y / left.y;
        double b = right.y / left.y;
        double denom = a * a - b;
        if (is_zero(denom)) {
            // is this state possible ? fall back to bisection
            flag = Flag::division_by_zero;
            return false;
        }

        proposal.x += d * a / std::sqrt(denom);
        proposal.y = fun(proposal.x);
        update_bracket();
        return true;
    }
};

/**
 * @brief Not a method. Illinois-type methods use the update bracket defined here.
 * See the documentation of the `illinois` method for details. Will act as
 * bisection method if called.
 */
struct illinois_type : public bisection {
    using bisection::bisection;
    void update_bracket(double gamma)
    {
        if (opposite_signs(proposal.y, right.y)) {
            left = right;
        } else {
            left.y *= gamma;
        }
        right = proposal;
    }
};

/**
 * @brief The "Illinois" method. Provide a function object and an initial
 * valid bracket. Bracket is valid if the sign of the function differs at
 * the end points.
 *
 * @details The "Illinois" method, so called because it was developed at
 * the University of Illinois in the 1950s, is effectively the same as
 * "regula falsi" except it catches the failure mode of the regula falsi.
 *
 * Regula falsi is slow if the same end point is retained twice in a row.
 * In the Illinois method, if the same end point is retained twice then
 * the value of the function is reduced by half for computing the next
 * iterate.
 *
 * Let \f$ [x_{n-1}, x_n]\f$ be a bracket, where the "right" side is most recent
 * endpoint. It is possible to have $x_n < x_{n-1}$ even though intervals are not
 * normally written that way. A new point is generated by the secant formula:
 *
 * \f[ x_{n+1} = \frac{f(x_n) x_{n-1} - f(x_{n-1}) x_{n}}{f(x_n) - f(x_{n-1})} \f]
 *
 * If the sign of \f$f(x_{n+1})\f$ is the same as \f$f(x_{n-1})\f$ then we form
 * a new bracket \f$ [x_n , x_{n+1} ] \f$ as normal in "regula falsi". Otherwise,
 * if the sign of \f$f(x_{n+1})\f$ is the same as \f$f(x_{n})\f$ then we build a
 * new bracket \f$[x_{n-1} , x_{n+1} ]\f$ but for the next iteration, we use the
 * value \f$f(x_{n-1})/2 \f$ instead of \f$f(x_{n-1})\f$. In practice, the saved
 * result `left.y = f(left.x)` is simply halved: ` left.y /= 2 `.
 *
 * The illinois method is the simplest of a family of methods, which
 * all rescale the value of the function `f` at the retained endpoint when
 * that endpoint is retained for a second iteration. `left.y *= gamma`.
 *
 * This method is basically always faster than "regula falsi" and has
 * robustness similar to the bisection method.
 *
 * References:
 * - Ford, J. A. (1995). "Improved Illinois-type methods for the solution
 *   of nonlinear equations." Technical Report, University of Essex Press.
 *
 * - Dowell, M.; Jarratt, P. (1971). "A modified regula falsi method for
 *   computing the root of an equation". BIT. 11 (2): 168–174.
 *   doi:10.1007/BF01934364
 */
struct illinois : public illinois_type {
    using illinois_type::illinois_type;
    // does not preserve left and right. Treats right as best guess.
    template <typename F>
    bool iterate(F&& fun)
    {
        proposal.x = get_secant_update(left, right);
        proposal.y = fun(proposal.x);
        update_bracket(0.5);
        return true;
    }
};

/**
 * @brief The "pegasus" method. An Illinois-type bracketing method. Provide a
 * valid bracket.
 *
 * @details The `pegasus` has the same general idea as the `illinois` method but
 * uses a different update for when the same endpoint is retained twice. See the
 * `illinois` method for details.
 *
 * The scaling factor is the ratio:
 *
 * \f[ \gamma = \frac{f(x_n)}{f(x_n) + f(x_{n+1})} \f]
 *
 * which is always positive because the update only occurs when \f$f(x_n)\f$ and
 * \f$f(x_{n+1})\f$ have the same same sign.
 *
 * It is slightly faster than the `illinois` method in numerical tests.
 * See references for additional details.
 *
 * References:
 * - Ford, J. A. (1995). "Improved Illinois-type methods for the solution
 *   of nonlinear equations." Technical Report, University of Essex Press.
 *
 * - Dowell, M., Jarratt, P. The “Pegasus” method for computing the root of an
 *   equation. BIT 12, 503–508 (1972). https://doi.org/10.1007/BF01932959
 */
struct pegasus : public illinois_type {
    using illinois_type::illinois_type;
    // does not preserve left and right. Treats right as best guess.
    template <typename F>
    bool iterate(F&& fun)
    {
        proposal.x = get_secant_update(left, right);
        proposal.y = fun(proposal.x);
        update_bracket(right.y / (right.y + proposal.y));
        return true;
    }
};

/**
 * @brief The "Anderson-Björck" method. An Illinois-type bracketing method. Provide a
 * valid bracket.
 *
 * @details The same general idea as the `illinois` method but using a different
 * update for when the same endpoint is retained twice. See the `illinois` method
 * for details.
 *
 *  The scaling factor is the ratio of the divided differences:
 * \f[ \gamma = \frac{f[x_{n+1}, x_n]}{f[x_n, x_{n-1}]} \f]
 * Whether the divided difference is:
 *
 * \f[ f[x_{n+1}, x_n] = \frac{f(x_{n+1}) - f(x_n)}{x_{n+1} - x_n}\f]
 *
 * See reference for details.
 *
 * References:
 * - Ford, J. A. (1995). "Improved Illinois-type methods for the solution
 *   of nonlinear equations." Technical Report, University of Essex Press.
 *
 */
struct anderson_bjorck : public illinois_type {
    using illinois_type::illinois_type;
    // does not preserve left and right. Treats right as best guess.
    template <typename F>
    bool iterate(F&& fun)
    {
        proposal.x = get_secant_update(left, right);
        proposal.y = fun(proposal.x);
        double m0 = (proposal.y - right.y) / (proposal.x - right.x);
        double m1 = (right.y - left.y) / (right.x - left.x);
        update_bracket(m0 / m1);
        return true;
    }
};
}  // namespace root_finding
#endif
