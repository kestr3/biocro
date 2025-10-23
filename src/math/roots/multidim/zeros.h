#ifndef ROOT_MULTIDIM_H
#define ROOT_MULTIDIM_H

#include <array>
#include <cmath>    // std::sqrt, std::isnan
#include <numeric>  // std::inner_product

// #include <iostream>

namespace root_multidim
{
// May be replaced later with BOOST's uBLAS code.
// vector type
template <size_t dim>
using vec_t = std::array<double, dim>;

// matrix type
template <size_t dim>
using mat_t = std::array<vec_t<dim>, dim>;

// termination states
enum class Flag {
    residual_zero,
    delta_x_zero,
    zero_is_nonfinite,
    max_iterations,

    function_is_nonfinite,

};

template <size_t dim>
struct result_t {
    vec_t<dim> zero;
    vec_t<dim> residual;
    size_t iteration;
    Flag flag;
};

template <size_t dim, typename Method>
struct zero_finding_method {
    // declare types for solving
    using vec_t = vec_t<dim>;
    using mat_t = mat_t<dim>;
    using result_t = result_t<dim>;

    zero_finding_method(size_t max_iter, double abs_tol, double rel_tol) : max_iterations{max_iter}, _abs_tol{abs_tol}, _rel_tol{rel_tol} {}
    // zero_finding_method(size_t max_iter, const double& abs_tol, const double& rel_tol) : max_iterations{max_iter}, _abs_tol{abs_tol}, _rel_tol{rel_tol} {}

    size_t max_iterations = 100;
    double _abs_tol = 1e-12;
    double _rel_tol = 1e-12;
    bool is_valid;
    Flag flag;

    template <typename F, typename... Args>
    result_t solve(F&& fun, Args&&... args)
    {
        is_valid = static_cast<Method*>(this)->initialize(std::forward<F>(fun), std::forward<Args>(args)...);
        for (size_t i = 0; i <= max_iterations; ++i) {
            if (!is_valid) {
                return make_result(i);
            }

            is_valid = static_cast<Method*>(this)->iterate(std::forward<F>(fun));

            if (is_valid) {
                is_valid = static_cast<Method*>(this)->has_converged();
            }
        }
        flag = Flag::max_iterations;
        return make_result(max_iterations);
    }

    template <typename F, typename... Args>
    inline result_t operator()(F&& fun, Args&&... args)
    {
        return solve(std::forward<F>(fun), std::forward<Args>(args)...);
    }

   protected:
    result_t make_result(size_t i)
    {
        result_t out;
        out.zero = static_cast<Method*>(this)->zero();
        out.residual = static_cast<Method*>(this)->residual();
        out.iteration = i;
        return out;
    }
    // All methods require tolerance-based floating point number equality tests.
    // Are two floating point numbers equal?
    // Equality is lax if both `x` and `y` are big.
    inline bool is_close(double x, double y) const
    {
        double norm = std::min(std::abs(x), std::abs(y));
        return std::abs(x - y) <= std::max(_abs_tol, _rel_tol * norm);
    }

    // Is a floating point number zero?
    inline bool is_zero(double x) const
    {
        return std::abs(x) <= _abs_tol;
    }

    // floating point errors can increase with dimension
    inline bool is_zero(const vec_t& y, const vec_t& x) const
    {
        double ysq = std::inner_product(y.cbegin(), y.cend(), y.cbegin(), 0.0);
        double xsq = std::inner_product(x.cbegin(), x.cend(), x.cbegin(), 0.0);
        return std::sqrt(ysq) < _abs_tol + _rel_tol * std::sqrt(xsq);
    }

    inline bool is_nan(const vec_t& x) const
    {
        for (const double& v : x) {
            if (std::isnan(v)) return true;
        }
        return false;
    }
};

// template <size_t dim>
// vec_t<dim>& operator+=(vec_t<dim>& x, const vec_t<dim>& y)
// {
//     for (size_t i = 0; i < dim; ++i) {
//         x[i] += y[i];
//     }
//     return x;
// }

// template <size_t dim>
// vec_t<dim>& operator-=(vec_t<dim>& x, const vec_t<dim>& y)
// {
//     for (size_t i = 0; i < dim; ++i) {
//         x[i] -= y[i];
//     }
//     return x;
// }

// template <size_t dim>
// vec_t<dim> operator+(const vec_t<dim>& x, const vec_t<dim>& y)
// {
//     vec_t<dim> out;
//     for (size_t i = 0; i < dim; ++i) {
//         out[i] = x[i] + y[i];
//     }
//     return x;
// }

// template <size_t dim>
// vec_t<dim> operator-(const vec_t<dim>& x, const vec_t<dim>& y)
// {
//     vec_t<dim> out;
//     for (size_t i = 0; i < dim; ++i) {
//         out[i] = x[i] - y[i];
//     }
//     return out;
// }

template <size_t dim>
inline double dot(const vec_t<dim>& u, const vec_t<dim>& v)
{
    return std::inner_product(u.cbegin(), u.cend(), v.cbegin(), 0.0);
}

template <size_t dim>
vec_t<dim> dot(const mat_t<dim>& A, const vec_t<dim>& v)
{
    vec_t<dim> out;
    for (size_t i = 0; i < v.size(); ++i) {
        out[i] = dot(A[i], v);
    }
    return out;
}

// // transpose

// template <size_t dim>
// vec_t<dim> dot(const vec_t<dim>& v, const std::array<vec_t<dim>, dim>& A)
// {
//     vec_t<dim> out;
//     out.fill(0);
//     for (size_t i = 0; i < dim; ++i) {
//         for (size_t j = 0; j < dim; ++j) {
//             out[j] += v[i] * A[i][j];
//         }
//     }
//     return out;
// }

}  // namespace root_multidim
#endif
