#ifndef ROOT_BROYDEN_H
#define ROOT_BROYDEN_H

#include <array>
#include <cmath>    // std::sqrt, std::isnan
#include <numeric>  // std::inner_product
#include "zeros.h"

namespace root_multidim
{

/** @brief Broyden's method.
 *
 *  @details
 */
template <size_t dim>
struct broyden : public zero_finding_method<dim, broyden<dim>> {
    using zero_finding_method<dim, broyden<dim>>::zero_finding_method;
    using vec_t = vec_t<dim>;
    using mat_t = mat_t<dim>;

    vec_t _zero;
    vec_t _residual;
    vec_t delta_x;
    vec_t delta_y;
    mat_t inv_jac;

    vec_t _tmp_a;
    vec_t _tmp_b;
    double _tmp_c;

    template <typename F>
    bool initialize(F&& fun, const vec_t& guess)
    {
        _zero = guess;
        _residual = fun(guess);
        // mat_t jac = approx_jac(fun, guess);
        inv_jac = identity();

        return true;
    }

    template <typename F>
    bool iterate(F&& fun)
    {
        /*
        UPDATE FORMULA
            x -> x + delta_x ;
            where
                delta_x = -inv_jac(y);
                y = fun(x);
                delta_y = fun(x + delta_x) - fun(x);
            inv_jac ->  inv_jac + outer(a, b) / c
            where
                a = delta_x - dot(inv_jac , delta_y)
                b = dot(delta_x, inv_jac)
                c = dot(delta_x, dot(inv_jac, delta_y))

            `update_x` does
                delta_x = dot(inv_jac, -y);
                x += delta_x;

            need to save current/old `residual` to compute new `delta_y`
            don't need `delta_y` so overwrite it with new `residual`
            delta_y = fun(_zero)
            `update_y`
                std::swap(_residual, delta_y);
                // _residual is now new `residual`
                // delta_y is now old `residual`
                delta_y = _residual - y_old

            Could check convergence now

            `update_inv_jac` computes
                a = delta_x - dot(inv_jac , delta_y)
                b = dot(delta_x, inv_jac)
                c = dot(delta_x, dot(inv_jac, delta_y))

        */
        update_x();
        delta_y = fun(_zero);
        update_y();
        update_inv_jac();
        return true;
    }

    bool has_converged()
    {
        if (this->is_zero(_residual, _zero)) {
            this->flag = Flag::residual_zero;
            return true;
        }

        if (this->is_zero(delta_x, _zero)) {
            this->flag = Flag::delta_x_zero;
            return true;
        }

        return false;
    }

    vec_t residual() const
    {
        return _residual;
    }

    vec_t zero() const
    {
        return _zero;
    }

    // vec_t& residual()
    // {
    //     return _residual;
    // }

    // vec_t& zero()
    // {
    //     return _zero;
    // }

    // template <typename F>
    // result_t solve(F&& fun, const vec_t& guess)
    // {
    //     result.zero = guess;
    //     result.residual = fun(result.zero);
    //     mat_t jac = approx_jac(fun, result.zero);
    //     mat_t inv_jac = invert(jac);

    //     vec_t delta_x;
    //     vec_t delta_y;

    //     // temporaries for updating `inv_jac`

    //     for (size_t i = 0; i <= max_iterations; ++i) {
    //         update_x(result.zero, delta_x, result.residual, inv_jac);
    //         // print(std::cout, result.zero);
    //         delta_y = fun(result.zero);
    //         update_y(result.residual, delta_y);
    //         if (is_zero(result.residual, result.zero) || is_nan(result.zero)) {
    //             result.iteration = i;
    //             return result;
    //         }
    //         inv_jac = update_inv_jac(inv_jac, delta_y, delta_x);
    //     }
    //     result.iteration = max_iterations;
    //     return result;
    // }

    inline void update_x()
    {
        // these formulae are computed in a single for-loop pass
        // delta_x = dot(inv_jac,  -residual)
        // zero += delta_x
        for (size_t i = 0; i < dim; ++i) {
            delta_x[i] = 0;
            for (size_t j = 0; j < dim; ++j) {
                delta_x[i] += inv_jac[i][j] * -_residual[j];
            }
            _zero[i] += delta_x[i];
        }
    }

    inline void update_y()
    {
        // this is computed in a single for-loop pass
        for (size_t i = 0; i < dim; ++i) {
            std::swap(delta_y[i], _residual[i]);     // swap values so that y[i] is new value, and delta_y[i] is old
            delta_y[i] = _residual[i] - delta_y[i];  //
        }
    }

    void update_inv_jac()
    {
        /*
        inv_jac ->  inv_jac + outer(a, b) / c
          where
            a = delta_x - dot(inv_jac , delta_y)
            b = dot(delta_x, inv_jac)
            c = dot(delta_x, dot(inv_jac, delta_y))


        */

        _tmp_c = 0;
        _tmp_b.fill(0);

        for (size_t i = 0; i < dim; ++i) {
            _tmp_a[i] = delta_x[i];
            for (size_t j = 0; j < dim; ++j) {
                _tmp_a[i] -= inv_jac[i][j] * delta_y[j];
                _tmp_b[j] += delta_x[i] * inv_jac[i][j];
                _tmp_c += delta_x[i] * inv_jac[i][j] * delta_y[j];
            }
        }

        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                inv_jac[i][j] += _tmp_a[i] * _tmp_b[j] / _tmp_c;
            }
        }
    }

    // void print(std::ostream& os, const vec_t& v)
    // {
    //     os << "[ ";
    //     for (size_t i = 0; i < dim; ++i) {
    //         os << v[i];
    //         if (i < dim - 1) {
    //             os << ", ";
    //         }
    //     }
    //     os << " ]\n";
    // }

    // void print(std::ostream& os, const mat_t& v)
    // {
    //     os << "[ ";
    //     for (size_t i = 0; i < dim; ++i) {
    //         os << v[i];
    //         if (i < dim - 1) {
    //             os << ", ";
    //         }
    //     }
    //     os << " ]\n";
    // }

   private:
    template <typename F>
    mat_t approx_jac(F&& fun, const vec_t& x)
    {
        vec_t x0 = x;
        mat_t jac;
        vec_t y0 = fun(x);
        vec_t y1;
        double constexpr eps = 1e-10;
        for (size_t i = 0; i < dim; ++i) {
            x0[i] += eps;
            if (i > 0) {
                x0[i - 1] -= eps;
            }
            y1 = fun(x0);
            for (size_t j = 0; j < dim; ++j) {
                jac[j][i] = (y1[j] - y0[j]) / eps;
            }
        }

        return jac;
    }

    inline mat_t invert(const mat_t& A)
    {
        mat_t out = identity();
        double b;
        // gauss seidel method
        for (size_t iteration = 0; iteration < 100; ++iteration) {
            for (size_t i = 0; i < dim; ++i) {
                for (size_t k = 0; k < dim; ++k) {
                    b = i == k ? 1 : 0;
                    out[i][k] = b;
                    for (size_t j = 0; j < dim; ++j) {
                        if (j != i)
                            out[i][k] -= A[i][j] * out[j][k];
                    }
                    out[i][k] /= A[i][i];

                    if (std::isnan(out[i][k])) {
                        return identity();
                    }
                }
            }
        }
        return out;
    }

    inline mat_t identity()
    {
        mat_t out;
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                out[i][j] = 0;
                if (i == j)
                    out[i][j] = 1;
            }
        }
        return out;
    }
};

}  // namespace root_multidim
#endif
