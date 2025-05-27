#ifndef BROYDEN_METHOD_H
#define BROYDEN_METHOD_H

// #include "root_onedim.h"
#include <array>
#include <cmath>
#include <functional>

template <typename T, size_t dim>
struct broyden {
    using vec_t = typename std::array<T, dim>;
    using mat_t = typename std::array<std::array<T, dim>, dim>;

    struct result_t {
        vec_t zero;
        vec_t residual;
        size_t iteration;
    };

    size_t max_iterations;
    double _abs_tol;
    double _rel_tol;

    template <typename F>
    result_t solve(F&& fun, vec_t x)
    {
        vec_t y;
        mat_t inv_jac;
        fun(y, x);
        set_identity(inv_jac);
        vec_t delta_x;
        vec_t delta_y;

        for (size_t i = 0; i <= max_iterations; ++i) {
            update_delta_x(delta_x, y, inv_jac);
            update_x(x, delta_x);
            fun(delta_y, x);
            update_delta_y(delta_y, y);
            if (is_zero(y)) {
                return result_t{x, y, i};
            }
            update_inv_jac(inv_jac, delta_y, delta_x);
        }
        return result_t{x, y, max_iterations};
    }

   private:
    void set_identity(mat_t& A)
    {
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                A[i][j] = 0;
                if (i == j)
                    A[i][j] = 1;
            }
        }
    }

    void update_delta_x(
        vec_t& delta_x, const vec_t& y, const mat_t& inv_jac)
    {
        // delta_x = - inv_jac * y
        for (size_t i = 0; i < dim; ++i) {
            delta_x[i] = -inv_jac[i][0] * y[0];
            for (size_t j = 1; j < dim; ++j) {
                delta_x[i] -= inv_jac[i][j] * y[j];
            }
        }
    }

    void update_x(
        vec_t& x, const vec_t& delta_x)
    {
        // x += delta_x
        for (size_t i = 0; i < dim; ++i) {
            x[i] += delta_x[i];
        }
    }

    void update_delta_y(vec_t& delta_y, vec_t& y)
    {
        for (size_t i = 0; i < dim; ++i) {
            std::swap(delta_y[i], y[i]);
            delta_y[i] = y[i] - delta_y[i];
        }
    }

    void update_inv_jac(mat_t& inv_jac, const vec_t& delta_y, const vec_t& delta_x)
    {
        // out = inv_jac + (delta_x - inv_jac * delta_y) / ( t(delta_x) * inv_jac * delta_y ) * t(delta_x) * inv_jac
        vec_t a;
        vec_t b;
        b.fill(0);
        T c = 0;
        for (size_t i = 0; i < dim; ++i) {
            a[i] = delta_x[i];
            for (size_t j = 0; j < dim; ++j) {
                a[i] -= inv_jac[i][j] * delta_y[j];
                b[j] += delta_x[i] * inv_jac[i][j];
                c += delta_x[i] * inv_jac[i][j] * delta_y[j];
            }
        }

        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                inv_jac[i][j] += a[i] * b[j] / c;
            }
        }
    }

    bool is_zero(const vec_t& y)
    {
        T sq = 0;
        for (size_t i = 0; i < dim; ++i) {
            sq += y[i] * y[i];
        }
        return std::sqrt(sq) < _abs_tol;
    }
};

#endif
