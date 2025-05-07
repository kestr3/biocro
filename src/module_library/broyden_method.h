#ifndef BROYDEN_METHOD_H
#define BROYDEN_METHOD_H

#include "root_onedim.h"
#include <array>

template <typename T, size_t dim>
struct broyden {
    using vec_t = typename std::array<T, dim>;

    struct result_t {
        vec_t zero;
        vec_t residual;
        size_t iteration;
        root_algorithm::Flag flag;
    };

    template <typename F>
    result_t solve(F&& fun, vec_t x0)
    {
    }
};

#endif
