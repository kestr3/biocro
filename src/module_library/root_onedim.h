#ifndef ROOTS_ONEDIM_H
#define ROOTS_ONEDIM_H




#include <cmath>
#include <limits>
#include <iterator>



template<typename F, typename T>
class function_iterator {
private:
    const F func;
    T current;
    size_t iteration;

public:
    // Iterator traits
    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = const T&;


    function_iterator(const F function, T initial)
        : func{function}, current{initial}, iteration{0} {}

    function_iterator(const F function, T initial, size_t max_iter) : func{function}, current{}, iteration{max_iter} {}

    reference operator*() const {
        return current;
    }

    function_iterator& operator++() {
        current = func(current);
        ++iteration;
        return *this;
    }

    // Post-increment == Pre-increment
    function_iterator operator++(int) {
        function_iterator tmp = *this;
        ++(*this);
        return tmp;
    }

    bool operator==(const function_iterator& other) const {
        return (iteration == other.iteration);
    }

    bool operator!=(const function_iterator& other) const {
        return !(*this == other);
    }

};

template<typename F, typename T>
class function_sequence {

private:
    const F& func;
    const T x0;
    size_t max_iterations;

public:
    function_sequence(const F& function, T initial,  size_t max_iter) :
        func{function}, x0{initial}, max_iterations{max_iter} {}

    function_iterator<F , T> begin() {
        return function_iterator<F , T>(func, x0);
    }

    function_iterator<F, T> end() {
        return function_iterator<F, T>(func, x0, max_iterations);
    }

};

template<typename F, typename T>
function_sequence<F, T> iterate(const F& function, T initial, size_t max_iter){
    return function_sequence<F, T>(function, initial, max_iter);
}

    template<typename Method>
    class iterator {

        Method::state current;
        size_t iteration;

    public:
        // Iterator traits
        using iterator_category = std::input_iterator_tag;
        using value_type = Method::state;
        using difference_type = std::ptrdiff_t;
        using pointer = const Method::state*;
        using reference = const Method::state&;


        iterator(Method::state& initial) : current{initial}, iteration{0} {}
        iterator(size_t max_iter) : current{}, iteration{max_iter} {}
        ~iterator()

        reference operator*() const {
            return current;
        }

        iterator& operator++() {
            Method::update(current);
            ++iteration;
            return *this;
        }

        // Post-increment == Pre-increment
        iterator operator++(int) {
            function_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const iterator& other) const {
            return (iteration == other.iteration);
        }

        bool operator!=(const iterator& other) const {
            return !(*this == other);
        }
    };
/*
methods
*/

// Helper function declarations
inline bool is_close(double x, double y, double tol, double rtol);
inline bool is_zero(double x, double tol);
inline bool same_signs(double x, double y);

template<typename F, typename Method>
class root_solver {

    const F& fun;
    Method method;

    size_t max_iterations;
    double _abs_tol;
    double _rel_tol;

    template<typename... Args>
    root_solver(
        const F& function, Args&&... args,
        size_t max_iter, double abs_tol, double rel_tol) :
        fun{function}, method{std::forward<Args>(args)...},
        max_iterations{max_iter}, _abs_tol{abs_tol}, _rel_tol{rel_tol} {}




};

template<typename Method>
inline bool has_converged(Method method, double&& tol){
    return is_zero(method.root(), tol);
}

template<typename Method>
struct root_finder {

    size_t max_iterations = 100;
    double _abs_tol = std::numeric_limits<double>::epsilon();
    double _rel_tol = 2*std::numeric_limits<double>::epsilon();

    root_finder(size_t max_iter, double abs_tol, double rel_tol) :
        max_iterations{max_iter}, _abs_tol{abs_tol}, _rel_tol{rel_tol} {}

    template<typename F, typename ...Args>
    root_result solve(const F& function, Args&& ...args){
        Method method{function, std::forward<Args>(args)...};
        for (size_t i = 0; i != max_iterations; ++i){
            if (has_converged(method, _abs_tol)){
                return root_result{method.root(), method.residual(), i};
            };
            method.iterate();
        }
        return root_result{method.root(), method.residual(), max_iterations};
    }
};


struct secant {

    struct state {
        double x0;
        double y0;
        double x1;
        double y1;
    };

    template<typename F>
    state initialize(F&& fun, double x0, double x1){
        double y0 = fun(x0);
        double y1 = fun(x1);
        return state{x0, y0, x1, y1};
    }

    template<typename F>
    void update(F&& fun , state& s){
        double secant_slope = (s.y1 - s.y0) / (s.x1 - s.x0);
        double x2 = s.x1 - s.y1 / secant_slope;
        double y2 = fun(x2);
        x0 = x1;
        x1 = x2;
        y0 = y1;
        y1 = y2;
    }

    inline bool has_converged(state& s){
        if is_close(s.x0, s.x1) {
            return true;
        }
    }

    inline double residual(){
        return y1;
    }

};
template<typename fun>
struct newton {

    const F& fun;

    struct state {
        double x;
        double y;
    };


    newton(const F& function) : fun{function} {}

    state initialize(double x){
        state{x, fun(x)};
    }

    void update(state& s){
        s.x = s.x - s.y / fun.jac(s.x);
        s.y = fun(s.x);
    }



};



template <typename F>
root_result find_root_secant_method(
    F const fun,
    double x0,
    double x1,
    size_t const max_iter,
    double const tol,
    double const rtol)
{

    // Check inputs for trivial solutions
    double y0 = fun(x0);
    if (is_zero(y0, tol)) {
        return root_result{x0, y0, 0};
    }

    // Guesses should be different
    if (is_close(x1, x0, tol, rtol)) {
        throw std::runtime_error(
            "Within tolerance, x0 == x1. Initial guesses should be distinct.");
    }

    // Check second solution
    double y1 = fun(x1);
    if (is_zero(y1, tol)) {
        return root_result{x1, y1, 0};
    }

    for (size_t n = 1; n <= max_iter; ++n) {


        double secant_slope = (y1 - y0) / (x1 - x0);

        // Terminate if division by zero encountered
        if (is_zero(secant_slope, tol)) {
            return root_result{x1, y1, n};
        }

        // Secant update formula
        double x2 = x1 - y1 / secant_slope;
        double y2 = fun(x2);


        // Test convergence; if y2 == 0, then x2 is a root.

        if (is_zero(y2, tol)) {
            return root_result{x2, y2, n};
        }
        // If no improvement in guesses, then terminate; prevents x1 == x0.
        if (is_close(x1, x2, tol, rtol)) {
            return root_result{x2, y2, n};
        }

        // Update
        x0 = x1;
        x1 = x2;
        y0 = y1;
        y1 = y2;
    }

    return root_result{x1, y1, max_iter};
}

template<typename F, typename G>
struct newton {

    const F& _fun;
    const G& _jac;

    double x_previous;
    double y_previous;
    double x_current;
    double y_current;

    newton(const F& fun, const G& jac) : _fun{fun}, _jac{jac} {}

    void initialize(double x){
        _x = x;
        _y = _fun(x);
    }

    void update(){
        _x = _x - _y/_jac(_x);
        _y = _fun(_x);
    }

};

template<typename F, typename G>
root_result newton(const F& fun, const G& jac, double x0, size_t max_iter, double tol){
    double x = x0;
    double y = fun(x);
    if (is_zero(y, tol)){
        return root_result{x, y, 0};
    }

    for (size_t i = 1; i <= max_iter; ++i){
        double x_new = x - y/jac(x);
        y = fun(x_new);

        if(is_close(x_new, x, tol, tol)){
            return root_result{x_new, y, i};
        }

        if (is_zero(y, tol)){
            return root_result{x_new, y, i};
        }

        x = x_new;
    }

    return root_result{x, y, max_iter};
}

template<typename F>
root_result bisection(const F& fun, double a, double b, size_t max_iter, double tol, double rtol){

    double ya = fun(a);
    if (is_zero(ya, tol)){
        return root_result{a, ya, 0};
    }

    double yb = fun(b);
    if (is_zero(yb, tol)){
        return root_result{b, yb, 0};
    }

    if(same_signs(ya, yb)){
        throw std::runtime_error("Signs at endpoints must be different.");
    }

    double x;
    double y;
    for (size_t i = 1; i <= max_iter; ++i){
        x = (a + b)/2;
        y = fun(x);

        if (is_zero(y, tol)){
            return root_result{x, y, i};
        }

        if (is_close(a, x, tol, rtol)){
            return root_result{x, y, i};
        }

        if (is_close(b, x, tol, rtol)){
            return root_result{x, y, i};
        }

        if (same_signs(y, ya)){
            a = x;
            ya = y;
        }
        else {
            b = x;
            yb = y;
        }
    }

    x = (a + b)/2;
    y = fun(x);

    return root_result{x, y, max_iter};
}



// Helper function definitions
inline bool is_close(double x, double y, double tol, double rtol)
{
    return std::abs(x - y) <=tol + rtol * std::abs(y);
}

inline bool is_zero(double x, double tol)
{
    return std::abs(x) <= tol;
}

inline bool same_signs(double x, double y){
    return x * y > 0;
}

inline bool opposite_signs(double x, double y){
    return x * y < 0;
}

// Test functions

double testfun(double x){
    if(x < 0 )  return -0.25;
    return  (std::pow(x, 4) - 1)/4;
}

double testfunjac(double x){
    if(x < 0 )  return 0.;
    return std::pow(x, 3);
}

template<size_t max_iter>
std::array<double, max_iter> test_sequence(double x0){
    std::array<double, max_iter> out;
    out[0] = x0;
    for (size_t i= 0; i != max_iter - 1; ++i){
        out[i+1] = out[i] - testfun(out[i])/testfunjac(out[i]);
    }
    return out;
}






#endif
