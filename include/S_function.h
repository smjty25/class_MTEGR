#include <math.h>

/*
static inline double smooth_step(double x, double x0, double x1, double f0, double f1) {
    double y = (x - x0) / (x1 - x0);  // y in [0,1]
    // Use (1 - cos(pi*y))/2 which smoothly goes from 0 to 1
    double poly = 0.5 * (1.0 - cos(M_PI * y));
    return f0 + (f1 - f0) * poly;
}

static inline double smooth_step_p(double x, double x0, double x1, double f0, double f1) {
    double y = (x - x0) / (x1 - x0);
    double factor = (f1 - f0) / (x1 - x0);
    // d/dx[(1-cos(pi*y))/2] = (pi/2)*sin(pi*y) * dy/dx
    return factor * 0.5 * M_PI * sin(M_PI * y);
}

static inline double smooth_step_pp(double x, double x0, double x1, double f0, double f1) {
    double y = (x - x0) / (x1 - x0);
    double factor = (f1 - f0) / ((x1 - x0) * (x1 - x0));
    // d²/dx² = (pi²/2)*cos(pi*y) * (dy/dx)²
    return factor * 0.5 * M_PI * M_PI * cos(M_PI * y);
}

inline double F(double x, double alpha, double beta, double gamma, double delta) {
    double left = gamma;
    double center = gamma + delta/2.0;
    double right = gamma + delta;
    if (x <= left) {
        return 1.0;
    } else if (x < center) {
        return smooth_step(x, left, center, 1.0, beta);
    } else if (x < right) {
        return smooth_step(x, center, right, beta, alpha);
    } else {
        return alpha;
    }
}

inline double Fp(double x, double alpha, double beta, double gamma, double delta) {
    double left = gamma;
    double center = gamma + delta/2.0;
    double right = gamma + delta;
    if (x <= left || x >= right) {
        return 0.0;               // constant regions
    } else if (x < center) {
        return smooth_step_p(x, left, center, 1.0, beta);
    } else {  // x between gamma and gamma+delta (including gamma)
        return smooth_step_p(x, center, right, beta, alpha);
    }
}

inline double Fpp(double x, double alpha, double beta, double gamma, double delta) {
    double left = gamma;
    double center = gamma + delta/2.0;
    double right = gamma + delta;
    if (x <= left || x >= right) {
        return 0.0;               // constant regions
    } else if (x < center) {
        return smooth_step_pp(x, left, center, 1.0, beta);
    } else {  // x between gamma and gamma+delta (including gamma)
        return smooth_step_pp(x, center, right, beta, alpha);
    }
}
*/

inline double F(double x, double alpha, double beta, double gamma, double delta) {
    // F(x) = 1.0 + (beta * x * x * exp(2.0 - 2.0 * x / delta)) / (delta * delta)
    return 1.0 + (beta * x * x * exp(2.0 - 2.0 * x / delta)) / (delta * delta);
}

/**
* First derivative of F(x).
*/
inline double Fp(double x, double alpha, double beta, double gamma, double delta) {
    // Let h(x) = beta * x^2 * exp(2 - 2x/delta) / delta^2
    // h'(x) = (beta/delta^2) * [2x * exp(2 - 2x/delta) + x^2 * exp(2 - 2x/delta) * (-2/delta)]
    //        = (beta * exp(2 - 2x/delta) / delta^2) * [2x - 2x^2/delta]
    //        = (2 * beta * x * exp(2 - 2x/delta) / delta^2) * (1 - x/delta)
    double exp_term = exp(2.0 - 2.0 * x / delta);
    return (2.0 * beta * x * exp_term / (delta * delta)) * (1.0 - x / delta);
}

/**
* Second derivative of F(x).
*/
inline double Fpp(double x, double alpha, double beta, double gamma, double delta) {
    // h''(x) = (2*beta/delta^2) * exp(2 - 2x/delta) * [1 - 4x/delta + 2x^2/delta^2]
    // (derived from product rule on h'(x))
    double exp_term = exp(2.0 - 2.0 * x / delta);
    double x_over_delta = x / delta;
    return (2.0 * beta * exp_term / (delta * delta)) * 
           (1.0 - 4.0 * x_over_delta + 2.0 * x_over_delta * x_over_delta);
}