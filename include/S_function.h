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
    // F(x) = 1.0 + (alpha * x + beta * x * x) * exp(-x/delta)
    return 1.0 + (alpha * x + beta * x * x) * exp(-x / delta);
}

/**
* First derivative of F(x).
*/
inline double Fp(double x, double alpha, double beta, double gamma, double delta) {
    // Let g(x) = alpha*x + beta*x^2
    // g'(x) = alpha + 2*beta*x
    // F'(x) = g'(x)*exp(-x/delta) + g(x)*exp(-x/delta)*(-1/delta)
    //        = exp(-x/delta) * [g'(x) - g(x)/delta]
    double g = alpha * x + beta * x * x;
    double gp = alpha + 2.0 * beta * x;
    return exp(-x / delta) * (gp - g / delta);
}

/**
* Second derivative of F(x).
*/
inline double Fpp(double x, double alpha, double beta, double gamma, double delta) {
    // F''(x) = exp(-x/delta) * [g''(x) - 2*g'(x)/delta + g(x)/delta^2]
    // where g''(x) = 2*beta
    double g = alpha * x + beta * x * x;
    double gp = alpha + 2.0 * beta * x;
    double gpp = 2.0 * beta;
    return exp(-x / delta) * (gpp - 2.0 * gp / delta + g / (delta * delta));
}
