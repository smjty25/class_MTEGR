#include <math.h>
/**
 * Smooth step function using cosine (C-infinity continuous)
 * Maps x from [x0,x1] to [f0,f1]
 */
static inline double smooth_step(double x, double x0, double x1, double f0, double f1) {
    double y = (x - x0) / (x1 - x0);  // y in [0,1]
    // Use (1 - cos(pi*y))/2 which smoothly goes from 0 to 1
    double poly = 0.5 * (1.0 - cos(M_PI * y));
    return f0 + (f1 - f0) * poly;
}

/**
 * First derivative of cosine smooth_step with respect to x
 */
static inline double smooth_step_p(double x, double x0, double x1, double f0, double f1) {
    double y = (x - x0) / (x1 - x0);
    double factor = (f1 - f0) / (x1 - x0);
    // d/dx[(1-cos(pi*y))/2] = (pi/2)*sin(pi*y) * dy/dx
    return factor * 0.5 * M_PI * sin(M_PI * y);
}

/**
 * Second derivative of cosine smooth_step with respect to x
 */
static inline double smooth_step_pp(double x, double x0, double x1, double f0, double f1) {
    double y = (x - x0) / (x1 - x0);
    double factor = (f1 - f0) / ((x1 - x0) * (x1 - x0));
    // d²/dx² = (pi²/2)*cos(pi*y) * (dy/dx)²
    return factor * 0.5 * M_PI * M_PI * cos(M_PI * y);
}
/**
 * Piecewise function F(x) defined as:
 *   F(x) = 1                         for x <= gamma - delta
 *          smooth_step(x, gamma-delta, gamma, 1, beta)   for gamma-delta < x < gamma
 *          smooth_step(x, gamma, gamma+delta, beta, alpha) for gamma <= x < gamma+delta
 *          alpha                      for x >= gamma+delta
 */
inline double F(double x, double alpha, double beta, double gamma, double delta) {
    double left = gamma - delta;
    double right = gamma + delta;

    if (x <= left) {
        return 1.0;
    } else if (x < gamma) {
        return smooth_step(x, left, gamma, 1.0, beta);
    } else if (x < right) {
        return smooth_step(x, gamma, right, beta, alpha);
    } else {
        return alpha;
    }
}

/**
 * First derivative of F(x).
 */
inline double Fp(double x, double alpha, double beta, double gamma, double delta) {
    double left = gamma - delta;
    double right = gamma + delta;

    if (x <= left || x >= right) {
        return 0.0;               // constant regions
    } else if (x < gamma) {
        return smooth_step_p(x, left, gamma, 1.0, beta);
    } else {  // x between gamma and gamma+delta (including gamma)
        return smooth_step_p(x, gamma, right, beta, alpha);
    }
}

/**
 * Second derivative of F(x).
 */
inline double Fpp(double x, double alpha, double beta, double gamma, double delta) {
    double left = gamma - delta;
    double right = gamma + delta;

    if (x <= left || x >= right) {
        return 0.0;               // constant regions
    } else if (x < gamma) {
        return smooth_step_pp(x, left, gamma, 1.0, beta);
    } else {  // x between gamma and gamma+delta (including gamma)
        return smooth_step_pp(x, gamma, right, beta, alpha);
    }
}