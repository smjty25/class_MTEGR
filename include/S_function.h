/**
 * Smooth step function (C2 continuous) mapping x from [x0,x1] to [f0,f1].
 * For x outside [x0,x1], the caller should ensure it's used only within the interval.
 * Formula: f0 + (f1-f0) * (6y^5 - 15y^4 + 10y^3), where y = (x-x0)/(x1-x0).
 */
static inline double smooth_step(double x, double x0, double x1, double f0, double f1) {
    double y = (x - x0) / (x1 - x0);
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y2 * y2;
    double y5 = y4 * y;
    double poly = 6.0 * y5 - 15.0 * y4 + 10.0 * y3;
    return f0 + (f1 - f0) * poly;
}

/**
 * First derivative of smooth_step with respect to x.
 */
static inline double smooth_step_p(double x, double x0, double x1, double f0, double f1) {
    double y = (x - x0) / (x1 - x0);
    double factor = (f1 - f0) / (x1 - x0);
    // Derivative = factor * 30 * y^2 * (1-y)^2
    double t = y * (1.0 - y);
    return factor * 30.0 * t * t;   // 30 * y^2 * (1-y)^2
}

/**
 * Second derivative of smooth_step with respect to x.
 */
static inline double smooth_step_pp(double x, double x0, double x1, double f0, double f1) {
    double y = (x - x0) / (x1 - x0);
    double factor = (f1 - f0) / ((x1 - x0) * (x1 - x0));
    // Second derivative = factor * 60 * y * (1-y) * (1-2y)
    return factor * 60.0 * y * (1.0 - y) * (1.0 - 2.0 * y);
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