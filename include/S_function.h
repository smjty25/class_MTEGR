#include <math.h>

// Helper: hyperbolic tangent and its derivatives (sech^2, -2 sech^2 tanh)
static inline double tanh_deriv(double u) {
    double sech2 = 1.0 / (cosh(u) * cosh(u));
    return sech2;
}
static inline double tanh_deriv2(double u) {
    double sech2 = 1.0 / (cosh(u) * cosh(u));
    return -2.0 * sech2 * tanh(u);
}

// Helper: Gaussian exp(-u^2) and its derivatives
static inline double gauss(double u) {
    return exp(-u * u);
}
static inline double gauss_deriv(double u) {
    return -2.0 * u * exp(-u * u);
}
static inline double gauss_deriv2(double u) {
    return (4.0 * u * u - 2.0) * exp(-u * u);
}

// F(x) = 1 + (alpha-1)*H(x) + (beta - (1+alpha)/2)*G(x)
// where H(x) = (1+tanh((x-gamma)/delta))/2, G(x) = exp(-((x-gamma)/delta)^2)
double F(double x, double alpha, double beta, double gamma, double delta) {
    double u = (x - gamma) / delta;
    double H = (1.0 + tanh(u)) / 2.0;
    double G = gauss(u);
    double offset = beta - (1.0 + alpha) / 2.0;
    return 1.0 + (alpha - 1.0) * H + offset * G;
}

// First derivative dF/dx
double Fp(double x, double alpha, double beta, double gamma, double delta) {
    double u = (x - gamma) / delta;
    double dH = tanh_deriv(u) / (2.0 * delta);          // dH/dx = (1/2)*sech^2(u)/delta
    double dG = gauss_deriv(u) / delta;                 // dG/dx = -2u*exp(-u^2)/delta
    double offset = beta - (1.0 + alpha) / 2.0;
    return (alpha - 1.0) * dH + offset * dG;
}

// Second derivative d²F/dx²
double Fpp(double x, double alpha, double beta, double gamma, double delta) {
    double u = (x - gamma) / delta;
    double d2H = tanh_deriv2(u) / (2.0 * delta * delta);   // d²H/dx² = (1/2)*d²(tanh)/du² / delta²
    double d2G = gauss_deriv2(u) / (delta * delta);        // d²G/dx² = (4u²-2)*exp(-u²) / delta²
    double offset = beta - (1.0 + alpha) / 2.0;
    return (alpha - 1.0) * d2H + offset * d2G;
}