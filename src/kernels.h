#pragma once
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef USE_GPU
#pragma omp declare target
#endif

namespace Kernels {

// --------------------------------------------------------------------------------
// Cubic Spline Kernel (Monaghan 1992)
// Support radius is 1h. Returns W (value) and dWdh
// --------------------------------------------------------------------------------
inline void cubic_spline(double r, double h, double& W, double& dWdh) {
    double q = r / h;
    double h3 = h * h * h;
    double norm = 8.0 / (M_PI * h3);  // 3D normalization for 1h support radius
    double dWdr = 0.0;

    if (q < 0.5) {
        W = norm * (1.0 - 6.0 * q * q + 6.0 * q * q * q);
        if (r > 1e-12) {
            dWdr = norm * (-12.0 * q + 18.0 * q * q) / h;
        } else {
            dWdr = 0.0;
        }
    } else if (q < 1.0) {
        double diff = 1.0 - q;
        W = norm * 2.0 * diff * diff * diff;
        dWdr = norm * -6.0 * diff * diff / h;
    } else {
        W = 0.0;
        dWdr = 0.0;
    }

    dWdh = -(3.0 / h) * W - q * dWdr;
}

// --------------------------------------------------------------------------------
// Simplified Cubic Spline (Value Only)
// --------------------------------------------------------------------------------
inline void cubic_spline_value(double r, double h, double& W) {
    double q = r / h;
    double norm = 8.0 / (M_PI * h * h * h);
    if (q < 0.5) {
        W = norm * (1.0 - 6.0 * q * q + 6.0 * q * q * q);
    } else if (q < 1.0) {
        double diff = 1.0 - q;
        W = norm * 2.0 * diff * diff * diff;
    } else {
        W = 0.0;
    }
}

// --------------------------------------------------------------------------------
// Adaptive Gravitational Softening Kernel (GIZMO / Price & Monaghan 2007)
// Returns the modified 1/r^2 force factor and the d(phi)/dh term.
// --------------------------------------------------------------------------------
inline void gravity_derivatives(double r, double h, double& dphi_dr,
                                double& dphi_dh) {
    double q = r / h;
    double q2 = q * q;
    double q3 = q2 * q;
    double q4 = q2 * q2;
    double q5 = q3 * q2;

    double h2 = h * h;

    if (q < 0.5) {
        dphi_dr = (1.0 / h2) * (10.6666666667 * q - 38.4 * q3 + 32.0 * q4);
        dphi_dh = (1.0 / h2) * (2.8 - 16.0 * q2 + 48.0 * q4 - 38.4 * q5);
    } else if (q < 1.0) {
        dphi_dr = (1.0 / h2) * (-(1.0 / 15.0) / q2 + (64.0 / 3.0) * q -
                                48.0 * q2 + 38.4 * q3 - (32.0 / 3.0) * q4);
        dphi_dh =
            (1.0 / h2) * (3.2 - 32.0 * q2 + 64.0 * q3 - 48.0 * q4 + 12.8 * q5);
    } else {
        dphi_dr = 1.0 / (r * r);
        dphi_dh = 0.0;
    }
    // OLD BUGGY POLYNOMIAL
    /*if (q < 0.5) {
        dphi_dr =
            (1.0 / h2) * (10.6666666667 * q - 19.2 * q3 + 10.6666666667 * q4);
        dphi_dh = (1.0 / h2) * (2.8 - (5.3333333333) * q2 + (11.52) * q4 -
                                (7.1111111111) * q5);
    } else if (q < 1.0) {
        dphi_dr =
            (1.0 / h2) * (10.6666666667 * q - 19.2 * q3 + 10.6666666667 * q4 -
                          0.0666666667 / q2 + 0.1 * q5 - 3.2 * q4 +
                          10.6666666667 * q3 - 16.0 * q2 + 11.2 * q - 3.2);
        dphi_dh = (1.0 / h2) * (3.2 - (10.6666666667) * q + (19.2) * q2 -
                                (11.52) * q3 + (1.0666666667) * q4 - 0.08 / q);
    } else {
        dphi_dr = 1.0 / (r * r);
        dphi_dh = 0.0;
    }*/
}

// --------------------------------------------------------------------------------
// Adaptive Gravitational Softening Kernel (Forces & Gradients)
// Returns the d(phi)/dr potential derivative and d(W)/dr smoothing derivative
// --------------------------------------------------------------------------------
inline void adaptive_gravity_terms(double r, double h, double& dphi_dr,
                                   double& dW_dr) {
    if (r >= h) {
        dphi_dr = 1.0 / (r * r);
        dW_dr = 0.0;
        return;
    }

    double q = r / h;
    double q2 = q * q;
    double q3 = q2 * q;
    double q4 = q3 * q;
    double q5 = q4 * q;  // the buggy polynomial needs it

    double h2 = h * h;
    double h3 = h2 * h;

    if (q < 0.5) {
        dphi_dr = (1.0 / h2) * (10.6666666667 * q - 38.4 * q3 + 32.0 * q4);
    } else {
        dphi_dr = (1.0 / h2) * (-(1.0 / 15.0) / q2 + (64.0 / 3.0) * q -
                                48.0 * q2 + 38.4 * q3 - (32.0 / 3.0) * q4);
    }
    // OLD BUGGY POLYNOMIAL
    /*if (q < 0.5) {
        dphi_dr =
            (1.0 / h2) * (10.6666666667 * q - 19.2 * q3 + 10.6666666667 * q4);
    } else {
        dphi_dr =
            (1.0 / h2) * (10.6666666667 * q - 19.2 * q3 + 10.6666666667 * q4 -
                          0.0666666667 / q2 + 0.1 * q5 - 3.2 * q4 +
                          10.6666666667 * q3 - 16.0 * q2 + 11.2 * q - 3.2);
    }*/

    double norm = 8.0 / (M_PI * h3);
    if (q < 0.5) {
        dW_dr = norm * (-12.0 * q + 18.0 * q2) / h;
    } else {
        double diff = 1.0 - q;
        dW_dr = norm * (-6.0 * diff * diff) / h;
    }
}

}  // namespace Kernels

#ifdef USE_GPU
#pragma omp end declare target
#endif