#include "cosmology.h"

#include "config.h"
#include "state.h"

namespace Cosmology {

void update(SimState& state, const Config& config) {
    if (!config.expanding_universe) {
        state.scale_factor = 1.0;
        state.hubble_param = 0.0;
        return;
    }

    if (config.omega_lambda == 0.0) {
        // Einstein-de Sitter
        double t_start = std::pow(config.a_start, 1.5);
        double current_t = t_start + state.total_time;
        state.scale_factor = std::pow(current_t, 2.0 / 3.0);
        state.hubble_param = (2.0 / 3.0) / current_t;
    } else {
        // Lambda-CDM
        const double Om = config.omega_m;
        const double Ol = config.omega_lambda;
        const double H0 = 2.0 / (3.0 * std::sqrt(Om));
        double factor = std::sqrt(Ol / Om);

        // Reverse-engineer the initial time from start_a
        double t_start = (2.0 / (3.0 * H0 * std::sqrt(Ol))) *
                         std::asinh(factor * std::pow(config.a_start, 1.5));
        double current_t = t_start + state.total_time;

        // Calculate current scale factor
        double sinh_term = std::sinh(1.5 * H0 * std::sqrt(Ol) * current_t);
        state.scale_factor =
            std::pow(Om / Ol, 1.0 / 3.0) * std::pow(sinh_term, 2.0 / 3.0);

        // Calculate current Hubble parameter H(a)
        double a3 = std::pow(state.scale_factor, 3.0);
        state.hubble_param = H0 * std::sqrt(Om / a3 + Ol);
    }
}

double get_time_from_scale_factor(double a, const Config& config) {
    if (!config.expanding_universe) {
        // Without expansion, 'a' doesn't change
        return std::numeric_limits<double>::infinity();
    }

    if (config.omega_lambda == 0.0) {
        // Einstein-de Sitter
        return std::pow(a, 1.5);
    } else {
        // Lambda-CDM
        const double Om = config.omega_m;
        const double Ol = config.omega_lambda;
        const double H0 = 2.0 / (3.0 * std::sqrt(Om));
        double factor = std::sqrt(Ol / Om);

        // This is the inverse of the Lambda-CDM solver
        return (2.0 / (3.0 * H0 * std::sqrt(Ol))) *
               std::asinh(factor * std::pow(a, 1.5));
    }
}

}  // namespace Cosmology
