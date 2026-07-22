#pragma once
#include <string>
#include <vector>

#include "config.h"

class Cooling {
   private:
    struct CoolingData {
        std::vector<double> redshift;
        std::vector<double> density;
        std::vector<double> temperature;
        std::vector<double> cooling_metal;
        std::vector<double> heating_metal;
        std::vector<double> cooling_primordial;
        std::vector<double> heating_primordial;
        size_t Nz, Nd, Nt;
    };

    // Encapsulated state
    CoolingData cooling_data;
    bool use_table = false;
    double mean_rho_code = 0.0;
    double subgrid_clumping_A = 0.0;

    // Private helper methods
    double compute_clumping_factor(double rho_code, double a,
                                   const Config& config) const;

    void load_cooling_tables(const Config& config);

    int get_index(double value, const std::vector<double>& axis) const;

    double trilinear_interpolate(const std::vector<double>& data, double z,
                                 double log_nH, double log_T) const;

    double compute_tabulated_cooling(double u_code, double rho_code, double Z,
                                     double a, double clumping_factor,
                                     const Config& config) const;

   public:
    static constexpr int MAX_ITER = 50;

    // Constructor replaces the old initialize() function
    Cooling(const Config& config);

    // Inline conversion helpers
    static inline double get_temp_from_internal_energy(double u_code, double a,
                                                       const Config& config) {
        return u_code * (a * a) * config.factor_u_to_t;
    }

    static inline double get_internal_energy_from_temp(double T_kelvin,
                                                       double a,
                                                       const Config& config) {
        return (T_kelvin * config.factor_t_to_u) / (a * a);
    }

    // Main public interfaces
    double compute_du_dt(double u_code, double rho_code, double Z, double a,
                         const Config& config) const;

    double solve_cooling_implicit(double u_old, double rho_code, double Z,
                                  double a, double dt, double u_rad_floor,
                                  const Config& config,
                                  int& iterations_taken) const;

    double get_u_rad_floor(double a, const Config& config) const;
};
