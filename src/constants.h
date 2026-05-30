#pragma once

namespace constants {
    // ---------------------------------------------------------
    // Fundamental Physical Constants
    // ---------------------------------------------------------

    // Masses
    constexpr double M_P_SI  = 1.67262192e-27; // Proton mass [kg]
    constexpr double M_P_CGS = 1.67262192e-24; // Proton mass [g]
    
    // Thermodynamics
    constexpr double K_B_SI  = 1.380649e-23;   // Boltzmann constant [J/K]
    constexpr double K_B_CGS = 1.380649e-16;   // Boltzmann constant [erg/K]

    // ---------------------------------------------------------
    // Astronomical Conversions
    // ---------------------------------------------------------
    
    constexpr double MSUN_CGS = 1.98847e33;    // Solar mass [g]
    constexpr double MPC_CGS  = 3.08567758e24; // Megaparsec [cm]
    constexpr double GYR_CGS  = 3.15576e16;    // Gigayear [s]
}