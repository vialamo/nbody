#include "gas.h"

HydroScratchPad::HydroScratchPad( int mesh_size ) :
    density( mesh_size, mesh_size ), mom_n( mesh_size, mesh_size ), mom_t( mesh_size, mesh_size ),
    energy( mesh_size, mesh_size ), v_n( mesh_size, mesh_size ), v_t( mesh_size, mesh_size ),
    pressure( mesh_size, mesh_size ), rho_L( mesh_size, mesh_size ), p_L( mesh_size, mesh_size ),
    vn_L( mesh_size, mesh_size ), vt_L( mesh_size, mesh_size ), E_L( mesh_size, mesh_size ),
    mom_n_L( mesh_size, mesh_size ), mom_t_L( mesh_size, mesh_size ), rho_R( mesh_size, mesh_size ),
    p_R( mesh_size, mesh_size ), vn_R( mesh_size, mesh_size ), vt_R( mesh_size, mesh_size ),
    E_R( mesh_size, mesh_size ), mom_n_R( mesh_size, mesh_size ), cs_L( mesh_size, mesh_size ),
    cs_R( mesh_size, mesh_size ), S_L( mesh_size, mesh_size ), S_R( mesh_size, mesh_size ),
    S_R_minus_S_L( mesh_size, mesh_size ), F_dens_L( mesh_size, mesh_size ), F_dens_R( mesh_size, mesh_size ),
    F_momn_L( mesh_size, mesh_size ), F_momn_R( mesh_size, mesh_size ), F_momt_L( mesh_size, mesh_size ),
    F_momt_R( mesh_size, mesh_size ), F_en_L( mesh_size, mesh_size ), F_en_R( mesh_size, mesh_size ),
    flux_density( mesh_size, mesh_size ), flux_mom_n( mesh_size, mesh_size ),
    flux_mom_t( mesh_size, mesh_size ), flux_energy( mesh_size, mesh_size ) {
}

Grid roll( const Grid& m, int shift, int axis ) {
    Grid result( m.rows(), m.cols() );
    if( axis == 0 ) {
        int rows = static_cast< int >( m.rows() );
        shift = ( shift % rows + rows ) % rows;
        result.bottomRows( rows - shift ) = m.topRows( rows - shift );
        result.topRows( shift ) = m.bottomRows( shift );
    }
    else {
        int cols = static_cast< int >( m.cols() );
        shift = ( shift % cols + cols ) % cols;
        result.rightCols( cols - shift ) = m.leftCols( cols - shift );
        result.leftCols( shift ) = m.rightCols( shift );
    }
    return result;
}

GasGrid::GasGrid( const Config& conf ) :
    density( conf.MESH_SIZE, conf.MESH_SIZE ),
    momentum_x( conf.MESH_SIZE, conf.MESH_SIZE ),
    momentum_y( conf.MESH_SIZE, conf.MESH_SIZE ),
    energy( conf.MESH_SIZE, conf.MESH_SIZE ),
    pressure( conf.MESH_SIZE, conf.MESH_SIZE ),
    velocity_x( conf.MESH_SIZE, conf.MESH_SIZE ),
    velocity_y( conf.MESH_SIZE, conf.MESH_SIZE ),
    accel_x( conf.MESH_SIZE, conf.MESH_SIZE ),
    accel_y( conf.MESH_SIZE, conf.MESH_SIZE ),
    scratch( conf.MESH_SIZE ),
    config( conf )
{
    density.setZero(); momentum_x.setZero(); momentum_y.setZero(); energy.setZero();
    pressure.setZero(); velocity_x.setZero(); velocity_y.setZero();
    accel_x.setZero(); accel_y.setZero();
}

void GasGrid::update_primitive_variables() {
    for( int i = 0; i < config.MESH_SIZE; ++i ) {
        for( int j = 0; j < config.MESH_SIZE; ++j ) {
            if( density( i, j ) > 1e-12 ) {
                velocity_x( i, j ) = momentum_x( i, j ) / density( i, j );
                velocity_y( i, j ) = momentum_y( i, j ) / density( i, j );
            }
            else {
                velocity_x( i, j ) = 0.0;
                velocity_y( i, j ) = 0.0;
            }
        }
    }

    Grid kinetic_energy_density = 0.5 * ( momentum_x.array().square() + momentum_y.array().square() );
    for( int i = 0; i < config.MESH_SIZE; ++i ) {
        for( int j = 0; j < config.MESH_SIZE; ++j ) {
            if( density( i, j ) > 1e-12 ) {
                kinetic_energy_density( i, j ) /= density( i, j );
            }
            else {
                kinetic_energy_density( i, j ) = 0.0;
            }
        }
    }

    Grid internal_energy_density = energy - kinetic_energy_density;
    pressure = ( config.GAMMA - 1.0 ) * internal_energy_density;

    // Apply pressure floor
    pressure = ( pressure.array() < 1e-12 ).select( 1e-12, pressure );
}

// Helper: Calculates HLL fluxes for one axis
void GasGrid::calculate_fluxes( int axis ) {
    if( axis == 1 ) { // Permute axes for y-direction
        scratch.density = density.transpose();
        scratch.mom_n = momentum_y.transpose();
        scratch.mom_t = momentum_x.transpose();
        scratch.energy = energy.transpose();
        scratch.v_n = velocity_y.transpose();
        scratch.v_t = velocity_x.transpose();
        scratch.pressure = pressure.transpose();
    }
    else {
        scratch.density = density;
        scratch.mom_n = momentum_x;
        scratch.mom_t = momentum_y;
        scratch.energy = energy;
        scratch.v_n = velocity_x;
        scratch.v_t = velocity_y;
        scratch.pressure = pressure;
    }

    const int rollDir = -1;
    scratch.rho_L = scratch.density;
    scratch.p_L = scratch.pressure;
    scratch.vn_L = scratch.v_n;
    scratch.vt_L = scratch.v_t;
    scratch.E_L = scratch.energy;
    scratch.mom_n_L = scratch.mom_n;
    scratch.mom_t_L = scratch.mom_t;
    scratch.rho_R = roll( scratch.rho_L, rollDir, 0 );
    scratch.p_R = roll( scratch.p_L, rollDir, 0 );
    scratch.vn_R = roll( scratch.vn_L, rollDir, 0 );
    scratch.vt_R = roll( scratch.vt_L, rollDir, 0 );
    scratch.E_R = roll( scratch.E_L, rollDir, 0 );
    scratch.mom_n_R = roll( scratch.mom_n_L, rollDir, 0 );
    scratch.mom_t_R = roll( scratch.mom_t_L, rollDir, 0 );

    scratch.cs_L = ( config.GAMMA * scratch.p_L.array() / scratch.rho_L.array() ).sqrt();
    scratch.cs_R = ( config.GAMMA * scratch.p_R.array() / scratch.rho_R.array() ).sqrt();

    scratch.cs_L = ( scratch.rho_L.array() > 1e-12 ).select( scratch.cs_L, 0.0 );
    scratch.cs_R = ( scratch.rho_R.array() > 1e-12 ).select( scratch.cs_R, 0.0 );
    // Handle NaNs from sqrt(negative pressure)
    scratch.cs_L = ( scratch.cs_L.array() == scratch.cs_L.array() ).select( scratch.cs_L, 0.0 );
    scratch.cs_R = ( scratch.cs_R.array() == scratch.cs_R.array() ).select( scratch.cs_R, 0.0 );


    scratch.S_L = ( scratch.vn_L.array() - scratch.cs_L.array() ).cwiseMin( scratch.vn_R.array() - scratch.cs_R.array() );
    scratch.S_R = ( scratch.vn_L.array() + scratch.cs_L.array() ).cwiseMax( scratch.vn_R.array() + scratch.cs_R.array() );

    scratch.F_dens_L = scratch.rho_L.array() * scratch.vn_L.array();
    scratch.F_dens_R = scratch.rho_R.array() * scratch.vn_R.array();
    scratch.F_momn_L = scratch.rho_L.array() * scratch.vn_L.array().square() + scratch.p_L.array();
    scratch.F_momn_R = scratch.rho_R.array() * scratch.vn_R.array().square() + scratch.p_R.array();
    scratch.F_momt_L = scratch.rho_L.array() * scratch.vn_L.array() * scratch.vt_L.array();
    scratch.F_momt_R = scratch.rho_R.array() * scratch.vn_R.array() * scratch.vt_R.array();
    scratch.F_en_L = ( scratch.E_L.array() + scratch.p_L.array() ) * scratch.vn_L.array();
    scratch.F_en_R = ( scratch.E_R.array() + scratch.p_R.array() ) * scratch.vn_R.array();

    scratch.S_R_minus_S_L = scratch.S_R - scratch.S_L;
    scratch.S_R_minus_S_L = ( scratch.S_R_minus_S_L.array().abs() < 1e-9 ).select( 1e-9, scratch.S_R_minus_S_L );

    scratch.flux_density = ( scratch.S_L.array() >= 0 ).select( scratch.F_dens_L,
        ( scratch.S_R.array() <= 0 ).select( scratch.F_dens_R,
            ( scratch.S_R.array() * scratch.F_dens_L.array() - scratch.S_L.array() * scratch.F_dens_R.array() + scratch.S_L.array() * scratch.S_R.array() * ( scratch.rho_R.array() - scratch.rho_L.array() ) ) / scratch.S_R_minus_S_L.array() ) );
    scratch.flux_mom_n = ( scratch.S_L.array() >= 0 ).select( scratch.F_momn_L,
        ( scratch.S_R.array() <= 0 ).select( scratch.F_momn_R,
            ( scratch.S_R.array() * scratch.F_momn_L.array() - scratch.S_L.array() * scratch.F_momn_R.array() + scratch.S_L.array() * scratch.S_R.array() * ( scratch.mom_n_R.array() - scratch.mom_n_L.array() ) ) / scratch.S_R_minus_S_L.array() ) );
    scratch.flux_mom_t = ( scratch.S_L.array() >= 0 ).select( scratch.F_momt_L,
        ( scratch.S_R.array() <= 0 ).select( scratch.F_momt_R,
            ( scratch.S_R.array() * scratch.F_momt_L.array() - scratch.S_L.array() * scratch.F_momt_R.array() + scratch.S_L.array() * scratch.S_R.array() * ( scratch.mom_t_R.array() - scratch.mom_t_L.array() ) ) / scratch.S_R_minus_S_L.array() ) );
    scratch.flux_energy = ( scratch.S_L.array() >= 0 ).select( scratch.F_en_L,
        ( scratch.S_R.array() <= 0 ).select( scratch.F_en_R,
            ( scratch.S_R.array() * scratch.F_en_L.array() - scratch.S_L.array() * scratch.F_en_R.array() + scratch.S_L.array() * scratch.S_R.array() * ( scratch.E_R.array() - scratch.E_L.array() ) ) / scratch.S_R_minus_S_L.array() ) );

    if( axis == 1 ) {
        scratch.flux_density.transposeInPlace();
        scratch.flux_mom_t.transposeInPlace();
        scratch.flux_mom_n.transposeInPlace();
        scratch.flux_energy.transposeInPlace();
    }
}

void GasGrid::hydro_step( double dt ) {
    update_primitive_variables();

    calculate_fluxes( 0 );
    double factor = dt / config.CELL_SIZE;
    density -= factor * ( scratch.flux_density - roll( scratch.flux_density, 1, 0 ) );
    momentum_x -= factor * ( scratch.flux_mom_n - roll( scratch.flux_mom_n, 1, 0 ) );
    momentum_y -= factor * ( scratch.flux_mom_t - roll( scratch.flux_mom_t, 1, 0 ) );
    energy -= factor * ( scratch.flux_energy - roll( scratch.flux_energy, 1, 0 ) );
    update_primitive_variables();

    calculate_fluxes( 1 );
    density -= factor * ( scratch.flux_density - roll( scratch.flux_density, 1, 1 ) );
    momentum_x -= factor * ( scratch.flux_mom_t - roll( scratch.flux_mom_t, 1, 1 ) );
    momentum_y -= factor * ( scratch.flux_mom_n - roll( scratch.flux_mom_n, 1, 1 ) );
    energy -= factor * ( scratch.flux_energy - roll( scratch.flux_energy, 1, 1 ) );
    update_primitive_variables();
}

double GasGrid::get_cfl_timestep() const {
    if( !config.USE_HYDRO ) return std::numeric_limits<double>::infinity();

    Grid cs_sq = ( config.GAMMA * pressure.array() / density.array() );
    cs_sq = ( density.array() > 1e-12 ).select( cs_sq, 0.0 );
    cs_sq = ( cs_sq.array() == cs_sq.array() ).select( cs_sq, 0.0 );

    Grid cs = cs_sq.array().sqrt();
    Grid v_mag = ( velocity_x.array().square() + velocity_y.array().square() ).sqrt();
    Grid signal_vel = v_mag + cs;

    double max_signal_vel = signal_vel.maxCoeff();
    if( max_signal_vel < 1e-9 ) return std::numeric_limits<double>::infinity();

    return ( config.CELL_SIZE / max_signal_vel ) * config.CFL_SAFETY_FACTOR;
}