#include "gas.h"

HydroScratchPad::HydroScratchPad( int mesh_size ) :
    density( mesh_size ), mom_n( mesh_size ), mom_t1( mesh_size ), mom_t2( mesh_size ),
    energy( mesh_size ), v_n( mesh_size ), v_t1( mesh_size ), v_t2( mesh_size ),
    pressure( mesh_size ), rho_L( mesh_size ), p_L( mesh_size ),
    vn_L( mesh_size ), vt1_L( mesh_size ), vt2_L( mesh_size ), E_L( mesh_size ),
    mom_n_L( mesh_size ), mom_t1_L( mesh_size ), mom_t2_L( mesh_size ), rho_R( mesh_size ),
    p_R( mesh_size ), vn_R( mesh_size ), vt1_R( mesh_size ), vt2_R( mesh_size ),
    E_R( mesh_size ), mom_n_R( mesh_size ), mom_t1_R( mesh_size ), mom_t2_R( mesh_size ),
    cs_L( mesh_size ), cs_R( mesh_size ), S_L( mesh_size ), S_R( mesh_size ),
    S_R_minus_S_L( mesh_size ), F_dens_L( mesh_size ), F_dens_R( mesh_size ),
    F_momn_L( mesh_size ), F_momn_R( mesh_size ), F_momt1_L( mesh_size ),
    F_momt1_R( mesh_size ), F_momt2_L( mesh_size ), F_momt2_R( mesh_size ),
    F_en_L( mesh_size ), F_en_R( mesh_size ),
    flux_density( mesh_size ), flux_mom_n( mesh_size ),
    flux_mom_t1( mesh_size ), flux_mom_t2( mesh_size ), flux_energy( mesh_size ) {
}

Grid3D roll( const Grid3D& m, int shift, int axis ) {
    int N = m.n;
    Grid3D res( N );
    for( int i = 0; i < N; ++i ) {
        for( int j = 0; j < N; ++j ) {
            for( int k = 0; k < N; ++k ) {
                if( axis == 0 ) res( ( i + shift + N ) % N, j, k ) = m( i, j, k );
                else if( axis == 1 ) res( i, ( j + shift + N ) % N, k ) = m( i, j, k );
                else if( axis == 2 ) res( i, j, ( k + shift + N ) % N ) = m( i, j, k );
            }
        }
    }
    return res;
}

GasGrid::GasGrid( const Config& conf ) :
    density( conf.MESH_SIZE ), momentum_x( conf.MESH_SIZE ), momentum_y( conf.MESH_SIZE ), momentum_z( conf.MESH_SIZE ),
    energy( conf.MESH_SIZE ), pressure( conf.MESH_SIZE ), velocity_x( conf.MESH_SIZE ), velocity_y( conf.MESH_SIZE ), velocity_z( conf.MESH_SIZE ),
    accel_x( conf.MESH_SIZE ), accel_y( conf.MESH_SIZE ), accel_z( conf.MESH_SIZE ),
    scratch( conf.MESH_SIZE ), config( conf )
{
}

void GasGrid::update_primitive_variables() {
    for( int i = 0; i < config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE; ++i ) {
        if( density.data[i] > 1e-12 ) {
            velocity_x.data[i] = momentum_x.data[i] / density.data[i];
            velocity_y.data[i] = momentum_y.data[i] / density.data[i];
            velocity_z.data[i] = momentum_z.data[i] / density.data[i];
        }
        else {
            velocity_x.data[i] = velocity_y.data[i] = velocity_z.data[i] = 0.0;
        }
    }

    Grid3D kin_energy( config.MESH_SIZE );
    kin_energy.data = 0.5 * ( momentum_x.array().square() + momentum_y.array().square() + momentum_z.array().square() );

    for( int i = 0; i < config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE; ++i ) {
        if( density.data[i] > 1e-12 ) kin_energy.data[i] /= density.data[i];
        else kin_energy.data[i] = 0.0;
    }

    pressure.data = ( config.GAMMA - 1.0 ) * ( energy.data - kin_energy.data );
    pressure.data = ( pressure.array() < 1e-12 ).select( 1e-12, pressure.data );
}

void GasGrid::calculate_fluxes( int axis ) {
    scratch.density = density;
    scratch.energy = energy;
    scratch.pressure = pressure;

    if( axis == 0 ) {
        scratch.mom_n = momentum_x; scratch.mom_t1 = momentum_y; scratch.mom_t2 = momentum_z;
        scratch.v_n = velocity_x; scratch.v_t1 = velocity_y; scratch.v_t2 = velocity_z;
    }
    else if( axis == 1 ) {
        scratch.mom_n = momentum_y; scratch.mom_t1 = momentum_z; scratch.mom_t2 = momentum_x;
        scratch.v_n = velocity_y; scratch.v_t1 = velocity_z; scratch.v_t2 = velocity_x;
    }
    else {
        scratch.mom_n = momentum_z; scratch.mom_t1 = momentum_x; scratch.mom_t2 = momentum_y;
        scratch.v_n = velocity_z; scratch.v_t1 = velocity_x; scratch.v_t2 = velocity_y;
    }

    const int rollDir = -1;
    scratch.rho_L = scratch.density;
    scratch.p_L = scratch.pressure;
    scratch.vn_L = scratch.v_n;
    scratch.vt1_L = scratch.v_t1;
    scratch.vt2_L = scratch.v_t2;
    scratch.E_L = scratch.energy;
    scratch.mom_n_L = scratch.mom_n;
    scratch.mom_t1_L = scratch.mom_t1;
    scratch.mom_t2_L = scratch.mom_t2;

    scratch.rho_R = roll( scratch.rho_L, rollDir, axis );
    scratch.p_R = roll( scratch.p_L, rollDir, axis );
    scratch.vn_R = roll( scratch.vn_L, rollDir, axis );
    scratch.vt1_R = roll( scratch.vt1_L, rollDir, axis );
    scratch.vt2_R = roll( scratch.vt2_L, rollDir, axis );
    scratch.E_R = roll( scratch.E_L, rollDir, axis );
    scratch.mom_n_R = roll( scratch.mom_n_L, rollDir, axis );
    scratch.mom_t1_R = roll( scratch.mom_t1_L, rollDir, axis );
    scratch.mom_t2_R = roll( scratch.mom_t2_L, rollDir, axis );

    scratch.cs_L.data = ( config.GAMMA * scratch.p_L.array() / scratch.rho_L.array() ).sqrt();
    scratch.cs_R.data = ( config.GAMMA * scratch.p_R.array() / scratch.rho_R.array() ).sqrt();
    scratch.cs_L.data = ( scratch.rho_L.array() > 1e-12 ).select( scratch.cs_L.array(), 0.0 );
    scratch.cs_R.data = ( scratch.rho_R.array() > 1e-12 ).select( scratch.cs_R.array(), 0.0 );
    scratch.cs_L.data = ( scratch.cs_L.array() == scratch.cs_L.array() ).select( scratch.cs_L.array(), 0.0 );
    scratch.cs_R.data = ( scratch.cs_R.array() == scratch.cs_R.array() ).select( scratch.cs_R.array(), 0.0 );

    scratch.S_L.data = ( scratch.vn_L.array() - scratch.cs_L.array() ).cwiseMin( scratch.vn_R.array() - scratch.cs_R.array() );
    scratch.S_R.data = ( scratch.vn_L.array() + scratch.cs_L.array() ).cwiseMax( scratch.vn_R.array() + scratch.cs_R.array() );

    scratch.F_dens_L.data = scratch.rho_L.array() * scratch.vn_L.array();
    scratch.F_dens_R.data = scratch.rho_R.array() * scratch.vn_R.array();
    scratch.F_momn_L.data = scratch.rho_L.array() * scratch.vn_L.array().square() + scratch.p_L.array();
    scratch.F_momn_R.data = scratch.rho_R.array() * scratch.vn_R.array().square() + scratch.p_R.array();
    scratch.F_momt1_L.data = scratch.rho_L.array() * scratch.vn_L.array() * scratch.vt1_L.array();
    scratch.F_momt1_R.data = scratch.rho_R.array() * scratch.vn_R.array() * scratch.vt1_R.array();
    scratch.F_momt2_L.data = scratch.rho_L.array() * scratch.vn_L.array() * scratch.vt2_L.array();
    scratch.F_momt2_R.data = scratch.rho_R.array() * scratch.vn_R.array() * scratch.vt2_R.array();
    scratch.F_en_L.data = ( scratch.E_L.array() + scratch.p_L.array() ) * scratch.vn_L.array();
    scratch.F_en_R.data = ( scratch.E_R.array() + scratch.p_R.array() ) * scratch.vn_R.array();

    scratch.S_R_minus_S_L.data = scratch.S_R.array() - scratch.S_L.array();
    scratch.S_R_minus_S_L.data = ( scratch.S_R_minus_S_L.array().abs() < 1e-9 ).select( 1e-9, scratch.S_R_minus_S_L.data );

    auto HLL = [&]( const Grid3D& FL, const Grid3D& FR, const Grid3D& UL, const Grid3D& UR ) {
        Grid3D flux( config.MESH_SIZE );
        flux.data = ( scratch.S_L.array() >= 0 ).select( FL.array(),
            ( scratch.S_R.array() <= 0 ).select( FR.array(),
                ( scratch.S_R.array() * FL.array() - scratch.S_L.array() * FR.array() +
                    scratch.S_L.array() * scratch.S_R.array() * ( UR.array() - UL.array() ) ) / scratch.S_R_minus_S_L.array() ) );
        return flux;
        };

    scratch.flux_density = HLL( scratch.F_dens_L, scratch.F_dens_R, scratch.rho_L, scratch.rho_R );
    scratch.flux_mom_n = HLL( scratch.F_momn_L, scratch.F_momn_R, scratch.mom_n_L, scratch.mom_n_R );
    scratch.flux_mom_t1 = HLL( scratch.F_momt1_L, scratch.F_momt1_R, scratch.mom_t1_L, scratch.mom_t1_R );
    scratch.flux_mom_t2 = HLL( scratch.F_momt2_L, scratch.F_momt2_R, scratch.mom_t2_L, scratch.mom_t2_R );
    scratch.flux_energy = HLL( scratch.F_en_L, scratch.F_en_R, scratch.E_L, scratch.E_R );
}

void GasGrid::hydro_step( double dt ) {
    update_primitive_variables();
    double factor = dt / config.CELL_SIZE;

    // X-Sweep
    calculate_fluxes( 0 );
    density.data = density.array() - factor * ( scratch.flux_density.array() - roll( scratch.flux_density, 1, 0 ).array() );
    momentum_x.data = momentum_x.array() - factor * ( scratch.flux_mom_n.array() - roll( scratch.flux_mom_n, 1, 0 ).array() );
    momentum_y.data = momentum_y.array() - factor * ( scratch.flux_mom_t1.array() - roll( scratch.flux_mom_t1, 1, 0 ).array() );
    momentum_z.data = momentum_z.array() - factor * ( scratch.flux_mom_t2.array() - roll( scratch.flux_mom_t2, 1, 0 ).array() );
    energy.data = energy.array() - factor * ( scratch.flux_energy.array() - roll( scratch.flux_energy, 1, 0 ).array() );
    update_primitive_variables();

    // Y-Sweep
    calculate_fluxes( 1 );
    density.data = density.array() - factor * ( scratch.flux_density.array() - roll( scratch.flux_density, 1, 1 ).array() );
    momentum_x.data = momentum_x.array() - factor * ( scratch.flux_mom_t2.array() - roll( scratch.flux_mom_t2, 1, 1 ).array() );
    momentum_y.data = momentum_y.array() - factor * ( scratch.flux_mom_n.array() - roll( scratch.flux_mom_n, 1, 1 ).array() );
    momentum_z.data = momentum_z.array() - factor * ( scratch.flux_mom_t1.array() - roll( scratch.flux_mom_t1, 1, 1 ).array() );
    energy.data = energy.array() - factor * ( scratch.flux_energy.array() - roll( scratch.flux_energy, 1, 1 ).array() );
    update_primitive_variables();

    // Z-Sweep
    calculate_fluxes( 2 );
    density.data = density.array() - factor * ( scratch.flux_density.array() - roll( scratch.flux_density, 1, 2 ).array() );
    momentum_x.data = momentum_x.array() - factor * ( scratch.flux_mom_t1.array() - roll( scratch.flux_mom_t1, 1, 2 ).array() );
    momentum_y.data = momentum_y.array() - factor * ( scratch.flux_mom_t2.array() - roll( scratch.flux_mom_t2, 1, 2 ).array() );
    momentum_z.data = momentum_z.array() - factor * ( scratch.flux_mom_n.array() - roll( scratch.flux_mom_n, 1, 2 ).array() );
    energy.data = energy.array() - factor * ( scratch.flux_energy.array() - roll( scratch.flux_energy, 1, 2 ).array() );
    update_primitive_variables();
}

double GasGrid::get_cfl_timestep() const {
    if( !config.USE_HYDRO ) return std::numeric_limits<double>::infinity();

    Grid3D cs_sq( config.MESH_SIZE );
    cs_sq.data = ( config.GAMMA * pressure.array() / density.array() );
    cs_sq.data = ( density.array() > 1e-12 ).select( cs_sq.array(), 0.0 );
    cs_sq.data = ( cs_sq.array() == cs_sq.array() ).select( cs_sq.array(), 0.0 );

    Grid3D v_mag( config.MESH_SIZE );
    v_mag.data = ( velocity_x.array().square() + velocity_y.array().square() + velocity_z.array().square() ).sqrt();

    Grid3D signal_vel( config.MESH_SIZE );
    signal_vel.data = v_mag.array() + cs_sq.array().sqrt();

    double max_signal_vel = signal_vel.maxCoeff();
    if( max_signal_vel < 1e-9 ) return std::numeric_limits<double>::infinity();

    return ( config.CELL_SIZE / max_signal_vel ) * config.CFL_SAFETY_FACTOR;
}