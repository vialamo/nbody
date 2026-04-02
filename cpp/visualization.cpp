#include "visualization.h"
#include <iostream>


static sf::Color lerpColor( sf::Color c1, sf::Color c2, double t ) {
    t = std::max( 0.0, std::min( 1.0, t ) );
    return sf::Color(
        static_cast< uint8_t >( c1.r + ( c2.r - c1.r ) * t ),
        static_cast< uint8_t >( c1.g + ( c2.g - c1.g ) * t ),
        static_cast< uint8_t >( c1.b + ( c2.b - c1.b ) * t )
    );
}

static sf::Color getPlasmaColor( double value ) {
    if( value < 0.25 ) return lerpColor( sf::Color( 0, 0, 0 ), sf::Color( 84, 2, 163 ), value / 0.25 );
    if( value < 0.5 ) return lerpColor( sf::Color( 84, 2, 163 ), sf::Color( 158, 28, 133 ), ( value - 0.25 ) / 0.25 );
    if( value < 0.75 ) return lerpColor( sf::Color( 158, 28, 133 ), sf::Color( 218, 107, 49 ), ( value - 0.5 ) / 0.25 );
    return lerpColor( sf::Color( 218, 107, 49 ), sf::Color( 255, 237, 200 ), ( value - 0.75 ) / 0.25 );
}


SimulationApp::SimulationApp( SimulationEngine& eng )
    : engine( eng ),
    window( sf::VideoMode( eng.config.RENDER_SIZE, eng.config.RENDER_SIZE ), "N-Body + Hydro Simulation" )
{
}

void SimulationApp::render_field( const Grid& field ) {
    const Config& cfg = engine.config;

    double min_f = field.minCoeff();
    double max_f = field.maxCoeff();
    double range_f = max_f - min_f;
    if( range_f < 1e-9 ) range_f = 1.0;

    sf::VertexArray field_mesh( sf::PrimitiveType::TriangleStrip, cfg.MESH_SIZE * cfg.MESH_SIZE * 4 );

    for( int i = 0; i < cfg.MESH_SIZE; ++i ) {
        for( int j = 0; j < cfg.MESH_SIZE; ++j ) {
            int i_next = ( i + 1 ); // No wrap for drawing
            int j_next = ( j + 1 ); // No wrap for drawing

            double val1 = ( field( i, j ) - min_f ) / range_f;
            sf::Color c1 = getPlasmaColor( val1 );
            sf::Vector2f pos1( static_cast< float >( i * cfg.CELL_SIZE * cfg.RENDER_SCALE ), static_cast< float >( j * cfg.CELL_SIZE * cfg.RENDER_SCALE ) );

            double val2 = ( j_next < cfg.MESH_SIZE ) ? ( field( i, j_next ) - min_f ) / range_f : val1;
            sf::Color c2 = ( j_next < cfg.MESH_SIZE ) ? getPlasmaColor( val2 ) : c1;
            sf::Vector2f pos2( static_cast< float >( i * cfg.CELL_SIZE * cfg.RENDER_SCALE ), static_cast< float >( j_next * cfg.CELL_SIZE * cfg.RENDER_SCALE ) );

            double val3 = ( i_next < cfg.MESH_SIZE && j_next < cfg.MESH_SIZE ) ? ( field( i_next, j_next ) - min_f ) / range_f : val1;
            sf::Color c3 = ( i_next < cfg.MESH_SIZE && j_next < cfg.MESH_SIZE ) ? getPlasmaColor( val3 ) : c1;
            sf::Vector2f pos3( static_cast< float >( i_next * cfg.CELL_SIZE * cfg.RENDER_SCALE ), static_cast< float >( j_next * cfg.CELL_SIZE * cfg.RENDER_SCALE ) );

            double val4 = ( i_next < cfg.MESH_SIZE ) ? ( field( i_next, j ) - min_f ) / range_f : val1;
            sf::Color c4 = ( i_next < cfg.MESH_SIZE ) ? getPlasmaColor( val4 ) : c1;
            sf::Vector2f pos4( static_cast< float >( i_next * cfg.CELL_SIZE * cfg.RENDER_SCALE ), static_cast< float >( j * cfg.CELL_SIZE * cfg.RENDER_SCALE ) );

            size_t idx = ( i * cfg.MESH_SIZE + j ) * 4;
            field_mesh[idx + 0].position = pos1; field_mesh[idx + 0].color = c1;
            field_mesh[idx + 1].position = pos2; field_mesh[idx + 1].color = c2;
            field_mesh[idx + 2].position = pos4; field_mesh[idx + 2].color = c4;
            field_mesh[idx + 3].position = pos3; field_mesh[idx + 3].color = c3;
        }
    }
    window.draw( field_mesh );
}

void SimulationApp::render() {
    const Config& cfg = engine.config;

    window.clear( sf::Color::Black );

    if( cfg.USE_HYDRO ) {
        render_field( engine.state.gas.pressure );
    }

    sf::CircleShape particle_shape( 1.0f );
    particle_shape.setFillColor( sf::Color::White );

    for( const auto& p : engine.state.dm.particles ) {
        particle_shape.setPosition( {
            static_cast< float >( p.pos.x * cfg.RENDER_SCALE ),
            static_cast< float >( p.pos.y * cfg.RENDER_SCALE )
            } );
        window.draw( particle_shape );
    }
    window.display();
}

void SimulationApp::save_frame( int frame_num ) {
    char filename[100];
    snprintf( filename, sizeof( filename ), "frame_%03d.png", frame_num );
    sf::Texture texture;
    texture.create( window.getSize().x, window.getSize().y );
    texture.update( window );
    if( texture.copyToImage().saveToFile( filename ) ) {
        std::cout << "Saved " << filename << std::endl;
    }
}

void SimulationApp::run() {
    const Config& cfg = engine.config;
    std::cout << "Starting graphical simulation loop..." << std::endl;

    while( window.isOpen() ) {
        sf::Event event;
        while( window.pollEvent( event ) ) {
            if( event.type == sf::Event::Closed ) {
                window.close();
            }
        }

        // Advance the physics engine
        engine.step();

        render();

        if( cfg.SAVE_RENDER_EVERY_CYCLES > 0 && engine.cycle_count % cfg.SAVE_RENDER_EVERY_CYCLES == 0 ) {
            save_frame( engine.cycle_count / cfg.SAVE_RENDER_EVERY_CYCLES );
        }
    }
}