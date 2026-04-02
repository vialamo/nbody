#pragma once
#include <SFML/Graphics.hpp>
#include "engine.h"

class SimulationApp {
    SimulationEngine& engine;
    sf::RenderWindow window;

public:
    SimulationApp( SimulationEngine& eng );
    void render_field( const Grid& field );
    void render();
    void save_frame( int frame_num );
    void run();
};
