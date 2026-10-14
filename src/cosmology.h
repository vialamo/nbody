struct Config;
struct SimState;

void update_cosmology(SimState& state, const Config& config);
double get_time_from_scale_factor(double a, const Config& config);
