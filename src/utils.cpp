#include "utils.h"

#include <ctime>
#include <filesystem>
#include <iomanip>
#include <sstream>

#include "config.h"

namespace Utils {
std::string get_timestamp() {
    std::time_t now = std::time(nullptr);
    std::tm ltm;
    localtime_r(&now, &ltm);
    std::stringstream ss;
    ss << std::put_time(&ltm, "%Y-%m-%d_%H-%M-%S");
    return ss.str();
}

std::string generate_directory_name(const Config& config,
                                    const std::string& timestamp) {
    std::string base = InitialConfig::to_string(config.initial_setup);
    std::filesystem::path base_path(config.outputs_folder_path);
    std::filesystem::path full_path = base_path / (base + "_" + timestamp);
    return full_path.string();
}
}  // namespace Utils
