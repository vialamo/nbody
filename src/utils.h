#pragma once
#include <string>

struct Config;

namespace Utils {
std::string get_timestamp();
std::string generate_directory_name(const Config& config,
                                    const std::string& timestamp);
}  // namespace Utils
