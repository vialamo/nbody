/*
 * A simple, single-header C++ INI parser.
 * Based on ini.h by rxi (https://github.com/rxi/ini)
 * Adapted for C++ std::map and std::string.
 */

#ifndef INI_H
#define INI_H

#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

namespace ini {

class Ini {
public:
    using Section = std::map<std::string, std::string>;
    std::map<std::string, Section> sections;

    // Load from file
    void load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("INI: Could not open file: " + filename);
        }
        std::string line, section;
        while (std::getline(file, line)) {
            line = trim(line);
            if (line.empty() || line[0] == ';') continue;
            if (line[0] == '[') {
                section = trim(line.substr(1, line.find(']') - 1));
            } else {
                size_t eq_pos = line.find('=');
                if (eq_pos != std::string::npos) {
                    std::string key = trim(line.substr(0, eq_pos));
                    std::string val = trim(line.substr(eq_pos + 1));
                    sections[section][key] = val;
                }
            }
        }
    }

    // Get string value
    std::string get(const std::string& section, const std::string& key, const std::string& default_val = "") {
        if (sections.find(section) == sections.end() || sections[section].find(key) == sections[section].end()) {
            return default_val;
        }
        return sections[section][key];
    }

    // Get integer value
    int get_int(const std::string& section, const std::string& key, int default_val = 0) {
        std::string val = get(section, key);
        if (val.empty()) return default_val;
        try {
            return std::stoi(val);
        } catch (...) {
            return default_val;
        }
    }

    // Get double value
    double get_double(const std::string& section, const std::string& key, double default_val = 0.0) {
        std::string val = get(section, key);
        if (val.empty()) return default_val;
        try {
            return std::stod(val);
        } catch (...) {
            return default_val;
        }
    }

    // Get boolean value
    bool get_bool(const std::string& section, const std::string& key, bool default_val = false) {
        std::string val = get(section, key);
        if (val.empty()) return default_val;
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        if (val == "true" || val == "yes" || val == "1") return true;
        if (val == "false" || val == "no" || val == "0") return false;
        return default_val;
    }

private:
    std::string trim(const std::string& str) {
        size_t first = str.find_first_not_of(" \t");
        if (std::string::npos == first) return "";
        size_t last = str.find_last_not_of(" \t");
        return str.substr(first, (last - first + 1));
    }
};

} // namespace ini

#endif
