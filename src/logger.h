#pragma once
#include <chrono>
#include <fstream>
#include <string>

// Forward declaration
class Diagnostics;

class Logger {
   private:
    std::ofstream log_file;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    bool header_written = false;

    static std::string format_double(double val, int precision,
                                     bool scientific = false);

   public:
    Logger(const std::string& run_dir);
    ~Logger();

    void write_header();
    void log(const Diagnostics& diag);
};