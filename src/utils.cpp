#include "utils.h"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

std::string utils::get_timestamp() {
    std::time_t now = std::time(nullptr);
    std::tm ltm;
#ifdef _MSC_VER
    localtime_s(&ltm, &now);
#else
    ltm = *std::localtime(&now);
#endif
    std::stringstream ss;
    ss << std::put_time(&ltm, "%Y-%m-%d_%H-%M-%S");
    return ss.str();
}
