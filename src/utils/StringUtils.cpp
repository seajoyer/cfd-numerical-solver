#include "utils/StringUtils.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <algorithm>

namespace utils {
auto DoubleWithoutDot(double value) -> std::string {
    std::ostringstream oss;

    // Handle zero as special case
    if (value == 0.0) {
        return "0e0";
    }

    // Handle negative numbers
    if (value < 0) {
        oss << '-';
        value = -value;
    }

    // Calculate base-10 exponent
    int exp = static_cast<int>(std::floor(std::log10(value)));

    // Maximum significant digits for double precision
    const int max_digits = std::numeric_limits<double>::digits10;

    // Scale value to get all significant digits as an integer
    int shift = max_digits - exp - 1;
    double scaled = value * std::pow(10.0, shift);

    // Round to handle floating-point imprecision
    auto mantissa = static_cast<int64_t>(std::round(scaled));

    // Remove trailing zeros, adjusting exponent accordingly
    while (mantissa != 0 && mantissa % 10 == 0) {
        mantissa /= 10;
        shift--;
    }

    // Final exponent is negation of shift
    int final_exp = -shift;

    oss << mantissa << 'e' << final_exp;
    return oss.str();
}

auto ToLower(std::string& str) -> std::string {
    std::string lower(str.size(), '\0');
    std::transform(str.begin(), str.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return lower;
}

auto GetTimestamp() -> std::string {
    using namespace std::chrono;
    
    // Get current time
    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
    
    // Convert to local time
    std::tm local_tm = *std::localtime(&now_time_t);
    
    // Format: DD-MM-YYYY_HH:MM:SS:mmm
    std::ostringstream oss;
    oss << std::setfill('0')
        << std::setw(2) << local_tm.tm_mday << '-'
        << std::setw(2) << (local_tm.tm_mon + 1) << '-'
        << (local_tm.tm_year + 1900) << '_'
        << std::setw(2) << local_tm.tm_hour << ':'
        << std::setw(2) << local_tm.tm_min << ':'
        << std::setw(2) << local_tm.tm_sec << ':'
        << std::setw(3) << now_ms.count();
    
    return oss.str();
}
}  // namespace utils
