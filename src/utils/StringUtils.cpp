#include "utils/StringUtils.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>

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

}  // namespace utils
