#ifndef STRINGUTILS_HPP
#define STRINGUTILS_HPP

#include <string>

namespace utils {

/**
 * @brief Converts a double to scientific notation string without decimal point.
 *
 * Formats a floating-point number as an integer mantissa followed by a
 * power-of-ten exponent, with no decimal point in the output.
 * Useful for creating filesystem-safe parameter strings.
 *
 * @param value The double value to convert.
 * @return std::string Formatted string (e.g., "123e-5" for 0.00123).
 *
 * @note Special cases:
 *       - Zero returns "0e0"
 *       - Negative values include leading minus sign
 *       - Trailing zeros in mantissa are removed
 *
 * @example
 * @code
 *     std::string s1 = utils::DoubleWithoutDot(0.00123);  // "123e-5"
 *     std::string s2 = utils::DoubleWithoutDot(1.5);      // "15e-1"
 *     std::string s3 = utils::DoubleWithoutDot(0.5);      // "5e-1"
 * @endcode
 */
auto DoubleWithoutDot(double value) -> std::string;

auto ToLower(std::string& str) -> std::string;

}  // namespace utils

#endif  // STRINGUTILS_HPP
