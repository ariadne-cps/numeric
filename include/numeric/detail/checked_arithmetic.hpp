/***************************************************************************
 *            numeric/detail/checked_arithmetic.hpp
 *
 *  Copyright  2026
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @file numeric/detail/checked_arithmetic.hpp
 *  @brief Declarations and utilities for checked arithmetic in Ariadne Numeric.
 */

#ifndef ARIADNE_NUMERIC_DETAIL_CHECKED_ARITHMETIC_HPP
#define ARIADNE_NUMERIC_DETAIL_CHECKED_ARITHMETIC_HPP

#include <concepts>
#include <limits>
#include <stdexcept>

namespace Ariadne::Detail {

template<std::unsigned_integral To, std::unsigned_integral From>
inline constexpr bool unsigned_narrowing =
    std::numeric_limits<To>::digits < std::numeric_limits<From>::digits;

template<std::unsigned_integral To, std::unsigned_integral From>
constexpr bool unsigned_fits(From value) noexcept {
    if constexpr (unsigned_narrowing<To, From>) {
        return value <= static_cast<From>(std::numeric_limits<To>::max());
    } else {
        (void)value;
        return true;
    }
}

// Compare the mathematical sum left+right with target, even when the
// sum cannot be represented by T. Subtraction is performed only in range.
template<std::unsigned_integral T>
constexpr int compare_unsigned_sum(T left, T right, T target) noexcept {
    if (left > target)
        return 1;
    const T remaining = target - left;
    return right < remaining ? -1 : right > remaining ? 1
                                                      : 0;
}

template<std::unsigned_integral T>
constexpr T checked_unsigned_multiply(T left, T right, const char* message) {
    if (right != 0 && left > std::numeric_limits<T>::max() / right)
        throw std::overflow_error(message);
    return static_cast<T>(left * right);
}

template<std::unsigned_integral T>
constexpr T checked_unsigned_add(T left, T right, const char* message) {
    if (left > std::numeric_limits<T>::max() - right)
        throw std::overflow_error(message);
    return left + right;
}

template<std::unsigned_integral To, std::unsigned_integral From>
constexpr To checked_unsigned_cast(From value, const char* message) {
    if constexpr (unsigned_narrowing<To, From>) {
        if (!unsigned_fits<To>(value))
            throw std::overflow_error(message);
    }
    return static_cast<To>(value);
}

} // namespace Ariadne::Detail
#endif
