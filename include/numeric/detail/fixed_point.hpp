/***************************************************************************
 *            numeric/detail/fixed_point.hpp
 *
 *  Copyright  2026
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @file numeric/detail/fixed_point.hpp
 *  @brief Declarations and utilities for fixed point in Ariadne Numeric.
 */

#ifndef ARIADNE_NUMERIC_DETAIL_FIXED_POINT_HPP
#define ARIADNE_NUMERIC_DETAIL_FIXED_POINT_HPP

#include <cstddef>

#include "numeric/big_unsigned.hpp"

namespace Ariadne::Detail {

// dropped>0 and its guard bit is inside the canonical source view.
// odd is the parity of the retained integer after truncation.
inline bool nearest_even_binary_increment(BigUnsigned::LimbView source,
                                          std::size_t dropped, bool odd) noexcept {
    const auto guard = dropped - 1u;
    return ((source[guard / 64u] >> (guard % 64u)) & 1u) != 0u &&
           (odd || BigUnsigned::any_bits_below(source, guard));
}

// The denominator is nonzero; callers validate scale metadata before
// entering this allocation boundary. Both exact and inexact ratios are valid.
inline BigUnsigned scaled_unsigned_quotient(const BigUnsigned& numerator,
                                            const BigUnsigned& denominator,
                                            std::size_t shift, bool round_upward) {
    auto [quotient, remainder] = divrem(shift_left(numerator, shift), denominator);
    if (round_upward && !remainder.is_zero())
        quotient = add(quotient, BigUnsigned(1u));
    return quotient;
}

} // namespace Ariadne::Detail
#endif
