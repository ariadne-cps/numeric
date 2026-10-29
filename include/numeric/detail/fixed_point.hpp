/***************************************************************************
 *            numeric/detail/fixed_point.hpp
 *
 *  Shared exact rounding primitives for native fixed-point arithmetic.
 ****************************************************************************/
#ifndef ARIADNE_NUMERIC_DETAIL_FIXED_POINT_HPP
#define ARIADNE_NUMERIC_DETAIL_FIXED_POINT_HPP

#include <cstddef>
#include "numeric/big_unsigned.hpp"

namespace Ariadne::Detail {

// dropped>0 and its guard bit is inside the canonical source view.
// odd is the parity of the retained integer after truncation.
inline bool nearest_even_binary_increment(BigUnsigned::LimbView source,
                                           std::size_t dropped,bool odd) noexcept {
    const auto guard=dropped-1u;
    return ((source[guard/64u]>>(guard%64u))&1u)!=0u &&
        (odd || BigUnsigned::any_bits_below(source,guard));
}

// The denominator is nonzero; callers validate scale metadata before
// entering this allocation boundary. Both exact and inexact ratios are valid.
inline BigUnsigned scaled_unsigned_quotient(BigUnsigned const& numerator,
                                            BigUnsigned const& denominator,
                                            std::size_t shift,bool round_upward) {
    auto [quotient,remainder]=divrem(shift_left(numerator,shift),denominator);
    if(round_upward && !remainder.is_zero()) quotient=add(quotient,BigUnsigned(1u));
    return quotient;
}

} // namespace Ariadne::Detail
#endif
