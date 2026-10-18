/***************************************************************************
 *            numeric/floatap_native_ratio.hpp
 *
 *  Copyright  2026  Ariadne contributors
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

/*! \file numeric/floatap_native_ratio.hpp
 *  \brief Correctly rounded positive unsigned ratios, without GMP.
 */

#ifndef ARIADNE_NUMERIC_FLOATAP_NATIVE_RATIO_HPP
#define ARIADNE_NUMERIC_FLOATAP_NATIVE_RATIO_HPP

#include "floatap_native_unsigned.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Ariadne::FloatAPNative {

// Signedness belongs to the floating-point value; inputs are magnitudes.
enum class RatioRounding { NearestEven, Downward, Upward, TowardZero };

struct RoundedRatio {
    Unsigned significand;
    std::int64_t exponent;
    // A carry out of the requested precision has already been shifted
    // away. Apply the exponent increment AFTER adding the caller's scale,
    // so that max-exponent carry can be converted into infinity.
    bool carry;
};

inline std::int64_t ratio_checked_add(std::int64_t a,std::int64_t b) {
    constexpr auto max=std::numeric_limits<std::int64_t>::max();
    constexpr auto min=std::numeric_limits<std::int64_t>::min();
    if ((b>0 && a>max-b) || (b<0 && a<min-b)) {
        throw std::overflow_error("FloatAP native ratio exponent overflow");
    }
    return a+b;
}

inline std::int64_t ratio_checked_sub(std::int64_t a,std::int64_t b) {
    if (b==std::numeric_limits<std::int64_t>::min()) {
        return ratio_checked_add(
            ratio_checked_add(a,std::numeric_limits<std::int64_t>::max()),1);
    }
    return ratio_checked_add(a,-b);
}

inline std::int64_t ratio_checked_bit_length(Unsigned const& value) {
    if (value.bit_length() >
        static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error("FloatAP native ratio size overflow");
    }
    return static_cast<std::int64_t>(value.bit_length());
}

inline std::size_t ratio_shift_amount(std::int64_t shift) {
    if (shift<0) {
        throw std::logic_error("negative FloatAP native ratio shift");
    }
    if constexpr (sizeof(std::size_t)<sizeof(std::uint64_t)) {
        if (static_cast<std::uint64_t>(shift)>
            std::numeric_limits<std::size_t>::max()) {
            throw std::overflow_error("FloatAP native ratio shift overflow");
        }
    }
    return static_cast<std::size_t>(shift);
}

// Return a normalized p-bit rounded significand and a binary exponent.
// The inputs must be strictly positive. The result represents
// n/d rounded at precision p, with the sign used only for directed
// rounding. The caller's external binary scale is not included.
inline RoundedRatio round_ratio(Unsigned const& n,Unsigned const& d,
                                unsigned long precision,bool negative,
                                RatioRounding mode)
{
    if (n.is_zero() || d.is_zero() || precision<2u ||
        static_cast<std::uint64_t>(precision)>
            static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
        throw std::invalid_argument("invalid FloatAP native ratio");
    }
    const std::int64_t p=static_cast<std::int64_t>(precision);
    std::int64_t k=ratio_checked_sub(ratio_checked_bit_length(n),
                                     ratio_checked_bit_length(d));
    if (k>=0) {
        if (compare(n,shift_left(d,ratio_shift_amount(k)))<0) { --k; }
    } else {
        if (compare(shift_left(n,ratio_shift_amount(ratio_checked_sub(0,k))),d)<0) {
            --k;
        }
    }
    const std::int64_t shift=ratio_checked_sub(p-1,k);
    const Unsigned numerator=shift>=0 ?
        shift_left(n,ratio_shift_amount(shift)) : n;
    const Unsigned denominator=shift<0 ?
        shift_left(d,ratio_shift_amount(ratio_checked_sub(0,shift))) : d;

    auto [quotient,remainder]=divrem(numerator,denominator);
    bool increment=false;
    if (!remainder.is_zero()) {
        switch (mode) {
            case RatioRounding::NearestEven: {
                const int twice_vs_divisor=
                    compare(shift_left(remainder,1u),denominator);
                increment=twice_vs_divisor>0 ||
                          (twice_vs_divisor==0 && quotient.bit(0u));
                break;
            }
            case RatioRounding::Downward:
                increment=negative;
                break;
            case RatioRounding::Upward:
                increment=!negative;
                break;
            case RatioRounding::TowardZero:
                break;
            default:
                throw std::invalid_argument("unknown FloatAP native rounding mode");
        }
    }
    if (increment) { quotient=add(quotient,Unsigned(1u)); }
    const bool carry=quotient.bit_length()>static_cast<std::size_t>(precision);
    if (carry) { quotient=shift_right(quotient,1u); }
    return {std::move(quotient),ratio_checked_sub(k,p-1),carry};
}

} // namespace Ariadne::FloatAPNative

#endif // ARIADNE_NUMERIC_FLOATAP_NATIVE_RATIO_HPP
