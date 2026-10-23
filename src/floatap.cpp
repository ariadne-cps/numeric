/***************************************************************************
 *            numeric/floatap.cpp
 *
 *  Copyright  2026  Luca Geretti
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

/*! \file numeric/floatap.cpp
 *  \brief Native arbitrary-precision floating-point operations.
 */

#include "numeric/floatap.hpp"
#include "numeric/integer.hpp"
#include "numeric/detail/integer_gmp_bridge.hpp"
#include "numeric/rational.hpp"
#include "numeric/dyadic.hpp"

#include "numeric/double.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/rounding.hpp"
#include "numeric/positive.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/decimal.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <istream>
#include <limits>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Ariadne {

// All unsigned arithmetic in FloatAP uses the repository's optimized
// 64-bit-limb BigUnsigned; the nonnegative significand is owned directly.
namespace {
bool floatap_test_bit(BigUnsigned const& n,std::size_t bit) noexcept {
    return ((n.limb(bit/64u)>>(bit%64u))&1u)!=0u;
}
bool floatap_any_bits_below(BigUnsigned const& n,std::size_t limit) noexcept {
    const auto digits=n.limbs();
    const std::size_t whole=std::min(limit/64u,digits.size());
    for(std::size_t i=0u;i<whole;++i) {
        if(digits[i]!=0u) return true;
    }
    const unsigned partial=static_cast<unsigned>(limit%64u);
    return partial!=0u && whole<digits.size() &&
        (digits[whole]&((BigUnsigned::Limb{1}<<partial)-1u))!=0u;
}
std::size_t floatap_popcount(BigUnsigned const& n) noexcept {
    std::size_t count=0u;
    for(auto word:n.limbs()) count+=static_cast<std::size_t>(std::popcount(word));
    return count;
}
BigUnsigned add_big_unsigned(BigUnsigned const&,BigUnsigned const&);
} // namespace
// Round a positive native magnitude in place, preserving tie-to-even.
void FloatAP::round_scaled_integer(FloatAP& result, Exponent scale, Rounding mode) {
    const auto precision = result._precision;
    const Exponent shift = checked_sub(bit_length(result._significand),
                                       static_cast<Exponent>(precision));
    result._exponent = checked_add(scale, shift);
    if (shift <= 0) {
        result._significand=shift_left(result._significand,
            static_cast<std::size_t>(shift_amount(checked_sub(0,shift))));
        return;
    }
    const std::size_t dropped=static_cast<std::size_t>(shift_amount(shift));
    const bool inexact=floatap_any_bits_below(result._significand,dropped);
    bool increment = false;
    if (inexact) {
        switch (mode.value()) {
            case Rounding::ToNearest:
                increment=floatap_test_bit(result._significand,dropped-1u) &&
                    (floatap_test_bit(result._significand,dropped) ||
                     floatap_any_bits_below(result._significand,dropped-1u));
                break;
            case Rounding::Downward:
                increment = result._negative;
                break;
            case Rounding::Upward:
                increment = !result._negative;
                break;
            case Rounding::TowardZero:
                break;
            default:
                throw std::invalid_argument("FloatAP unknown rounding mode");
        }
    }
    result._significand=shift_right(result._significand,dropped);
    if (increment) result._significand=add_big_unsigned(
        result._significand,BigUnsigned(1u));
    if (bit_length(result._significand)>static_cast<Exponent>(precision)) {
        result._significand=shift_right(result._significand,1u);
        // Crossing the maximal finite exponent by a rounded carry yields infinity.
        if (result._exponent == std::numeric_limits<Exponent>::max()) {
            result = infinity(precision, result._negative);
            return;
        }
        result._exponent = checked_add(result._exponent, 1);
    }
}

namespace {
// BigUnsigned's hidden-friend add relies on ADL, which is suppressed
// by FloatAP::add during lookup inside FloatAP member functions.
BigUnsigned add_big_unsigned(BigUnsigned const& a, BigUnsigned const& b) {
    return add(a,b);
}
// ADL for BigUnsigned's hidden friends is suppressed inside FloatAP::mul.
BigUnsigned mul_big_unsigned(BigUnsigned const& a, BigUnsigned const& b) {
    return mul(a,b);
}
BigUnsigned sub_big_unsigned(BigUnsigned const& a, BigUnsigned const& b) {
    return sub(a,b);
}
} // namespace

// D64: exact division and remainder now use the shared native
// BigUnsigned kernels, including for the old <=64-bit fast path.
// GMP is temporarily confined to the boundary wrapper and FloatAP
// significand storage until the remaining storage migration.
FloatAP FloatAP::from_native_ratio(Integer const& numerator,
                                  Integer const& denominator,
                                  Exponent scale, unsigned long precision,
                                  Rounding mode, bool negative_zero,
                                  bool invert_sign) {
    if (denominator.magnitude().is_zero()) {
        throw std::domain_error("FloatAP zero denominator");
    }
    FloatAP result(precision);
    auto const& n = numerator.magnitude();
    auto const& d = denominator.magnitude();
    if (n.is_zero()) {
        result.set_zero(negative_zero);
        return result;
    }
    result._negative =
        (numerator.negative() != denominator.negative()) != invert_sign;

    // Power-of-two denominators are exact binary scalings. Preserve
    // the fast path used by the old mpq_t constructor without creating a
    // GMP rational or executing a general multi-limb division.
    bool denominator_is_power_of_two = true;
    bool saw_nonzero_limb = false;
    for (auto limb : d.limbs()) {
        if (limb == 0u) continue;
        if (saw_nonzero_limb || !std::has_single_bit(limb)) {
            denominator_is_power_of_two = false;
            break;
        }
        saw_nonzero_limb = true;
    }
    if (denominator_is_power_of_two) {
        const Exponent denominator_shift =
            floatap_exponent_from_bit_length(d.bit_length()-1u);
        result._significand=n;
        round_scaled_integer(result,checked_sub(scale,denominator_shift),mode);
        return result;
    }

    Exponent k=checked_sub(floatap_exponent_from_bit_length(n.bit_length()),
                           floatap_exponent_from_bit_length(d.bit_length()));
    const auto nonnegative_shift=[](Exponent e) -> std::size_t {
        if (e<0 || static_cast<std::uintmax_t>(e)>
                     std::numeric_limits<std::size_t>::max()) {
            throw std::overflow_error("FloatAP native shift out of range");
        }
        return static_cast<std::size_t>(e);
    };
    if (k>=0) {
        if (compare(n,shift_left(d,nonnegative_shift(k)))<0) {
            k=checked_sub(k,1);
        }
    } else {
        if (compare(shift_left(n,nonnegative_shift(checked_sub(0,k))),d)<0) {
            k=checked_sub(k,1);
        }
    }

    Exponent const shift=checked_sub(static_cast<Exponent>(precision-1u),k);
    BigUnsigned dividend;
    BigUnsigned divisor;
    if (shift>=0) {
        dividend=shift_left(n,nonnegative_shift(shift));
        divisor=d;
    } else {
        dividend=n;
        divisor=shift_left(d,nonnegative_shift(checked_sub(0,shift)));
    }
    auto [quotient,remainder]=divrem(dividend,divisor);
    bool increment=false;
    if (!remainder.is_zero()) {
        switch (mode.value()) {
            case Rounding::ToNearest: {
                auto const twice=add_big_unsigned(remainder,remainder);
                int const order=compare(twice,divisor);
                increment=order>0 || (order==0 && (quotient.limb(0u)&1u));
                break;
            }
            case Rounding::Downward: increment=result._negative; break;
            case Rounding::Upward: increment=!result._negative; break;
            case Rounding::TowardZero: break;
            default: throw std::invalid_argument("FloatAP unknown rounding mode");
        }
    }
    if (increment) quotient=add_big_unsigned(quotient,BigUnsigned(1u));
    result._exponent=checked_add(
        scale,checked_sub(k,static_cast<Exponent>(precision-1u)));
    if (quotient.bit_length()>precision) {
        quotient=shift_right(quotient,1u);
        if (result._exponent==std::numeric_limits<Exponent>::max()) {
            return infinity(precision,result._negative);
        }
        result._exponent=checked_add(result._exponent,1);
    }
    result._significand=std::move(quotient);
    return result;
}

FloatAP FloatAP::negate(FloatAP x) {
    if (!x.is_nan()) x._negative = !x._negative;
    return x;
}

FloatAP FloatAP::add(FloatAP const& a, FloatAP const& b,
unsigned long p, Rounding mode) {
    return add_sub(a, b, p, mode, false);
}

FloatAP FloatAP::sub(FloatAP const& a, FloatAP const& b,
unsigned long p, Rounding mode) {
    return add_sub(a, b, p, mode, true);
}

FloatAP FloatAP::add_sub(FloatAP const& a, FloatAP const& b,
unsigned long p, Rounding mode, bool subtract) {
    const bool b_negative = b._negative != subtract;
    if (a.is_nan() || b.is_nan()) return nan(p);
    if (a.is_inf() || b.is_inf()) {
        if (a.is_inf() && b.is_inf() && a._negative != b_negative) return nan(p);
        return infinity(p, a.is_inf() ? a._negative : b_negative);
    }
    if (a.is_zero() && b.is_zero()) {
        return zero(p, a._negative == b_negative ? a._negative : mode == Rounding::Downward);
    }
    if (a.is_zero()) return FloatAP(subtract ? negate(b) : b, mode, PrecisionType(p));
    if (b.is_zero()) return FloatAP(a, mode, PrecisionType(p));

    // Align the exact integers for ordinary gaps; for immense gaps use
    // a rigorously small, signed sticky term before a single rounding.
    // Unlike the old approximation shortcut, this preserves directed
    // subtraction and ties, including when the tiny operand opposes the
    // large operand.
    const bool a_high = a._exponent >= b._exponent;
    FloatAP const& high = a_high ? a : b;
    FloatAP const& low = a_high ? b : a;
    const bool high_negative = a_high ? a._negative : b_negative;
    const bool low_negative = a_high ? b_negative : a._negative;
    // Unsigned subtraction handles even INT64_MAX - INT64_MIN exactly.
    const std::uint64_t gap = static_cast<std::uint64_t>(high._exponent) -
                              static_cast<std::uint64_t>(low._exponent);
    // This conservative gap bound places the smaller operand below all
    // noncoincident rounding boundaries even when p > high.precision.
    const std::uint64_t high_bits = static_cast<std::uint64_t>(high._precision);
    const std::uint64_t low_bits = static_cast<std::uint64_t>(low._precision);
    const bool far_apart = gap > high_bits && gap - high_bits > low_bits &&
        gap - high_bits - low_bits > static_cast<std::uint64_t>(p) + 2u;
    if (far_apart) {
        // Preserve the certified signed-sticky approximation using native
        // Integer for signed arithmetic and BigUnsigned for exact shifts.
        const Exponent padding=checked_add(
            checked_add(static_cast<Exponent>(high._precision),
                        static_cast<Exponent>(p)),2);
        const auto shift=static_cast<std::size_t>(shift_amount(padding));
        Integer exact=Integer::from_magnitude(
            shift_left(high._significand,shift),high_negative);
        exact=low_negative ? exact-Integer(1) : exact+Integer(1);
        const Integer denominator=Integer::from_magnitude(
            shift_left(BigUnsigned(1u),shift),false);
        return from_native_ratio(exact,denominator,high._exponent,p,mode);
    }

    // For ordinary gaps, align and combine positive magnitudes natively.
    // Preserve the independent sign metadata and a single final rounding.
    FloatAP result(p);
    result._negative=high_negative;
    BigUnsigned aligned=shift_left(high._significand,
        static_cast<std::size_t>(
            shift_amount(checked_sub(high._exponent,low._exponent))));
    if(a._negative==b_negative) {
        aligned=add_big_unsigned(aligned,low._significand);
    } else {
        if(compare(aligned,low._significand)<0) {
            aligned=sub_big_unsigned(low._significand,aligned);
            result._negative=!result._negative;
        } else {
            aligned=sub_big_unsigned(aligned,low._significand);
        }
    }
    if(aligned.is_zero()) {
        result.set_zero(mode==Rounding::Downward);
        return result;
    }
    result._significand=std::move(aligned);
    round_scaled_integer(result,low._exponent,mode);
    return result;
}

FloatAP FloatAP::mul(FloatAP const& a, FloatAP const& b,
unsigned long p, Rounding mode) {
    if (a.is_nan() || b.is_nan() || (a.is_inf() && b.is_zero()) || (a.is_zero() && b.is_inf())) return nan(p);
    const bool neg = a._negative != b._negative;
    if (a.is_inf() || b.is_inf()) return infinity(p, neg);
    if (a.is_zero() || b.is_zero()) return zero(p, neg);
    FloatAP result(p);
    result._negative = neg;
    result._significand=mul_big_unsigned(a._significand,b._significand);
    round_scaled_integer(result,checked_add(a._exponent,b._exponent),mode);
    return result;
}

FloatAP FloatAP::div(FloatAP const& a, FloatAP const& b,
unsigned long p, Rounding mode) {
    if (a.is_nan() || b.is_nan() || (a.is_zero() && b.is_zero()) || (a.is_inf() && b.is_inf())) return nan(p);
    const bool neg = a._negative != b._negative;
    if (a.is_inf() || b.is_zero()) return infinity(p, neg);
    if (a.is_zero() || b.is_inf()) return zero(p, neg);
    // Sign remains metadata, and inputs already own unsigned magnitudes.
    return from_native_ratio(
                             Integer::from_magnitude(a._significand,false),
                             Integer::from_magnitude(b._significand,false),
                             checked_sub(a._exponent,b._exponent),
                             p,mode,false,neg);
}

FloatAP FloatAP::fma(FloatAP const& a, FloatAP const& b,
FloatAP const& c, unsigned long p, Rounding mode) {
    if (a.is_nan() || b.is_nan() || c.is_nan() || (a.is_zero() && b.is_inf()) || (a.is_inf() && b.is_zero())) return nan(p);
    const bool neg = a._negative != b._negative;
    if (a.is_inf() || b.is_inf()) {
        if (c.is_inf() && c._negative != neg) return nan(p);
        return infinity(p, neg);
    }
    if (c.is_inf()) return infinity(p, c._negative);
    if (a.is_zero() || b.is_zero()) return add(zero(p, neg), c, p, mode);
    BigUnsigned product=mul_big_unsigned(a._significand,b._significand);
    FloatAP exact(static_cast<unsigned long>(product.bit_length()));
    exact._significand=std::move(product);
    exact._negative = neg;
    exact._exponent = checked_add(a._exponent, b._exponent);
    return add(exact, c, p, mode);
}

// Integer owns a native signed magnitude; no mpz_t import is required,
// including for INT64_MIN.
FloatAP::Float(std::int64_t n, PrecisionType pr)
    : Float(Integer(n), ROUND_TO_NEAREST, pr) { }

namespace {
thread_local FloatAP::RoundingModeType floatap_current_rounding =
    FloatAP::ROUND_TO_NEAREST;
thread_local std::size_t floatap_default_precision = 53;
}

FloatAP::RoundingModeType FloatAP::get_rounding_mode() noexcept {
    return floatap_current_rounding;
}

FloatAP::PrecisionType FloatAP::get_default_precision() {
    return PrecisionType(floatap_default_precision);
}

void FloatAP::set_default_precision(PrecisionType precision) {
    const FloatAP validated(precision);
    floatap_default_precision = validated.precision().bits();
}

void FloatAP::set_rounding_mode(RoundingModeType mode) {
    switch (mode.value()) {
        case Rounding::ToNearest:
        case Rounding::Downward:
        case Rounding::Upward:
        case Rounding::TowardZero:
            floatap_current_rounding = mode;
            return;
        default:
            throw std::invalid_argument("FloatAP unknown rounding mode");
    }
}

FloatAP::Rounding::Rounding(RoundToNearest const&) noexcept : _value(ToNearest) { }
FloatAP::Rounding::Rounding(RoundDownward const&) noexcept : _value(Downward) { }
FloatAP::Rounding::Rounding(RoundUpward const&) noexcept : _value(Upward) { }
FloatAP::Rounding::Rounding(RoundTowardZero const&) noexcept : _value(TowardZero) { }
FloatAP::Rounding::Rounding(RoundApproximately const&) noexcept : _value(ToNearest) { }
FloatAP::Rounding::Rounding(CurrentRoundingMode const&) noexcept
    : _value(FloatAP::get_rounding_mode().value()) { }

// Forward current-rounding operations to the existing native implementation.
FloatAP sqr(CurrentRoundingMode, FloatAP const& x) { return sqr(FloatAP::get_rounding_mode(),x); }
FloatAP rec(CurrentRoundingMode, FloatAP const& x) { return rec(FloatAP::get_rounding_mode(),x); }
FloatAP add(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2) { return add(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP sub(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2) { return sub(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP mul(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2) { return mul(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP div(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2) { return div(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP fma(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2, FloatAP const& x3) { return fma(FloatAP::get_rounding_mode(),x1,x2,x3); }
FloatAP pow(CurrentRoundingMode, FloatAP const& x, Nat m) { return pow(FloatAP::get_rounding_mode(),x,m); }
FloatAP pow(CurrentRoundingMode, FloatAP const& x, Int n) { return pow(FloatAP::get_rounding_mode(),x,n); }
FloatAP sqrt(CurrentRoundingMode, FloatAP const& x) { return sqrt(FloatAP::get_rounding_mode(),x); }
FloatAP exp(CurrentRoundingMode, FloatAP const& x) { return exp(FloatAP::get_rounding_mode(),x); }
FloatAP log(CurrentRoundingMode, FloatAP const& x) { return log(FloatAP::get_rounding_mode(),x); }
FloatAP sin(CurrentRoundingMode, FloatAP const& x) { return sin(FloatAP::get_rounding_mode(),x); }
FloatAP cos(CurrentRoundingMode, FloatAP const& x) { return cos(FloatAP::get_rounding_mode(),x); }
FloatAP tan(CurrentRoundingMode, FloatAP const& x) { return tan(FloatAP::get_rounding_mode(),x); }
FloatAP asin(CurrentRoundingMode, FloatAP const& x) { return asin(FloatAP::get_rounding_mode(),x); }
FloatAP acos(CurrentRoundingMode, FloatAP const& x) { return acos(FloatAP::get_rounding_mode(),x); }
FloatAP atan(CurrentRoundingMode, FloatAP const& x) { return atan(FloatAP::get_rounding_mode(),x); }
FloatAP FloatAP::pi(CurrentRoundingMode, ArbitraryPrecision pr) { return pi(FloatAP::get_rounding_mode(),pr); }
FloatAP add(CurrentRoundingMode, FloatAP const& x1, FloatDP const& x2) { return add(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP sub(CurrentRoundingMode, FloatAP const& x1, FloatDP const& x2) { return sub(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP mul(CurrentRoundingMode, FloatAP const& x1, FloatDP const& x2) { return mul(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP div(CurrentRoundingMode, FloatAP const& x1, FloatDP const& x2) { return div(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP add(CurrentRoundingMode, FloatDP const& x1, FloatAP const& x2) { return add(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP sub(CurrentRoundingMode, FloatDP const& x1, FloatAP const& x2) { return sub(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP mul(CurrentRoundingMode, FloatDP const& x1, FloatAP const& x2) { return mul(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP div(CurrentRoundingMode, FloatDP const& x1, FloatAP const& x2) { return div(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP add(CurrentRoundingMode, FloatAP const& x1, ExactDouble x2) { return add(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP sub(CurrentRoundingMode, FloatAP const& x1, ExactDouble x2) { return sub(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP mul(CurrentRoundingMode, FloatAP const& x1, ExactDouble x2) { return mul(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP div(CurrentRoundingMode, FloatAP const& x1, ExactDouble x2) { return div(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP add(CurrentRoundingMode, ExactDouble x1, FloatAP const& x2) { return add(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP sub(CurrentRoundingMode, ExactDouble x1, FloatAP const& x2) { return sub(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP mul(CurrentRoundingMode, ExactDouble x1, FloatAP const& x2) { return mul(FloatAP::get_rounding_mode(),x1,x2); }
FloatAP div(CurrentRoundingMode, ExactDouble x1, FloatAP const& x2) { return div(FloatAP::get_rounding_mode(),x1,x2); }

FloatAP::Float(FloatDP const& value, PrecisionType precision)
    : Float(Rational(value), ROUND_TO_NEAREST, precision) { }

FloatAP::Float(FloatDP const& value, RoundingModeType mode, PrecisionType precision)
    : Float(Rational(value), mode, precision) { }

FloatAP::Float(FloatAP const& value, PrecisionType precision)
    : Float(value, ROUND_TO_NEAREST, precision) { }

FloatAP::Float(FloatAP const& value, RoundingModeType mode, PrecisionType precision)
    : Float(precision) {
    if (value._precision == _precision) {
        *this = value;
        return;
    }
    if (value.is_nan()) {
        *this = nan(precision);
    } else if (value.is_inf()) {
        *this = infinity(static_cast<unsigned long>(precision.bits()), value.negative());
    } else if (value.is_zero()) {
        *this = zero(static_cast<unsigned long>(precision.bits()), value.negative());
    } else {
        // Rescale the existing binary significand without materializing
        // significand * 2^exponent as an enormous Rational. This is
        // essential at the native exponent limits.
        const Integer numerator=Integer::from_magnitude(
            value._significand,value._negative);
        *this=from_native_ratio(numerator,Integer(1u),value._exponent,
                                static_cast<unsigned long>(precision.bits()),
                                mode);
    }
}

FloatAP& FloatAP::operator=(Integer const& n) {
    return *this = FloatAP(n, get_rounding_mode(), precision());
}

FloatAP& FloatAP::operator=(Dyadic const& n) {
    // Preserve the exact-assignment contract of FloatMP. Constructing the
    // replacement first also leaves *this unchanged if it is inexact.
    return *this = FloatAP(n, precision());
}

FloatAP& FloatAP::operator=(TwoExp const& n) {
    return *this = FloatAP(n, precision());
}

FloatAP& FloatAP::operator=(ExactDouble const& n) {
    return *this = FloatAP(n, get_rounding_mode(), precision());
}

FloatAP& FloatAP::assign_unsigned(std::uint64_t n) {
    return *this = FloatAP(Integer(n), get_rounding_mode(), precision());
}

FloatAP::operator Dyadic() const {
    if (is_nan()) {
        return Dyadic::nan();
    }
    if (is_inf()) {
        return Dyadic::inf(negative() ? Sign::NEGATIVE : Sign::POSITIVE);
    }
    // Export the exact native binary significand/exponent.
    BigUnsigned magnitude=_significand;
    if (magnitude.is_zero()) { return Dyadic(0); }
    if (_exponent>=0) {
        return Dyadic(Integer::from_magnitude(
            shift_left(magnitude,shift_amount(_exponent)),_negative));
    }
    const auto shift=shift_amount(checked_sub(0,_exponent));
    return Dyadic(Integer::from_magnitude(std::move(magnitude),_negative),shift);
}

Integer cast_integer(FloatAP const& value) {
    if (!value.is_finite()) {
        throw std::domain_error("FloatAP cannot cast a nonfinite value to Integer");
    }
    const Rational exact = static_cast<Rational>(value);
    if (exact.denominator()!=1) {
        throw std::domain_error("FloatAP is not an integer");
    }
    return exact.numerator();
}

Comparison cmp(FloatAP const& a, Rational const& b) {
    if (a.is_nan() || is_nan(b)) {
        return Comparison::INCOMPARABLE;
    }
    return cmp(static_cast<Rational>(a), b);
}

Comparison cmp(Rational const& a, FloatAP const& b) {
    if (b.is_nan() || is_nan(a)) {
        return Comparison::INCOMPARABLE;
    }
    return cmp(a, static_cast<Rational>(b));
}

FloatAP::Float(ExactDouble const& value, PrecisionType precision)
    : Float(Rational(value), ROUND_TO_NEAREST, precision) { }

FloatAP::Float(Dyadic const& value, PrecisionType precision)
    : Float(Rational(value), ROUND_TO_NEAREST, precision) {
    // Unlike the explicitly rounded constructor, this conversion must be
    // exact. In particular, a Dyadic with more significant bits than the
    // requested precision cannot be silently rounded.
    if (!this->is_nan() && cmp(*this, Rational(value)) != Comparison::EQUAL) {
        throw std::runtime_error("FloatAP: Dyadic is not exactly representable at the requested precision");
    }
}

FloatAP::Float(TwoExp const& value, PrecisionType precision)
    : Float(Dyadic(value), precision) { }

FloatAP::Float(std::string const& value, PrecisionType precision)
    : Float(Dyadic(value), precision) { }

FloatAP::Float(ExactDouble const& value, RoundingModeType mode,
               PrecisionType precision)
    : Float(Rational(value), mode, precision) { }

FloatAP::Float(Dyadic const& value, RoundingModeType mode,
               PrecisionType precision)
    : Float(Rational(value), mode, precision) { }

void FloatAP::set_precision(PrecisionType precision) {
    if (precision.bits() == _precision) {
        return;
    }
    if (is_nan()) {
        *this = nan(precision);
    } else if (is_inf()) {
        *this = infinity(static_cast<unsigned long>(precision.bits()), _negative);
    } else if (is_zero()) {
        *this = zero(static_cast<unsigned long>(precision.bits()), _negative);
    } else {
        *this = FloatAP(static_cast<Rational>(*this), ROUND_TO_NEAREST, precision);
    }
}


Nat FloatAP::output_places=20;
void FloatAP::set_output_places(Nat places) { output_places=places; }

namespace {
// Format an exact binary rational to fixed decimal places using integer
// division and a single correctly directed decimal rounding step.
std::string decimal_fixed(FloatAP const& x, unsigned int places,
                          FloatAP::RoundingModeType mode) {
    if (x.is_nan()) return "nan";
    if (x.is_inf()) return x.negative() ? "-inf" : "inf";
    // Preserve the exact dyadic value rather than constructing a
    // canonical Rational and a temporary GMP rational. Decimal digits
    // are obtained by one native unsigned division and a single exact
    // directed rounding decision.
    BigUnsigned numerator=x.significand();
    BigUnsigned denominator(1u);
    const auto binary_exponent=x.exponent();
    if(binary_exponent>0) {
        numerator=shift_left(numerator,
            static_cast<std::size_t>(binary_exponent));
    } else if(binary_exponent<0) {
        // Avoid negating INT64_MIN; preserve the checked shift domain.
        const auto distance=static_cast<std::uint64_t>(
            -(binary_exponent+1))+1u;
        if(distance>std::numeric_limits<std::size_t>::max())
            throw std::overflow_error("FloatAP decimal exponent shift out of range");
        denominator=shift_left(denominator,static_cast<std::size_t>(distance));
    }
    const bool neg=x.negative();
    const BigUnsigned scale=pow(BigUnsigned(10u),
        static_cast<std::uint64_t>(places));
    auto [quotient,remainder]=divrem(mul(numerator,scale),denominator);
    if (!remainder.is_zero()) {
        bool increase=false;
        switch(mode.value()) {
            case FloatAP::Rounding::ToNearest: {
                const int order=compare(shift_left(remainder,1u),denominator);
                increase=order>0 || (order==0 &&
                    (quotient.limb(0u)&BigUnsigned::Limb{1})!=0u);
                break;
            }
            case FloatAP::Rounding::Downward: increase=neg; break;
            case FloatAP::Rounding::Upward: increase=!neg; break;
            case FloatAP::Rounding::TowardZero: break;
            default: throw std::invalid_argument("Invalid FloatAP decimal rounding mode");
        }
        if(increase) quotient=add(quotient,BigUnsigned(1u));
    }
    std::string digits=to_decimal(quotient);
    if(places>0) {
        if(digits.size()<=places) digits.insert(0,places+1-digits.size(),'0');
        digits.insert(digits.size()-places,1,'.');
    } else digits.push_back('.');
    if(neg) digits.insert(digits.begin(),'-');
    return digits;
}
}

std::string print(FloatAP const& x, DecimalPlaces places, FloatAP::RoundingModeType mode) {
    return decimal_fixed(x,static_cast<unsigned int>(places),mode);
}

std::string print(FloatAP const& x, DecimalPrecision figures, FloatAP::RoundingModeType mode) {
    if(!x.is_finite() || x.is_zero()) return decimal_fixed(x,0,mode);
    const double absolute=std::fabs(x.get_d());
    // DecimalPrecision is a count of significant figures rather than places.
    // For values outside the range of double, derive the base-10 exponent
    // from the full integer significand and binary exponent.
    int magnitude=0;
    if(std::isfinite(absolute) && absolute>0) {
        magnitude=static_cast<int>(std::floor(std::log10(absolute)));
    } else {
        const double estimated=(static_cast<double>(x.significand().bit_length())-1.0
                               +static_cast<double>(x.exponent()))*0.3010299956639812;
        magnitude=static_cast<int>(std::floor(estimated));
    }
    const long long places=static_cast<long long>(static_cast<unsigned int>(figures))-1-magnitude;
    return decimal_fixed(x,places>0?static_cast<unsigned int>(places):0u,mode);
}

OutputStream& write(OutputStream& out, FloatAP const& x, DecimalPlaces places, FloatAP::RoundingModeType mode) {
    return out<<print(x,places,mode);
}

OutputStream& write(OutputStream& out, FloatAP const& x, DecimalPrecision figures, FloatAP::RoundingModeType mode) {
    return out<<print(x,figures,mode);
}

OutputStream& repr(OutputStream& out, FloatAP const& x, FloatAP::RoundingModeType mode) {
    return write(out,x,DecimalPrecision(static_cast<Nat>(x.precision().bits()/3+4)),mode);
}

OutputStream& repr(OutputStream& out, FloatAP const& x) {
    return repr(out,x,FloatAP::get_rounding_mode());
}

std::string FloatAP::literal() const {
    return literal(get_rounding_mode());
}

std::string FloatAP::literal(RoundingModeType mode) const {
    std::ostringstream out;
    repr(out,*this,mode);
    return out.str();
}

OutputStream& operator<<(OutputStream& os, FloatAP const& value) {
    return os << static_cast<Rational>(value);
}

InputStream& operator>>(InputStream& is, FloatAP& value) {
    // Follow FloatMP's token boundaries, leaving punctuation for the caller.
    std::string token;
    int c = is.get();
    while (c == ' ' || c == '\t' || c == '\n') {
        c = is.get();
    }
    while ((c >= '0' && c <= '9') || c == '+' || c == '-' ||
           c == '.' || c == 'e' || c == 'E') {
        token.push_back(static_cast<char>(c));
        c = is.get();
    }
    const bool at_eof = c == std::char_traits<char>::eof();
    if (!at_eof) {
        is.putback(static_cast<char>(c));
    }

    try {
        const std::size_t exponent_start = token.find_first_of("eE");
        const String mantissa = token.substr(0, exponent_start);
        bool has_digit = false;
        for (char digit : mantissa) {
            if (digit >= '0' && digit <= '9') {
                has_digit = true;
                break;
            }
        }
        if (!has_digit) {
            throw std::invalid_argument("FloatAP invalid decimal literal");
        }
        Rational exact = static_cast<Rational>(Decimal(mantissa));
        if (exponent_start != std::string::npos) {
            const std::string suffix = token.substr(exponent_start + 1);
            std::size_t parsed = 0;
            const int exponent = std::stoi(suffix, &parsed);
            if (parsed != suffix.size()) {
                throw std::invalid_argument("FloatAP invalid decimal exponent");
            }
            exact *= pow(Rational(10), static_cast<Int>(exponent));
        }
        // Input is rounded once from its exact decimal value, independently
        // of the mutable arithmetic rounding mode.
        FloatAP parsed(exact, FloatAP::ROUND_TO_NEAREST, value.precision());
        value = std::move(parsed);
    } catch (const std::exception&) {
        is.setstate(std::ios::failbit);
        return is;
    }
    if (at_eof) {
        is.clear(is.rdstate() & ~std::ios::failbit);
    }
    return is;
}

Positive<FloatAP> abs(FloatAP const& x) {
    return cast_positive(x.negative() ? neg(x) : x);
}

PositiveUpperBound<FloatAP> mag(FloatAP const& x) {
    return PositiveUpperBound<FloatAP>(abs(x));
}

PositiveLowerBound<FloatAP> mig(FloatAP const& x) {
    return PositiveLowerBound<FloatAP>(abs(x));
}

// The immediate neighbors are computed on the p-bit integer significand;
// the exponent only changes on a carry or at a power-of-two boundary.
FloatAP next(RoundUpward, FloatAP const& x) {
    if (x.is_nan() || (x.is_inf() && !x.negative())) return x;
    if (x.is_inf()) return neg(FloatAP::max(x.precision()));
    if (x.is_zero()) return FloatAP::min(x.precision());
    FloatAP r(x);
    const bool all_ones = floatap_popcount(r._significand) == r._precision;
    const bool first = floatap_popcount(r._significand) == 1;
    if (r._negative) {
        if (first) {
            if (r._exponent == std::numeric_limits<FloatAP::Exponent>::min()) {
                return FloatAP::zero(r._precision, true);
            }
            r._significand=sub_big_unsigned(
                shift_left(BigUnsigned(1u),r._precision),BigUnsigned(1u));
            --r._exponent;
        } else {
            r._significand=sub_big_unsigned(r._significand,BigUnsigned(1u));
        }
    } else if (all_ones) {
        if (r._exponent == std::numeric_limits<FloatAP::Exponent>::max()) {
            return FloatAP::inf(r.precision());
        }
        r._significand=shift_left(BigUnsigned(1u),r._precision-1u);
        ++r._exponent;
    } else {
        r._significand=add_big_unsigned(r._significand,BigUnsigned(1u));
    }
    return r;
}

FloatAP next(RoundDownward, FloatAP const& x) {
    if (x.is_nan() || (x.is_inf() && x.negative())) return x;
    if (x.is_inf()) return FloatAP::max(x.precision());
    if (x.is_zero()) return neg(FloatAP::min(x.precision()));
    FloatAP r(x);
    const bool all_ones = floatap_popcount(r._significand) == r._precision;
    const bool first = floatap_popcount(r._significand) == 1;
    if (r._negative) {
        if (all_ones) {
            if (r._exponent == std::numeric_limits<FloatAP::Exponent>::max()) {
                return FloatAP::inf(Sign::NEGATIVE, r.precision());
            }
            r._significand=shift_left(BigUnsigned(1u),r._precision-1u);
            ++r._exponent;
        } else {
            r._significand=add_big_unsigned(r._significand,BigUnsigned(1u));
        }
    } else if (first) {
        if (r._exponent == std::numeric_limits<FloatAP::Exponent>::min()) {
            return FloatAP::zero(r._precision);
        }
        r._significand=sub_big_unsigned(
                shift_left(BigUnsigned(1u),r._precision),BigUnsigned(1u));
        --r._exponent;
    } else {
        r._significand=sub_big_unsigned(r._significand,BigUnsigned(1u));
    }
    return r;
}

// Integral rounding works directly on the BigUnsigned significand.
// In particular, signed minimum exponents never trigger huge 2^k allocations.
namespace {

bool has_discarded_bits(BigUnsigned const& magnitude,std::size_t shift) {
    const std::size_t whole=shift/64u;
    const std::size_t partial=shift%64u;
    const auto limbs=magnitude.limbs();
    for(std::size_t i=0;i<std::min(whole,limbs.size());++i) {
        if(limbs[i]!=0u) return true;
    }
    if(partial!=0u && whole<limbs.size()) {
        const std::uint64_t mask=(std::uint64_t{1}<<partial)-1u;
        return (limbs[whole]&mask)!=0u;
    }
    return false;
}

bool has_bit(BigUnsigned const& magnitude,std::size_t index) {
    return (magnitude.limb(index/64u) &
            (std::uint64_t{1}<<(index%64u)))!=0u;
}

FloatAP signed_integral_float(BigUnsigned magnitude,bool negative,
                              ArbitraryPrecision precision) {
    if(magnitude.is_zero()) {
        return FloatAP::zero(static_cast<unsigned long>(precision.bits()),negative);
    }
    return FloatAP(Integer::from_magnitude(std::move(magnitude),negative),
                   FloatAP::ROUND_TO_NEAREST,precision);
}

} // namespace

FloatAP floor(FloatAP const& x) {
    if(!x.is_finite() || x.is_zero() || x._exponent>=0) return x;
    const std::uint64_t shift=std::uint64_t{0}-
        static_cast<std::uint64_t>(x._exponent);
    const BigUnsigned& magnitude=x._significand;
    if(shift>=magnitude.bit_length()) {
        return FloatAP(x._negative ? std::int64_t{-1} : std::int64_t{0},x.precision());
    }
    const std::size_t dropped=static_cast<std::size_t>(shift);
    BigUnsigned integral=shift_right(magnitude,dropped);
    if(x._negative && has_discarded_bits(magnitude,dropped)) {
        integral=add_big_unsigned(integral,BigUnsigned(1u));
    }
    return signed_integral_float(std::move(integral),x._negative,x.precision());
}

FloatAP ceil(FloatAP const& x) {
    if(!x.is_finite() || x.is_zero() || x._exponent>=0) return x;
    const std::uint64_t shift=std::uint64_t{0}-
        static_cast<std::uint64_t>(x._exponent);
    const BigUnsigned& magnitude=x._significand;
    if(shift>=magnitude.bit_length()) {
        return FloatAP(x._negative ? std::int64_t{0} : std::int64_t{1},x.precision());
    }
    const std::size_t dropped=static_cast<std::size_t>(shift);
    BigUnsigned integral=shift_right(magnitude,dropped);
    if(!x._negative && has_discarded_bits(magnitude,dropped)) {
        integral=add_big_unsigned(integral,BigUnsigned(1u));
    }
    return signed_integral_float(std::move(integral),x._negative,x.precision());
}

FloatAP round(FloatAP const& x) {
    if(!x.is_finite() || x.is_zero() || x._exponent>=0) return x;
    const std::uint64_t shift=std::uint64_t{0}-
        static_cast<std::uint64_t>(x._exponent);
    const BigUnsigned& magnitude=x._significand;
    const std::size_t bits=magnitude.bit_length();
    if(shift>bits) {
        return FloatAP::zero(x._precision,x._negative);
    }
    const std::size_t dropped=static_cast<std::size_t>(shift);
    BigUnsigned integral=shift_right(magnitude,dropped);
    // Preserve the historical nearest-integer/ties-away-from-zero contract.
    if(has_bit(magnitude,dropped-1u)) {
        integral=add_big_unsigned(integral,BigUnsigned(1u));
    }
    return signed_integral_float(std::move(integral),x._negative,x.precision());
}

FloatAP abs(FloatAP::RoundingModeType, FloatAP const& x) {
    return x.negative() ? neg(x) : x;
}

FloatAP mag(FloatAP::RoundingModeType mode, FloatAP const& x) {
    return abs(mode, x);
}

FloatAP& iadd(FloatAP::RoundingModeType mode, FloatAP& a, FloatAP const& b) {
    return a = FloatAP::add(a, b, static_cast<unsigned long>(a.precision().bits()), mode);
}

FloatAP& isub(FloatAP::RoundingModeType mode, FloatAP& a, FloatAP const& b) {
    return a = FloatAP::sub(a, b, static_cast<unsigned long>(a.precision().bits()), mode);
}

FloatAP& imul(FloatAP::RoundingModeType mode, FloatAP& a, FloatAP const& b) {
    return a = FloatAP::mul(a, b, static_cast<unsigned long>(a.precision().bits()), mode);
}

FloatAP& idiv(FloatAP::RoundingModeType mode, FloatAP& a, FloatAP const& b) {
    return a = FloatAP::div(a, b, static_cast<unsigned long>(a.precision().bits()), mode);
}

FloatAP add(FloatAP::RoundingModeType mode, FloatAP const& a, FloatDP const& b) {
    return FloatAP(Rational(a) + Rational(b), mode, a.precision());
}

FloatAP add(FloatAP::RoundingModeType mode, FloatDP const& a, FloatAP const& b) {
    return FloatAP(Rational(a) + Rational(b), mode, b.precision());
}

FloatAP sub(FloatAP::RoundingModeType mode, FloatAP const& a, FloatDP const& b) {
    return FloatAP(Rational(a) - Rational(b), mode, a.precision());
}

FloatAP sub(FloatAP::RoundingModeType mode, FloatDP const& a, FloatAP const& b) {
    return FloatAP(Rational(a) - Rational(b), mode, b.precision());
}

FloatAP mul(FloatAP::RoundingModeType mode, FloatAP const& a, FloatDP const& b) {
    return FloatAP(Rational(a) * Rational(b), mode, a.precision());
}

FloatAP mul(FloatAP::RoundingModeType mode, FloatDP const& a, FloatAP const& b) {
    return FloatAP(Rational(a) * Rational(b), mode, b.precision());
}

FloatAP div(FloatAP::RoundingModeType mode, FloatAP const& a, FloatDP const& b) {
    return FloatAP(Rational(a) / Rational(b), mode, a.precision());
}

FloatAP div(FloatAP::RoundingModeType mode, FloatDP const& a, FloatAP const& b) {
    return FloatAP(Rational(a) / Rational(b), mode, b.precision());
}

FloatAP add(FloatAP::RoundingModeType mode, FloatAP const& a, ExactDouble const& b) {
    return FloatAP(Rational(a) + Rational(b), mode, a.precision());
}

FloatAP add(FloatAP::RoundingModeType mode, ExactDouble const& a, FloatAP const& b) {
    return FloatAP(Rational(a) + Rational(b), mode, b.precision());
}

FloatAP sub(FloatAP::RoundingModeType mode, FloatAP const& a, ExactDouble const& b) {
    return FloatAP(Rational(a) - Rational(b), mode, a.precision());
}

FloatAP sub(FloatAP::RoundingModeType mode, ExactDouble const& a, FloatAP const& b) {
    return FloatAP(Rational(a) - Rational(b), mode, b.precision());
}

FloatAP mul(FloatAP::RoundingModeType mode, FloatAP const& a, ExactDouble const& b) {
    return FloatAP(Rational(a) * Rational(b), mode, a.precision());
}

FloatAP mul(FloatAP::RoundingModeType mode, ExactDouble const& a, FloatAP const& b) {
    return FloatAP(Rational(a) * Rational(b), mode, b.precision());
}

FloatAP div(FloatAP::RoundingModeType mode, FloatAP const& a, ExactDouble const& b) {
    return FloatAP(Rational(a) / Rational(b), mode, a.precision());
}

FloatAP div(FloatAP::RoundingModeType mode, ExactDouble const& a, FloatAP const& b) {
    return FloatAP(Rational(a) / Rational(b), mode, b.precision());
}


FloatAP operator*(TwoExp const& t, FloatAP const& x) {
    return shft(x, t.exponent());
}

FloatAP operator*(FloatAP const& x, TwoExp const& t) {
    return shft(x, t.exponent());
}

FloatAP operator/(FloatAP const& x, TwoExp const& t) {
    return shft(x, -t.exponent());
}

FloatAP mul(FloatAP const& x, TwoExp const& t) {
    return shft(x, t.exponent());
}

FloatAP div(FloatAP const& x, TwoExp const& t) {
    return shft(x, -t.exponent());
}

FloatAP& operator*=(FloatAP& x, TwoExp const& t) {
    return x = x * t;
}

FloatAP& operator/=(FloatAP& x, TwoExp const& t) {
    return x = x / t;
}

// Native fused multiply-add: one outward rounding per endpoint.
Bounds<FloatAP> fma(FloatAP const& a, FloatAP const& b, FloatAP const& c) {
    return Bounds<FloatAP>(fma(FloatAP::ROUND_DOWNWARD,a,b,c),
                           fma(FloatAP::ROUND_UPWARD,a,b,c));
}

// Other outward-rounded FloatAP operations are defined by Bounds<F> in
// float_bounds.hpp; do not define competing overloads here.
FloatAP pow(FloatAP::RoundingModeType mode, FloatAP const& x, Int exponent) {
    if (exponent >= 0) {
        return pow(mode, x, static_cast<Nat>(exponent));
    }
    const auto bits = static_cast<unsigned long>(x.precision().bits());
    const bool odd = (exponent % 2) != 0;
    if (x.is_nan()) return FloatAP::nan(bits);
    if (x.is_inf()) return FloatAP::zero(bits, x.negative() && odd);
    if (x.is_zero()) return FloatAP::infinity(bits, x.negative() && odd);
    // The reciprocal of a rounded power need not be correctly rounded:
    // compute the exact rational power, then round once.
    const Rational exact = pow(static_cast<Rational>(x), exponent);
    return FloatAP(exact, mode, x.precision());
}

Comparison cmp(FloatAP const& a, FloatDP const& b) {
    return cmp(a, Rational(b));
}

Comparison cmp(FloatDP const& a, FloatAP const& b) {
    return cmp(Rational(a), b);
}

bool operator==(FloatAP const& a, FloatDP const& b) { return cmp(a,b)==Comparison::EQUAL; }
bool operator!=(FloatAP const& a, FloatDP const& b) { return cmp(a,b)!=Comparison::EQUAL; }
bool operator< (FloatAP const& a, FloatDP const& b) { return cmp(a,b)==Comparison::LESS; }
bool operator> (FloatAP const& a, FloatDP const& b) { return cmp(a,b)==Comparison::GREATER; }
bool operator<=(FloatAP const& a, FloatDP const& b) {
    auto c=cmp(a,b); return c==Comparison::LESS || c==Comparison::EQUAL;
}

bool operator>=(FloatAP const& a, FloatDP const& b) {
    auto c=cmp(a,b); return c==Comparison::GREATER || c==Comparison::EQUAL;
}

bool operator==(FloatDP const& a, FloatAP const& b) { return cmp(a,b)==Comparison::EQUAL; }
bool operator!=(FloatDP const& a, FloatAP const& b) { return cmp(a,b)!=Comparison::EQUAL; }
bool operator< (FloatDP const& a, FloatAP const& b) { return cmp(a,b)==Comparison::LESS; }
bool operator> (FloatDP const& a, FloatAP const& b) { return cmp(a,b)==Comparison::GREATER; }
bool operator<=(FloatDP const& a, FloatAP const& b) {
    auto c=cmp(a,b); return c==Comparison::LESS || c==Comparison::EQUAL;
}

bool operator>=(FloatDP const& a, FloatAP const& b) {
    auto c=cmp(a,b); return c==Comparison::GREATER || c==Comparison::EQUAL;
}

FloatAP shft(FloatAP const& value, Int amount) {
    FloatAP result(value);
    if (result.is_finite() && !result.is_zero()) {
        result._exponent = FloatAP::checked_add(
            result._exponent, static_cast<FloatAP::Exponent>(amount));
    }
    return result;
}

FloatAP::Float(Integer const& n, RoundingModeType mode, PrecisionType pr)
: Float(from_native_ratio(n,Integer(1u),0,
                          static_cast<unsigned long>(pr.bits()),mode)) {
}

namespace {
// Rational declares is_nan/is_inf as friend functions: ADL finds them,
// whereas qualified lookup (Ariadne::is_nan) does not.
bool rational_is_nan(Rational const& q) {
    return is_nan(q);
}

bool rational_is_inf(Rational const& q) {
    return is_inf(q);
}
}

FloatAP FloatAP::from_ratio(Integer const& numerator,
                            Integer const& denominator,
                            RoundingModeType mode,PrecisionType precision) {
    return from_native_ratio(numerator,denominator,0,
                             static_cast<unsigned long>(precision.bits()),mode);
}

FloatAP::Float(Rational const& q, RoundingModeType mode, PrecisionType pr)
: Float(pr) {
    if (rational_is_nan(q)) {
        *this=nan(_precision);
    }
    else if (rational_is_inf(q)) {
        *this=infinity(_precision,q.numerator().negative());
    } else {
        // Rational already owns canonical Integer/BigUnsigned values.
        // Avoid exporting them to GMP merely to import them back.
        *this=from_native_ratio(q.numerator(),q.denominator(),0,_precision,mode);
    }
}

FloatAP::operator Rational() const {
    if (is_nan()) return Rational::nan();
    if (is_inf()) return Rational::inf(negative() ? Sign::NEGATIVE : Sign::POSITIVE);
    BigUnsigned n=_significand;
    BigUnsigned d(1u);
    if (_exponent>0) { n=shift_left(n,shift_amount(_exponent)); }
    if (_exponent<0) { d=shift_left(d,shift_amount(checked_sub(0,_exponent))); }
    return Rational(Integer::from_magnitude(std::move(n),_negative),Integer::from_magnitude(std::move(d),false));
}

double FloatAP::get_d() const {
    // FloatAP stores an exact binary significand and an exponent. Converting
    // through Rational builds/canonicalizes a potentially enormous 2^k
    // denominator and also silently truncates rather than respecting the
    // active rounding mode. Extract 53 leading bits directly from GMP.
    if (is_nan()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (is_inf()) {
        const double infinity = std::numeric_limits<double>::infinity();
        return _negative ? -infinity : infinity;
    }
    if (is_zero()) {
        return _negative ? -0.0 : 0.0;
    }

    constexpr Exponent digits = std::numeric_limits<double>::digits;
    constexpr Exponent largest = std::numeric_limits<double>::max_exponent - 1;
    constexpr Exponent min_normal = std::numeric_limits<double>::min_exponent - 1;
    constexpr Exponent min_subnormal = min_normal - digits + 1;
    const Rounding::Value mode = get_rounding_mode().value();
    const auto signed_value = [this](double value) {
        return _negative ? -value : value;
    };
    const auto rounds_away = [this,mode]() {
        return mode == Rounding::ToNearest ||
               (mode == Rounding::Upward && !_negative) ||
               (mode == Rounding::Downward && _negative);
    };

    // The finite significand is positive and has at most _precision bits,
    // which itself fits Exponent by construction.
    const Exponent top=static_cast<Exponent>(_significand.bit_length()-1u);
    if (_exponent > 0 && top > std::numeric_limits<Exponent>::max() - _exponent) {
        return signed_value(rounds_away() ?
                            std::numeric_limits<double>::infinity() :
                            std::numeric_limits<double>::max());
    }
    const Exponent highest = _exponent + top;
    if (highest > largest) {
        return signed_value(rounds_away() ?
                            std::numeric_limits<double>::infinity() :
                            std::numeric_limits<double>::max());
    }
    if (highest < min_subnormal - 1) {
        const bool directed_away = (mode == Rounding::Upward && !_negative) ||
                                   (mode == Rounding::Downward && _negative);
        return signed_value(directed_away ?
                            std::numeric_limits<double>::denorm_min() : 0.0);
    }

    // Extract the 53 leading bits exactly, truncating any low bits.
    const std::size_t width=_significand.bit_length();
    const BigUnsigned top_bits=width>static_cast<std::size_t>(digits)
        ? shift_right(_significand,width-static_cast<std::size_t>(digits))
        : shift_left(_significand,static_cast<std::size_t>(digits)-width);
    const std::uint64_t leading=top_bits.limb(0u);
    const Exponent subnormal_shift = highest < min_normal
                                    ? min_normal - highest : 0;
    std::uint64_t mantissa = leading >> static_cast<unsigned>(subnormal_shift);

    // The original significand bits below this cut decide every rounding
    // mode, including IEEE halfway ties and the transition to subnormals.
    const Exponent discarded = top - (digits - 1) + subnormal_shift;
    if (discarded > 0) {
        const std::size_t dropped=static_cast<std::size_t>(discarded);
        bool increment = false;
        switch (mode) {
            case Rounding::ToNearest:
                increment=floatap_test_bit(_significand,dropped-1u) &&
                    ((mantissa&1u)!=0u ||
                     floatap_any_bits_below(_significand,dropped-1u));
                break;
            case Rounding::Downward:
                increment=_negative &&
                    floatap_any_bits_below(_significand,dropped);
                break;
            case Rounding::Upward:
                increment=!_negative &&
                    floatap_any_bits_below(_significand,dropped);
                break;
            case Rounding::TowardZero:
                break;
            default:
                throw std::invalid_argument("FloatAP unknown double conversion mode");
        }
        if (increment) {
            ++mantissa;
        }
    }

    if (highest >= min_normal) {
        // A carry may move the leading bit into the next exponent, or
        // overflow the largest finite double after rounding.
        if (mantissa == (std::uint64_t{1} << digits)) {
            mantissa >>= 1u;
            if (highest == largest) {
                return signed_value(rounds_away() ?
                                    std::numeric_limits<double>::infinity() :
                                    std::numeric_limits<double>::max());
            }
            return signed_value(std::ldexp(
                static_cast<double>(mantissa),
                static_cast<int>(highest - digits + 2)));
        }
        return signed_value(std::ldexp(
            static_cast<double>(mantissa),
            static_cast<int>(highest - digits + 1)));
    }
    // Integer multiples of the least subnormal are exactly representable.
    // mantissa == 2^52 also correctly produces the minimum normal double.
    return signed_value(std::ldexp(
        static_cast<double>(mantissa), static_cast<int>(min_subnormal)));
}

Comparison cmp(FloatAP const& a, FloatAP const& b) {
    if (a.is_nan() || b.is_nan()) {
        return Comparison::INCOMPARABLE;
    }
    if (a.is_inf() || b.is_inf()) {
        if (a.is_inf() && b.is_inf()) {
            return a.negative() == b.negative() ? Comparison::EQUAL :
                   (a.negative() ? Comparison::LESS : Comparison::GREATER);
        }
        return a.is_inf() ? (a.negative() ? Comparison::LESS : Comparison::GREATER) :
               (b.negative() ? Comparison::GREATER : Comparison::LESS);
    }

    // Signed zeros compare equal, independent of their stored sign.
    if (a.is_zero() || b.is_zero()) {
        if (a.is_zero() && b.is_zero()) {
            return Comparison::EQUAL;
        }
        if (a.is_zero()) {
            return b.negative() ? Comparison::GREATER : Comparison::LESS;
        }
        return a.negative() ? Comparison::LESS : Comparison::GREATER;
    }

    if (a._negative != b._negative) {
        return a._negative ? Comparison::LESS : Comparison::GREATER;
    }

    // Nonzero normalized operands of equal precision have the same
    // significand bit-length. Their exponents alone decide magnitude
    // whenever they differ, regardless of the significand values.
    // This avoids GMP size/shift work on the common equal-precision case.
    if (a._precision == b._precision && a._exponent != b._exponent) {
        const bool greater = (a._exponent > b._exponent) != a._negative;
        return greater ? Comparison::GREATER : Comparison::LESS;
    }

    // Compare significand * 2^exponent directly. A different most
    // significant binary position decides the order without allocating
    // any GMP temporaries or constructing Rational values.
    int magnitude = 0;
    if (a._exponent == b._exponent) {
        magnitude = compare(a._significand,b._significand);
    } else {
        const std::size_t a_bits=a._significand.bit_length();
        const std::size_t b_bits=b._significand.bit_length();

        if (a._exponent > b._exponent) {
            // Unsigned subtraction also handles a signed exponent gap
            // larger than INT64_MAX without undefined behavior.
            const std::uint64_t gap =
                static_cast<std::uint64_t>(a._exponent) -
                static_cast<std::uint64_t>(b._exponent);

            if (a_bits >= b_bits || gap > b_bits - a_bits) {
                magnitude = 1;
            } else if (gap < b_bits - a_bits) {
                magnitude = -1;
            } else {
                // The leading positions coincide. Shift the smaller
                // representation right; the discarded bits break ties.
                const auto shift=static_cast<std::size_t>(gap);
                const BigUnsigned truncated=shift_right(b._significand,shift);
                magnitude=compare(a._significand,truncated);
                if(magnitude==0 && floatap_any_bits_below(b._significand,shift)) {
                    magnitude = -1;
                }
            }
        } else {
            const std::uint64_t gap =
                static_cast<std::uint64_t>(b._exponent) -
                static_cast<std::uint64_t>(a._exponent);

            if (b_bits >= a_bits || gap > a_bits - b_bits) {
                magnitude = -1;
            } else if (gap < a_bits - b_bits) {
                magnitude = 1;
            } else {
                const auto shift=static_cast<std::size_t>(gap);
                const BigUnsigned truncated=shift_right(a._significand,shift);
                magnitude=compare(truncated,b._significand);
                if(magnitude==0 && floatap_any_bits_below(a._significand,shift)) {
                    magnitude = 1;
                }
            }
        }
    }

    if (a._negative) {
        magnitude = -magnitude;
    }
    return magnitude < 0 ? Comparison::LESS :
           (magnitude > 0 ? Comparison::GREATER : Comparison::EQUAL);
}

// The normalized significand has precisely p bits. After scaling it
// by p-1 or p binary positions, its integer square root also has p
// bits. The integer remainder determines the exact rounding direction.
// This avoids adaptive approximations and any artificial iteration cap.
FloatAP sqrt(FloatAP::RoundingModeType mode, FloatAP const& value) {
    const unsigned long bits = static_cast<unsigned long>(value.precision().bits());
    if (value.is_nan() || (value.negative() && !value.is_zero())) {
        return FloatAP::nan(bits);
    }
    if (value.is_inf() || value.is_zero()) {
        return value;
    }

    const unsigned long offset = bits - 1;
    const bool odd_exponent = value._exponent % 2 != 0;
    const bool odd_magnitude = (offset % 2 != 0) != odd_exponent;
    const unsigned long shift = offset + (odd_magnitude ? 1ul : 0ul);

    // Compute floor((exponent + p - 1)/2) without overflowing int64_t
    // and with mathematical floor rather than truncation toward zero.
    const FloatAP::Exponent half_exponent =
        value._exponent / 2 - ((value._exponent < 0 && odd_exponent) ? 1 : 0);
    FloatAP::Exponent result_exponent = FloatAP::checked_add(
        half_exponent, static_cast<FloatAP::Exponent>(offset / 2));
    if (odd_exponent && offset % 2 != 0) {
        result_exponent = FloatAP::checked_add(result_exponent, 1);
    }
    result_exponent = FloatAP::checked_sub(
        result_exponent, static_cast<FloatAP::Exponent>(offset));

    // Compute the exact integer square root and remainder with native
    // unsigned limbs. GMP is not part of the root/rounding calculation.
    // Only the final owned significand export remains transitional.
    const BigUnsigned& magnitude=value._significand;
    const BigUnsigned scaled=shift_left(magnitude,static_cast<std::size_t>(
        FloatAP::shift_amount(static_cast<FloatAP::Exponent>(shift))));
    auto [root,remainder]=sqrtrem(scaled);

    bool increment = false;
    if (!remainder.is_zero()) {
        switch (mode.value()) {
            case FloatAP::Rounding::ToNearest:
                // scaled=root*root+remainder. The exact nearest criterion
                // is remainder>root; a halfway tie is impossible.
                increment = compare(remainder,root)>0;
                break;
            case FloatAP::Rounding::Upward:
                increment = true;
                break;
            case FloatAP::Rounding::Downward:
            case FloatAP::Rounding::TowardZero:
                break;
            default:
                throw std::invalid_argument("FloatAP unknown square-root rounding mode");
        }
    }
    if (increment) {
        root=add_big_unsigned(root,BigUnsigned(1u));
    }
    if (root.bit_length()>bits) {
        root=shift_right(root,1u);
        result_exponent = FloatAP::checked_add(result_exponent, 1);
    }

    FloatAP result(bits);
    result._significand=std::move(root);
    result._exponent = result_exponent;
    return result;
}

FloatAP pow(FloatAP::RoundingModeType mode, FloatAP const& value, Nat exponent) {
    const unsigned long bits = static_cast<unsigned long>(value.precision().bits());
    const unsigned long degree = static_cast<unsigned long>(exponent);
    if (degree == 0) {
        return FloatAP(std::int64_t{1}, value.precision());
    }
    if (value.is_nan()) {
        return FloatAP::nan(bits);
    }
    if (value.is_inf()) {
        return FloatAP::infinity(bits, value.negative() && (degree % 2 != 0));
    }
    if (value.is_zero()) {
        return FloatAP::zero(bits, value.negative() && (degree % 2 != 0));
    }

    // The operand is exactly sign * significand * 2^exponent.
    // A nonnegative integer power therefore has an exact integer
    // significand and a binary scale. Round ONCE, after powering:
    // multiplying already-rounded Floats would introduce double rounding.
    if (degree == 1ul) {
        return value;
    }

    // Calculate exponent*degree with checked binary multiplication.
    // Exception: the unnormalized product may exceed Exponent even when
    // normalization cancels it (e.g. very large powers of exact one).
    // Preserve the old general rational path in that rare situation.
    FloatAP::Exponent product_scale = 0;
    FloatAP::Exponent doubled_scale = value._exponent;
    unsigned long remaining_degree = degree;
    bool scale_overflow = false;
    try {
        while (remaining_degree != 0ul) {
            if ((remaining_degree & 1ul) != 0ul) {
                product_scale = FloatAP::checked_add(
                    product_scale, doubled_scale);
            }
            remaining_degree >>= 1u;
            if (remaining_degree != 0ul) {
                doubled_scale = FloatAP::checked_add(
                    doubled_scale, doubled_scale);
            }
        }
    } catch (std::overflow_error const&) {
        scale_overflow = true;
    }

    if (!scale_overflow) {
        FloatAP result(bits);
        result._negative = value._negative && ((degree & 1ul) != 0ul);
        // Exact power on the already-owned unsigned magnitude.
        result._significand=pow(value._significand,
            static_cast<std::uint64_t>(degree));
        FloatAP::round_scaled_integer(result, product_scale, mode);
        return result;
    }

    // Rare normalization-cancellation/scale-overflow fallback.
    // Keep an exact rational result and one final rounding, without
    // roundtrips through mpq_t/mpz_t. Rational::pow is native.
    const Rational exact=static_cast<Rational>(value);
    return FloatAP(pow(exact,exponent),mode,value.precision());
}

template<> String class_name<ArbitraryPrecision>() { return "ArbitraryPrecision"; }
template<> String class_name<FloatAP>() { return "FloatAP"; }
template<> String class_name<Rounded<FloatAP>>() { return "Rounded<FloatAP>"; }


// Exact rational intervals shared by certified transcendental operations.
namespace {

class Fraction {
    mpq_t _value;

public:
    Fraction() {
        mpq_init(_value);
    }

    explicit Fraction(long integer) : Fraction() {
        mpq_set_si(_value, integer, 1);
    }

    explicit Fraction(mpq_srcptr rational) : Fraction() {
        mpq_set(_value, rational);
    }

    // Certified interval scratch is still GMP-backed until the
    // separate native transcendental migration. Export the exact
    // binary dyadic directly: constructing a canonical Rational here
    // needlessly invokes a potentially large GCD at every precision.
    explicit Fraction(FloatAP const& binary) : Fraction() {
        if (!binary.is_finite()) {
            throw std::domain_error("FloatAP nonfinite to rational interval");
        }
        integer_to_gmp(mpq_numref(_value),
            Integer::from_magnitude(binary.significand(),binary.negative()));
        mpz_set_ui(mpq_denref(_value),1u);
        const auto exponent=binary.exponent();
        if (exponent==std::numeric_limits<FloatAP::Exponent>::min()) {
            throw std::overflow_error("FloatAP exponent shift overflow");
        }
        const std::uint64_t distance=exponent<0
            ? static_cast<std::uint64_t>(-exponent)
            : static_cast<std::uint64_t>(exponent);
        if (distance>std::numeric_limits<mp_bitcnt_t>::max()) {
            throw std::overflow_error("FloatAP GMP interval shift overflow");
        }
        const mp_bitcnt_t shift=static_cast<mp_bitcnt_t>(distance);
        if (exponent>0) {
            mpz_mul_2exp(mpq_numref(_value),mpq_numref(_value),shift);
        } else if (exponent<0) {
            mpz_mul_2exp(mpq_denref(_value),mpq_denref(_value),shift);
        }
        mpq_canonicalize(_value);
    }

    // Import a native fixed-point endpoint once; no GMP in the series loop.
    Fraction(BigUnsigned const& numerator, unsigned long denominator_bits) : Fraction() {
        const auto digits=numerator.limbs();
        mpz_import(mpq_numref(_value),digits.size(),-1,
                   sizeof(BigUnsigned::Limb),0,0,digits.data());
        mpz_set_ui(mpq_denref(_value),1u);
        mpz_mul_2exp(mpq_denref(_value),mpq_denref(_value),denominator_bits);
        mpq_canonicalize(_value);
    }

    // A certified dyadic endpoint may straddle zero; attach its sign
    // only once at the Fraction boundary, not inside the native series.
    Fraction(BigUnsigned const& magnitude, unsigned long denominator_bits,
             bool negative) : Fraction(magnitude,denominator_bits) {
        if (negative && !magnitude.is_zero()) mpq_neg(_value,_value);
    }

    Fraction(mpz_srcptr numerator, unsigned long denominator_bits) : Fraction() {
        mpz_set(mpq_numref(_value), numerator);
        mpz_set_ui(mpq_denref(_value), 1);
        mpz_mul_2exp(mpq_denref(_value), mpq_denref(_value), denominator_bits);
        mpq_canonicalize(_value);
    }

    Fraction(Fraction const& other) : Fraction(other._value) { }

    Fraction(Fraction&& other) noexcept : Fraction() {
        mpq_swap(_value, other._value);
    }

    ~Fraction() {
        mpq_clear(_value);
    }

    Fraction& operator=(Fraction other) noexcept {
        mpq_swap(_value, other._value);
        return *this;
    }

    mpq_srcptr raw() const noexcept {
        return _value;
    }

    // Round exact interval endpoints directly: do not construct a
    // canonical Rational. Rational(n,d) performs a multiprecision GCD
    // even though the certified interval is already normalized.
    // It was introduced into every transcendental certification
    // attempt in D120 and caused a substantial regression.
    FloatAP rounded(unsigned long bits,FloatAP::RoundingModeType mode) const {
        return FloatAP::from_ratio(
            integer_from_gmp(mpq_numref(_value)),
            integer_from_gmp(mpq_denref(_value)),mode,
            FloatAP::PrecisionType(bits));
    }

    // Avoid a temporary Fraction and general rational division when
    // reducing an argument by an exact factor of two.
    void halve() {
        mpq_div_2exp(_value, _value, 1u);
    }

    int sign() const noexcept {
        return mpq_sgn(_value);
    }

    friend int compare(Fraction const& a, Fraction const& b) {
        return mpq_cmp(a._value, b._value);
    }

    friend Fraction operator+(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_add(result._value, a._value, b._value);
        return result;
    }

    friend Fraction operator-(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_sub(result._value, a._value, b._value);
        return result;
    }

    friend Fraction operator*(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_mul(result._value, a._value, b._value);
        return result;
    }

    friend Fraction operator/(Fraction const& a, Fraction const& b) {
        if (b.sign() == 0) {
            throw std::domain_error("FloatAP zero rational denominator");
        }
        Fraction result;
        mpq_div(result._value, a._value, b._value);
        return result;
    }
};

struct Interval {
    Fraction lower;
    Fraction upper;
};

Fraction absolute(Fraction const& x) {
    return x.sign() < 0 ? Fraction{} - x : x;
}

Fraction rational_argument(FloatAP const& x) {
    return Fraction(x);
}

namespace FloatAPExpLogDetail {

// Evaluate exp on |x| <= 1/2 with a single downward-rounded
// fixed-point Taylor recurrence. Let Q=2^B, a=floor(Q*|x|),
// t_0=Q, t_n=floor(t_(n-1)*a/(Q*n)), and T_n=Q*|x|^n/n!.
// As |x| <= 1/2, T_n - t_n < 6 for every n >= 1:
// the recurrence error satisfies e_n < e_(n-1)/(2*n)+2/n+1.
// After summing through n, the omitted exponential tail is at
// most T_n <= t_n+6, so the complete result is enclosed by
// [sum(t_i), sum(t_i)+t_n+6*(n+1)] / Q. This avoids a second
// rounded recurrence without relying on unproved negligible errors.
// For negative x take the reciprocal, reversing endpoints.
Interval exp_series(Fraction const& x, unsigned long work) {
    if (x.sign() == 0) {
        return {Fraction(1L), Fraction(1L)};
    }

    const Fraction magnitude = absolute(x);
    const unsigned long max_terms = work / 2ul + 64ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(max_terms));
    if (work > std::numeric_limits<unsigned long>::max() - guard) {
        throw std::overflow_error("FloatAP exponential work precision overflow");
    }
    const unsigned long fixed_bits = work + guard;
    if (static_cast<std::uintmax_t>(fixed_bits) >
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max()) ||
        static_cast<std::uintmax_t>(fixed_bits) >
        static_cast<std::uintmax_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("FloatAP exponential fixed-point shift overflow");
    }

    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const auto fixed_shift=static_cast<std::size_t>(fixed_bits);
    const std::size_t whole=fixed_shift/64u;
    const unsigned partial=static_cast<unsigned>(fixed_shift%64u);
    const Integer native_numerator=integer_from_gmp(mpq_numref(magnitude.raw()));
    const Integer native_denominator=integer_from_gmp(mpq_denref(magnitude.raw()));
    const BigUnsigned argument=quot(
        shift_left(native_numerator.magnitude(),fixed_shift),
        native_denominator.magnitude());

    // The five buffers retain their capacities across all Taylor terms.
    // In contrast to chained owning arithmetic, no temporary BigUnsigned
    // or GMP integer is constructed for multiplication, shifting, division,
    // or accumulation in the recurrence.
    Words term(whole+1u,0u);
    term[whole]=Limb{1}<<partial;
    Words sum_lower=term;
    Words product,shifted,next_term,next_sum;
    const BigUnsigned threshold=shift_left(BigUnsigned(1u),
                                           static_cast<std::size_t>(guard));
    const auto trim=[](Words& value) {
        while (!value.empty() && value.back()==0u) value.pop_back();
    };
    const auto below_threshold=[&](Words const& value) {
        const auto bound=threshold.limbs();
        if (value.size()!=bound.size()) return value.size()<bound.size();
        for (std::size_t i=value.size();i!=0u;--i) {
            if (value[i-1u]!=bound[i-1u]) return value[i-1u]<bound[i-1u];
        }
        return true;
    };

    for (unsigned long n = 1ul; n <= max_terms; ++n) {
        // t_n=floor(floor(t_(n-1)*a/2^B)/n). Reuse the same
        // multiplication, shifting and short-division buffers.
        if (term.empty() || argument.is_zero()) {
            product.clear();
        } else {
            if (term.size()>std::numeric_limits<std::size_t>::max()-argument.size()) {
                throw std::overflow_error("FloatAP exponential limb count overflow");
            }
            product.resize(term.size()+argument.size());
            BigUnsigned::mul_into(term,argument.limbs(),product);
            trim(product);
        }

        BigUnsigned::LimbView right_shifted;
        if (product.size()>whole) {
            const std::size_t length=product.size()-whole;
            if (partial==0u) {
                right_shifted=BigUnsigned::LimbView(product.data()+whole,length);
            } else {
                shifted.resize(length);
                for (std::size_t i=0u;i<length;++i) {
                    const Limb low=product[whole+i]>>partial;
                    const Limb high=i+1u<length
                        ? product[whole+i+1u]<<(64u-partial) : Limb{0};
                    shifted[i]=low|high;
                }
                trim(shifted);
                right_shifted=shifted;
            }
        }
        next_term.resize(right_shifted.size());
        (void)BigUnsigned::divrem_limb_into(
            right_shifted,static_cast<Limb>(n),next_term);
        trim(next_term);
        term.swap(next_term);

        const std::size_t count=std::max(sum_lower.size(),term.size());
        next_sum.resize(count+1u);
        const Limb carry=BigUnsigned::add_into(
            sum_lower,term,BigUnsigned::MutableLimbView(next_sum.data(),count));
        next_sum[count]=carry;
        trim(next_sum);
        sum_lower.swap(next_sum);

        if (below_threshold(term) || n == max_terms) {
            // The same proven enclosure: sum + term + 6*(n+1).
            const BigUnsigned lower(BigUnsigned::LimbView(
                sum_lower.data(),sum_lower.size()));
            const BigUnsigned last_term(BigUnsigned::LimbView(
                term.data(),term.size()));
            const BigUnsigned tail=mul(
                add(BigUnsigned(static_cast<Limb>(n)),BigUnsigned(1u)),
                BigUnsigned(6u));
            const BigUnsigned upper=add(add(lower,last_term),tail);
            const Interval positive{Fraction(lower,fixed_bits),
                                    Fraction(upper,fixed_bits)};
            if (x.sign() > 0) {
                return positive;
            }
            return {Fraction(1L) / positive.upper,
                    Fraction(1L) / positive.lower};
        }
    }
    throw std::logic_error("FloatAP exponential series did not terminate");
}

Interval exp_interval(Fraction x, unsigned long work) {
    // Preserve argument reduction, Taylor enclosure, and the exact
    // outward-directed repeated-squaring reconstruction.
    const unsigned long reduction_bits = std::min(
        16ul, std::max(6ul, work / 32ul));
    const Fraction reduction_limit =
        Fraction(1L) / Fraction(static_cast<long>(1ul << reduction_bits));

    unsigned halvings = 0;
    while (compare(absolute(x), reduction_limit) > 0) {
        if (++halvings > 2048u) {
            throw std::overflow_error("FloatAP exp argument range");
        }
        x.halve();
    }
    if (work > std::numeric_limits<unsigned long>::max() -
                    static_cast<unsigned long>(halvings) - 64ul) {
        throw std::overflow_error("FloatAP exponential reconstruction precision overflow");
    }
    const unsigned long reduced_work = work + halvings + 8ul;
    Interval interval = exp_series(x, reduced_work);
    if (halvings == 0u) {
        return interval;
    }

    const unsigned long fixed_bits = work + halvings + 64ul;
    if (static_cast<std::uintmax_t>(fixed_bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max()) ||
        static_cast<std::uintmax_t>(fixed_bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("FloatAP exponential reconstruction shift overflow");
    }
    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const std::size_t shift=static_cast<std::size_t>(fixed_bits);
    const std::size_t whole=shift/64u;
    const unsigned partial=static_cast<unsigned>(shift%64u);

    // Convert only the two initial certified rational bounds. The
    // lower is rounded down and the upper rounded up at the same
    // dyadic scale. Avoid canonical Rational/GCD round trips.
    const auto to_scaled=[&](Fraction const& bound,bool round_upward) {
        const Integer numerator=integer_from_gmp(mpq_numref(bound.raw()));
        const Integer denominator=integer_from_gmp(mpq_denref(bound.raw()));
        auto [value,remainder]=divrem(
            shift_left(numerator.magnitude(),shift),
            denominator.magnitude());
        if (round_upward && !remainder.is_zero()) {
            value=add(value,BigUnsigned(1u));
        }
        const auto digits=value.limbs();
        return Words(digits.begin(),digits.end());
    };
    Words lower=to_scaled(interval.lower,false);
    Words upper=to_scaled(interval.upper,true);
    Words product,next_lower,next_upper;
    const auto trim=[](Words& digits) {
        while (!digits.empty() && digits.back()==0u) digits.pop_back();
    };
    const auto increment=[](Words& digits) {
        for (std::size_t i=0u;i<digits.size();++i) {
            if (++digits[i]!=0u) return;
        }
        digits.push_back(Limb{1});
    };

    // floor(y*y/2^B) or ceil(y*y/2^B) as required. The discarded
    // full limbs and low partial-limb bits determine whether ceil
    // increments. The same product and output scratch is reused.
    const auto square_outward=[&](Words const& value,Words& result,bool round_upward) {
        if (value.empty()) {
            result.clear();
            return;
        }
        if (value.size()>std::numeric_limits<std::size_t>::max()/2u) {
            throw std::overflow_error("FloatAP exponential square size overflow");
        }
        product.resize(2u*value.size());
        BigUnsigned::mul_into(value,value,product);
        trim(product);

        bool discarded=false;
        if (round_upward) {
            for (std::size_t j=0u;j<std::min(whole,product.size());++j) {
                discarded = discarded || product[j]!=0u;
            }
            if (partial!=0u && product.size()>whole) {
                discarded = discarded ||
                    (product[whole] & ((Limb{1}<<partial)-1u))!=0u;
            }
        }

        if (product.size()<=whole) {
            result.clear();
        } else {
            const std::size_t count=product.size()-whole;
            result.resize(count);
            for (std::size_t j=0u;j<count;++j) {
                const Limb low=product[whole+j]>>partial;
                const Limb high=partial!=0u && j+1u<count
                    ? product[whole+j+1u]<<(64u-partial) : Limb{0};
                result[j]=low|high;
            }
            trim(result);
        }
        if (round_upward && discarded) increment(result);
    };

    // Square only positive endpoint magnitudes, preserving the
    // certified enclosure through monotonicity of squaring.
    for (unsigned i=0u;i<halvings;++i) {
        square_outward(lower,next_lower,false);
        lower.swap(next_lower);
        square_outward(upper,next_upper,true);
        upper.swap(next_upper);
    }
    return {
        Fraction(BigUnsigned(BigUnsigned::LimbView(
                     lower.data(),lower.size())),fixed_bits),
        Fraction(BigUnsigned(BigUnsigned::LimbView(
                     upper.data(),upper.size())),fixed_bits)
    };
}

// For |z| <= 1/3, log((1+z)/(1-z)) = 2 sum z^(2n+1)/(2n+1).
// Instead of growing exact rational denominators at every term, evaluate
// the series in fixed-point arithmetic with outward-rounded GMP integers.
// The first N terms are enclosed term by term. Since |z| <= 1/3,
// the remaining tail is at most (9/4)*3^(-2N-1), strictly below
// 2^(-work) for N = floor(work/3)+4.
Interval log_series(Fraction const& z, unsigned long work) {
    if (z.sign() == 0) {
        return {Fraction{}, Fraction{}};
    }

    const unsigned long terms = work / 3ul + 4ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(terms));
    if (work > std::numeric_limits<unsigned long>::max() - guard) {
        throw std::overflow_error("FloatAP logarithm work precision overflow");
    }
    const unsigned long fixed_bits = work + guard;
    if (static_cast<std::uintmax_t>(fixed_bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max()) ||
        static_cast<std::uintmax_t>(fixed_bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("FloatAP logarithm fixed-point shift overflow");
    }

    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const std::size_t shift=static_cast<std::size_t>(fixed_bits);
    const std::size_t whole=shift/64u;
    const unsigned partial=static_cast<unsigned>(shift%64u);
    const Fraction magnitude=absolute(z);
    const Integer num=integer_from_gmp(mpq_numref(magnitude.raw()));
    const Integer den=integer_from_gmp(mpq_denref(magnitude.raw()));
    auto [z_lower,remainder]=divrem(
        shift_left(num.magnitude(),shift),den.magnitude());
    const BigUnsigned z_upper=remainder.is_zero()
        ? z_lower : add(z_lower,BigUnsigned(1u));

    // Match GMP's exact floor/ceil endpoints. Only the initial fraction
    // import and final Fraction construction use the legacy GMP carrier.
    // Every fixed-point product, power, integer division and sum below
    // operates on reusable native caller-owned limb buffers.
    Words power_lo(z_lower.limbs().begin(),z_lower.limbs().end());
    Words power_hi(z_upper.limbs().begin(),z_upper.limbs().end());
    Words square_lo,square_hi,sum_lo,sum_hi;
    Words product,next_power,quotient,tail;
    const BigUnsigned threshold=shift_left(BigUnsigned(1u),
                                            static_cast<std::size_t>(guard));
    const auto trim=[](Words& v) {
        while (!v.empty() && v.back()==0u) v.pop_back();
    };
    const auto increment=[](Words& v) {
        std::size_t index=0u;
        while (index<v.size()) {
            if (++v[index]!=0u) return;
            ++index;
        }
        v.push_back(Limb{1});
    };

    // Given a,b >= 0, compute floor(a*b/2^B), or ceil when outward=true.
    // The full product is borrowed from one scratch buffer: discarded
    // bits, including the partial limb, determine the upward correction.
    const auto scaled_product=[&](Words const& left,Words const& right,
                                   bool outward,Words& output) {
        if (left.empty() || right.empty()) {
            output.clear();
            return;
        }
        if (left.size()>std::numeric_limits<std::size_t>::max()-right.size()) {
            throw std::overflow_error("FloatAP logarithm product size overflow");
        }
        product.resize(left.size()+right.size());
        BigUnsigned::mul_into(left,right,product);
        trim(product);

        bool discarded=false;
        if (outward) {
            for (std::size_t i=0u;i<std::min(whole,product.size());++i) {
                discarded = discarded || product[i]!=0u;
            }
            if (partial!=0u && product.size()>whole) {
                discarded = discarded ||
                    (product[whole] & ((Limb{1}<<partial)-1u))!=0u;
            }
        }

        if (product.size()<=whole) {
            output.clear();
        } else {
            const std::size_t count=product.size()-whole;
            output.resize(count);
            for (std::size_t i=0u;i<count;++i) {
                const Limb low=product[whole+i]>>partial;
                const Limb high=partial!=0u && i+1u<count
                    ? product[whole+i+1u]<<(64u-partial) : Limb{0};
                output[i]=low|high;
            }
            trim(output);
        }
        if (outward && discarded) increment(output);
    };

    // In-place nonnegative accumulator; carry propagation is bounded
    // by the current active sum and never creates a BigUnsigned.
    const auto accumulate=[&](Words& sum,Words const& addend) {
        if (sum.size()<addend.size()) sum.resize(addend.size(),0u);
        Limb carry=0u;
        for (std::size_t i=0u;i<addend.size();++i) {
            const Limb old=sum[i];
            const Limb value=old+addend[i];
            const Limb result=value+carry;
            carry=static_cast<Limb>(value<old || result<value);
            sum[i]=result;
        }
        for (std::size_t i=addend.size();carry!=0u && i<sum.size();++i) {
            ++sum[i];
            carry=static_cast<Limb>(sum[i]==0u);
        }
        if (carry!=0u) sum.push_back(carry);
    };
    const auto within_threshold=[&](Words const& value) {
        const auto bound=threshold.limbs();
        if (value.size()!=bound.size()) return value.size()<bound.size();
        for (std::size_t i=value.size();i!=0u;--i) {
            if (value[i-1u]!=bound[i-1u]) return value[i-1u]<bound[i-1u];
        }
        return true;
    };

    scaled_product(power_lo,power_lo,false,square_lo);
    scaled_product(power_hi,power_hi,true,square_hi);

    // The original proof and work-dependent bound are unchanged.
    // The Taylor tail after each term is <= 3 * upper odd power.
    for (unsigned long n=0u;n<terms;++n) {
        if (n>(std::numeric_limits<unsigned long>::max()-1ul)/2ul) {
            throw std::overflow_error("FloatAP logarithm series index overflow");
        }
        const Limb divisor=static_cast<Limb>(2ul*n+1ul);
        quotient.resize(power_lo.size());
        (void)BigUnsigned::divrem_limb_into(power_lo,divisor,quotient);
        trim(quotient);
        accumulate(sum_lo,quotient);

        quotient.resize(power_hi.size());
        const Limb upper_remainder=
            BigUnsigned::divrem_limb_into(power_hi,divisor,quotient);
        trim(quotient);
        if (upper_remainder!=0u) increment(quotient);
        accumulate(sum_hi,quotient);

        scaled_product(power_lo,square_lo,false,next_power);
        power_lo.swap(next_power);
        scaled_product(power_hi,square_hi,true,next_power);
        power_hi.swap(next_power);

        if (power_hi.empty()) {
            tail.clear();
        } else {
            const Limb three=3u;
            tail.resize(power_hi.size()+1u);
            BigUnsigned::mul_into(power_hi,BigUnsigned::LimbView(&three,1u),
                                  tail);
            trim(tail);
        }
        if (within_threshold(tail) || n+1ul==terms) {
            break;
        }
    }

    const BigUnsigned lower=shift_left(
        BigUnsigned(BigUnsigned::LimbView(sum_lo.data(),sum_lo.size())),1u);
    const BigUnsigned upper=add(
        shift_left(BigUnsigned(BigUnsigned::LimbView(
            sum_hi.data(),sum_hi.size())),1u),
        BigUnsigned(BigUnsigned::LimbView(tail.data(),tail.size())));
    const Interval positive{Fraction(lower,fixed_bits),
                            Fraction(upper,fixed_bits)};
    return z.sign()<0
        ? Interval{Fraction{}-positive.upper,Fraction{}-positive.lower}
        : positive;
}

Interval log2_interval(unsigned long work) {
    // An enclosure computed at higher working precision remains valid
    // for all lower-precision requests. This cache contains an interval,
    // not a rounded scalar, and is independent of the rounding mode.
    struct Cache {
        unsigned long work = 0;
        std::optional<Interval> value;
    };
    static thread_local Cache cache;
    if (cache.value && cache.work >= work) {
        return *cache.value;
    }

    const Interval result = log_series(Fraction(1L) / Fraction(3L), work);
    cache.value = result;
    cache.work = work;
    return result;
}

Interval log_interval(Fraction x, unsigned long work) {
    long power = 0;
    const Fraction lower_limit = Fraction(3L) / Fraction(4L);
    const Fraction upper_limit = Fraction(3L) / Fraction(2L);
    // Reduce to [3/4,3/2], where |(x-1)/(x+1)| <= 1/5.
    // This also moves numbers near a power of two close to 1,
    // permitting the fixed-point series to terminate very early.
    while (compare(x, upper_limit) > 0) {
        if (++power > 2048) {
            throw std::overflow_error("FloatAP log argument range");
        }
        x.halve();
    }
    while (compare(x, lower_limit) < 0) {
        if (--power < -2048) {
            throw std::overflow_error("FloatAP log argument range");
        }
        x = x * Fraction(2L);
    }
    const Fraction z = (x - Fraction(1L)) / (x + Fraction(1L));
    Interval result = log_series(z, work);
    if (power != 0) {
        const Interval ln2 = log2_interval(work);
        const Fraction count(power < 0 ? -power : power);
        if (power > 0) {
            result = {result.lower + count * ln2.lower,
                      result.upper + count * ln2.upper};
        } else {
            result = {result.lower - count * ln2.upper,
                      result.upper - count * ln2.lower};
        }
    }
    return result;
}

// Certified half-centered exponential shortcut. The sign-dependent
// constants are cached as rigorous outward DYADIC enclosures, not
// approximations rounded in a caller-selected mode. The general
// exp_interval remains the fallback everywhere outside the window.
bool exp_inside_half_window(Fraction const& delta, unsigned long work) {
    // |delta| <= 2^(-ceil(work/2)), with work>=4 guaranteeing the
    // entire segment around +/-1/2 lies inside [-3/4,3/4].
    if (work < 4ul) return false;
    const unsigned long cutoff_bits = work / 2ul + work % 2ul;
    if (static_cast<std::uintmax_t>(cutoff_bits) >
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        return false;
    }
    mpz_t scaled_absolute_offset;
    mpz_init(scaled_absolute_offset);
    mpz_abs(scaled_absolute_offset, mpq_numref(delta.raw()));
    mpz_mul_2exp(scaled_absolute_offset, scaled_absolute_offset,
                 static_cast<mp_bitcnt_t>(cutoff_bits));
    const bool inside_center_window =
        mpz_cmp(scaled_absolute_offset, mpq_denref(delta.raw())) <= 0;
    mpz_clear(scaled_absolute_offset);
    return inside_center_window;
}

// Convert a signed certified rational endpoint to an exact outward
// dyadic quotient. Positive floor/ceil use quotient/remainder directly;
// for negative endpoints the direction of magnitude rounding reverses.
// This is a temporary boundary import until Fraction itself is native.
Fraction outward_dyadic_endpoint(Fraction const& exact,
                                 unsigned long bits,
                                 bool round_toward_positive) {
    if (static_cast<std::uintmax_t>(bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max()) ||
        static_cast<std::uintmax_t>(bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("FloatAP dyadic endpoint shift overflow");
    }
    const Integer numerator=integer_from_gmp(mpq_numref(exact.raw()));
    const Integer denominator=integer_from_gmp(mpq_denref(exact.raw()));
    auto [magnitude,remainder]=divrem(
        shift_left(numerator.magnitude(),static_cast<std::size_t>(bits)),
        denominator.magnitude());
    const bool negative=numerator.negative();
    const bool round_magnitude_up=negative
        ? !round_toward_positive : round_toward_positive;
    if (round_magnitude_up && !remainder.is_zero()) {
        magnitude=add(magnitude,BigUnsigned(1u));
    }
    return Fraction(magnitude,bits,negative);
}

Interval exp_outward_dyadic(Interval const& bounds,
                            unsigned long denominator_bits) {
    return {
        outward_dyadic_endpoint(bounds.lower,denominator_bits,false),
        outward_dyadic_endpoint(bounds.upper,denominator_bits,true)
    };
}

struct ExponentialHalfCenters {
    Interval positive;
    Interval negative;
};

ExponentialHalfCenters const& exponential_half_centers(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<ExponentialHalfCenters> enclosure;
    };
    static thread_local Cache cache;
    if (!cache.enclosure || cache.work < work) {
        if (work > std::numeric_limits<unsigned long>::max() - 40ul) {
            throw std::overflow_error("FloatAP exponential half precision overflow");
        }
        const Fraction half = Fraction(1L) / Fraction(2L);
        const Interval positive = exp_interval(half, work + 40ul);
        // Both bounds are strictly positive, so reciprocity reverses
        // their order and yields an exact interval for exp(-1/2).
        const Interval negative{
            Fraction(1L) / positive.upper,
            Fraction(1L) / positive.lower
        };
        ExponentialHalfCenters fresh{
            exp_outward_dyadic(positive, work + 24ul),
            exp_outward_dyadic(negative, work + 24ul)
        };
        cache.enclosure = std::move(fresh);
        cache.work = work;
    }
    return *cache.enclosure;
}

Interval exp_half_centered_or_general(Fraction const& input,
                                     unsigned long work) {
    const Fraction half = Fraction(1L) / Fraction(2L);
    const Fraction center = input.sign() < 0 ? Fraction{} - half : half;
    const Fraction delta = input - center;
    if (!exp_inside_half_window(delta, work)) {
        return exp_interval(input, work);
    }

    // exp(center+delta)=exp(center)*(1+delta)+R.
    // On the segment [-3/4,3/4], |exp''(t)| < e < 3,
    // hence Taylor gives |R| < 3*delta^2/2 <= 2*delta^2.
    // 1+delta >= 3/4 is positive, so endpoint multiplication
    // preserves order. All operations below use exact rationals.
    const auto& cache = exponential_half_centers(work);
    const Interval& base = input.sign() < 0
        ? cache.negative : cache.positive;
    const Fraction multiplier = Fraction(1L) + delta;
    const Fraction remainder = Fraction(2L) * delta * delta;
    return {base.lower * multiplier - remainder,
            base.upper * multiplier + remainder};
}

// Certified log close to 1/2, analogous to the centered exp route.
// For x=1/2+delta with |delta|<=2^-ceil(work/2), work>=4,
// x lies in [1/4,3/4] and |log''(x)|=1/x^2<=16.
// Consequently log(x)=-log(2)+2*delta+R, |R|<=8*delta^2.
// Cache an outward interval of -log(2); the exact rational slope
// costs no extra certified constant. Other arguments use log_interval.
Interval log_half_centered_or_general(Fraction const& input,
                                     unsigned long work) {
    const Fraction delta = input - Fraction(1L) / Fraction(2L);
    if (!exp_inside_half_window(delta, work)) {
        return log_interval(input, work);
    }

    struct Cache {
        unsigned long work = 0;
        std::optional<Interval> negative_log_two;
    };
    static thread_local Cache cache;
    if (!cache.negative_log_two || cache.work < work) {
        if (work > std::numeric_limits<unsigned long>::max() - 40ul) {
            throw std::overflow_error("FloatAP centered log precision overflow");
        }
        const Interval certified_log_two = log2_interval(work + 40ul);
        const Interval negative_center{
            Fraction{} - certified_log_two.upper,
            Fraction{} - certified_log_two.lower
        };
        // The helper rounds interval endpoints outward by exact GMP
        // floor/ceiling, including for negative values.
        cache.negative_log_two = exp_outward_dyadic(
            negative_center, work + 24ul);
        cache.work = work;
    }

    const Fraction linear_term = Fraction(2L) * delta;
    const Fraction maximum_remainder = Fraction(8L) * delta * delta;
    return {cache.negative_log_two->lower + linear_term - maximum_remainder,
            cache.negative_log_two->upper + linear_term + maximum_remainder};
}

template<class Evaluate>
FloatAP certified(Evaluate evaluate, FloatAP::RoundingModeType mode,
                  FloatAP::PrecisionType precision) {
    const unsigned long bits = static_cast<unsigned long>(precision.bits());
    unsigned long extra = 32;
    for (;;) {
        if (extra > std::numeric_limits<unsigned long>::max() - bits - 16) {
            throw std::overflow_error("FloatAP certified precision overflow");
        }
        const Interval interval = evaluate(bits + extra);
        const FloatAP lower = interval.lower.rounded(bits,mode);
        const FloatAP upper = interval.upper.rounded(bits,mode);
        if (lower == upper) {
            return lower;
        }
        if (extra > std::numeric_limits<unsigned long>::max() / 2) {
            throw std::overflow_error("FloatAP certified refinement overflow");
        }
        extra *= 2;
    }
}

} // namespace FloatAPExpLogDetail

namespace FloatAPTrigonometricDetail {

// Arithmetic on these fractions is exact. The approximation error is
// represented only by the endpoints of Interval, never by floating point.

Interval point(Fraction const& value) {
    return {value, value};
}

Interval add(Interval const& a, Interval const& b) {
    return {a.lower + b.lower, a.upper + b.upper};
}

Interval subtract(Interval const& a, Interval const& b) {
    return {a.lower - b.upper, a.upper - b.lower};
}

Interval negative(Interval const& value) {
    return {Fraction{} - value.upper, Fraction{} - value.lower};
}

Interval multiply(Interval const& a, Interval const& b) {
    const Fraction candidates[4] = {
        a.lower * b.lower, a.lower * b.upper,
        a.upper * b.lower, a.upper * b.upper
    };
    Fraction lower(candidates[0]);
    Fraction upper(candidates[0]);
    for (unsigned i = 1; i < 4; ++i) {
        if (compare(candidates[i], lower) < 0) {
            lower = candidates[i];
        }
        if (compare(candidates[i], upper) > 0) {
            upper = candidates[i];
        }
    }
    return {std::move(lower), std::move(upper)};
}

Interval divide(Interval const& numerator, Interval const& denominator) {
    if (denominator.lower.sign() <= 0 && denominator.upper.sign() >= 0) {
        throw std::domain_error("FloatAP trigonometric interval crosses a pole");
    }

    // Monotonicity and endpoint signs identify the extrema. Only two
    // exact rational divisions are needed instead of two reciprocals
    // and four interval cross-products.
    if (denominator.lower.sign() > 0) {
        if (numerator.lower.sign() >= 0) {
            return {numerator.lower / denominator.upper,
                    numerator.upper / denominator.lower};
        }
        if (numerator.upper.sign() <= 0) {
            return {numerator.lower / denominator.lower,
                    numerator.upper / denominator.upper};
        }
        return {numerator.lower / denominator.lower,
                numerator.upper / denominator.lower};
    }

    if (numerator.lower.sign() >= 0) {
        return {numerator.upper / denominator.upper,
                numerator.lower / denominator.lower};
    }
    if (numerator.upper.sign() <= 0) {
        return {numerator.upper / denominator.lower,
                numerator.lower / denominator.upper};
    }
    return {numerator.upper / denominator.upper,
            numerator.lower / denominator.upper};
}

// For |x| <= 1, sine and cosine have alternating, decreasing Taylor
// terms. Instead of two directed fixed-point recurrences, calculate one
// downward recurrence and add a proven bound for every discarded bit.
// The cosine needed by tangent is still derived from sine by a directed
// integer square root, so only a single Taylor series is evaluated.
enum class TrigOutputs { sine, cosine, both };

std::pair<Interval, Interval> small_sine_cosine(Fraction const& x,
                                               unsigned long work,
                                               TrigOutputs outputs) {
    if (x.sign() == 0) {
        return {point(Fraction{}), point(Fraction(1L))};
    }

    const unsigned long max_terms = work / 2ul + 48ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(max_terms));
    if (work > std::numeric_limits<unsigned long>::max() - guard) {
        throw std::overflow_error("FloatAP sine/cosine work precision overflow");
    }
    const unsigned long fixed_bits = work + guard;
    if (static_cast<std::uintmax_t>(fixed_bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max()) ||
        static_cast<std::uintmax_t>(fixed_bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("FloatAP sine/cosine fixed-point shift overflow");
    }

    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const std::size_t shift=static_cast<std::size_t>(fixed_bits);
    const std::size_t whole=shift/64u;
    const unsigned partial=static_cast<unsigned>(shift%64u);
    const Fraction magnitude=absolute(x);
    const Integer numerator=integer_from_gmp(mpq_numref(magnitude.raw()));
    const Integer denominator=integer_from_gmp(mpq_denref(magnitude.raw()));
    const BigUnsigned scale=shift_left(BigUnsigned(1u),shift);
    const BigUnsigned argument=quot(
        shift_left(numerator.magnitude(),shift),denominator.magnitude());

    Words square,term,sum,product,next_term,quotient;
    const BigUnsigned threshold=shift_left(BigUnsigned(1u),
                                            static_cast<std::size_t>(guard));
    const auto trim=[](Words& value) {
        while (!value.empty() && value.back()==0u) value.pop_back();
    };
    const auto within_threshold=[&](Words const& value) {
        const auto bound=threshold.limbs();
        if (value.size()!=bound.size()) return value.size()<bound.size();
        for (std::size_t i=value.size();i!=0u;--i) {
            if (value[i-1u]!=bound[i-1u]) return value[i-1u]<bound[i-1u];
        }
        return true;
    };

    // The series is alternating with decreasing positive terms for |x|<=1.
    // Its partial sum is therefore nonnegative, so update one limb vector
    // in place rather than allocating signed Integer/Rational temporaries.
    const auto accumulate=[&](Words const& addend,bool subtract) {
        if (!subtract) {
            if (sum.size()<addend.size()) sum.resize(addend.size(),0u);
            Limb carry=0u;
            for (std::size_t i=0u;i<addend.size();++i) {
                const Limb old=sum[i];
                const Limb value=old+addend[i];
                const Limb result=value+carry;
                carry=static_cast<Limb>(value<old || result<value);
                sum[i]=result;
            }
            for (std::size_t i=addend.size();carry!=0u && i<sum.size();++i) {
                ++sum[i];
                carry=static_cast<Limb>(sum[i]==0u);
            }
            if (carry!=0u) sum.push_back(carry);
        } else {
            if (addend.size()>sum.size()) {
                throw std::logic_error("FloatAP sine/cosine series sum underflow");
            }
            Limb borrow=0u;
            for (std::size_t i=0u;i<sum.size();++i) {
                const Limb right=i<addend.size()?addend[i]:Limb{0};
                const Limb old=sum[i];
                const Limb difference=old-right;
                const Limb result=difference-borrow;
                borrow=static_cast<Limb>(old<right || difference<borrow);
                sum[i]=result;
            }
            if (borrow!=0u) {
                throw std::logic_error("FloatAP sine/cosine series sum underflow");
            }
            trim(sum);
        }
    };

    // floor(a*b/2^B), with disjoint scratch and output limb storage.
    // Reuse product and result capacity across Taylor terms.
    const auto scaled_product=[&](Words const& left,Words const& right,
                                   Words& output) {
        if (left.empty() || right.empty()) {
            output.clear();
            return;
        }
        if (left.size()>std::numeric_limits<std::size_t>::max()-right.size()) {
            throw std::overflow_error("FloatAP sine/cosine product size overflow");
        }
        product.resize(left.size()+right.size());
        BigUnsigned::mul_into(left,right,product);
        trim(product);
        if (product.size()<=whole) {
            output.clear();
            return;
        }
        const std::size_t count=product.size()-whole;
        output.resize(count);
        for (std::size_t i=0u;i<count;++i) {
            const Limb low=product[whole+i]>>partial;
            const Limb high=partial!=0u && i+1u<count
                ? product[whole+i+1u]<<(64u-partial) : Limb{0};
            output[i]=low|high;
        }
        trim(output);
    };
    Words argument_words(argument.limbs().begin(),argument.limbs().end());
    scaled_product(argument_words,argument_words,square);

    // When both functions are requested, derive cosine from the sine
    // enclosure; retain its native endpoints to avoid GMP round-trips.
    BigUnsigned sine_lower_nonnegative,sine_upper_nonnegative;
    const auto series=[&](bool sine) -> Interval {
        const auto start=sine ? argument.limbs() : scale.limbs();
        term.assign(start.begin(),start.end());
        sum.clear();

        for (unsigned long n=0u;n<max_terms;++n) {
            accumulate(term,n%2ul!=0ul);

            if (n>(std::numeric_limits<unsigned long>::max()-3ul)/2ul) {
                throw std::overflow_error("FloatAP sine/cosine series index");
            }
            const unsigned long twice_n=2ul*n;
            const unsigned long a=sine ? twice_n+2ul : twice_n+1ul;
            const unsigned long b=sine ? twice_n+3ul : twice_n+2ul;
            if (a>std::numeric_limits<unsigned long>::max()/b) {
                throw std::overflow_error("FloatAP sine/cosine factorial overflow");
            }
            const Limb divisor=static_cast<Limb>(a*b);

            scaled_product(term,square,next_term);
            quotient.resize(next_term.size());
            (void)BigUnsigned::divrem_limb_into(next_term,divisor,quotient);
            trim(quotient);
            term.swap(quotient);

            if (within_threshold(term) || n+1ul==max_terms) {
                // Identical certificate: 8*(n+1) accumulated units,
                // followed by the first omitted term plus 8 units.
                const BigUnsigned result(BigUnsigned::LimbView(
                    sum.data(),sum.size()));
                const BigUnsigned error=mul(
                    BigUnsigned(static_cast<Limb>(n+1ul)),BigUnsigned(8u));
                const BigUnsigned tail=add(
                    BigUnsigned(BigUnsigned::LimbView(
                        term.data(),term.size())),BigUnsigned(8u));
                const BigUnsigned lower_error=n%2ul==0ul
                    ? add(error,tail) : error;
                const BigUnsigned upper_error=n%2ul==0ul
                    ? error : add(error,tail);
                const bool lower_negative=compare(result,lower_error)<0;
                const BigUnsigned lower=lower_negative
                    ? sub(lower_error,result) : sub(result,lower_error);
                const BigUnsigned upper=add(result,upper_error);
                if (sine && outputs==TrigOutputs::both) {
                    sine_lower_nonnegative=lower_negative
                        ? BigUnsigned{} : lower;
                    sine_upper_nonnegative=upper;
                }
                return {Fraction(lower,fixed_bits,lower_negative),
                        Fraction(upper,fixed_bits)};
            }
        }
        throw std::logic_error("FloatAP sine/cosine series did not terminate");
    };

    const bool want_sine=outputs!=TrigOutputs::cosine;
    Interval sine=want_sine ? series(true) : point(Fraction{});
    Interval cosine=point(Fraction{});

    if (outputs==TrigOutputs::cosine) {
        cosine=series(false);
    } else if (outputs==TrigOutputs::both) {
        // The sine interval must be intersected with [0,1].
        // For |x|<=1, cosine is nonnegative and sin(|x|)<=sin(1).
        const BigUnsigned& lower=sine_lower_nonnegative;
        const BigUnsigned upper=compare(sine_upper_nonnegative,scale)>0
            ? scale : sine_upper_nonnegative;
        if (compare(lower,upper)>0) {
            throw std::logic_error("FloatAP invalid certified sine bounds");
        }

        // When U<=8Q/9, sqrt(Q^2-s^2) is 2-Lipschitz on
        // the whole [L,U]. One integer square root is sufficient;
        // otherwise compute exact outward roots for both endpoints.
        const BigUnsigned squared_scale=mul(scale,scale);
        BigUnsigned cosine_lower,cosine_upper;
        if (compare(mul(upper,BigUnsigned(9u)),
                    mul(scale,BigUnsigned(8u)))<=0) {
            cosine_lower=isqrt(sub(squared_scale,mul(upper,upper)));
            const BigUnsigned width=sub(upper,lower);
            cosine_upper=add(cosine_lower,
                add(shift_left(width,1u),BigUnsigned(1u)));
        } else {
            auto [ceil_root,root_remainder]=sqrtrem(
                sub(squared_scale,mul(lower,lower)));
            if (!root_remainder.is_zero()) {
                ceil_root=add(ceil_root,BigUnsigned(1u));
            }
            cosine_lower=isqrt(sub(squared_scale,mul(upper,upper)));
            cosine_upper=std::move(ceil_root);
        }
        cosine={Fraction(cosine_lower,fixed_bits),
                Fraction(cosine_upper,fixed_bits)};
    }

    if (want_sine && x.sign()<0) {
        sine=negative(sine);
    }
    return {std::move(sine),std::move(cosine)};
}

// Multiply an interval by a signed exact scalar without evaluating four
// products. This is used in the centered trigonometric enclosure.
Interval multiply_by_fraction(Interval const& input, Fraction const& scalar) {
    if (scalar.sign() >= 0) {
        return {input.lower * scalar, input.upper * scalar};
    }
    return {input.upper * scalar, input.lower * scalar};
}

// Cache certified enclosures of sin(1/2) and cos(1/2), rather than
// rounded scalars. The same (or tighter) intervals are valid for every
// later request of no greater working precision on this thread.
std::pair<Interval, Interval> const& half_trigonometric_intervals(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<std::pair<Interval, Interval>> bounds;
    };
    static thread_local Cache cache;
    if (!cache.bounds || cache.work < work) {
        cache.bounds = small_sine_cosine(
            Fraction(1L) / Fraction(2L), work, TrigOutputs::both);
        cache.work = work;
    }
    return *cache.bounds;
}

// The cached half-angle sine/cosine endpoints are certified rationals.
// Convert each one OUTWARD to a shared dyadic scale once per thread and
// work-precision upgrade, then perform the local Taylor enclosure using
// integer significands rather than generic rational arithmetic.
std::pair<Interval, Interval> centered_trigonometric_dyadic(
    Fraction const& argument, Fraction const& offset,
    unsigned long work, TrigOutputs output) {
    struct Cache {
        unsigned long work = 0;
        unsigned long bits = 0;
        mpz_t sine_lower, sine_upper, cosine_lower, cosine_upper;
        Cache() {
            mpz_inits(sine_lower, sine_upper, cosine_lower, cosine_upper, nullptr);
        }
        ~Cache() {
            mpz_clears(sine_lower, sine_upper, cosine_lower, cosine_upper, nullptr);
        }
        Cache(Cache const&) = delete;
        Cache& operator=(Cache const&) = delete;
    };
    static thread_local Cache fixed_center;

    if (fixed_center.work < work) {
        if (work > std::numeric_limits<unsigned long>::max() - 32ul) {
            throw std::overflow_error("FloatAP centered fixed-point precision overflow");
        }
        const unsigned long target_bits = work + 16ul;
        const auto& [sine, cosine] =
            half_trigonometric_intervals(work + 32ul);

        mpz_t numerator_shifted;
        mpz_init(numerator_shifted);
        struct Cleanup {
            mpz_ptr value;
            ~Cleanup() { mpz_clear(value); }
        } cleanup{numerator_shifted};
        const auto assign_outward = [&](mpz_ptr destination, Fraction const& value,
                                        bool round_toward_positive_infinity) {
            mpz_mul_2exp(numerator_shifted, mpq_numref(value.raw()),
                         static_cast<mp_bitcnt_t>(target_bits));
            if (round_toward_positive_infinity) {
                mpz_cdiv_q(destination, numerator_shifted, mpq_denref(value.raw()));
            } else {
                mpz_fdiv_q(destination, numerator_shifted, mpq_denref(value.raw()));
            }
        };
        assign_outward(fixed_center.sine_lower, sine.lower, false);
        assign_outward(fixed_center.sine_upper, sine.upper, true);
        assign_outward(fixed_center.cosine_lower, cosine.lower, false);
        assign_outward(fixed_center.cosine_upper, cosine.upper, true);
        fixed_center.bits = target_bits;
        fixed_center.work = work;
    }

    // The caller established a dyadic offset, with positive canonical
    // denominator 2^k, and |offset| <= 2^-ceil(work/2).
    const mp_bitcnt_t shift = mpz_scan1(mpq_denref(offset.raw()), 0u);
    mpz_srcptr delta_numerator = mpq_numref(offset.raw());

    mpz_t term_lower, term_upper, result_lower, result_upper, remainder;
    mpz_inits(term_lower, term_upper, result_lower, result_upper, remainder, nullptr);
    struct TermsCleanup {
        mpz_ptr a, b, c, d, e;
        ~TermsCleanup() { mpz_clears(a, b, c, d, e, nullptr); }
    } terms_cleanup{term_lower,term_upper,result_lower,result_upper,remainder};

    // Taylor: |R|<=offset^2/2<=2^(-work-1). At the cached scale
    // Q=2^bits this is exactly <=2^(bits-work-1) integer units.
    const unsigned long remaining_bits = fixed_center.bits - work - 1ul;
    mpz_set_ui(remainder, 1ul);
    mpz_mul_2exp(remainder, remainder,
                 static_cast<mp_bitcnt_t>(remaining_bits));

    const auto mul_outward = [&](mpz_srcptr coefficient_lower,
                                 mpz_srcptr coefficient_upper) {
        const bool nonnegative = mpz_sgn(delta_numerator) >= 0;
        mpz_mul(term_lower,
                nonnegative ? coefficient_lower : coefficient_upper,
                delta_numerator);
        mpz_fdiv_q_2exp(term_lower, term_lower, shift);
        mpz_mul(term_upper,
                nonnegative ? coefficient_upper : coefficient_lower,
                delta_numerator);
        mpz_cdiv_q_2exp(term_upper, term_upper, shift);
    };

    Interval sine_result = point(Fraction{});
    Interval cosine_result = point(Fraction{});
    if (output != TrigOutputs::cosine) {
        // sin(1/2+d) = S + C*d + R, |R|<=d^2/2.
        mul_outward(fixed_center.cosine_lower, fixed_center.cosine_upper);
        mpz_add(result_lower, fixed_center.sine_lower, term_lower);
        mpz_sub(result_lower, result_lower, remainder);
        mpz_add(result_upper, fixed_center.sine_upper, term_upper);
        mpz_add(result_upper, result_upper, remainder);
        sine_result = {Fraction(result_lower, fixed_center.bits),
                       Fraction(result_upper, fixed_center.bits)};
        if (argument.sign() < 0) sine_result = negative(sine_result);
    }
    if (output != TrigOutputs::sine) {
        // cos(1/2+d) = C - S*d + R, |R|<=d^2/2.
        mul_outward(fixed_center.sine_lower, fixed_center.sine_upper);
        mpz_sub(result_lower, fixed_center.cosine_lower, term_upper);
        mpz_sub(result_lower, result_lower, remainder);
        mpz_sub(result_upper, fixed_center.cosine_upper, term_lower);
        mpz_add(result_upper, result_upper, remainder);
        cosine_result = {Fraction(result_lower, fixed_center.bits),
                         Fraction(result_upper, fixed_center.bits)};
    }
    return {std::move(sine_result), std::move(cosine_result)};
}

// For finite delta, test |delta| <= 2^(-ceil(work/2)) using exact
// integer arithmetic. Then delta^2/2 <= 2^(-work-1) (work even)
// or even smaller (work odd), so a first-order Taylor enclosure
// has a certified remainder below one bit at the requested work.
bool sufficiently_near_half(Fraction const& delta, unsigned long work) {
    const unsigned long half_work = work / 2ul + work % 2ul;
    if (static_cast<std::uintmax_t>(half_work) >
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        return false;
    }
    mpz_t scaled;
    mpz_init(scaled);
    mpz_abs(scaled, mpq_numref(delta.raw()));
    mpz_mul_2exp(scaled, scaled, static_cast<mp_bitcnt_t>(half_work));
    const bool within_half_threshold =
        mpz_cmp(scaled, mpq_denref(delta.raw())) <= 0;
    mpz_clear(scaled);
    return within_half_threshold;
}

// Evaluate atan(x) for |x| <= 1/2 as an alternating Taylor series.
// All intermediates are GMP fixed-point integers: rounding each positive
// power and each term outward yields an enclosure of the partial sum.
// The alternating remainder has the sign of the first omitted term and
// magnitude bounded by that term. No GMP rational is allocated per term.
Interval small_atan(Fraction const& x, unsigned long work) {
    if (x.sign() == 0) {
        return point(Fraction{});
    }

    const unsigned long max_terms = work / 2ul + 4ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(max_terms));
    if (work > std::numeric_limits<unsigned long>::max() - guard) {
        throw std::overflow_error("FloatAP arctangent work precision overflow");
    }
    const unsigned long fixed_bits = work + guard;
    if (static_cast<std::uintmax_t>(fixed_bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max()) ||
        static_cast<std::uintmax_t>(fixed_bits) >
            static_cast<std::uintmax_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("FloatAP arctangent fixed-point shift overflow");
    }

    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const std::size_t shift=static_cast<std::size_t>(fixed_bits);
    const std::size_t whole=shift/64u;
    const unsigned partial=static_cast<unsigned>(shift%64u);
    const Fraction magnitude=absolute(x);
    const Integer numerator=integer_from_gmp(mpq_numref(magnitude.raw()));
    const Integer denominator=integer_from_gmp(mpq_denref(magnitude.raw()));
    const BigUnsigned argument=quot(
        shift_left(numerator.magnitude(),shift),denominator.magnitude());

    Words power(argument.limbs().begin(),argument.limbs().end());
    Words sum,square,term,product,next_power;
    const BigUnsigned limit=shift_left(BigUnsigned(1u),
                                       static_cast<std::size_t>(guard));
    const auto trim=[](Words& digits) {
        while (!digits.empty() && digits.back()==0u) digits.pop_back();
    };

    // All Taylor powers are nonnegative. Output and product retain their
    // capacity across terms; no owning BigUnsigned is made in the loop.
    const auto scaled_product=[&](Words const& left,Words const& right,
                                   Words& output) {
        if (left.empty() || right.empty()) {
            output.clear();
            return;
        }
        if (left.size()>std::numeric_limits<std::size_t>::max()-right.size()) {
            throw std::overflow_error("FloatAP arctangent product size overflow");
        }
        product.resize(left.size()+right.size());
        BigUnsigned::mul_into(left,right,product);
        trim(product);
        if (product.size()<=whole) {
            output.clear();
            return;
        }
        const std::size_t count=product.size()-whole;
        output.resize(count);
        for (std::size_t i=0u;i<count;++i) {
            const Limb low=product[whole+i]>>partial;
            const Limb high=partial!=0u && i+1u<count
                ? product[whole+i+1u]<<(64u-partial) : Limb{0};
            output[i]=low|high;
        }
        trim(output);
    };

    // Alternating partial sums are nonnegative: terms decrease for
    // |x|<=1/2, and the very first term is positive.
    const auto accumulate=[&](Words const& addend,bool subtract) {
        if (!subtract) {
            if (sum.size()<addend.size()) sum.resize(addend.size(),0u);
            Limb carry=0u;
            for (std::size_t i=0u;i<addend.size();++i) {
                const Limb old=sum[i];
                const Limb value=old+addend[i];
                const Limb result=value+carry;
                carry=static_cast<Limb>(value<old || result<value);
                sum[i]=result;
            }
            for (std::size_t i=addend.size();carry!=0u && i<sum.size();++i) {
                ++sum[i];
                carry=static_cast<Limb>(sum[i]==0u);
            }
            if (carry!=0u) sum.push_back(carry);
        } else {
            Limb borrow=0u;
            for (std::size_t i=0u;i<sum.size();++i) {
                const Limb part=i<addend.size()?addend[i]:Limb{0};
                const Limb old=sum[i];
                const Limb difference=old-part;
                const Limb result=difference-borrow;
                borrow=static_cast<Limb>(old<part || difference<borrow);
                sum[i]=result;
            }
            if (borrow!=0u || addend.size()>sum.size()) {
                throw std::logic_error("FloatAP arctangent series sum underflow");
            }
            trim(sum);
        }
    };
    const auto under_limit=[&](Words const& digits) {
        const auto bound=limit.limbs();
        if (digits.size()!=bound.size()) return digits.size()<bound.size();
        for (std::size_t i=digits.size();i!=0u;--i) {
            if (digits[i-1u]!=bound[i-1u]) return digits[i-1u]<bound[i-1u];
        }
        return true;
    };

    // Q=2^B, a=floor(Q*|x|), b=floor(a*a/Q), then
    // p_(n+1)=floor(p_n*b/Q). Keep the original 4-unit power
    // and 5-unit per-term error certificates.
    scaled_product(power,power,square);
    for (unsigned long n=0u;n<max_terms;++n) {
        if (n>(std::numeric_limits<unsigned long>::max()-1ul)/2ul) {
            throw std::overflow_error("FloatAP arctangent series index overflow");
        }
        const Limb divisor=static_cast<Limb>(2ul*n+1ul);
        term.resize(power.size());
        (void)BigUnsigned::divrem_limb_into(power,divisor,term);
        trim(term);
        accumulate(term,n%2ul!=0ul);

        scaled_product(power,square,next_power);
        power.swap(next_power);

        if (under_limit(power) || n+1ul==max_terms) {
            // The first omitted alternating term is <= power+4.
            // n+1 included terms introduce at most 5*(n+1) units.
            const BigUnsigned result(BigUnsigned::LimbView(
                sum.data(),sum.size()));
            const BigUnsigned error=mul(
                BigUnsigned(static_cast<Limb>(n+1ul)),BigUnsigned(5u));
            const BigUnsigned tail=add(
                BigUnsigned(BigUnsigned::LimbView(
                    power.data(),power.size())),BigUnsigned(4u));
            const BigUnsigned lower_error=n%2ul==0ul
                ? add(error,tail) : error;
            const BigUnsigned upper_error=n%2ul==0ul
                ? error : add(error,tail);

            const bool lower_negative=compare(result,lower_error)<0;
            const BigUnsigned lower=lower_negative
                ? sub(lower_error,result) : sub(result,lower_error);
            const BigUnsigned upper=add(result,upper_error);
            const Interval positive{
                Fraction(lower,fixed_bits,lower_negative),
                Fraction(upper,fixed_bits)
            };
            return x.sign()<0 ? negative(positive) : positive;
        }
    }
    throw std::logic_error("FloatAP arctangent series did not terminate");
}

Interval pi_interval(unsigned long work) {
    // Keep a single certified enclosure per thread, at the greatest
    // requested work precision. A tighter enclosure is always valid
    // for a smaller work request, independently of rounding mode.
    struct Cache {
        unsigned long work = 0;
        std::optional<Interval> value;
    };
    static thread_local Cache cache;
    if (cache.value && cache.work >= work) {
        return *cache.value;
    }

    // Machin: pi = 16 atan(1/5) - 4 atan(1/239).
    const Interval a = small_atan(Fraction(1L) / Fraction(5L), work);
    const Interval b = small_atan(Fraction(1L) / Fraction(239L), work);
    const Interval result = subtract(multiply(point(Fraction(16L)), a),
                                     multiply(point(Fraction(4L)), b));
    cache.value = result;
    cache.work = work;
    return result;
}

// An interval for atan(1/2) computed at greater working precision
// remains certified for lower precision requests. Cache the enclosure,
// never a rounded scalar, so all four rounding modes share the result.
Interval atan_half_interval(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<Interval> value;
    };
    static thread_local Cache cache;
    if (cache.value && cache.work >= work) {
        return *cache.value;
    }
    const Interval result = small_atan(Fraction(1L) / Fraction(2L), work);
    cache.value = result;
    cache.work = work;
    return result;
}

Interval atan_interval(Fraction const& x, unsigned long work) {
    // The centered branch has additional certified-constant and rational
    // reduction costs. At moderate precision the direct Taylor kernel is
    // faster for small arguments (notably asin's argument near 2-sqrt(3)).
    // At high precision the saved Taylor iterations justify that overhead.
    const Fraction direct_limit = work <= 512ul
        ? Fraction(1L) / Fraction(3L) : Fraction(1L) / Fraction(4L);
    if (compare(absolute(x), direct_limit) < 0) {
        return small_atan(x, work);
    }
    if (x.sign() < 0) {
        return negative(atan_interval(Fraction{} - x, work));
    }

    const Fraction half_center_limit = Fraction(5L) / Fraction(7L);
    if (compare(x, half_center_limit) <= 0) {
        // atan(x) = atan(1/2) + atan((2*x-1)/(x+2)).
        // On [1/4,5/7] the second argument lies in [-2/9,3/19].
        // For x near 1/2 it is tiny, so one Taylor iteration often
        // replaces dozens or hundreds. Both intervals are certified.
        const Fraction reduced = (Fraction(2L) * x - Fraction(1L)) /
                                 (x + Fraction(2L));
        return add(atan_half_interval(work), small_atan(reduced, work));
    }

    const Interval quarter_pi = divide(pi_interval(work), point(Fraction(4L)));
    if (compare(x, Fraction(2L)) <= 0) {
        // atan(x) = pi/4 + atan((x-1)/(x+1)), |argument| <= 1/3.
        const Fraction reduced = (x - Fraction(1L)) / (x + Fraction(1L));
        return add(quarter_pi, small_atan(reduced, work));
    }
    // atan(x) = pi/2 - atan(1/x) for x > 2.
    return subtract(multiply(point(Fraction(2L)), quarter_pi),
                    small_atan(Fraction(1L) / x, work));
}

std::pair<Interval, Interval> sine_cosine_intervals(Fraction x,
                                                   unsigned long work,
                                                   TrigOutputs outputs) {
    unsigned halvings = 0;
    while (compare(absolute(x), Fraction(1L)) > 0) {
        if (++halvings > 4096u) {
            throw std::overflow_error("FloatAP trigonometric argument too large");
        }
        x.halve();
    }

    // For an argument very close to 1/2, use a cached certified
    // sin(1/2), cos(1/2) and the Taylor remainder bound
    // |f(1/2+d)-f(1/2)-f'(1/2)*d| <= d*d/2 for f=sin,cos.
    // This is valid for EVERY real argument, not merely binary inputs.
    // If d is not tiny relative to work, use the original Taylor kernel.
    // The *dyadic fixed-point* centered path is cheap even at 53-bit
    // target precision (85-bit initial work), once its rigorous cached
    // center is initialized. Keep the original exact neighborhood check
    // and Taylor fallback. Tangent needs BOTH intervals and retains its
    // separate certified path: the cached path is for standalone
    // sine/cosine only.
    if (halvings == 0u && work >= 64ul &&
        outputs != TrigOutputs::both) {
        const Fraction half = Fraction(1L) / Fraction(2L);
        const Fraction delta = absolute(x) - half;
        if (sufficiently_near_half(delta, work)) {
            // A FloatAP input is dyadic. Use a single cached integer
            // scale and outward-directed GMP shifts instead of creating
            // general rational intermediate products every time.
            if (mpz_popcount(mpq_denref(delta.raw())) == 1u) {
                return centered_trigonometric_dyadic(x, delta, work, outputs);
            }
            if (work > std::numeric_limits<unsigned long>::max() - 16ul) {
                throw std::overflow_error("FloatAP centered trigonometric precision overflow");
            }
            const auto& [center_sine, center_cosine] =
                half_trigonometric_intervals(work + 16ul);
            const Fraction remainder = (delta * delta) / Fraction(2L);
            Interval sine = point(Fraction{});
            Interval cosine = point(Fraction{});

            if (outputs != TrigOutputs::cosine) {
                sine = add(center_sine, multiply_by_fraction(center_cosine, delta));
                sine.lower = sine.lower - remainder;
                sine.upper = sine.upper + remainder;
                if (x.sign() < 0) sine = negative(sine);
            }
            if (outputs != TrigOutputs::sine) {
                cosine = subtract(center_cosine,
                                  multiply_by_fraction(center_sine, delta));
                cosine.lower = cosine.lower - remainder;
                cosine.upper = cosine.upper + remainder;
            }
            return {std::move(sine), std::move(cosine)};
        }
    }

    // Double-angle reconstruction needs BOTH values, regardless of
    // which function the caller eventually requests.
    const TrigOutputs kernel_outputs = halvings == 0u ? outputs : TrigOutputs::both;
    auto results = small_sine_cosine(x, work + halvings + 16u, kernel_outputs);
    for (unsigned n = 0; n < halvings; ++n) {
        const Interval sine = multiply(point(Fraction(2L)),
                                       multiply(results.first, results.second));
        const Interval cosine = subtract(multiply(results.second, results.second),
                                          multiply(results.first, results.first));
        results = {sine, cosine};
    }
    return results;
}

Interval square_root_interval(Fraction const& x, unsigned long work) {
    if (x.sign() < 0) {
        throw std::domain_error("FloatAP inverse trigonometric domain");
    }
    if (x.sign() == 0) {
        return point(Fraction{});
    }
    // The argument is an exact positive fraction. Its integer quotient
    // and exact square-root remainder are native; only the two final
    // interval endpoints cross the legacy Fraction boundary.
    if (work > std::numeric_limits<mp_bitcnt_t>::max() / 2u ||
        work > std::numeric_limits<std::size_t>::max() / 2u) {
        throw std::overflow_error("FloatAP radical precision too large");
    }
    const Integer numerator=integer_from_gmp(mpq_numref(x.raw()));
    const Integer denominator=integer_from_gmp(mpq_denref(x.raw()));
    auto [quotient,division_remainder]=divrem(
        shift_left(numerator.magnitude(),static_cast<std::size_t>(2u*work)),
        denominator.magnitude());
    auto [root,square_remainder]=sqrtrem(quotient);
    const Fraction lower(root,work);
    // Equality to the original rational requires BOTH exact integer
    // division and an exact square. No extra root^2 * denominator
    // multiply or canonical Rational/GCD construction is necessary.
    if (division_remainder.is_zero() && square_remainder.is_zero()) {
        return point(lower);
    }
    return {lower,Fraction(add(root,BigUnsigned(1u)),work)};
}

Interval asin_interval(Fraction const& x, unsigned long work) {
    // asin(x) = 2 atan(x / (1 + sqrt(1-x*x))).
    // Unlike atan(x/sqrt(1-x*x)), the denominator remains >= 1.
    const Interval radical = square_root_interval(
        Fraction(1L) - x * x, work);
    const Interval denominator = add(point(Fraction(1L)), radical);
    const Interval argument = divide(point(x), denominator);

    // atan is monotone and globally 1-Lipschitz:
    // 0 <= atan(b)-atan(a) <= b-a for a<=b, because
    // atan'(t)=1/(1+t*t) is between zero and one.
    // Evaluating one endpoint, plus the exact input interval width,
    // avoids a second certified arctangent evaluation. Choose the
    // endpoint nearer to zero to reduce its Taylor-series work.
    const Fraction width = argument.upper - argument.lower;
    if (compare(absolute(argument.lower), absolute(argument.upper)) <= 0) {
        const Interval left = atan_interval(argument.lower, work);
        return {Fraction(2L) * left.lower,
                Fraction(2L) * (left.upper + width)};
    }
    const Interval right = atan_interval(argument.upper, work);
    return {Fraction(2L) * (right.lower - width),
            Fraction(2L) * right.upper};
}

enum class Function {
    sine, cosine, tangent, arctangent, arcsine, arccosine,
    pi, half_pi, negative_half_pi
};

// Replace exact division in repeated near-half tangent evaluations with
// dyadic outward endpoints. The initial general division is computed only
// on a cache miss or when the requested work precision increases.
Interval outward_dyadic_interval(Interval const& input, unsigned long bits) {
    return {
        outward_dyadic_endpoint(input.lower,bits,false),
        outward_dyadic_endpoint(input.upper,bits,true)
    };
}

// Certified linearization of atan around 1/2, for signed inputs.
// Cache an OUTWARD interval of atan(1/2), never a rounded scalar.
// For |delta|<=2^-ceil(work/2) and work>=4, y=|x|=1/2+delta
// lies inside [1/4,3/4]. There
// |atan''(y)|=2*y/(1+y*y)^2 <= 3/2, hence |R|<=delta^2.
// The derivative 4/5 is an EXACT rational, not a second cache.
Interval centered_atan_half_interval(Fraction const& input,
                                      Fraction const& delta,
                                      unsigned long work) {
    // All FloatAP arguments are dyadic. For a sufficiently small delta,
    // evaluate the linear correction (4/5)*delta and the certified
    // quadratic error with integer significands on a cached binary
    // scale. No rational multiplication/division is needed per call.
    // The generic exact-rational branch remains as a fallback.
    const auto rational_fallback = [&]() {
        struct Cache {
            unsigned long work = 0;
            std::optional<Interval> center;
        };
        static thread_local Cache cached_rational;
        if (!cached_rational.center || cached_rational.work < work) {
            if (work > std::numeric_limits<unsigned long>::max() - 40ul) {
                throw std::overflow_error("FloatAP centered atan precision overflow");
            }
            cached_rational.center = outward_dyadic_interval(
                atan_half_interval(work + 40ul), work + 24ul);
            cached_rational.work = work;
        }
        const Fraction linear = Fraction(4L) * delta / Fraction(5L);
        const Fraction remainder = delta * delta;
        const Interval positive{
            cached_rational.center->lower + linear - remainder,
            cached_rational.center->upper + linear + remainder
        };
        return input.sign() < 0 ? negative(positive) : positive;
    };

    if (work > std::numeric_limits<unsigned long>::max() - 40ul) {
        throw std::overflow_error("FloatAP centered atan precision overflow");
    }
    const unsigned long bits = work + 24ul;
    if (static_cast<std::uintmax_t>(bits + 2ul) >
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max()) ||
        mpz_popcount(mpq_denref(delta.raw())) != 1u) {
        return rational_fallback();
    }
    const mp_bitcnt_t delta_shift =
        mpz_scan1(mpq_denref(delta.raw()), 0u);
    if (delta_shift > static_cast<mp_bitcnt_t>(bits + 2ul)) {
        return rational_fallback();
    }

    struct Cache {
        unsigned long work = 0;
        unsigned long bits = 0;
        mpz_t lower, upper;
        Cache() { mpz_inits(lower, upper, nullptr); }
        ~Cache() { mpz_clears(lower, upper, nullptr); }
        Cache(Cache const&) = delete;
        Cache& operator=(Cache const&) = delete;
    };
    static thread_local Cache centered;
    if (centered.work < work) {
        const Interval certified_center = atan_half_interval(work + 40ul);
        mpz_t scaled;
        mpz_init(scaled);
        struct ScaledCleanup {
            mpz_ptr value;
            ~ScaledCleanup() { mpz_clear(value); }
        } scaled_cleanup{scaled};
        mpz_mul_2exp(scaled, mpq_numref(certified_center.lower.raw()),
                     static_cast<mp_bitcnt_t>(bits));
        mpz_fdiv_q(centered.lower, scaled,
                   mpq_denref(certified_center.lower.raw()));
        mpz_mul_2exp(scaled, mpq_numref(certified_center.upper.raw()),
                     static_cast<mp_bitcnt_t>(bits));
        mpz_cdiv_q(centered.upper, scaled,
                   mpq_denref(certified_center.upper.raw()));
        centered.bits = bits;
        centered.work = work;
    }

    mpz_t scaled, linear_lower, linear_upper, lower, upper, remainder;
    mpz_inits(scaled, linear_lower, linear_upper, lower, upper, remainder,
              nullptr);
    struct Cleanup {
        mpz_ptr a,b,c,d,e,f;
        ~Cleanup() { mpz_clears(a,b,c,d,e,f,nullptr); }
    } cleanup{scaled,linear_lower,linear_upper,lower,upper,remainder};

    // delta = n / 2^k. At scale Q=2^centered.bits,
    // (4/5)*delta*Q = n*2^(centered.bits+2-k)/5 exactly.
    // Directed integer quotients give rigorous linear endpoints,
    // including negative deltas.
    const mp_bitcnt_t shift = static_cast<mp_bitcnt_t>(
        centered.bits + 2ul) - delta_shift;
    mpz_mul_2exp(scaled, mpq_numref(delta.raw()), shift);
    mpz_fdiv_q_ui(linear_lower, scaled, 5ul);
    mpz_cdiv_q_ui(linear_upper, scaled, 5ul);

    // |atan''(y)|<=3/2 on [1/4,3/4], so
    // |R| <= delta^2 <= 2^-work for the checked neighborhood.
    // On the cached (possibly tighter) binary scale, 2^-work
    // corresponds exactly to 2^(centered.bits-work) integer units.
    mpz_set_ui(remainder, 1ul);
    mpz_mul_2exp(remainder, remainder, static_cast<mp_bitcnt_t>(
        centered.bits - work));
    mpz_add(lower, centered.lower, linear_lower);
    mpz_sub(lower, lower, remainder);
    mpz_add(upper, centered.upper, linear_upper);
    mpz_add(upper, upper, remainder);
    const Interval positive{Fraction(lower, centered.bits),
                            Fraction(upper, centered.bits)};
    return input.sign() < 0 ? negative(positive) : positive;
}

struct TangentHalfEnclosures {
    Interval tangent;
    Interval derivative;
};

// Compute an interval for tan(1/2) from the existing certified
// sine/cosine kernel and round its endpoints OUTWARD to dyadics.
// Because tan(1/2)>0, the derivative 1+tan(1/2)^2 is monotone
// in the cached interval. Quantize the derivative outward too.
TangentHalfEnclosures const& tangent_half_enclosures(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<TangentHalfEnclosures> enclosures;
    };
    static thread_local Cache cache;
    if (!cache.enclosures || cache.work < work) {
        if (work > std::numeric_limits<unsigned long>::max() - 48ul) {
            throw std::overflow_error("FloatAP tangent cache work precision overflow");
        }
        const unsigned long target_bits = work + 32ul;
        const auto& [sine, cosine] = half_trigonometric_intervals(work + 48ul);
        const Interval tangent_ratio_bounds = divide(sine, cosine);
        const Interval cached_tangent_bounds = outward_dyadic_interval(tangent_ratio_bounds, target_bits);
        if (cached_tangent_bounds.lower.sign() <= 0) {
            throw std::logic_error("FloatAP invalid cached tangent enclosure");
        }
        const Interval raw_derivative{
            Fraction(1L) + cached_tangent_bounds.lower * cached_tangent_bounds.lower,
            Fraction(1L) + cached_tangent_bounds.upper * cached_tangent_bounds.upper
        };
        TangentHalfEnclosures fresh_enclosures{
            cached_tangent_bounds,
            outward_dyadic_interval(raw_derivative, target_bits)
        };
        cache.enclosures = std::move(fresh_enclosures);
        cache.work = work;
    }
    return *cache.enclosures;
}

// On |x-1/2| <= 2^-ceil(work/2), Taylor's theorem gives
// tan(x) = T + (1+T^2)*delta + R, |R|<=6*delta^2.
// Indeed x stays in [1/4,3/4]; here |sin(x)|<=3/4 and
// cos(x)>=1-x*x/2 >= 23/32 > 1/2, hence |tan(x)|<=3/2,
// sec(x)^2<=4, |tan''(x)|=|2*tan(x)*sec(x)^2|<=12.
// Interval arithmetic with cached dyadic T and T' therefore
// gives certified bounds for any sign and rounding mode.
Interval centered_tangent_interval(Fraction const& x,
                                   Fraction const& delta,
                                   unsigned long work) {
    const auto& cache = tangent_half_enclosures(work);
    Interval result = add(cache.tangent,
                          multiply_by_fraction(cache.derivative, delta));
    const Fraction remainder = Fraction(6L) * delta * delta;
    result.lower = result.lower - remainder;
    result.upper = result.upper + remainder;
    return x.sign() < 0 ? negative(result) : result;
}

// For y=|x| close to 1/2, the inverse-trigonometric functions admit
// a certified local enclosure: asin(1/2+d)=pi/6 + (2/sqrt(3))*d + R,
// |R|<=2*d*d. This avoids the radical, rational quotient and atan
// sequence of the general inverse-trigonometric implementation.
// Cache OUTWARD enclosures, never rounded center values.
struct InverseHalfEnclosures {
    Interval arcsine_half;
    Interval arccosine_positive_half;
    Interval arccosine_negative_half;
    Interval slope_half;
};

InverseHalfEnclosures const& inverse_half_enclosures(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<InverseHalfEnclosures> bounds;
    };
    static thread_local Cache center_cache;
    if (!center_cache.bounds || center_cache.work < work) {
        if (work > std::numeric_limits<unsigned long>::max() - 40ul) {
            throw std::overflow_error("FloatAP inverse-half precision overflow");
        }
        const unsigned long dyadic_bits = work + 24ul;
        const Interval certified_pi = pi_interval(work + 40ul);
        const Interval asin_half_bounds =
            divide(certified_pi, point(Fraction(6L)));
        const Interval acos_positive_bounds =
            divide(certified_pi, point(Fraction(3L)));
        const Interval acos_negative_bounds =
            multiply(point(Fraction(2L)), acos_positive_bounds);
        const Interval certified_sqrt_three =
            square_root_interval(Fraction(3L), work + 40ul);
        const Interval slope_bounds =
            divide(point(Fraction(2L)), certified_sqrt_three);
        center_cache.bounds = InverseHalfEnclosures{
            outward_dyadic_interval(asin_half_bounds, dyadic_bits),
            outward_dyadic_interval(acos_positive_bounds, dyadic_bits),
            outward_dyadic_interval(acos_negative_bounds, dyadic_bits),
            outward_dyadic_interval(slope_bounds, dyadic_bits)
        };
        center_cache.work = work;
    }
    return *center_cache.bounds;
}

// When |d|<=2^-ceil(work/2) and work>=4, y=1/2+d is in
// [1/4,3/4]. There asin''(y)=y/(1-y*y)^(3/2) <24/7<4:
// y<=3/4, 1-y*y>=7/16 and sqrt(1-y*y)>1/2.
// Taylor's theorem gives |R| <= 2*d*d, for BOTH asin and
// acos. Use oddness of asin and acos(-y)=pi-acos(y).
Interval centered_inverse_half_interval(Fraction const& input,
                                        Fraction const& offset,
                                        unsigned long work,
                                        bool evaluate_arccosine) {
    const auto& cached = inverse_half_enclosures(work);
    const Interval linear_term_interval =
        multiply_by_fraction(cached.slope_half, offset);

    Interval enclosure;
    if (!evaluate_arccosine) {
        enclosure = add(cached.arcsine_half, linear_term_interval);
    } else if (input.sign() >= 0) {
        enclosure = subtract(cached.arccosine_positive_half,
                             linear_term_interval);
    } else {
        enclosure = add(cached.arccosine_negative_half,
                        linear_term_interval);
    }

    const Fraction maximum_remainder = Fraction(2L) * offset * offset;
    enclosure.lower = enclosure.lower - maximum_remainder;
    enclosure.upper = enclosure.upper + maximum_remainder;
    return (!evaluate_arccosine && input.sign() < 0)
        ? negative(enclosure) : enclosure;
}

Interval evaluate(Function function, Fraction const& x, unsigned long work) {
    switch (function) {
        case Function::sine:
            return sine_cosine_intervals(x, work, TrigOutputs::sine).first;
        case Function::cosine:
            return sine_cosine_intervals(x, work, TrigOutputs::cosine).second;
        case Function::tangent: {
            // Short-circuit repeated near-half evaluations without
            // recomputing sin/cos and dividing two non-dyadic intervals.
            // Outside this exact precision-dependent neighborhood, retain
            // the complete general-purpose certified path and pole checks.
            const Fraction delta = absolute(x) - Fraction(1L) / Fraction(2L);
            if (sufficiently_near_half(delta, work)) {
                return centered_tangent_interval(x, delta, work);
            }
            const auto [sine, cosine] =
                sine_cosine_intervals(x, work, TrigOutputs::both);
            return divide(sine, cosine);
        }
        case Function::arctangent: {
            const Fraction delta = absolute(x) -
                Fraction(1L) / Fraction(2L);
            if (work >= 4ul && sufficiently_near_half(delta, work)) {
                return centered_atan_half_interval(x, delta, work);
            }
            return atan_interval(x, work);
        }
        case Function::arcsine: {
            const Fraction offset = absolute(x) -
                Fraction(1L) / Fraction(2L);
            if (work >= 4ul && sufficiently_near_half(offset, work)) {
                return centered_inverse_half_interval(x, offset, work, false);
            }
            return asin_interval(x, work);
        }
        case Function::arccosine: {
            const Fraction offset = absolute(x) -
                Fraction(1L) / Fraction(2L);
            if (work >= 4ul && sufficiently_near_half(offset, work)) {
                return centered_inverse_half_interval(x, offset, work, true);
            }
            return subtract(divide(pi_interval(work), point(Fraction(2L))),
                            asin_interval(x, work));
        }
        case Function::pi:
            return pi_interval(work);
        case Function::half_pi:
            return divide(pi_interval(work), point(Fraction(2L)));
        case Function::negative_half_pi:
            return negative(divide(pi_interval(work), point(Fraction(2L))));
        default:
            throw std::invalid_argument("FloatAP unknown trigonometric function");
    }
}

// A result is returned only when all real values in its certified
// enclosure round to exactly the same FloatAP. This includes directed
// rounding and correctly handles the nearest-even halfway boundaries.
FloatAP certified(Function function, Fraction const& x,
                  FloatAP::RoundingModeType mode, FloatAP::PrecisionType precision) {
    const unsigned long bits = static_cast<unsigned long>(precision.bits());
    unsigned long extra = 32u;
    for (;;) {
        if (extra > std::numeric_limits<unsigned long>::max() - bits - 32u) {
            throw std::overflow_error("FloatAP trigonometric precision too large");
        }

        const Interval bounds = evaluate(function, x, bits + extra);
        const FloatAP low = bounds.lower.rounded(bits,mode);
        const FloatAP high = bounds.upper.rounded(bits,mode);
        if (low == high) {
            return low;
        }

        if (extra > std::numeric_limits<unsigned long>::max() / 2u) {
            throw std::overflow_error("FloatAP trigonometric refinement overflow");
        }
        extra *= 2u;
    }
}

FloatAP unary(Function function, FloatAP::RoundingModeType mode,
              FloatAP const& x) {
    const unsigned long bits = static_cast<unsigned long>(x.precision().bits());
    if (x.is_nan()) {
        return FloatAP::nan(bits);
    }

    if (x.is_inf()) {
        if (function == Function::arctangent) {
            return certified(x.negative() ? Function::negative_half_pi :
                           Function::half_pi, Fraction{}, mode, x.precision());
        }
        return FloatAP::nan(bits);
    }

    if (x.is_zero()) {
        switch (function) {
            case Function::sine:
            case Function::tangent:
            case Function::arctangent:
            case Function::arcsine:
                return FloatAP::zero(bits, x.negative());
            case Function::cosine:
                return FloatAP(std::int64_t{1}, x.precision());
            case Function::arccosine:
                return certified(Function::arccosine, Fraction{}, mode, x.precision());
            case Function::pi:
            case Function::half_pi:
            case Function::negative_half_pi:
            default:
                break;
        }
    }

    const Fraction argument = rational_argument(x);
    if (function == Function::arcsine || function == Function::arccosine) {
        if (compare(absolute(argument), Fraction(1L)) > 0) {
            return FloatAP::nan(bits);
        }
        if (compare(argument, Fraction(1L)) == 0) {
            if (function == Function::arccosine) {
                return FloatAP::zero(bits);
            }
            return certified(Function::half_pi, Fraction{}, mode, x.precision());
        }
        if (compare(argument, Fraction(-1L)) == 0) {
            return certified(function == Function::arcsine ?
                             Function::negative_half_pi : Function::pi,
                             Fraction{}, mode, x.precision());
        }
    }
    return certified(function, argument, mode, x.precision());
}

} // namespace FloatAPTrigonometricDetail

} // namespace

FloatAP exp(FloatAP::RoundingModeType mode, FloatAP const& x) {
    using namespace FloatAPExpLogDetail;
    const unsigned long bits = static_cast<unsigned long>(x.precision().bits());
    if (x.is_nan()) {
        return FloatAP::nan(bits);
    }
    if (x.is_inf()) {
        return x.negative() ? FloatAP::zero(bits) : FloatAP::infinity(bits);
    }
    if (x.is_zero()) {
        return FloatAP(std::int64_t{1}, x.precision());
    }
    const Fraction input = rational_argument(x);
    return certified([&](unsigned long work) {
        return exp_half_centered_or_general(input, work);
    }, mode, x.precision());
}

FloatAP log(FloatAP::RoundingModeType mode, FloatAP const& x) {
    using namespace FloatAPExpLogDetail;
    const unsigned long bits = static_cast<unsigned long>(x.precision().bits());
    if (x.is_nan() || (x.negative() && !x.is_zero())) {
        return FloatAP::nan(bits);
    }
    if (x.is_zero()) {
        return FloatAP::infinity(bits, true);
    }
    if (x.is_inf()) {
        return FloatAP::infinity(bits);
    }
    const Fraction input = rational_argument(x);
    if (compare(input, Fraction(1)) == 0) {
        return FloatAP::zero(bits);
    }
    return certified([&](unsigned long work) {
        return log_half_centered_or_general(input, work);
    }, mode, x.precision());
}

FloatAP sin(FloatAP::RoundingModeType mode, FloatAP const& x) {
    using namespace FloatAPTrigonometricDetail;
    return unary(Function::sine, mode, x);
}

FloatAP cos(FloatAP::RoundingModeType mode, FloatAP const& x) {
    using namespace FloatAPTrigonometricDetail;
    return unary(Function::cosine, mode, x);
}

FloatAP tan(FloatAP::RoundingModeType mode, FloatAP const& x) {
    using namespace FloatAPTrigonometricDetail;
    return unary(Function::tangent, mode, x);
}

FloatAP atan(FloatAP::RoundingModeType mode, FloatAP const& x) {
    using namespace FloatAPTrigonometricDetail;
    return unary(Function::arctangent, mode, x);
}

FloatAP asin(FloatAP::RoundingModeType mode, FloatAP const& x) {
    using namespace FloatAPTrigonometricDetail;
    return unary(Function::arcsine, mode, x);
}

FloatAP acos(FloatAP::RoundingModeType mode, FloatAP const& x) {
    using namespace FloatAPTrigonometricDetail;
    return unary(Function::arccosine, mode, x);
}

FloatAP Float<ArbitraryPrecision>::pi(RoundingModeType mode, PrecisionType precision) {
    using namespace FloatAPTrigonometricDetail;
    unsigned index;
    switch (mode.value()) {
        case Rounding::ToNearest: index = 0u; break;
        case Rounding::Downward: index = 1u; break;
        case Rounding::Upward: index = 2u; break;
        case Rounding::TowardZero: index = 3u; break;
        default: throw std::invalid_argument("Invalid FloatAP pi rounding mode");
    }

    // Unlike the certified interval cache, the rounded result is specific
    // to BOTH target precision and rounding mode. Keep four values for
    // the most recently used precision per thread; no global lock, no
    // dependence on a mutable default rounding mode.
    struct Cache {
        unsigned long bits = 0;
        std::optional<FloatAP> rounded[4];
    };
    static thread_local Cache cache;
    const unsigned long bits = static_cast<unsigned long>(precision.bits());
    if (cache.bits != bits) {
        for (auto& entry : cache.rounded) {
            entry.reset();
        }
        cache.bits = bits;
    }
    auto& entry = cache.rounded[index];
    if (!entry) {
        entry.emplace(certified(Function::pi, Fraction{}, mode, precision));
    }
    return *entry;
}

} // namespace Ariadne
