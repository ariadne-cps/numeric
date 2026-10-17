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
#include "numeric/rational.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/double.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/rounding.hpp"
#include "numeric/positive.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/decimal.hpp"

#include <algorithm>
#include <bit>
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

namespace Ariadne {

// Round an already-populated, nonzero positive GMP significand in-place.
// Keeping the destination as the accumulator avoids extra mpz allocations
// and signed-integer copies in add_sub() and mul(). The sign is stored
// separately, and the rounding logic is shared with from_ratio(n,1,...).
void FloatAP::round_scaled_integer(FloatAP& result, Exponent scale, Rounding mode) {
    const auto precision = result._precision;
    const Exponent shift = checked_sub(bit_length(result._significand),
                                       static_cast<Exponent>(precision));
    result._exponent = checked_add(scale, shift);
    if (shift <= 0) {
        mpz_mul_2exp(result._significand, result._significand,
                     shift_amount(checked_sub(0, shift)));
        return;
    }

    const mp_bitcnt_t dropped = shift_amount(shift);
    const bool inexact = mpz_scan1(result._significand, 0u) < dropped;
    bool increment = false;
    if (inexact) {
        switch (mode.value()) {
            case Rounding::ToNearest:
                increment = mpz_tstbit(result._significand, dropped - 1u) != 0 &&
                    (mpz_tstbit(result._significand, dropped) != 0 ||
                     mpz_scan1(result._significand, 0u) < dropped - 1u);
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
    mpz_tdiv_q_2exp(result._significand, result._significand, dropped);
    if (increment) mpz_add_ui(result._significand, result._significand, 1u);
    if (bit_length(result._significand) > static_cast<Exponent>(precision)) {
        mpz_tdiv_q_2exp(result._significand, result._significand, 1u);
        // Crossing the maximal finite exponent by a rounded carry yields infinity.
        if (result._exponent == std::numeric_limits<Exponent>::max()) {
            result = infinity(precision, result._negative);
            return;
        }
        result._exponent = checked_add(result._exponent, 1);
    }
}

FloatAP FloatAP::from_ratio(mpz_srcptr numerator, mpz_srcptr denominator,
Exponent scale, unsigned long precision, Rounding mode,
bool negative_zero, bool invert_sign) {
    FloatAP result(precision);
    if (mpz_sgn(denominator) == 0) throw std::domain_error("FloatAP zero denominator");
    if (mpz_sgn(numerator) == 0) {
        result.set_zero(negative_zero);
        return result;
    }
    result._negative = ((mpz_sgn(numerator) < 0) !=
                        (mpz_sgn(denominator) < 0)) != invert_sign;

    // General rational normalization and division are unnecessary when
    // the denominator is one. Round the existing magnitude in-place.
    if (mpz_cmp_ui(denominator, 1ul) == 0) {
        mpz_abs(result._significand, numerator);
        round_scaled_integer(result, scale, mode);
        return result;
    }

    // Non-unit denominator: reuse input magnitudes whenever already
    // positive. The quotient is produced directly in the final FloatAP
    // significand, while GMP needs only the remainder and aligned scratch.
    mpz_t abs_n, abs_d, remainder, scratch;
    mpz_inits(abs_n, abs_d, remainder, scratch, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c, d;
        ~Cleanup() { mpz_clears(a, b, c, d, nullptr); }
    } cleanup{abs_n, abs_d, remainder, scratch};
    mpz_srcptr n = numerator;
    mpz_srcptr d = denominator;
    if (mpz_sgn(numerator) < 0) {
        mpz_abs(abs_n, numerator);
        n = abs_n;
    }
    if (mpz_sgn(denominator) < 0) {
        mpz_abs(abs_d, denominator);
        d = abs_d;
    }

    Exponent k = checked_sub(bit_length(n), bit_length(d));
    if (k >= 0) {
        mpz_mul_2exp(scratch, d, shift_amount(k));
        if (mpz_cmp(n, scratch) < 0) k = checked_sub(k, 1);
    } else {
        mpz_mul_2exp(scratch, n, shift_amount(checked_sub(0, k)));
        if (mpz_cmp(scratch, d) < 0) k = checked_sub(k, 1);
    }

    const Exponent shift = checked_sub(static_cast<Exponent>(precision - 1), k);
    mpz_srcptr divisor = d;
    if (shift >= 0) {
        mpz_mul_2exp(scratch, n, shift_amount(shift));
        mpz_tdiv_qr(result._significand, remainder, scratch, d);
    } else {
        mpz_mul_2exp(scratch, d, shift_amount(checked_sub(0, shift)));
        divisor = scratch;
        mpz_tdiv_qr(result._significand, remainder, n, scratch);
    }

    bool increment = false;
    if (mpz_size(remainder) != 0) {
        switch (mode.value()) {
            case Rounding::ToNearest:
                mpz_mul_2exp(remainder, remainder, 1);
                increment = mpz_cmp(remainder, divisor) > 0 ||
                    (mpz_cmp(remainder, divisor) == 0 &&
                     mpz_odd_p(result._significand));
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
    if (increment) mpz_add_ui(result._significand, result._significand, 1);
    result._exponent = checked_add(scale, checked_sub(k, static_cast<Exponent>(precision - 1)));
    if (bit_length(result._significand) > static_cast<Exponent>(precision)) {
        mpz_tdiv_q_2exp(result._significand, result._significand, 1);
        if (result._exponent == std::numeric_limits<Exponent>::max()) {
            return infinity(precision, result._negative);
        }
        result._exponent = checked_add(result._exponent, 1);
    }
    return result;
}

FloatAP FloatAP::from_integer(mpz_srcptr n, unsigned long p, Rounding mode) {
    mpz_t one;
    mpz_init_set_ui(one, 1);
    struct Cleanup {
        mpz_ptr z;
        ~Cleanup() {
            mpz_clear(z);
        }
    }
    cleanup { one };
    return from_ratio(n, one, 0, p, mode);
}

FloatAP FloatAP::from_rational(mpq_srcptr r, unsigned long p, Rounding mode) {
    // Dyadic rationals can be rounded as integer significands without a
    // general GMP quotient/remainder division. Both paths are exact.
    mpz_srcptr denominator = mpq_denref(r);
    if (mpz_popcount(denominator) == 1u) {
        const mp_bitcnt_t power = mpz_scan1(denominator, 0u);
        if (static_cast<std::uintmax_t>(power) <=
            static_cast<std::uintmax_t>(std::numeric_limits<Exponent>::max())) {
            FloatAP result(p);
            mpz_srcptr numerator = mpq_numref(r);
            if (mpz_sgn(numerator) == 0) return result;
            result._negative = mpz_sgn(numerator) < 0;
            mpz_abs(result._significand, numerator);
            round_scaled_integer(result, -static_cast<Exponent>(power), mode);
            return result;
        }
    }
    return from_ratio(mpq_numref(r), denominator, 0, p, mode);
}

void FloatAP::to_rational(mpq_ptr out) const {
    if (!is_finite()) throw std::domain_error("FloatAP nonfinite to Rational");
    mpz_set(mpq_numref(out), _significand);
    if (_negative) mpz_neg(mpq_numref(out), mpq_numref(out));
    mpz_set_ui(mpq_denref(out), 1);
    if (_exponent > 0) mpz_mul_2exp(mpq_numref(out), mpq_numref(out), shift_amount(_exponent));
    if (_exponent < 0) mpz_mul_2exp(mpq_denref(out), mpq_denref(out), shift_amount(checked_sub(0, _exponent)));
    mpq_canonicalize(out);
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
        // Keep the certified signed-sticky path for enormous exponent gaps.
        mpz_t exact, one;
        mpz_inits(exact, one, nullptr);
        struct Cleanup {
            mpz_ptr a, b;
            ~Cleanup() { mpz_clears(a, b, nullptr); }
        } cleanup{exact, one};
        const Exponent padding = checked_add(
            checked_add(static_cast<Exponent>(high._precision),
                        static_cast<Exponent>(p)), 2);
        mpz_mul_2exp(exact, high._significand, shift_amount(padding));
        if (high_negative) mpz_neg(exact, exact);
        if (low_negative) mpz_sub_ui(exact, exact, 1u);
        else mpz_add_ui(exact, exact, 1u);
        mpz_set_ui(one, 1u);
        mpz_mul_2exp(one, one, shift_amount(padding));
        return from_ratio(exact, one, high._exponent, p, mode);
    }

    // For ordinary gaps, accumulate directly into the returned value.
    // The unshifted input significands are positive; sign is metadata.
    FloatAP result(p);
    result._negative = high_negative;
    mpz_mul_2exp(result._significand, high._significand,
                 shift_amount(checked_sub(high._exponent, low._exponent)));
    if (a._negative == b_negative) {
        mpz_add(result._significand, result._significand, low._significand);
    } else {
        mpz_sub(result._significand, result._significand, low._significand);
        if (mpz_sgn(result._significand) < 0) {
            mpz_neg(result._significand, result._significand);
            result._negative = !result._negative;
        }
    }
    if (mpz_size(result._significand) == 0) {
        result.set_zero(mode == Rounding::Downward);
        return result;
    }
    round_scaled_integer(result, low._exponent, mode);
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
    mpz_mul(result._significand, a._significand, b._significand);
    round_scaled_integer(result, checked_add(a._exponent, b._exponent), mode);
    return result;
}

FloatAP FloatAP::div(FloatAP const& a, FloatAP const& b,
unsigned long p, Rounding mode) {
    if (a.is_nan() || b.is_nan() || (a.is_zero() && b.is_zero()) || (a.is_inf() && b.is_inf())) return nan(p);
    const bool neg = a._negative != b._negative;
    if (a.is_inf() || b.is_zero()) return infinity(p, neg);
    if (a.is_zero() || b.is_inf()) return zero(p, neg);
    // Sign is metadata. The two significands are positive, immutable
    // GMP integers; do not allocate, copy or negate the numerator merely
    // to pass it to the exact quotient-and-remainder rounding kernel.
    return from_ratio(a._significand, b._significand,
                      checked_sub(a._exponent, b._exponent), p, mode,
                      false, neg);
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
    mpz_t product;
    mpz_init(product);
    struct Cleanup {
        mpz_ptr z;
        ~Cleanup() {
            mpz_clear(z);
        }
    }
    cleanup { product };
    mpz_mul(product, a._significand, b._significand);
    FloatAP exact(static_cast<unsigned long>(mpz_sizeinbase(product, 2)));
    mpz_set(exact._significand, product);
    exact._negative = neg;
    exact._exponent = checked_add(a._exponent, b._exponent);
    return add(exact, c, p, mode);
}

namespace {
// A signed 64-bit integer is imported directly, without a decimal
// intermediate. In particular, INT64_MIN is handled without negating
// a signed value outside its range.
FloatAP float_from_int64(std::int64_t n, ArbitraryPrecision precision) {
    mpz_t value;
    mpz_init(value);
    struct Cleanup {
        mpz_ptr value;
        ~Cleanup() {
            mpz_clear(value);
        }
    }
    cleanup { value };

    std::uint64_t magnitude = static_cast<std::uint64_t>(n);
    if (n < 0) {
        magnitude = std::uint64_t { 0 } - magnitude;
    }
    mpz_import(value, 1, 1, sizeof(magnitude), 0, 0,& magnitude);
    if (n < 0) {
        mpz_neg(value, value);
    }
    return FloatAP::from_integer(value, static_cast<unsigned long>(precision.bits()),
    FloatAP::ROUND_TO_NEAREST);
}
}

FloatAP::Float(std::int64_t n, PrecisionType pr)
: Float(float_from_int64(n, pr)) {
}

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
        mpz_t numerator, one;
        mpz_init_set(numerator, value._significand);
        mpz_init_set_ui(one, 1u);
        struct Cleanup {
            mpz_ptr a, b;
            ~Cleanup() { mpz_clears(a, b, nullptr); }
        } cleanup{numerator, one};
        if (value._negative) mpz_neg(numerator, numerator);
        *this = from_ratio(numerator, one, value._exponent,
                           static_cast<unsigned long>(precision.bits()), mode);
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
    mpf_t exact;
    mpf_init2(exact, static_cast<mp_bitcnt_t>(_precision));
    struct Cleanup {
        mpf_ptr value;
        ~Cleanup() { mpf_clear(value); }
    } cleanup{exact};
    mpf_set_z(exact, _significand);
    if (_exponent >= 0) {
        mpf_mul_2exp(exact, exact, shift_amount(_exponent));
    } else {
        mpf_div_2exp(exact, exact, shift_amount(checked_sub(0, _exponent)));
    }
    if (_negative) {
        mpf_neg(exact, exact);
    }
    return Dyadic(exact);
}

Integer cast_integer(FloatAP const& value) {
    if (!value.is_finite()) {
        throw std::domain_error("FloatAP cannot cast a nonfinite value to Integer");
    }
    const Rational exact = static_cast<Rational>(value);
    mpq_srcptr q = exact.get_mpq();
    if (mpz_cmp_ui(mpq_denref(q), 1u) != 0) {
        throw std::domain_error("FloatAP is not an integer");
    }
    return Integer(mpq_numref(q));
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
    const Rational rational=static_cast<Rational>(x);
    mpq_srcptr q=rational.get_mpq();
    mpz_t numerator,denominator,scale,quotient,remainder;
    mpz_inits(numerator,denominator,scale,quotient,remainder,nullptr);
    struct Clear {
        mpz_ptr a,b,c,d,e;
        ~Clear(){mpz_clears(a,b,c,d,e,nullptr);}
    } clear{numerator,denominator,scale,quotient,remainder};
    const bool neg=x.negative();
    mpz_abs(numerator,mpq_numref(q));
    mpz_set(denominator,mpq_denref(q));
    mpz_ui_pow_ui(scale,10ul,static_cast<unsigned long>(places));
    mpz_mul(numerator,numerator,scale);
    mpz_tdiv_qr(quotient,remainder,numerator,denominator);
    if (mpz_sgn(remainder)!=0) {
        bool increase=false;
        switch(mode.value()) {
            case FloatAP::Rounding::ToNearest: {
                mpz_mul_ui(remainder,remainder,2u);
                const int order=mpz_cmp(remainder,denominator);
                increase=order>0 || (order==0 && mpz_odd_p(quotient));
                break;
            }
            case FloatAP::Rounding::Downward: increase=neg; break;
            case FloatAP::Rounding::Upward: increase=!neg; break;
            case FloatAP::Rounding::TowardZero: break;
            default: throw std::invalid_argument("Invalid FloatAP decimal rounding mode");
        }
        if(increase) mpz_add_ui(quotient,quotient,1);
    }
    const std::size_t capacity=mpz_sizeinbase(quotient,10)+2;
    std::string digits(capacity,'\0');
    mpz_get_str(digits.data(),10,quotient);
    digits.resize(std::char_traits<char>::length(digits.c_str()));
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
        const double estimated=(static_cast<double>(mpz_sizeinbase(x.significand(),2))-1.0
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
    const bool all_ones = mpz_popcount(r._significand) == r._precision;
    const bool first = mpz_popcount(r._significand) == 1;
    if (r._negative) {
        if (first) {
            if (r._exponent == std::numeric_limits<FloatAP::Exponent>::min()) {
                return FloatAP::zero(r._precision, true);
            }
            mpz_setbit(r._significand, static_cast<mp_bitcnt_t>(r._precision));
            mpz_sub_ui(r._significand, r._significand, 1);
            --r._exponent;
        } else {
            mpz_sub_ui(r._significand, r._significand, 1);
        }
    } else if (all_ones) {
        if (r._exponent == std::numeric_limits<FloatAP::Exponent>::max()) {
            return FloatAP::inf(r.precision());
        }
        mpz_set_ui(r._significand, 0);
        mpz_setbit(r._significand, static_cast<mp_bitcnt_t>(r._precision-1));
        ++r._exponent;
    } else {
        mpz_add_ui(r._significand, r._significand, 1);
    }
    return r;
}

FloatAP next(RoundDownward, FloatAP const& x) {
    if (x.is_nan() || (x.is_inf() && x.negative())) return x;
    if (x.is_inf()) return FloatAP::max(x.precision());
    if (x.is_zero()) return neg(FloatAP::min(x.precision()));
    FloatAP r(x);
    const bool all_ones = mpz_popcount(r._significand) == r._precision;
    const bool first = mpz_popcount(r._significand) == 1;
    if (r._negative) {
        if (all_ones) {
            if (r._exponent == std::numeric_limits<FloatAP::Exponent>::max()) {
                return FloatAP::inf(Sign::NEGATIVE, r.precision());
            }
            mpz_set_ui(r._significand, 0);
            mpz_setbit(r._significand, static_cast<mp_bitcnt_t>(r._precision-1));
            ++r._exponent;
        } else {
            mpz_add_ui(r._significand, r._significand, 1);
        }
    } else if (first) {
        if (r._exponent == std::numeric_limits<FloatAP::Exponent>::min()) {
            return FloatAP::zero(r._precision);
        }
        mpz_setbit(r._significand, static_cast<mp_bitcnt_t>(r._precision));
        mpz_sub_ui(r._significand, r._significand, 1);
        --r._exponent;
    } else {
        mpz_sub_ui(r._significand, r._significand, 1);
    }
    return r;
}

// Integral rounding uses the exact GMP significand and binary exponent.
// It does not convert an extreme exponent into an enormous Rational.
FloatAP floor(FloatAP const& x) {
    if (!x.is_finite() || x.is_zero() || x._exponent >= 0) return x;
    const std::uint64_t shift = std::uint64_t{0} -
        static_cast<std::uint64_t>(x._exponent);
    if (shift >= mpz_sizeinbase(x._significand, 2)) {
        return FloatAP(x._negative ? std::int64_t{-1} : std::int64_t{0}, x.precision());
    }
    mpz_t integer;
    mpz_init(integer);
    struct Cleanup { mpz_ptr p; ~Cleanup() { mpz_clear(p); } } cleanup{integer};
    const auto k = FloatAP::shift_amount(static_cast<FloatAP::Exponent>(shift));
    if (x._negative) {
        mpz_cdiv_q_2exp(integer, x._significand, k);
        mpz_neg(integer, integer);
    } else {
        mpz_fdiv_q_2exp(integer, x._significand, k);
    }
    return FloatAP::from_integer(integer, x._precision, FloatAP::ROUND_TO_NEAREST);
}

FloatAP ceil(FloatAP const& x) {
    if (!x.is_finite() || x.is_zero() || x._exponent >= 0) return x;
    const std::uint64_t shift = std::uint64_t{0} -
        static_cast<std::uint64_t>(x._exponent);
    if (shift >= mpz_sizeinbase(x._significand, 2)) {
        return FloatAP(x._negative ? std::int64_t{0} : std::int64_t{1}, x.precision());
    }
    mpz_t integer;
    mpz_init(integer);
    struct Cleanup { mpz_ptr p; ~Cleanup() { mpz_clear(p); } } cleanup{integer};
    const auto k = FloatAP::shift_amount(static_cast<FloatAP::Exponent>(shift));
    if (x._negative) {
        mpz_fdiv_q_2exp(integer, x._significand, k);
        mpz_neg(integer, integer);
    } else {
        mpz_cdiv_q_2exp(integer, x._significand, k);
    }
    return FloatAP::from_integer(integer, x._precision, FloatAP::ROUND_TO_NEAREST);
}

FloatAP round(FloatAP const& x) {
    if (!x.is_finite() || x.is_zero() || x._exponent >= 0) return x;
    const std::uint64_t shift = std::uint64_t{0} -
        static_cast<std::uint64_t>(x._exponent);
    const std::size_t bits = mpz_sizeinbase(x._significand, 2);
    if (shift > bits) return FloatAP::zero(x._precision, x._negative);
    mpz_t integer;
    mpz_init(integer);
    struct Cleanup { mpz_ptr p; ~Cleanup() { mpz_clear(p); } } cleanup{integer};
    const auto k = FloatAP::shift_amount(static_cast<FloatAP::Exponent>(shift));
    mpz_tdiv_q_2exp(integer, x._significand, k);
    // mpfr_round uses nearest integer with halfway cases away from zero.
    if (mpz_tstbit(x._significand, k-1) != 0) {
        mpz_add_ui(integer, integer, 1);
    }
    if (x._negative) mpz_neg(integer, integer);
    if (mpz_sgn(integer) == 0) return FloatAP::zero(x._precision, x._negative);
    return FloatAP::from_integer(integer, x._precision, FloatAP::ROUND_TO_NEAREST);
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

// Outward-rounded FloatAP operations are defined by Bounds<F> in
// float_bounds.hpp; do not define an independent competing overload here.
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
: Float(from_integer(n.get_mpz(), static_cast<unsigned long>(pr.bits()), mode)) {
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

FloatAP::Float(Rational const& q, RoundingModeType mode, PrecisionType pr)
: Float(pr) {
    if (rational_is_nan(q)) {
        *this=nan(_precision);
    }
    else if (rational_is_inf(q)) {
        *this=infinity(_precision,mpz_sgn(mpq_numref(q.get_mpq()))<0);
    } else {
        *this=from_rational(q.get_mpq(),_precision,mode);
    }
}

FloatAP::operator Rational() const {
    if (is_nan()) return Rational::nan();
    if (is_inf()) return Rational::inf(negative() ? Sign::NEGATIVE : Sign::POSITIVE);
    Rational q;
    to_rational(q._mpq);
    return q;
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
    const Exponent top = static_cast<Exponent>(mpz_sizeinbase(_significand, 2) - 1u);
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

    // GMP truncates this normalized 53-bit significand, regardless of the
    // stored precision. The conversion and scaling to an integer are exact.
    long unused_exponent = 0;
    const double leading_fraction = mpz_get_d_2exp(&unused_exponent, _significand);
    const std::uint64_t leading = static_cast<std::uint64_t>(
        std::ldexp(leading_fraction, static_cast<int>(digits)));
    const Exponent subnormal_shift = highest < min_normal
                                    ? min_normal - highest : 0;
    std::uint64_t mantissa = leading >> static_cast<unsigned>(subnormal_shift);

    // The original significand bits below this cut decide every rounding
    // mode, including IEEE halfway ties and the transition to subnormals.
    const Exponent discarded = top - (digits - 1) + subnormal_shift;
    if (discarded > 0) {
        if (static_cast<std::uintmax_t>(discarded) >
            static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
            throw std::overflow_error("FloatAP conversion shift exceeds GMP limit");
        }
        const mp_bitcnt_t dropped = static_cast<mp_bitcnt_t>(discarded);
        bool increment = false;
        switch (mode) {
            case Rounding::ToNearest:
                increment = mpz_tstbit(_significand, dropped - 1u) != 0 &&
                    ((mantissa & 1u) != 0u ||
                     mpz_scan1(_significand, 0u) < dropped - 1u);
                break;
            case Rounding::Downward:
                increment = _negative &&
                    mpz_scan1(_significand, 0u) < dropped;
                break;
            case Rounding::Upward:
                increment = !_negative &&
                    mpz_scan1(_significand, 0u) < dropped;
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

    // Compare significand * 2^exponent directly. A different most
    // significant binary position decides the order without allocating
    // any GMP temporaries or constructing Rational values.
    int magnitude = 0;
    if (a._exponent == b._exponent) {
        magnitude = mpz_cmp(a._significand, b._significand);
    } else {
        const std::size_t a_bits = mpz_sizeinbase(a._significand, 2);
        const std::size_t b_bits = mpz_sizeinbase(b._significand, 2);

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
                mpz_t truncated;
                mpz_init(truncated);
                struct Cleanup {
                    mpz_ptr value;
                    ~Cleanup() { mpz_clear(value); }
                } cleanup{truncated};
                const auto shift = static_cast<mp_bitcnt_t>(gap);
                mpz_tdiv_q_2exp(truncated, b._significand, shift);
                magnitude = mpz_cmp(a._significand, truncated);
                if (magnitude == 0 && !mpz_divisible_2exp_p(b._significand, shift)) {
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
                mpz_t truncated;
                mpz_init(truncated);
                struct Cleanup {
                    mpz_ptr value;
                    ~Cleanup() { mpz_clear(value); }
                } cleanup{truncated};
                const auto shift = static_cast<mp_bitcnt_t>(gap);
                mpz_tdiv_q_2exp(truncated, a._significand, shift);
                magnitude = mpz_cmp(truncated, b._significand);
                if (magnitude == 0 && !mpz_divisible_2exp_p(a._significand, shift)) {
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

    mpz_t scaled, root, remainder, threshold;
    mpz_inits(scaled, root, remainder, threshold, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c, d;
        ~Cleanup() { mpz_clears(a, b, c, d, nullptr); }
    } cleanup{scaled, root, remainder, threshold};

    mpz_mul_2exp(scaled, value._significand,
                 FloatAP::shift_amount(static_cast<FloatAP::Exponent>(shift)));
    mpz_sqrtrem(root, remainder, scaled);

    bool increment = false;
    if (mpz_size(remainder) != 0) {
        switch (mode.value()) {
            case FloatAP::Rounding::ToNearest:
                // sqrt(scaled) > root+1/2 iff 4*scaled > (2*root+1)^2.
                // Equality is impossible: the left side is divisible
                // by four, while the right side is odd.
                mpz_mul_2exp(threshold, root, 1);
                mpz_add_ui(threshold, threshold, 1);
                mpz_mul(threshold, threshold, threshold);
                mpz_mul_2exp(scaled, scaled, 2);
                increment = mpz_cmp(scaled, threshold) > 0;
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
        mpz_add_ui(root, root, 1);
    }
    if (mpz_sizeinbase(root, 2) > bits) {
        mpz_tdiv_q_2exp(root, root, 1);
        result_exponent = FloatAP::checked_add(result_exponent, 1);
    }

    FloatAP result(bits);
    mpz_set(result._significand, root);
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
        mpz_pow_ui(result._significand, value._significand, degree);
        FloatAP::round_scaled_integer(result, product_scale, mode);
        return result;
    }

    // Rare normalization-cancellation/scale-overflow fallback:
    // retain the original exact rational computation, including its
    // exception and rounding behavior.
    mpq_t exact;
    mpq_init(exact);
    struct RationalCleanup {
        mpq_ptr value;
        ~RationalCleanup() { mpq_clear(value); }
    } cleanup{exact};
    value.to_rational(exact);

    mpz_t numerator, denominator;
    mpz_inits(numerator, denominator, nullptr);
    struct IntegerCleanup {
        mpz_ptr a, b;
        ~IntegerCleanup() { mpz_clears(a, b, nullptr); }
    } ints{numerator, denominator};

    mpz_pow_ui(numerator, mpq_numref(exact), degree);
    mpz_pow_ui(denominator, mpq_denref(exact), degree);
    return FloatAP::from_ratio(numerator, denominator, 0, bits, mode);
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
    const Rational rational = static_cast<Rational>(x);
    return Fraction(rational.get_mpq());
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
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        throw std::overflow_error("FloatAP exponential GMP shift overflow");
    }
    const mp_bitcnt_t shift = static_cast<mp_bitcnt_t>(fixed_bits);

    mpz_t numerator, argument, term, sum_lower, sum_upper, threshold;
    mpz_inits(numerator, argument, term, sum_lower, sum_upper, threshold, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c, d, e, f;
        ~Cleanup() { mpz_clears(a, b, c, d, e, f, nullptr); }
    } cleanup{numerator, argument, term, sum_lower, sum_upper, threshold};

    mpz_mul_2exp(numerator, mpq_numref(magnitude.raw()), shift);
    mpz_fdiv_q(argument, numerator, mpq_denref(magnitude.raw()));

    mpz_set_ui(term, 1ul);
    mpz_mul_2exp(term, term, shift);
    mpz_set(sum_lower, term);

    mpz_set_ui(threshold, 1ul);
    mpz_mul_2exp(threshold, threshold, static_cast<mp_bitcnt_t>(guard));

    for (unsigned long n = 1ul; n <= max_terms; ++n) {
        // floor(floor(y/Q)/n) == floor(y/(Q*n)) for nonnegative y.
        mpz_mul(term, term, argument);
        mpz_fdiv_q_2exp(term, term, shift);
        mpz_fdiv_q_ui(term, term, n);
        mpz_add(sum_lower, sum_lower, term);

        if (mpz_cmp(term, threshold) <= 0 || n == max_terms) {
            // 6*n bounds the accumulated downward-rounding error;
            // t_n+6 bounds the complete positive omitted tail.
            mpz_set(sum_upper, sum_lower);
            mpz_add(sum_upper, sum_upper, term);
            mpz_set_ui(numerator, n);
            mpz_add_ui(numerator, numerator, 1ul);
            mpz_mul_ui(numerator, numerator, 6ul);
            mpz_add(sum_upper, sum_upper, numerator);

            const Interval positive{Fraction(sum_lower, fixed_bits),
                                    Fraction(sum_upper, fixed_bits)};
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
    // Reduce more aggressively when greater precision makes the
    // exponential Taylor series expensive. Fixed-point reconstruction
    // keeps the size of subsequent square intermediates bounded.
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

    // The series endpoints are certified positive rationals. Convert them
    // to directed B-bit fixed-point integers once, then square without
    // letting their exact rational denominators double at each step.
    const unsigned long fixed_bits = work + halvings + 64ul;
    if (static_cast<std::uintmax_t>(fixed_bits) >
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        throw std::overflow_error("FloatAP exponential reconstruction GMP shift overflow");
    }
    const mp_bitcnt_t shift = static_cast<mp_bitcnt_t>(fixed_bits);
    mpz_t numerator, lower, upper;
    mpz_inits(numerator, lower, upper, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c;
        ~Cleanup() { mpz_clears(a, b, c, nullptr); }
    } cleanup{numerator, lower, upper};

    mpz_mul_2exp(numerator, mpq_numref(interval.lower.raw()), shift);
    mpz_fdiv_q(lower, numerator, mpq_denref(interval.lower.raw()));
    mpz_mul_2exp(numerator, mpq_numref(interval.upper.raw()), shift);
    mpz_cdiv_q(upper, numerator, mpq_denref(interval.upper.raw()));

    // Squaring is increasing on positive reals: floor lower, ceil upper.
    // With B >= work + halvings + 64, accumulated rounding errors are
    // tiny, but correctness relies only on outward bounds, not that size.
    for (unsigned i = 0u; i < halvings; ++i) {
        mpz_mul(lower, lower, lower);
        mpz_fdiv_q_2exp(lower, lower, shift);
        mpz_mul(upper, upper, upper);
        mpz_cdiv_q_2exp(upper, upper, shift);
    }
    return {Fraction(lower, fixed_bits), Fraction(upper, fixed_bits)};
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
    const unsigned long guard = 32ul + static_cast<unsigned long>(std::bit_width(terms));
    if (work > std::numeric_limits<unsigned long>::max() - guard) {
        throw std::overflow_error("FloatAP logarithm work precision overflow");
    }
    const unsigned long fixed_bits = work + guard;

    mpz_t numerator, zlo, zhi, square_lo, square_hi, power_lo, power_hi;
    mpz_t sum_lo, sum_hi, quotient, tail;
    mpz_inits(numerator, zlo, zhi, square_lo, square_hi, power_lo, power_hi,
              sum_lo, sum_hi, quotient, tail, nullptr);
    struct Cleanup {
        mpz_ptr a,b,c,d,e,f,g,h,i,j,k;
        ~Cleanup() { mpz_clears(a,b,c,d,e,f,g,h,i,j,k,nullptr); }
    } cleanup{numerator,zlo,zhi,square_lo,square_hi,power_lo,power_hi,
              sum_lo,sum_hi,quotient,tail};

    const Fraction magnitude = absolute(z);
    mpz_mul_2exp(numerator, mpq_numref(magnitude.raw()), fixed_bits);
    mpz_fdiv_q(zlo, numerator, mpq_denref(magnitude.raw()));
    mpz_cdiv_q(zhi, numerator, mpq_denref(magnitude.raw()));

    mpz_mul(square_lo, zlo, zlo);
    mpz_fdiv_q_2exp(square_lo, square_lo, fixed_bits);
    mpz_mul(square_hi, zhi, zhi);
    mpz_cdiv_q_2exp(square_hi, square_hi, fixed_bits);
    mpz_set(power_lo, zlo);
    mpz_set(power_hi, zhi);

    // For |z| <= 1/3 the remaining doubled series after N terms is
    // at most (9/4)*|z|^(2N+1). The upper fixed-point odd power is
    // outward rounded, so three times that upper bound is a safe
    // integer enclosure of the entire tail. Stop as soon as the
    // *actual* argument converges, not after a worst-case term count.
    mpz_t threshold;
    mpz_init_set_ui(threshold, 1ul);
    struct ThresholdCleanup {
        mpz_ptr value;
        ~ThresholdCleanup() { mpz_clear(value); }
    } threshold_cleanup{threshold};
    mpz_mul_2exp(threshold, threshold, fixed_bits - work);

    for (unsigned long n = 0; n < terms; ++n) {
        const unsigned long divisor = 2ul * n + 1ul;
        mpz_fdiv_q_ui(quotient, power_lo, divisor);
        mpz_add(sum_lo, sum_lo, quotient);
        mpz_cdiv_q_ui(quotient, power_hi, divisor);
        mpz_add(sum_hi, sum_hi, quotient);

        mpz_mul(power_lo, power_lo, square_lo);
        mpz_fdiv_q_2exp(power_lo, power_lo, fixed_bits);
        mpz_mul(power_hi, power_hi, square_hi);
        mpz_cdiv_q_2exp(power_hi, power_hi, fixed_bits);

        mpz_mul_ui(tail, power_hi, 3ul);
        if (mpz_cmp(tail, threshold) <= 0 || n + 1ul == terms) {
            break;
        }
    }

    // The partial sum contains only positive magnitudes; negate and
    // reverse the interval for an originally negative argument.
    mpz_mul_2exp(sum_lo, sum_lo, 1);
    mpz_mul_2exp(sum_hi, sum_hi, 1);
    mpz_add(sum_hi, sum_hi, tail);

    const Interval positive{Fraction(sum_lo, fixed_bits),
                            Fraction(sum_hi, fixed_bits)};
    return z.sign() < 0
        ? Interval{Fraction{} - positive.upper, Fraction{} - positive.lower}
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

Interval exp_outward_dyadic(Interval const& bounds,
                            unsigned long denominator_bits) {
    if (static_cast<std::uintmax_t>(denominator_bits) >
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        throw std::overflow_error("FloatAP cached exponential shift overflow");
    }
    const mp_bitcnt_t scale = static_cast<mp_bitcnt_t>(denominator_bits);
    mpz_t scaled_numerator, endpoint_lower, endpoint_upper;
    mpz_inits(scaled_numerator, endpoint_lower, endpoint_upper, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c;
        ~Cleanup() { mpz_clears(a, b, c, nullptr); }
    } cleanup{scaled_numerator, endpoint_lower, endpoint_upper};

    mpz_mul_2exp(scaled_numerator, mpq_numref(bounds.lower.raw()), scale);
    mpz_fdiv_q(endpoint_lower, scaled_numerator,
               mpq_denref(bounds.lower.raw()));
    mpz_mul_2exp(scaled_numerator, mpq_numref(bounds.upper.raw()), scale);
    mpz_cdiv_q(endpoint_upper, scaled_numerator,
               mpq_denref(bounds.upper.raw()));
    return {Fraction(endpoint_lower, denominator_bits),
            Fraction(endpoint_upper, denominator_bits)};
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
        const FloatAP lower = FloatAP::from_rational(interval.lower.raw(), bits, mode);
        const FloatAP upper = FloatAP::from_rational(interval.upper.raw(), bits, mode);
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

    const Fraction magnitude = absolute(x);
    const unsigned long max_terms = work / 2ul + 48ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(max_terms));
    if (work > std::numeric_limits<unsigned long>::max() - guard) {
        throw std::overflow_error("FloatAP sine/cosine work precision overflow");
    }
    const unsigned long fixed_bits = work + guard;
    if (static_cast<std::uintmax_t>(fixed_bits) >
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        throw std::overflow_error("FloatAP sine/cosine GMP shift overflow");
    }
    const mp_bitcnt_t shift = static_cast<mp_bitcnt_t>(fixed_bits);

    mpz_t scaled, argument, square, term, sum, lower, upper;
    mpz_t error, tail, scratch;
    mpz_inits(scaled, argument, square, term, sum, lower, upper,
              error, tail, scratch, nullptr);
    struct Cleanup {
        mpz_ptr a,b,c,d,e,f,g,h,i,j;
        ~Cleanup() { mpz_clears(a,b,c,d,e,f,g,h,i,j,nullptr); }
    } cleanup{scaled,argument,square,term,sum,lower,upper,
              error,tail,scratch};

    mpz_set_ui(scaled, 1ul);
    mpz_mul_2exp(scaled, scaled, shift);
    mpz_mul_2exp(scratch, mpq_numref(magnitude.raw()), shift);
    mpz_fdiv_q(argument, scratch, mpq_denref(magnitude.raw()));
    mpz_mul(square, argument, argument);
    mpz_fdiv_q_2exp(square, square, shift);

    mpz_t threshold;
    mpz_init_set_ui(threshold, 1ul);
    struct ThresholdCleanup {
        mpz_ptr value;
        ~ThresholdCleanup() { mpz_clear(value); }
    } threshold_cleanup{threshold};
    mpz_mul_2exp(threshold, threshold, static_cast<mp_bitcnt_t>(guard));

    // Let Q=2^B and a=floor(Q*|x|). The fixed-point square
    // b=floor(a*a/Q) underapproximates Q*x*x by less than 3.
    // For each nonnegative Taylor term T_n, the recurrence
    // t_(n+1)=floor(t_n*b/(Q*d_n)) satisfies 0<=T_n-t_n<8:
    // the error obeys e_(n+1)<(e_n+3)/d_n+1, d_n>=2.
    // This applies to both sin and cos for every |x|<=1.
    const auto series = [&](bool sine) -> Interval {
        if (sine) {
            mpz_set(term, argument);
        } else {
            mpz_set(term, scaled);
        }
        mpz_set_ui(sum, 0ul);

        for (unsigned long n = 0; n < max_terms; ++n) {
            if (n % 2ul == 0ul) {
                mpz_add(sum, sum, term);
            } else {
                mpz_sub(sum, sum, term);
            }

            if (n > (std::numeric_limits<unsigned long>::max() - 3ul) / 2ul) {
                throw std::overflow_error("FloatAP sine/cosine series index");
            }
            const unsigned long twice_n = 2ul * n;
            const unsigned long a = sine ? twice_n + 2ul : twice_n + 1ul;
            const unsigned long b = sine ? twice_n + 3ul : twice_n + 2ul;
            if (a > std::numeric_limits<unsigned long>::max() / b) {
                throw std::overflow_error("FloatAP sine/cosine factorial overflow");
            }
            const unsigned long divisor = a * b;

            mpz_mul(term, term, square);
            mpz_fdiv_q_2exp(term, term, shift);
            mpz_fdiv_q_ui(term, term, divisor);

            if (mpz_cmp(term, threshold) <= 0 || n + 1ul == max_terms) {
                // The n+1 included terms contribute at most 8*(n+1)
                // total units of downward-rounding error. The first
                // omitted exact term is <= t_(n+1)+8. Its alternating
                // sign identifies which endpoint needs the tail.
                mpz_set_ui(error, n);
                mpz_add_ui(error, error, 1ul);
                mpz_mul_ui(error, error, 8ul);
                mpz_set(lower, sum);
                mpz_sub(lower, lower, error);
                mpz_set(upper, sum);
                mpz_add(upper, upper, error);
                mpz_set(tail, term);
                mpz_add_ui(tail, tail, 8ul);
                if (n % 2ul == 0ul) {
                    mpz_sub(lower, lower, tail);
                } else {
                    mpz_add(upper, upper, tail);
                }

                return {Fraction(lower, fixed_bits),
                        Fraction(upper, fixed_bits)};
            }
        }
        throw std::logic_error("FloatAP sine/cosine series did not terminate");
    };

    const bool want_sine = outputs != TrigOutputs::cosine;
    Interval sine = want_sine ? series(true) : point(Fraction{});
    Interval cosine = point(Fraction{});

    if (outputs == TrigOutputs::cosine) {
        cosine = series(false);
    } else if (outputs == TrigOutputs::both) {
        // On |x|<=1, the sine magnitude is in [0,1] and cosine is
        // nonnegative. The positive-sine interval is still represented
        // in the shared `lower` and `upper` scaled integer endpoints.
        if (mpz_sgn(lower) < 0) {
            mpz_set_ui(lower, 0ul);
        }
        if (mpz_cmp(upper, scaled) > 0) {
            mpz_set(upper, scaled);
        }
        if (mpz_cmp(lower, upper) > 0) {
            throw std::logic_error("FloatAP invalid certified sine bounds");
        }

        // For 0<=L<=sin(|x|)*Q<=U<=8Q/9, let
        // f(s)=sqrt(Q^2-s^2). Because |f'(s)|=s/f(s)<=8/sqrt(17)<2
        // on [0,8Q/9], f(L)-f(U)<=2(U-L).
        // Thus just ONE downward integer square root of Q^2-U^2
        // certifies both cosine endpoints: C_lo=floor(f(U)),
        // C_hi=C_lo+1+2(U-L). An alternative exact two-root path
        // handles any overly wide sine interval outside the bound.
        mpz_mul_ui(scratch, upper, 9ul);
        mpz_mul_ui(term, scaled, 8ul);
        if (mpz_cmp(scratch, term) <= 0) {
            mpz_mul(scratch, scaled, scaled);
            mpz_mul(term, upper, upper);
            mpz_sub(term, scratch, term);
            mpz_sqrt(sum, term);
            mpz_sub(tail, upper, lower);
            mpz_mul_2exp(tail, tail, 1u);
            mpz_add_ui(tail, tail, 1u);
            mpz_add(scratch, sum, tail);
            mpz_set(lower, sum);
            mpz_set(upper, scratch);
        } else {
            // Conservative general fallback: two exact integer roots.
            mpz_mul(scratch, scaled, scaled);
            mpz_mul(term, lower, lower);
            mpz_sub(term, scratch, term);
            mpz_sqrtrem(sum, tail, term);
            if (mpz_sgn(tail) != 0) {
                mpz_add_ui(sum, sum, 1ul);
            }
            mpz_mul(term, upper, upper);
            mpz_sub(term, scratch, term);
            mpz_sqrt(lower, term);
            mpz_set(upper, sum);
        }

        cosine = {Fraction(lower, fixed_bits),
                  Fraction(upper, fixed_bits)};
    }

    if (want_sine && x.sign() < 0) {
        sine = negative(sine);
    }
    return {std::move(sine), std::move(cosine)};
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

    const Fraction magnitude = absolute(x);
    const unsigned long max_terms = work / 2ul + 4ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(max_terms));
    if (work > std::numeric_limits<unsigned long>::max() - guard) {
        throw std::overflow_error("FloatAP arctangent work precision overflow");
    }
    const unsigned long fixed_bits = work + guard;
    if (static_cast<std::uintmax_t>(fixed_bits) >
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        throw std::overflow_error("FloatAP arctangent GMP shift overflow");
    }
    const mp_bitcnt_t shift = static_cast<mp_bitcnt_t>(fixed_bits);

    mpz_t scaled, argument, square, power, sum, term, tail, error;
    mpz_inits(scaled, argument, square, power, sum, term, tail, error, nullptr);
    struct Cleanup {
        mpz_ptr a,b,c,d,e,f,g,h;
        ~Cleanup() { mpz_clears(a,b,c,d,e,f,g,h,nullptr); }
    } cleanup{scaled,argument,square,power,sum,term,tail,error};

    // Let Q=2^fixed_bits, a=floor(Q*|x|) and b=floor(a*a/Q).
    // For |x|<=1/2, 0<=Q*x*x-b<2 and a<=Q*x.
    mpz_mul_2exp(scaled, mpq_numref(magnitude.raw()), shift);
    mpz_fdiv_q(argument, scaled, mpq_denref(magnitude.raw()));
    mpz_mul(square, argument, argument);
    mpz_fdiv_q_2exp(square, square, shift);
    mpz_set(power, argument);
    mpz_set_ui(sum, 0ul);

    // p_n=floor(p_(n-1)*b/Q), p_0=a. If P_n=Q*|x|^(2n+1),
    // e_n=P_n-p_n obeys e_0<1 and e_(n+1)<e_n/4+2,
    // so 0<=e_n<4 for every n. Consequently each downward-rounded
    // Taylor term p_n/(2n+1) differs from the exact scaled term by
    // less than 5 integer units, independently of n. The first
    // omitted alternating term is bounded by p_(n+1)+4 units.
    // One recurrence therefore certifies BOTH endpoints: no need for
    // separate upward and downward powers, squares and divisions.
    mpz_t limit;
    mpz_init_set_ui(limit, 1ul);
    struct LimitCleanup {
        mpz_ptr value;
        ~LimitCleanup() { mpz_clear(value); }
    } limit_cleanup{limit};
    mpz_mul_2exp(limit, limit, static_cast<mp_bitcnt_t>(guard));

    for (unsigned long n = 0; n < max_terms; ++n) {
        const unsigned long divisor = 2ul * n + 1ul;
        mpz_fdiv_q_ui(term, power, divisor);
        if (n % 2ul == 0ul) mpz_add(sum, sum, term);
        else mpz_sub(sum, sum, term);

        mpz_mul(power, power, square);
        mpz_fdiv_q_2exp(power, power, shift);

        if (mpz_cmp(power, limit) <= 0 || n + 1ul == max_terms) {
            // All n+1 term errors together have magnitude <5*(n+1).
            // The next term is the alternating-series remainder bound.
            mpz_set_ui(error, n);
            mpz_add_ui(error, error, 1ul);
            mpz_mul_ui(error, error, 5ul);
            mpz_t lower, upper;
            mpz_inits(lower, upper, nullptr);
            struct EndpointCleanup {
                mpz_ptr a, b;
                ~EndpointCleanup() { mpz_clears(a, b, nullptr); }
            } endpoint_cleanup{lower, upper};
            mpz_sub(lower, sum, error);
            mpz_add(upper, sum, error);
            mpz_add_ui(tail, power, 4ul);
            if (n % 2ul == 0ul) mpz_sub(lower, lower, tail);
            else mpz_add(upper, upper, tail);

            const Interval positive{Fraction(lower, fixed_bits),
                                    Fraction(upper, fixed_bits)};
            return x.sign() < 0 ? negative(positive) : positive;
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
    // At modest precision the extra exact-rational operations cost more
    // than direct Taylor. For tangent, evaluating BOTH local intervals
    // leaves large rational denominators for the final division and
    // loses much more than it saves. Keep the dyadic Taylor/one-root
    // path for that case at every precision.
    if (halvings == 0u && work >= 128ul &&
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
    if (work > std::numeric_limits<mp_bitcnt_t>::max() / 2u) {
        throw std::overflow_error("FloatAP radical precision too large");
    }

    mpz_t scaled, quotient, root, squared;
    mpz_inits(scaled, quotient, root, squared, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c, d;
        ~Cleanup() { mpz_clears(a, b, c, d, nullptr); }
    } cleanup{scaled, quotient, root, squared};

    mpz_mul_2exp(scaled, mpq_numref(x.raw()), 2u * work);
    mpz_fdiv_q(quotient, scaled, mpq_denref(x.raw()));
    mpz_sqrt(root, quotient);
    const Fraction lower(root, work);

    mpz_mul(squared, root, root);
    mpz_mul(squared, squared, mpq_denref(x.raw()));
    if (mpz_cmp(squared, scaled) == 0) {
        return point(lower);
    }

    mpz_add_ui(root, root, 1u);
    return {lower, Fraction(root, work)};
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
    if (static_cast<std::uintmax_t>(bits) >
        static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        throw std::overflow_error("FloatAP dyadic tangent interval shift overflow");
    }
    const mp_bitcnt_t shift = static_cast<mp_bitcnt_t>(bits);
    mpz_t scaled, lower, upper;
    mpz_inits(scaled, lower, upper, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c;
        ~Cleanup() { mpz_clears(a, b, c, nullptr); }
    } cleanup{scaled, lower, upper};

    mpz_mul_2exp(scaled, mpq_numref(input.lower.raw()), shift);
    mpz_fdiv_q(lower, scaled, mpq_denref(input.lower.raw()));
    mpz_mul_2exp(scaled, mpq_numref(input.upper.raw()), shift);
    mpz_cdiv_q(upper, scaled, mpq_denref(input.upper.raw()));
    return {Fraction(lower, bits), Fraction(upper, bits)};
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
        case Function::arctangent:
            return atan_interval(x, work);
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
        const FloatAP low = FloatAP::from_rational(bounds.lower.raw(), bits, mode);
        const FloatAP high = FloatAP::from_rational(bounds.upper.raw(), bits, mode);
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
        return log_interval(input, work);
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
