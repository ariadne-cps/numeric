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

FloatAP FloatAP::from_ratio(mpz_srcptr numerator, mpz_srcptr denominator,
Exponent scale, unsigned long precision, Rounding mode, bool negative_zero) {
    FloatAP result(precision);
    if (mpz_sgn(denominator) == 0) throw std::domain_error("FloatAP zero denominator");
    if (mpz_sgn(numerator) == 0) {
        result.set_zero(negative_zero);
        return result;
    }
    result._negative =(mpz_sgn(numerator) < 0) !=(mpz_sgn(denominator) < 0);
    mpz_t n, d, q, r, w;
    mpz_inits(n, d, q, r, w, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c, d, e;
        ~Cleanup() {
            mpz_clears(a, b, c, d, e, nullptr);
        }
    }
    cleanup { n, d, q, r, w };
    mpz_abs(n, numerator);
    mpz_abs(d, denominator);
    Exponent k = checked_sub(bit_length(n), bit_length(d));
    if (k >= 0) {
        mpz_mul_2exp(w, d, shift_amount(k));
        if (mpz_cmp(n, w) < 0) k = checked_sub(k, 1);
    } else {
        mpz_mul_2exp(w, n, shift_amount(checked_sub(0, k)));
        if (mpz_cmp(w, d) < 0) k = checked_sub(k, 1);
    }
    const Exponent shift = checked_sub(static_cast<Exponent>(precision - 1), k);
    if (shift >= 0) {
        mpz_mul_2exp(w, n, shift_amount(shift));
        mpz_tdiv_qr(q, r, w, d);
        mpz_set(w, d);
    } else {
        mpz_mul_2exp(w, d, shift_amount(checked_sub(0, shift)));
        mpz_tdiv_qr(q, r, n, w);
    }
    bool increment = false;
    if (mpz_size(r) != 0) {
        switch (mode.value()) {
            case Rounding::ToNearest:
            mpz_mul_2exp(r, r, 1);
            increment = mpz_cmp(r, w) > 0 || (mpz_cmp(r, w) == 0 && mpz_odd_p(q));
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
    if (increment) mpz_add_ui(q, q, 1);
    result._exponent = checked_add(scale, checked_sub(k, static_cast<Exponent>(precision - 1)));
    if (bit_length(q) > static_cast<Exponent>(precision)) {
        mpz_tdiv_q_2exp(q, q, 1);
        result._exponent = checked_add(result._exponent, 1);
    }
    mpz_set(result._significand, q);
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
    return from_ratio(mpq_numref(r), mpq_denref(r), 0, p, mode);
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
    if (a.is_nan() || b.is_nan()) return nan(p);
    if (a.is_inf() || b.is_inf()) {
        if (a.is_inf() && b.is_inf() && a._negative != b._negative) return nan(p);
        return infinity(p, a.is_inf() ? a._negative : b._negative);
    }
    if (a.is_zero() && b.is_zero()) return zero(p, a._negative == b._negative ? a._negative : mode == Rounding::Downward);
    if (a.is_zero()) return FloatAP(b, mode, PrecisionType(p));
    if (b.is_zero()) return FloatAP(a, mode, PrecisionType(p));
    // For equal-precision normalized operands separated by more than
    // p+2 binary positions, the small addend cannot affect nearest-even.
    // A directed result only differs by one representable step. This
    // avoids allocating an impossible exponent-gap-sized GMP integer.
    if (a._precision == p && b._precision == p && a._exponent != b._exponent) {
        FloatAP const& dominant = a._exponent > b._exponent ? a : b;
        FloatAP const& smaller = a._exponent > b._exponent ? b : a;
        const std::uint64_t gap =
            static_cast<std::uint64_t>(dominant._exponent) -
            static_cast<std::uint64_t>(smaller._exponent);
        if (gap > static_cast<std::uint64_t>(p) + 2u) {
            if (mode == Rounding::Upward && !smaller._negative) return next(up, dominant);
            if (mode == Rounding::Downward && smaller._negative) return next(down, dominant);
            if (mode == Rounding::TowardZero && dominant._negative != smaller._negative) {
                return dominant._negative ? next(up, dominant) : next(down, dominant);
            }
            return dominant;
        }
    }
    const Exponent low = a._exponent < b._exponent ? a._exponent : b._exponent;
    mpz_t x, y, sum, one;
    mpz_inits(x, y, sum, one, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c, d;
        ~Cleanup() {
            mpz_clears(a, b, c, d, nullptr);
        }
    }
    cleanup { x, y, sum, one };
    mpz_mul_2exp(x, a._significand, shift_amount(checked_sub(a._exponent, low)));
    mpz_mul_2exp(y, b._significand, shift_amount(checked_sub(b._exponent, low)));
    if (a._negative) mpz_neg(x, x);
    if (b._negative) mpz_neg(y, y);
    mpz_add(sum, x, y);
    mpz_set_ui(one, 1);
    return from_ratio(sum, one, low, p, mode, mode == Rounding::Downward);
}

FloatAP FloatAP::sub(FloatAP const& a, FloatAP const& b,
unsigned long p, Rounding mode) {
    return add(a, negate(b), p, mode);
}

FloatAP FloatAP::mul(FloatAP const& a, FloatAP const& b,
unsigned long p, Rounding mode) {
    if (a.is_nan() || b.is_nan() || (a.is_inf() && b.is_zero()) || (a.is_zero() && b.is_inf())) return nan(p);
    const bool neg = a._negative != b._negative;
    if (a.is_inf() || b.is_inf()) return infinity(p, neg);
    if (a.is_zero() || b.is_zero()) return zero(p, neg);
    mpz_t n, one;
    mpz_inits(n, one, nullptr);
    struct Cleanup {
        mpz_ptr a, b;
        ~Cleanup() {
            mpz_clears(a, b, nullptr);
        }
    }
    cleanup { n, one };
    mpz_mul(n, a._significand, b._significand);
    if (neg) mpz_neg(n, n);
    mpz_set_ui(one, 1);
    return from_ratio(n, one, checked_add(a._exponent, b._exponent), p, mode);
}

FloatAP FloatAP::div(FloatAP const& a, FloatAP const& b,
unsigned long p, Rounding mode) {
    if (a.is_nan() || b.is_nan() || (a.is_zero() && b.is_zero()) || (a.is_inf() && b.is_inf())) return nan(p);
    const bool neg = a._negative != b._negative;
    if (a.is_inf() || b.is_zero()) return infinity(p, neg);
    if (a.is_zero() || b.is_inf()) return zero(p, neg);
    mpz_t n;
    mpz_init_set(n, a._significand);
    struct Cleanup {
        mpz_ptr z;
        ~Cleanup() {
            mpz_clear(z);
        }
    }
    cleanup { n };
    if (neg) mpz_neg(n, n);
    return from_ratio(n, b._significand, checked_sub(a._exponent, b._exponent), p, mode);
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
    : Float(static_cast<Rational>(value), mode, precision) {
    if (value.is_nan()) {
        *this = nan(precision);
    } else if (value.is_inf()) {
        *this = infinity(static_cast<unsigned long>(precision.bits()), value.negative());
    } else if (value.is_zero()) {
        *this = zero(static_cast<unsigned long>(precision.bits()), value.negative());
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
    return static_cast<Rational>(* this).get_d();
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

// Evaluate exp on |x| <= 1/2 with outward-rounded fixed-point GMP
// integers. Taylor terms for exp(|x|) are positive; the omitted tail is
// at most twice the first omitted term. For x < 0, take the reciprocal
// of the positive enclosure, exchanging its endpoints.
Interval exp_series(Fraction const& x, unsigned long work) {
    if (x.sign() == 0) {
        return {Fraction(1), Fraction(1)};
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

    mpz_t scaled, argument_lo, argument_hi, term_lo, term_hi;
    mpz_t sum_lo, sum_hi, next_hi, remainder, threshold;
    mpz_inits(scaled, argument_lo, argument_hi, term_lo, term_hi,
              sum_lo, sum_hi, next_hi, remainder, threshold, nullptr);
    struct Cleanup {
        mpz_ptr a,b,c,d,e,f,g,h,i,j;
        ~Cleanup() { mpz_clears(a,b,c,d,e,f,g,h,i,j,nullptr); }
    } cleanup{scaled,argument_lo,argument_hi,term_lo,term_hi,
              sum_lo,sum_hi,next_hi,remainder,threshold};

    mpz_mul_2exp(scaled, mpq_numref(magnitude.raw()), shift);
    mpz_fdiv_q(argument_lo, scaled, mpq_denref(magnitude.raw()));
    mpz_cdiv_q(argument_hi, scaled, mpq_denref(magnitude.raw()));

    mpz_set_ui(term_lo, 1ul);
    mpz_mul_2exp(term_lo, term_lo, shift);
    mpz_set(term_hi, term_lo);
    mpz_set(sum_lo, term_lo);
    mpz_set(sum_hi, term_hi);

    mpz_set_ui(threshold, 1ul);
    mpz_mul_2exp(threshold, threshold, static_cast<mp_bitcnt_t>(guard));

    for (unsigned long n = 1ul; n <= max_terms; ++n) {
        // Division by 2^B and then by n, with floor/ceiling on the
        // respective endpoints, encloses the exact nth Taylor term.
        mpz_mul(term_lo, term_lo, argument_lo);
        mpz_fdiv_q_2exp(term_lo, term_lo, shift);
        mpz_fdiv_q_ui(term_lo, term_lo, n);

        mpz_mul(term_hi, term_hi, argument_hi);
        mpz_cdiv_q_2exp(term_hi, term_hi, shift);
        mpz_cdiv_q_ui(term_hi, term_hi, n);

        mpz_add(sum_lo, sum_lo, term_lo);
        mpz_add(sum_hi, sum_hi, term_hi);

        // The omitted positive tail is bounded by twice term (n+1).
        mpz_mul(next_hi, term_hi, argument_hi);
        mpz_cdiv_q_2exp(next_hi, next_hi, shift);
        mpz_cdiv_q_ui(next_hi, next_hi, n + 1ul);
        mpz_mul_2exp(remainder, next_hi, 1ul);

        if (mpz_cmp(remainder, threshold) <= 0 || n == max_terms) {
            mpz_add(sum_hi, sum_hi, remainder);
            const Interval positive{Fraction(sum_lo, fixed_bits),
                                    Fraction(sum_hi, fixed_bits)};
            if (x.sign() > 0) {
                return positive;
            }
            return {Fraction(1) / positive.upper,
                    Fraction(1) / positive.lower};
        }
    }
    throw std::logic_error("FloatAP exponential series did not terminate");
}

Interval exp_interval(Fraction x, unsigned long work) {
    unsigned halvings = 0;
    while (compare(absolute(x), Fraction(1) / Fraction(2)) > 0) {
        if (++halvings > 2048) {
            throw std::overflow_error("FloatAP exp argument range");
        }
        x = x / Fraction(2);
    }
    Interval interval = exp_series(x, work + halvings + 8);
    for (unsigned i = 0; i < halvings; ++i) {
        // exp(x) is positive, so endpoint multiplication is monotone.
        interval = {interval.lower * interval.lower,
                    interval.upper * interval.upper};
    }
    return interval;
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
    }

    // Multiply by two, then add the certified bound for the omitted tail.
    mpz_mul_2exp(sum_lo, sum_lo, 1);
    mpz_mul_2exp(sum_hi, sum_hi, 1);
    mpz_set_ui(tail, 1);
    mpz_mul_2exp(tail, tail, fixed_bits - work);
    mpz_add(sum_hi, sum_hi, tail);

    const Interval positive{Fraction(sum_lo, fixed_bits),
                            Fraction(sum_hi, fixed_bits)};
    return z.sign() < 0
        ? Interval{Fraction{} - positive.upper, Fraction{} - positive.lower}
        : positive;
}

Interval log_interval(Fraction x, unsigned long work) {
    long power = 0;
    while (compare(x, Fraction(2)) > 0) {
        if (++power > 2048) {
            throw std::overflow_error("FloatAP log argument range");
        }
        x = x / Fraction(2);
    }
    while (compare(x, Fraction(1) / Fraction(2)) < 0) {
        if (--power < -2048) {
            throw std::overflow_error("FloatAP log argument range");
        }
        x = x * Fraction(2);
    }
    const Fraction z = (x - Fraction(1)) / (x + Fraction(1));
    Interval result = log_series(z, work);
    if (power != 0) {
        const Interval ln2 = log_series(Fraction(1) / Fraction(3), work);
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
    const Interval reciprocal = {
        Fraction(1L) / denominator.upper,
        Fraction(1L) / denominator.lower
    };
    return multiply(numerator, reciprocal);
}

// For |x| <= 1, the positive magnitudes of both sine and cosine
// Taylor terms decrease monotonically. Evaluate both series together
// using outward-rounded fixed-point GMP integers: no mpq_t intermediates
// are allocated in the loop. The first omitted alternating term bounds
// the exact remainder on one side of the partial sum.
std::pair<Interval, Interval> small_sine_cosine(Fraction const& x,
                                               unsigned long work) {
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

    mpz_t scaled, argument_lo, argument_hi, square_lo, square_hi;
    mpz_t sine_term_lo, sine_term_hi, cosine_term_lo, cosine_term_hi;
    mpz_t sine_sum_lo, sine_sum_hi, cosine_sum_lo, cosine_sum_hi, threshold;
    mpz_inits(scaled, argument_lo, argument_hi, square_lo, square_hi,
              sine_term_lo, sine_term_hi, cosine_term_lo, cosine_term_hi,
              sine_sum_lo, sine_sum_hi, cosine_sum_lo, cosine_sum_hi,
              threshold, nullptr);
    struct Cleanup {
        mpz_ptr a,b,c,d,e,f,g,h,i,j,k,l,m,n;
        ~Cleanup() { mpz_clears(a,b,c,d,e,f,g,h,i,j,k,l,m,n,nullptr); }
    } cleanup{scaled,argument_lo,argument_hi,square_lo,square_hi,
              sine_term_lo,sine_term_hi,cosine_term_lo,cosine_term_hi,
              sine_sum_lo,sine_sum_hi,cosine_sum_lo,cosine_sum_hi,threshold};

    mpz_mul_2exp(scaled, mpq_numref(magnitude.raw()), shift);
    mpz_fdiv_q(argument_lo, scaled, mpq_denref(magnitude.raw()));
    mpz_cdiv_q(argument_hi, scaled, mpq_denref(magnitude.raw()));
    mpz_mul(square_lo, argument_lo, argument_lo);
    mpz_fdiv_q_2exp(square_lo, square_lo, shift);
    mpz_mul(square_hi, argument_hi, argument_hi);
    mpz_cdiv_q_2exp(square_hi, square_hi, shift);

    mpz_set(sine_term_lo, argument_lo);
    mpz_set(sine_term_hi, argument_hi);
    mpz_set_ui(cosine_term_lo, 1ul);
    mpz_mul_2exp(cosine_term_lo, cosine_term_lo, shift);
    mpz_set(cosine_term_hi, cosine_term_lo);

    mpz_set_ui(threshold, 1ul);
    mpz_mul_2exp(threshold, threshold, static_cast<mp_bitcnt_t>(guard));

    for (unsigned long n = 0; n < max_terms; ++n) {
        if (n % 2ul == 0ul) {
            mpz_add(sine_sum_lo, sine_sum_lo, sine_term_lo);
            mpz_add(sine_sum_hi, sine_sum_hi, sine_term_hi);
            mpz_add(cosine_sum_lo, cosine_sum_lo, cosine_term_lo);
            mpz_add(cosine_sum_hi, cosine_sum_hi, cosine_term_hi);
        } else {
            mpz_sub(sine_sum_lo, sine_sum_lo, sine_term_hi);
            mpz_sub(sine_sum_hi, sine_sum_hi, sine_term_lo);
            mpz_sub(cosine_sum_lo, cosine_sum_lo, cosine_term_hi);
            mpz_sub(cosine_sum_hi, cosine_sum_hi, cosine_term_lo);
        }

        if (n > (std::numeric_limits<unsigned long>::max() - 3ul) / 2ul) {
            throw std::overflow_error("FloatAP sine/cosine series index");
        }
        const unsigned long twice_n = 2ul * n;
        const unsigned long a = twice_n + 1ul;
        const unsigned long b = twice_n + 2ul;
        const unsigned long c = twice_n + 3ul;
        if (a > std::numeric_limits<unsigned long>::max() / b ||
            b > std::numeric_limits<unsigned long>::max() / c) {
            throw std::overflow_error("FloatAP sine/cosine factorial overflow");
        }
        const unsigned long cosine_divisor = a * b;
        const unsigned long sine_divisor = b * c;

        mpz_mul(sine_term_lo, sine_term_lo, square_lo);
        mpz_fdiv_q_2exp(sine_term_lo, sine_term_lo, shift);
        mpz_fdiv_q_ui(sine_term_lo, sine_term_lo, sine_divisor);
        mpz_mul(sine_term_hi, sine_term_hi, square_hi);
        mpz_cdiv_q_2exp(sine_term_hi, sine_term_hi, shift);
        mpz_cdiv_q_ui(sine_term_hi, sine_term_hi, sine_divisor);

        mpz_mul(cosine_term_lo, cosine_term_lo, square_lo);
        mpz_fdiv_q_2exp(cosine_term_lo, cosine_term_lo, shift);
        mpz_fdiv_q_ui(cosine_term_lo, cosine_term_lo, cosine_divisor);
        mpz_mul(cosine_term_hi, cosine_term_hi, square_hi);
        mpz_cdiv_q_2exp(cosine_term_hi, cosine_term_hi, shift);
        mpz_cdiv_q_ui(cosine_term_hi, cosine_term_hi, cosine_divisor);

        if ((mpz_cmp(sine_term_hi, threshold) <= 0 &&
             mpz_cmp(cosine_term_hi, threshold) <= 0) ||
            n + 1ul == max_terms) {
            // If the last included term was positive, the omitted
            // alternating tail is negative, and conversely.
            if (n % 2ul == 0ul) {
                mpz_sub(sine_sum_lo, sine_sum_lo, sine_term_hi);
                mpz_sub(cosine_sum_lo, cosine_sum_lo, cosine_term_hi);
            } else {
                mpz_add(sine_sum_hi, sine_sum_hi, sine_term_hi);
                mpz_add(cosine_sum_hi, cosine_sum_hi, cosine_term_hi);
            }

            Interval sine{Fraction(sine_sum_lo, fixed_bits),
                          Fraction(sine_sum_hi, fixed_bits)};
            Interval cosine{Fraction(cosine_sum_lo, fixed_bits),
                            Fraction(cosine_sum_hi, fixed_bits)};
            if (x.sign() < 0) {
                sine = negative(sine);
            }
            return {std::move(sine), std::move(cosine)};
        }
    }
    throw std::logic_error("FloatAP sine/cosine series did not terminate");
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

    mpz_t scaled, xlo, xhi, square_lo, square_hi, power_lo, power_hi;
    mpz_t sum_lo, sum_hi, term, tail, limit;
    mpz_inits(scaled, xlo, xhi, square_lo, square_hi, power_lo, power_hi,
              sum_lo, sum_hi, term, tail, limit, nullptr);
    struct Cleanup {
        mpz_ptr a,b,c,d,e,f,g,h,i,j,k,l;
        ~Cleanup() { mpz_clears(a,b,c,d,e,f,g,h,i,j,k,l,nullptr); }
    } cleanup{scaled,xlo,xhi,square_lo,square_hi,power_lo,power_hi,
              sum_lo,sum_hi,term,tail,limit};

    mpz_mul_2exp(scaled, mpq_numref(magnitude.raw()),
                 static_cast<mp_bitcnt_t>(fixed_bits));
    mpz_fdiv_q(xlo, scaled, mpq_denref(magnitude.raw()));
    mpz_cdiv_q(xhi, scaled, mpq_denref(magnitude.raw()));

    mpz_mul(square_lo, xlo, xlo);
    mpz_fdiv_q_2exp(square_lo, square_lo,
                    static_cast<mp_bitcnt_t>(fixed_bits));
    mpz_mul(square_hi, xhi, xhi);
    mpz_cdiv_q_2exp(square_hi, square_hi,
                    static_cast<mp_bitcnt_t>(fixed_bits));
    mpz_set(power_lo, xlo);
    mpz_set(power_hi, xhi);

    mpz_set_ui(limit, 1ul);
    mpz_mul_2exp(limit, limit, static_cast<mp_bitcnt_t>(guard));

    for (unsigned long n = 0; n < max_terms; ++n) {
        const unsigned long divisor = 2ul * n + 1ul;
        if (n % 2ul == 0ul) {
            mpz_fdiv_q_ui(term, power_lo, divisor);
            mpz_add(sum_lo, sum_lo, term);
            mpz_cdiv_q_ui(term, power_hi, divisor);
            mpz_add(sum_hi, sum_hi, term);
        } else {
            mpz_cdiv_q_ui(term, power_hi, divisor);
            mpz_sub(sum_lo, sum_lo, term);
            mpz_fdiv_q_ui(term, power_lo, divisor);
            mpz_sub(sum_hi, sum_hi, term);
        }

        mpz_mul(power_lo, power_lo, square_lo);
        mpz_fdiv_q_2exp(power_lo, power_lo,
                        static_cast<mp_bitcnt_t>(fixed_bits));
        mpz_mul(power_hi, power_hi, square_hi);
        mpz_cdiv_q_2exp(power_hi, power_hi,
                        static_cast<mp_bitcnt_t>(fixed_bits));

        // The exact remainder has the sign of term n+1 and is bounded
        // by |x|^(2n+3)/(2n+3). power_hi encloses its numerator.
        mpz_cdiv_q_ui(tail, power_hi, 2ul * n + 3ul);
        if (mpz_cmp(tail, limit) <= 0 || n + 1ul == max_terms) {
            if (n % 2ul == 0ul) {
                mpz_sub(sum_lo, sum_lo, tail);
            } else {
                mpz_add(sum_hi, sum_hi, tail);
            }
            const Interval positive{Fraction(sum_lo, fixed_bits),
                                    Fraction(sum_hi, fixed_bits)};
            return x.sign() < 0
                ? negative(positive) : positive;
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

Interval atan_interval(Fraction const& x, unsigned long work) {
    const Fraction half = Fraction(1L) / Fraction(2L);
    if (compare(absolute(x), half) <= 0) {
        return small_atan(x, work);
    }
    if (x.sign() < 0) {
        return negative(atan_interval(Fraction{} - x, work));
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
                                                   unsigned long work) {
    unsigned halvings = 0;
    while (compare(absolute(x), Fraction(1L)) > 0) {
        if (++halvings > 4096u) {
            throw std::overflow_error("FloatAP trigonometric argument too large");
        }
        x = x / Fraction(2L);
    }

    auto results = small_sine_cosine(x, work + halvings + 16u);
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
    const Interval left = atan_interval(argument.lower, work);
    const Interval right = atan_interval(argument.upper, work);
    return {Fraction(2L) * left.lower, Fraction(2L) * right.upper};
}

enum class Function {
    sine, cosine, tangent, arctangent, arcsine, arccosine,
    pi, half_pi, negative_half_pi
};

Interval evaluate(Function function, Fraction const& x, unsigned long work) {
    switch (function) {
        case Function::sine:
            return sine_cosine_intervals(x, work).first;
        case Function::cosine:
            return sine_cosine_intervals(x, work).second;
        case Function::tangent: {
            const auto [sine, cosine] = sine_cosine_intervals(x, work);
            return divide(sine, cosine);
        }
        case Function::arctangent:
            return atan_interval(x, work);
        case Function::arcsine:
            return asin_interval(x, work);
        case Function::arccosine:
            return subtract(divide(pi_interval(work), point(Fraction(2L))),
                            asin_interval(x, work));
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
        return exp_interval(input, work);
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
