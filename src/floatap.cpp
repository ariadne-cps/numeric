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
#include <cmath>
#include <iomanip>
#include <sstream>
#include <istream>
#include <ostream>
#include <string>
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

}
// namespace Ariadne

