/***************************************************************************
 *            numeric/floatap.hpp
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

/*! \file numeric/floatap.hpp
 *  \brief Native arbitrary-precision floating-point type and operations.
 */

#ifndef ARIADNE_FLOATAP_HPP
#define ARIADNE_FLOATAP_HPP

#include "numeric/bits.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/sign.hpp"
#include <cstddef>
#include <ostream>
#include <gmp.h>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <type_traits>

namespace Ariadne {

class Integer;
class Rational;
class Dyadic;
class ExactDouble;
class TwoExp;
struct RoundDownward;
struct RoundUpward;
struct RoundToNearest;
struct RoundTowardZero;
struct RoundApproximately;
struct CurrentRoundingMode;


// A size_t bit count is not necessarily representable by the signed
// exponent. Keep the narrowing check testable without allocating an
// infeasibly large GMP integer.
inline std::int64_t floatap_exponent_from_bit_length(std::size_t bit_count) {
    if constexpr (std::numeric_limits<std::size_t>::digits >
                  std::numeric_limits<std::int64_t>::digits) {
        if (bit_count > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
            throw std::overflow_error("FloatAP operand too large");
        }
    }
    return static_cast<std::int64_t>(bit_count);
}

// The precision tag is deliberately independent of MPFR types.
class ArbitraryPrecision {
    std::size_t _bits;
  public:
    explicit ArbitraryPrecision(Bits count) : _bits(static_cast<unsigned long>(count)) {
    }
    explicit ArbitraryPrecision(std::size_t count) : _bits(count) {
    }
    std::size_t bits() const noexcept {
        return _bits;
    }
    friend OutputStream& operator<<(OutputStream& os, ArbitraryPrecision p) {
        return os << "ap(" << p.bits() << ")";
    }
    friend OutputStream& repr(OutputStream& os, ArbitraryPrecision p) {
        return os << "ArbitraryPrecision(" << p.bits() << ")";
    }
    friend bool operator==(ArbitraryPrecision a, ArbitraryPrecision b) noexcept {
        return a._bits == b._bits;
    }
    friend bool operator<(ArbitraryPrecision a, ArbitraryPrecision b) noexcept {
        return a._bits < b._bits;
    }
    friend bool operator<=(ArbitraryPrecision a, ArbitraryPrecision b) noexcept {
        return a._bits <= b._bits;
    }
    friend bool operator>(ArbitraryPrecision a, ArbitraryPrecision b) noexcept {
        return b < a;
    }
    friend bool operator>=(ArbitraryPrecision a, ArbitraryPrecision b) noexcept {
        return b <= a;
    }
};

// The public native arbitrary-precision type. The API is growing; elementary
// arithmetic is implemented, while transcendental functions are pending.
using FloatAP = Float<ArbitraryPrecision>;

// MPFR-independent floating-point representation and arithmetic.
// A finite nonzero value is (-1)^negative * significand * 2^exponent,
// with precisely precision significant bits. GMP is temporary limb storage.
template<> class Float<ArbitraryPrecision> {
  public:
    // A value type rather than an enum: Ariadne's rounding tags convert
    // directly to it, without making the GMP kernel depend on MPFR.
    class Rounding {
      public:
        enum Value { ToNearest, Downward, Upward, TowardZero };
      private:
        Value _value;
      public:
        constexpr Rounding(Value value) noexcept : _value(value) { }
        Rounding(RoundToNearest const&) noexcept;
        Rounding(RoundDownward const&) noexcept;
        Rounding(RoundUpward const&) noexcept;
        Rounding(RoundTowardZero const&) noexcept;
        Rounding(RoundApproximately const&) noexcept;
        Rounding(CurrentRoundingMode const&) noexcept;
        constexpr Value value() const noexcept { return _value; }
        friend constexpr bool operator==(Rounding const&, Rounding const&) = default;
    };
    using Paradigm = ExactTag;
    using NumericType = FloatAP;
    using PrecisionType = ArbitraryPrecision;
    using RoundingModeType = Rounding;
    static const RoundingModeType ROUND_TO_NEAREST;
    static const RoundingModeType ROUND_DOWNWARD;
    static const RoundingModeType ROUND_UPWARD;
    static const RoundingModeType ROUND_TOWARD_ZERO;
    enum class Kind {
        Finite, Infinity, NaN
    };
    using Exponent = std::int64_t;

  private:
    mpz_t _significand;
    unsigned long _precision;
    Exponent _exponent = 0;
    Kind _kind = Kind::Finite;
    bool _negative = false;

    static Exponent checked_add(Exponent a, Exponent b) {
        constexpr auto max = std::numeric_limits<Exponent>::max();
        constexpr auto min = std::numeric_limits<Exponent>::min();
        if ((b > 0 && a > max - b) || (b < 0 && a < min - b)) {
            throw std::overflow_error("FloatAP exponent overflow");
        }
        return a + b;
    }
    static Exponent checked_sub(Exponent a, Exponent b) {
        if (b == std::numeric_limits<Exponent>::min()) {
            return checked_add(checked_add(a, std::numeric_limits<Exponent>::max()), 1);
        }
        return checked_add(a,- b);
    }
    static mp_bitcnt_t shift_amount(Exponent n) {
        // All callers pass a nonnegative shift, enforced by their branches
        // and by checked_sub(). Only narrower GMP bit-count types need a
        // representability check (notably possible on Windows).
        if constexpr (std::numeric_limits<mp_bitcnt_t>::digits <
                      std::numeric_limits<Exponent>::digits) {
            if (static_cast<std::uint64_t>(n) > std::numeric_limits<mp_bitcnt_t>::max()) {
                throw std::overflow_error("FloatAP shift out of range");
            }
        }
        return static_cast<mp_bitcnt_t>(n);
    }
    static Exponent bit_length(mpz_srcptr x) {
        return floatap_exponent_from_bit_length(mpz_sizeinbase(x, 2));
    }
    void set_zero(bool negative) {
        _kind = Kind::Finite;
        _negative = negative;
        _exponent = 0;
        mpz_set_ui(_significand, 0);
    }

  public:
    explicit Float(unsigned long precision) : _precision(precision) {
        if (precision < 2 || precision > static_cast<std::uint64_t>(std::numeric_limits<Exponent>::max())) {
            throw std::invalid_argument("FloatAP precision must be >=2 and fit in exponent range");
        }
        mpz_init(_significand);
    }
    explicit Float(PrecisionType precision)
    : Float(static_cast<unsigned long>(precision.bits())) {
    }
    Float(std::int64_t n, PrecisionType precision);
    Float(ExactDouble const& n, PrecisionType precision);
    Float(Dyadic const& n, PrecisionType precision);
    explicit Float(String const& value, PrecisionType precision);
    Float(TwoExp const& n, PrecisionType precision);
    Float(FloatDP const& n, PrecisionType precision);
    Float(FloatAP const& n, PrecisionType precision);
    Float(FloatAP const& n, RoundingModeType mode, PrecisionType precision);
    Float(FloatDP const& n, RoundingModeType mode, PrecisionType precision);
    Float(ExactDouble const& n, RoundingModeType mode, PrecisionType precision);
    Float(Dyadic const& n, RoundingModeType mode, PrecisionType precision);
    Float(Integer const& n, PrecisionType precision)
        : Float(n, ROUND_TO_NEAREST, precision) { }

    Float(Integer const& n, RoundingModeType mode, PrecisionType precision);
    Float(Rational const& n, RoundingModeType mode, PrecisionType precision);
    explicit Float(Rational const& n, PrecisionType precision)
    : Float(n, ROUND_TO_NEAREST, precision) {
    }
    explicit operator Rational() const;
    explicit operator Dyadic() const;
    double get_d() const;
    ~Float() {
        mpz_clear(_significand);
    }
    Float(FloatAP const& x)
    : _precision(x._precision), _exponent(x._exponent), _kind(x._kind), _negative(x._negative) {
        mpz_init_set(_significand, x._significand);
    }
    Float(FloatAP && x) noexcept
    : _precision(x._precision), _exponent(x._exponent), _kind(x._kind), _negative(x._negative) {
        mpz_init(_significand);
        mpz_swap(_significand, x._significand);
    }
    FloatAP& operator=(FloatAP x) noexcept {
        swap(x);
        return *this;
    }
    FloatAP& operator=(Integer const& n);
    FloatAP& operator=(Dyadic const& n);
    FloatAP& operator=(TwoExp const& n);
    FloatAP& operator=(ExactDouble const& n);
    FloatAP& assign_unsigned(std::uint64_t n);
    template<BuiltinIntegral N> FloatAP& operator=(N n) {
        if constexpr (std::is_signed_v<N>) {
            return *this = FloatAP(static_cast<std::int64_t>(n), precision());
        } else {
            return assign_unsigned(static_cast<std::uint64_t>(n));
        }
    }
    void swap(FloatAP& x) noexcept {
        mpz_swap(_significand, x._significand);
        std::swap(_precision, x._precision);
        std::swap(_exponent, x._exponent);
        std::swap(_kind, x._kind);
        std::swap(_negative, x._negative);
    }

    static RoundingModeType get_rounding_mode() noexcept;
    static void set_rounding_mode(RoundingModeType mode);
    static void set_rounding_to_nearest() { set_rounding_mode(ROUND_TO_NEAREST); }
    static void set_rounding_downward() { set_rounding_mode(ROUND_DOWNWARD); }
    static void set_rounding_upward() { set_rounding_mode(ROUND_UPWARD); }
    static void set_rounding_toward_zero() { set_rounding_mode(ROUND_TOWARD_ZERO); }
    static PrecisionType get_default_precision();
    static void set_default_precision(PrecisionType precision);
    void set_precision(PrecisionType precision);
    FloatAP const& raw() const noexcept { return *this; }
    FloatAP& raw() noexcept { return *this; }

    PrecisionType precision() const noexcept {
        return PrecisionType(_precision);
    }
    Exponent exponent() const noexcept {
        return _exponent;
    }
    mpz_srcptr significand() const noexcept {
        return _significand;
    }
    Kind kind() const noexcept {
        return _kind;
    }
    bool negative() const noexcept {
        return _negative;
    }
    bool is_zero() const noexcept {
        return _kind == Kind::Finite && mpz_size(_significand) == 0;
    }
    bool is_finite() const noexcept {
        return _kind == Kind::Finite;
    }
    bool is_inf() const noexcept {
        return _kind == Kind::Infinity;
    }
    bool is_nan() const noexcept {
        return _kind == Kind::NaN;
    }

    static FloatAP zero(unsigned long p, bool negative = false) {
        FloatAP z(p);
        z.set_zero(negative);
        return z;
    }
    static FloatAP infinity(unsigned long p, bool negative = false) {
        FloatAP z(p);
        z._kind = Kind::Infinity;
        z._negative = negative;
        return z;
    }
    static FloatAP nan(unsigned long p) {
        FloatAP z(p);
        z._kind = Kind::NaN;
        return z;
    }
    static FloatAP inf(Sign sign, PrecisionType p) {
        return sign == Sign::ZERO ? nan(p) : infinity(static_cast<unsigned long>(p.bits()), sign == Sign::NEGATIVE);
    }
    static FloatAP inf(PrecisionType p) {
        return infinity(static_cast<unsigned long>(p.bits()));
    }
    static FloatAP nan(PrecisionType p) {
        return nan(static_cast<unsigned long>(p.bits()));
    }

    friend bool is_nan(FloatAP const& x) noexcept {
        return x.is_nan();
    }
    friend bool is_inf(FloatAP const& x) noexcept {
        return x.is_inf();
    }
    friend bool is_finite(FloatAP const& x) noexcept {
        return x.is_finite();
    }
    friend bool is_zero(FloatAP const& x) noexcept {
        return x.is_zero();
    }
    friend Comparison cmp(FloatAP const& a, FloatAP const& b);
    friend Comparison cmp(FloatAP const& a, Rational const& b);
    friend Comparison cmp(Rational const& a, FloatAP const& b);
    friend Comparison cmp(FloatAP const& a, FloatDP const& b);
    friend Comparison cmp(FloatDP const& a, FloatAP const& b);
    friend bool operator==(FloatAP const& a, FloatDP const& b);
    friend bool operator!=(FloatAP const& a, FloatDP const& b);
    friend bool operator<(FloatAP const& a, FloatDP const& b);
    friend bool operator>(FloatAP const& a, FloatDP const& b);
    friend bool operator<=(FloatAP const& a, FloatDP const& b);
    friend bool operator>=(FloatAP const& a, FloatDP const& b);
    friend bool operator==(FloatDP const& a, FloatAP const& b);
    friend bool operator!=(FloatDP const& a, FloatAP const& b);
    friend bool operator<(FloatDP const& a, FloatAP const& b);
    friend bool operator>(FloatDP const& a, FloatAP const& b);
    friend bool operator<=(FloatDP const& a, FloatAP const& b);
    friend bool operator>=(FloatDP const& a, FloatAP const& b);
    friend Integer cast_integer(FloatAP const& value);
    friend bool operator==(FloatAP const& a, Rational const& b) {
        return cmp(a, b) == Comparison::EQUAL;
    }
    friend bool operator!=(FloatAP const& a, Rational const& b) {
        return cmp(a, b) != Comparison::EQUAL;
    }
    friend bool operator<(FloatAP const& a, Rational const& b) {
        return cmp(a, b) == Comparison::LESS;
    }
    friend bool operator>(FloatAP const& a, Rational const& b) {
        return cmp(a, b) == Comparison::GREATER;
    }
    friend bool operator<=(FloatAP const& a, Rational const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::LESS || c == Comparison::EQUAL;
    }
    friend bool operator>=(FloatAP const& a, Rational const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::GREATER || c == Comparison::EQUAL;
    }
    friend bool operator==(Rational const& a, FloatAP const& b) {
        return cmp(a, b) == Comparison::EQUAL;
    }
    friend bool operator!=(Rational const& a, FloatAP const& b) {
        return cmp(a, b) != Comparison::EQUAL;
    }
    friend bool operator<(Rational const& a, FloatAP const& b) {
        return cmp(a, b) == Comparison::LESS;
    }
    friend bool operator>(Rational const& a, FloatAP const& b) {
        return cmp(a, b) == Comparison::GREATER;
    }
    friend bool operator<=(Rational const& a, FloatAP const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::LESS || c == Comparison::EQUAL;
    }
    friend bool operator>=(Rational const& a, FloatAP const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::GREATER || c == Comparison::EQUAL;
    }
    friend OutputStream& operator<<(OutputStream& os, FloatAP const& x);
    friend InputStream& operator>>(InputStream& is, FloatAP& x);
    friend bool same(FloatAP const& a, FloatAP const& b) { return a == b; }
    friend Positive<FloatAP> abs(FloatAP const& x);
    friend FloatAP abs(RoundingModeType, FloatAP const& x);
    friend FloatAP operator*(TwoExp const& t, FloatAP const& x);
    friend FloatAP operator*(FloatAP const& x, TwoExp const& t);
    friend FloatAP operator/(FloatAP const& x, TwoExp const& t);
    friend FloatAP& operator*=(FloatAP& x, TwoExp const& t);
    friend FloatAP& operator/=(FloatAP& x, TwoExp const& t);
    friend Bounds<FloatAP> operator+(FloatAP const& a, FloatAP const& b);
    friend Bounds<FloatAP> operator-(FloatAP const& a, FloatAP const& b);
    friend Bounds<FloatAP> operator*(FloatAP const& a, FloatAP const& b);
    friend Bounds<FloatAP> operator/(FloatAP const& a, FloatAP const& b);
    friend Bounds<FloatAP> add(FloatAP const& a, FloatAP const& b);
    friend Bounds<FloatAP> sub(FloatAP const& a, FloatAP const& b);
    friend Bounds<FloatAP> mul(FloatAP const& a, FloatAP const& b);
    friend Bounds<FloatAP> div(FloatAP const& a, FloatAP const& b);
    friend Bounds<FloatAP> sqr(FloatAP const& x);
    friend Bounds<FloatAP> rec(FloatAP const& x);
    friend Bounds<FloatAP> pow(FloatAP const& x, Nat n);
    friend Bounds<FloatAP> pow(FloatAP const& x, Int n);
    friend FloatAP pow(RoundingModeType mode, FloatAP const& x, Int exponent);
    friend FloatAP operator+(FloatAP const& x) { return x; }
    friend FloatAP operator-(FloatAP const& x) { return negate(x); }
    friend FloatAP nul(FloatAP const& x) { return zero(x._precision); }
    friend FloatAP pos(FloatAP const& x) { return x; }
    friend FloatAP neg(FloatAP const& x) { return negate(x); }
    friend FloatAP shft(FloatAP const& x, Int amount);
    friend FloatAP hlf(FloatAP const& x) { return shft(x, Int(-1)); }
    friend FloatAP max(FloatAP const& a, FloatAP const& b) { return a < b ? b : a; }
    friend FloatAP min(FloatAP const& a, FloatAP const& b) { return a < b ? a : b; }
    friend FloatAP sqr(RoundingModeType mode, FloatAP const& x) {
        return FloatAP::mul(x, x, x._precision, mode);
    }
    friend FloatAP rec(RoundingModeType mode, FloatAP const& x) {
        const FloatAP one(std::int64_t{1}, x.precision());
        return FloatAP::div(one, x, x._precision, mode);
    }
    friend FloatAP med(RoundingModeType mode, FloatAP const& a, FloatAP const& b) {
        const unsigned long precision = a._precision < b._precision
            ? a._precision : b._precision;
        return hlf(FloatAP::add(a, b, precision, mode));
    }
    friend FloatAP rad(RoundingModeType mode, FloatAP const& a, FloatAP const& b) {
        const unsigned long precision = a._precision < b._precision
            ? a._precision : b._precision;
        return hlf(FloatAP::sub(b, a, precision, mode));
    }
    friend bool operator==(FloatAP const& a, FloatAP const& b) {
        return cmp(a, b) == Comparison::EQUAL;
    }
    friend bool operator!=(FloatAP const& a, FloatAP const& b) {
        return !(a == b);
    }
    friend bool operator<(FloatAP const& a, FloatAP const& b) {
        return cmp(a, b) == Comparison::LESS;
    }
    friend bool operator>(FloatAP const& a, FloatAP const& b) {
        return cmp(a, b) == Comparison::GREATER;
    }
    friend bool operator<=(FloatAP const& a, FloatAP const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::LESS || c == Comparison::EQUAL;
    }
    friend bool operator>=(FloatAP const& a, FloatAP const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::GREATER || c == Comparison::EQUAL;
    }

    friend FloatAP add(RoundingModeType mode, FloatAP const& a, FloatAP const& b) {
        return add(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    friend FloatAP sub(RoundingModeType mode, FloatAP const& a, FloatAP const& b) {
        return sub(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    friend FloatAP mul(RoundingModeType mode, FloatAP const& a, FloatAP const& b) {
        return mul(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    friend FloatAP div(RoundingModeType mode, FloatAP const& a, FloatAP const& b) {
        return div(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    friend FloatAP fma(RoundingModeType mode, FloatAP const& a,
    FloatAP const& b, FloatAP const& c) {
        const auto p = a._precision < b._precision ? a._precision : b._precision;
        return fma(a, b, c, p < c._precision ? p : c._precision, mode);
    }

    // Square root is rounded from an exact integer-square-root enclosure.
    // Integer powers are computed exactly, then rounded only once.
    friend FloatAP sqrt(RoundingModeType mode, FloatAP const& x);
    friend FloatAP pow(RoundingModeType mode, FloatAP const& x, Nat exponent);
    friend FloatAP exp(RoundingModeType mode, FloatAP const& x);
    friend FloatAP log(RoundingModeType mode, FloatAP const& x);
    friend FloatAP sin(RoundingModeType mode, FloatAP const& x);
    friend FloatAP cos(RoundingModeType mode, FloatAP const& x);
    friend FloatAP tan(RoundingModeType mode, FloatAP const& x);
    friend FloatAP atan(RoundingModeType mode, FloatAP const& x);
    friend FloatAP asin(RoundingModeType mode, FloatAP const& x);
    friend FloatAP acos(RoundingModeType mode, FloatAP const& x);
    static FloatAP pi(RoundingModeType mode, PrecisionType precision);

    // Round signed numerator/positive denominator * 2^scale with one rounding.
    static FloatAP from_ratio(mpz_srcptr numerator, mpz_srcptr denominator,
    Exponent scale, unsigned long precision,
    Rounding mode, bool negative_zero = false);
    static FloatAP from_integer(mpz_srcptr n, unsigned long precision, Rounding mode);
    static FloatAP from_rational(mpq_srcptr value, unsigned long precision, Rounding mode);
    // Exact conversion, not an approximation to double.
    void to_rational(mpq_ptr out) const;
    static FloatAP negate(FloatAP x);
    static FloatAP add(FloatAP const& a, FloatAP const& b,
    unsigned long precision, Rounding mode);
    static FloatAP sub(FloatAP const& a, FloatAP const& b,
    unsigned long precision, Rounding mode);
    static FloatAP mul(FloatAP const& a, FloatAP const& b,
    unsigned long precision, Rounding mode);
    static FloatAP div(FloatAP const& a, FloatAP const& b,
    unsigned long precision, Rounding mode);
    // Exact product followed by a single rounding after addition.
    static FloatAP fma(FloatAP const& a, FloatAP const& b,
    FloatAP const& c, unsigned long precision, Rounding mode);
};

inline constexpr FloatAP::RoundingModeType FloatAP::ROUND_TO_NEAREST =
    FloatAP::Rounding::ToNearest;
inline constexpr FloatAP::RoundingModeType FloatAP::ROUND_DOWNWARD =
    FloatAP::Rounding::Downward;
inline constexpr FloatAP::RoundingModeType FloatAP::ROUND_UPWARD =
    FloatAP::Rounding::Upward;
inline constexpr FloatAP::RoundingModeType FloatAP::ROUND_TOWARD_ZERO =
    FloatAP::Rounding::TowardZero;

}
#endif

