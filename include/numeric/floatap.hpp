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
#include "numeric/number.decl.hpp"
#include "numeric/number.hpp"
#include "numeric/rounding.hpp"
#include "numeric/sign.hpp"
#include "numeric/big_unsigned.hpp"
#include <bit>
#include <cstddef>
#include <ostream>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <type_traits>

namespace Ariadne {

class Integer;
class Rational;
class Decimal;
class Dyadic;
class ExactDouble;
class TwoExp;
namespace Detail {
    enum class FloatTranscendental {
        sine,cosine,tangent,arctangent,arcsine,arccosine,
        pi,half_pi,negative_half_pi,exponential,logarithm,hyperbolic_tangent
    };
    struct FloatAPShortTranscendentals;
}
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
    friend ArbitraryPrecision min(ArbitraryPrecision a, ArbitraryPrecision b) noexcept {
        return a._bits <= b._bits ? a : b;
    }
    friend ArbitraryPrecision max(ArbitraryPrecision a, ArbitraryPrecision b) noexcept {
        return a._bits >= b._bits ? a : b;
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
// with precisely precision significant bits and a native BigUnsigned
// significand. GMP remains only at external legacy/interval boundaries.
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
        Rounding(Ariadne::Rounding const&) noexcept;
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
    friend struct Detail::FloatAPShortTranscendentals;
    // Exact nonnegative magnitude; sign and special-value kind are separate.
    static FloatAP from_native_ratio(Integer const& numerator,
                                     Integer const& denominator,
                                     Exponent scale,unsigned long precision,
                                     Rounding mode,bool negative_zero=false,
                                     bool invert_sign=false);
    static FloatAP from_unsigned_ratio(BigUnsigned const& numerator,
                                       BigUnsigned const& denominator,
                                       Exponent scale,unsigned long precision,
                                       Rounding mode,bool negative,
                                       bool negative_zero=false);
    BigUnsigned _significand;
    unsigned long _precision;
    Exponent _exponent = 0;
    Kind _kind = Kind::Finite;
    bool _negative = false;
    double get_d_general() const;

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
    static std::size_t shift_amount(Exponent n) {
        // Every native shift must fit size_t, including on 32-bit targets.
        // Call sites provide a nonnegative exponent or a checked distance.
        if (n<0 || static_cast<std::uintmax_t>(n)>
                   std::numeric_limits<std::size_t>::max()) {
            throw std::overflow_error("FloatAP shift out of range");
        }
        return static_cast<std::size_t>(n);
    }
    static Exponent bit_length(BigUnsigned const& x) {
        return floatap_exponent_from_bit_length(x.bit_length());
    }
    // Normalize and round an already-owned nonzero positive significand.
    static void round_scaled_integer(FloatAP& result, Exponent scale, Rounding mode);
    static void round_short(FloatAP& result,BigUnsigned::LimbView magnitude,
                            Exponent scale,Rounding mode);
    void set_zero(bool negative) {
        _kind = Kind::Finite;
        _negative = negative;
        _exponent = 0;
        _significand=BigUnsigned{};
    }

  public:
    explicit Float(unsigned long precision) : _precision(precision) {
        if (precision < 2 || precision > static_cast<std::uint64_t>(std::numeric_limits<Exponent>::max())) {
            throw std::invalid_argument("FloatAP precision must be >=2 and fit in exponent range");
        }
    }
    explicit Float(PrecisionType precision)
    : Float(static_cast<unsigned long>(precision.bits())) {
    }
    Float(std::int64_t n, PrecisionType precision);
    Float(ExactDouble const& n, PrecisionType precision);
    Float(Dyadic const& n, PrecisionType precision);
    explicit Float(std::string const& value, PrecisionType precision);
    Float(TwoExp const& n, PrecisionType precision);
    Float(FloatDP const& n, PrecisionType precision);
    Float(FloatAP const& n, PrecisionType precision);
    Float(FloatAP const& n, RoundingModeType mode, PrecisionType precision);
    Float(FloatDP const& n, RoundingModeType mode, PrecisionType precision);
    Float(ExactDouble const& n, RoundingModeType mode, PrecisionType precision);
    Float(Dyadic const& n, RoundingModeType mode, PrecisionType precision);
    Float(Decimal const& n, RoundingModeType mode, PrecisionType precision);
    Float(Integer const& n, PrecisionType precision)
        : Float(n, ROUND_TO_NEAREST, precision) { }

    Float(Integer const& n, RoundingModeType mode, PrecisionType precision);
    template<BuiltinIntegral N>
    Float(N n, RoundingModeType mode, PrecisionType precision)
        : Float(Integer(n), mode, precision) { }
    Float(Rational const& n, RoundingModeType mode, PrecisionType precision);
    //! Round an exact signed native numerator/denominator once.
    //! No canonical Rational intermediate (which would incur GCD).
    //! A zero denominator throws std::domain_error.
    static FloatAP from_ratio(Integer const& numerator,
                              Integer const& denominator,
                              RoundingModeType mode,PrecisionType precision);
    explicit Float(Rational const& n, PrecisionType precision)
    : Float(n, ROUND_TO_NEAREST, precision) {
    }
    explicit operator Rational() const;
    explicit operator Dyadic() const;
    double get_d() const {
        if constexpr (sizeof(double)==sizeof(std::uint64_t) &&
                      std::numeric_limits<double>::is_iec559 &&
                      std::numeric_limits<double>::digits==53) {
            // A normal binary64 with <=53 significant bits is exact in
            // every rounding mode. Keep this frequent conversion inline.
            if (_kind==Kind::Finite && _precision<=53u) {
                if (_significand.is_zero()) return _negative ? -0.0 : 0.0;
                if (_exponent>=-1022-static_cast<Exponent>(_precision-1u) &&
                    _exponent<=1023-static_cast<Exponent>(_precision-1u)) {
                    const auto highest=_exponent+static_cast<Exponent>(_precision-1u);
                    const std::uint64_t fraction=(_significand.limb(0u)<<(53u-_precision)) &
                        ((std::uint64_t{1}<<52u)-1u);
                    return std::bit_cast<double>((static_cast<std::uint64_t>(_negative)<<63u) |
                        (static_cast<std::uint64_t>(highest+1023)<<52u) | fraction);
                }
            }
        }
        return get_d_general();
    }
    std::string literal() const;
    std::string literal(RoundingModeType mode) const;
    static Nat output_places;
    static void set_output_places(Nat places);
    friend std::string print(FloatAP const& x, DecimalPlaces places, RoundingModeType mode);
    friend std::string print(FloatAP const& x, DecimalPrecision figures, RoundingModeType mode);
    friend OutputStream& write(OutputStream& out, FloatAP const& x, DecimalPlaces places, RoundingModeType mode);
    friend OutputStream& write(OutputStream& out, FloatAP const& x, DecimalPrecision figures, RoundingModeType mode);
    friend OutputStream& repr(OutputStream& out, FloatAP const& x, RoundingModeType mode);
    friend OutputStream& repr(OutputStream& out, FloatAP const& x);
    ~Float() = default;
    Float(FloatAP const& x)
    : _significand(x._significand), _precision(x._precision),
      _exponent(x._exponent), _kind(x._kind), _negative(x._negative) { }
    Float(FloatAP&& x) noexcept
    : _significand(std::move(x._significand)), _precision(x._precision),
      _exponent(x._exponent), _kind(x._kind), _negative(x._negative) { }
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
        _significand.swap(x._significand);
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
    operator ExactNumber() const;
    Ball<FloatAP,FloatAP> pm(Error<FloatAP> const& error) const;
    Ball<FloatAP,FloatDP> pm(Error<FloatDP> const& error) const;

    PrecisionType precision() const noexcept {
        return PrecisionType(_precision);
    }
    Exponent exponent() const noexcept {
        return _exponent;
    }
    BigUnsigned const& significand() const noexcept {
        return _significand;
    }
    Kind kind() const noexcept {
        return _kind;
    }
    bool negative() const noexcept {
        return _negative;
    }
    bool is_zero() const noexcept {
        return _kind == Kind::Finite && _significand.is_zero();
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
    // Unit in the last place at 1, with the selected significand width.
    static FloatAP eps(PrecisionType p) {
        FloatAP value(std::int64_t{1}, p);
        value._exponent = checked_sub(value._exponent,
                                      static_cast<Exponent>(p.bits() - 1));
        return value;
    }
    static FloatAP nan(PrecisionType p) {
        return nan(static_cast<unsigned long>(p.bits()));
    }
    // Exponent range comes from the signed native exponent representation.
    // No subnormal significands are supported by this kernel.
    static FloatAP min(PrecisionType p) {
        FloatAP result(p);
        result._significand=shift_left(BigUnsigned(1u),p.bits()-1u);
        result._exponent = std::numeric_limits<Exponent>::min();
        return result;
    }
    static FloatAP max(PrecisionType p) {
        FloatAP result(p);
        const BigUnsigned power=shift_left(BigUnsigned(1u),p.bits());
        const BigUnsigned one(1u);
        result._significand=BigUnsigned::sub_shifted(
            power.limbs(),0u,one.limbs(),0u);
        result._exponent = std::numeric_limits<Exponent>::max();
        return result;
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
    friend bool eq(FloatAP const& a, FloatAP const& b) { return a == b; }
    friend bool lt(FloatAP const& a, FloatAP const& b) { return a < b; }
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
    friend PositiveUpperBound<FloatAP> mag(FloatAP const& x);
    friend PositiveLowerBound<FloatAP> mig(FloatAP const& x);
    friend FloatAP abs(RoundingModeType, FloatAP const& x);
    friend FloatAP mag(RoundingModeType mode, FloatAP const& x);
    friend FloatAP& iadd(RoundingModeType mode, FloatAP& a, FloatAP const& b);
    friend FloatAP& isub(RoundingModeType mode, FloatAP& a, FloatAP const& b);
    friend FloatAP& imul(RoundingModeType mode, FloatAP& a, FloatAP const& b);
    friend FloatAP& idiv(RoundingModeType mode, FloatAP& a, FloatAP const& b);
    friend FloatAP add(RoundingModeType mode, FloatAP const& a, FloatDP const& b);
    friend FloatAP add(RoundingModeType mode, FloatDP const& a, FloatAP const& b);
    friend FloatAP sub(RoundingModeType mode, FloatAP const& a, FloatDP const& b);
    friend FloatAP sub(RoundingModeType mode, FloatDP const& a, FloatAP const& b);
    friend FloatAP mul(RoundingModeType mode, FloatAP const& a, FloatDP const& b);
    friend FloatAP mul(RoundingModeType mode, FloatDP const& a, FloatAP const& b);
    friend FloatAP div(RoundingModeType mode, FloatAP const& a, FloatDP const& b);
    friend FloatAP div(RoundingModeType mode, FloatDP const& a, FloatAP const& b);
    friend FloatAP add(RoundingModeType mode, FloatAP const& a, ExactDouble const& b);
    friend FloatAP add(RoundingModeType mode, ExactDouble const& a, FloatAP const& b);
    friend FloatAP sub(RoundingModeType mode, FloatAP const& a, ExactDouble const& b);
    friend FloatAP sub(RoundingModeType mode, ExactDouble const& a, FloatAP const& b);
    friend FloatAP mul(RoundingModeType mode, FloatAP const& a, ExactDouble const& b);
    friend FloatAP mul(RoundingModeType mode, ExactDouble const& a, FloatAP const& b);
    friend FloatAP div(RoundingModeType mode, FloatAP const& a, ExactDouble const& b);
    friend FloatAP div(RoundingModeType mode, ExactDouble const& a, FloatAP const& b);
    friend FloatAP next(RoundUpward, FloatAP const& x);
    friend FloatAP next(RoundDownward, FloatAP const& x);
    friend FloatAP floor(FloatAP const& x);
    friend FloatAP ceil(FloatAP const& x);
    friend FloatAP round(FloatAP const& x);
    friend FloatAP operator*(TwoExp const& t, FloatAP const& x);
    friend FloatAP operator*(FloatAP const& x, TwoExp const& t);
    friend FloatAP operator/(FloatAP const& x, TwoExp const& t);
    friend FloatAP mul(FloatAP const& x, TwoExp const& t);
    friend FloatAP div(FloatAP const& x, TwoExp const& t);
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
    friend Bounds<FloatAP> fma(FloatAP const& a, FloatAP const& b, FloatAP const& c);
    friend Bounds<FloatAP> sqr(FloatAP const& x);
    friend Bounds<FloatAP> rec(FloatAP const& x);
    friend Bounds<FloatAP> pow(FloatAP const& x, Nat n);
    friend Bounds<FloatAP> pow(FloatAP const& x, Int n);
    // Definitions of the outward-rounded overloads are provided by
    // Bounds<FloatAP> (float_bounds.hpp), not duplicated in FloatAP.
    friend Bounds<FloatAP> sqrt(FloatAP const& x);
    friend Bounds<FloatAP> exp(FloatAP const& x);
    friend Bounds<FloatAP> log(FloatAP const& x);
    friend Bounds<FloatAP> sin(FloatAP const& x);
    friend Bounds<FloatAP> cos(FloatAP const& x);
    friend Bounds<FloatAP> tan(FloatAP const& x);
    friend Bounds<FloatAP> asin(FloatAP const& x);
    friend Bounds<FloatAP> acos(FloatAP const& x);
    friend Bounds<FloatAP> atan(FloatAP const& x);
    // Current-rounding overloads, using native FloatAP.
    friend FloatAP sqr(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP rec(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP add(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2);
    friend FloatAP sub(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2);
    friend FloatAP mul(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2);
    friend FloatAP div(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2);
    friend FloatAP fma(CurrentRoundingMode, FloatAP const& x1, FloatAP const& x2, FloatAP const& x3);
    friend FloatAP pow(CurrentRoundingMode, FloatAP const& x, Nat m);
    friend FloatAP pow(CurrentRoundingMode, FloatAP const& x, Int n);
    friend FloatAP sqrt(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP exp(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP log(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP sin(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP cos(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP tan(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP tanh(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP asin(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP acos(CurrentRoundingMode, FloatAP const& x);
    friend FloatAP atan(CurrentRoundingMode, FloatAP const& x);
    static FloatAP pi(CurrentRoundingMode, ArbitraryPrecision pr);
    friend FloatAP add(CurrentRoundingMode, FloatAP const& x1, FloatDP const& x2);
    friend FloatAP sub(CurrentRoundingMode, FloatAP const& x1, FloatDP const& x2);
    friend FloatAP mul(CurrentRoundingMode, FloatAP const& x1, FloatDP const& x2);
    friend FloatAP div(CurrentRoundingMode, FloatAP const& x1, FloatDP const& x2);
    friend FloatAP add(CurrentRoundingMode, FloatDP const& x1, FloatAP const& x2);
    friend FloatAP sub(CurrentRoundingMode, FloatDP const& x1, FloatAP const& x2);
    friend FloatAP mul(CurrentRoundingMode, FloatDP const& x1, FloatAP const& x2);
    friend FloatAP div(CurrentRoundingMode, FloatDP const& x1, FloatAP const& x2);
    friend FloatAP add(CurrentRoundingMode, FloatAP const& x1, ExactDouble x2);
    friend FloatAP sub(CurrentRoundingMode, FloatAP const& x1, ExactDouble x2);
    friend FloatAP mul(CurrentRoundingMode, FloatAP const& x1, ExactDouble x2);
    friend FloatAP div(CurrentRoundingMode, FloatAP const& x1, ExactDouble x2);
    friend FloatAP add(CurrentRoundingMode, ExactDouble x1, FloatAP const& x2);
    friend FloatAP sub(CurrentRoundingMode, ExactDouble x1, FloatAP const& x2);
    friend FloatAP mul(CurrentRoundingMode, ExactDouble x1, FloatAP const& x2);
    friend FloatAP div(CurrentRoundingMode, ExactDouble x1, FloatAP const& x2);
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
    friend FloatAP max(RoundingModeType mode, FloatAP const& a, FloatAP const& b) {
        const PrecisionType p(a._precision < b._precision ? a._precision : b._precision);
        if (a.is_nan() && b.is_nan()) return nan(p);
        if (a.is_nan()) return FloatAP(b,mode,p);
        if (b.is_nan()) return FloatAP(a,mode,p);
        if (a.is_zero() && b.is_zero()) return zero(static_cast<unsigned long>(p.bits()), a._negative && b._negative);
        return FloatAP(a < b ? b : a,mode,p);
    }
    friend FloatAP min(RoundingModeType mode, FloatAP const& a, FloatAP const& b) {
        const PrecisionType p(a._precision < b._precision ? a._precision : b._precision);
        if (a.is_nan() && b.is_nan()) return nan(p);
        if (a.is_nan()) return FloatAP(b,mode,p);
        if (b.is_nan()) return FloatAP(a,mode,p);
        if (a.is_zero() && b.is_zero()) return zero(static_cast<unsigned long>(p.bits()), a._negative || b._negative);
        return FloatAP(a < b ? a : b,mode,p);
    }
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
    friend FloatAP tanh(RoundingModeType mode, FloatAP const& x);
    friend FloatAP atan(RoundingModeType mode, FloatAP const& x);
    friend FloatAP asin(RoundingModeType mode, FloatAP const& x);
    friend FloatAP acos(RoundingModeType mode, FloatAP const& x);
    static FloatAP pi(RoundingModeType mode, PrecisionType precision);

    static FloatAP negate(FloatAP x);
    static FloatAP add_sub(FloatAP const& a, FloatAP const& b,
    unsigned long precision, Rounding mode, bool subtract);
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

// Shared certified short kernel; the binary64 adapter preserves the
// same enclosure/rounding contract without constructing Rational inputs.
namespace Detail {
    bool try_floatdp_transcendental(FloatTranscendental function,
        FloatAP::Rounding mode,double input,double& output);
}

template<class R, class A> R integer_cast(A const& value);
template<> inline Int integer_cast(const FloatAP& x) { return static_cast<Int>(x.get_d()); }
template<> inline Nat integer_cast(const FloatAP& x) { return static_cast<Nat>(x.get_d()); }

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

