/***************************************************************************
 *            numeric/floatmp.hpp
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

/*! \file numeric/floatmp.hpp
 *  \brief Native multiple-precision floating-point type and operations.
 */

#ifndef ARIADNE_FLOATMP_HPP
#define ARIADNE_FLOATMP_HPP

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
    struct FloatMPShortTranscendentals;
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
inline std::int64_t floatmp_exponent_from_bit_length(std::size_t bit_count) {
    if constexpr (std::numeric_limits<std::size_t>::digits >
                  std::numeric_limits<std::int64_t>::digits) {
        if (bit_count > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
            throw std::overflow_error("FloatMP operand too large");
        }
    }
    return static_cast<std::int64_t>(bit_count);
}

// The precision tag is deliberately independent of independent exact arithmetic types.
class MultiplePrecision {
    std::size_t _bits;
  public:
    explicit MultiplePrecision(Bits count) : _bits(static_cast<unsigned long>(count)) {
    }
    explicit MultiplePrecision(std::size_t count) : _bits(count) {
    }
    std::size_t bits() const noexcept {
        return _bits;
    }
    friend OutputStream& operator<<(OutputStream& os, MultiplePrecision p) {
        return os << "mp(" << p.bits() << ")";
    }
    friend OutputStream& repr(OutputStream& os, MultiplePrecision p) {
        return os << "MultiplePrecision(" << p.bits() << ")";
    }
    friend MultiplePrecision min(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits <= b._bits ? a : b;
    }
    friend MultiplePrecision max(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits >= b._bits ? a : b;
    }
    friend bool operator==(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits == b._bits;
    }
    friend bool operator<(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits < b._bits;
    }
    friend bool operator<=(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits <= b._bits;
    }
    friend bool operator>(MultiplePrecision a, MultiplePrecision b) noexcept {
        return b < a;
    }
    friend bool operator>=(MultiplePrecision a, MultiplePrecision b) noexcept {
        return b <= a;
    }
};

// The public native multiple-precision type. The API is growing; elementary
// arithmetic is implemented, while transcendental functions are pending.
using FloatMP = Float<MultiplePrecision>;

// independent exact arithmetic-independent floating-point representation and arithmetic.
// A finite nonzero value is (-1)^negative * significand * 2^exponent,
// with precisely precision significant bits and a native BigUnsigned
// significand. GMP remains only at external legacy/interval boundaries.
template<> class Float<MultiplePrecision> {
  public:
    // A value type rather than an enum: Ariadne's rounding tags convert
    // directly to it, without making the GMP kernel depend on independent exact arithmetic.
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
    using NumericType = FloatMP;
    using PrecisionType = MultiplePrecision;
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
    friend struct Detail::FloatMPShortTranscendentals;
    // Exact nonnegative magnitude; sign and special-value kind are separate.
    static FloatMP from_native_ratio(Integer const& numerator,
                                     Integer const& denominator,
                                     Exponent scale,unsigned long precision,
                                     Rounding mode,bool negative_zero=false,
                                     bool invert_sign=false);
    static FloatMP from_unsigned_ratio(BigUnsigned const& numerator,
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
            throw std::overflow_error("FloatMP exponent overflow");
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
        // Every caller supplies a nonnegative shift: a guarded exponent,
        // a checked distance between ordered exponents, precision padding,
        // or the square-root shift (at most the validated precision).
        // Only targets narrower than Exponent require a range check.
        if constexpr (std::numeric_limits<std::size_t>::digits <
                      std::numeric_limits<Exponent>::digits) {
            if (static_cast<std::uintmax_t>(n)>
                std::numeric_limits<std::size_t>::max()) {
                throw std::overflow_error("FloatMP shift out of range");
            }
        }
        return static_cast<std::size_t>(n);
    }
    static Exponent bit_length(BigUnsigned const& x) {
        return floatmp_exponent_from_bit_length(x.bit_length());
    }
    // Normalize and round an already-owned nonzero positive significand.
    static void round_scaled_integer(FloatMP& result, Exponent scale, Rounding mode);
    static void round_short(FloatMP& result,BigUnsigned::LimbView magnitude,
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
            throw std::invalid_argument("FloatMP precision must be >=2 and fit in exponent range");
        }
    }
    explicit Float(PrecisionType precision)
    : Float(static_cast<unsigned long>(precision.bits())) {
    }
    Float(std::int64_t n, PrecisionType precision);
    template<BuiltinUnsignedIntegral M>
    Float(M n, PrecisionType precision)
        : Float(Integer(n), ROUND_TO_NEAREST, precision) { }
    Float(ExactDouble const& n, PrecisionType precision);
    Float(Dyadic const& n, PrecisionType precision);
    explicit Float(std::string const& value, PrecisionType precision);
    Float(TwoExp const& n, PrecisionType precision);
    Float(FloatDP const& n, PrecisionType precision);
    Float(FloatMP const& n, PrecisionType precision);
    Float(FloatMP const& n, RoundingModeType mode, PrecisionType precision);
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
    static FloatMP from_ratio(Integer const& numerator,
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
    friend std::string print(FloatMP const& x, DecimalPlaces places, RoundingModeType mode);
    friend std::string print(FloatMP const& x, DecimalPrecision figures, RoundingModeType mode);
    friend OutputStream& write(OutputStream& out, FloatMP const& x, DecimalPlaces places, RoundingModeType mode);
    friend OutputStream& write(OutputStream& out, FloatMP const& x, DecimalPrecision figures, RoundingModeType mode);
    friend OutputStream& repr(OutputStream& out, FloatMP const& x, RoundingModeType mode);
    friend OutputStream& repr(OutputStream& out, FloatMP const& x);
    ~Float() = default;
    Float(FloatMP const& x)
    : _significand(x._significand), _precision(x._precision),
      _exponent(x._exponent), _kind(x._kind), _negative(x._negative) { }
    Float(FloatMP&& x) noexcept
    : _significand(std::move(x._significand)), _precision(x._precision),
      _exponent(x._exponent), _kind(x._kind), _negative(x._negative) { }
    FloatMP& operator=(FloatMP const& x) {
        if (this != &x) {
            FloatMP copy(x);
            swap(copy);
        }
        return *this;
    }
    FloatMP& operator=(FloatMP&& x) noexcept {
        if (this != &x) {
            swap(x);
        }
        return *this;
    }
    FloatMP& operator=(Integer const& n);
    FloatMP& operator=(Dyadic const& n);
    FloatMP& operator=(TwoExp const& n);
    FloatMP& operator=(ExactDouble const& n);
    FloatMP& assign_unsigned(std::uint64_t n);
    template<BuiltinIntegral N> FloatMP& operator=(N n) {
        if constexpr (std::is_signed_v<N>) {
            return *this = FloatMP(static_cast<std::int64_t>(n), precision());
        } else {
            return assign_unsigned(static_cast<std::uint64_t>(n));
        }
    }
    void swap(FloatMP& x) noexcept {
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
    FloatMP const& raw() const noexcept { return *this; }
    FloatMP& raw() noexcept { return *this; }
    operator ExactNumber() const;
    Ball<FloatMP,FloatMP> pm(Error<FloatMP> const& error) const;
    Ball<FloatMP,FloatDP> pm(Error<FloatDP> const& error) const;

    PrecisionType precision() const noexcept {
        return PrecisionType(_precision);
    }
    PrecisionType characteristics() const noexcept {
        return precision();
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

    static FloatMP zero(unsigned long p, bool negative = false) {
        FloatMP z(p);
        z.set_zero(negative);
        return z;
    }
    static FloatMP infinity(unsigned long p, bool negative = false) {
        FloatMP z(p);
        z._kind = Kind::Infinity;
        z._negative = negative;
        return z;
    }
    static FloatMP nan(unsigned long p) {
        FloatMP z(p);
        z._kind = Kind::NaN;
        return z;
    }
    static FloatMP inf(Sign sign, PrecisionType p) {
        return sign == Sign::ZERO ? nan(p) : infinity(static_cast<unsigned long>(p.bits()), sign == Sign::NEGATIVE);
    }
    static FloatMP inf(PrecisionType p) {
        return infinity(static_cast<unsigned long>(p.bits()));
    }
    // Unit in the last place at 1, with the selected significand width.
    static FloatMP eps(PrecisionType p) {
        FloatMP value(std::int64_t{1}, p);
        value._exponent = checked_sub(value._exponent,
                                      static_cast<Exponent>(p.bits() - 1));
        return value;
    }
    static FloatMP nan(PrecisionType p) {
        return nan(static_cast<unsigned long>(p.bits()));
    }
    // Exponent range comes from the signed native exponent representation.
    // No subnormal significands are supported by this kernel.
    static FloatMP min(PrecisionType p) {
        FloatMP result(p);
        result._significand=shift_left(BigUnsigned(1u),p.bits()-1u);
        result._exponent = std::numeric_limits<Exponent>::min();
        return result;
    }
    static FloatMP max(PrecisionType p) {
        FloatMP result(p);
        const BigUnsigned power=shift_left(BigUnsigned(1u),p.bits());
        const BigUnsigned one(1u);
        result._significand=BigUnsigned::sub_shifted(
            power.limbs(),0u,one.limbs(),0u);
        result._exponent = std::numeric_limits<Exponent>::max();
        return result;
    }

    friend bool is_nan(FloatMP const& x) noexcept {
        return x.is_nan();
    }
    friend bool is_inf(FloatMP const& x) noexcept {
        return x.is_inf();
    }
    friend bool is_finite(FloatMP const& x) noexcept {
        return x.is_finite();
    }
    friend bool is_zero(FloatMP const& x) noexcept {
        return x.is_zero();
    }
    friend Comparison cmp(FloatMP const& a, FloatMP const& b);
    friend bool eq(FloatMP const& a, FloatMP const& b) { return a == b; }
    friend bool lt(FloatMP const& a, FloatMP const& b) { return a < b; }
    friend Comparison cmp(FloatMP const& a, Rational const& b);
    friend Comparison cmp(Rational const& a, FloatMP const& b);
    friend Comparison cmp(FloatMP const& a, FloatDP const& b);
    friend Comparison cmp(FloatDP const& a, FloatMP const& b);
    friend bool operator==(FloatMP const& a, FloatDP const& b);
    friend bool operator!=(FloatMP const& a, FloatDP const& b);
    friend bool operator<(FloatMP const& a, FloatDP const& b);
    friend bool operator>(FloatMP const& a, FloatDP const& b);
    friend bool operator<=(FloatMP const& a, FloatDP const& b);
    friend bool operator>=(FloatMP const& a, FloatDP const& b);
    friend bool operator==(FloatDP const& a, FloatMP const& b);
    friend bool operator!=(FloatDP const& a, FloatMP const& b);
    friend bool operator<(FloatDP const& a, FloatMP const& b);
    friend bool operator>(FloatDP const& a, FloatMP const& b);
    friend bool operator<=(FloatDP const& a, FloatMP const& b);
    friend bool operator>=(FloatDP const& a, FloatMP const& b);
    friend Integer cast_integer(FloatMP const& value);
    friend bool operator==(FloatMP const& a, Rational const& b) {
        return cmp(a, b) == Comparison::EQUAL;
    }
    friend bool operator!=(FloatMP const& a, Rational const& b) {
        return cmp(a, b) != Comparison::EQUAL;
    }
    friend bool operator<(FloatMP const& a, Rational const& b) {
        return cmp(a, b) == Comparison::LESS;
    }
    friend bool operator>(FloatMP const& a, Rational const& b) {
        return cmp(a, b) == Comparison::GREATER;
    }
    friend bool operator<=(FloatMP const& a, Rational const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::LESS || c == Comparison::EQUAL;
    }
    friend bool operator>=(FloatMP const& a, Rational const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::GREATER || c == Comparison::EQUAL;
    }
    friend bool operator==(Rational const& a, FloatMP const& b) {
        return cmp(a, b) == Comparison::EQUAL;
    }
    friend bool operator!=(Rational const& a, FloatMP const& b) {
        return cmp(a, b) != Comparison::EQUAL;
    }
    friend bool operator<(Rational const& a, FloatMP const& b) {
        return cmp(a, b) == Comparison::LESS;
    }
    friend bool operator>(Rational const& a, FloatMP const& b) {
        return cmp(a, b) == Comparison::GREATER;
    }
    friend bool operator<=(Rational const& a, FloatMP const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::LESS || c == Comparison::EQUAL;
    }
    friend bool operator>=(Rational const& a, FloatMP const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::GREATER || c == Comparison::EQUAL;
    }
    friend OutputStream& operator<<(OutputStream& os, FloatMP const& x);
    friend InputStream& operator>>(InputStream& is, FloatMP& x);
    friend bool same(FloatMP const& a, FloatMP const& b) { return a == b; }
    friend Positive<FloatMP> abs(FloatMP const& x);
    friend PositiveUpperBound<FloatMP> mag(FloatMP const& x);
    friend PositiveLowerBound<FloatMP> mig(FloatMP const& x);
    friend FloatMP abs(RoundingModeType, FloatMP const& x);
    friend FloatMP mag(RoundingModeType mode, FloatMP const& x);
    friend FloatMP& iadd(RoundingModeType mode, FloatMP& a, FloatMP const& b);
    friend FloatMP& isub(RoundingModeType mode, FloatMP& a, FloatMP const& b);
    friend FloatMP& imul(RoundingModeType mode, FloatMP& a, FloatMP const& b);
    friend FloatMP& idiv(RoundingModeType mode, FloatMP& a, FloatMP const& b);
    friend FloatMP add(RoundingModeType mode, FloatMP const& a, FloatDP const& b);
    friend FloatMP add(RoundingModeType mode, FloatDP const& a, FloatMP const& b);
    friend FloatMP sub(RoundingModeType mode, FloatMP const& a, FloatDP const& b);
    friend FloatMP sub(RoundingModeType mode, FloatDP const& a, FloatMP const& b);
    friend FloatMP mul(RoundingModeType mode, FloatMP const& a, FloatDP const& b);
    friend FloatMP mul(RoundingModeType mode, FloatDP const& a, FloatMP const& b);
    friend FloatMP div(RoundingModeType mode, FloatMP const& a, FloatDP const& b);
    friend FloatMP div(RoundingModeType mode, FloatDP const& a, FloatMP const& b);
    friend FloatMP add(RoundingModeType mode, FloatMP const& a, ExactDouble const& b);
    friend FloatMP add(RoundingModeType mode, ExactDouble const& a, FloatMP const& b);
    friend FloatMP sub(RoundingModeType mode, FloatMP const& a, ExactDouble const& b);
    friend FloatMP sub(RoundingModeType mode, ExactDouble const& a, FloatMP const& b);
    friend FloatMP mul(RoundingModeType mode, FloatMP const& a, ExactDouble const& b);
    friend FloatMP mul(RoundingModeType mode, ExactDouble const& a, FloatMP const& b);
    friend FloatMP div(RoundingModeType mode, FloatMP const& a, ExactDouble const& b);
    friend FloatMP div(RoundingModeType mode, ExactDouble const& a, FloatMP const& b);
    friend FloatMP next(RoundUpward, FloatMP const& x);
    friend FloatMP next(RoundDownward, FloatMP const& x);
    friend FloatMP floor(FloatMP const& x);
    friend FloatMP ceil(FloatMP const& x);
    friend FloatMP round(FloatMP const& x);
    friend FloatMP operator*(TwoExp const& t, FloatMP const& x);
    friend FloatMP operator*(FloatMP const& x, TwoExp const& t);
    friend FloatMP operator/(FloatMP const& x, TwoExp const& t);
    friend FloatMP mul(FloatMP const& x, TwoExp const& t);
    friend FloatMP div(FloatMP const& x, TwoExp const& t);
    friend FloatMP& operator*=(FloatMP& x, TwoExp const& t);
    friend FloatMP& operator/=(FloatMP& x, TwoExp const& t);
    friend Bounds<FloatMP> operator+(FloatMP const& a, FloatMP const& b);
    friend Bounds<FloatMP> operator-(FloatMP const& a, FloatMP const& b);
    friend Bounds<FloatMP> operator*(FloatMP const& a, FloatMP const& b);
    friend Bounds<FloatMP> operator/(FloatMP const& a, FloatMP const& b);
    friend Bounds<FloatMP> add(FloatMP const& a, FloatMP const& b);
    friend Bounds<FloatMP> sub(FloatMP const& a, FloatMP const& b);
    friend Bounds<FloatMP> mul(FloatMP const& a, FloatMP const& b);
    friend Bounds<FloatMP> div(FloatMP const& a, FloatMP const& b);
    friend Bounds<FloatMP> fma(FloatMP const& a, FloatMP const& b, FloatMP const& c);
    friend Bounds<FloatMP> sqr(FloatMP const& x);
    friend Bounds<FloatMP> rec(FloatMP const& x);
    friend Bounds<FloatMP> pow(FloatMP const& x, Nat n);
    friend Bounds<FloatMP> pow(FloatMP const& x, Int n);
    // Definitions of the outward-rounded overloads are provided by
    // Bounds<FloatMP> (float_bounds.hpp), not duplicated in FloatMP.
    friend Bounds<FloatMP> sqrt(FloatMP const& x);
    friend Bounds<FloatMP> exp(FloatMP const& x);
    friend Bounds<FloatMP> log(FloatMP const& x);
    friend Bounds<FloatMP> sin(FloatMP const& x);
    friend Bounds<FloatMP> cos(FloatMP const& x);
    friend Bounds<FloatMP> tan(FloatMP const& x);
    friend Bounds<FloatMP> asin(FloatMP const& x);
    friend Bounds<FloatMP> acos(FloatMP const& x);
    friend Bounds<FloatMP> atan(FloatMP const& x);
    // Current-rounding overloads, using native FloatMP.
    friend FloatMP sqr(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP rec(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP add(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP sub(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP mul(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP div(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2);
    friend FloatMP fma(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2, FloatMP const& x3);
    friend FloatMP pow(CurrentRoundingMode, FloatMP const& x, Nat m);
    friend FloatMP pow(CurrentRoundingMode, FloatMP const& x, Int n);
    friend FloatMP sqrt(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP exp(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP log(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP sin(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP cos(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP tan(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP tanh(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP asin(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP acos(CurrentRoundingMode, FloatMP const& x);
    friend FloatMP atan(CurrentRoundingMode, FloatMP const& x);
    static FloatMP pi(CurrentRoundingMode, MultiplePrecision pr);
    friend FloatMP add(CurrentRoundingMode, FloatMP const& x1, FloatDP const& x2);
    friend FloatMP sub(CurrentRoundingMode, FloatMP const& x1, FloatDP const& x2);
    friend FloatMP mul(CurrentRoundingMode, FloatMP const& x1, FloatDP const& x2);
    friend FloatMP div(CurrentRoundingMode, FloatMP const& x1, FloatDP const& x2);
    friend FloatMP add(CurrentRoundingMode, FloatDP const& x1, FloatMP const& x2);
    friend FloatMP sub(CurrentRoundingMode, FloatDP const& x1, FloatMP const& x2);
    friend FloatMP mul(CurrentRoundingMode, FloatDP const& x1, FloatMP const& x2);
    friend FloatMP div(CurrentRoundingMode, FloatDP const& x1, FloatMP const& x2);
    friend FloatMP add(CurrentRoundingMode, FloatMP const& x1, ExactDouble x2);
    friend FloatMP sub(CurrentRoundingMode, FloatMP const& x1, ExactDouble x2);
    friend FloatMP mul(CurrentRoundingMode, FloatMP const& x1, ExactDouble x2);
    friend FloatMP div(CurrentRoundingMode, FloatMP const& x1, ExactDouble x2);
    friend FloatMP add(CurrentRoundingMode, ExactDouble x1, FloatMP const& x2);
    friend FloatMP sub(CurrentRoundingMode, ExactDouble x1, FloatMP const& x2);
    friend FloatMP mul(CurrentRoundingMode, ExactDouble x1, FloatMP const& x2);
    friend FloatMP div(CurrentRoundingMode, ExactDouble x1, FloatMP const& x2);
    friend FloatMP pow(RoundingModeType mode, FloatMP const& x, Int exponent);
    friend FloatMP operator+(FloatMP const& x) { return x; }
    friend FloatMP operator-(FloatMP const& x) { return negate(x); }
    friend FloatMP nul(FloatMP const& x) { return zero(x._precision); }
    friend FloatMP pos(FloatMP const& x) { return x; }
    friend FloatMP neg(FloatMP const& x) { return negate(x); }
    friend FloatMP shft(FloatMP const& x, Int amount);
    friend FloatMP hlf(FloatMP const& x) { return shft(x, Int(-1)); }
    friend FloatMP max(FloatMP const& a, FloatMP const& b) { return a < b ? b : a; }
    friend FloatMP min(FloatMP const& a, FloatMP const& b) { return a < b ? a : b; }
    friend FloatMP max(RoundingModeType mode, FloatMP const& a, FloatMP const& b) {
        const PrecisionType p(a._precision < b._precision ? a._precision : b._precision);
        if (a.is_nan() && b.is_nan()) return nan(p);
        if (a.is_nan()) return FloatMP(b,mode,p);
        if (b.is_nan()) return FloatMP(a,mode,p);
        if (a.is_zero() && b.is_zero()) return zero(static_cast<unsigned long>(p.bits()), a._negative && b._negative);
        return FloatMP(a < b ? b : a,mode,p);
    }
    friend FloatMP min(RoundingModeType mode, FloatMP const& a, FloatMP const& b) {
        const PrecisionType p(a._precision < b._precision ? a._precision : b._precision);
        if (a.is_nan() && b.is_nan()) return nan(p);
        if (a.is_nan()) return FloatMP(b,mode,p);
        if (b.is_nan()) return FloatMP(a,mode,p);
        if (a.is_zero() && b.is_zero()) return zero(static_cast<unsigned long>(p.bits()), a._negative || b._negative);
        return FloatMP(a < b ? a : b,mode,p);
    }
    friend FloatMP sqr(RoundingModeType mode, FloatMP const& x) {
        return FloatMP::mul(x, x, x._precision, mode);
    }
    friend FloatMP rec(RoundingModeType mode, FloatMP const& x) {
        const FloatMP one(std::int64_t{1}, x.precision());
        return FloatMP::div(one, x, x._precision, mode);
    }
    friend FloatMP med(RoundingModeType mode, FloatMP const& a, FloatMP const& b) {
        const unsigned long precision = a._precision < b._precision
            ? a._precision : b._precision;
        return hlf(FloatMP::add(a, b, precision, mode));
    }
    friend FloatMP rad(RoundingModeType mode, FloatMP const& a, FloatMP const& b) {
        const unsigned long precision = a._precision < b._precision
            ? a._precision : b._precision;
        return hlf(FloatMP::sub(b, a, precision, mode));
    }
    friend bool operator==(FloatMP const& a, FloatMP const& b) {
        return cmp(a, b) == Comparison::EQUAL;
    }
    friend bool operator!=(FloatMP const& a, FloatMP const& b) {
        return !(a == b);
    }
    friend bool operator<(FloatMP const& a, FloatMP const& b) {
        return cmp(a, b) == Comparison::LESS;
    }
    friend bool operator>(FloatMP const& a, FloatMP const& b) {
        return cmp(a, b) == Comparison::GREATER;
    }
    friend bool operator<=(FloatMP const& a, FloatMP const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::LESS || c == Comparison::EQUAL;
    }
    friend bool operator>=(FloatMP const& a, FloatMP const& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::GREATER || c == Comparison::EQUAL;
    }

    friend FloatMP add(RoundingModeType mode, FloatMP const& a, FloatMP const& b) {
        return add(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    friend FloatMP sub(RoundingModeType mode, FloatMP const& a, FloatMP const& b) {
        return sub(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    friend FloatMP mul(RoundingModeType mode, FloatMP const& a, FloatMP const& b) {
        return mul(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    friend FloatMP div(RoundingModeType mode, FloatMP const& a, FloatMP const& b) {
        return div(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    friend FloatMP fma(RoundingModeType mode, FloatMP const& a,
    FloatMP const& b, FloatMP const& c) {
        const auto p = a._precision < b._precision ? a._precision : b._precision;
        return fma(a, b, c, p < c._precision ? p : c._precision, mode);
    }

    // Square root is rounded from an exact integer-square-root enclosure.
    // Integer powers are computed exactly, then rounded only once.
    friend FloatMP sqrt(RoundingModeType mode, FloatMP const& x);
    friend FloatMP pow(RoundingModeType mode, FloatMP const& x, Nat exponent);
    friend FloatMP exp(RoundingModeType mode, FloatMP const& x);
    friend FloatMP log(RoundingModeType mode, FloatMP const& x);
    friend FloatMP sin(RoundingModeType mode, FloatMP const& x);
    friend FloatMP cos(RoundingModeType mode, FloatMP const& x);
    friend FloatMP tan(RoundingModeType mode, FloatMP const& x);
    friend FloatMP tanh(RoundingModeType mode, FloatMP const& x);
    friend FloatMP atan(RoundingModeType mode, FloatMP const& x);
    friend FloatMP asin(RoundingModeType mode, FloatMP const& x);
    friend FloatMP acos(RoundingModeType mode, FloatMP const& x);
    static FloatMP pi(RoundingModeType mode, PrecisionType precision);

    static FloatMP negate(FloatMP x);
    static FloatMP add_sub(FloatMP const& a, FloatMP const& b,
    unsigned long precision, Rounding mode, bool subtract);
    static FloatMP add(FloatMP const& a, FloatMP const& b,
    unsigned long precision, Rounding mode);
    static FloatMP sub(FloatMP const& a, FloatMP const& b,
    unsigned long precision, Rounding mode);
    static FloatMP mul(FloatMP const& a, FloatMP const& b,
    unsigned long precision, Rounding mode);
    static FloatMP div(FloatMP const& a, FloatMP const& b,
    unsigned long precision, Rounding mode);
    // Exact product followed by a single rounding after addition.
    static FloatMP fma(FloatMP const& a, FloatMP const& b,
    FloatMP const& c, unsigned long precision, Rounding mode);
};


template<class R, class A> R integer_cast(A const& value);
template<> inline Int integer_cast(const FloatMP& x) { return static_cast<Int>(x.get_d()); }
template<> inline Nat integer_cast(const FloatMP& x) { return static_cast<Nat>(x.get_d()); }

inline constexpr FloatMP::RoundingModeType FloatMP::ROUND_TO_NEAREST =
    FloatMP::Rounding::ToNearest;
inline constexpr FloatMP::RoundingModeType FloatMP::ROUND_DOWNWARD =
    FloatMP::Rounding::Downward;
inline constexpr FloatMP::RoundingModeType FloatMP::ROUND_UPWARD =
    FloatMP::Rounding::Upward;
inline constexpr FloatMP::RoundingModeType FloatMP::ROUND_TOWARD_ZERO =
    FloatMP::Rounding::TowardZero;

}
#endif

