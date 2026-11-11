/***************************************************************************
 *            numeric/floatmp.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @file numeric/floatmp.hpp
 *  @brief Native multiple-precision floating-point type and operations.
 */

#ifndef ARIADNE_FLOATMP_HPP
#define ARIADNE_FLOATMP_HPP

#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "numeric/big_unsigned.hpp"
#include "numeric/bits.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/number.hpp"
#include "numeric/rounding.hpp"
#include "numeric/sign.hpp"

namespace Ariadne {

class Integer;
class Rational;
class Decimal;
class Dyadic;
class ExactDouble;
class TwoExp;
namespace Detail {
enum class FloatTranscendental {
    sine,
    cosine,
    tangent,
    arctangent,
    arcsine,
    arccosine,
    pi,
    half_pi,
    negative_half_pi,
    exponential,
    logarithm,
    hyperbolic_tangent,
    logistic_sigmoid,
    softplus
};
struct FloatMPShortTranscendentals;
} // namespace Detail
struct RoundDownward;
struct RoundUpward;
struct RoundToNearest;
struct RoundTowardZero;
struct RoundApproximately;
struct CurrentRoundingMode;

// A size_t bit count is not necessarily representable by the signed
// exponent. Keep the narrowing check testable without allocating an
// infeasibly large GMP integer.
/** @brief Compute the multiprecision exponent from the significant bit count. */
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
/** @brief Runtime-selected binary precision for multiprecision arithmetic. */
class MultiplePrecision {
    std::size_t _bits;

  public:
    /** @brief Construct a \c MultiplePrecision from count. */
    explicit MultiplePrecision(Bits count) : _bits(static_cast<unsigned long>(count)) {
    }
    /** @brief Construct a \c MultiplePrecision from count. */
    explicit MultiplePrecision(std::size_t count) : _bits(count) {
    }
    /** @brief Return the number or representation of bits. */
    std::size_t bits() const noexcept {
        return _bits;
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, MultiplePrecision p) {
        return os << "mp(" << p.bits() << ")";
    }
    /** @brief Write the representation of the object. */
    friend OutputStream& repr(OutputStream& os, MultiplePrecision p) {
        return os << "MultiplePrecision(" << p.bits() << ")";
    }
    /** @brief Return the minimum of the operands. */
    friend MultiplePrecision min(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits <= b._bits ? a : b;
    }
    /** @brief Return the maximum of the operands. */
    friend MultiplePrecision max(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits >= b._bits ? a : b;
    }
    /** @brief Test equality of the operands. */
    friend bool operator==(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits == b._bits;
    }
    /** @brief Test whether the first operand is less than the second. */
    friend bool operator<(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits < b._bits;
    }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend bool operator<=(MultiplePrecision a, MultiplePrecision b) noexcept {
        return a._bits <= b._bits;
    }
    /** @brief Test whether the first operand is greater than the second. */
    friend bool operator>(MultiplePrecision a, MultiplePrecision b) noexcept {
        return b < a;
    }
    /** @brief Test whether the first operand is greater than or equal to the second. */
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
/** @brief Floating-point number parameterized by its precision representation. */
template<>
class Float<MultiplePrecision> {
  public:
    // A value type rather than an enum: Ariadne's rounding tags convert
    // directly to it, without making the GMP kernel depend on independent exact arithmetic.
    /** @brief Floating-point number parameterized by its precision representation. */
    class Rounding {
      public:
        // Without a fixed underlying type, the value range is exactly 0..3.
        // Every defined value therefore names one of these four modes.
        enum Value { ToNearest,
                     Downward,
                     Upward,
                     TowardZero };

      private:
        Value _value;

      public:
        /** @brief Construct the selected rounding policy. */
        constexpr Rounding(Value value) noexcept : _value(value) {}
        /** @brief Construct a \c Rounding with the specified characteristics. */
        Rounding(const RoundToNearest&) noexcept;
        /** @brief Construct a \c Rounding with the specified characteristics. */
        Rounding(const RoundDownward&) noexcept;
        /** @brief Construct a \c Rounding with the specified characteristics. */
        Rounding(const RoundUpward&) noexcept;
        /** @brief Construct a \c Rounding with the specified characteristics. */
        Rounding(const RoundTowardZero&) noexcept;
        /** @brief Construct a \c Rounding with the specified characteristics. */
        Rounding(const RoundApproximately&) noexcept;
        /** @brief Construct a \c Rounding with the specified characteristics. */
        Rounding(const CurrentRoundingMode&) noexcept;
        /** @brief Construct a \c Rounding with the specified characteristics. */
        Rounding(const Ariadne::Rounding&) noexcept;
        /** @brief Return the stored value. */
        constexpr Value value() const noexcept { return _value; }
        /** @brief Test equality of the operands. */
        friend constexpr bool operator==(const Rounding&, const Rounding&) = default;
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
        Finite,
        Infinity,
        NaN
    };
    using Exponent = std::int64_t;

  private:
    friend struct Detail::FloatMPShortTranscendentals;
    // Exact nonnegative magnitude; sign and special-value kind are separate.
    /** @brief Construct the exact rational from native numerator and denominator values. */
    static FloatMP from_native_ratio(const Integer& numerator,
                                     const Integer& denominator,
                                     Exponent scale, unsigned long precision,
                                     Rounding mode, bool negative_zero = false,
                                     bool invert_sign = false);
    /** @brief Construct the exact rational from unsigned numerator and denominator values. */
    static FloatMP from_unsigned_ratio(const BigUnsigned& numerator,
                                       const BigUnsigned& denominator,
                                       Exponent scale, unsigned long precision,
                                       Rounding mode, bool negative,
                                       bool negative_zero = false);
    BigUnsigned _significand;
    unsigned long _precision;
    Exponent _exponent = 0;
    Kind _kind = Kind::Finite;
    bool _negative = false;
    /** @brief Return the d general. */
    double get_d_general() const;

    /** @brief Add the operands with overflow checking. */
    static Exponent checked_add(Exponent a, Exponent b) {
        constexpr auto max = std::numeric_limits<Exponent>::max();
        constexpr auto min = std::numeric_limits<Exponent>::min();
        if ((b > 0 && a > max - b) || (b < 0 && a < min - b)) {
            throw std::overflow_error("FloatMP exponent overflow");
        }
        return a + b;
    }
    /** @brief Subtract the operands with overflow checking. */
    static Exponent checked_sub(Exponent a, Exponent b) {
        if (b == std::numeric_limits<Exponent>::min()) {
            return checked_add(checked_add(a, std::numeric_limits<Exponent>::max()), 1);
        }
        return checked_add(a, -b);
    }
    /** @brief Return the binary shift count. */
    static std::size_t shift_amount(Exponent n) {
        // Every caller supplies a nonnegative shift: a guarded exponent,
        // a checked distance between ordered exponents, precision padding,
        // or the square-root shift (at most the validated precision).
        // Only targets narrower than Exponent require a range check.
        if constexpr (std::numeric_limits<std::size_t>::digits <
                      std::numeric_limits<Exponent>::digits) {
            if (static_cast<std::uintmax_t>(n) >
                std::numeric_limits<std::size_t>::max()) {
                throw std::overflow_error("FloatMP shift out of range");
            }
        }
        return static_cast<std::size_t>(n);
    }
    /** @brief Return the number of significant binary digits. */
    static Exponent bit_length(const BigUnsigned& x) {
        return floatmp_exponent_from_bit_length(x.bit_length());
    }
    // Normalize and round an already-owned nonzero positive significand.
    /** @brief Round the binary-scaled integer to the target precision. */
    static void round_scaled_integer(FloatMP& result, Exponent scale, Rounding mode);
    /** @brief Round the short quotient using the remainder. */
    static void round_short(FloatMP& result, BigUnsigned::LimbView magnitude,
                            Exponent scale, Rounding mode);
    /** @brief Set the zero. */
    void set_zero(bool negative) {
        _kind = Kind::Finite;
        _negative = negative;
        _exponent = 0;
        _significand = BigUnsigned{};
    }

  public:
    /** @brief Construct a \c Float from precision. */
    explicit Float(unsigned long precision) : _precision(precision) {
        if (precision < 2 || precision > static_cast<std::uint64_t>(std::numeric_limits<Exponent>::max())) {
            throw std::invalid_argument("FloatMP precision must be >=2 and fit in exponent range");
        }
    }
    /** @brief Construct a \c Float from precision. */
    explicit Float(PrecisionType precision)
        : Float(static_cast<unsigned long>(precision.bits())) {
    }
    /** @brief Construct a \c Float from n, precision. */
    Float(std::int64_t n, PrecisionType precision);
    /** @brief Construct a \c Float from n, precision. */
    template<BuiltinUnsignedIntegral M>
    Float(M n, PrecisionType precision)
        : Float(Integer(n), ROUND_TO_NEAREST, precision) {}
    /** @brief Construct a \c Float from n, precision. */
    Float(const ExactDouble& n, PrecisionType precision);
    /** @brief Construct a \c Float from n, precision. */
    Float(const Dyadic& n, PrecisionType precision);
    /** @brief Construct a \c Float from value, precision. */
    explicit Float(const std::string& value, PrecisionType precision);
    /** @brief Construct a \c Float from n, precision. */
    Float(const TwoExp& n, PrecisionType precision);
    /** @brief Construct a \c Float from n, precision. */
    Float(const FloatDP& n, PrecisionType precision);
    /** @brief Construct a \c Float from n, precision. */
    Float(const FloatMP& n, PrecisionType precision);
    /** @brief Construct a \c Float from n, mode, precision. */
    Float(const FloatMP& n, RoundingModeType mode, PrecisionType precision);
    /** @brief Construct a \c Float from n, mode, precision. */
    Float(const FloatDP& n, RoundingModeType mode, PrecisionType precision);
    /** @brief Construct a \c Float from n, mode, precision. */
    Float(const ExactDouble& n, RoundingModeType mode, PrecisionType precision);
    /** @brief Construct a \c Float from n, mode, precision. */
    Float(const Dyadic& n, RoundingModeType mode, PrecisionType precision);
    /** @brief Construct a \c Float from n, mode, precision. */
    Float(const Decimal& n, RoundingModeType mode, PrecisionType precision);
    /** @brief Construct a \c Float from n, precision. */
    Float(const Integer& n, PrecisionType precision)
        : Float(n, ROUND_TO_NEAREST, precision) {}

    /** @brief Construct a \c Float from n, mode, precision. */
    Float(const Integer& n, RoundingModeType mode, PrecisionType precision);
    /** @brief Construct a \c Float from n, mode, precision. */
    template<BuiltinIntegral N>
    Float(N n, RoundingModeType mode, PrecisionType precision)
        : Float(Integer(n), mode, precision) {}
    /** @brief Construct a \c Float from n, mode, precision. */
    Float(const Rational& n, RoundingModeType mode, PrecisionType precision);
    //! Round an exact signed native numerator/denominator once.
    //! No canonical Rational intermediate (which would incur GCD).
    //! A zero denominator throws std::domain_error.
    static FloatMP from_ratio(const Integer& numerator,
                              const Integer& denominator,
                              RoundingModeType mode, PrecisionType precision);
    /** @brief Construct a \c Float from n, precision. */
    explicit Float(const Rational& n, PrecisionType precision)
        : Float(n, ROUND_TO_NEAREST, precision) {
    }
    /** @brief Convert to rational. */
    explicit operator Rational() const;
    /** @brief Convert to dyadic. */
    explicit operator Dyadic() const;
    /** @brief Convert to a built-in double. */
    double get_d() const {
        if constexpr (sizeof(double) == sizeof(std::uint64_t) &&
                      std::numeric_limits<double>::is_iec559 &&
                      std::numeric_limits<double>::digits == 53) {
            // A normal binary64 with <=53 significant bits is exact in
            // every rounding mode. Keep this frequent conversion inline.
            if (_kind == Kind::Finite && _precision <= 53u) {
                if (_significand.is_zero())
                    return _negative ? -0.0 : 0.0;
                if (_exponent >= -1022 - static_cast<Exponent>(_precision - 1u) &&
                    _exponent <= 1023 - static_cast<Exponent>(_precision - 1u)) {
                    const auto highest = _exponent + static_cast<Exponent>(_precision - 1u);
                    const std::uint64_t fraction = (_significand.limb(0u) << (53u - _precision)) &
                                                   ((std::uint64_t{1} << 52u) - 1u);
                    return std::bit_cast<double>((static_cast<std::uint64_t>(_negative) << 63u) |
                                                 (static_cast<std::uint64_t>(highest + 1023) << 52u) | fraction);
                }
            }
        }
        return get_d_general();
    }
    /** @brief Return the literal representation. */
    std::string literal() const;
    /** @brief Return the literal representation. */
    std::string literal(RoundingModeType mode) const;
    static Nat output_places;
    /** @brief Set the output places. */
    static void set_output_places(Nat places);
    /** @brief Print the value. */
    friend std::string print(const FloatMP& x, DecimalPlaces places, RoundingModeType mode);
    /** @brief Print the value. */
    friend std::string print(const FloatMP& x, DecimalPrecision figures, RoundingModeType mode);
    /** @brief Write the value to the output stream. */
    friend OutputStream& write(OutputStream& out, const FloatMP& x, DecimalPlaces places, RoundingModeType mode);
    /** @brief Write the value to the output stream. */
    friend OutputStream& write(OutputStream& out, const FloatMP& x, DecimalPrecision figures, RoundingModeType mode);
    /** @brief Write the representation of the object. */
    friend OutputStream& repr(OutputStream& out, const FloatMP& x, RoundingModeType mode);
    /** @brief Write the representation of the object. */
    friend OutputStream& repr(OutputStream& out, const FloatMP& x);
    /** @brief Destroy the object and release its owned resources. */
    ~Float() = default;
    /** @brief Copy-construct the object from the supplied instance. */
    Float(const FloatMP& x)
        : _significand(x._significand), _precision(x._precision),
          _exponent(x._exponent), _kind(x._kind), _negative(x._negative) {}
    /** @brief Move-construct the object from the supplied instance. */
    Float(FloatMP&& x) noexcept
        : _significand(std::move(x._significand)), _precision(x._precision),
          _exponent(x._exponent), _kind(x._kind), _negative(x._negative) {}
    /** @brief Assign the supplied value and return this object. */
    FloatMP& operator=(const FloatMP& x) {
        if (this != &x) {
            FloatMP copy(x);
            swap(copy);
        }
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    FloatMP& operator=(FloatMP&& x) noexcept {
        if (this != &x) {
            swap(x);
        }
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    FloatMP& operator=(const Integer& n);
    /** @brief Assign the supplied value and return this object. */
    FloatMP& operator=(const Dyadic& n);
    /** @brief Assign the supplied value and return this object. */
    FloatMP& operator=(const TwoExp& n);
    /** @brief Assign the supplied value and return this object. */
    FloatMP& operator=(const ExactDouble& n);
    /** @brief Assign an unsigned integer without narrowing through a signed type. */
    FloatMP& assign_unsigned(std::uint64_t n);
    /** @brief Assign the supplied value and return this object. */
    template<BuiltinIntegral N>
    FloatMP& operator=(N n) {
        if constexpr (std::is_signed_v<N>) {
            return *this = FloatMP(static_cast<std::int64_t>(n), precision());
        } else {
            return assign_unsigned(static_cast<std::uint64_t>(n));
        }
    }
    /** @brief Exchange the contents of the operands. */
    void swap(FloatMP& x) noexcept {
        _significand.swap(x._significand);
        std::swap(_precision, x._precision);
        std::swap(_exponent, x._exponent);
        std::swap(_kind, x._kind);
        std::swap(_negative, x._negative);
    }

    /** @brief Return the rounding mode. */
    static RoundingModeType get_rounding_mode() noexcept;
    /** @brief Set the rounding mode. */
    static void set_rounding_mode(RoundingModeType mode);
    /** @brief Set the rounding to nearest. */
    static void set_rounding_to_nearest() { set_rounding_mode(ROUND_TO_NEAREST); }
    /** @brief Set the rounding downward. */
    static void set_rounding_downward() { set_rounding_mode(ROUND_DOWNWARD); }
    /** @brief Set the rounding upward. */
    static void set_rounding_upward() { set_rounding_mode(ROUND_UPWARD); }
    /** @brief Set the rounding toward zero. */
    static void set_rounding_toward_zero() { set_rounding_mode(ROUND_TOWARD_ZERO); }
    /** @brief Return the default precision. */
    static PrecisionType get_default_precision();
    /** @brief Set the default precision. */
    static void set_default_precision(PrecisionType precision);
    /** @brief Set the precision. */
    void set_precision(PrecisionType precision);
    /** @brief Return the underlying numeric representation. */
    const FloatMP& raw() const noexcept { return *this; }
    /** @brief Return the underlying numeric representation. */
    FloatMP& raw() noexcept { return *this; }
    /** @brief Convert to exact number. */
    operator ExactNumber() const;
    /** @brief Return an enclosure enlarged by the supplied symmetric error. */
    Ball<FloatMP, FloatMP> pm(const Error<FloatMP>& error) const;
    /** @brief Return an enclosure enlarged by the supplied symmetric error. */
    Ball<FloatMP, FloatDP> pm(const Error<FloatDP>& error) const;

    /** @brief Return the working precision. */
    PrecisionType precision() const noexcept {
        return PrecisionType(_precision);
    }
    /** @brief Return the characteristics required to construct compatible values. */
    PrecisionType characteristics() const noexcept {
        return precision();
    }
    /** @brief Return the exponent. */
    Exponent exponent() const noexcept {
        return _exponent;
    }
    /** @brief Return the significand. */
    const BigUnsigned& significand() const noexcept {
        return _significand;
    }
    /** @brief Return the category discriminator. */
    Kind kind() const noexcept {
        return _kind;
    }
    /** @brief Test whether the value is negative. */
    bool negative() const noexcept {
        return _negative;
    }
    /** @brief Test whether the value is zero. */
    bool is_zero() const noexcept {
        return _kind == Kind::Finite && _significand.is_zero();
    }
    /** @brief Test whether the value is finite. */
    bool is_finite() const noexcept {
        return _kind == Kind::Finite;
    }
    /** @brief Test whether the value is infinite. */
    bool is_inf() const noexcept {
        return _kind == Kind::Infinity;
    }
    /** @brief Test whether the value is not a number. */
    bool is_nan() const noexcept {
        return _kind == Kind::NaN;
    }

    /** @brief Return the additive zero. */
    static FloatMP zero(unsigned long p, bool negative = false) {
        FloatMP z(p);
        z.set_zero(negative);
        return z;
    }
    /** @brief Return an infinite value. */
    static FloatMP infinity(unsigned long p, bool negative = false) {
        FloatMP z(p);
        z._kind = Kind::Infinity;
        z._negative = negative;
        return z;
    }
    /** @brief Return a not-a-number value. */
    static FloatMP nan(unsigned long p) {
        FloatMP z(p);
        z._kind = Kind::NaN;
        return z;
    }
    /** @brief Return an infinite value. */
    static FloatMP inf(Sign sign, PrecisionType p) {
        return sign == Sign::ZERO ? nan(p) : infinity(static_cast<unsigned long>(p.bits()), sign == Sign::NEGATIVE);
    }
    /** @brief Return an infinite value. */
    static FloatMP inf(PrecisionType p) {
        return infinity(static_cast<unsigned long>(p.bits()));
    }
    // Unit in the last place at 1, with the selected significand width.
    /** @brief Return the precision-dependent unit roundoff or tolerance. */
    static FloatMP eps(PrecisionType p) {
        FloatMP value(std::int64_t{1}, p);
        value._exponent = checked_sub(value._exponent,
                                      static_cast<Exponent>(p.bits() - 1));
        return value;
    }
    /** @brief Return a not-a-number value. */
    static FloatMP nan(PrecisionType p) {
        return nan(static_cast<unsigned long>(p.bits()));
    }
    // Exponent range comes from the signed native exponent representation.
    // No subnormal significands are supported by this kernel.
    /** @brief Return the minimum of the operands. */
    static FloatMP min(PrecisionType p) {
        FloatMP result(p);
        result._significand = shift_left(BigUnsigned(1u), p.bits() - 1u);
        result._exponent = std::numeric_limits<Exponent>::min();
        return result;
    }
    /** @brief Return the maximum of the operands. */
    static FloatMP max(PrecisionType p) {
        FloatMP result(p);
        const BigUnsigned power = shift_left(BigUnsigned(1u), p.bits());
        const BigUnsigned one(1u);
        result._significand = BigUnsigned::sub_shifted(
            power.limbs(), 0u, one.limbs(), 0u);
        result._exponent = std::numeric_limits<Exponent>::max();
        return result;
    }

    /** @brief Test whether the value is not a number. */
    friend bool is_nan(const FloatMP& x) noexcept {
        return x.is_nan();
    }
    /** @brief Test whether the value is infinite. */
    friend bool is_inf(const FloatMP& x) noexcept {
        return x.is_inf();
    }
    /** @brief Test whether the value is finite. */
    friend bool is_finite(const FloatMP& x) noexcept {
        return x.is_finite();
    }
    /** @brief Test whether the value is zero. */
    friend bool is_zero(const FloatMP& x) noexcept {
        return x.is_zero();
    }
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const FloatMP& a, const FloatMP& b);
    /** @brief Test equality of the operands. */
    friend bool eq(const FloatMP& a, const FloatMP& b) { return a == b; }
    /** @brief Test whether the first operand is less than the second. */
    friend bool lt(const FloatMP& a, const FloatMP& b) { return a < b; }
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const FloatMP& a, const Rational& b);
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const Rational& a, const FloatMP& b);
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const FloatMP& a, const FloatDP& b);
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const FloatDP& a, const FloatMP& b);
    /** @brief Test equality of the operands. */
    friend bool operator==(const FloatMP& a, const FloatDP& b);
    /** @brief Test inequality of the operands. */
    friend bool operator!=(const FloatMP& a, const FloatDP& b);
    /** @brief Test whether the first operand is less than the second. */
    friend bool operator<(const FloatMP& a, const FloatDP& b);
    /** @brief Test whether the first operand is greater than the second. */
    friend bool operator>(const FloatMP& a, const FloatDP& b);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend bool operator<=(const FloatMP& a, const FloatDP& b);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend bool operator>=(const FloatMP& a, const FloatDP& b);
    /** @brief Test equality of the operands. */
    friend bool operator==(const FloatDP& a, const FloatMP& b);
    /** @brief Test inequality of the operands. */
    friend bool operator!=(const FloatDP& a, const FloatMP& b);
    /** @brief Test whether the first operand is less than the second. */
    friend bool operator<(const FloatDP& a, const FloatMP& b);
    /** @brief Test whether the first operand is greater than the second. */
    friend bool operator>(const FloatDP& a, const FloatMP& b);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend bool operator<=(const FloatDP& a, const FloatMP& b);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend bool operator>=(const FloatDP& a, const FloatMP& b);
    /** @brief Convert to an integer using the selected integer representation. */
    friend Integer cast_integer(const FloatMP& value);
    /** @brief Test equality of the operands. */
    friend bool operator==(const FloatMP& a, const Rational& b) {
        return cmp(a, b) == Comparison::EQUAL;
    }
    /** @brief Test inequality of the operands. */
    friend bool operator!=(const FloatMP& a, const Rational& b) {
        return cmp(a, b) != Comparison::EQUAL;
    }
    /** @brief Test whether the first operand is less than the second. */
    friend bool operator<(const FloatMP& a, const Rational& b) {
        return cmp(a, b) == Comparison::LESS;
    }
    /** @brief Test whether the first operand is greater than the second. */
    friend bool operator>(const FloatMP& a, const Rational& b) {
        return cmp(a, b) == Comparison::GREATER;
    }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend bool operator<=(const FloatMP& a, const Rational& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::LESS || c == Comparison::EQUAL;
    }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend bool operator>=(const FloatMP& a, const Rational& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::GREATER || c == Comparison::EQUAL;
    }
    /** @brief Test equality of the operands. */
    friend bool operator==(const Rational& a, const FloatMP& b) {
        return cmp(a, b) == Comparison::EQUAL;
    }
    /** @brief Test inequality of the operands. */
    friend bool operator!=(const Rational& a, const FloatMP& b) {
        return cmp(a, b) != Comparison::EQUAL;
    }
    /** @brief Test whether the first operand is less than the second. */
    friend bool operator<(const Rational& a, const FloatMP& b) {
        return cmp(a, b) == Comparison::LESS;
    }
    /** @brief Test whether the first operand is greater than the second. */
    friend bool operator>(const Rational& a, const FloatMP& b) {
        return cmp(a, b) == Comparison::GREATER;
    }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend bool operator<=(const Rational& a, const FloatMP& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::LESS || c == Comparison::EQUAL;
    }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend bool operator>=(const Rational& a, const FloatMP& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::GREATER || c == Comparison::EQUAL;
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FloatMP& x);
    /** @brief Read the value from the input stream. */
    friend InputStream& operator>>(InputStream& is, FloatMP& x);
    /** @brief Test whether the operands have the same representation. */
    friend bool same(const FloatMP& a, const FloatMP& b) { return a == b; }
    /** @brief Return the absolute value. */
    friend Positive<FloatMP> abs(const FloatMP& x);
    /** @brief Return an upper bound on the absolute value. */
    friend PositiveUpperBound<FloatMP> mag(const FloatMP& x);
    /** @brief Return a lower bound on the absolute value. */
    friend PositiveLowerBound<FloatMP> mig(const FloatMP& x);
    /** @brief Return the absolute value. */
    friend FloatMP abs(RoundingModeType, const FloatMP& x);
    /** @brief Return an upper bound on the absolute value. */
    friend FloatMP mag(RoundingModeType mode, const FloatMP& x);
    /** @brief Add the operand in place. */
    friend FloatMP& iadd(RoundingModeType mode, FloatMP& a, const FloatMP& b);
    /** @brief Subtract the operand in place. */
    friend FloatMP& isub(RoundingModeType mode, FloatMP& a, const FloatMP& b);
    /** @brief Multiply by the operand in place. */
    friend FloatMP& imul(RoundingModeType mode, FloatMP& a, const FloatMP& b);
    /** @brief Divide by the operand in place. */
    friend FloatMP& idiv(RoundingModeType mode, FloatMP& a, const FloatMP& b);
    /** @brief Return the sum of the operands. */
    friend FloatMP add(RoundingModeType mode, const FloatMP& a, const FloatDP& b);
    /** @brief Return the sum of the operands. */
    friend FloatMP add(RoundingModeType mode, const FloatDP& a, const FloatMP& b);
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(RoundingModeType mode, const FloatMP& a, const FloatDP& b);
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(RoundingModeType mode, const FloatDP& a, const FloatMP& b);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(RoundingModeType mode, const FloatMP& a, const FloatDP& b);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(RoundingModeType mode, const FloatDP& a, const FloatMP& b);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(RoundingModeType mode, const FloatMP& a, const FloatDP& b);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(RoundingModeType mode, const FloatDP& a, const FloatMP& b);
    /** @brief Return the sum of the operands. */
    friend FloatMP add(RoundingModeType mode, const FloatMP& a, const ExactDouble& b);
    /** @brief Return the sum of the operands. */
    friend FloatMP add(RoundingModeType mode, const ExactDouble& a, const FloatMP& b);
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(RoundingModeType mode, const FloatMP& a, const ExactDouble& b);
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(RoundingModeType mode, const ExactDouble& a, const FloatMP& b);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(RoundingModeType mode, const FloatMP& a, const ExactDouble& b);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(RoundingModeType mode, const ExactDouble& a, const FloatMP& b);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(RoundingModeType mode, const FloatMP& a, const ExactDouble& b);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(RoundingModeType mode, const ExactDouble& a, const FloatMP& b);
    /** @brief Return the next term or representable value. */
    friend FloatMP next(RoundUpward, const FloatMP& x);
    /** @brief Return the next term or representable value. */
    friend FloatMP next(RoundDownward, const FloatMP& x);
    /** @brief Return the greatest integer not exceeding the value. */
    friend FloatMP floor(const FloatMP& x);
    /** @brief Return the least integer not less than the value. */
    friend FloatMP ceil(const FloatMP& x);
    /** @brief Round the value to an integer. */
    friend FloatMP round(const FloatMP& x);
    /** @brief Return the product. */
    friend FloatMP operator*(const TwoExp& t, const FloatMP& x);
    /** @brief Return the product. */
    friend FloatMP operator*(const FloatMP& x, const TwoExp& t);
    /** @brief Return the quotient. */
    friend FloatMP operator/(const FloatMP& x, const TwoExp& t);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(const FloatMP& x, const TwoExp& t);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(const FloatMP& x, const TwoExp& t);
    /** @brief Multiply by the operand in place. */
    friend FloatMP& operator*=(FloatMP& x, const TwoExp& t);
    /** @brief Divide by the operand in place. */
    friend FloatMP& operator/=(FloatMP& x, const TwoExp& t);
    /** @brief Return the sum or unary positive value. */
    friend Bounds<FloatMP> operator+(const FloatMP& a, const FloatMP& b);
    /** @brief Return the difference or additive inverse. */
    friend Bounds<FloatMP> operator-(const FloatMP& a, const FloatMP& b);
    /** @brief Return the product. */
    friend Bounds<FloatMP> operator*(const FloatMP& a, const FloatMP& b);
    /** @brief Return the quotient. */
    friend Bounds<FloatMP> operator/(const FloatMP& a, const FloatMP& b);
    /** @brief Return the sum of the operands. */
    friend Bounds<FloatMP> add(const FloatMP& a, const FloatMP& b);
    /** @brief Return the difference of the operands. */
    friend Bounds<FloatMP> sub(const FloatMP& a, const FloatMP& b);
    /** @brief Return the product of the operands. */
    friend Bounds<FloatMP> mul(const FloatMP& a, const FloatMP& b);
    /** @brief Return the quotient of the operands. */
    friend Bounds<FloatMP> div(const FloatMP& a, const FloatMP& b);
    /** @brief Compute the product of the first two operands plus the third. */
    friend Bounds<FloatMP> fma(const FloatMP& a, const FloatMP& b, const FloatMP& c);
    /** @brief Return the square. */
    friend Bounds<FloatMP> sqr(const FloatMP& x);
    /** @brief Return the multiplicative inverse. */
    friend Bounds<FloatMP> rec(const FloatMP& x);
    /** @brief Raise the value to the specified power. */
    friend Bounds<FloatMP> pow(const FloatMP& x, Nat n);
    /** @brief Raise the value to the specified power. */
    friend Bounds<FloatMP> pow(const FloatMP& x, Int n);
    // Definitions of the outward-rounded overloads are provided by
    // Bounds<FloatMP> (float_bounds.hpp), not duplicated in FloatMP.
    /** @brief Return the square root. */
    friend Bounds<FloatMP> sqrt(const FloatMP& x);
    /** @brief Return the exponential. */
    friend Bounds<FloatMP> exp(const FloatMP& x);
    /** @brief Return the natural logarithm. */
    friend Bounds<FloatMP> log(const FloatMP& x);
    /** @brief Return the sine. */
    friend Bounds<FloatMP> sin(const FloatMP& x);
    /** @brief Return the cosine. */
    friend Bounds<FloatMP> cos(const FloatMP& x);
    /** @brief Return the tangent. */
    friend Bounds<FloatMP> tan(const FloatMP& x);
    /** @brief Return the inverse sine. */
    friend Bounds<FloatMP> asin(const FloatMP& x);
    /** @brief Return the inverse cosine. */
    friend Bounds<FloatMP> acos(const FloatMP& x);
    /** @brief Return the inverse tangent. */
    friend Bounds<FloatMP> atan(const FloatMP& x);
    // Current-rounding overloads, using native FloatMP.
    /** @brief Return the square. */
    friend FloatMP sqr(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the multiplicative inverse. */
    friend FloatMP rec(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the sum of the operands. */
    friend FloatMP add(CurrentRoundingMode, const FloatMP& x1, const FloatMP& x2);
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(CurrentRoundingMode, const FloatMP& x1, const FloatMP& x2);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(CurrentRoundingMode, const FloatMP& x1, const FloatMP& x2);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(CurrentRoundingMode, const FloatMP& x1, const FloatMP& x2);
    /** @brief Compute the product of the first two operands plus the third. */
    friend FloatMP fma(CurrentRoundingMode, const FloatMP& x1, const FloatMP& x2, const FloatMP& x3);
    /** @brief Raise the value to the specified power. */
    friend FloatMP pow(CurrentRoundingMode, const FloatMP& x, Nat m);
    /** @brief Raise the value to the specified power. */
    friend FloatMP pow(CurrentRoundingMode, const FloatMP& x, Int n);
    /** @brief Return the square root. */
    friend FloatMP sqrt(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the exponential. */
    friend FloatMP exp(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the natural logarithm. */
    friend FloatMP log(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the sine. */
    friend FloatMP sin(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the cosine. */
    friend FloatMP cos(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the tangent. */
    friend FloatMP tan(CurrentRoundingMode, const FloatMP& x);
    //! \brief Correctly rounded logistic sigmoid using the active rounding mode.
    friend FloatMP sigmoid(CurrentRoundingMode, const FloatMP& x);
    //! \brief Correctly rounded softplus using the active rounding mode.
    friend FloatMP softplus(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the hyperbolic tangent. */
    friend FloatMP tanh(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the inverse sine. */
    friend FloatMP asin(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the inverse cosine. */
    friend FloatMP acos(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return the inverse tangent. */
    friend FloatMP atan(CurrentRoundingMode, const FloatMP& x);
    /** @brief Return pi at the requested precision. */
    static FloatMP pi(CurrentRoundingMode, MultiplePrecision pr);
    /** @brief Return the sum of the operands. */
    friend FloatMP add(CurrentRoundingMode, const FloatMP& x1, const FloatDP& x2);
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(CurrentRoundingMode, const FloatMP& x1, const FloatDP& x2);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(CurrentRoundingMode, const FloatMP& x1, const FloatDP& x2);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(CurrentRoundingMode, const FloatMP& x1, const FloatDP& x2);
    /** @brief Return the sum of the operands. */
    friend FloatMP add(CurrentRoundingMode, const FloatDP& x1, const FloatMP& x2);
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(CurrentRoundingMode, const FloatDP& x1, const FloatMP& x2);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(CurrentRoundingMode, const FloatDP& x1, const FloatMP& x2);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(CurrentRoundingMode, const FloatDP& x1, const FloatMP& x2);
    /** @brief Return the sum of the operands. */
    friend FloatMP add(CurrentRoundingMode, const FloatMP& x1, ExactDouble x2);
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(CurrentRoundingMode, const FloatMP& x1, ExactDouble x2);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(CurrentRoundingMode, const FloatMP& x1, ExactDouble x2);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(CurrentRoundingMode, const FloatMP& x1, ExactDouble x2);
    /** @brief Return the sum of the operands. */
    friend FloatMP add(CurrentRoundingMode, ExactDouble x1, const FloatMP& x2);
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(CurrentRoundingMode, ExactDouble x1, const FloatMP& x2);
    /** @brief Return the product of the operands. */
    friend FloatMP mul(CurrentRoundingMode, ExactDouble x1, const FloatMP& x2);
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(CurrentRoundingMode, ExactDouble x1, const FloatMP& x2);
    /** @brief Raise the value to the specified power. */
    friend FloatMP pow(RoundingModeType mode, const FloatMP& x, Int exponent);
    /** @brief Return the sum or unary positive value. */
    friend FloatMP operator+(const FloatMP& x) { return x; }
    /** @brief Return the difference or additive inverse. */
    friend FloatMP operator-(const FloatMP& x) { return negate(x); }
    /** @brief Return the additive zero with the operand characteristics. */
    friend FloatMP nul(const FloatMP& x) { return zero(x._precision); }
    /** @brief Return the unary positive value. */
    friend FloatMP pos(const FloatMP& x) { return x; }
    /** @brief Return the additive inverse. */
    friend FloatMP neg(const FloatMP& x) { return negate(x); }
    /** @brief Shift the value by the supplied number of binary digits. */
    friend FloatMP shft(const FloatMP& x, Int amount);
    /** @brief Divide the value by two. */
    friend FloatMP hlf(const FloatMP& x) { return shft(x, Int(-1)); }
    /** @brief Return the maximum of the operands. */
    friend FloatMP max(const FloatMP& a, const FloatMP& b) { return a < b ? b : a; }
    /** @brief Return the minimum of the operands. */
    friend FloatMP min(const FloatMP& a, const FloatMP& b) { return a < b ? a : b; }
    /** @brief Return the maximum of the operands. */
    friend FloatMP max(RoundingModeType mode, const FloatMP& a, const FloatMP& b) {
        const PrecisionType p(a._precision < b._precision ? a._precision : b._precision);
        if (a.is_nan() && b.is_nan())
            return nan(p);
        if (a.is_nan())
            return FloatMP(b, mode, p);
        if (b.is_nan())
            return FloatMP(a, mode, p);
        if (a.is_zero() && b.is_zero())
            return zero(static_cast<unsigned long>(p.bits()), a._negative && b._negative);
        return FloatMP(a < b ? b : a, mode, p);
    }
    /** @brief Return the minimum of the operands. */
    friend FloatMP min(RoundingModeType mode, const FloatMP& a, const FloatMP& b) {
        const PrecisionType p(a._precision < b._precision ? a._precision : b._precision);
        if (a.is_nan() && b.is_nan())
            return nan(p);
        if (a.is_nan())
            return FloatMP(b, mode, p);
        if (b.is_nan())
            return FloatMP(a, mode, p);
        if (a.is_zero() && b.is_zero())
            return zero(static_cast<unsigned long>(p.bits()), a._negative || b._negative);
        return FloatMP(a < b ? a : b, mode, p);
    }
    /** @brief Return the square. */
    friend FloatMP sqr(RoundingModeType mode, const FloatMP& x) {
        return FloatMP::mul(x, x, x._precision, mode);
    }
    /** @brief Return the multiplicative inverse. */
    friend FloatMP rec(RoundingModeType mode, const FloatMP& x) {
        const FloatMP one(std::int64_t{1}, x.precision());
        return FloatMP::div(one, x, x._precision, mode);
    }
    /** @brief Return the midpoint value. */
    friend FloatMP med(RoundingModeType mode, const FloatMP& a, const FloatMP& b) {
        const unsigned long precision = a._precision < b._precision
                                          ? a._precision
                                          : b._precision;
        return hlf(FloatMP::add(a, b, precision, mode));
    }
    /** @brief Return the radius. */
    friend FloatMP rad(RoundingModeType mode, const FloatMP& a, const FloatMP& b) {
        const unsigned long precision = a._precision < b._precision
                                          ? a._precision
                                          : b._precision;
        return hlf(FloatMP::sub(b, a, precision, mode));
    }
    /** @brief Test equality of the operands. */
    friend bool operator==(const FloatMP& a, const FloatMP& b) {
        return cmp(a, b) == Comparison::EQUAL;
    }
    /** @brief Test inequality of the operands. */
    friend bool operator!=(const FloatMP& a, const FloatMP& b) {
        return !(a == b);
    }
    /** @brief Test whether the first operand is less than the second. */
    friend bool operator<(const FloatMP& a, const FloatMP& b) {
        return cmp(a, b) == Comparison::LESS;
    }
    /** @brief Test whether the first operand is greater than the second. */
    friend bool operator>(const FloatMP& a, const FloatMP& b) {
        return cmp(a, b) == Comparison::GREATER;
    }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend bool operator<=(const FloatMP& a, const FloatMP& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::LESS || c == Comparison::EQUAL;
    }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend bool operator>=(const FloatMP& a, const FloatMP& b) {
        const Comparison c = cmp(a, b);
        return c == Comparison::GREATER || c == Comparison::EQUAL;
    }

    /** @brief Return the sum of the operands. */
    friend FloatMP add(RoundingModeType mode, const FloatMP& a, const FloatMP& b) {
        return add(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    /** @brief Return the difference of the operands. */
    friend FloatMP sub(RoundingModeType mode, const FloatMP& a, const FloatMP& b) {
        return sub(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    /** @brief Return the product of the operands. */
    friend FloatMP mul(RoundingModeType mode, const FloatMP& a, const FloatMP& b) {
        return mul(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    /** @brief Return the quotient of the operands. */
    friend FloatMP div(RoundingModeType mode, const FloatMP& a, const FloatMP& b) {
        return div(a, b, a._precision < b._precision ? a._precision : b._precision, mode);
    }
    /** @brief Compute the product of the first two operands plus the third. */
    friend FloatMP fma(RoundingModeType mode, const FloatMP& a,
                       const FloatMP& b, const FloatMP& c) {
        const auto p = a._precision < b._precision ? a._precision : b._precision;
        return fma(a, b, c, p < c._precision ? p : c._precision, mode);
    }

    // Square root is rounded from an exact integer-square-root enclosure.
    // Integer powers are computed exactly, then rounded only once.
    /** @brief Return the square root. */
    friend FloatMP sqrt(RoundingModeType mode, const FloatMP& x);
    /** @brief Raise the value to the specified power. */
    friend FloatMP pow(RoundingModeType mode, const FloatMP& x, Nat exponent);
    /** @brief Return the exponential. */
    friend FloatMP exp(RoundingModeType mode, const FloatMP& x);
    /** @brief Return the natural logarithm. */
    friend FloatMP log(RoundingModeType mode, const FloatMP& x);
    /** @brief Return the sine. */
    friend FloatMP sin(RoundingModeType mode, const FloatMP& x);
    /** @brief Return the cosine. */
    friend FloatMP cos(RoundingModeType mode, const FloatMP& x);
    /** @brief Return the tangent. */
    friend FloatMP tan(RoundingModeType mode, const FloatMP& x);
    //! \brief Correctly rounded logistic sigmoid.
    friend FloatMP sigmoid(RoundingModeType mode, const FloatMP& x);
    //! \brief Correctly rounded softplus.
    friend FloatMP softplus(RoundingModeType mode, const FloatMP& x);
    /** @brief Return the hyperbolic tangent. */
    friend FloatMP tanh(RoundingModeType mode, const FloatMP& x);
    /** @brief Return the inverse tangent. */
    friend FloatMP atan(RoundingModeType mode, const FloatMP& x);
    /** @brief Return the inverse sine. */
    friend FloatMP asin(RoundingModeType mode, const FloatMP& x);
    /** @brief Return the inverse cosine. */
    friend FloatMP acos(RoundingModeType mode, const FloatMP& x);
    /** @brief Return pi at the requested precision. */
    static FloatMP pi(RoundingModeType mode, PrecisionType precision);

    /** @brief Negate the value. */
    static FloatMP negate(FloatMP x);
    /** @brief Add or subtract according to the selected sign. */
    static FloatMP add_sub(const FloatMP& a, const FloatMP& b,
                           unsigned long precision, Rounding mode, bool subtract);
    /** @brief Return the sum of the operands. */
    static FloatMP add(const FloatMP& a, const FloatMP& b,
                       unsigned long precision, Rounding mode);
    /** @brief Return the difference of the operands. */
    static FloatMP sub(const FloatMP& a, const FloatMP& b,
                       unsigned long precision, Rounding mode);
    /** @brief Return the product of the operands. */
    static FloatMP mul(const FloatMP& a, const FloatMP& b,
                       unsigned long precision, Rounding mode);
    /** @brief Return the quotient of the operands. */
    static FloatMP div(const FloatMP& a, const FloatMP& b,
                       unsigned long precision, Rounding mode);
    // Exact product followed by a single rounding after addition.
    /** @brief Compute the product of the first two operands plus the third. */
    static FloatMP fma(const FloatMP& a, const FloatMP& b,
                       const FloatMP& c, unsigned long precision, Rounding mode);
};

template<class R, class A>
R integer_cast(const A& value);
/** @brief Convert the floating-point value to the requested integer type. */
template<>
inline Int integer_cast(const FloatMP& x) { return static_cast<Int>(x.get_d()); }
/** @brief Convert the floating-point value to the requested integer type. */
template<>
inline Nat integer_cast(const FloatMP& x) { return static_cast<Nat>(x.get_d()); }

inline constexpr FloatMP::RoundingModeType FloatMP::ROUND_TO_NEAREST =
    FloatMP::Rounding::ToNearest;
inline constexpr FloatMP::RoundingModeType FloatMP::ROUND_DOWNWARD =
    FloatMP::Rounding::Downward;
inline constexpr FloatMP::RoundingModeType FloatMP::ROUND_UPWARD =
    FloatMP::Rounding::Upward;
inline constexpr FloatMP::RoundingModeType FloatMP::ROUND_TOWARD_ZERO =
    FloatMP::Rounding::TowardZero;

} // namespace Ariadne
#endif

