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

namespace Ariadne {

class Integer;
class Rational;

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
    friend bool operator==(ArbitraryPrecision a, ArbitraryPrecision b) noexcept {
        return a._bits == b._bits;
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
    enum class Rounding {
        nearest, downward, upward, toward_zero
    };
    using PrecisionType = ArbitraryPrecision;
    using RoundingModeType = Rounding;
    static constexpr RoundingModeType ROUND_TO_NEAREST = Rounding::nearest;
    static constexpr RoundingModeType ROUND_DOWNWARD = Rounding::downward;
    static constexpr RoundingModeType ROUND_UPWARD = Rounding::upward;
    static constexpr RoundingModeType ROUND_TOWARD_ZERO = Rounding::toward_zero;
    enum class Kind {
        finite, infinity, nan
    };
    using Exponent = std::int64_t;

  private:
    mpz_t _significand;
    unsigned long _precision;
    Exponent _exponent = 0;
    Kind _kind = Kind::finite;
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
        _kind = Kind::finite;
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
    Float(Integer const& n, RoundingModeType mode, PrecisionType precision);
    Float(Rational const& n, RoundingModeType mode, PrecisionType precision);
    explicit Float(Rational const& n, PrecisionType precision)
    : Float(n, ROUND_TO_NEAREST, precision) {
    }
    explicit operator Rational() const;
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
        return* this;
    }
    void swap(FloatAP& x) noexcept {
        mpz_swap(_significand, x._significand);
        std::swap(_precision, x._precision);
        std::swap(_exponent, x._exponent);
        std::swap(_kind, x._kind);
        std::swap(_negative, x._negative);
    }

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
        return _kind == Kind::finite && mpz_size(_significand) == 0;
    }
    bool is_finite() const noexcept {
        return _kind == Kind::finite;
    }
    bool is_inf() const noexcept {
        return _kind == Kind::infinity;
    }
    bool is_nan() const noexcept {
        return _kind == Kind::nan;
    }

    static FloatAP zero(unsigned long p, bool negative = false) {
        FloatAP z(p);
        z.set_zero(negative);
        return z;
    }
    static FloatAP infinity(unsigned long p, bool negative = false) {
        FloatAP z(p);
        z._kind = Kind::infinity;
        z._negative = negative;
        return z;
    }
    static FloatAP nan(unsigned long p) {
        FloatAP z(p);
        z._kind = Kind::nan;
        return z;
    }
    static FloatAP inf(Sign sign, PrecisionType p) {
        return infinity(static_cast<unsigned long>(p.bits()), sign == Sign::NEGATIVE);
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

}
#endif

