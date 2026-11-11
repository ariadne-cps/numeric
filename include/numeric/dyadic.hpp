/***************************************************************************
 *            numeric/dyadic.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file numeric/dyadic.hpp
 *  @brief Dyadic numbers.
 */

#ifndef ARIADNE_DYADIC_HPP
#define ARIADNE_DYADIC_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "foundation/logical.hpp"
#include "numeric/builtin.hpp"
#include "numeric/decimal.hpp"
#include "numeric/integer.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/rational.hpp"
#include "numeric/twoexp.hpp"
#include "utility/handle.hpp"
#include "utility/writable.hpp"

namespace Ariadne {
class Dyadic;
extern const Dyadic infty;
/** @brief Exception reporting infinity. */
class InfinityException : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};
// Exact finite value: coefficient / 2^denominator_exponent.
// Canonical form has no factor of two in the coefficient if exponent>0.
// GMP is deliberately absent from representation and arithmetic.
/** @brief Exact rational number whose denominator is a power of two. */
class Dyadic {
    static Writer<Dyadic> _default_writer;

  public:
    typedef ExactTag Paradigm;
    /** @brief Set the default writer. */
    static Void set_default_writer(Writer<Dyadic> writer);
    /** @brief Return the current numeric output-format policy. */
    static Writer<Dyadic> default_writer();
    /** @brief Convert to exact number. */
    operator ExactNumber() const;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream&, const Dyadic&);

    enum class Kind : unsigned char { Finite,
                                      PositiveInfinity,
                                      NegativeInfinity,
                                      NaN };

  private:
    Integer _coefficient;
    std::size_t _denominator_exponent = 0u;
    Kind _kind = Kind::Finite;

    /** @brief Validate and convert a nonnegative exponent. */
    static std::size_t checked_natural_exponent(const Natural& exp) {
        const BigUnsigned& magnitude = exp.magnitude();
        if (magnitude.bit_length() > std::numeric_limits<std::size_t>::digits) {
            throw std::overflow_error("Dyadic denominator exponent too large");
        }
        const std::size_t value = exp.template get<std::size_t>();
        // Match the original Natural overload's signed Int exponent range.
        if (!std::in_range<Int>(value)) {
            throw std::overflow_error("Dyadic denominator exponent too large");
        }
        return value;
    }
    /** @brief Return a value shifted by the supplied binary exponent. */
    static Integer shifted(const Integer& x, std::size_t bits) {
        if (bits == 0u || x.magnitude().is_zero()) {
            return x;
        }
        return Integer::from_magnitude(shift_left(x.magnitude(), bits), x.negative());
    }

    /** @brief Normalize the numeric representation. */
    void normalize() {
        // All callers construct a finite coefficient or check is_finite first.
        if (_coefficient.magnitude().is_zero()) {
            _denominator_exponent = 0u;
            return;
        }
        const auto& magnitude = _coefficient.magnitude();
        std::size_t trailing = 0u;
        // A nonzero canonical magnitude contains a nonzero limb, so this
        // scan always stops before the end of its storage.
        std::size_t i = 0u;
        while (magnitude.limb(i) == 0u) {
            trailing += 64u;
            ++i;
        }
        trailing += static_cast<std::size_t>(std::countr_zero(magnitude.limb(i)));
        const std::size_t reduce = std::min(trailing, _denominator_exponent);
        if (reduce != 0u) {
            _coefficient = Integer::from_magnitude(
                shift_right(magnitude, reduce), _coefficient.negative());
            _denominator_exponent -= reduce;
        }
    }

    /** @brief Return the truncated integer or dyadic representation. */
    static Integer truncated(const Dyadic& x) {
        return Integer::from_magnitude(
            shift_right(x._coefficient.magnitude(), x._denominator_exponent),
            x._coefficient.negative());
    }

    // Exact comparison of positive n/2^en and m/2^em without
    // constructing shifted operands.
    /** @brief Compare the values with the supplied binary scaling. */
    static int compare_scaled(const BigUnsigned& n, std::size_t en,
                              const BigUnsigned& m, std::size_t em) noexcept {
        if (n.is_zero() || m.is_zero()) {
            // cmp() has already handled two zeros; the owned arithmetic
            // caller requires at least one large, nonzero magnitude.
            return n.is_zero() ? -1 : 1;
        }
        const std::size_t nb = n.bit_length(), mb = m.bit_length();
        const bool ng = nb >= en, mg = mb >= em;
        if (ng != mg) {
            return ng ? 1 : -1;
        }
        if (ng) {
            if (nb - en != mb - em) {
                return nb - en > mb - em ? 1 : -1;
            }
        } else {
            if (en - nb != em - mb) {
                return en - nb < em - mb ? 1 : -1;
            }
        }
        return nb >= mb ? compare_aligned(n, m) : -compare_aligned(m, n);
    }

    // Align shorter positive magnitude virtually to the longer
    // magnitude's most significant bit. No owning BigUnsigned.
    /** @brief Compare the values after aligning their exponents. */
    static int compare_aligned(const BigUnsigned& a,
                               const BigUnsigned& b) noexcept {
        const std::size_t shift = a.bit_length() - b.bit_length();
        const std::size_t whole = shift / 64u;
        const unsigned part = static_cast<unsigned>(shift % 64u);
        for (std::size_t i = a.size(); i != 0u; --i) {
            const std::size_t index = i - 1u;
            std::uint64_t word = 0u;
            if (index >= whole) {
                const std::size_t j = index - whole;
                if (j < b.size()) {
                    word = part == 0u ? b.limb(j) : b.limb(j) << part;
                }
                // Aligning the leading bits gives j<=b.size(); with j>0
                // the carry source j-1 is therefore always in bounds.
                if (part != 0u && j != 0u) {
                    word |= b.limb(j - 1u) >> (64u - part);
                }
            }
            const std::uint64_t actual = a.limb(index);
            if (actual != word) {
                return actual < word ? -1 : 1;
            }
        }
        return 0;
    }

    /** @brief Reject a nonfinite operand. */
    static void require_finite(const Dyadic& x) {
        if (x._kind != Kind::Finite) {
            throw std::domain_error("non-finite Dyadic cannot be rounded");
        }
    }

    /** @brief Return the sign of the value. */
    static int signum(const Dyadic& x) noexcept {
        if (x._kind == Kind::PositiveInfinity) {
            return 1;
        }
        if (x._kind == Kind::NegativeInfinity) {
            return -1;
        }
        if (x._kind == Kind::NaN || x._coefficient.magnitude().is_zero()) {
            return 0;
        }
        return x._coefficient.negative() ? -1 : 1;
    }

  public:
    /** @brief Construct a \c Dyadic with the specified characteristics. */
    Dyadic() = default;
    /** @brief Construct a \c Dyadic from numerator, exponent. */
    Dyadic(Integer numerator, std::size_t exponent = 0u)
        : _coefficient(std::move(numerator)), _denominator_exponent(exponent) {
        normalize();
    }
    /** @brief Construct a \c Dyadic from numerator, exponent. */
    Dyadic(const Integer& numerator, Natural exponent)
        : Dyadic(numerator, checked_natural_exponent(exponent)) {}
    /** @brief Disallow this construction to preserve the ownership contract. */
    template<BuiltinSignedIntegral N>
    Dyadic(const Integer&, N) = delete;
    /** @brief Construct a \c Dyadic from n. */
    template<BuiltinIntegral N>
    Dyadic(N n) : Dyadic(Integer(n)) {}

    // Same value-conversion entry points as Dyadic; all arithmetic
    // remains native, and only exact dyadic Rational values are accepted.
    /** @brief Construct a \c Dyadic from q. */
    explicit Dyadic(const Rational& q) {
        if (is_nan(q)) {
            *this = nan();
            return;
        }
        if (is_inf(q)) {
            *this = inf(sgn(q) == Sign::POSITIVE ? 1 : -1);
            return;
        }
        const auto denominator = q.denominator();
        const auto& magnitude = denominator.magnitude();
        // A finite Rational has a nonzero denominator. Nonfinite kinds
        // were handled above, so bit_length()-1 is well-defined here.
        const std::size_t exp = magnitude.bit_length() - 1u;
        if (compare(magnitude, shift_left(BigUnsigned(1u), exp)) != 0) {
            throw std::invalid_argument("Rational is not a Dyadic");
        }
        _coefficient = q.numerator();
        _denominator_exponent = exp;
        normalize();
    }
    /** @brief Construct a \c Dyadic from d. */
    Dyadic(const ExactDouble& d) : Dyadic(Rational(d)) {}
    /** @brief Construct a \c Dyadic from d. */
    explicit Dyadic(double d) : Dyadic(Rational(ExactDouble(d))) {}
    /** @brief Construct a \c Dyadic from s. */
    explicit Dyadic(const String& s) : Dyadic(Rational(Decimal(s))) {}
    /** @brief Construct a \c Dyadic from d. */
    explicit Dyadic(const FloatDP& d);
    /** @brief Construct a \c Dyadic from d. */
    explicit Dyadic(const FloatMP& d);
    /** @brief Construct a \c Dyadic from power. */
    Dyadic(const TwoExp& power) {
        const Int e = power.exponent();
        if (e >= 0) {
            _coefficient = shifted(Integer(1), static_cast<std::size_t>(e));
        } else {
            // Avoid signed overflow for INT_MIN.
            _coefficient = Integer(1);
            _denominator_exponent = static_cast<std::size_t>(-(e + 1)) + 1u;
        }
    }
    /** @brief Return the integer mantissa. */
    Integer mantissa() const { return _coefficient; }
    /** @brief Return the exponent. */
    Int exponent() const {
        if (_denominator_exponent >
            static_cast<std::size_t>(std::numeric_limits<Int>::max())) {
            throw std::overflow_error("Dyadic exponent out of range");
        }
        return static_cast<Int>(_denominator_exponent);
    }
    /** @brief Return the literal representation. */
    String literal() const {
        if (_kind == Kind::NaN) {
            return String("NaN");
        }
        if (_kind == Kind::PositiveInfinity) {
            return String("inf");
        }
        if (_kind == Kind::NegativeInfinity) {
            return String("-inf");
        }
        if (_denominator_exponent == 0u) {
            return _coefficient.literal();
        }
        if (_denominator_exponent > std::numeric_limits<Nat>::max()) {
            throw std::overflow_error("Dyadic literal exponent out of range");
        }
        return Decimal(_coefficient * pow(Integer(5),
                                          static_cast<Nat>(_denominator_exponent)),
                       static_cast<Nat>(_denominator_exponent))
            .literal();
    }

    /** @brief Return an infinite value. */
    static Dyadic inf() { return inf(1); }
    /** @brief Return an infinite value. */
    static Dyadic inf(int sign) {
        Dyadic result;
        result._kind = sign == 0 ? Kind::NaN : sign > 0 ? Kind::PositiveInfinity
                                                        : Kind::NegativeInfinity;
        return result;
    }
    /** @brief Return a not-a-number value. */
    static Dyadic nan() { return inf(0); }
    /** @brief Return an infinite value. */
    static Dyadic inf(Sign sign) {
        return inf(sign == Sign::ZERO ? 0 : sign == Sign::POSITIVE ? 1
                                                                   : -1);
    }

    /** @brief Return the category discriminator. */
    Kind kind() const noexcept { return _kind; }
    /** @brief Return the selected coefficient. */
    const Integer& coefficient() const noexcept { return _coefficient; }
    /** @brief Return the power of two in the denominator. */
    std::size_t denominator_exponent() const noexcept {
        return _denominator_exponent;
    }
    /** @brief Test whether the value is finite. */
    bool is_finite() const noexcept { return _kind == Kind::Finite; }
    /** @brief Test whether the value is zero. */
    bool is_zero() const noexcept {
        return is_finite() && _coefficient.magnitude().is_zero();
    }

    /** @brief Convert to an exact rational. */
    Rational to_rational() const {
        if (_kind == Kind::NaN) {
            return Rational::nan();
        }
        if (_kind == Kind::PositiveInfinity) {
            return Rational::inf(Sign::POSITIVE);
        }
        if (_kind == Kind::NegativeInfinity) {
            return Rational::inf(Sign::NEGATIVE);
        }
        if (_denominator_exponent == 0u) {
            return Rational(_coefficient);
        }
        const Integer denominator = Integer::from_magnitude(
            shift_left(BigUnsigned(1u), _denominator_exponent), false);
        return Rational(_coefficient, denominator);
    }
    // Truncate to IEEE binary64 from the native magnitude.
    // GMP mpf_get_d returns positive zero on underflow even
    // for negative operands. No Rational allocation.
    /** @brief Convert to a built-in double. */
    double get_d() const {
        if (_kind == Kind::NaN) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (_kind == Kind::PositiveInfinity) {
            return std::numeric_limits<double>::infinity();
        }
        if (_kind == Kind::NegativeInfinity) {
            return -std::numeric_limits<double>::infinity();
        }
        const BigUnsigned& n = _coefficient.magnitude();
        if (n.is_zero()) {
            return 0.0;
        }
        const bool negative = _coefficient.negative();
        const std::size_t bits = n.bit_length();
        int exponent;
        if (bits > _denominator_exponent) {
            const std::size_t positive = bits - 1u - _denominator_exponent;
            if (positive > 1023u) {
                return negative ? -std::numeric_limits<double>::infinity()
                                : std::numeric_limits<double>::infinity();
            }
            exponent = static_cast<int>(positive);
        } else {
            const std::size_t negative_power = _denominator_exponent - bits + 1u;
            if (negative_power > 1074u) {
                return 0.0;
            }
            exponent = -static_cast<int>(negative_power);
        }
        std::uint64_t significand;
        if (bits <= 53u) {
            significand = n.limb(0u) << (53u - bits);
        } else {
            const std::size_t first = bits - 53u;
            const std::size_t index = first / 64u;
            const unsigned shift = static_cast<unsigned>(first % 64u);
            significand = n.limb(index) >> shift;
            if (shift != 0u && index + 1u < n.size()) {
                significand |= n.limb(index + 1u) << (64u - shift);
            }
            significand &= (std::uint64_t{1} << 53u) - 1u;
        }
        double value;
        if (exponent >= -1022) {
            value = std::ldexp(static_cast<double>(significand), exponent - 52);
        } else {
            const unsigned drop = static_cast<unsigned>(-1022 - exponent);
            value = std::ldexp(static_cast<double>(significand >> drop), -1074);
        }
        // Underflow returned above. The retained leading bit yields at
        // least one subnormal unit, so value is strictly positive here.
        return negative ? -value : value;
    }

    /** @brief Return the difference or additive inverse. */
    friend Dyadic operator-(const Dyadic& x) {
        if (x._kind == Kind::NaN) {
            return nan();
        }
        if (x._kind == Kind::PositiveInfinity) {
            return inf(-1);
        }
        if (x._kind == Kind::NegativeInfinity) {
            return inf(1);
        }
        return Dyadic(-x._coefficient, x._denominator_exponent);
    }
    /** @brief Return the sum or unary positive value. */
    friend Dyadic operator+(const Dyadic& x) { return x; }
    /** @brief Return the additive zero with the operand characteristics. */
    friend Dyadic nul(const Dyadic&) { return Dyadic(0); }
    /** @brief Return the unary positive value. */
    friend Dyadic pos(const Dyadic& x) { return x; }
    /** @brief Return the additive inverse. */
    friend Dyadic neg(const Dyadic& x) { return -x; }
    /** @brief Return the sum of the operands. */
    friend Dyadic add(const Dyadic& a, const Dyadic& b) {
        return a + b;
    }
    /** @brief Return the difference of the operands. */
    friend Dyadic sub(const Dyadic& a, const Dyadic& b) {
        return a - b;
    }
    /** @brief Return the product of the operands. */
    friend Dyadic mul(const Dyadic& a, const Dyadic& b) {
        return a * b;
    }
    /** @brief Compute the product of the first two operands plus the third. */
    friend Dyadic fma(const Dyadic& a, const Dyadic& b,
                      const Dyadic& c) { return a * b + c; }
    /** @brief Add the operand in place. */
    Dyadic& operator+=(const Dyadic& rhs) {
        *this = *this + rhs;
        return *this;
    }
    /** @brief Subtract the operand in place. */
    Dyadic& operator-=(const Dyadic& rhs) {
        *this = *this - rhs;
        return *this;
    }
    /** @brief Multiply by the operand in place. */
    Dyadic& operator*=(const Dyadic& rhs) {
        *this = *this * rhs;
        return *this;
    }

  private:
    // R36: on ARM64 the single-buffer accumulator improves add/sub
    // across 5..64-limb magnitudes; x86-64 benefits through 16 limbs
    // but regresses at 32..64 limbs. Keep a single exact algorithm,
    // changing only the architecture-specific crossover to fallback.
    // All other architectures conservatively use the shorter crossover.
#if defined(__aarch64__) || defined(_M_ARM64)
    static constexpr std::size_t aligned_add_sub_max_limbs = 64u;
#else
    static constexpr std::size_t aligned_add_sub_max_limbs = 16u;
#endif

    /** @brief Test whether aligned addition or subtraction is preferable for these operands. */
    static bool prefer_aligned_add_sub(const Dyadic& a,
                                       const Dyadic& b) noexcept {
        const std::size_t as = a._coefficient.magnitude().size();
        const std::size_t bs = b._coefficient.magnitude().size();
        const std::size_t limit = aligned_add_sub_max_limbs;
        const std::size_t max_exponent = std::max(
            a._denominator_exponent, b._denominator_exponent);
        const std::size_t min_exponent = std::min(
            a._denominator_exponent, b._denominator_exponent);
        return (as >= 5u || bs >= 5u) && as <= limit && bs <= limit &&
               max_exponent - min_exponent <= 63u;
    }

    // R38: use one owned magnitude allocation directly in BigUnsigned.
    // Large inputs bypass both the two-shift fallback and the R35
    // stack-buffer/owning-buffer double materialization. The equality
    // and sign decision use virtual magnitude comparison only.
    /** @brief Test whether owned shifted addition or subtraction is preferable for these operands. */
    static bool prefer_owned_shifted_add_sub(const Dyadic& a,
                                             const Dyadic& b) noexcept {
#if defined(__aarch64__) || defined(_M_ARM64)
        constexpr std::size_t min_limbs = 32u;
#else
        constexpr std::size_t min_limbs = 17u;
#endif
        const std::size_t largest = std::max(
            a._coefficient.magnitude().size(),
            b._coefficient.magnitude().size());
        const std::size_t high = std::max(
            a._denominator_exponent, b._denominator_exponent);
        const std::size_t low = std::min(
            a._denominator_exponent, b._denominator_exponent);
        return largest >= min_limbs && high - low <= 63u;
    }

    /** @brief Add or subtract using owned shifted integer storage. */
    static Dyadic owned_shifted_add_sub(
        const Dyadic& a, const Dyadic& b, bool subtract_second) {
        const BigUnsigned& an = a._coefficient.magnitude();
        const BigUnsigned& bn = b._coefficient.magnitude();
        const std::size_t exponent = std::max(
            a._denominator_exponent, b._denominator_exponent);
        const unsigned shift_a = static_cast<unsigned>(
            exponent - a._denominator_exponent);
        const unsigned shift_b = static_cast<unsigned>(
            exponent - b._denominator_exponent);
        const bool negative_a = a._coefficient.negative();
        const bool negative_b = b._coefficient.negative() != subtract_second;
        const bool same_sign = negative_a == negative_b;

        BigUnsigned magnitude;
        bool result_negative = negative_a;
        if (same_sign) {
            magnitude = shift_a == 0u && shift_b == 0u
                          ? add(an, bn)
                          : BigUnsigned::add_shifted(
                                an.limbs(), shift_a, bn.limbs(), shift_b);
        } else {
            const int order = compare_scaled(
                an, a._denominator_exponent,
                bn, b._denominator_exponent);
            if (order == 0) {
                return Dyadic();
            }
            if (order > 0) {
                magnitude = shift_a == 0u && shift_b == 0u
                              ? sub(an, bn)
                              : BigUnsigned::sub_shifted(
                                    an.limbs(), shift_a, bn.limbs(), shift_b);
            } else {
                magnitude = shift_a == 0u && shift_b == 0u
                              ? sub(bn, an)
                              : BigUnsigned::sub_shifted(
                                    bn.limbs(), shift_b, an.limbs(), shift_a);
                result_negative = negative_b;
            }
        }
        return Dyadic(Integer::from_magnitude(
                          std::move(magnitude), result_negative),
                      exponent);
    }

    // R35: for aligned finite signed arithmetic on 5..64-limb
    // magnitudes with an exponent gap <=63, scan borrowed virtual
    // shifted limbs and produce one canonical owned magnitude.
    // No shifted Integer operands or temporary negation.
    /** @brief Add or subtract after aligning the dyadic exponents. */
    static Dyadic aligned_add_sub(const Dyadic& a, const Dyadic& b,
                                  bool subtract_second) {
        using Limb = BigUnsigned::Limb;
        const BigUnsigned& an = a._coefficient.magnitude();
        const BigUnsigned& bn = b._coefficient.magnitude();
        if (an.is_zero()) {
            return subtract_second ? -b : b;
        }
        if (bn.is_zero()) {
            return a;
        }

        const std::size_t exponent = std::max(a._denominator_exponent,
                                              b._denominator_exponent);
        const unsigned shift_a = static_cast<unsigned>(
            exponent - a._denominator_exponent);
        const unsigned shift_b = static_cast<unsigned>(
            exponent - b._denominator_exponent);
        const std::size_t a_bits = an.bit_length() + shift_a;
        const std::size_t b_bits = bn.bit_length() + shift_b;
        const std::size_t width = (std::max(a_bits, b_bits) + 63u) / 64u;
        // Up to 65 shifted limbs, plus one final addition carry.
        std::array<Limb, 66u> output{};

        const auto digit = [](const BigUnsigned& n, std::size_t i,
                              unsigned shift) -> Limb {
            if (shift == 0u) {
                return n.limb(i);
            }
            const Limb low = n.limb(i) << shift;
            return i == 0u ? low : low | (n.limb(i - 1u) >> (64u - shift));
        };

        const bool a_negative = a._coefficient.negative();
        const bool b_negative = b._coefficient.negative() != subtract_second;
        const bool adding = a_negative == b_negative;
        bool negative = a_negative, reverse = false;
        if (!adding) {
            int order = 0;
            if (a_bits != b_bits) {
                order = a_bits > b_bits ? 1 : -1;
            } else {
                for (std::size_t i = width; i != 0u; --i) {
                    const Limb av = digit(an, i - 1u, shift_a);
                    const Limb bv = digit(bn, i - 1u, shift_b);
                    if (av != bv) {
                        order = av > bv ? 1 : -1;
                        break;
                    }
                }
            }
            if (order == 0) {
                return Dyadic();
            }
            reverse = order < 0;
            negative = reverse ? b_negative : a_negative;
        }

        Limb carry_or_borrow = 0u;
        for (std::size_t i = 0u; i < width; ++i) {
            Limb av = digit(an, i, shift_a), bv = digit(bn, i, shift_b);
            if (reverse) {
                std::swap(av, bv);
            }
            if (adding) {
                const Limb partial = av + bv;
                const Limb result = partial + carry_or_borrow;
                carry_or_borrow = static_cast<Limb>(
                    partial < av || result < partial);
                output[i] = result;
            } else {
                const Limb partial = av - bv;
                const Limb result = partial - carry_or_borrow;
                carry_or_borrow = static_cast<Limb>(
                    av < bv || partial < carry_or_borrow);
                output[i] = result;
            }
        }
        std::size_t used = width;
        if (adding && carry_or_borrow != 0u) {
            output[used++] = carry_or_borrow;
        }
        return Dyadic(Integer::from_magnitude(
                          BigUnsigned(BigUnsigned::LimbView(output.data(), used)), negative),
                      exponent);
    }

  public:
    /** @brief Return the sum or unary positive value. */
    friend Dyadic operator+(const Dyadic& a, const Dyadic& b) {
        if (!a.is_finite() || !b.is_finite()) {
            if (a._kind == Kind::NaN || b._kind == Kind::NaN) {
                return nan();
            }
            if (!a.is_finite() && !b.is_finite()) {
                return a._kind == b._kind ? a : nan();
            }
            return a.is_finite() ? b : a;
        }
        const std::size_t e = std::max(a._denominator_exponent,
                                       b._denominator_exponent);
        if (prefer_owned_shifted_add_sub(a, b)) {
            return owned_shifted_add_sub(a, b, false);
        }
        if (prefer_aligned_add_sub(a, b)) {
            return aligned_add_sub(a, b, false);
        }
        return Dyadic(shifted(a._coefficient, e - a._denominator_exponent) + shifted(b._coefficient, e - b._denominator_exponent), e);
    }
    /** @brief Return the difference or additive inverse. */
    friend Dyadic operator-(const Dyadic& a, const Dyadic& b) {
        if (a.is_finite() && b.is_finite()) {
            if (prefer_owned_shifted_add_sub(a, b)) {
                return owned_shifted_add_sub(a, b, true);
            }
            if (prefer_aligned_add_sub(a, b)) {
                return aligned_add_sub(a, b, true);
            }
        }
        return a + (-b);
    }
    /** @brief Return the product. */
    friend Dyadic operator*(const Dyadic& a, const Dyadic& b) {
        if (!a.is_finite() || !b.is_finite()) {
            if (a._kind == Kind::NaN || b._kind == Kind::NaN ||
                signum(a) == 0 || signum(b) == 0) {
                return nan();
            }
            return inf(signum(a) == signum(b) ? 1 : -1);
        }
        if (a._denominator_exponent >
            std::numeric_limits<std::size_t>::max() - b._denominator_exponent) {
            throw std::overflow_error("Dyadic binary exponent overflow");
        }
        return Dyadic(a._coefficient * b._coefficient,
                      a._denominator_exponent + b._denominator_exponent);
    }
    /** @brief Return the quotient. */
    friend Dyadic operator/(Dyadic x, int denominator) {
        if (denominator != 2) {
            throw std::domain_error("only divide by two is a Dyadic operation");
        }
        if (x.is_finite() && !x.is_zero()) {
            if (x._denominator_exponent == std::numeric_limits<std::size_t>::max()) {
                throw std::overflow_error("Dyadic binary exponent overflow");
            }
            ++x._denominator_exponent;
            x.normalize();
        }
        return x;
    }
    /** @brief Divide the value by two. */
    friend Dyadic hlf(const Dyadic& x) { return x / 2; }
    /** @brief Return the quotient. */
    friend Rational operator/(const Dyadic& a, const Dyadic& b) {
        return div(a.to_rational(), b.to_rational());
    }
    // Quotients and reciprocals of dyadics are exact rational numbers.
    /** @brief Return the multiplicative inverse. */
    friend Rational rec(const Dyadic& x) {
        return rec(x.to_rational());
    }
    /** @brief Return the quotient of the operands. */
    friend Rational div(const Dyadic& a, const Dyadic& b) {
        return div(a.to_rational(), b.to_rational());
    }
    /** @brief Return the quotient of the operands. */
    friend Rational div(const Dyadic& a, const Integer& b) {
        return div(a.to_rational(), Rational(b));
    }
    /** @brief Return the quotient of the operands. */
    friend Rational div(const Integer& a, const Dyadic& b) {
        return div(Rational(a), b.to_rational());
    }

    /** @brief Raise the value to the specified power. */
    friend Dyadic pow(const Dyadic& x, Nat exponent) {
        Dyadic result(1);
        Dyadic factor(x);
        for (unsigned n = exponent; n != 0u; n >>= 1u) {
            if ((n & 1u) != 0u) {
                result = result * factor;
            }
            if (n > 1u) {
                factor = factor * factor;
            }
        }
        return result;
    }
    /** @brief Test equality of the operands. */
    friend Boolean eq(const Dyadic& a, const Dyadic& b) {
        if (a._kind == Kind::NaN || b._kind == Kind::NaN) {
            return false;
        }
        if (a._kind != b._kind) {
            return false;
        }
        if (!a.is_finite()) {
            return true;
        }
        return a._denominator_exponent == b._denominator_exponent &&
               static_cast<bool>(a._coefficient == b._coefficient);
    }
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const Dyadic& a, const Dyadic& b) {
        if (a._kind == Kind::NaN || b._kind == Kind::NaN) {
            return Comparison::INCOMPARABLE;
        }
        if (a._kind == b._kind && !a.is_finite()) {
            return Comparison::EQUAL;
        }
        if (a._kind == Kind::NegativeInfinity ||
            b._kind == Kind::PositiveInfinity) {
            return Comparison::LESS;
        }
        if (a._kind == Kind::PositiveInfinity ||
            b._kind == Kind::NegativeInfinity) {
            return Comparison::GREATER;
        }
        const int sa = signum(a), sb = signum(b);
        if (sa != sb) {
            return sa < sb ? Comparison::LESS : Comparison::GREATER;
        }
        if (sa == 0) {
            return Comparison::EQUAL;
        }
        const int magnitude = compare_scaled(
            a._coefficient.magnitude(), a._denominator_exponent,
            b._coefficient.magnitude(), b._denominator_exponent);
        const int order = sa < 0 ? -magnitude : magnitude;
        return order < 0 ? Comparison::LESS : order > 0 ? Comparison::GREATER
                                                        : Comparison::EQUAL;
    }
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const Dyadic& a, Int b) { return cmp(a, Dyadic(b)); }
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(Int a, const Dyadic& b) { return cmp(Dyadic(a), b); }
    /** @brief Test whether the first operand is less than the second. */
    friend Boolean lt(const Dyadic& a, const Dyadic& b) {
        if (a._kind == Kind::NaN || b._kind == Kind::NaN) {
            return false;
        }
        return cmp(a, b) == Comparison::LESS;
    }
    /** @brief Return the sign of the value. */
    friend Sign sgn(const Dyadic& x) {
        return signum(x) > 0 ? Sign::POSITIVE : signum(x) < 0 ? Sign::NEGATIVE
                                                              : Sign::ZERO;
    }
    /** @brief Test whether the value is not a number. */
    friend Bool is_nan(const Dyadic& x) {
        return x._kind == Kind::NaN;
    }
    /** @brief Test whether the value is infinite. */
    friend Bool is_inf(const Dyadic& x) {
        return x._kind == Kind::PositiveInfinity ||
               x._kind == Kind::NegativeInfinity;
    }
    /** @brief Test whether the value is finite. */
    friend Bool is_finite(const Dyadic& x) { return x.is_finite(); }
    /** @brief Test whether the value is zero. */
    friend Bool is_zero(const Dyadic& x) { return x.is_zero(); }
    /** @brief Test equality of the operands. */
    friend Boolean operator==(const Dyadic& a, const Dyadic& b) { return eq(a, b); }
    /** @brief Test inequality of the operands. */
    friend Boolean operator!=(const Dyadic& a, const Dyadic& b) { return !eq(a, b); }

    /** @brief Test whether the first operand is less than the second. */
    friend Boolean operator<(const Dyadic& a, const Dyadic& b) { return lt(a, b); }
    /** @brief Test whether the first operand is greater than the second. */
    friend Boolean operator>(const Dyadic& a, const Dyadic& b) { return lt(b, a); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend Boolean operator<=(const Dyadic& a, const Dyadic& b) { return !lt(b, a); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend Boolean operator>=(const Dyadic& a, const Dyadic& b) { return !lt(a, b); }
    /** @brief Raise the value to the specified power. */
    friend Rational pow(const Dyadic& x, Int exponent) {
        if (exponent >= 0) {
            return pow(x, static_cast<Nat>(exponent)).to_rational();
        }
        const Nat positive = static_cast<Nat>(-(exponent + 1)) + Nat(1);
        return Rational(1) / pow(x, positive).to_rational();
    }
    /** @brief Return the minimum of the operands. */
    friend Dyadic min(const Dyadic& a, const Dyadic& b) {
        if (a._kind == Kind::NaN || b._kind == Kind::NaN) {
            return nan();
        }
        return b < a ? b : a;
    }
    /** @brief Return the maximum of the operands. */
    friend Dyadic max(const Dyadic& a, const Dyadic& b) {
        if (a._kind == Kind::NaN || b._kind == Kind::NaN) {
            return nan();
        }
        return a < b ? b : a;
    }
    /** @brief Return the absolute value. */
    friend Positive<Dyadic> abs(const Dyadic& x);
    /** @brief Return the square. */
    friend Positive<Dyadic> sqr(const Dyadic& x);
    /** @brief Return the greatest integer not exceeding the value. */
    friend Integer floor(const Dyadic& x) {
        require_finite(x);
        Integer z = truncated(x);
        if (x._coefficient.negative() &&
            x._denominator_exponent != 0u) {
            z = z - Integer(1);
        }
        return z;
    }
    /** @brief Return the least integer not less than the value. */
    friend Integer ceil(const Dyadic& x) {
        require_finite(x);
        Integer z = truncated(x);
        if (!x._coefficient.negative() &&
            x._denominator_exponent != 0u) {
            z = z + Integer(1);
        }
        return z;
    }
    /** @brief Round the value to an integer. */
    friend Integer round(const Dyadic& x) {
        require_finite(x);
        const auto& magnitude = x._coefficient.magnitude();
        BigUnsigned integer = shift_right(magnitude, x._denominator_exponent);
        // Round ties away from zero, matching current Dyadic.
        if (x._denominator_exponent != 0u) {
            const auto half = shift_right(magnitude, x._denominator_exponent - 1u);
            if ((half.limb(0u) & 1u) != 0u) {
                integer = add(integer, BigUnsigned(1u));
            }
        }
        return Integer::from_magnitude(std::move(integer), x._coefficient.negative());
    }
};

/** @brief Output policy formatting numbers in decimal notation. */
class DecimalWriter : public WriterInterface<Dyadic> {
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os, const Dyadic& x) const final override;
};
/** @brief Output policy formatting numbers in scientific notation. */
class ScientificWriter : public WriterInterface<Dyadic> {
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os, const Dyadic& x) const final override;
};
/** @brief Output policy formatting exact numbers as fractions. */
class FractionWriter : public WriterInterface<Dyadic> {
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os, const Dyadic& x) const final override;
};
/** @brief Output policy selecting the stored value representation. */
template<>
class RepresentationWriter<Dyadic> : public WriterInterface<Dyadic> {
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os, const Dyadic& x) const final override;
};
/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<>
class Positive<Dyadic> : public Dyadic {
  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : Dyadic() {}
    /** @brief Construct a \c Positive from n. */
    template<BuiltinUnsignedIntegral M>
    Positive(M n) : Dyadic(n) {}
    /** @brief Disallow this construction to preserve the ownership contract. */
    Positive(int) = delete;
    /** @brief Construct a \c Positive from x. */
    explicit Positive(const Dyadic& x) : Dyadic(x) {
        if (sgn(x) == Sign::NEGATIVE) {
            throw std::domain_error("negative Positive<Dyadic>");
        }
    }
    /** @brief Move-construct the object from the supplied instance. */
    explicit Positive(Dyadic&& x) : Dyadic(std::move(x)) {
        if (sgn(*this) == Sign::NEGATIVE) {
            throw std::domain_error("negative Positive<Dyadic>");
        }
    }
};
/** @brief Convert to the positive representation. */
inline Positive<Dyadic> cast_positive(const Dyadic& x) { return Positive<Dyadic>(x); }
/** @brief Return the absolute value. */
inline Positive<Dyadic> abs(const Dyadic& x) {
    // D42: copy the canonical magnitude exactly once, without forming
    // a negated Dyadic or renormalizing an already canonical value.
    Dyadic result;
    result._coefficient = Integer::from_magnitude(
        x._coefficient.magnitude(), false);
    result._denominator_exponent = x._denominator_exponent;
    result._kind = x._kind == Dyadic::Kind::NegativeInfinity
                     ? Dyadic::Kind::PositiveInfinity
                     : x._kind;
    return Positive<Dyadic>(std::move(result));
}
/** @brief Return the square. */
inline Positive<Dyadic> sqr(const Dyadic& x) {
    return Positive<Dyadic>(x * x);
}

// Existing native core implementation ends here. Wrapper specializations
// live in this same header, matching the final numeric/dyadic.hpp layout.
} // namespace Ariadne

#include "foundation/logical.hpp"
#include "numeric/arithmetic.hpp"
#include "utility/handle.hpp"
#include "utility/writable.hpp"

namespace Ariadne {

using PositiveDyadic = Positive<Dyadic>;

//! \relates Dyadic
//! \name Type synonyms
//!@{
//! \ingroup NumericModule
using DyadicBounds = Bounds<Dyadic>;               //!< Alias for dyadic bounds on a number.
using DyadicBall = Ball<Dyadic, Dyadic>;           //!< Alias for ball about a number with dyadic value and error.
using DyadicUpperBound = UpperBound<Dyadic>;       //!< Alias for dyadic upper bound for a number.
using DyadicLowerBound = LowerBound<Dyadic>;       //!< Alias for dyadic lower bound for a number.
using DyadicApproximation = Approximation<Dyadic>; //!< Alias for dyadic approximateion to a number.

using PositiveDyadicApproximation = Positive<Approximation<Dyadic>>; //!< <p/>
using PositiveDyadicLowerBound = Positive<LowerBound<Dyadic>>;       //!< <p/>
using PositiveDyadicUpperBound = Positive<UpperBound<Dyadic>>;       //!< <p/>
using PositiveDyadicBounds = Positive<Bounds<Dyadic>>;               //!< <p/>
//!@}

/** @brief Validated numeric enclosure with lower and upper endpoints. */
template<>
class Bounds<Dyadic> {
    Dyadic _l, _u;

  public:
    typedef ValidatedTag Paradigm;
    /** @brief Construct a \c Bounds from w. */
    Bounds(Dyadic w) : _l(w), _u(w) {}
    /** @brief Construct a \c Bounds from x. */
    template<ConvertibleTo<Dyadic> X>
    Bounds(const X& x)
        : Bounds(Dyadic(x)) {}
    /** @brief Construct a \c Bounds from l, u. */
    Bounds(Dyadic l, Dyadic u) : _l(l), _u(u) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<class X>
        requires Constructible<Dyadic, X>
    Bounds(const Bounds<X>& x)
        : Bounds(Dyadic(x.lower_raw()), Dyadic(x.upper_raw())) {}
    /** @brief Convert to validated number. */
    operator ValidatedNumber() const;
    /** @brief Return an enclosure enlarged by the supplied symmetric error. */
    Bounds<Dyadic> pm(Dyadic e) { return Bounds<Dyadic>(_l - e, _u + e); }
    /** @brief Return the lower endpoint. */
    Dyadic lower() const { return _l; }
    /** @brief Return the upper endpoint. */
    Dyadic upper() const { return _u; }
    /** @brief Return the underlying lower-endpoint representation. */
    Dyadic lower_raw() const { return _l; }
    /** @brief Return the underlying upper-endpoint representation. */
    Dyadic upper_raw() const { return _u; }
    /** @brief Return the requested element or stored value. */
    Bounds<FloatDP> get(DoublePrecision pr) const;
    /** @brief Return the requested element or stored value. */
    Bounds<FloatMP> get(MultiplePrecision pr) const;
    /** @brief Return the sum or unary positive value. */
    friend Bounds<Dyadic> operator+(const Bounds<Dyadic>& w) { return Bounds<Dyadic>(+w._l, w._u); }
    /** @brief Return the difference or additive inverse. */
    friend Bounds<Dyadic> operator-(const Bounds<Dyadic>& w) { return Bounds<Dyadic>(-w._u, -w._l); }
    /** @brief Return the sum or unary positive value. */
    friend Bounds<Dyadic> operator+(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        return Bounds<Dyadic>(w1._l + w2._l, w1._u + w2._u);
    }
    /** @brief Return the difference or additive inverse. */
    friend Bounds<Dyadic> operator-(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        return Bounds<Dyadic>(w1._l - w2._u, w1._u - w2._l);
    }
    /** @brief Return the product. */
    friend Bounds<Dyadic> operator*(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        // Exact multiplication enclosure: all four vertices, no rounding.
        Dyadic p = w1._l * w2._l, q = w1._l * w2._u;
        Dyadic r = w1._u * w2._l, s = w1._u * w2._u;
        return Bounds<Dyadic>(min(min(p, q), min(r, s)),
                              max(max(p, q), max(r, s)));
    }
    /** @brief Return the quotient. */
    friend Bounds<Rational> operator/(const Bounds<Rational>& w1, const Bounds<Rational>& w2);
    /** @brief Return the sum of the operands. */
    friend Bounds<Dyadic> add(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        return Bounds<Dyadic>(w1._l + w2._l, w1._u + w2._u);
    }
    /** @brief Return the difference of the operands. */
    friend Bounds<Dyadic> sub(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        return Bounds<Dyadic>(w1._l - w2._u, w1._u - w2._l);
    }
    /** @brief Return the product of the operands. */
    friend Bounds<Dyadic> mul(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) { return w1 * w2; }
    /** @brief Return the quotient of the operands. */
    friend Bounds<Rational> div(const Bounds<Rational>& w1, const Bounds<Rational>& w2);
    /** @brief Raise the value to the specified power. */
    friend Bounds<Dyadic> pow(const Bounds<Dyadic>& w, Nat m) {
        if (m == 0u) {
            return Bounds<Dyadic>(Dyadic(1));
        }
        if (m % 2u != 0u) {
            return Bounds<Dyadic>(pow(w._l, m), pow(w._u, m));
        }
        if (w._l >= 0) {
            return Bounds<Dyadic>(pow(w._l, m), pow(w._u, m));
        }
        if (w._u <= 0) {
            return Bounds<Dyadic>(pow(w._u, m), pow(w._l, m));
        }
        return Bounds<Dyadic>(Dyadic(0),
                              max(pow(w._l, m), pow(w._u, m)));
    }
    /** @brief Raise the value to the specified power. */
    friend Bounds<Dyadic> pow(const Bounds<Dyadic>& w, Int m) {
        if (m < 0) {
            throw std::domain_error("negative Dyadic interval power");
        }
        return pow(w, static_cast<Nat>(m));
    }
    /** @brief Return the additive zero with the operand characteristics. */
    friend Bounds<Dyadic> nul(const Bounds<Dyadic>& w) {
        return Bounds<Dyadic>(nul(w._l), nul(w._u));
    }
    /** @brief Return the unary positive value. */
    friend Bounds<Dyadic> pos(const Bounds<Dyadic>& w) {
        return Bounds<Dyadic>(pos(w._l), pos(w._u));
    }
    /** @brief Return the additive inverse. */
    friend Bounds<Dyadic> neg(const Bounds<Dyadic>& w) {
        return Bounds<Dyadic>(neg(w._u), neg(w._l));
    }
    /** @brief Return the square. */
    friend Bounds<Dyadic> sqr(const Bounds<Dyadic>& w) {
        if (w._l > 0) {
            return Bounds<Dyadic>(sqr(w._l), sqr(w._u));
        } else if (w._u < 0) {
            return Bounds<Dyadic>(sqr(w._u), sqr(w._l));
        } else {
            return Bounds<Dyadic>(nul(w._l), max(sqr(w._l), sqr(w._u)));
        }
    }
    /** @brief Divide the value by two. */
    friend Bounds<Dyadic> hlf(const Bounds<Dyadic>& w) { return Bounds<Dyadic>(hlf(w._l), hlf(w._u)); }
    /** @brief Return the multiplicative inverse. */
    friend Bounds<Rational> rec(const Bounds<Rational>& w);
    /** @brief Return the absolute value. */
    friend Bounds<Dyadic> abs(const Bounds<Dyadic>& w) { return Bounds<Dyadic>(max(min(w._l, -w._u), 0), max(-w._l, w._u)); }
    /** @brief Return the maximum of the operands. */
    friend Bounds<Dyadic> max(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) { return Bounds<Dyadic>(max(w1._l, w2._l), max(w1._u, w2._u)); }
    /** @brief Return the minimum of the operands. */
    friend Bounds<Dyadic> min(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) { return Bounds<Dyadic>(min(w1._l, w2._l), min(w1._u, w2._u)); }

    /** @brief Return the sign of the value. */
    friend ValidatedKleenean sgn(const Bounds<Dyadic>& w) {
        if (w._l > 0) {
            return true;
        } else if (w._u < 0) {
            return false;
        } else {
            return indeterminate;
        }
    }
    /** @brief Test equality of the operands. */
    friend ValidatedKleenean operator==(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        if (w1._l >= w2._u && w1._u <= w2._l) {
            return true;
        } else if (w1._u < w2._l || w1._l > w2._u) {
            return false;
        } else {
            return indeterminate;
        }
    }
    /** @brief Test inequality of the operands. */
    friend ValidatedKleenean operator!=(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        if (w1._u < w2._l || w1._l > w2._u) {
            return true;
        } else if (w1._l >= w2._u && w1._u <= w2._l) {
            return false;
        } else {
            return indeterminate;
        }
    }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedKleenean operator<=(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        if (w1._u <= w2._l) {
            return true;
        } else if (w1._l > w2._u) {
            return false;
        } else {
            return indeterminate;
        }
    }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedKleenean operator>=(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        if (w1._l >= w2._u) {
            return true;
        } else if (w1._u < w2._l) {
            return false;
        } else {
            return indeterminate;
        }
    }
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        if (w1._u < w2._l) {
            return true;
        } else if (w1._l >= w2._u) {
            return false;
        } else {
            return indeterminate;
        }
    }
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        if (w1._l > w2._u) {
            return true;
        } else if (w1._u <= w2._l) {
            return false;
        } else {
            return indeterminate;
        }
    }
    /** @brief Test whether the enclosures admit a common value. */
    friend Bool consistent(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        return w1._l <= w2._u && w1._u >= w2._l;
    }
    /** @brief Test whether the enclosures cannot represent a common value. */
    friend Boolean inconsistent(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        return w1._l > w2._u || w1._u < w2._l;
    }
    /** @brief Test whether the first enclosure is contained in the second. */
    friend Boolean refines(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        return w1._l >= w2._l and w1._u <= w2._u;
    }
    /** @brief Return the common refinement of the enclosures. */
    friend Bounds<Dyadic> refinement(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        return Bounds<Dyadic>(max(w1._l, w2._l), min(w1._u, w2._u));
    }
    /** @brief Return an enclosure containing both operands. */
    friend Bounds<Dyadic> coarsening(const Bounds<Dyadic>& w1, const Bounds<Dyadic>& w2) {
        return Bounds<Dyadic>(min(w1._l, w2._l), max(w1._u, w2._u));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, Bounds<Dyadic> y) { return os << "[" << y._l << ":" << y._u << "]"; }
};

/** @brief Validated numeric enclosure represented by a centre and error radius. */
template<>
class Ball<Dyadic, Dyadic> {
    Dyadic _v, _e;

  public:
    /** @brief Construct a \c Ball from w. */
    Ball(Dyadic w) : _v(w), _e(0) {}
    /** @brief Construct a \c Ball from v, e. */
    Ball(Dyadic v, Dyadic e) : _v(v), _e(e) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<class X, class XE>
        requires Constructible<Dyadic, X> and Constructible<Dyadic, XE>
    Ball(const Ball<X, XE>& x)
        : Ball(Dyadic(x.value_raw()), Dyadic(x.error_raw())) {}
    /** @brief Construct a \c Ball from w. */
    explicit Ball(const Bounds<Dyadic>& w) : _v(hlf(w.lower() + w.upper())), _e(hlf(w.upper() - w.lower())) {}
    /** @brief Convert to bounds< dyadic >. */
    explicit operator Bounds<Dyadic>() const { return Bounds<Dyadic>(_v - _e, _v + _e); }
    /** @brief Return the stored value. */
    Dyadic value() const { return _v; }
    /** @brief Return the error bound. */
    Dyadic error() const { return _e; }
    /** @brief Return the underlying value representation. */
    Dyadic value_raw() const { return _v; }
    /** @brief Return the underlying error representation. */
    Dyadic error_raw() const { return _e; }
    /** @brief Return the sum or unary positive value. */
    friend DyadicBall operator+(const DyadicBall& w1, const DyadicBall& w2) { return DyadicBall(w1._v + w2._v, w1._e + w2._e); }
    /** @brief Return the difference or additive inverse. */
    friend DyadicBall operator-(const DyadicBall& w1, const DyadicBall& w2) { return DyadicBall(w1._v - w2._v, w1._e + w2._e); }
    /** @brief Return the product. */
    friend DyadicBall operator*(const DyadicBall& w1, const DyadicBall& w2) { return DyadicBall(w1._v * w2._v, abs(w1._v) * w2._e + w1._e * abs(w2._v) + w1._e * w2._e); }
    /** @brief Return the absolute value. */
    friend DyadicBall abs(const DyadicBall& w) {
        if (abs(w._v) >= w._e) {
            return DyadicBall(abs(w._v), w._e);
        } else {
            Dyadic av = hlf(max(w._e - w._v, w._v + w._e));
            return DyadicBall(av, av);
        }
    }
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const DyadicBall& w1, const DyadicBall& w2) {
        if (w1._v + w1._e < w2._v - w2._e) {
            return true;
        } else if (w1._v - w1._e >= w2._v + w2._e) {
            return false;
        } else {
            return indeterminate;
        }
    }
    /** @brief Test whether the enclosures cannot represent a common value. */
    friend Boolean inconsistent(const DyadicBall& w1, const DyadicBall& w2) { return abs(w1._v - w2._v) > w1._e + w2._e; }
    /** @brief Test whether the first enclosure is contained in the second. */
    friend Boolean refines(const DyadicBall& w1, const DyadicBall& w2) { return abs(w1._v - w2._v) + w1._e <= w2._e; }
    /** @brief Return the common refinement of the enclosures. */
    friend DyadicBall refinement(const DyadicBall& w1, const DyadicBall& w2) { return DyadicBall(refinement(DyadicBounds(w1), DyadicBounds(w2))); }
    /** @brief Return an enclosure containing both operands. */
    friend DyadicBall coarsening(const DyadicBall& w1, const DyadicBall& w2) { return DyadicBall(coarsening(DyadicBounds(w1), DyadicBounds(w2))); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, DyadicBall y) { return os << "[" << y._v << ":" << y._e << "]"; }
};

/** @brief One-sided lower bound on a numeric value. */
template<>
class LowerBound<Dyadic> {
    Dyadic _l;

  public:
    /** @brief Construct a \c LowerBound from l. */
    LowerBound(Dyadic l) : _l(l) {}
    /** @brief Construct a \c LowerBound from lu. */
    LowerBound(Bounds<Dyadic> lu) : _l(lu.lower_raw()) {}
    /** @brief Copy-construct the object from the supplied instance. */
    LowerBound(const LowerBound<FloatDP>& x);
    /** @brief Copy-construct the object from the supplied instance. */
    LowerBound(const LowerBound<FloatMP>& x);
    /** @brief Return the underlying numeric representation. */
    Dyadic raw() const { return _l; }
    /** @brief Return the requested element or stored value. */
    LowerBound<FloatDP> get(DoublePrecision pr) const;
    /** @brief Return the requested element or stored value. */
    LowerBound<FloatMP> get(MultiplePrecision pr) const;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const LowerBound<Dyadic>& y) { return os << y._l << ":"; }
};

/** @brief One-sided upper bound on a numeric value. */
template<>
class UpperBound<Dyadic> {
    Dyadic _u;

  public:
    /** @brief Construct a \c UpperBound from u. */
    explicit UpperBound(Dyadic u) : _u(u) {}
    /** @brief Construct a \c UpperBound from lu. */
    UpperBound(Bounds<Dyadic> lu) : _u(lu.upper_raw()) {}
    /** @brief Copy-construct the object from the supplied instance. */
    UpperBound(const UpperBound<FloatDP>& x);
    /** @brief Copy-construct the object from the supplied instance. */
    UpperBound(const UpperBound<FloatMP>& x);
    /** @brief Return the underlying numeric representation. */
    Dyadic raw() const { return _u; }
    /** @brief Return the requested element or stored value. */
    UpperBound<FloatDP> get(DoublePrecision pr) const;
    /** @brief Return the requested element or stored value. */
    UpperBound<FloatMP> get(MultiplePrecision pr) const;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const UpperBound<Dyadic>& y) { return os << ":" << y._u; }
};

/** @brief Numeric approximation without a guaranteed error enclosure. */
template<>
class Approximation<Dyadic> {
    Dyadic _a;

  public:
    /** @brief Construct a \c Approximation from a. */
    Approximation(Dyadic a) : _a(a) {}
    /** @brief Copy-construct the object from the supplied instance. */
    Approximation(const Approximation<FloatDP>& x);
    /** @brief Return the underlying numeric representation. */
    Dyadic raw() const { return _a; }
    /** @brief Return the requested element or stored value. */
    Approximation<FloatDP> get(DoublePrecision pr) const;
    /** @brief Return the requested element or stored value. */
    Approximation<FloatMP> get(MultiplePrecision pr) const;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Approximation<Dyadic>& y) {
        return os << y._a;
    }
};

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<>
class Positive<Bounds<Dyadic>> : public Bounds<Dyadic> {
  public:
    /** @brief Construct a \c Positive from w. */
    Positive(Bounds<Dyadic> w) : Bounds<Dyadic>(w) {};
};
/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<>
class Positive<LowerBound<Dyadic>> : public LowerBound<Dyadic> {
  public:
    /** @brief Construct a \c Positive from w. */
    Positive(LowerBound<Dyadic> w) : LowerBound<Dyadic>(w) {};
};
/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<>
class Positive<UpperBound<Dyadic>> : public UpperBound<Dyadic> {
  public:
    /** @brief Construct a \c Positive from w. */
    Positive(UpperBound<Dyadic> w) : UpperBound<Dyadic>(w) {};
};
/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<>
class Positive<Approximation<Dyadic>> : public Approximation<Dyadic> {
  public:
    /** @brief Construct a \c Positive from w. */
    Positive(Approximation<Dyadic> w) : Approximation<Dyadic>(w) {}
};

/** @brief Construct the typed value from the literal. */
inline Dyadic operator""_dyadic(long double x) { return Dyadic(static_cast<double>(x)); }
/** @brief Construct the typed value from the literal. */
inline Dyadic operator""_dy(long double x) { return operator""_dyadic(x); }
/** @brief Construct the typed value from the literal. */
inline Dyadic operator""_q2(long double x) { return operator""_dyadic(x); }
/** @brief Construct the typed value from the literal. */
inline Dyadic operator""_bin(long double x) { return operator""_dyadic(x); }
/** @brief Construct a dyadic from the supplied literal components. */
inline Dyadic dy_(long double x) { return operator""_dyadic(x); }
/** @brief Construct dyadic from the supplied arguments. */
Dyadic make_dyadic(unsigned long long int);
/** @brief Divide the value by two. */
Dyadic hlf(const Integer& n);
// Exact mixed dyadic/rational enclosure operations. Overloads prevent
// ambiguous implicit promotions between distinct interval specializations.
/** @brief Compute a rational enclosure of the supplied value. */
inline RationalBounds rational_enclosure(const DyadicBounds& x) {
    return RationalBounds(Rational(x.lower_raw()), Rational(x.upper_raw()));
}
/** @brief Return the quotient of the operands. */
inline RationalBounds div(const DyadicBounds& a, const DyadicBounds& b) {
    return div(rational_enclosure(a), rational_enclosure(b));
}
#define ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(NAME)                              \
    inline RationalBounds NAME(RationalBounds const& q, DyadicBounds const& d) { \
        return NAME(q, rational_enclosure(d));                                   \
    }                                                                            \
    inline RationalBounds NAME(DyadicBounds const& d, RationalBounds const& q) { \
        return NAME(rational_enclosure(d), q);                                   \
    }
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(add)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(sub)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(mul)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(div)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(max)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(min)
#undef ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY

#define ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(OP)                                         \
    inline ValidatedKleenean operator OP(RationalBounds const& q, DyadicBounds const& d) { \
        return q OP rational_enclosure(d);                                                 \
    }                                                                                      \
    inline ValidatedKleenean operator OP(DyadicBounds const& d, RationalBounds const& q) { \
        return rational_enclosure(d) OP q;                                                 \
    }
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(==)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(!=)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(<)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(>)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(<=)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(>=)
#undef ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE

} // namespace Ariadne
#endif // ARIADNE_DYADIC_HPP
