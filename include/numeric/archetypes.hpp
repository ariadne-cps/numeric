/***************************************************************************
 *            numeric/archetypes.hpp
 *
 *  Copyright  2020  Pieter Collins
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

/** @file numeric/archetypes.hpp
 *  @brief Declarations and utilities for archetypes in Ariadne Numeric.
 */

#include "foundation/paradigm.hpp"
#include "numeric/concepts.hpp"
#include "utility/typedefs.hpp"

#ifndef ARIADNE_NUMERIC_ARCHETYPES_HPP
#define ARIADNE_NUMERIC_ARCHETYPES_HPP

namespace Ariadne {

class Real;
class TwoExp;

class DoublePrecision;

struct DecimalPlaces;
struct DecimalPrecision;

template<class X>
class Positive;

class ApproximateDouble;
class ExactDouble;
class Integer;
class Dyadic;
class Decimal;
class Rational;
class DoublePrecision;
struct RoundDownward;
struct RoundUpward;
struct RoundToNearest;
struct RoundApproximately;

template<class P>
class Number;
template<class P>
class LowerNumber;
template<class P>
class UpperNumber;
using ValidatedLowerNumber = LowerNumber<ValidatedTag>;
using ValidatedUpperNumber = UpperNumber<ValidatedTag>;

/** @brief Floating-point number parameterized by its precision representation. */
template<class... PRS>
class Float;
using FloatDP = Float<DoublePrecision>;

/** @brief Reference type specifying the requirements of rounded. */
class RoundedArchetype {
  public:
    /** @brief Representation of a floating-point rounding direction. */
    class RoundingModeType {
      public:
        /** @brief Construct a \c RoundingModeType with the specified characteristics. */
        RoundingModeType(RoundDownward);
        /** @brief Construct a \c RoundingModeType with the specified characteristics. */
        RoundingModeType(RoundUpward);
        /** @brief Construct a \c RoundingModeType with the specified characteristics. */
        RoundingModeType(RoundToNearest);
        /** @brief Construct a \c RoundingModeType with the specified characteristics. */
        RoundingModeType(RoundApproximately);
    };
    /** @brief Descriptor of numeric precision. */
    class PrecisionType {
      public:
        /** @brief Test equality of the operands. */
        friend bool operator==(PrecisionType, PrecisionType);
        /** @brief Return the maximum of the operands. */
        friend PrecisionType max(PrecisionType, PrecisionType);
    };

  private:
    using X = RoundedArchetype;
    using RND = RoundingModeType;
    using PR = PrecisionType;

  public:
    /** @brief Return the working precision. */
    PrecisionType precision() const;
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    RoundedArchetype();
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(int);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(int, PrecisionType);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(ApproximateDouble, RoundingModeType, PrecisionType);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(ExactDouble, PrecisionType);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(TwoExp, PrecisionType);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(ExactDouble, RoundingModeType, PrecisionType);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(Integer, RoundingModeType, PrecisionType);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(Dyadic, RoundingModeType, PrecisionType);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(const Decimal&, RoundingModeType, PrecisionType);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(Rational, RoundingModeType, PrecisionType);
    /** @brief Construct a \c RoundedArchetype with the specified characteristics. */
    explicit RoundedArchetype(RoundedArchetype, RoundingModeType, PrecisionType);

    /** @brief Return the rounding mode. */
    static RoundingModeType get_rounding_mode();
    /** @brief Set the rounding mode. */
    static Void set_rounding_mode(RoundingModeType);

    /** @brief Return an infinite value. */
    static X inf(PR);
    /** @brief Return pi at the requested precision. */
    static X pi(RND, PR);

    /** @brief Round the value to an integer. */
    friend X round(X);

    /** @brief Return the sum or unary positive value. */
    friend X operator+(X);
    /** @brief Return the difference or additive inverse. */
    friend X operator-(X);

    /** @brief Return the additive zero with the operand characteristics. */
    friend auto nul(X) -> X;
    /** @brief Return the unary positive value. */
    friend auto pos(X) -> X;
    /** @brief Return the additive inverse. */
    friend auto neg(X) -> X;
    /** @brief Divide the value by two. */
    friend auto hlf(X) -> X;
    /** @brief Return the multiplicative inverse. */
    friend auto rec(RoundingModeType, X) -> X;
    /** @brief Return the sum of the operands. */
    friend auto add(RoundingModeType, X, X) -> X;
    /** @brief Return the difference of the operands. */
    friend auto sub(RoundingModeType, X, X) -> X;
    /** @brief Return the product of the operands. */
    friend auto mul(RoundingModeType, X, X) -> X;
    /** @brief Return the quotient of the operands. */
    friend auto div(RoundingModeType, X, X) -> X;
    /** @brief Compute the product of the first two operands plus the third. */
    friend auto fma(RoundingModeType, X, X, X) -> X;
    /** @brief Raise the value to the specified power. */
    friend auto pow(RoundingModeType, X, Int) -> X;

    /** @brief Return the midpoint value. */
    friend auto med(RoundToNearest, X, X) -> X;
    /*
        friend auto sqrt(RND, X) -> X;
        friend auto exp(RND, X) -> X;
        friend auto log(RND, X) -> X;
        friend auto sin(RND, X) -> X;
        friend auto cos(RND, X) -> X;
        friend auto tan(RND, X) -> X;
    */
    /** @brief Return the inverse tangent. */
    friend auto atan(RND, X) -> X;
    /** @brief Return the absolute value. */
    friend auto abs(X) -> X;
    /** @brief Return the maximum of the operands. */
    friend auto max(X, X) -> X;
    /** @brief Return the minimum of the operands. */
    friend auto min(X, X) -> X;

    /** @brief Test whether the first operand is less than or equal to the second. */
    friend bool operator<=(X, X);
    /** @brief Test whether the first operand is greater than the second. */
    friend bool operator>(X, X);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend bool operator>=(X, X);

    /** @brief Test equality of the operands. */
    friend bool operator==(X, Int);
    /** @brief Test whether the first operand is less than the second. */
    friend bool operator<(X, Int);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend bool operator<=(X, Int);
    /** @brief Test whether the first operand is greater than the second. */
    friend bool operator>(X, Int);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend bool operator>=(X, Int);

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream&, const X&);
    /** @brief Convert to a built-in double. */
    double get_d() const;
};
/** @brief Write the value to the output stream. */
OutputStream& write(OutputStream&, RoundedArchetype, DecimalPrecision, RoundedArchetype::RoundingModeType);
/** @brief Write the value to the output stream. */
OutputStream& write(OutputStream&, RoundedArchetype, DecimalPlaces, RoundedArchetype::RoundingModeType);

/** @brief Return the lower. */
auto get_lower(const Real& r, RoundedArchetype::PrecisionType pr) -> RoundedArchetype;
/** @brief Return the lower. */
auto get_lower(const ValidatedLowerNumber& y, RoundedArchetype::PrecisionType pr) -> RoundedArchetype;
/** @brief Return the upper. */
auto get_upper(const Real& r, RoundedArchetype::PrecisionType pr) -> RoundedArchetype;
/** @brief Return the upper. */
auto get_upper(const ValidatedUpperNumber& y, RoundedArchetype::PrecisionType pr) -> RoundedArchetype;
/** @brief Return the approx. */
auto get_approx(const Real& r, RoundedArchetype::PrecisionType pr) -> RoundedArchetype;
/** @brief Return the approx. */
auto get_approx(const ValidatedLowerNumber& y, RoundedArchetype::PrecisionType pr) -> RoundedArchetype;

/** @brief Convert to the raw float representation. */
RoundedArchetype cast_raw_float(RoundedArchetype::PrecisionType);

/** @brief Associated type definitions for raw float. */
template<class PR>
struct RawFloatTypedef;
/** @brief Associated type definitions for raw float. */
template<>
struct RawFloatTypedef<RoundedArchetype::PrecisionType> {
    typedef RoundedArchetype Type;
};

static_assert(IsRoundedField<RoundedArchetype>);
// template<class X> concept Rounded = RoundedMode<X,typename X::RoundingModeType>;

} // namespace Ariadne

#endif
