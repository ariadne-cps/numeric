/***************************************************************************
 *            numeric/flt64.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/** @file numeric/flt64.hpp
 *  @brief Declarations and utilities for flt64 in Ariadne Numeric.
 */

#ifndef ARIADNE_FLT64_HPP
#define ARIADNE_FLT64_HPP

#include <cmath>
#include <iostream>

#include "foundation/logical.decl.hpp"
#include "foundation/logical.hpp"
#include "foundation/paradigm.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/is_number.hpp"
#include "numeric/number.hpp"
#include "numeric/real.hpp"
#include "numeric/sign.hpp"
#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

/************ Flt64 ********************************************************/

/** @brief Set the default rounding. */
Void set_default_rounding();

Void FloatDP::set_rounding_to_nearest();
Void FloatDP::set_rounding_downward();
Void FloatDP::set_rounding_upward();
Void FloatDP::set_rounding_toward_zero();

typedef unsigned short RoundingModeType;
RoundingModeType FloatDP::get_rounding_mode();
Void FloatDP::set_rounding_mode(RoundingModeType);

extern const RoundingModeType ROUND_NEAR;
extern const RoundingModeType ROUND_DOWN;
extern const RoundingModeType ROUND_UP;
extern const RoundingModeType ROUND_ZERO;

// ----------------- RawTag FloatDP class ---------------------------------------------------------

/** @brief Convert the floating-point value to the requested integer type. */
template<BuiltinIntegral Z>
Z integer_cast(Flt64);

// ----------------- RawTag FloatDP class ---------------------------------------------------------

//! \ingroup FloatDPSubModule
//! \brief Wrapper class for double-precision floating-point numbers.
//! \details
//! The main purpose of this class is as a base implementation for the \em safe floating-point number types
//! FloatDPApproximation, FloatDPLowerBound, FloatDPUpperBound,FloatDPBounds, FloatDPBall, FloatDPError and FloatDP.
//! These classes
//!
//! Default arithmetic operations are approximate, and comparisons are exact, so this class is \em unsafe.
//! %Rounded operations on \b FloatDP classes are provided using the \c _down
//! and \c _up suffixes. %ApproximateTag arithmetic without control of the
//! rounding mode can be specified by the \c _approx rounding suffix, or \c _near if round-to-nearest is available.
//! %ExactTag arithmetic can be specified explicitly using the \c _exact rounding suffix.
//!
//! Since the raw \b FloatDP classes represent unsafe values, they cannot be converted to other numbers.
//! The exception is conversion to \b FloatDPApproximation classes, since these represent approximate values with no control of the error.
//! Since the raw \b FloatDP classes are used as the implementation of the safe \c FloatDP classes,
//! they are valid arguments to constructors, but these constructors are all declared \c explicit.
//! Comparison of raw \b %FloatDP data and other numbers is performed as if the %FloatDP object were an \b %FloatDPApproximation.
//!
//! When testing, it is often useful to perform comparisons with \c %Flt64 or \c double values.
//! Although care must be taken, since compiler rounding may change the truth of certain comparisons.
//! comparisons with \c double are performed as if the value were \c FloatDP.

class Flt64 {
  public:
    volatile double d;

  public:
    typedef RawTag Paradigm;
    typedef Flt64 NumericType;
    /** @brief Construct a \c Flt64 with the specified characteristics. */
    Flt64() : d(0.0) {}
    /** @brief Construct a \c Flt64 from dbl. */
    Flt64(double dbl) : d(dbl) {}
    /** @brief Convert to double. */
    explicit operator double() const { return d; }
    /** @brief Return the underlying numeric representation. */
    Flt64 raw() const { return *this; }
    /** @brief Return the underlying floating-point value. */
    Flt64 get_flt() const { return *this; }
    /** @brief Convert to a built-in double. */
    double get_d() const { return d; }

    /** @brief Return the precision-dependent unit roundoff or tolerance. */
    static Flt64 eps();
    /** @brief Return an infinite value. */
    static Flt64 inf();

    /** @brief Set the rounding mode. */
    static Void FloatDP::set_rounding_mode(RoundingModeType);
    /** @brief Return the rounding mode. */
    static RoundingModeType FloatDP::get_rounding_mode();

    /** @brief Set the rounding to nearest. */
    static Void FloatDP::set_rounding_to_nearest();
    /** @brief Set the rounding downward. */
    static Void FloatDP::set_rounding_downward();
    /** @brief Set the rounding upward. */
    static Void FloatDP::set_rounding_upward();
    /** @brief Set the rounding toward zero. */
    static Void FloatDP::set_rounding_toward_zero();

    /** @brief Return the sum or unary positive value. */
    friend inline Flt64 operator+(Flt64 x) { return Flt64{+x.d}; }
    /** @brief Return the difference or additive inverse. */
    friend inline Flt64 operator-(Flt64 x) { return Flt64{-x.d}; }
    /** @brief Return the sum or unary positive value. */
    friend inline Flt64 operator+(Flt64 x1, Flt64 x2) { return Flt64{x1.d + x2.d}; }
    /** @brief Return the difference or additive inverse. */
    friend inline Flt64 operator-(Flt64 x1, Flt64 x2) { return Flt64{x1.d - x2.d}; }
    /** @brief Return the product. */
    friend inline Flt64 operator*(Flt64 x1, Flt64 x2) { return Flt64{x1.d * x2.d}; }
    /** @brief Return the quotient. */
    friend inline Flt64 operator/(Flt64 x1, Flt64 x2) { return Flt64{x1.d / x2.d}; }
    /** @brief Return the quotient. */
    friend inline Flt64 operator/(Flt64 x1, Int n2) { return Flt64{x1.d / n2}; }
    /** @brief Add the operand in place. */
    friend inline Flt64& operator+=(Flt64& x1, Flt64 x2) {
        x1.d += x2.d;
        return x1;
    }
    /** @brief Subtract the operand in place. */
    friend inline Flt64& operator-=(Flt64& x1, Flt64 x2) {
        x1.d -= x2.d;
        return x1;
    }
    /** @brief Multiply by the operand in place. */
    friend inline Flt64& operator*=(Flt64& x1, Flt64 x2) {
        x1.d *= x2.d;
        return x1;
    }
    /** @brief Divide by the operand in place. */
    friend inline Flt64& operator/=(Flt64& x1, Flt64 x2) {
        x1.d /= x2.d;
        return x1;
    }

    /** @brief Return the minimum of the operands. */
    friend inline Flt64 min(Flt64 x1, Flt64 x2) { return Flt64{x2.d < x1.d ? x2.d : x1.d}; }
    /** @brief Return the maximum of the operands. */
    friend inline Flt64 max(Flt64 x1, Flt64 x2) { return Flt64{x2.d > x1.d ? x2.d : x1.d}; }
    /** @brief Return the absolute value. */
    friend inline Flt64 abs(Flt64 x) { return Flt64{std::fabs(x.d)}; }
    /** @brief Return a lower bound on the absolute value. */
    friend inline Flt64 mig(Flt64 x) { return abs(x); }
    /** @brief Return an upper bound on the absolute value. */
    friend inline Flt64 mag(Flt64 x) { return abs(x); }

    /** @brief Return the sum of the operands. */
    friend inline Flt64 add(Flt64 x1, Flt64 x2) { return Flt64{x1.d + x2.d}; }
    /** @brief Return the difference of the operands. */
    friend inline Flt64 sub(Flt64 x1, Flt64 x2) { return Flt64{x1.d - x2.d}; }
    /** @brief Return the product of the operands. */
    friend inline Flt64 mul(Flt64 x1, Flt64 x2) { return Flt64{x1.d * x2.d}; }
    /** @brief Return the quotient of the operands. */
    friend inline Flt64 div(Flt64 x1, Flt64 x2) { return Flt64{x1.d / x2.d}; }

    /** @brief Return the unary positive value. */
    friend inline Flt64 pos(Flt64 x) { return Flt64{x.d}; }
    /** @brief Divide the value by two. */
    friend inline Flt64 hlf(Flt64 x) { return Flt64{x.d / 2}; }
    /** @brief Return the square. */
    friend inline Flt64 sqr(Flt64 x) { return Flt64{x.d * x.d}; }

    /** @brief Return the additive inverse. */
    friend inline Flt64 neg(Flt64 x) { return Flt64{-x.d}; }
    /** @brief Return the multiplicative inverse. */
    friend inline Flt64 rec(Flt64 x) { return Flt64{1.0 / x.d}; }

    /** @brief Raise the value to the specified power. */
    friend inline Flt64 pow(Flt64 x, Nat m) { return std::pow(x.d, m); }
    /** @brief Raise the value to the specified power. */
    friend inline Flt64 pow(Flt64 x, Int n) { return std::pow(x.d, n); }
    /** @brief Return the square root. */
    friend inline Flt64 sqrt(Flt64 x) { return std::sqrt(x.d); }
    /** @brief Return the exponential. */
    friend inline Flt64 exp(Flt64 x) { return std::exp(x.d); }
    /** @brief Return the natural logarithm. */
    friend inline Flt64 log(Flt64 x) { return std::log(x.d); }
    /** @brief Return the sine. */
    friend inline Flt64 sin(Flt64 x) { return std::sin(x.d); }
    /** @brief Return the cosine. */
    friend inline Flt64 cos(Flt64 x) { return std::cos(x.d); }
    /** @brief Return the tangent. */
    friend inline Flt64 tan(Flt64 x) { return std::tan(x.d); }
    /** @brief Return the inverse sine. */
    friend inline Flt64 asin(Flt64 x) { return std::asin(x.d); }
    /** @brief Return the inverse cosine. */
    friend inline Flt64 acos(Flt64 x) { return std::acos(x.d); }
    /** @brief Return the inverse tangent. */
    friend inline Flt64 atan(Flt64 x) { return std::atan(x.d); }

    /** @brief Return the greatest integer not exceeding the value. */
    friend inline Int32Type floor(Flt64 x) { return int32_t(std::floor(x.d)); }
    /** @brief Return the least integer not less than the value. */
    friend inline Int32Type ceil(Flt64 x) { return int32_t(std::ceil(x.d)); }

    /** @brief Test whether the value is not a number. */
    friend Bool isnan(Flt64 x);

    /** @brief Write the value to the output stream. */
    friend inline OutputStream& operator<<(OutputStream& os, Flt64 x) { return os << x.get_d(); }
    /** @brief Read the value from the input stream. */
    friend inline InputStream& operator>>(InputStream& is, Flt64& x) { return is >> reinterpret_cast<double&>(x); }

    /** @brief Compare the operands using the comparison result type. */
    friend inline Comparison cmp(Flt64 x1, Flt64 x2) {
        return (x1.d == x2.d) ? Comparison::EQUAL : (x1.d > x2.d) ? Comparison::GREATER
                                                                  : Comparison::LESS;
    }
    /** @brief Test equality of the operands. */
    friend inline Bool operator==(Flt64 x1, Flt64 x2) { return x1.d == x2.d; }
    /** @brief Test inequality of the operands. */
    friend inline Bool operator!=(Flt64 x1, Flt64 x2) { return !(x1 == x2); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend inline Bool operator<=(Flt64 x1, Flt64 x2) { return x1.d <= x2.d; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend inline Bool operator>=(Flt64 x1, Flt64 x2) { return (x2 <= x1); }
    /** @brief Test whether the first operand is less than the second. */
    friend inline Bool operator<(Flt64 x1, Flt64 x2) { return !(x2 <= x1); }
    /** @brief Test whether the first operand is greater than the second. */
    friend inline Bool operator>(Flt64 x1, Flt64 x2) { return !(x1 <= x2); }

    /** @brief Test equality of the operands. */
    friend inline Bool operator==(Flt64 x1, double x2) { return x1.d == x2; }
    /** @brief Test inequality of the operands. */
    friend inline Bool operator!=(Flt64 x1, double x2) { return !(x1 == x2); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend inline Bool operator<=(Flt64 x1, double x2) { return x1.d <= x2; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend inline Bool operator>=(Flt64 x1, double x2) { return (x2 <= x1); }
    /** @brief Test whether the first operand is less than the second. */
    friend inline Bool operator<(Flt64 x1, double x2) { return !(x2 <= x1); }
    /** @brief Test whether the first operand is greater than the second. */
    friend inline Bool operator>(Flt64 x1, double x2) { return !(x1 <= x2); }
};

inline Flt64 Flt64::inf() { return Flt64{1.0 / 0.0}; }

/************ Rounded arithmetic of Flt64 ***************************************************/

/** @brief Divide the value by two. */
Flt64 hlf(exact, Flt64);

/** @brief Return the next term or representable value. */
Flt64 next(down, Flt64);
/** @brief Return the next term or representable value. */
Flt64 next(up, Flt64);

/** @brief Return the sum of the operands. */
Flt64 add(near, Flt64, Flt64);
/** @brief Return the sum of the operands. */
Flt64 add(down, Flt64, Flt64);
/** @brief Return the sum of the operands. */
Flt64 add(up, Flt64, Flt64);

/** @brief Return the difference of the operands. */
Flt64 sub(near, Flt64, Flt64);
/** @brief Return the difference of the operands. */
Flt64 sub(down, Flt64, Flt64);
/** @brief Return the difference of the operands. */
Flt64 sub(up, Flt64, Flt64);

/** @brief Return the product of the operands. */
Flt64 mul(near, Flt64, Flt64);
/** @brief Return the product of the operands. */
Flt64 mul(down, Flt64, Flt64);
/** @brief Return the product of the operands. */
Flt64 mul(up, Flt64, Flt64);

/** @brief Return the quotient of the operands. */
Flt64 div(near, Flt64, Flt64);
/** @brief Return the quotient of the operands. */
Flt64 div(down, Flt64, Flt64);
/** @brief Return the quotient of the operands. */
Flt64 div(up, Flt64, Flt64);

/** @brief Return the radius. */
Flt64 rad(up, Flt64, Flt64);
/** @brief Return the midpoint value. */
Flt64 med(near, Flt64, Flt64);

/** @brief Raise the value to the specified power. */
Flt64 pow(approx, Flt64, Int);
/** @brief Raise the value to the specified power. */
Flt64 pow(down, Flt64, Int);
/** @brief Raise the value to the specified power. */
Flt64 pow(up, Flt64, Int);

/** @brief Return the square root. */
Flt64 sqrt(approx, Flt64);
/** @brief Return the square root. */
Flt64 sqrt(down, Flt64);
/** @brief Return the square root. */
Flt64 sqrt(up, Flt64);
/** @brief Return the exponential. */
Flt64 exp(approx, Flt64);
/** @brief Return the exponential. */
Flt64 exp(down, Flt64);
/** @brief Return the exponential. */
Flt64 exp(up, Flt64);
/** @brief Return the natural logarithm. */
Flt64 log(approx, Flt64);
/** @brief Return the natural logarithm. */
Flt64 log(down, Flt64);
/** @brief Return the natural logarithm. */
Flt64 log(up, Flt64);
/** @brief Return the sine. */
Flt64 sin(approx, Flt64);
/** @brief Return the sine. */
Flt64 sin(down, Flt64);
/** @brief Return the sine. */
Flt64 sin(up, Flt64);
/** @brief Return the cosine. */
Flt64 cos(approx, Flt64);
/** @brief Return the cosine. */
Flt64 cos(down, Flt64);
/** @brief Return the cosine. */
Flt64 cos(up, Flt64);
/** @brief Return the tangent. */
Flt64 tan(approx, Flt64);
/** @brief Return the tangent. */
Flt64 tan(down, Flt64);
/** @brief Return the tangent. */
Flt64 tan(up, Flt64);
/** @brief Return the inverse tangent. */
Flt64 atan(approx, Flt64);
/** @brief Return the inverse tangent. */
Flt64 atan(down, Flt64);
/** @brief Return the inverse tangent. */
Flt64 atan(up, Flt64);

// Arithmetic respecting the rounding mode
/** @brief Return the sum of the operands with the specified rounding mode. */
Flt64 add_rnd(Flt64 x, Flt64 y);
/** @brief Return the difference of the operands with the specified rounding mode. */
Flt64 sub_rnd(Flt64 x, Flt64 y);
/** @brief Return the product of the operands with the specified rounding mode. */
Flt64 mul_rnd(Flt64 x, Flt64 y);
/** @brief Return the quotient of the operands with the specified rounding mode. */
Flt64 div_rnd(Flt64 x, Flt64 y);
/** @brief Return the sum of the operands with opposite directed rounding. */
Flt64 add_opp(Flt64 x, Flt64 y);
/** @brief Return the difference of the operands with opposite directed rounding. */
Flt64 sub_opp(Flt64 x, Flt64 y);
/** @brief Return the product of the operands with opposite directed rounding. */
Flt64 mul_opp(Flt64 x, Flt64 y);
/** @brief Return the quotient of the operands with opposite directed rounding. */
Flt64 div_opp(Flt64 x, Flt64 y);

} // namespace Ariadne

#endif
