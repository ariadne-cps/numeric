/***************************************************************************
 *            numeric/float64-crtp.hpp
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

/** @file numeric/float64-crtp.hpp
 *  @brief Declarations and utilities for float64 crtp in Ariadne Numeric.
 */

#ifndef ARIADNE_FLOAT64_HPP
#define ARIADNE_FLOAT64_HPP

#include "foundation/paradigm.hpp"
#include "numeric/sign.hpp"
#include "utility/metaprogramming.hpp"
#include "utility/module.hpp"

namespace Ariadne {

/************ FloatDP ********************************************************/

struct ExactTag {};

class Boolean;
class ValidatedKleenean;
class ValidatedSierpinskian;
class Fuzzy;

class FloatDPApproximation;
class FloatDPLowerBound;
class FloatDPUpperBound;
class ValidFloatDP;
class FloatDPBounds;
class FloatDPError;
class FloatDP;

class Rational;

FloatDP operator""_x(long double);
/** @brief Construct the typed value from the literal. */
FloatDPError operator""_e(long double);
/** @brief Construct the typed value from the literal. */
FloatDPLowerBound operator""_l(long double);
/** @brief Construct the typed value from the literal. */
FloatDPUpperBound operator""_u(long double);
FloatDPApproximation operator""_a(long double);

Void FloatDP::set_rounding_to_nearest();
Void FloatDP::set_rounding_downward();
Void FloatDP::set_rounding_upward();
Void FloatDP::set_rounding_toward_zero();

/** @brief Add the operand in place. */
template<class N>
inline N& operator+=(NumberObject<N>& n1, const NumberObject<N>& n2) {
    n1.upcast() = n1.upcast() + n2.upcast();
    return n1.upcast();
}
/** @brief Multiply by the operand in place. */
template<class N>
inline N& operator*=(NumberObject<N>& n1, const NumberObject<N>& n2) {
    return n1.upcast() = n1.upcast() * n2.upcast();
    return n1.upcast();
}

/** @brief Double-precision floating-point value with explicit rounding operations. */
class FloatDP : public ScalarObject<FloatDP> {
    volatile double d;

  public:
    typedef Void Paradigm;
    /** @brief Construct a \c FloatDP with the specified characteristics. */
    FloatDP() : d(0.0) {}
    /** @brief Construct a \c FloatDP from dbl. */
    FloatDP(double dbl) : d(dbl) {}
    /** @brief Convert to double. */
    explicit operator double() const { return d; }
    /** @brief Return the underlying floating-point value. */
    FloatDP get_flt() const { return *this; }
    /** @brief Convert to a built-in double. */
    double get_d() const { return d; }
    /** @brief Return the sum or unary positive value. */
    friend FloatDP operator+(FloatDP x);
    /** @brief Return the difference or additive inverse. */
    friend FloatDP operator-(FloatDP x);
    /** @brief Return the sum or unary positive value. */
    friend FloatDP operator+(FloatDP x1, FloatDP x2);
    /** @brief Return the difference or additive inverse. */
    friend FloatDP operator-(FloatDP x1, FloatDP x2);
    /** @brief Return the product. */
    friend FloatDP operator*(FloatDP x1, FloatDP x2);
    /** @brief Return the quotient. */
    friend FloatDP operator/(FloatDP x1, FloatDP x2);
    /** @brief Return the quotient. */
    friend FloatDP operator/(FloatDP x1, Int n2);
    /** @brief Add the operand in place. */
    friend FloatDP& operator+=(FloatDP& x1, FloatDP x2);
    /** @brief Subtract the operand in place. */
    friend FloatDP& operator-=(FloatDP& x1, FloatDP x2);
    /** @brief Multiply by the operand in place. */
    friend FloatDP& operator*=(FloatDP& x1, FloatDP x2);
    /** @brief Divide by the operand in place. */
    friend FloatDP& operator/=(FloatDP& x1, FloatDP x2);
    /** @brief Return the maximum of the operands. */
    friend FloatDP max(FloatDP x1, FloatDP x2);
    /** @brief Return the additive inverse. */
    friend FloatDP neg(FloatDP x);
    /** @brief Return the absolute value. */
    friend FloatDP abs(FloatDP x);
    /** @brief Return the greatest integer not exceeding the value. */
    friend Int32 floor(FloatDP x);
    /** @brief Return the least integer not less than the value. */
    friend Int32 ceil(FloatDP x);
    /** @brief Return the sign of the floating-point value. */
    friend Sign sng(FloatDP x);
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(FloatDP x1, FloatDP x2);
    /** @brief Test equality of the operands. */
    friend Bool operator==(FloatDP x1, FloatDP x2);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend Bool operator<=(FloatDP x1, FloatDP x2);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, FloatDP x);
};
/** @brief Return the sum or unary positive value. */
inline FloatDP operator+(FloatDP x) { return FloatDP{+x.d}; }
/** @brief Return the difference or additive inverse. */
inline FloatDP operator-(FloatDP x) { return FloatDP{-x.d}; }
/** @brief Return the sum or unary positive value. */
inline FloatDP operator+(FloatDP x1, FloatDP x2) { return FloatDP{x1.d + x2.d}; }
/** @brief Return the difference or additive inverse. */
inline FloatDP operator-(FloatDP x1, FloatDP x2) { return FloatDP{x1.d - x2.d}; }
/** @brief Return the product. */
inline FloatDP operator*(FloatDP x1, FloatDP x2) { return FloatDP{x1.d * x2.d}; }
/** @brief Return the quotient. */
inline FloatDP operator/(FloatDP x1, FloatDP x2) { return FloatDP{x1.d / x2.d}; }
/** @brief Return the quotient. */
inline FloatDP operator/(FloatDP x1, Int n2) { return FloatDP{x1.d / n2}; }
/** @brief Add the operand in place. */
inline FloatDP& operator+=(FloatDP& x1, FloatDP x2) {
    x1.d += x2.d;
    return x1;
}
/** @brief Subtract the operand in place. */
inline FloatDP& operator-=(FloatDP& x1, FloatDP x2) {
    x1.d -= x2.d;
    return x1;
}
/** @brief Multiply by the operand in place. */
inline FloatDP& operator*=(FloatDP& x1, FloatDP x2) {
    x1.d *= x2.d;
    return x1;
}
/** @brief Divide by the operand in place. */
inline FloatDP& operator/=(FloatDP& x1, FloatDP x2) {
    x1.d /= x2.d;
    return x1;
}
/** @brief Compare the operands using the comparison result type. */
inline Comparison cmp(FloatDP x1, FloatDP x2) { return (x1.d == x2.d) ? Comparison::EQUAL : (x1.d > x2.d) ? Comparison::GREATER
                                                                                                          : Comparison::LESS; }
/** @brief Test whether the value is not a number. */
Bool isnan(FloatDP x);

/** @brief Return the absolute value. */
FloatDP abs(FloatDP x);
/** @brief Return a lower bound on the absolute value. */
inline FloatDP mig(FloatDP x) { return abs(x); }
/** @brief Return an upper bound on the absolute value. */
inline FloatDP mag(FloatDP x) { return abs(x); }

/** @brief Return the additive inverse. */
inline FloatDP neg(FloatDP x) { return FloatDP{-x.d}; }

/** @brief Test equality of the operands. */
inline Bool operator==(FloatDP x1, FloatDP x2) { return x1.d == x2.d; }
/** @brief Test inequality of the operands. */
inline Bool operator!=(FloatDP x1, FloatDP x2) { return !(x1 == x2); }
/** @brief Test whether the first operand is less than or equal to the second. */
inline Bool operator<=(FloatDP x1, FloatDP x2) { return x1.d <= x2.d; }
/** @brief Test whether the first operand is greater than or equal to the second. */
inline Bool operator>=(FloatDP x1, FloatDP x2) { return (x2 <= x1); }
/** @brief Test whether the first operand is less than the second. */
inline Bool operator<(FloatDP x1, FloatDP x2) { return !(x2 <= x1); }
/** @brief Test whether the first operand is greater than the second. */
inline Bool operator>(FloatDP x1, FloatDP x2) { return !(x1 <= x2); }

/************ FloatDP ***************************************************/

extern const FloatDP inf;

class FloatDP
    : public NumberObject<FloatDP> {
    friend class FloatDPBounds;
    friend class FloatDPApproximation;
    volatile double v;
    static Int output_places;

  public:
    /** @brief Set the output places. */
    static Void set_output_places(Int);

  public:
    typedef ExactTag Paradigm;
    /** @brief Construct a \c FloatDP with the specified characteristics. */
    FloatDP();
    /** @brief Construct a \c FloatDP with the specified characteristics. */
    template<BuiltinIntegral N>
    FloatDP(N);
    /** @brief Construct a \c FloatDP with the specified characteristics. */
    FloatDP(Nat);
    /** @brief Construct a \c FloatDP with the specified characteristics. */
    FloatDP(Int);
    /** @brief Construct a \c FloatDP with the specified characteristics. */
    explicit FloatDP(double);
    /** @brief Construct a \c FloatDP with the specified characteristics. */
    explicit FloatDP(FloatDP);
    /** @brief Convert to float dp bounds. */
    operator FloatDPBounds() const;
    /** @brief Convert to valid float dp. */
    operator ValidFloatDP() const;
    /** @brief Convert to float dp lower bound. */
    operator FloatDPLowerBound() const;
    /** @brief Convert to float dp upper bound. */
    operator FloatDPUpperBound() const;
    /** @brief Convert to float dp approximation. */
    operator FloatDPApproximation() const;
    /** @brief Return the underlying floating-point value. */
    FloatDP get_flt() const;
    /** @brief Convert to a built-in double. */
    double get_d() const;
    /** @brief Return the sum or unary positive value. */
    friend FloatDP operator+(FloatDP);
    /** @brief Return the difference or additive inverse. */
    friend FloatDP operator-(FloatDP);
    /** @brief Return the sum or unary positive value. */
    friend FloatDPBounds operator+(FloatDP x1, FloatDP x2);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPBounds operator-(FloatDP x1, FloatDP x2);
    /** @brief Return the product. */
    friend FloatDPBounds operator*(FloatDP x1, FloatDP x2);
    /** @brief Return the quotient. */
    friend FloatDPBounds operator/(FloatDP x1, FloatDP x2);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FloatDP&);

  private:
    /** @brief Construct a \c FloatDP from n. */
    FloatDP(long long int n, std::nullptr_t);
};
Boolean operator==(FloatDP x1, FloatDP x2);
Boolean operator!=(FloatDP x1, FloatDP x2);
Boolean operator<=(FloatDP x1, FloatDP x2);
Boolean operator>=(FloatDP x1, FloatDP x2);
Boolean operator<(FloatDP x1, FloatDP x2);
Boolean operator>(FloatDP x1, FloatDP x2);
/** @brief Return the minimum of the operands. */
FloatDP min(FloatDP x1, FloatDP x2);
/** @brief Return the maximum of the operands. */
FloatDP max(FloatDP x1, FloatDP x2);

/** @brief Exception reporting float dp. */
class FloatDPError : public NumberObject<FloatDPUpperBound> {
    volatile double e;
    static Int output_places;

  public:
    /** @brief Set the output places. */
    static Void set_output_places(Int);

  public:
    /** @brief Construct a \c FloatDPError with the specified characteristics. */
    FloatDPError();
    template<BuiltinUnsignedIntegral M>
    = dummy > FloatDPError(M m);
    /** @brief Construct a \c FloatDPError from x. */
    template<BuiltinFloatingPoint X>
    explicit FloatDPError(X x);
    /** @brief Construct a \c FloatDPError with the specified characteristics. */
    explicit FloatDPError(double);
    /** @brief Construct a \c FloatDPError with the specified characteristics. */
    explicit FloatDPError(FloatDP);
    /** @brief Convert to float dp upper bound. */
    explicit operator FloatDPUpperBound() const;
    /** @brief Return the underlying floating-point value. */
    FloatDP get_flt() const;
    /** @brief Return the underlying floating-point value. */
    FloatDP& get_flt();
    /** @brief Convert to a built-in double. */
    double get_d() const;
    /** @brief Convert to double. */
    explicit operator double() const;
    /** @brief Return the sum or unary positive value. */
    friend FloatDPError operator+(FloatDPError);
    /** @brief Return the sum or unary positive value. */
    friend FloatDPError operator+(FloatDPError, FloatDPError);
    /** @brief Return the product. */
    friend FloatDPError operator*(FloatDPError, FloatDPError);
    /** @brief Raise the value to the specified power. */
    friend FloatDPError pow(FloatDPError, Nat);
    /** @brief Return the quotient. */
    friend FloatDPUpperBound operator/(FloatDPError x1, FloatDP x2);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FloatDPError&);
    /** @brief Test equality of the operands. */
    friend Fuzzy operator==(FloatDPError, FloatDP);
};

/** @brief Return the maximum of the operands. */
FloatDPError max(FloatDPError x1, FloatDPError x2);

/** @brief Return the quotient. */
template<BuiltinUnsignedIntegral M>
inline FloatDPError operator/(FloatDPError x1, M m2) {
    return FloatDPError(x1.get_d() / m2);
}

/** @brief Validated double-precision value representation. */
class ValidFloatDP : public NumberObject<ValidFloatDP> {
    volatile double v;
    volatile double e;

  private:
    /** @brief Construct a \c ValidFloatDP from d. */
    ValidFloatDP(double d, ExactTag);

  public:
    typedef ValidatedTag Paradigm;
    /** @brief Construct a \c ValidFloatDP with the specified characteristics. */
    ValidFloatDP();
    /** @brief Construct a \c ValidFloatDP from n. */
    template<BuiltinIntegral N>
    ValidFloatDP(N n);
    /** @brief Construct a \c ValidFloatDP with the specified characteristics. */
    explicit ValidFloatDP(double);
    /** @brief Construct a \c ValidFloatDP with the specified characteristics. */
    explicit ValidFloatDP(double, double);
    /** @brief Construct a \c ValidFloatDP with the specified characteristics. */
    ValidFloatDP(FloatDP, FloatDPError);
    /** @brief Convert to float dp approximation. */
    operator FloatDPApproximation() const;
    /** @brief Return the stored value. */
    FloatDP value() const;
    /** @brief Return the error bound. */
    FloatDPError error() const;
    /** @brief Return the upper endpoint. */
    FloatDPUpperBound upper() const;
    /** @brief Return the lower endpoint. */
    FloatDPLowerBound lower() const;
    /** @brief Return the sum or unary positive value. */
    friend ValidFloatDP operator+(ValidFloatDP);
    /** @brief Return the difference or additive inverse. */
    friend ValidFloatDP operator-(ValidFloatDP);
    /** @brief Return the sum or unary positive value. */
    friend ValidFloatDP operator+(ValidFloatDP, ValidFloatDP);
    /** @brief Return the difference or additive inverse. */
    friend ValidFloatDP operator-(ValidFloatDP, ValidFloatDP);
    /** @brief Return the product. */
    friend ValidFloatDP operator*(ValidFloatDP, ValidFloatDP);
    /** @brief Return the quotient. */
    friend ValidFloatDP operator/(ValidFloatDP, ValidFloatDP);
    /** @brief Return the square. */
    friend ValidFloatDP sqr(ValidFloatDP);
    /** @brief Return the multiplicative inverse. */
    friend ValidFloatDP rec(ValidFloatDP x);
    /** @brief Test equality of the operands. */
    friend Bool operator==(ValidFloatDP, Int);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ValidFloatDP&);
};
/** @brief Return the minimum of the operands. */
ValidFloatDP min(ValidFloatDP x1, ValidFloatDP x2);
/** @brief Return the maximum of the operands. */
ValidFloatDP max(ValidFloatDP x1, ValidFloatDP x2);

/** @brief Test whether the operands have the same representation. */
Bool same(ValidFloatDP x1, ValidFloatDP x2);

template<BuiltinIntegral N>
inline ValidFloatDP::ValidFloatDP(N n) : ValidFloatDP(double(n), ExactTag()) {}

/** @brief Validated enclosure with double-precision endpoints. */
class FloatDPBounds : public NumberObject<FloatDPBounds> {
    friend class FloatDPApproximation;
    friend class ValidFloatDP;
    volatile double l;
    volatile double u;
    static Int output_places;

  public:
    /** @brief Set the output places. */
    static Void set_output_places(Int);

  public:
    typedef ValidatedTag Paradigm;
    /** @brief Construct a \c FloatDPBounds with the specified characteristics. */
    FloatDPBounds();
    /** @brief Construct a \c FloatDPBounds with the specified characteristics. */
    FloatDPBounds(Int);
    /** @brief Construct a \c FloatDPBounds with the specified characteristics. */
    FloatDPBounds(const Rational&);
    /** @brief Construct a \c FloatDPBounds with the specified characteristics. */
    explicit FloatDPBounds(double);
    /** @brief Construct a \c FloatDPBounds with the specified characteristics. */
    explicit FloatDPBounds(double, double);
    /** @brief Construct a \c FloatDPBounds with the specified characteristics. */
    FloatDPBounds(FloatDPLowerBound, FloatDPUpperBound);
    /** @brief Convert to float dp approximation. */
    operator FloatDPApproximation() const;
    /** @brief Convert to float dp upper bound. */
    operator FloatDPUpperBound() const;
    /** @brief Convert to float dp lower bound. */
    operator FloatDPLowerBound() const;
    /** @brief Convert to valid float dp. */
    operator ValidFloatDP() const;
    /** @brief Return the sum or unary positive value. */
    friend FloatDPBounds operator+(FloatDPBounds);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPBounds operator-(FloatDPBounds);
    /** @brief Return the sum or unary positive value. */
    friend FloatDPBounds operator+(FloatDPBounds, FloatDPBounds);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPBounds operator-(FloatDPBounds, FloatDPBounds);
    /** @brief Return the product. */
    friend FloatDPBounds operator*(FloatDPBounds, FloatDPBounds);
    /** @brief Return the quotient. */
    friend FloatDPBounds operator/(FloatDPBounds, FloatDPBounds);
    /** @brief Return the square. */
    friend FloatDPBounds sqr(FloatDPBounds);
    /** @brief Return the maximum of the operands. */
    friend FloatDPBounds max(FloatDPBounds, FloatDPBounds);
    /** @brief Return the minimum of the operands. */
    friend FloatDPBounds min(FloatDPBounds, FloatDPBounds);
    /** @brief Test equality of the operands. */
    friend Bool operator==(FloatDPBounds, Int);
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(FloatDPBounds, FloatDPBounds);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FloatDPBounds&);

    /** @brief Return the lower endpoint. */
    FloatDPLowerBound lower() const;
    /** @brief Return the upper endpoint. */
    FloatDPUpperBound upper() const;
    /** @brief Return the error bound. */
    FloatDPError error() const;
    /** @brief Return the difference between the endpoints. */
    FloatDPError width() const;
    /** @brief Return the enclosure radius. */
    FloatDPError radius() const;
};
/** @brief Return the minimum of the operands. */
FloatDPBounds min(FloatDPBounds x1, FloatDPBounds x2);
/** @brief Return the maximum of the operands. */
FloatDPBounds max(FloatDPBounds x1, FloatDPBounds x2);
/** @brief Return an upper bound on the absolute value. */
FloatDPError mag(FloatDPBounds x);
/** @brief Return the square. */
FloatDPBounds sqr(FloatDPBounds);

/** @brief Test whether the operands have the same representation. */
Bool same(FloatDPBounds x1, FloatDPBounds x2);

/** @brief Test equality of the operands. */
ValidatedKleenean operator==(FloatDPBounds x1, FloatDPBounds x2);
/** @brief Test inequality of the operands. */
ValidatedKleenean operator!=(FloatDPBounds x1, FloatDPBounds x2);
/** @brief Test whether the first operand is less than or equal to the second. */
ValidatedKleenean operator<=(FloatDPBounds x1, FloatDPBounds x2);
/** @brief Test whether the first operand is greater than or equal to the second. */
ValidatedKleenean operator>=(FloatDPBounds x1, FloatDPBounds x2);
/** @brief Test whether the first operand is less than the second. */
ValidatedKleenean operator<(FloatDPBounds x1, FloatDPBounds x2);
/** @brief Test whether the first operand is greater than the second. */
ValidatedKleenean operator>(FloatDPBounds x1, FloatDPBounds x2);

/** @brief Validated double-precision lower bound. */
class FloatDPLowerBound : public NumberObject<FloatDPLowerBound> {
    volatile double l;

  public:
    /** @brief Construct a \c FloatDPLowerBound with the specified characteristics. */
    FloatDPLowerBound();
    /** @brief Construct a \c FloatDPLowerBound with the specified characteristics. */
    FloatDPLowerBound(Int);
    /** @brief Construct a \c FloatDPLowerBound with the specified characteristics. */
    explicit FloatDPLowerBound(double);
    /** @brief Convert to float dp approximation. */
    operator FloatDPApproximation() const;
    /** @brief Return the underlying floating-point value. */
    FloatDP get_flt() const;
    /** @brief Return the sum or unary positive value. */
    friend FloatDPLowerBound operator+(FloatDPLowerBound);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPUpperBound operator-(FloatDPLowerBound);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPLowerBound operator-(FloatDPUpperBound);
    /** @brief Return the sum or unary positive value. */
    friend FloatDPLowerBound operator+(FloatDPLowerBound, FloatDPLowerBound);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPLowerBound operator-(FloatDPLowerBound, FloatDPUpperBound);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPUpperBound operator-(FloatDPUpperBound, FloatDPLowerBound);
    /** @brief Return the minimum of the operands. */
    friend FloatDPUpperBound min(FloatDPUpperBound, FloatDPUpperBound);
    /** @brief Return the maximum of the operands. */
    friend FloatDPUpperBound max(FloatDPUpperBound, FloatDPUpperBound);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FloatDPLowerBound&);
};
/** @brief Return the minimum of the operands. */
FloatDPLowerBound min(FloatDPLowerBound, FloatDPLowerBound);
/** @brief Return the maximum of the operands. */
FloatDPLowerBound max(FloatDPLowerBound, FloatDPLowerBound);

/** @brief Return the multiplicative inverse. */
FloatDPUpperBound rec(FloatDPLowerBound);

/** @brief Validated double-precision upper bound. */
class FloatDPUpperBound : public NumberObject<FloatDPUpperBound> {
    volatile double u;

  public:
    /** @brief Construct a \c FloatDPUpperBound with the specified characteristics. */
    FloatDPUpperBound();
    /** @brief Construct a \c FloatDPUpperBound with the specified characteristics. */
    FloatDPUpperBound(Int);
    /** @brief Construct a \c FloatDPUpperBound with the specified characteristics. */
    explicit FloatDPUpperBound(double);
    /** @brief Convert to float dp approximation. */
    operator FloatDPApproximation() const;
    /** @brief Return the underlying floating-point value. */
    FloatDP get_flt() const;
    /** @brief Return the sum or unary positive value. */
    friend FloatDPUpperBound operator+(FloatDPUpperBound);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPUpperBound operator-(FloatDPLowerBound);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPLowerBound operator-(FloatDPUpperBound);
    /** @brief Return the sum or unary positive value. */
    friend FloatDPUpperBound operator+(FloatDPUpperBound, FloatDPUpperBound);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPUpperBound operator-(FloatDPUpperBound, FloatDPLowerBound);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPLowerBound operator-(FloatDPLowerBound, FloatDPUpperBound);
    /** @brief Return the product. */
    friend FloatDPUpperBound operator*(FloatDPUpperBound, FloatDPUpperBound);
    /** @brief Return the minimum of the operands. */
    friend FloatDPLowerBound min(FloatDPLowerBound, FloatDPLowerBound);
    /** @brief Return the maximum of the operands. */
    friend FloatDPLowerBound max(FloatDPLowerBound, FloatDPLowerBound);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FloatDPUpperBound&);
};
/** @brief Return the minimum of the operands. */
FloatDPUpperBound min(FloatDPUpperBound, FloatDPUpperBound);
/** @brief Return the maximum of the operands. */
FloatDPUpperBound max(FloatDPUpperBound, FloatDPUpperBound);
/** @brief Return the multiplicative inverse. */
FloatDPLowerBound rec(FloatDPUpperBound);

/** @brief Return the sum or unary positive value. */
FloatDPUpperBound operator+(FloatDPUpperBound, FloatDPError);
/** @brief Return the difference or additive inverse. */
FloatDPLowerBound operator-(FloatDPLowerBound, FloatDPError);

/** @brief Approximate double-precision numeric value. */
class FloatDPApproximation : public NumberObject<FloatDPApproximation> {
    volatile double a;

  public:
    typedef ApproximateTag Paradigm;
    /** @brief Construct a \c FloatDPApproximation with the specified characteristics. */
    FloatDPApproximation();
    /** @brief Construct a \c FloatDPApproximation with the specified characteristics. */
    FloatDPApproximation(double);
    /** @brief Construct a \c FloatDPApproximation with the specified characteristics. */
    FloatDPApproximation(const Rational&);
    /** @brief Return the underlying floating-point value. */
    FloatDP get_flt() const;
    /** @brief Convert to a built-in double. */
    double get_d() const;
    /** @brief Return the sum or unary positive value. */
    friend FloatDPApproximation operator+(FloatDPApproximation);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPApproximation operator-(FloatDPApproximation);
    /** @brief Return the sum or unary positive value. */
    friend FloatDPApproximation operator+(FloatDPApproximation, FloatDPApproximation);
    /** @brief Return the difference or additive inverse. */
    friend FloatDPApproximation operator-(FloatDPApproximation, FloatDPApproximation);
    /** @brief Return the product. */
    friend FloatDPApproximation operator*(FloatDPApproximation, FloatDPApproximation);
    /** @brief Return the quotient. */
    friend FloatDPApproximation operator/(FloatDPApproximation, FloatDPApproximation);
    /** @brief Add the operand in place. */
    friend FloatDPApproximation& operator+=(FloatDPApproximation&, FloatDPApproximation);
    /** @brief Return the additive inverse. */
    friend FloatDPApproximation neg(FloatDPApproximation);
    /** @brief Return the square. */
    friend FloatDPApproximation sqr(FloatDPApproximation);
    /** @brief Return the minimum of the operands. */
    friend FloatDPApproximation min(FloatDPApproximation, FloatDPApproximation);
    /** @brief Return the maximum of the operands. */
    friend FloatDPApproximation max(FloatDPApproximation, FloatDPApproximation);
    /** @brief Test equality of the operands. */
    friend Fuzzy operator==(FloatDPApproximation, FloatDPApproximation);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FloatDPApproximation&);
};

/** @brief Test equality of the operands. */
Fuzzy operator==(FloatDPApproximation, FloatDPApproximation);
/** @brief Test inequality of the operands. */
Fuzzy operator!=(FloatDPApproximation, FloatDPApproximation);
/** @brief Test whether the first operand is less than or equal to the second. */
Fuzzy operator<=(FloatDPApproximation, FloatDPApproximation);
/** @brief Test whether the first operand is greater than or equal to the second. */
Fuzzy operator>=(FloatDPApproximation, FloatDPApproximation);
/** @brief Test whether the first operand is less than the second. */
Fuzzy operator<(FloatDPApproximation, FloatDPApproximation);
/** @brief Test whether the first operand is greater than the second. */
Fuzzy operator>(FloatDPApproximation, FloatDPApproximation);
/** @brief Return the additive inverse. */
FloatDPApproximation neg(FloatDPApproximation);
/** @brief Return the absolute value. */
FloatDPApproximation abs(FloatDPApproximation);
/** @brief Return the square. */
FloatDPApproximation sqr(FloatDPApproximation);

class TwoExp {
    Int _n;

  public:
    /** @brief Construct a \c TwoExp from n. */
    TwoExp(Int n) : _n(n) {}
    /** @brief Convert to a built-in double. */
    double get_d() const {
        if (_n >= 0) {
            return 1 << _n;
        } else {
            return 1.0 / (1 << (-_n));
        }
    }
    /** @brief Convert to float dp. */
    operator FloatDP() const { return FloatDP(this->get_d()); }
    /** @brief Convert to float dp error. */
    operator FloatDPError() const { return FloatDPError(this->get_d()); }
    /** @brief Convert to valid float dp. */
    operator ValidFloatDP() const { return ValidFloatDP(this->get_d()); }
    /** @brief Convert to float dp bounds. */
    operator FloatDPBounds() const { return FloatDPBounds(this->get_d()); }
};
/** @brief Construct the specified power of two. */
inline TwoExp two_exp(Int n) { return TwoExp(n); }

/** @brief Convert to the exact representation. */
FloatDP cast_exact(const FloatDP&);
/** @brief Construct error from the supplied arguments. */
FloatDPError make_error(const FloatDP&);
/** @brief Construct valid from the supplied arguments. */
ValidFloatDP make_valid(const FloatDP&);
FloatDP cast_raw(const FloatDPApproximation&);

} // namespace Ariadne

#endif
