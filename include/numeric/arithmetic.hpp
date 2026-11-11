/***************************************************************************
 *            numeric/arithmetic.hpp
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

/** @file numeric/arithmetic.hpp
 *  @brief Arithmetic declarations
 */

#ifndef ARIADNE_ARITHMETIC_HPP
#define ARIADNE_ARITHMETIC_HPP

#include "foundation/logical.decl.hpp"
#include "foundation/logical.hpp" // TODO: Try to remove; needed for specialisation of Boolean DefineMixedComparisonOperators
#include "number.decl.hpp"
#include "sign.hpp"
#include "utility/metaprogramming.hpp"

namespace Ariadne {

class Natural;
class Integer;

/** @brief Compute a power using the generic arithmetic operations. */
template<class X>
X generic_pow(const X& x, Nat m) {
    X r = x;
    r *= 0;
    r += 1;
    X p = x;
    while (m != 0) {
        if (m % 2 == 1) {
            r = r * p;
        }
        p = p * p;
        m = m / 2;
    }
    return r;
}
/** @brief Compute a power using the generic arithmetic operations. */
template<class X>
X generic_pow(const X& x, Int n) {
    return n >= 0 ? generic_pow(x, Nat(n)) : rec(generic_pow(x, Nat(-n)));
}

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a mutable field.
template<class X, class Y, class NY, class R = X>
struct DeclareInplaceMixedDirectedGroupOperators {
};
/** @brief Declare inplace mixed directed group operators for the supported operand types. */
template<class X, class Y, class NY>
struct DeclareInplaceMixedDirectedGroupOperators<X, Y, NY, X> {
    /** @brief Add the operand in place. */
    friend X& operator+=(X& x1, const Y& x2);
    /** @brief Subtract the operand in place. */
    friend X& operator-=(X& x1, NY const& x2);
};
/** @brief Declare inplace directed group operators for the supported operand types. */
template<class X, class NX, class R = X>
struct DeclareInplaceDirectedGroupOperators : DeclareInplaceMixedDirectedGroupOperators<X, X, NX, R> {
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a mutable semifield.
template<class X, class Y, class QY, class R = X>
struct DeclareInplaceMixedDirectedSemifieldOperators {
};
/** @brief Declare inplace mixed directed semifield operators for the supported operand types. */
template<class X, class Y, class QY>
struct DeclareInplaceMixedDirectedSemifieldOperators<X, Y, QY, X> {
    /** @brief Add the operand in place. */
    friend X& operator+=(X& x1, const Y& x2);
    /** @brief Multiply by the operand in place. */
    friend X& operator*=(X& x1, const Y& x2);
    /** @brief Divide by the operand in place. */
    friend X& operator/=(X& x1, QY const& x2);
};
/** @brief Declare inplace directed semifield operators for the supported operand types. */
template<class X, class QX, class R = X>
struct DeclareInplaceDirectedSemifieldOperators : DeclareInplaceMixedDirectedSemifieldOperators<X, X, QX, R> {
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a mutable field.
template<class X, class Y, class R = X>
struct DeclareInplaceMixedFieldOperators {
};
/** @brief Declare inplace mixed field operators for the supported operand types. */
template<class X, class Y>
struct DeclareInplaceMixedFieldOperators<X, Y, X> {
    //! \brief Inplace addition. May be implemented as \c x1=x1+x2.
    friend X& operator+=(X& x1, const Y& x2);
    //! \brief Inplace subtraction. May be implemented as \c x1=x1-x2.
    friend X& operator-=(X& x1, const Y& x2);
    //! \brief Inplace multiplication. May be implemented as \c x1=x1*x2.
    friend X& operator*=(X& x1, const Y& x2);
    //! \brief Inplace division. May be implemented as \c x1=x1/x2.
    friend X& operator/=(X& x1, const Y& x2);
};
/** @brief Declare inplace field operators for the supported operand types. */
template<class X, class R = X>
struct DeclareInplaceFieldOperators : DeclareInplaceMixedFieldOperators<X, X, R> {
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a (mathematical) ring.
//! \brief \sa DeclareOrderedOperations, DeclareComparisons
//! \details Inheriting from this class declares operators which are found by argument-dependent lookup.
template<class X, class NX, class R = X, class NR = NX>
struct DeclareDirectedGroupOperations
    : DeclareInplaceDirectedGroupOperators<X, NX, R> {
    //!@{
    //! \name Named operations.

    //! \brief Zero \a 0. Preserves accuracy parameters.
    friend X nul(const X& x);
    //! \brief Positive \a +x. Usually constructs a direct copy.
    friend X pos(const X& x);
    //! \brief Negative \a -x. Usually constructs a direct copy with reversed sign.
    friend NX neg(const X& x);
    /** @brief Return the additive inverse. */
    friend X neg(NX const& x);

    //! \brief Sum \a x1+x2.
    friend R add(const X& x1, const X& x2);
    //! \brief Difference \a x1-x2. May be implemented in terms of add() and neg() as <code>add(x1,neg(x2))</code>.
    friend R sub(const X& x1, NX const& x2);
    /** @brief Return the difference of the operands. */
    friend NR sub(NX const& x1, const X& x2);
    //!@}

    //!@{
    //! \name  Standard overloadable operators.

    //! \brief Positive \a +x. Usually dispatches to <code>pos(x)</code>
    friend X operator+(const X& x);
    //! \brief Negative \a -x. Usually dispatches to <code>neg(x)</code>
    friend NX operator-(const X& x);
    /** @brief Return the difference or additive inverse. */
    friend X operator-(NX const& x);
    //! \brief Sum \a x1+x2. Usually dispatches to <code>sum(x1,x2)</code>
    friend R operator+(const X& x1, const X& x2);
    //! \brief Difference \a x1-x2. Usually dispatches to <code>sub(x1,x2)</code>
    friend R operator-(const X& x1, NX const& x2);
    /** @brief Return the difference or additive inverse. */
    friend NR operator-(NX const& x1, const X& x2);
    //!@}
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a (mathematical) group.
template<class X, class R = X>
struct DeclareGroupOperations
    : DeclareDirectedGroupOperations<X, X, R> {};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a (mathematical) ring.
//! \brief \sa DeclareOrderedOperations, DeclareComparisons
//! \details Inheriting from this class declares operators which are found by argument-dependent lookup.
template<class X, class R = X, class PR = R>
struct DeclareRingOperations
    : DeclareGroupOperations<X, R> {
    //!@{
    //! \name Named operations.

    //! \brief Square \a x^2.
    friend PR sqr(const X& x);
    //! \brief Product \a x1*x2.
    friend R mul(const X& x1, const X& x2);
    //! \brief Power \a x<sup>m</sup>. May be implemented in terms of mul() and/or sqr().
    friend R pow(const X& x, Nat m);
    //!@}

    //!@{
    //! \name  Standard overloadable operators.

    //! \brief Product \a x1*x2. Usually dispatches to <code>mul(x1,x2)</code>
    friend R operator*(const X& x1, const X& x2);

    //!@}
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being a (mathematical) field.
//! Includes all operations supported in DeclareRingOperations<X>, and those listed below.
template<class X, class R = X>
struct DeclareFieldOperations
    : DeclareRingOperations<X, R>,
      DeclareInplaceFieldOperators<X, R> {
    //! \brief Half \a x/2.
    friend R hlf(const X& x);
    //! \brief Reciprocal \a 1/x.
    friend R rec(const X& x);
    //! \brief Quotient \a x1/x2. May be implemented in terms of mul() and rec() as <code>mul(x1,rec(x2))</code>.
    friend R div(const X& x1, const X& x2);
    //! \brief Power \a x<sup>n</sup>. May be implemented using rec() and pow(X,Nat) as  <code> n>=0 ? pow(x,Nat(n)) : rec(pow(x,Nat(-n)))</code>
    friend R pow(const X& x, Int n);

    //! \brief Quotient \a x1*x2. Usually dispatches to <code>div(x1,x2)</code>
    friend R operator/(const X&, const X&);
};

/** @brief Declare field operators for the supported operand types. */
template<class X, class R = X>
struct DeclareFieldOperators : DeclareInplaceFieldOperators<X, R> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X&, const X&);
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const X&, const X&);
    /** @brief Return the product. */
    friend R operator*(const X&, const X&);
    /** @brief Return the quotient. */
    friend R operator/(const X&, const X&);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare elementary algebraic and transcendental operations.
template<class X, class R = X>
class DeclareMonotoneOperations {
    //! \brief The square root of \a x, √\a x. Requires \c x>=0.
    friend X sqrt(const X& x);
    //! \brief The natural exponent of \a x, \em e<sup>x</sup>.
    friend X exp(const X& x);
    //! \brief The natural logarithm of \a x. Requires \c x>=0.
    friend X log(const X& x);
    //! \brief The arc-tangent of \a x.
    friend X atan(const X& x);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare elementary algebraic and transcendental operations.
template<class X, class R = X>
class DeclareTranscendentalOperations
    : DeclareMonotoneOperations<X, R> {
    //! \brief The sine of \a x.
    friend X sin(const X& x);
    //! \brief The cosine of \a x.
    friend X cos(const X& x);
    //! \brief The tangent of \a x, sin(\a x)/cos(\a x).
    friend X tan(const X& x);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare elementary arithmetic, algebraic and transcendental operations.
template<class X, class R = X>
class DeclareAnalyticFieldOperations
    : DeclareFieldOperations<X, R>,
      DeclareTranscendentalOperations<X, R> {
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare operators corresponding to \a X being an ordered mutable field.
//! \sa DeclareRingOperations, DeclareFieldOperations, DeclareComparisons
template<class X, class PX = Void>
class DeclareLatticeOperations {
    //! \brief The absolute value of \a x.
    friend PX abs(const X& x);
    //! \brief The mimimum of \a x1 and \a x2.
    friend X min(const X& x1, const X& x2);
    /** @brief Return the minimum of the operands. */
    friend PX min(PX const& x1, PX const& x2);
    //! \brief The maximum of \a x1 and \a x2.
    friend X max(const X& x1, const X& x2);
    /** @brief Return the maximum of the operands. */
    friend PX max(PX const& x1, const X& x2);
    /** @brief Return the maximum of the operands. */
    friend PX max(const X& x1, PX const& x2);
    /** @brief Return the maximum of the operands. */
    friend PX max(PX const& x1, PX const& x2);
    // friend U mag(X const& x);
    // friend L mig(X const& x);
};
/** @brief Declare lattice operations for the supported operand types. */
template<class X>
class DeclareLatticeOperations<X, Void> {
    /** @brief Return the minimum of the operands. */
    friend X min(const X& x1, const X& x2);
    /** @brief Return the maximum of the operands. */
    friend X max(const X& x1, const X& x2);
};

//! \ingroup NumericAlgebraSubModule
//! \brief Declare comparison operations.
//! \sa DeclareRingOperations, DeclareComparisons
template<class X, class LT, class EQ = LT>
struct DeclareComparisonOperations {
    typedef LT GT;
    typedef decltype(not declval<LT>()) GEQ;
    typedef decltype(not declval<GT>()) LEQ;
    typedef decltype(not declval<EQ>()) NEQ;

    //! \brief Tests if \a x1 is equal to \a x2.
    friend EQ eq(const X& x1, const X& x2);
    //! \brief Tests if \a x1 is less than \a x2.
    friend LT lt(const X& x1, const X& x2);
    //! \brief Tests if \a x1 is equal to \a x2. If equality is undecidable, may not return \a true even if values are equal.
    friend EQ operator==(const X& x1, const X& x2);
    //! \brief Tests if \a x1 is strictly less than \a x2.
    friend LT operator<(const X& x1, const X& x2);
    //! \brief Tests if \a x1 is not equal to \a x2.
    friend NEQ operator!=(const X& x1, const X& x2);
    //! \brief Tests if \a x1 is strictly greater than \a x2.
    friend GT operator>(const X& x1, const X& x2);
    //! \brief Tests if \a x1 is less than or equal to \a x2.
    friend LEQ operator<=(const X& x1, const X& x2);
    //! \brief Tests if \a x1 is greater than or equal to \a x2.
    friend GEQ operator>=(const X& x1, const X& x2);
};

/** @brief Declare directed comparison operations for the supported operand types. */
template<class X, class NX, class LT, class EQ>
struct DeclareDirectedComparisonOperations {
    typedef decltype(not declval<EQ>()) NEQ;
    typedef decltype(not declval<LT>()) GT;
    typedef LT LEQ;
    typedef GT GEQ;
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, NX const& x2);
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, NX const& x2);
    /** @brief Test equality of the operands. */
    friend EQ operator==(const X& x1, NX const& x2);
    /** @brief Test equality of the operands. */
    friend EQ operator==(NX const& x1, const X& x2);
    /** @brief Test inequality of the operands. */
    friend NEQ operator!=(const X& x1, NX const& x2);
    /** @brief Test inequality of the operands. */
    friend NEQ operator!=(NX const& x1, const X& x2);
    /** @brief Test whether the first operand is less than the second. */
    friend LT operator<(const X& x1, NX const& x2);
    /** @brief Test whether the first operand is greater than the second. */
    friend GT operator>(const X& x1, NX const& x2);
    /** @brief Test whether the first operand is greater than the second. */
    friend LT operator>(NX const& x1, const X& x2);
    /** @brief Test whether the first operand is less than the second. */
    friend GT operator<(NX const& x1, const X& x2);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend LEQ operator<=(const X& x1, NX const& x2);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend GEQ operator>=(const X& x1, NX const& x2);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend LEQ operator>=(NX const& x1, const X& x2);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend GEQ operator<=(NX const& x1, const X& x2);
};

/** @brief Declare ordered analytic field operations for the supported operand types. */
template<class X, class LT, class EQ = LT>
class DeclareOrderedAnalyticFieldOperations
    : DeclareFieldOperations<X>,
      DeclareTranscendentalOperations<X>,
      DeclareLatticeOperations<X>,
      DeclareComparisonOperations<X, LT, EQ> {
};

/** @brief Declare mixed directed group operators for the supported operand types. */
template<class X, class NX, class Y, class NY, class R = X, class NR = NX>
struct DeclareMixedDirectedGroupOperators
    : DeclareInplaceMixedDirectedGroupOperators<X, Y, NY, R>,
      DeclareInplaceMixedDirectedGroupOperators<NX, NY, Y, NR> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const Y& y2);
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const X& x1, NY const& y2);
    /** @brief Return the sum or unary positive value. */
    friend NR operator+(NX const& x1, NY const& y2);
    /** @brief Return the difference or additive inverse. */
    friend NR operator-(NX const& x1, const Y& y2);
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const Y& y1, const X& x2);
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const Y& y1, NX const& x2);
    /** @brief Return the sum or unary positive value. */
    friend NR operator+(NY const& y1, NX const& x2);
    /** @brief Return the difference or additive inverse. */
    friend NR operator-(NY const& y1, const X& x2);
};

/** @brief Declare mixed directed semifield operators for the supported operand types. */
template<class X, class QX, class Y, class QY, class R = X, class QR = QX>
struct DeclareMixedDirectedSemifieldOperators
    : DeclareInplaceMixedDirectedSemifieldOperators<X, Y, QY, R>,
      DeclareInplaceMixedDirectedGroupOperators<QX, QY, Y, QR> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const Y& y2);
    /** @brief Return the product. */
    friend R operator*(const X& x1, const Y& y2);
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, QY const& y2);
    /** @brief Return the sum or unary positive value. */
    friend QR operator+(QX const& x1, QY const& y2);
    /** @brief Return the product. */
    friend QR operator*(QX const& x1, QY const& y2);
    /** @brief Return the quotient. */
    friend QR operator/(QX const& x1, const Y& y2);
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const Y& y1, const X& x2);
    /** @brief Return the product. */
    friend R operator*(const Y& y1, const X& x2);
    /** @brief Return the quotient. */
    friend R operator/(const Y& y1, QX const& x2);
    /** @brief Return the sum or unary positive value. */
    friend QR operator+(QY const& y1, QX const& x2);
    /** @brief Return the product. */
    friend QR operator*(QY const& y1, QX const& x2);
    /** @brief Return the quotient. */
    friend QR operator/(QY const& y1, const X& x2);
};

/** @brief Declare mixed field operators for the supported operand types. */
template<class X, class Y, class R = X>
struct DeclareMixedFieldOperators
    : DeclareInplaceMixedFieldOperators<X, Y, R> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const Y& y2);
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const X& x1, const Y& y2);
    /** @brief Return the product. */
    friend R operator*(const X& x1, const Y& y2);
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, const Y& y2);
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const Y& y1, const X& x2);
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const Y& y1, const X& x2);
    /** @brief Return the product. */
    friend R operator*(const Y& y1, const X& x2);
    /** @brief Return the quotient. */
    friend R operator/(const Y& y1, const X& x2);
};
/** @brief Declare mixed field operations for the supported operand types. */
template<class X, class Y, class R = X>
struct DeclareMixedFieldOperations
    : DeclareMixedFieldOperators<X, Y, R> {
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const Y& y2);
    /** @brief Return the difference of the operands. */
    friend R sub(const X& x1, const Y& y2);
    /** @brief Return the product of the operands. */
    friend R mul(const X& x1, const Y& y2);
    /** @brief Return the quotient of the operands. */
    friend R div(const X& x1, const Y& y2);
    /** @brief Return the sum of the operands. */
    friend R add(const Y& y1, const X& x2);
    /** @brief Return the difference of the operands. */
    friend R sub(const Y& y1, const X& x2);
    /** @brief Return the product of the operands. */
    friend R mul(const Y& y1, const X& x2);
    /** @brief Return the quotient of the operands. */
    friend R div(const Y& y1, const X& x2);
};

/** @brief Declare mixed arithmetic operators for the supported operand types. */
template<class X, class Y, class R = X>
struct DeclareMixedArithmeticOperators : DeclareMixedFieldOperators<X, Y, R> {};
/** @brief Declare mixed arithmetic operations for the supported operand types. */
template<class X, class Y, class R = X>
struct DeclareMixedArithmeticOperations : DeclareMixedFieldOperations<X, Y, R> {};

/** @brief Declare directed numeric operations for the supported operand types. */
template<class X, class NX, class R = X>
struct DeclareDirectedNumericOperations
    : DeclareDirectedGroupOperations<X, NX, R>,
      DeclareMonotoneOperations<X, R>,
      DeclareLatticeOperations<X> {
    /** @brief Return the additive zero with the operand characteristics. */
    friend X nul(const X& x);
    /** @brief Return the unary positive value. */
    friend X pos(const X& x);
    /** @brief Return the additive inverse. */
    friend X neg(NX const& x);
    /** @brief Divide the value by two. */
    friend X hlf(const X& x);
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const X& x2);
    /** @brief Return the difference of the operands. */
    friend R sub(const X& x1, NX const& x2);
    /** @brief Return the square root. */
    friend R sqrt(const X& x);
    /** @brief Return the exponential. */
    friend R exp(const X& x);
    /** @brief Return the natural logarithm. */
    friend R log(const X& x);
    /** @brief Return the inverse tangent. */
    friend R atan(const X& x);
    /** @brief Return the maximum of the operands. */
    friend X max(const X& x1, const X& x2);
    /** @brief Return the minimum of the operands. */
    friend X min(const X& x1, const X& x2);
};

/** @brief Declare numeric operations for the supported operand types. */
template<class X, class R = X>
struct DeclareNumericOperations
    : DeclareAnalyticFieldOperations<X, R> {
    /** @brief Return the additive zero with the operand characteristics. */
    friend X nul(const X& x);
    /** @brief Return the unary positive value. */
    friend X pos(const X& x);
    /** @brief Return the additive inverse. */
    friend X neg(const X& x);
    /** @brief Divide the value by two. */
    friend X hlf(const X& x);
    /** @brief Return the square. */
    friend R sqr(const X& x);
    /** @brief Return the multiplicative inverse. */
    friend R rec(const X& x);
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const X& x2);
    /** @brief Return the difference of the operands. */
    friend R sub(const X& x1, const X& x2);
    /** @brief Return the product of the operands. */
    friend R mul(const X& x1, const X& x2);
    /** @brief Return the quotient of the operands. */
    friend R div(const X& x1, const X& x2);
    /** @brief Raise the value to the specified power. */
    friend R pow(const X& x, Nat m);
    /** @brief Raise the value to the specified power. */
    friend R pow(const X& x, Int n);
    /** @brief Return the square root. */
    friend R sqrt(const X& x);
    /** @brief Return the exponential. */
    friend R exp(const X& x);
    /** @brief Return the natural logarithm. */
    friend R log(const X& x);
    /** @brief Return the sine. */
    friend R sin(const X& x);
    /** @brief Return the cosine. */
    friend R cos(const X& x);
    /** @brief Return the tangent. */
    friend R tan(const X& x);
    /** @brief Return the inverse sine. */
    friend R asin(const X& x);
    /** @brief Return the inverse cosine. */
    friend R acos(const X& x);
    /** @brief Return the inverse tangent. */
    friend R atan(const X& x);
    /** @brief Return the maximum of the operands. */
    friend X max(const X& x1, const X& x2);
    /** @brief Return the minimum of the operands. */
    friend X min(const X& x1, const X& x2);
    /** @brief Return the absolute value. */
    friend X abs(const X& x);
};

/** @brief Declare directed numeric operators for the supported operand types. */
template<class X, class NX = X>
struct DeclareDirectedNumericOperators {
    /** @brief Return the sum or unary positive value. */
    friend X operator+(const X& x);
    /** @brief Return the difference or additive inverse. */
    friend X operator-(NX const& x);
    /** @brief Return the difference or additive inverse. */
    friend NX operator-(const X& x);
    /** @brief Return the sum or unary positive value. */
    friend X operator+(const X& x1, const X& x2);
    /** @brief Return the difference or additive inverse. */
    friend X operator-(const X& x1, NX const& x2);
    /** @brief Return the difference or additive inverse. */
    friend NX operator-(NX const& x1, const X& x2);
    /** @brief Add the operand in place. */
    friend X& operator+=(X& x1, const X& x2);
    /** @brief Subtract the operand in place. */
    friend X& operator-=(X& x1, NX const& x2);

    /** @brief Return the additive zero with the operand characteristics. */
    friend X nul(const X& x);
    /** @brief Return the unary positive value. */
    friend X pos(const X& x);
    /** @brief Return the additive inverse. */
    friend X neg(NX const& x);
    /** @brief Return the additive inverse. */
    friend NX neg(const X& x);
    /** @brief Divide the value by two. */
    friend X hlf(const X& x);
    /** @brief Return the sum of the operands. */
    friend X add(const X& x1, const X& x2);
    /** @brief Return the difference of the operands. */
    friend X sub(const X& x1, NX const& x2);
    /** @brief Return the difference of the operands. */
    friend NX sub(NX const& x1, const X& x2);

    /** @brief Return the maximum of the operands. */
    friend X max(const X& x1, const X& x2);
    /** @brief Return the minimum of the operands. */
    friend X min(const X& x1, const X& x2);
};

/** @brief Declare positive directed numeric operations for the supported operand types. */
template<class X, class QX>
struct DeclarePositiveDirectedNumericOperations {
    /** @brief Return the sum of the operands. */
    friend X add(const X&, const X&);
    /** @brief Return the product of the operands. */
    friend X mul(const X&, const X&);
    /** @brief Raise the value to the specified power. */
    friend X pow(const X&, const Natural&);
    /** @brief Return the quotient of the operands. */
    friend X div(const X&, QX const&);
    /** @brief Return the multiplicative inverse. */
    friend QX rec(const X&);
    /** @brief Return the square root. */
    friend X sqrt(const X&);
    /** @brief Return the inverse tangent. */
    friend X atan(const X&);

    /** @brief Return the maximum of the operands. */
    friend X max(const X&, const X&);
    /** @brief Return the minimum of the operands. */
    friend X min(const X&, const X&);
};

/** @brief Declare real operations for the supported operand types. */
template<class X, class PX = X>
struct DeclareRealOperations
    : DeclareFieldOperators<X> {
    /** @brief Return the additive inverse. */
    friend X neg(const X&);
    /** @brief Divide the value by two. */
    friend X hlf(const X&);
    /** @brief Return the sum of the operands. */
    friend X add(const X&, const X&);
    /** @brief Return the difference of the operands. */
    friend X sub(const X&, const X&);
    /** @brief Return the product of the operands. */
    friend X mul(const X&, const X&);
    /** @brief Raise the value to the specified power. */
    friend X pow(const X&, const Integer&);
    /** @brief Return the multiplicative inverse. */
    friend X rec(const X&);
    /** @brief Return the square root. */
    friend PX sqrt(PX const&);
    /** @brief Return the square root. */
    friend X sqrt(const X&);
    /** @brief Return the exponential. */
    friend X exp(const X&);
    /** @brief Return the natural logarithm. */
    friend X log(const X&);
    /** @brief Return the sine. */
    friend X sin(const X&);
    /** @brief Return the cosine. */
    friend X cos(const X&);
    /** @brief Return the tangent. */
    friend X tan(const X&);
    /** @brief Return the inverse tangent. */
    friend X atan(const X&);
    /** @brief Return the inverse tangent. */
    friend PX atan(PX const&);

    /** @brief Return the maximum of the operands. */
    friend X max(const X&, const X&);
    /** @brief Return the maximum of the operands. */
    friend PX max(const X&, PX const&);
    /** @brief Return the maximum of the operands. */
    friend PX max(PX const&, const X&);
    /** @brief Return the minimum of the operands. */
    friend X min(const X&, const X&);
    /** @brief Return the minimum of the operands. */
    friend PX min(PX const&, PX const&);
    /** @brief Return the absolute value. */
    friend PX abs(const X&);
    /** @brief Return the distance between the operands. */
    friend PX dist(const X&, const X&);
};

/** @brief Declare directed real operations for the supported operand types. */
template<class X, class NX, class PX>
struct DeclareDirectedRealOperations {
    /** @brief Return the additive inverse. */
    friend NX neg(const X&);
    /** @brief Return the sum of the operands. */
    friend X add(const X&, const X&);
    /** @brief Return the difference of the operands. */
    friend X sub(const X&, NX const&);
    /** @brief Return the square root. */
    friend PX sqrt(PX const&);
    /** @brief Return the exponential. */
    friend PX exp(const X&);
    /** @brief Return the natural logarithm. */
    friend X log(PX const&);
    /** @brief Return the inverse tangent. */
    friend X atan(const X&);

    /** @brief Return the maximum of the operands. */
    friend X max(const X&, const X&);
    /** @brief Return the maximum of the operands. */
    friend PX max(const X&, PX const&);
    /** @brief Return the maximum of the operands. */
    friend PX max(PX const&, const X&);
    /** @brief Return the minimum of the operands. */
    friend X min(const X&, const X&);
};

/** @brief Declare positive real operations for the supported operand types. */
template<class PX>
struct DeclarePositiveRealOperations {
    /** @brief Return the sum of the operands. */
    friend PX add(PX const&, PX const&);
    /** @brief Return the product of the operands. */
    friend PX mul(PX const&, PX const&);
    /** @brief Return the quotient of the operands. */
    friend PX div(PX const&, PX const&);
    /** @brief Return the multiplicative inverse. */
    friend PX rec(PX const&);
    /** @brief Return the square root. */
    friend PX sqrt(PX const&);
    /** @brief Return the inverse tangent. */
    friend PX atan(PX const&);

    /** @brief Return the maximum of the operands. */
    friend PX max(PX const&, PX const&);
    /** @brief Return the minimum of the operands. */
    friend PX min(PX const&, PX const&);
};

/** @brief Declare positive directed real operations for the supported operand types. */
template<class X, class QX>
struct DeclarePositiveDirectedRealOperations {
    /** @brief Return the sum of the operands. */
    friend X add(const X&, const X&);
    /** @brief Return the product of the operands. */
    friend X mul(const X&, const X&);
    /** @brief Raise the value to the specified power. */
    friend X pow(const X&, const Natural&);
    /** @brief Return the quotient of the operands. */
    friend X div(const X&, QX const&);
    /** @brief Return the multiplicative inverse. */
    friend QX rec(const X&);
    /** @brief Return the square root. */
    friend X sqrt(const X&);
    /** @brief Return the inverse tangent. */
    friend X atan(const X&);

    /** @brief Return the maximum of the operands. */
    friend X max(const X&, const X&);
    /** @brief Return the minimum of the operands. */
    friend X min(const X&, const X&);
};

/** @brief Declare arithmetic operators for the supported operand types. */
template<class T, class NT = T, class QT = NT>
struct DeclareArithmeticOperators {
    /** @brief Return the sum or unary positive value. */
    friend T operator+(const T& t);
    /** @brief Return the difference or additive inverse. */
    friend T operator-(NT const& t);
    /** @brief Return the difference or additive inverse. */
    friend NT operator-(const T& t);
    /** @brief Return the sum or unary positive value. */
    friend T operator+(const T& t1, const T& t2);
    /** @brief Return the difference or additive inverse. */
    friend T operator-(const T& t1, NT const& t2);
    /** @brief Return the product. */
    friend T operator*(const T& t1, const T& t2);
    /** @brief Return the quotient. */
    friend T operator/(const T& t1, QT const& t2);
    /** @brief Add the operand in place. */
    friend T& operator+=(T& t1, const T& t2);
    /** @brief Subtract the operand in place. */
    friend T& operator-=(T& t1, NT const& t2);
    /** @brief Multiply by the operand in place. */
    friend T& operator*=(T& t1, const T& t2);
    /** @brief Divide by the operand in place. */
    friend T& operator/=(T& t1, QT const& t2);
};
/** @brief Implement arithmetic operators for the supported operand types. */
template<class T, class NT = T, class QT = NT>
struct DefineArithmeticOperators {
    /** @brief Return the sum or unary positive value. */
    friend T operator+(const T& t) { return pos(t); }
    /** @brief Return the difference or additive inverse. */
    friend NT operator-(const T& t) { return neg(t); }
    /** @brief Return the sum or unary positive value. */
    friend T operator+(const T& t1, const T& t2) { return add(t1, t2); }
    /** @brief Return the difference or additive inverse. */
    friend T operator-(const T& t1, NT const& t2) { return sub(t1, t2); }
    /** @brief Return the product. */
    friend T operator*(const T& t1, const T& t2) { return mul(t1, t2); }
    /** @brief Return the quotient. */
    friend T operator/(const T& t1, NT const& t2) { return div(t1, t2); }
    /** @brief Add the operand in place. */
    friend T& operator+=(T& t1, const T& t2) { return t1 = add(t1, t2); }
    /** @brief Subtract the operand in place. */
    friend T& operator-=(T& t1, NT const& t2) { return t1 = sub(t1, t2); }
    /** @brief Multiply by the operand in place. */
    friend T& operator*=(T& t1, const T& t2) { return t1 = mul(t1, t2); }
    /** @brief Divide by the operand in place. */
    friend T& operator/=(T& t1, NT const& t2) { return t1 = div(t1, t2); }
};

/** @brief Implement inplace directed group operators for the supported operand types. */
template<class X, class Y, class NY, class R = X>
struct DefineInplaceDirectedGroupOperators {
};
/** @brief Implement inplace directed group operators for the supported operand types. */
template<class X, class Y, class NY>
struct DefineInplaceDirectedGroupOperators<X, Y, NY, X> {
    /** @brief Add the operand in place. */
    friend X& operator+=(X& x1, const Y& y2) { return x1 = operator+(x1, y2); }
    /** @brief Subtract the operand in place. */
    friend X& operator-=(X& x1, NY const& y2) { return x1 = operator-(x1, y2); }
};
/** @brief Implement directed group operators for the supported operand types. */
template<class X, class NX>
class DefineDirectedGroupOperators {
    // BUG: Cannot declare these due to bug in Clang
    // friend X operator-(NX const&);
    // friend NX operator-(NX const&, X const&);
    /** @brief Return the sum or unary positive value. */
    friend X operator+(const X& x) { return pos(x); }
    /** @brief Return the difference or additive inverse. */
    friend NX operator-(const X& x) { return neg(x); }
    /** @brief Return the sum or unary positive value. */
    friend X operator+(const X& x1, const X& x2) { return add(x1, x2); }
    /** @brief Return the difference or additive inverse. */
    friend X operator-(const X& x1, NX const& x2) { return sub(x1, x2); }
    /** @brief Add the operand in place. */
    friend X& operator+=(X& x1, const X& x2) { return x1 = add(x1, x2); }
    /** @brief Subtract the operand in place. */
    friend X& operator-=(X& x1, NX const& x2) { return x1 = sub(x1, x2); }
};

/** @brief Implement inplace directed semi field operators for the supported operand types. */
template<class X, class Y, class QY, class R = X>
struct DefineInplaceDirectedSemiFieldOperators {
};
/** @brief Implement inplace directed semi field operators for the supported operand types. */
template<class X, class Y, class QY>
struct DefineInplaceDirectedSemiFieldOperators<X, Y, QY, X> {
    /** @brief Add the operand in place. */
    friend X operator+=(X& x1, const Y& y2) { return x1 = operator+(x1, y2); }
    /** @brief Multiply by the operand in place. */
    friend X operator*=(X& x1, const Y& y2) { return x1 = operator*(x1, y2); }
    /** @brief Divide by the operand in place. */
    friend X operator/=(X& x1, const QY& y2) { return x1 = operator/(x1, y2); }
};
/** @brief Implement directed semi field operators for the supported operand types. */
template<class X, class QX, class R = X, class QR = QX>
struct DefineDirectedSemiFieldOperators
    : DefineInplaceDirectedSemiFieldOperators<X, X, QX, R> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const X& x2) { return add(x1, x2); }
    /** @brief Return the product. */
    friend R operator*(const X& x1, const X& x2) { return mul(x1, x2); }
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, const QX& qx2) { return div(x1, qx2); }
    // BUG: Cannot declare this due to bug in Clang
    // friend QR operator/(const QX& qx1, const X& x2);
};
/** @brief Implement mixed directed semi field operators for the supported operand types. */
template<class X, class QX, class Y, class QY, class R = X, class QR = QX>
struct DefineMixedDirectedSemiFieldOperators
    : DefineInplaceDirectedSemiFieldOperators<X, Y, QY, R> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const Y& y2) { return operator+(x1, factory(x1).create(y2)); }
    /** @brief Return the product. */
    friend R operator*(const X& x1, const Y& y2) { return operator*(x1, factory(x1).create(y2)); }
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, const QY& qy2) { return operator/(x1, create(qy2, x1)); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const Y& y1, const X& x2) { return operator+(factory(x2).create(y1), x2); }
    /** @brief Return the product. */
    friend R operator*(const Y& y1, const X& x2) { return operator*(factory(x2).create(y1), x2); }
    /** @brief Return the quotient. */
    friend QR operator/(const QY& y1, const X& x2) { return operator/(factory(x2).create(y1), x2); }
};

/** @brief Implement semi field operators for the supported operand types. */
template<class X, class R = X>
struct DefineSemiFieldOperators
    : DefineDirectedSemiFieldOperators<X, X, R, R> {
};
/** @brief Implement mixed semi field operators for the supported operand types. */
template<class X, class Y, class R = X>
struct DefineMixedSemiFieldOperators
    : DefineMixedDirectedSemiFieldOperators<X, X, Y, Y, R, R> {
};

/** @brief Implement inplace ring operators for the supported operand types. */
template<class X, class Y = X, class R = X>
struct DefineInplaceRingOperators {
};
/** @brief Implement inplace ring operators for the supported operand types. */
template<class X, class Y>
struct DefineInplaceRingOperators<X, Y, X> {
    /** @brief Add the operand in place. */
    friend X& operator+=(X& x1, const Y& y2) { return x1 = operator+(x1, y2); }
    /** @brief Subtract the operand in place. */
    friend X& operator-=(X& x1, const Y& y2) { return x1 = operator-(x1, y2); }
    /** @brief Multiply by the operand in place. */
    friend X& operator*=(X& x1, const Y& y2) { return x1 = operator*(x1, y2); }
};
/** @brief Implement ring operators for the supported operand types. */
template<class X, class R = X>
struct DefineRingOperators : DefineInplaceRingOperators<X, X, R> {
    /** @brief Return the sum or unary positive value. */
    friend X operator+(const X& x) { return pos(x); }
    /** @brief Return the difference or additive inverse. */
    friend X operator-(const X& x) { return neg(x); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const X& x2) { return add(x1, x2); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const X& x1, const X& x2) { return sub(x1, x2); }
    /** @brief Return the product. */
    friend R operator*(const X& x1, const X& x2) { return mul(x1, x2); }
};

/** @brief Implement inplace field operators for the supported operand types. */
template<class X, class Y = X, class R = X>
struct DefineInplaceFieldOperators {
};
/** @brief Implement inplace field operators for the supported operand types. */
template<class X, class Y>
struct DefineInplaceFieldOperators<X, Y, X> {
    /** @brief Add the operand in place. */
    friend X& operator+=(X& x1, const Y& y2) { return x1 = operator+(x1, y2); }
    /** @brief Subtract the operand in place. */
    friend X& operator-=(X& x1, const Y& y2) { return x1 = operator-(x1, y2); }
    /** @brief Multiply by the operand in place. */
    friend X& operator*=(X& x1, const Y& y2) { return x1 = operator*(x1, y2); }
    /** @brief Divide by the operand in place. */
    friend X& operator/=(X& x1, const Y& y2) { return x1 = operator/(x1, y2); }
};
/** @brief Implement field operators for the supported operand types. */
template<class X, class R = X>
struct DefineFieldOperators : DefineInplaceFieldOperators<X, X, R> {
    /** @brief Return the sum or unary positive value. */
    friend X operator+(const X& x) { return pos(x); }
    /** @brief Return the difference or additive inverse. */
    friend X operator-(const X& x) { return neg(x); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const X& x2) { return add(x1, x2); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const X& x1, const X& x2) { return sub(x1, x2); }
    /** @brief Return the product. */
    friend R operator*(const X& x1, const X& x2) { return mul(x1, x2); }
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, const X& x2) { return div(x1, x2); }
};

/** @brief Implement mixed field operators for the supported operand types. */
template<class X, class Y, class R = X>
struct DefineMixedFieldOperators : DefineInplaceFieldOperators<X, Y, R> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const Y& y2) { return add(x1, y2); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const X& x1, const Y& y2) { return sub(x1, y2); }
    /** @brief Return the product. */
    friend R operator*(const X& x1, const Y& y2) { return mul(x1, y2); }
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, const Y& y2) { return div(x1, y2); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const Y& y1, const X& x2) { return add(y1, x2); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const Y& y1, const X& x2) { return sub(y1, x2); }
    /** @brief Return the product. */
    friend R operator*(const Y& y1, const X& x2) { return mul(y1, x2); }
    /** @brief Return the quotient. */
    friend R operator/(const Y& y1, const X& x2) { return div(y1, x2); }
};

/** @brief Implement comparison operators for the supported operand types. */
template<class X, class LT, class EQ = LT>
struct DefineComparisonOperators {
    typedef LT GT;
    typedef decltype(not declval<EQ>()) NEQ;
    typedef decltype(not declval<LT>()) GEQ;
    typedef decltype(not declval<GT>()) LEQ;
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, const X& x2);
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, const X& x2);
    /** @brief Test equality of the operands. */
    friend EQ operator==(const X& x1, const X& x2) { return eq(x1, x2); }
    /** @brief Test inequality of the operands. */
    friend NEQ operator!=(const X& x1, const X& x2) { return not eq(x1, x2); }
    /** @brief Test whether the first operand is less than the second. */
    friend LT operator<(const X& x1, const X& x2) { return lt(x1, x2); }
    /** @brief Test whether the first operand is greater than the second. */
    friend GT operator>(const X& x1, const X& x2) { return lt(x2, x1); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend LEQ operator<=(const X& x1, const X& x2) { return not lt(x2, x1); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend GEQ operator>=(const X& x1, const X& x2) { return not lt(x1, x2); }
};

/** @brief Implement mixed comparison operators for the supported operand types. */
template<class X, class Y, class LT, class EQ = LT>
struct DefineMixedComparisonOperators {
    typedef LT GT;
    typedef decltype(not declval<EQ>()) NEQ;
    typedef decltype(not declval<LT>()) GEQ;
    typedef decltype(not declval<GT>()) LEQ;
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, const Y& y2);
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, const Y& y2);
    /** @brief Test whether the first operand is greater than the second. */
    friend GT gt(const X& x1, const Y& y2);
    /** @brief Test equality of the operands. */
    friend EQ operator==(const X& x1, const Y& y2) { return eq(x1, y2); }
    /** @brief Test inequality of the operands. */
    friend NEQ operator!=(const X& x1, const Y& y2) { return not eq(x1, y2); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend LEQ operator<=(const X& x1, const Y& y2) { return not gt(x1, y2); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend GEQ operator>=(const X& x1, const Y& y2) { return not lt(x1, y2); }
    /** @brief Test whether the first operand is less than the second. */
    friend LT operator<(const X& x1, const Y& y2) { return lt(x1, y2); }
    /** @brief Test whether the first operand is greater than the second. */
    friend GT operator>(const X& x1, const Y& y2) { return gt(x1, y2); }
    /** @brief Test equality of the operands. */
    friend EQ operator==(const Y& y1, const X& x2) { return eq(x2, y1); }
    /** @brief Test inequality of the operands. */
    friend NEQ operator!=(const Y& y1, const X& x2) { return not eq(x2, y1); }
    /** @brief Test whether the first operand is less than the second. */
    friend GT operator<(const Y& y1, const X& x2) { return gt(x2, y1); }
    /** @brief Test whether the first operand is greater than the second. */
    friend LT operator>(const Y& y1, const X& x2) { return lt(x2, y1); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend GEQ operator<=(const Y& y1, const X& x2) { return not lt(x2, y1); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend LEQ operator>=(const Y& y1, const X& x2) { return not gt(x2, y1); }
};

/** @brief Implement mixed comparison operators for the supported operand types. */
template<class X, class Y>
struct DefineMixedComparisonOperators<X, Y, Boolean> {
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const X& x1, const Y& y2);
    /** @brief Test equality of the operands. */
    friend Boolean eq(const X& x1, const Y& y2) { return cmp(x1, y2) == Comparison::EQUAL; }
    /** @brief Test whether the first operand is less than the second. */
    friend Boolean lt(const X& x1, const Y& y2) { return cmp(x1, y2) == Comparison::LESS; }
    /** @brief Test whether the first operand is greater than the second. */
    friend Boolean gt(const X& x1, const Y& y2) { return cmp(x1, y2) == Comparison::GREATER; }
    /** @brief Test equality of the operands. */
    friend Boolean operator==(const X& x1, const Y& y2) { return eq(x1, y2); }
    /** @brief Test inequality of the operands. */
    friend Boolean operator!=(const X& x1, const Y& y2) { return not eq(x1, y2); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend Boolean operator<=(const X& x1, const Y& y2) { return not gt(x1, y2); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend Boolean operator>=(const X& x1, const Y& y2) { return not lt(x1, y2); }
    /** @brief Test whether the first operand is less than the second. */
    friend Boolean operator<(const X& x1, const Y& y2) { return lt(x1, y2); }
    /** @brief Test whether the first operand is greater than the second. */
    friend Boolean operator>(const X& x1, const Y& y2) { return gt(x1, y2); }
    /** @brief Test equality of the operands. */
    friend Boolean operator==(const Y& y1, const X& x2) { return eq(x2, y1); }
    /** @brief Test inequality of the operands. */
    friend Boolean operator!=(const Y& y1, const X& x2) { return not eq(x2, y1); }
    /** @brief Test whether the first operand is less than the second. */
    friend Boolean operator<(const Y& y1, const X& x2) { return gt(x2, y1); }
    /** @brief Test whether the first operand is greater than the second. */
    friend Boolean operator>(const Y& y1, const X& x2) { return lt(x2, y1); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend Boolean operator<=(const Y& y1, const X& x2) { return not lt(x2, y1); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend Boolean operator>=(const Y& y1, const X& x2) { return not gt(x2, y1); }
};

/** @brief Implement directed comparison operators for the supported operand types. */
template<class X, class NX, class LT, class EQ>
struct DefineDirectedComparisonOperators {
    typedef decltype(not declval<EQ>()) NEQ;
    typedef decltype(not declval<LT>()) GT;
    typedef LT LEQ;
    typedef GT GEQ;
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, NX const& nx2);
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, NX const& nx2);
    /** @brief Test equality of the operands. */
    friend EQ operator==(const X& x1, NX const& nx2) { return eq(x1, nx2); }
    /** @brief Test inequality of the operands. */
    friend NEQ operator!=(const X& x1, NX const& nx2) { return not eq(x1, nx2); }
    /** @brief Test whether the first operand is less than the second. */
    friend LT operator<(const X& x1, NX const& nx2) { return lt(x1, nx2); }
    /** @brief Test whether the first operand is greater than the second. */
    friend LT operator>(NX const& nx1, const X& x2) { return lt(x2, nx1); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend GEQ operator>=(const X& x1, NX const& nx2) { return not lt(x1, nx2); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend GEQ operator<=(NX const& nx1, const X& x2) { return not lt(x2, nx1); }

    // BUG: Cannot declare these due to bug in Clang
    // friend EQ  operator==(NX const& nx1, X const& x2);
    // friend NEQ operator!=(NX const& nx1, X const& x2);
    // friend GT  operator< (NX const& nx1, X const& x2);
    // friend GT  operator> (X const& x1, NX const& nx2);
    // friend LEQ operator>=(NX const& nx1, X const& x2);
    // friend LEQ operator<=(X const& x1, NX const& nx2);
};

/** @brief Implement mixed directed comparison operators for the supported operand types. */
template<class X, class NY, class LT, class EQ>
struct DefineMixedDirectedComparisonOperators {
    typedef decltype(not declval<EQ>()) NEQ;
    typedef decltype(not declval<LT>()) GT;
    typedef LT LEQ;
    typedef GT GEQ;
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, NY const& y2);
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, NY const& y2);
    /** @brief Test whether the first operand is greater than the second. */
    friend GT gt(const X& x1, NY const& y2);
    /** @brief Test equality of the operands. */
    friend EQ operator==(const X& x1, NY const& ny2) { return eq(x1, ny2); }
    /** @brief Test inequality of the operands. */
    friend NEQ operator!=(const X& x1, NY const& ny2) { return not eq(x1, ny2); }
    /** @brief Test whether the first operand is less than the second. */
    friend LT operator<(const X& x1, NY const& ny2) { return lt(x1, ny2); }
    /** @brief Test whether the first operand is greater than the second. */
    friend GT operator>(const X& x1, NY const& ny2) { return gt(x1, ny2); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend LEQ operator<=(const X& x1, NY const& ny2) { return not gt(x1, ny2); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend GEQ operator>=(const X& x1, NY const& ny2) { return not lt(x1, ny2); }
    /** @brief Test equality of the operands. */
    friend EQ operator==(NY const& ny1, const X& x2) { return eq(x2, ny1); }
    /** @brief Test inequality of the operands. */
    friend NEQ operator!=(NY const& ny1, const X& x2) { return not eq(x2, ny1); }
    /** @brief Test whether the first operand is less than the second. */
    friend GT operator<(NY const& ny1, const X& x2) { return gt(x2, ny1); }
    /** @brief Test whether the first operand is greater than the second. */
    friend LT operator>(NY const& ny1, const X& x2) { return lt(x2, ny1); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend GEQ operator<=(NY const& ny1, const X& x2) { return not lt(x2, ny1); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend LEQ operator>=(NY const& ny1, const X& x2) { return not gt(x2, ny1); }
};

template<class X, class Y, class R = X>
struct ProvideConvertedFieldOperations;

/** @brief Provide converted field operations for the supported operand types. */
template<class X, class Y, class R>
struct ProvideConvertedFieldOperations {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(X x1, Y y2) { return operator+(R(x1), R(y2)); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(X x1, Y y2) { return operator-(R(x1), R(y2)); }
    /** @brief Return the product. */
    friend R operator*(X x1, Y y2) { return operator*(R(x1), R(y2)); }
    /** @brief Return the quotient. */
    friend R operator/(X x1, Y y2) { return operator/(R(x1), R(y2)); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(Y y1, X x2) { return operator+(R(y1), R(x2)); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(Y y1, X x2) { return operator-(R(y1), R(x2)); }
    /** @brief Return the product. */
    friend R operator*(Y y1, X x2) { return operator*(R(y1), R(x2)); }
    /** @brief Return the quotient. */
    friend R operator/(Y y1, X x2) { return operator/(R(y1), R(x2)); }
    /** @brief Return the sum of the operands. */
    friend R add(X x1, Y y2) { return add(R(x1), R(y2)); }
    /** @brief Return the difference of the operands. */
    friend R sub(X x1, Y y2) { return sub(R(x1), R(y2)); }
    /** @brief Return the product of the operands. */
    friend R mul(X x1, Y y2) { return mul(R(x1), R(y2)); }
    /** @brief Return the quotient of the operands. */
    friend R div(X x1, Y y2) { return div(R(x1), R(y2)); }
    /** @brief Return the sum of the operands. */
    friend R add(Y y1, X x2) { return add(R(y1), R(x2)); }
    /** @brief Return the difference of the operands. */
    friend R sub(Y y1, X x2) { return sub(R(y1), R(x2)); }
    /** @brief Return the product of the operands. */
    friend R mul(Y y1, X x2) { return mul(R(y1), R(x2)); }
    /** @brief Return the quotient of the operands. */
    friend R div(Y y1, X x2) { return div(R(y1), R(x2)); }
};

/** @brief Provide converted field operations for the supported operand types. */
template<class X, class Y>
struct ProvideConvertedFieldOperations<X, Y, X>
    : DefineInplaceFieldOperators<X, Y, X> {
    typedef X R;
    /** @brief Return the sum or unary positive value. */
    friend R operator+(X x1, Y y2) { return operator+(x1, R(y2)); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(X x1, Y y2) { return operator-(x1, R(y2)); }
    /** @brief Return the product. */
    friend R operator*(X x1, Y y2) { return operator*(x1, R(y2)); }
    /** @brief Return the quotient. */
    friend R operator/(X x1, Y y2) { return operator/(x1, R(y2)); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(Y y1, X x2) { return operator+(R(y1), x2); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(Y y1, X x2) { return operator-(R(y1), x2); }
    /** @brief Return the product. */
    friend R operator*(Y y1, X x2) { return operator*(R(y1), x2); }
    /** @brief Return the quotient. */
    friend R operator/(Y y1, X x2) { return operator/(R(y1), x2); }
    /** @brief Return the sum of the operands. */
    friend R add(X x1, Y y2) { return add(x1, R(y2)); }
    /** @brief Return the difference of the operands. */
    friend R sub(X x1, Y y2) { return sub(x1, R(y2)); }
    /** @brief Return the product of the operands. */
    friend R mul(X x1, Y y2) { return mul(x1, R(y2)); }
    /** @brief Return the quotient of the operands. */
    friend R div(X x1, Y y2) { return div(x1, R(y2)); }
    /** @brief Return the sum of the operands. */
    friend R add(Y y1, X x2) { return add(R(y1), x2); }
    /** @brief Return the difference of the operands. */
    friend R sub(Y y1, X x2) { return sub(R(y1), x2); }
    /** @brief Return the product of the operands. */
    friend R mul(Y y1, X x2) { return mul(R(y1), x2); }
    /** @brief Return the quotient of the operands. */
    friend R div(Y y1, X x2) { return div(R(y1), x2); }
};

/** @brief Provide converted field operations for the supported operand types. */
template<class X, class R>
struct ProvideConvertedFieldOperations<X, X, R> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const X& x2) { return operator+(R(x1), R(x2)); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const X& x1, const X& x2) { return operator-(R(x1), R(x2)); }
    /** @brief Return the product. */
    friend R operator*(const X& x1, const X& x2) { return operator*(R(x1), R(x2)); }
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, const X& x2) { return operator/(R(x1), R(x2)); }
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const X& x2) { return add(R(x1), R(x2)); }
    /** @brief Return the difference of the operands. */
    friend R sub(const X& x1, const X& x2) { return sub(R(x1), R(x2)); }
    /** @brief Return the product of the operands. */
    friend R mul(const X& x1, const X& x2) { return mul(R(x1), R(x2)); }
    /** @brief Return the quotient of the operands. */
    friend R div(const X& x1, const X& x2) { return div(R(x1), R(x2)); }
};

/** @brief Provide concrete generic directed group operators for the supported operand types. */
template<class X, class NX, class Y, class NY, class R = X, class NR = NX>
struct ProvideConcreteGenericDirectedGroupOperators
    : DefineInplaceDirectedGroupOperators<X, Y, NY, R> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const Y& y2) { return operator+(x1, x1.create(y2)); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const X& x1, const NY& ny2) { return operator-(x1, create(ny2, x1)); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const Y& y1, const X& x2) { return operator+(x2.create(y1), x2); }
    /** @brief Return the difference or additive inverse. */
    friend NR operator-(const NY& ny1, const X& x2) { return operator-(x2.create(ny1), x2); }
};
/** @brief Provide concrete generic directed group operations for the supported operand types. */
template<class X, class NX, class Y, class NY, class R = X, class NR = NX>
struct ProvideConcreteGenericDirectedGroupOperations
    : ProvideConcreteGenericDirectedGroupOperators<X, NX, Y, NY, R, NR> {
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const Y& y2) { return add(x1, x1.create(y2)); }
    /** @brief Return the difference of the operands. */
    friend R sub(const X& x1, const NY& ny2) { return sub(x1, x1.create(ny2)); }
    /** @brief Return the sum of the operands. */
    friend R add(const Y& y1, const X& x2) { return add(x2.create(y1), x2); }
    /** @brief Return the difference of the operands. */
    friend NR sub(const NY& ny1, const X& x2) { return sub(x2.create(ny1), x2); }
};

/** @brief Provide concrete generic directed semi field operators for the supported operand types. */
template<class X, class QX, class Y, class QY, class R = X, class QR = QX>
struct ProvideConcreteGenericDirectedSemiFieldOperators
    : DefineInplaceDirectedSemiFieldOperators<X, Y, QY, R> {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const Y& y2) { return operator+(x1, x1.create(y2)); }
    /** @brief Return the product. */
    friend R operator*(const X& x1, const Y& y2) { return operator*(x1, x1.create(y2)); }
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, const QY& qy2) { return operator/(x1, x1.create(qy2)); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const Y& y1, const X& x2) { return operator+(x2.create(y1), x2); }
    /** @brief Return the product. */
    friend R operator*(const Y& y1, const X& x2) { return operator*(x2.create(y1), x2); }
    /** @brief Return the quotient. */
    friend QR operator/(const QY& y1, const X& x2) { return operator/(x2.create(y1), x2); }
};
/** @brief Provide concrete generic directed semi field operations for the supported operand types. */
template<class X, class QX, class Y, class QY, class R = X, class QR = QX>
struct ProvideConcreteGenericDirectedSemiFieldOperations
    : ProvideConcreteGenericDirectedSemiFieldOperators<X, QX, Y, QY, R, QR> {
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const Y& y2) { return add(x1, x1.create(y2)); }
    /** @brief Return the product of the operands. */
    friend R mul(const X& x1, const Y& y2) { return mul(x1, x1.create(y2)); }
    /** @brief Return the quotient of the operands. */
    friend R div(const X& x1, const QY& y2) { return div(x1, x1.create(y2)); }
    /** @brief Return the sum of the operands. */
    friend R add(const Y& y1, const X& x2) { return add(x2.create(y1), x2); }
    /** @brief Return the product of the operands. */
    friend R mul(const Y& y1, const X& x2) { return mul(x2.create(y1), x2); }
    /** @brief Return the quotient of the operands. */
    friend QR div(const QY& y1, const X& x2) { return div(x2.create(y1), x2); }
};

/** @brief Provide concrete generic field operators for the supported operand types. */
template<class X, class Y, class R = X>
struct ProvideConcreteGenericFieldOperators
    : DefineInplaceFieldOperators<X, Y, R> {
    // static R create(Y const& y, X const& x) { return X(y,x.precision()); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    friend decltype(auto) _create(const Y& y, const X& x) { return factory(x).create(y); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const Y& y2) { return operator+(x1, _create(y2, x1)); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const X& x1, const Y& y2) { return operator-(x1, _create(y2, x1)); }
    /** @brief Return the product. */
    friend R operator*(const X& x1, const Y& y2) { return operator*(x1, _create(y2, x1)); }
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, const Y& y2) { return operator/(x1, _create(y2, x1)); }
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const Y& y1, const X& x2) { return operator+(_create(y1, x2), x2); }
    /** @brief Return the difference or additive inverse. */
    friend R operator-(const Y& y1, const X& x2) { return operator-(_create(y1, x2), x2); }
    /** @brief Return the product. */
    friend R operator*(const Y& y1, const X& x2) { return operator*(_create(y1, x2), x2); }
    /** @brief Return the quotient. */
    friend R operator/(const Y& y1, const X& x2) { return operator/(_create(y1, x2), x2); }
};
/** @brief Provide concrete generic field operations for the supported operand types. */
template<class X, class Y, class R = X>
struct ProvideConcreteGenericFieldOperations : ProvideConcreteGenericFieldOperators<X, Y, R> {
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const Y& y2) { return add(x1, _create(y2, x1)); }
    /** @brief Return the difference of the operands. */
    friend R sub(const X& x1, const Y& y2) { return sub(x1, _create(y2, x1)); }
    /** @brief Return the product of the operands. */
    friend R mul(const X& x1, const Y& y2) { return mul(x1, _create(y2, x1)); }
    /** @brief Return the quotient of the operands. */
    friend R div(const X& x1, const Y& y2) { return div(x1, _create(y2, x1)); }
    /** @brief Return the sum of the operands. */
    friend R add(const Y& y1, const X& x2) { return add(_create(y1, x2), x2); }
    /** @brief Return the difference of the operands. */
    friend R sub(const Y& y1, const X& x2) { return sub(_create(y1, x2), x2); }
    /** @brief Return the product of the operands. */
    friend R mul(const Y& y1, const X& x2) { return mul(_create(y1, x2), x2); }
    /** @brief Return the quotient of the operands. */
    friend R div(const Y& y1, const X& x2) { return div(_create(y1, x2), x2); }
};

/** @brief Provide concrete generic field operators for the supported operand types. */
template<class X, class Y>
struct ProvideConcreteGenericFieldOperators<X, Y, Void> {
    /** @brief Return the sum or unary positive value. */
    friend decltype(auto) operator+(const X& x, const Y& y) { return add(x, factory(x).create(y)); }
    /** @brief Return the difference or additive inverse. */
    friend decltype(auto) operator-(const X& x, const Y& y) { return sub(x, factory(x).create(y)); }
    /** @brief Return the product. */
    friend decltype(auto) operator*(const X& x, const Y& y) { return mul(x, factory(x).create(y)); }
    /** @brief Return the quotient. */
    friend decltype(auto) operator/(const X& x, const Y& y) { return div(x, factory(x).create(y)); }
    /** @brief Return the sum or unary positive value. */
    friend decltype(auto) operator+(const Y& y, const X& x) { return add(factory(x).create(y), x); }
    /** @brief Return the difference or additive inverse. */
    friend decltype(auto) operator-(const Y& y, const X& x) { return sub(factory(x).create(y), x); }
    /** @brief Return the product. */
    friend decltype(auto) operator*(const Y& y, const X& x) { return mul(factory(x).create(y), x); }
    /** @brief Return the quotient. */
    friend decltype(auto) operator/(const Y& y, const X& x) { return div(factory(x).create(y), x); }
    /** @brief Add the operand in place. */
    friend decltype(auto) operator+=(X& x, const Y& y) { return x = add(x, factory(x).create(y)); }
    /** @brief Subtract the operand in place. */
    friend decltype(auto) operator-=(X& x, const Y& y) { return x = sub(x, factory(x).create(y)); }
    /** @brief Multiply by the operand in place. */
    friend decltype(auto) operator*=(X& x, const Y& y) { return x = mul(x, factory(x).create(y)); }
    /** @brief Divide by the operand in place. */
    friend decltype(auto) operator/=(X& x, const Y& y) { return x = div(x, factory(x).create(y)); }
};
/** @brief Provide concrete generic field operations for the supported operand types. */
template<class X, class Y>
struct ProvideConcreteGenericFieldOperations<X, Y, Void> : ProvideConcreteGenericFieldOperators<X, Y, Void> {
    /** @brief Return the sum of the operands. */
    friend decltype(auto) add(const X& x, const Y& y) { return add(x, factory(x).create(y)); }
    /** @brief Return the difference of the operands. */
    friend decltype(auto) sub(const X& x, const Y& y) { return sub(x, factory(x).create(y)); }
    /** @brief Return the product of the operands. */
    friend decltype(auto) mul(const X& x, const Y& y) { return mul(x, factory(x).create(y)); }
    /** @brief Return the quotient of the operands. */
    friend decltype(auto) div(const X& x, const Y& y) { return div(x, factory(x).create(y)); }
    /** @brief Return the sum of the operands. */
    friend decltype(auto) add(const Y& y, const X& x) { return add(factory(x).create(y), x); }
    /** @brief Return the difference of the operands. */
    friend decltype(auto) sub(const Y& y, const X& x) { return sub(factory(x).create(y), x); }
    /** @brief Return the product of the operands. */
    friend decltype(auto) mul(const Y& y, const X& x) { return mul(factory(x).create(y), x); }
    /** @brief Return the quotient of the operands. */
    friend decltype(auto) div(const Y& y, const X& x) { return div(factory(x).create(y), x); }
};

/** @brief Provide concrete generic field operators for the supported operand types. */
template<class X>
struct ProvideConcreteGenericFieldOperators<X, Void> {
    /** @brief Return the sum or unary positive value. */
    template<GenericNumber Y>
    friend decltype(auto) operator+(const X& x, const Y& y) { return add(x, factory(x).create(y)); }
    /** @brief Return the difference or additive inverse. */
    template<GenericNumber Y>
    friend decltype(auto) operator-(const X& x, const Y& y) { return sub(x, factory(x).create(y)); }
    /** @brief Return the product. */
    template<GenericNumber Y>
    friend decltype(auto) operator*(const X& x, const Y& y) { return mul(x, factory(x).create(y)); }
    /** @brief Return the quotient. */
    template<GenericNumber Y>
    friend decltype(auto) operator/(const X& x, const Y& y) { return div(x, factory(x).create(y)); }
    /** @brief Return the sum or unary positive value. */
    template<GenericNumber Y>
    friend decltype(auto) operator+(const Y& y, const X& x) { return add(factory(x).create(y), x); }
    /** @brief Return the difference or additive inverse. */
    template<GenericNumber Y>
    friend decltype(auto) operator-(const Y& y, const X& x) { return sub(factory(x).create(y), x); }
    /** @brief Return the product. */
    template<GenericNumber Y>
    friend decltype(auto) operator*(const Y& y, const X& x) { return mul(factory(x).create(y), x); }
    /** @brief Return the quotient. */
    template<GenericNumber Y>
    friend decltype(auto) operator/(const Y& y, const X& x) { return div(factory(x).create(y), x); }
    /** @brief Add the operand in place. */
    template<GenericNumber Y>
    friend decltype(auto) operator+=(X& x, const Y& y) { return x = add(x, factory(x).create(y)); }
    /** @brief Subtract the operand in place. */
    template<GenericNumber Y>
    friend decltype(auto) operator-=(X& x, const Y& y) { return x = sub(x, factory(x).create(y)); }
    /** @brief Multiply by the operand in place. */
    template<GenericNumber Y>
    friend decltype(auto) operator*=(X& x, const Y& y) { return x = mul(x, factory(x).create(y)); }
    /** @brief Divide by the operand in place. */
    template<GenericNumber Y>
    friend decltype(auto) operator/=(X& x, const Y& y) { return x = div(x, factory(x).create(y)); }
};
/** @brief Provide concrete generic field operations for the supported operand types. */
template<class X>
struct ProvideConcreteGenericFieldOperations<X, Void> : ProvideConcreteGenericFieldOperators<X, Void> {
    /** @brief Return the sum of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) add(const X& x, const Y& y) { return add(x, factory(x).create(y)); }
    /** @brief Return the difference of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) sub(const X& x, const Y& y) { return sub(x, factory(x).create(y)); }
    /** @brief Return the product of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) mul(const X& x, const Y& y) { return mul(x, factory(x).create(y)); }
    /** @brief Return the quotient of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) div(const X& x, const Y& y) { return div(x, factory(x).create(y)); }
    /** @brief Return the sum of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) add(const Y& y, const X& x) { return add(factory(x).create(y), x); }
    /** @brief Return the difference of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) sub(const Y& y, const X& x) { return sub(factory(x).create(y), x); }
    /** @brief Return the product of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) mul(const Y& y, const X& x) { return mul(factory(x).create(y), x); }
    /** @brief Return the quotient of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) div(const Y& y, const X& x) { return div(factory(x).create(y), x); }
};

/** @brief Provide concrete generic arithmetic operators for the supported operand types. */
template<class X, class Y = Void, class R = X>
struct ProvideConcreteGenericArithmeticOperators : ProvideConcreteGenericFieldOperators<X, Y, R> {};
/** @brief Provide concrete generic arithmetic operations for the supported operand types. */
template<class X, class Y = Void, class R = X>
struct ProvideConcreteGenericArithmeticOperations : ProvideConcreteGenericFieldOperations<X, Y, R> {};

/** @brief Provide concrete generic lattice operations for the supported operand types. */
template<class X, class Y = Void, class R = X>
struct ProvideConcreteGenericLatticeOperations {
    /** @brief Return the maximum of the operands. */
    friend R max(const X& x, const Y& y) { return max(x, factory(x).create(y)); }
    /** @brief Return the minimum of the operands. */
    friend R min(const X& x, const Y& y) { return min(x, factory(x).create(y)); }
    /** @brief Return the maximum of the operands. */
    friend R max(const Y& y, const X& x) { return max(factory(x).create(y), x); }
    /** @brief Return the minimum of the operands. */
    friend R min(const Y& y, const X& x) { return min(factory(x).create(y), x); }
};

/** @brief Provide concrete generic lattice operations for the supported operand types. */
template<class X, class Y>
struct ProvideConcreteGenericLatticeOperations<X, Y, Void> {
    /** @brief Return the maximum of the operands. */
    friend decltype(auto) max(const X& x, const Y& y) { return max(x, factory(x).create(y)); }
    /** @brief Return the minimum of the operands. */
    friend decltype(auto) min(const X& x, const Y& y) { return min(x, factory(x).create(y)); }
    /** @brief Return the maximum of the operands. */
    friend decltype(auto) max(const Y& y, const X& x) { return max(factory(x).create(y), x); }
    /** @brief Return the minimum of the operands. */
    friend decltype(auto) min(const Y& y, const X& x) { return min(factory(x).create(y), x); }
};

/** @brief Provide concrete generic elementary operations for the supported operand types. */
template<class X, class Y = Void, class R = X>
struct ProvideConcreteGenericElementaryOperations
    : ProvideConcreteGenericArithmeticOperations<X, Y, R>,
      ProvideConcreteGenericLatticeOperations<X, Y, R> {
};

/** @brief Provide converted comparison operations for the supported operand types. */
template<class X, class Y, class R, class LT, class EQ = LT>
struct ProvideConvertedComparisonOperations : DefineMixedComparisonOperators<X, Y, LT, EQ> {
    typedef LT GT;
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, const Y& y2) { return lt(R(x1), R(y2)); }
    /** @brief Test whether the first operand is greater than the second. */
    friend GT gt(const X& x1, const Y& y2) { return lt(R(y2), R(x1)); }
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, const Y& y2) { return eq(R(x1), R(y2)); }
};

/** @brief Provide concrete generic comparison operations for the supported operand types. */
template<class X, class Y, class LT, class EQ = LT>
struct ProvideConcreteGenericComparisonOperations : DefineMixedComparisonOperators<X, Y, LT, EQ> {
    typedef LT GT;
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, const Y& y2) { return lt(x1, factory(x1).create(y2)); }
    /** @brief Test whether the first operand is greater than the second. */
    friend GT gt(const X& x1, const Y& y2) { return lt(factory(x1).create(y2), x1); }
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, const Y& y2) { return eq(x1, factory(x1).create(y2)); }
};

/** @brief Provide concrete generic directed comparison operations for the supported operand types. */
template<class X, class NY, class LT, class EQ>
struct ProvideConcreteGenericDirectedComparisonOperations : DefineMixedDirectedComparisonOperators<X, NY, LT, EQ> {
    typedef decltype(not declval<LT>()) GT;
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, const NY& ny2) { return eq(x1, factory(x1).create(ny2)); }
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, const NY& ny2) { return lt(x1, factory(x1).create(ny2)); }
    /** @brief Test whether the first operand is greater than the second. */
    friend GT gt(const X& x1, const NY& ny2) { return lt(factory(x1).create(ny2), x1); }
};

/** @brief Provide concrete generic operations for the supported operand types. */
template<class X, class Y, class R, class LT, class EQ = LT>
struct ProvideConcreteGenericOperations
    : ProvideConcreteGenericFieldOperations<X, Y, R>,
      ProvideConcreteGenericComparisonOperations<X, Y, LT, EQ> {};

/** @brief Provide concrete generic directed operations for the supported operand types. */
template<class X, class NX, class Y, class NY, class LT, class EQ>
struct ProvideConcreteGenericDirectedOperations
    : ProvideConcreteGenericDirectedGroupOperations<X, NX, Y, NY>,
      ProvideConcreteGenericDirectedComparisonOperations<X, NY, LT, EQ> {};

/** @brief Provide directed semi field operators for the supported operand types. */
template<class X, class QX = X, class R = X, class QR = QX>
class ProvideDirectedSemiFieldOperators {
    /** @brief Return the sum or unary positive value. */
    friend R operator+(const X& x1, const X& x2) { return add(x1, x2); }
    /** @brief Return the product. */
    friend R operator*(const X& x1, const X& x2) { return mul(x1, x2); }
    /** @brief Return the quotient. */
    friend R operator/(const X& x1, QX const& x2) { return div(x1, x2); }
    /** @brief Add the operand in place. */
    friend X& operator+=(X& x1, const X& x2) { return x1 = add(x1, x2); }
    /** @brief Multiply by the operand in place. */
    friend X& operator*=(X& x1, const X& x2) { return x1 = mul(x1, x2); }
    /** @brief Divide by the operand in place. */
    friend X& operator/=(X& x1, QX const& x2) { return x1 = div(x1, x2); }

    // BUG: Cannot declare these due to bug in Clang
    // friend R rec(QX const&);
    // friend QR rec(X const&);
    // friend QR operator/(QX const&, X const&);
};

// Evaluate generic division viability only when both operand types are
// complete. Checking the conversion in a trailing return type eagerly
// instantiates CRTP bases and breaks positive/directed numeric classes.
template<class X, class Y>
concept CanForwardGenericDivisionLeft =
    requires(const X& x, const Y& y) { div(x, factory(x).create(y)); };
template<class X, class Y>
concept CanForwardGenericDivisionRight =
    requires(const X& x, const Y& y) { div(factory(x).create(y), x); };

/** @brief Implement concrete generic arithmetic operations for the supported operand types. */
template<class X>
struct DefineConcreteGenericArithmeticOperations {
    /** @brief Return the sum of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) add(const X& x, const Y& y) { return add(x, factory(x).create(y)); }
    /** @brief Return the difference of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) sub(const X& x, const Y& y) { return sub(x, factory(x).create(y)); }
    /** @brief Return the product of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) mul(const X& x, const Y& y) { return mul(x, factory(x).create(y)); }
    /** @brief Return the quotient of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) div(const X& x, const Y& y)
        requires CanForwardGenericDivisionLeft<X, Y>
    { return div(x, factory(x).create(y)); }

    /** @brief Return the sum of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) add(const Y& y, const X& x) { return add(factory(x).create(y), x); }
    /** @brief Return the difference of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) sub(const Y& y, const X& x) { return sub(factory(x).create(y), x); }
    /** @brief Return the product of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) mul(const Y& y, const X& x) { return mul(factory(x).create(y), x); }
    /** @brief Return the quotient of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) div(const Y& y, const X& x)
        requires CanForwardGenericDivisionRight<X, Y>
    { return div(factory(x).create(y), x); }

    /** @brief Return the maximum of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) max(const X& x, const Y& y) { return max(x, factory(x).create(y)); }
    /** @brief Return the minimum of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) min(const X& x, const Y& y) { return min(x, factory(x).create(y)); }
    /** @brief Return the maximum of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) max(const Y& y, const X& x) { return max(factory(x).create(y), x); }
    /** @brief Return the minimum of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) min(const Y& y, const X& x) { return min(factory(x).create(y), x); }
};

/** @brief Implement concrete generic arithmetic operators for the supported operand types. */
template<class X>
struct DefineConcreteGenericArithmeticOperators
    : DefineConcreteGenericArithmeticOperations<X> {
    /** @brief Return the sum or unary positive value. */
    template<GenericNumber Y>
    friend decltype(auto) operator+(const X& x, const Y& y) { return x + factory(x).create(y); }
    /** @brief Return the difference or additive inverse. */
    template<GenericNumber Y>
    friend decltype(auto) operator-(const X& x, const Y& y) { return x - factory(x).create(y); }
    /** @brief Return the product. */
    template<GenericNumber Y>
    friend decltype(auto) operator*(const X& x, const Y& y) { return x * factory(x).create(y); }
    /** @brief Return the quotient. */
    template<GenericNumber Y>
    friend decltype(auto) operator/(const X& x, const Y& y) { return x / factory(x).create(y); }

    /** @brief Return the sum or unary positive value. */
    template<GenericNumber Y>
    friend decltype(auto) operator+(const Y& y, const X& x) { return factory(x).create(y) + x; }
    /** @brief Return the difference or additive inverse. */
    template<GenericNumber Y>
    friend decltype(auto) operator-(const Y& y, const X& x) { return factory(x).create(y) - x; }
    /** @brief Return the product. */
    template<GenericNumber Y>
    friend decltype(auto) operator*(const Y& y, const X& x) { return factory(x).create(y) * x; }
    /** @brief Return the quotient. */
    template<GenericNumber Y>
    friend decltype(auto) operator/(const Y& y, const X& x) { return factory(x).create(y) / x; }

    /** @brief Add the operand in place. */
    template<GenericNumber Y>
    friend decltype(auto) operator+=(X& x, const Y& y) { return x += factory(x).create(y); }
    /** @brief Subtract the operand in place. */
    template<GenericNumber Y>
    friend decltype(auto) operator-=(X& x, const Y& y) { return x -= factory(x).create(y); }
    /** @brief Multiply by the operand in place. */
    template<GenericNumber Y>
    friend decltype(auto) operator*=(X& x, const Y& y) { return x *= factory(x).create(y); }
    /** @brief Divide by the operand in place. */
    template<GenericNumber Y>
    friend decltype(auto) operator/=(X& x, const Y& y) { return x /= factory(x).create(y); }
};

/** @brief Implement concrete generic comparison operators for the supported operand types. */
template<class X>
struct DefineConcreteGenericComparisonOperators {
    /** @brief Test equality of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) operator==(const X& x, const Y& y) { return x == factory(x).create(y); }
    /** @brief Test inequality of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) operator!=(const X& x, const Y& y) { return x != factory(x).create(y); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    template<GenericNumber Y>
    friend decltype(auto) operator<=(const X& x, const Y& y) { return x <= factory(x).create(y); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    template<GenericNumber Y>
    friend decltype(auto) operator>=(const X& x, const Y& y) { return x >= factory(x).create(y); }
    /** @brief Test whether the first operand is less than the second. */
    template<GenericNumber Y>
    friend decltype(auto) operator<(const X& x, const Y& y) { return x < factory(x).create(y); }
    /** @brief Test whether the first operand is greater than the second. */
    template<GenericNumber Y>
    friend decltype(auto) operator>(const X& x, const Y& y) { return x > factory(x).create(y); }

    /** @brief Test equality of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) operator==(const Y& y, const X& x) { return factory(x).create(y) == x; }
    /** @brief Test inequality of the operands. */
    template<GenericNumber Y>
    friend decltype(auto) operator!=(const Y& y, const X& x) { return factory(x).create(y) != x; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    template<GenericNumber Y>
    friend decltype(auto) operator<=(const Y& y, const X& x) { return factory(x).create(y) <= x; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    template<GenericNumber Y>
    friend decltype(auto) operator>=(const Y& y, const X& x) { return factory(x).create(y) >= x; }
    /** @brief Test whether the first operand is less than the second. */
    template<GenericNumber Y>
    friend decltype(auto) operator<(const Y& y, const X& x) { return factory(x).create(y) < x; }
    /** @brief Test whether the first operand is greater than the second. */
    template<GenericNumber Y>
    friend decltype(auto) operator>(const Y& y, const X& x) { return factory(x).create(y) > x; }
};

/** @brief Implement concrete generic operators for the supported operand types. */
template<class X>
struct DefineConcreteGenericOperators
    : DefineConcreteGenericArithmeticOperators<X>,
      DefineConcreteGenericComparisonOperators<X> {};

/** @brief Compile-time traits for numeric. */
template<class T>
struct NumericTraits;

/** @brief Operations supported by  objects. */
template<class X>
class Operations {
    typedef decltype(add(declval<X>(), declval<X>())) R;
    typedef decltype(neg(declval<X>())) NX;
    typedef typename NumericTraits<X>::OppositeType QX;
    typedef decltype(rec(declval<X>())) QR;
    typedef R PR;
    typedef decltype(abs(declval<X>())) PX;
    typedef decltype(eq(declval<X>(), declval<NX>())) EQ;
    typedef decltype(lt(declval<X>(), declval<NX>())) LT;

  public:
    /** @brief Return the additive zero with the operand characteristics. */
    static X _nul(const X& x);
    /** @brief Return the unary positive value. */
    static X _pos(const X& x);
    /** @brief Return the additive inverse. */
    static NX _neg(const X& x);
    /** @brief Divide the value by two. */
    static X _hlf(const X& x);
    /** @brief Return the square. */
    static PR _sqr(const X& x);
    /** @brief Return the multiplicative inverse. */
    static QR _rec(const X& x);
    /** @brief Return the sum of the operands. */
    static R _add(const X& x1, const X& x2);
    /** @brief Return the difference of the operands. */
    static R _sub(const X& x1, NX const& x2);
    /** @brief Return the product of the operands. */
    static R _mul(const X& x1, const X& x2);
    /** @brief Return the quotient of the operands. */
    static R _div(const X& x1, QX const& x2);
    /** @brief Raise the value to the specified power. */
    static R _pow(const X& x, Nat m);
    /** @brief Raise the value to the specified power. */
    static R _pow(const X& x, Int n);
    /** @brief Return the square root. */
    static R _sqrt(const X& x);
    /** @brief Return the exponential. */
    static R _exp(const X& x);
    /** @brief Return the natural logarithm. */
    static R _log(const X& x);
    /** @brief Return the sine. */
    static R _sin(const X& x);
    /** @brief Return the cosine. */
    static R _cos(const X& x);
    /** @brief Return the tangent. */
    static R _tan(const X& x);
    /** @brief Return the inverse sine. */
    static R _asin(const X& x);
    /** @brief Return the inverse cosine. */
    static R _acos(const X& x);
    /** @brief Return the inverse tangent. */
    static R _atan(const X& x);
    /** @brief Return the maximum of the operands. */
    static X _max(const X& x1, const X& x2);
    /** @brief Return the minimum of the operands. */
    static X _min(const X& x1, const X& x2);
    /** @brief Return the absolute value. */
    static PX _abs(const X& x);
    /** @brief Test equality of the operands. */
    static EQ _eq(const X& x1, NX const& x2);
    /** @brief Test whether the first operand is less than the second. */
    static LT _lt(const X& x1, NX const& x2);
    /** @brief Write the value to the output stream. */
    static OutputStream& _write(OutputStream& os, const X& x);
    /** @brief Read the value from the input stream. */
    static InputStream& _read(InputStream& is, X& x);
};

/** @brief Dispatch numeric operations for the supported operand types. */
template<class X, class R = X>
struct DispatchNumericOperations
    : DefineFieldOperators<X, R> {
    typedef Operations<X> OperationsType;
    /** @brief Return the additive zero with the operand characteristics. */
    friend X nul(const X& x) { return OperationsType::_nul(x); }
    /** @brief Return the unary positive value. */
    friend X pos(const X& x) { return OperationsType::_pos(x); }
    /** @brief Return the additive inverse. */
    friend X neg(const X& x) { return OperationsType::_neg(x); }
    /** @brief Divide the value by two. */
    friend X hlf(const X& x) { return OperationsType::_hlf(x); }
    /** @brief Return the square. */
    friend R sqr(const X& x) { return OperationsType::_sqr(x); }
    /** @brief Return the multiplicative inverse. */
    friend R rec(const X& x) { return OperationsType::_rec(x); }
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const X& x2) { return OperationsType::_add(x1, x2); }
    /** @brief Return the difference of the operands. */
    friend R sub(const X& x1, const X& x2) { return OperationsType::_sub(x1, x2); }
    /** @brief Return the product of the operands. */
    friend R mul(const X& x1, const X& x2) { return OperationsType::_mul(x1, x2); }
    /** @brief Return the quotient of the operands. */
    friend R div(const X& x1, const X& x2) { return OperationsType::_div(x1, x2); }
    /** @brief Raise the value to the specified power. */
    friend R pow(const X& x, Nat m) { return OperationsType::_pow(x, m); }
    /** @brief Raise the value to the specified power. */
    friend R pow(const X& x, Int n) { return OperationsType::_pow(x, n); }
    /** @brief Return the square root. */
    friend R sqrt(const X& x) { return OperationsType::_sqrt(x); }
    /** @brief Return the exponential. */
    friend R exp(const X& x) { return OperationsType::_exp(x); }
    /** @brief Return the natural logarithm. */
    friend R log(const X& x) { return OperationsType::_log(x); }
    /** @brief Return the sine. */
    friend R sin(const X& x) { return OperationsType::_sin(x); }
    /** @brief Return the cosine. */
    friend R cos(const X& x) { return OperationsType::_cos(x); }
    /** @brief Return the tangent. */
    friend R tan(const X& x) { return OperationsType::_tan(x); }
    /** @brief Return the inverse sine. */
    friend R asin(const X& x) { return OperationsType::_asin(x); }
    /** @brief Return the inverse cosine. */
    friend R acos(const X& x) { return OperationsType::_acos(x); }
    /** @brief Return the inverse tangent. */
    friend R atan(const X& x) { return OperationsType::_atan(x); }
    /** @brief Return the maximum of the operands. */
    friend X max(const X& x1, const X& x2) { return OperationsType::_max(x1, x2); }
    /** @brief Return the minimum of the operands. */
    friend X min(const X& x1, const X& x2) { return OperationsType::_min(x1, x2); }
    /** @brief Return the absolute value. */
    friend X abs(const X& x) { return OperationsType::_abs(x); }
};

/** @brief Dispatch directed numeric operations for the supported operand types. */
template<class X, class NX, class R = X, class NR = NX>
struct DispatchDirectedNumericOperations
    : DefineDirectedGroupOperators<X, NX> {
    typedef Operations<X> OperationsType;
    /** @brief Return the additive zero with the operand characteristics. */
    friend X nul(const X& x) { return OperationsType::_nul(x); }
    /** @brief Return the unary positive value. */
    friend X pos(const X& x) { return OperationsType::_pos(x); }
    /** @brief Return the additive inverse. */
    friend NX neg(const X& x) { return OperationsType::_neg(x); }
    /** @brief Divide the value by two. */
    friend X hlf(const X& x) { return OperationsType::_hlf(x); }
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const X& x2) { return OperationsType::_add(x1, x2); }
    /** @brief Return the difference of the operands. */
    friend R sub(const X& x1, NX const& x2) { return OperationsType::_sub(x1, x2); }
    /** @brief Return the square root. */
    friend R sqrt(const X& x) { return OperationsType::_sqrt(x); }
    /** @brief Return the exponential. */
    friend R exp(const X& x) { return OperationsType::_exp(x); }
    /** @brief Return the natural logarithm. */
    friend R log(const X& x) { return OperationsType::_log(x); }
    /** @brief Return the inverse tangent. */
    friend R atan(const X& x) { return OperationsType::_atan(x); }
    /** @brief Return the maximum of the operands. */
    friend X max(const X& x1, const X& x2) { return OperationsType::_max(x1, x2); }
    /** @brief Return the minimum of the operands. */
    friend X min(const X& x1, const X& x2) { return OperationsType::_min(x1, x2); }
};

/** @brief Dispatch positive directed numeric operations for the supported operand types. */
template<class X, class QX, class R = X, class QR = QX>
struct DispatchPositiveDirectedNumericOperations
    : ProvideDirectedSemiFieldOperators<X, QX, R, QR> {
    typedef Operations<X> OperationsType;
    /** @brief Return the additive zero with the operand characteristics. */
    friend X nul(const X& x) { return OperationsType::_nul(x); }
    /** @brief Divide the value by two. */
    friend X hlf(const X& x) { return OperationsType::_hlf(x); }
    /** @brief Return the square. */
    friend R sqr(const X& x) { return OperationsType::_sqr(x); }
    /** @brief Return the multiplicative inverse. */
    friend QR rec(const X& x) { return OperationsType::_rec(x); }
    /** @brief Return the sum of the operands. */
    friend R add(const X& x1, const X& x2) { return OperationsType::_add(x1, x2); }
    /** @brief Return the product of the operands. */
    friend R mul(const X& x1, const X& x2) { return OperationsType::_mul(x1, x2); }
    /** @brief Return the quotient of the operands. */
    friend R div(const X& x1, QX const& x2) { return OperationsType::_div(x1, x2); }
    /** @brief Raise the value to the specified power. */
    friend R pow(const X& x, Nat m) { return OperationsType::_pow(x, m); }
    /** @brief Return the square root. */
    friend R sqrt(const X& x) { return OperationsType::_sqrt(x); }
    /** @brief Return the inverse tangent. */
    friend R atan(const X& x) { return OperationsType::_atan(x); }
    /** @brief Return the maximum of the operands. */
    friend X max(const X& x1, const X& x2) { return OperationsType::_max(x1, x2); }
    /** @brief Return the minimum of the operands. */
    friend X min(const X& x1, const X& x2) { return OperationsType::_min(x1, x2); }
    /** @brief Return the absolute value. */
    friend X abs(const X& x) { return OperationsType::_abs(x); }
};

/** @brief Dispatch positive numeric operations for the supported operand types. */
template<class X, class R = X>
struct DispatchPositiveNumericOperations
    : DispatchPositiveDirectedNumericOperations<X, X, R, R> {
};

/** @brief Dispatch comparison operations for the supported operand types. */
template<class X, class LT, class EQ = LT>
struct DispatchComparisonOperations
    : DefineComparisonOperators<X, LT, EQ> {
    typedef Operations<X> OperationsType;
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, const X& x2) { return OperationsType::_eq(x1, x2); }
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, const X& x2) { return OperationsType::_lt(x1, x2); }
};

/** @brief Dispatch directed comparison operations for the supported operand types. */
template<class X, class NX, class LT, class EQ>
struct DispatchDirectedComparisonOperations
    : DefineDirectedComparisonOperators<X, NX, LT, EQ> {
    typedef Operations<X> OperationsType;
    /** @brief Test equality of the operands. */
    friend EQ eq(const X& x1, NX const& nx2) { return OperationsType::_eq(x1, nx2); }
    /** @brief Test whether the first operand is less than the second. */
    friend LT lt(const X& x1, NX const& nx2) { return OperationsType::_lt(x1, nx2); }
};

} // namespace Ariadne

#endif
