/***************************************************************************
 *            numeric/floatdp.hpp
 *
 *  Copyright  2008-22  Pieter Collins
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

/** @file numeric/floatdp.hpp
 *  @brief Raw floating-point number class based on double-precision floats.
 */

#ifndef ARIADNE_FLOAT64_HPP
#define ARIADNE_FLOAT64_HPP

#include <algorithm> // For std::max, std::min
#include <cmath>     // For std::floor std::ceil etc
#include <iosfwd>    // For std::floor std::ceil etc
#include <limits>    // For std::numeric_limits<double>

#include "numeric/builtin.hpp"
#include "numeric/declarations.hpp"
#include "numeric/double.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/float_operations.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/operators.hpp"
#include "numeric/rounding.hpp"
#include "numeric/sign.hpp"

namespace Ariadne {

/** @brief Tag requesting an uninitialized numeric object. */
struct NoInit {};

class Rational;
enum class Comparison : ComparableEnumerationType;

// FIXME: Returns EQUAL if either argument is NaN
/** @brief Compare values through their generic numeric interfaces. */
template<class X1, class X2>
Comparison _generic_cmp(X1 x1, X2 x2) {
    return x1 > x2 ? Comparison::GREATER : (x1 < x2 ? Comparison::LESS : Comparison::EQUAL);
}
/** @brief Compare values through their generic numeric interfaces. */
template<class OP, class X1, class X2>
Comparison _generic_cmp(OP op, X1 x1, X2 x2) {
    return op(cmp(x1, x2), Comparison::EQUAL);
}

//! \ingroup NumericModule
//! \brief The precision of a FloatDP object. Since this is fixed, the class is only a tag; all objects are equal.
//! \relates FloatDP
class DoublePrecision {
  public:
    //! \brief Default constructor
    constexpr DoublePrecision() {}
    //! \brief <p/>
    /** @brief Return the maximum of the operands. */
    friend constexpr DoublePrecision max(DoublePrecision, DoublePrecision) { return DoublePrecision(); }
    //! \brief <p/>
    /** @brief Return the minimum of the operands. */
    friend constexpr DoublePrecision min(DoublePrecision, DoublePrecision) { return DoublePrecision(); }
    //! \brief <p/>
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend constexpr Bool operator<=(DoublePrecision, DoublePrecision) { return true; }
    //! \brief Compare the two precisions for equality. Always \c true.
    friend constexpr Bool operator==(DoublePrecision, DoublePrecision) { return true; }
    //! \brief Write to an output stream. <p/>
    friend OutputStream& operator<<(OutputStream& os, DoublePrecision) { return os << "dp"; }
    //! \brief Write a full representation to an output stream. <p/>
    friend OutputStream& repr(OutputStream& os, DoublePrecision) { return os << "DoublePrecision()"; }
};
//! \brief Shorthand for DoublePrecision.
//! \relates DoublePrecision
using DP = DoublePrecision;
static const DoublePrecision double_precision = DoublePrecision();
static const DoublePrecision dp = DP();

template<class FLT>
class Rounded;
using RoundedFloatDP = Rounded<FloatDP>;

//! \brief Shorthand for \ref Float<DP>.
typedef Float<DP> FloatDP;

//! \ingroup NumericModule
//! \brief Floating point numbers (double precision) using rounded arithmetic. Abbreviation \ref Ariadne::FloatDP "Flobber::FloatDP".
//! \details
//! The \c %FloatDP class represents floating-point numbers.
//! Unless otherwise mentioned, operations on floating-point numbers are performed approximately, with no guarantees
//! on the output.
//!
//! To implement <em>rounded arithmetic</em>, arithmetical operations of \c %FloatDP can be performed with guaranteed rounding by
//! specifying \c _up and \c _down suffixes to arithmetical functions \c add, \c sub, \c mul and \c div.
//! Additionally, operations can be performed in the current <em>rounding mode</em> by using the \c _rnd suffix,
//! or with rounding reversed using the \c _opp suffix.
//! The \c _approx suffix is provided to specifically indicate that the operation is computed approximately.
//!
//! %Ariadne floating-point numbers can be constructed by conversion from built-in C++ types.
//! Note that the value of a built-in floating-point value may differ from the mathematical value of the literal.
//! For example, while <c>%FloatDP(3.25)</c> is represented exactly, <c>%FloatDP(3.3)</c> has a value of \f$3.2999999999999998224\ldots\f$.
//! \note In the future, the construction of a \c %FloatDP from a string literal may be supported.
//! \sa Float<MP>, Real, FloatBall, FloatBounds, FloatUpperBound, FloatLowerBound, FloatApproximation, FloatError
template<>
class Float<DP>
    : public DeclareRoundedOperations<Float<DP>> {
  public:
    volatile double dbl;

  public:
    typedef ExactTag Paradigm;
    typedef FloatDP NumericType;
    //! \brief The abstract number class providing the same kind of information.
    typedef ExactNumber GenericType;
    //! \brief The type representing the precision of the number.
    typedef DoublePrecision PrecisionType;
    //! \brief The type representing the rounding modes available for arithmetic.
    typedef BuiltinRoundingModeType RoundingModeType;

  public:
    static const RoundingModeType ROUND_TO_NEAREST = Ariadne::ROUND_TO_NEAREST;
    static const RoundingModeType ROUND_DOWNWARD = Ariadne::ROUND_DOWNWARD;
    static const RoundingModeType ROUND_UPWARD = Ariadne::ROUND_UPWARD;
    static const RoundingModeType ROUND_TOWARD_ZERO = Ariadne::ROUND_TOWARD_ZERO;

    /** @brief Return the rounding mode. */
    static RoundingModeType get_rounding_mode(); //!< <p/>
    /** @brief Set the rounding mode. */
    static Void set_rounding_mode(RoundingModeType); //!< <p/>
    /** @brief Set the rounding downward. */
    static Void set_rounding_downward(); //!< <p/>
    /** @brief Set the rounding upward. */
    static Void set_rounding_upward(); //!< <p/>
    /** @brief Set the rounding to nearest. */
    static Void set_rounding_to_nearest(); //!< <p/>
    /** @brief Set the rounding toward zero. */
    static Void set_rounding_toward_zero(); //!< <p/>

  public:
    /** @brief Return the default precision. */
    static DoublePrecision get_default_precision(); //!< <p/>
    /** @brief Return the working precision. */
    DoublePrecision precision() const; //!< <p/>
    /** @brief Set the precision. */
    Void set_precision(DoublePrecision); //!< <p/>

  public:
    /** @brief Return a not-a-number value. */
    static FloatDP nan(DoublePrecision pr); //!< <p/>
    /** @brief Return an infinite value. */
    static FloatDP inf(Sign sgn, DoublePrecision pr); //!< <p/>
    /** @brief Return an infinite value. */
    static FloatDP inf(DoublePrecision pr); //!< <p/>

    /** @brief Return the maximum of the operands. */
    static FloatDP max(DoublePrecision pr); //!< <p/>
    /** @brief Return the precision-dependent unit roundoff or tolerance. */
    static FloatDP eps(DoublePrecision pr); //!< <p/>
    /** @brief Return the minimum of the operands. */
    static FloatDP min(DoublePrecision pr); //!< <p/>

  private:
    //! \brief Default constructor creates an uninitialised number.
    Float() : dbl() {}
    //! \brief Convert from a built-in double-precision floating-point number.
    explicit Float(double x) : dbl(x) {}

  public:
    /** @brief Construct a \c Float with the specified characteristics. */
    explicit Float(NoInit) : dbl() {}
    /** @brief Construct a \c Float with the specified characteristics. */
    explicit Float(DoublePrecision) : dbl() {}
    //! \brief Convert from a built-in double-precision floating-point number.
    //    explicit Float(ExactDouble const& x);
    //    explicit Float(double x, DoublePrecision) : dbl(x) { }
    template<BuiltinIntegral N>
    Float(N n, DoublePrecision) : dbl(static_cast<double>(n)) {}
    /** @brief Construct a \c Float from x. */
    Float(const ExactDouble& x, DoublePrecision);
    /** @brief Construct a \c Float from x. */
    Float(const TwoExp& x, DoublePrecision);
    /** @brief Construct a \c Float from w. */
    Float(const Dyadic& w, DoublePrecision);
    //! \brief Copy constructor.
    Float(const FloatDP& x) : dbl(x.dbl) {}
    /** @brief Construct a \c Float from x. */
    Float(const FloatDP& x, DoublePrecision) : dbl(x.dbl) {}
    /** @brief Assign the supplied value and return this object. */
    template<BuiltinIntegral N>
    FloatDP& operator=(N n) {
        return this->operator=(ExactDouble(n));
    }
    //! \brief Construct from an exact string literal.
    explicit Float(const String& str, DoublePrecision);
    /** @brief Assign the supplied value and return this object. */
    FloatDP& operator=(const ExactDouble& x);
    //! \brief Assignment from a dyadic number \a w. Requires \a w to be exactly representable as a double-precision floating-point number.
    FloatDP& operator=(const Dyadic& w);
    //! \brief Copy assignment.
    FloatDP& operator=(const FloatDP& x) {
        this->dbl = x.dbl;
        return *this;
    }

    //! \brief Construct from another FloatDP using given rounding
    explicit Float(const FloatDP& d, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from an integer using given rounding
    template<BuiltinIntegral N>
    explicit Float(N n, RoundingModeType rnd, PrecisionType pr)
        : Float(ExactDouble(n), rnd, pr) {}
    //! \brief Construct from an integer using given rounding
    explicit Float(ExactDouble x, RoundingModeType, PrecisionType)
        : Float(x.get_d()) {}
    //! \brief Construct from an integer number using given rounding
    explicit Float(const Integer&, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from a dyadic number with given rounding
    explicit Float(const Dyadic& w, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from a decimal number with given rounding
    explicit Float(const Decimal& d, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from a FloatMP using given rounding
    explicit Float(const FloatMP& d, RoundingModeType rnd, PrecisionType pr);
    //! \brief Construct from a rational number with given rounding
    explicit Float(const Rational& q, RoundingModeType rnd, PrecisionType pr);
    //! \brief Convert to a dyadic number.
    explicit operator Dyadic() const;
    //! \brief Convert to a rational number.
    explicit operator Rational() const;

    //! \brief Convert to a generic number.
    operator ExactNumber() const;

    //! \brief Construct a ball with radius \a e.
    Ball<FloatDP> pm(const Error<FloatDP>& e) const;

  public:
    /** @brief Return the underlying numeric representation. */
    const FloatDP& raw() const { return *this; }
    /** @brief Return the underlying numeric representation. */
    FloatDP& raw() { return *this; }
    //! \brief An approximation by a built-in double-precision floating-point number.
    double get_d() const { return this->dbl; }
    //! \brief The exact value as a decimal string.
    String literal() const;
    //! \brief An approximate value as a decimal string, rounded by \a rnd.
    //! The exact value is guaranteed to be the nearest value to the value's precision.
    String literal(RoundingModeType rnd) const;

  public:
    //! \name Special value tests
    //!@{
    /** @brief Test whether the value is not a number. */
    friend Bool is_nan(FloatDP x) { return std::isnan(x.dbl); } //!< <p/>
    /** @brief Test whether the value is infinite. */
    friend Bool is_inf(FloatDP x) { return std::isinf(x.dbl); } //!< <p/>
    /** @brief Test whether the value is finite. */
    friend Bool is_finite(FloatDP x) { return std::isfinite(x.dbl); } //!< <p/>
    /** @brief Test whether the value is zero. */
    friend Bool is_zero(FloatDP x) { return x.dbl == 0.0; } //!< <p/>
    //!@}

    //! \name Rounding operations
    //!@{
    friend FloatDP next(RoundUpward rnd, FloatDP x) { return add(rnd, x, FloatDP::min(x.precision())); }   //!< Round \a x up by 1 ulp.
    friend FloatDP next(RoundDownward rnd, FloatDP x) { return sub(rnd, x, FloatDP::min(x.precision())); } //!< Round \a x up by 1 ulp.

    friend FloatDP floor(FloatDP x) { return FloatDP(std::floor(x.dbl)); } //!< Round \a x down to the nearest lower integer.
    friend FloatDP ceil(FloatDP x) { return FloatDP(std::ceil(x.dbl)); }   //!< Round \a x to the nearest integer. %Rounding of halves is implementation-dependent.
    friend FloatDP round(FloatDP x) { return FloatDP(std::round(x.dbl)); } //!< Round \a x up to the nearest higher integer.

    friend Integer cast_integer(FloatDP x); //!< Round \a x to the nearest integer
                                            //!@}

  public:
    //! \name Arithmetic operators
    //!@{
    /** @brief Return the sum or unary positive value. */
    friend FloatDP operator+(FloatDP x) {
        volatile double xv = x.dbl;
        return FloatDP(+xv);
    } //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend FloatDP operator-(FloatDP x) {
        volatile double xv = x.dbl;
        return FloatDP(-xv);
    } //!< <p/>
    /** @brief Return the sum or unary positive value. */
    friend Bounds<FloatDP> operator+(const FloatDP& x1, const FloatDP& x2); //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend Bounds<FloatDP> operator-(const FloatDP& x1, const FloatDP& x2); //!< <p/>
    /** @brief Return the product. */
    friend Bounds<FloatDP> operator*(const FloatDP& x1, const FloatDP& x2); //!< <p/>
    /** @brief Return the quotient. */
    friend Bounds<FloatDP> operator/(const FloatDP& x1, const FloatDP& x2); //!< <p/>
    /** @brief Return the product. */
    friend FloatDP operator*(TwoExp e1, FloatDP x2) { return shft(x2, e1.exponent()); } //!< <p/>
    /** @brief Return the product. */
    friend FloatDP operator*(FloatDP x1, TwoExp e2) { return shft(x1, e2.exponent()); } //!< <p/>
    /** @brief Return the quotient. */
    friend FloatDP operator/(FloatDP x1, TwoExp e2) { return shft(x1, -e2.exponent()); } //!< <p/>
    /** @brief Multiply by the operand in place. */
    friend FloatDP& operator*=(FloatDP& v1, TwoExp e2) { return v1 = v1 * e2; } //!< <p/>
    /** @brief Divide by the operand in place. */
    friend FloatDP& operator/=(FloatDP& v1, TwoExp e2) { return v1 = v1 / e2; } //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operators
    /** @brief Test equality of the operands. */
    friend Boolean operator==(FloatDP x1, FloatDP x2) { return _cmp(Eq(), x1, x2); } //!< <p/>
    /** @brief Test inequality of the operands. */
    friend Boolean operator!=(FloatDP x1, FloatDP x2) { return _cmp(Ne(), x1, x2); } //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend Boolean operator<(FloatDP x1, FloatDP x2) { return _cmp(Lt(), x1, x2); } //!< <p/>
    /** @brief Test whether the first operand is greater than the second. */
    friend Boolean operator>(FloatDP x1, FloatDP x2) { return _cmp(Gt(), x1, x2); } //!< <p/>
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend Boolean operator<=(FloatDP x1, FloatDP x2) { return _cmp(Le(), x1, x2); } //!< <p/>
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend Boolean operator>=(FloatDP x1, FloatDP x2) { return _cmp(Ge(), x1, x2); } //!< <p/>
    //!@}

    //! \name Exact arithmetical functions
    //!@{
    friend FloatDP nul(FloatDP) { return FloatDP(+0.0); } //!< Zero \a 0.
    friend FloatDP pos(FloatDP x) {
        volatile double xv = x.dbl;
        return FloatDP(+xv);
    } //!< Identity \a +x.
    friend FloatDP neg(FloatDP x) {
        volatile double xv = x.dbl;
        return FloatDP(-xv);
    } //!< Negation \a -x.
    friend FloatDP hlf(FloatDP x) {
        volatile double xv = x.dbl;
        return FloatDP(xv / 2);
    } //!< Half \a x÷2.
    friend FloatDP shft(FloatDP x, Int n) {
        volatile double xv = x.dbl;
        return FloatDP(ldexp(xv, n));
    } //!< Bit-shift \a x×2<sup>n</sup>.
    friend FloatDP mul(FloatDP x1, TwoExp e2) { return shft(x1, e2.exponent()); }  //!< Multiplication by a power of two.
    friend FloatDP div(FloatDP x1, TwoExp e2) { return shft(x1, -e2.exponent()); } //!< Division by a power of two.
    //!@}

    //! \name Outward-rounded arithmetical functions
    //!@{
    friend Bounds<FloatDP> sqr(const FloatDP& x);                                        //!< Square \a x<sup>2</sup>.
    friend Bounds<FloatDP> rec(const FloatDP& x);                                        //!< Reciprocal \a 1/x.
    friend Bounds<FloatDP> add(const FloatDP& x1, const FloatDP& x2);                    //!< \brief Add \a x1+x2.
    friend Bounds<FloatDP> sub(const FloatDP& x1, const FloatDP& x2);                    //!< \brief Subtract \a x1-x2.
    friend Bounds<FloatDP> mul(const FloatDP& x1, const FloatDP& x2);                    //!< \brief Multiply \a x1×x2.
    friend Bounds<FloatDP> div(const FloatDP& x1, const FloatDP& x2);                    //!< \brief Divide \a x1÷x2.
    friend Bounds<FloatDP> fma(const FloatDP& x1, const FloatDP& x2, const FloatDP& x3); //!< \brief Fused multiply-and-add \a x1×x2+x3.
    friend Bounds<FloatDP> pow(const FloatDP& x, Nat m);                                 //!< \brief Power \a x<sup>m</sup>.
    friend Bounds<FloatDP> pow(const FloatDP& x, Int n);                                 //!< \brief Power \a x<sup>n</sup>.
    //!@}

    //! \name Algebraic and transcendental functions
    //!@{
    friend Bounds<FloatDP> sqrt(const FloatDP& x); //!< The square root of \a x, √\a x. Requires \a x ≥ 0.
    friend Bounds<FloatDP> exp(const FloatDP& x);  //!< The natural exponent of \a x, \em e<sup>x</sup>.
    friend Bounds<FloatDP> log(const FloatDP& x);  //!< The natural logarithm of \a x, log<em><sub>e</sub> x</em> or ln <em>x</em>. Requires \a x ≥ 0.
    friend Bounds<FloatDP> sin(const FloatDP& x);  //!< The sine of \a x.
    friend Bounds<FloatDP> cos(const FloatDP& x);  //!< The cosine of \a x.
    friend Bounds<FloatDP> tan(const FloatDP& x);  //!< The tangent of \a x, sin(\a x)/cos(\a x).
    friend Bounds<FloatDP> asin(const FloatDP& x); //!< The arc-sine of \a x.
    friend Bounds<FloatDP> acos(const FloatDP& x); //!< The arc-cosine of \a x.
    friend Bounds<FloatDP> atan(const FloatDP& x); //!< The arc-tangent of \a x.
    //!@}

    //! \name Lattice operations
    //!@{
    friend Positive<FloatDP> abs(FloatDP x);                                                 //!< Absolute value \a |x|.
    friend FloatDP max(FloatDP x1, FloatDP x2) { return FloatDP(std::max(x1.dbl, x2.dbl)); } //!< Maximum \a x1∨x2.
    friend FloatDP min(FloatDP x1, FloatDP x2) { return FloatDP(std::min(x1.dbl, x2.dbl)); } //!< Mimimum \a x1∧x2.

    friend PositiveUpperBound<FloatDP> mag(FloatDP x); //!< An upper bound for the absolute value \a |x|.
    friend PositiveLowerBound<FloatDP> mig(FloatDP x); //!< A lower bound for the absolute value \a |x|.
    //!@}

    //! \name Comparison operations
    //!@{
    /** @brief Test equality of the operands. */
    friend Boolean eq(const FloatDP& x1, const FloatDP& x2) { return x1.dbl == x2.dbl; } //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend Boolean lt(const FloatDP& x1, const FloatDP& x2) { return x1.dbl < x2.dbl; } //!< <p/>

    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(FloatDP x1, FloatDP x2) { return _generic_cmp(x1.dbl, x2.dbl); } //!< <p/>
    /** @brief Compare the operands using the comparison result type. */
    template<BuiltinIntegral N>
    friend Comparison cmp(const FloatDP& x1, N n2) { return _generic_cmp(x1.dbl, n2); }
    /** @brief Compare the operands using the comparison result type. */
    template<BuiltinIntegral N>
    friend Comparison cmp(N n1, const FloatDP& x2) { return _generic_cmp(n1, x2.dbl); }
#ifdef DOXYGEN
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const FloatDP& x1, Nat m2); //!< <p/>
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(Nat m1, const FloatDP& x2); //!< <p/>
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const FloatDP& x1, Int n2); //!< <p/>
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(Int n1, const FloatDP& x2); //!< <p/>
#endif                                                // DOXYGEN
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(FloatDP x1, ExactDouble x2) { return _generic_cmp(x1.dbl, x2.get_d()); } //!< <p/>
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(ExactDouble x1, FloatDP x2) { return _generic_cmp(x1.get_d(), x2.dbl); } //!< <p/>
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(FloatDP x1, const Rational& x2); //!< <p/>
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const Rational& x1, FloatDP x2); //!< <p/>
                                                           //!@}

  public:
    //! \name Correctly rounded operations
    //!@{
    /** @brief Return the square. */
    friend FloatDP sqr(CurrentRoundingMode, FloatDP x) { return FloatDP(sqr_rnd(x.dbl)); } //!< <p/>
    /** @brief Return the multiplicative inverse. */
    friend FloatDP rec(CurrentRoundingMode, FloatDP x) { return FloatDP(rec_rnd(x.dbl)); } //!< <p/>
    /** @brief Return the sum of the operands. */
    friend FloatDP add(CurrentRoundingMode, FloatDP x1, FloatDP x2) { return FloatDP(add_rnd(x1.dbl, x2.dbl)); } //!< <p/>
    /** @brief Return the difference of the operands. */
    friend FloatDP sub(CurrentRoundingMode, FloatDP x1, FloatDP x2) { return FloatDP(sub_rnd(x1.dbl, x2.dbl)); } //!< <p/>
    /** @brief Return the product of the operands. */
    friend FloatDP mul(CurrentRoundingMode, FloatDP x1, FloatDP x2) { return FloatDP(mul_rnd(x1.dbl, x2.dbl)); } //!< <p/>
    /** @brief Return the quotient of the operands. */
    friend FloatDP div(CurrentRoundingMode, FloatDP x1, FloatDP x2) { return FloatDP(div_rnd(x1.dbl, x2.dbl)); } //!< <p/>
    /** @brief Compute the product of the first two operands plus the third. */
    friend FloatDP fma(CurrentRoundingMode, FloatDP x1, FloatDP x2, FloatDP x3) { return FloatDP(fma_rnd(x1.dbl, x2.dbl, x3.dbl)); } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend FloatDP pow(CurrentRoundingMode, FloatDP x, Nat m) { return FloatDP(pow_rnd(x.dbl, m)); } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend FloatDP pow(CurrentRoundingMode, FloatDP x, Int n) { return FloatDP(pow_rnd(x.dbl, n)); } //!< <p/>
    /** @brief Return the square root. */
    friend FloatDP sqrt(CurrentRoundingMode, FloatDP x) { return FloatDP(sqrt_rnd(x.dbl)); } //!< <p/>
    /** @brief Return the exponential. */
    friend FloatDP exp(CurrentRoundingMode, FloatDP x) { return exp(FloatDP::get_rounding_mode(), x); }
    /** @brief Return the natural logarithm. */
    friend FloatDP log(CurrentRoundingMode, FloatDP x) { return log(FloatDP::get_rounding_mode(), x); }
    /** @brief Return the sine. */
    friend FloatDP sin(CurrentRoundingMode, FloatDP x) { return sin(FloatDP::get_rounding_mode(), x); }
    /** @brief Return the cosine. */
    friend FloatDP cos(CurrentRoundingMode, FloatDP x) { return cos(FloatDP::get_rounding_mode(), x); }
    /** @brief Return the tangent. */
    friend FloatDP tan(CurrentRoundingMode, FloatDP x) { return tan(FloatDP::get_rounding_mode(), x); }
    //! \brief Correctly rounded logistic sigmoid (1/(1+exp(-x))).
    friend FloatDP sigmoid(CurrentRoundingMode, FloatDP x) { return sigmoid(FloatDP::get_rounding_mode(), x); }
    //! \brief Correctly rounded softplus (log(1+exp(x))).
    friend FloatDP softplus(CurrentRoundingMode, FloatDP x) { return softplus(FloatDP::get_rounding_mode(), x); }
    /** @brief Return the hyperbolic tangent. */
    friend FloatDP tanh(CurrentRoundingMode, FloatDP x) { return tanh(FloatDP::get_rounding_mode(), x); }
    /** @brief Return the inverse sine. */
    friend FloatDP asin(CurrentRoundingMode, FloatDP x) { return asin(FloatDP::get_rounding_mode(), x); }
    /** @brief Return the inverse cosine. */
    friend FloatDP acos(CurrentRoundingMode, FloatDP x) { return acos(FloatDP::get_rounding_mode(), x); }
    /** @brief Return the inverse tangent. */
    friend FloatDP atan(CurrentRoundingMode, FloatDP x) { return atan(FloatDP::get_rounding_mode(), x); }
    /** @brief Return pi at the requested precision. */
    static FloatDP pi(CurrentRoundingMode, PrecisionType) { return FloatDP(pi_rnd()); } //!< <p/>
                                                                                        //!@}

  private:
    /** @brief Apply the supplied operation to the operands with the specified rounding mode. */
    template<class OP>
    static FloatDP _apply_rnd(OP op, RoundingModeType rnd, FloatDP x1, FloatDP x2, FloatDP x3) {
        RoundingSession<FloatDP> session;
        session.set(rnd);
        return FloatDP(op(x1.dbl, x2.dbl, x3.dbl));
    }
    /** @brief Apply the supplied operation to the operands with the specified rounding mode. */
    template<class OP>
    static FloatDP _apply_rnd(OP op, RoundingModeType rnd, FloatDP x1, FloatDP x2) {
        RoundingSession<FloatDP> session;
        session.set(rnd);
        return FloatDP(op(x1.dbl, x2.dbl));
    }
    /** @brief Apply the supplied operation to the operands with the specified rounding mode. */
    template<class OP>
    static FloatDP _apply_rnd(OP op, RoundingModeType rnd, FloatDP x) {
        RoundingSession<FloatDP> session;
        session.set(rnd);
        return FloatDP(op(x.dbl));
    }
    /** @brief Apply the supplied operation to the operands with the specified rounding mode. */
    template<class OP, BuiltinIntegral N>
    static FloatDP _apply_rnd(OP op, RoundingModeType rnd, FloatDP x, N n) {
        RoundingSession<FloatDP> session;
        session.set(rnd);
        return FloatDP(op(x.dbl, n));
    }
    // D7: certified nearest-mode residual path for directed IEEE-754 arithmetic.
    /** @brief Return the sum of the operands with a certified rounding bound. */
    static FloatDP _certified_add(RoundingModeType, FloatDP, FloatDP);
    /** @brief Return the difference of the operands with a certified rounding bound. */
    static FloatDP _certified_sub(RoundingModeType, FloatDP, FloatDP);
    /** @brief Return the product of the operands with a certified rounding bound. */
    static FloatDP _certified_mul(RoundingModeType, FloatDP, FloatDP);
    /** @brief Return the quotient of the operands with a certified rounding bound. */
    static FloatDP _certified_div(RoundingModeType, FloatDP, FloatDP);
    /** @brief Compute the product of the first two operands plus the third with a certified rounding bound. */
    static FloatDP _certified_fma(RoundingModeType, FloatDP, FloatDP, FloatDP);
    /** @brief Return the square root with a certified rounding bound. */
    static FloatDP _certified_sqrt(RoundingModeType, FloatDP);
    /** @brief Raise the value to the specified power with a certified rounding bound. */
    static FloatDP _certified_pow(RoundingModeType, FloatDP, Nat);

  public:
    //! \name Explicitly rounded operations
    //!@{
    /** @brief Return the maximum of the operands. */
    friend FloatDP max(RoundingModeType, FloatDP x1, FloatDP x2) { return FloatDP(std::max(x1.dbl, x2.dbl)); } //!< <p/>
    /** @brief Return the minimum of the operands. */
    friend FloatDP min(RoundingModeType, FloatDP x1, FloatDP x2) { return FloatDP(std::min(x1.dbl, x2.dbl)); } //!< <p/>
    /** @brief Return the absolute value. */
    friend FloatDP abs(RoundingModeType, FloatDP x) { return FloatDP(std::fabs(x.dbl)); } //!< <p/>
    /** @brief Return an upper bound on the absolute value. */
    friend FloatDP mag(RoundingModeType, FloatDP x) { return FloatDP(std::fabs(x.dbl)); } //!< <p/>

    /** @brief Add the operand in place. */
    friend FloatDP& iadd(RoundingModeType rnd, FloatDP& x1, FloatDP x2); //!< <p/>
    /** @brief Subtract the operand in place. */
    friend FloatDP& isub(RoundingModeType rnd, FloatDP& x1, FloatDP x2); //!< <p/>
    /** @brief Multiply by the operand in place. */
    friend FloatDP& imul(RoundingModeType rnd, FloatDP& x1, FloatDP x2); //!< <p/>
    /** @brief Divide by the operand in place. */
    friend FloatDP& idiv(RoundingModeType rnd, FloatDP& x1, FloatDP x2); //!< <p/>

    /** @brief Return pi at the requested precision. */
    static FloatDP pi(RoundingModeType rnd, PrecisionType pr); //!< <p/>

    /** @brief Return the additive zero with the operand characteristics. */
    friend FloatDP nul(RoundingModeType rnd, FloatDP x) { return _apply_rnd(&nul_rnd, rnd, x); } //!< <p/>
    /** @brief Return the unary positive value. */
    friend FloatDP pos(RoundingModeType rnd, FloatDP x) { return _apply_rnd(&pos_rnd, rnd, x); } //!< <p/>
    /** @brief Return the additive inverse. */
    friend FloatDP neg(RoundingModeType rnd, FloatDP x) { return _apply_rnd(&neg_rnd, rnd, x); } //!< <p/>
    /** @brief Divide the value by two. */
    friend FloatDP hlf(RoundingModeType rnd, FloatDP x) { return _apply_rnd(&hlf_rnd, rnd, x); } //!< <p/>
    /** @brief Return the sum of the operands. */
    friend FloatDP add(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return (rnd == ROUND_DOWNWARD || rnd == ROUND_UPWARD) ? _certified_add(rnd, x1, x2) : _apply_rnd(&add_rnd, rnd, x1, x2); } //!< <p/>
    /** @brief Return the difference of the operands. */
    friend FloatDP sub(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return (rnd == ROUND_DOWNWARD || rnd == ROUND_UPWARD) ? _certified_sub(rnd, x1, x2) : _apply_rnd(&sub_rnd, rnd, x1, x2); } //!< <p/>
    /** @brief Return the product of the operands. */
    friend FloatDP mul(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return (rnd == ROUND_DOWNWARD || rnd == ROUND_UPWARD) ? _certified_mul(rnd, x1, x2) : _apply_rnd(&mul_rnd, rnd, x1, x2); } //!< <p/>
    /** @brief Return the quotient of the operands. */
    friend FloatDP div(RoundingModeType rnd, FloatDP x1, FloatDP x2) { return (rnd == ROUND_DOWNWARD || rnd == ROUND_UPWARD) ? _certified_div(rnd, x1, x2) : _apply_rnd(&div_rnd, rnd, x1, x2); } //!< <p/>
    /** @brief Compute the product of the first two operands plus the third. */
    friend FloatDP fma(RoundingModeType rnd, FloatDP x1, FloatDP x2, FloatDP x3) { return (rnd == ROUND_DOWNWARD || rnd == ROUND_UPWARD) ? _certified_fma(rnd, x1, x2, x3) : _apply_rnd(&fma_rnd, rnd, x1, x2, x3); } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend FloatDP pow(RoundingModeType rnd, FloatDP x, Nat m) { return (rnd == ROUND_DOWNWARD || rnd == ROUND_UPWARD) ? _certified_pow(rnd, x, m) : _apply_rnd((double (*)(double, unsigned int))&pow_rnd, rnd, x, m); } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend FloatDP pow(RoundingModeType rnd, FloatDP x, Int n) { return _apply_rnd((double (*)(double, int))&pow_rnd, rnd, x, n); } //!< <p/>
    /** @brief Return the square. */
    friend FloatDP sqr(RoundingModeType rnd, FloatDP x) { return _apply_rnd(&sqr_rnd, rnd, x); } //!< <p/>
    /** @brief Return the multiplicative inverse. */
    friend FloatDP rec(RoundingModeType rnd, FloatDP x) { return _apply_rnd(&rec_rnd, rnd, x); } //!< <p/>
    /** @brief Return the square root. */
    friend FloatDP sqrt(RoundingModeType rnd, FloatDP x) { return (rnd == ROUND_DOWNWARD || rnd == ROUND_UPWARD) ? _certified_sqrt(rnd, x) : _apply_rnd(&sqrt_rnd, rnd, x); } //!< <p/>
    /** @brief Return the exponential. */
    friend FloatDP exp(RoundingModeType rnd, FloatDP x);
    /** @brief Return the natural logarithm. */
    friend FloatDP log(RoundingModeType rnd, FloatDP x);
    /** @brief Return the sine. */
    friend FloatDP sin(RoundingModeType rnd, FloatDP x);
    /** @brief Return the cosine. */
    friend FloatDP cos(RoundingModeType rnd, FloatDP x);
    /** @brief Return the tangent. */
    friend FloatDP tan(RoundingModeType rnd, FloatDP x);
    //! \brief Correctly rounded logistic sigmoid in the given direction.
    friend FloatDP sigmoid(RoundingModeType rnd, FloatDP x);
    //! \brief Correctly rounded softplus in the given direction.
    friend FloatDP softplus(RoundingModeType rnd, FloatDP x);
    /** @brief Return the hyperbolic tangent. */
    friend FloatDP tanh(RoundingModeType rnd, FloatDP x);
    /** @brief Return the inverse sine. */
    friend FloatDP asin(RoundingModeType rnd, FloatDP x);
    /** @brief Return the inverse cosine. */
    friend FloatDP acos(RoundingModeType rnd, FloatDP x);
    /** @brief Return the inverse tangent. */
    friend FloatDP atan(RoundingModeType rnd, FloatDP x);

    /** @brief Return the midpoint value. */
    friend FloatDP med(RoundingModeType rnd, FloatDP x1, FloatDP x2) {
        rounding_mode_t rounding_mode = get_rounding_mode();
        set_rounding_mode(rnd); //!< <p/>
        FloatDP r = hlf(add(rnd, x1, x2));
        set_rounding_mode(rounding_mode);
        return r;
    } //!< <p/>
    /** @brief Return the radius. */
    friend FloatDP rad(RoundingModeType rnd, FloatDP x1, FloatDP x2) {
        rounding_mode_t rounding_mode = get_rounding_mode();
        set_rounding_mode(rnd); //!< <p/>
        FloatDP r = hlf(sub(rnd, x2, x1));
        set_rounding_mode(rounding_mode);
        return r;
    } //!< <p/>
      //!@}

  private:
    /** @brief Compare the operands using the comparison result type. */
    template<class OP>
    static Boolean _cmp(OP op, FloatDP x1, FloatDP x2) {
        return op(x1.dbl, x2.dbl);
    }
    /*
        template<class OP, AnExactNumber Y> static Boolean _cmp(OP op, FloatDP const&  x1, Y const& y2) {
            if constexpr (BuiltinIntegral<Y>) { return op(x1.dbl,y2); }
            else if constexpr (Same<Y,ExactDouble>) { return op(x1.dbl,y2.get_d()); }
            else { return op(cmp(x1,y2),Comparison::EQUAL); }
        }
        template<class OP, AnExactNumber Y> static Boolean _cmp(OP op, Y const& y1, FloatDP const&  x2) {
            if constexpr (BuiltinIntegral<Y>) { return op(y1,x2.dbl); }
            else if constexpr (Same<Y,ExactDouble>) { return op(y1.get_d(),x2.dbl); }
            else { return op(cmp(y1,x2),Comparison::EQUAL); }
        }
    */

  public:
    //!@{
    //! \name Exact information tests
    friend Bool same(const FloatDP& x1, const FloatDP& x2) { return x1 == x2; } //!< Test equality of \a x1 and x2.
                                                                                //!@}

  public:
    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const FloatDP& x); //!< Write to an output stream.
    friend InputStream& operator>>(InputStream& os, FloatDP& x);         //!< Read from an input stream.
                                                                         //!@}

  private:
    /** @brief Write the value to the output stream. */
    friend OutputStream& write(OutputStream& os, const FloatDP& x, DecimalPlaces dgts, RoundingModeType rnd);
    /** @brief Write the value to the output stream. */
    friend OutputStream& write(OutputStream& os, const FloatDP& x, DecimalPrecision dgts, RoundingModeType rnd);
    /** @brief Write the representation of the object. */
    friend OutputStream& repr(OutputStream& os, const FloatDP& x);
    /** @brief Write the representation of the object. */
    friend OutputStream& repr(OutputStream& os, const FloatDP& x, RoundingModeType rnd);

  private:
    // Rounded arithmetic
    /** @brief Raise the value to the specified power with the specified rounding mode. */
    friend FloatDP pow_rnd(FloatDP x, Nat m);
    /** @brief Raise the value to the specified power with the specified rounding mode. */
    friend FloatDP pow_rnd(FloatDP x, Int n);
    /** @brief Return the square root with the specified rounding mode. */
    friend FloatDP sqrt_rnd(FloatDP x);
    /** @brief Return the exponential with the specified rounding mode. */
    friend FloatDP exp_rnd(FloatDP x);
    /** @brief Return the natural logarithm with the specified rounding mode. */
    friend FloatDP log_rnd(FloatDP x);
    /** @brief Return the sine with the specified rounding mode. */
    friend FloatDP sin_rnd(FloatDP x);
    /** @brief Return the cosine with the specified rounding mode. */
    friend FloatDP cos_rnd(FloatDP x);
    /** @brief Return the tangent with the specified rounding mode. */
    friend FloatDP tan_rnd(FloatDP x);
    /** @brief Return the hyperbolic tangent with the specified rounding mode. */
    friend FloatDP tanh_rnd(FloatDP x);
    /** @brief Return the inverse tangent with the specified rounding mode. */
    friend FloatDP atan_rnd(FloatDP x);

    // Opposite rounded arithmetic
    /** @brief Return the unary positive value with opposite directed rounding. */
    friend FloatDP pos_opp(FloatDP x) {
        volatile double t = -x.dbl;
        return FloatDP(-t);
    }
    /** @brief Return the additive inverse with opposite directed rounding. */
    friend FloatDP neg_opp(FloatDP x) {
        volatile double t = x.dbl;
        return FloatDP(-t);
    }
    /** @brief Return the square with opposite directed rounding. */
    friend FloatDP sqr_opp(FloatDP x) {
        volatile double t = -x.dbl;
        t = t * x.dbl;
        return FloatDP(-t);
    }
    /** @brief Return the multiplicative inverse with opposite directed rounding. */
    friend FloatDP rec_opp(FloatDP x) {
        volatile double t = -1.0 / (volatile double&)x.dbl;
        return FloatDP(-t);
    }
    /** @brief Return the sum of the operands with opposite directed rounding. */
    friend FloatDP add_opp(FloatDP x, FloatDP y) {
        volatile double t = -x.dbl;
        t = t - y.dbl;
        return FloatDP(-t);
    }
    /** @brief Return the difference of the operands with opposite directed rounding. */
    friend FloatDP sub_opp(FloatDP x, FloatDP y) {
        volatile double t = -x.dbl;
        t = t + y.dbl;
        return FloatDP(-t);
    }
    /** @brief Return the product of the operands with opposite directed rounding. */
    friend FloatDP mul_opp(FloatDP x, FloatDP y) {
        volatile double t = -x.dbl;
        t = t * y.dbl;
        return FloatDP(-t);
    }
    /** @brief Return the quotient of the operands with opposite directed rounding. */
    friend FloatDP div_opp(FloatDP x, FloatDP y) {
        volatile double t = x.dbl;
        t = t / y.dbl;
        return FloatDP(-t);
    }
    /** @brief Raise the value to the specified power with opposite directed rounding. */
    friend FloatDP pow_opp(FloatDP x, int n);

  public:
    static Nat output_places;
    /** @brief Set the output places. */
    static Void set_output_places(Nat pl);

  private:
    friend class Rounded<FloatDP>;
};

template<class R, class A>
R integer_cast(const A& a);
/** @brief Convert the floating-point value to the requested integer type. */
template<>
Nat integer_cast<Nat, FloatDP>(const FloatDP&);
/** @brief Convert the floating-point value to the requested integer type. */
template<>
Int integer_cast<Int, FloatDP>(const FloatDP&);

/** @brief Single-precision floating-point value. */
struct Float32 {
    float flt;

  public:
    /** @brief Construct a \c Float32 from x, rnd. */
    explicit Float32(FloatDP x, BuiltinRoundingModeType rnd) {
        set_builtin_rounding_mode(rnd);
        volatile double value = x.dbl;
        flt = static_cast<float>(value);
    }
    /** @brief Convert to float dp. */
    explicit operator FloatDP() const;
};

} // namespace Ariadne

#endif
