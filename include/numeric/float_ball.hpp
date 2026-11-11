/***************************************************************************
 *            numeric/float_ball.hpp
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

/** @file numeric/float_ball.hpp
 *  @brief Floating-point approximations to real numbers.
 */

#ifndef ARIADNE_FLOAT_BALL_HPP
#define ARIADNE_FLOAT_BALL_HPP

#include "float.decl.hpp"
#include "float_bounds.hpp"
#include "float_factory.hpp"
#include "float_operations.hpp"
#include "float_traits.hpp"
#include "foundation/logical.decl.hpp"
#include "number.decl.hpp"
#include "positive.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

/** @brief Return the precision used for error bounds. */
template<class PRE, class PR>
    requires DefaultConstructible<PRE>
inline PRE _error_precision(PR const&) { return PRE(); }
/** @brief Return the precision used for error bounds. */
template<class PRE, class PR>
    requires(not DefaultConstructible<PRE>) and Constructible<PRE, PR>
inline PRE _error_precision(PR const& pr) { return PRE(pr); }

//! \ingroup NumericModule
//! \brief Floating point approximations to a real number with guaranteed error bounds.
//! \sa Real, FloatDP, Error, Float, Bounds, Approximation.
template<class F, class FE>
class Ball
    : public DefineFloatOperations<Ball<F, FE>>,
      public ProvideConvertedFieldOperations<Bounds<F>, Ball<F, FE>> {
    typedef ValidatedTag P;
    typedef typename F::PrecisionType PR;
    typedef typename FE::PrecisionType PRE;
    static_assert(Constructible<PR, PRE> or DefaultConstructible<PRE>);

  public:
    //! <p/>
    typedef P Paradigm;
    //! <p/>
    typedef Ball<F, FE> NumericType;
    //! <p/>
    typedef Number<P> GenericType;
    //! <p/>
    typedef F RawType;
    //! <p/>
    typedef FE RawErrorType;
    //! <p/>
    typedef PR PrecisionType;
    //! <p/>
    typedef PRE ErrorPrecisionType;
    //! <p/>
    typedef Pair<PR, PRE> CharacteristicsType;

  public:
    //! Construct a ball of radius \f$0\f$ about zero, using precision \a pr for the centre (value),
    //! and using either value \a pr for the precision of the radius (error), or the default value.
    explicit Ball(PrecisionType pr) : _v(0.0_x, pr), _e(0.0_x, _error_precision<PRE>(pr)) {}
    //! Construct a ball of radius \f$0\f$ about zero,
    //! using precision \a pr for the centre (value), and \a pre for the radius (error).
    explicit Ball(PrecisionType pr, ErrorPrecisionType pre) : _v(0.0_x, pr), _e(0.0_x, pre) {}
    /** @brief Construct a \c Ball from prs. */
    explicit Ball(CharacteristicsType prs) : _v(0.0_x, std::get<0>(prs)), _e(0.0_x, std::get<1>(prs)) {}
    /** @brief Construct a \c Ball from v. */
    Ball(const F& v) : _v(v), _e(0.0_x, _error_precision<PRE>(v.precision())) {}
    /** @brief Construct a \c Ball from v, pre. */
    explicit Ball(const F& v, PRE pre) : _v(v), _e(0.0_x, pre) {}
    //! Construct a ball of radius \a e about \a v.
    explicit Ball(const F& v, FE const& e) : _v(v), _e(e) {}
    //! Construct a ball of radius \a error about \a v.
    Ball(const F& v, const Error<FE>& error);
    //! Construct a ball containing the bounds \a x, with error bound of precision \a pre.
    Ball(const Bounds<F>& x, PRE pre);
    /** @brief Disallow this construction to preserve the ownership contract. */
    Ball(const LowerBound<F>& lower, const UpperBound<F>& upper) = delete;

    /** @brief Construct a \c Ball from dv, de, pr. */
    Ball(const Decimal& dv, const Decimal& de, PrecisionType pr);

    /** @brief Construct a \c Ball from d, pr. */
    Ball(const ExactDouble& d, PR pr);
    /** @brief Construct a \c Ball from t, pr. */
    Ball(const TwoExp& t, PR pr);
    /** @brief Construct a \c Ball from z, pr. */
    Ball(const Integer& z, PR pr);
    /** @brief Construct a \c Ball from w, pr. */
    Ball(const Dyadic& w, PR pr);
    /** @brief Construct a \c Ball from d, pr. */
    Ball(const Decimal& d, PR pr);
    /** @brief Construct a \c Ball from q, pr. */
    Ball(const Rational& q, PR pr);
    /** @brief Construct a \c Ball from r, pr. */
    Ball(const Real& r, PR pr);
    /** @brief Construct a \c Ball from x, pr. */
    Ball(const Ball<F, FE>& x, PR pr);
    /** @brief Construct a \c Ball from x, pr, pre. */
    Ball(const Ball<F, FE>& x, PR pr, PRE pre);
    /** @brief Construct a \c Ball from y, pr. */
    Ball(const ValidatedNumber& y, PR pr);

    // FIXME: Constructors for other types
    /** @brief Construct a \c Ball from z, pr, pre. */
    Ball(const Integer& z, PR pr, PRE pre);
    /** @brief Construct a \c Ball from w, pr, pre. */
    Ball(const Dyadic& w, PR pr, PRE pre);
    /** @brief Construct a \c Ball from q, pr, pre. */
    Ball(const Rational& q, PR pr, PRE pre);
    /** @brief Construct a \c Ball from q, pr, pre. */
    Ball(const Real& q, PR pr, PRE pre);
    //! Construct a ball guaranteed to contain the generic validated number \a y,
    //! using precision \a pr for the centre (value), and \a pre for the radius (error).
    Ball(const ValidatedNumber& y, PR pr, PRE pre);

    //! Construct a ball containing the bounds \a x.
    //! The precision of the error is either that of \a F, or the default value.
    explicit Ball(const Bounds<F>& x);

    //! Assign from generic validated bounds \a y, keeping the same characteristics.
    Ball<F, FE>& operator=(const ValidatedNumber& y) { return *this = Ball<F, FE>(y, this->precision(), this->error_precision()); }

    //! Downcast to a generic validated value.
    operator ValidatedNumber() const;

    //! Create a ball from the generic validated number \a y with the same characteristics as \a this.
    Ball<F, FE> create(const ValidatedNumber& y) const { return Ball<F, FE>(y, this->precision(), this->error_precision()); }

    //! A lower bound.
    const LowerBound<F> lower() const;
    //! An upper bound.
    const UpperBound<F> upper() const;
    //! The centre value of the ball.
    const F value() const;
    //! The radius of the ball, giving a bound on the error of \ref value().
    const Error<FE> error() const;

    /** @brief Return the stored value. */
    friend F value(const Ball<F, FE>& x) { return x.value(); }
    /** @brief Return the error bound. */
    friend Error<FE> error(const Ball<F, FE>& x) { return x.error(); }

    /** @brief Return the underlying lower-endpoint representation. */
    const RawType lower_raw() const { return sub(down, _v, _e); }
    /** @brief Return the underlying upper-endpoint representation. */
    const RawType upper_raw() const { return add(up, _v, _e); }
    /** @brief Return the underlying value representation. */
    const RawType& value_raw() const { return _v; }
    /** @brief Return the underlying error representation. */
    const RawErrorType& error_raw() const { return _e; }
    //! Approximate by a builtin double-precision value. DEPRECATED
    double get_d() const { return _v.get_d(); }

    //! The precision of the floating-point type used for the centre value.
    PrecisionType precision() const { return _v.precision(); }
    //! The precision of the floating-point type used for the error bound.
    ErrorPrecisionType error_precision() const { return _e.precision(); }
    //! The compuational characteristics needed to create the ball; equivalent to the precision and error-precision.
    CharacteristicsType characteristics() const { return std::make_pair(_v.precision(), _e.precision()); }
    //! Downcast to a generic validated ball.
    GenericType generic() const { return this->operator GenericType(); }
    //! Add \a e to error bound.
    Ball<F, FE> pm(const Error<FE>& e) const;
    /** @brief Round the value to an integer. */
    friend Approximation<F> round(const Approximation<F>& x);

  public:
#ifdef DOXYGEN
    //!@{
    //! \name Arithmetic operators
    /** @brief Return the sum or unary positive value. */
    friend Ball<F, FE> operator+(const Ball<F, FE>& x);                         //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend Ball<F, FE> operator-(const Ball<F, FE>& x);                         //!< <p/>
    /** @brief Return the sum or unary positive value. */
    friend Ball<F, FE> operator+(const Ball<F, FE>& x1, const Ball<F, FE>& x2); //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend Ball<F, FE> operator-(const Ball<F, FE>& x1, const Ball<F, FE>& x2); //!< <p/>
    /** @brief Return the product. */
    friend Ball<F, FE> operator*(const Ball<F, FE>& x1, const Ball<F, FE>& x2); //!< <p/>
    /** @brief Return the quotient. */
    friend Ball<F, FE> operator/(const Ball<F, FE>& x1, const Ball<F, FE>& x2); //!< <p/>
    /** @brief Add the operand in place. */
    friend Ball<F, FE>& operator+=(Ball<F, FE>& x1, const Ball<F, FE>& x2);     //!< <p/>
    /** @brief Subtract the operand in place. */
    friend Ball<F, FE>& operator-=(Ball<F, FE>& x1, const Ball<F, FE>& x2);     //!< <p/>
    /** @brief Multiply by the operand in place. */
    friend Ball<F, FE>& operator*=(Ball<F, FE>& x1, const Ball<F, FE>& x2);     //!< <p/>
    /** @brief Divide by the operand in place. */
    friend Ball<F, FE>& operator/=(Ball<F, FE>& x1, const Ball<F, FE>& x2);     //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operators
    /** @brief Test equality of the operands. */
    friend ValidatedKleenean operator==(const Ball<F, FE>& x1, const Ball<F, FE>& x2); //!< <p/>
    /** @brief Test inequality of the operands. */
    friend ValidatedKleenean operator!=(const Ball<F, FE>& x1, const Ball<F, FE>& x2); //!< <p/>
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedKleenean operator<=(const Ball<F, FE>& x1, const Ball<F, FE>& x2); //!< <p/>
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedKleenean operator>=(const Ball<F, FE>& x1, const Ball<F, FE>& x2); //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const Ball<F, FE>& x1, const Ball<F, FE>& x2);  //!< <p/>
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(const Ball<F, FE>& x1, const Ball<F, FE>& x2);  //!< <p/>
    //!@}
#endif // DOXYGEN

    //!@{
    //! \name Arithmetic operations
    /** @brief Return the additive zero with the operand characteristics. */
    friend Ball<F, FE> nul(const Ball<F, FE>& x) {
        return Ball<F, FE>(nul(x._v), nul(x._e));
    } //!< <p/>
    /** @brief Return the unary positive value. */
    friend Ball<F, FE> pos(const Ball<F, FE>& x) {
        return Ball<F, FE>(pos(x._v), x._e);
    } //!< <p/>
    /** @brief Return the additive inverse. */
    friend Ball<F, FE> neg(const Ball<F, FE>& x) {
        return Ball<F, FE>(neg(x._v), x._e);
    } //!< <p/>
    /** @brief Divide the value by two. */
    friend Ball<F, FE> hlf(const Ball<F, FE>& x) {
        return Ball<F, FE>(hlf(x._v), hlf(x._e));
    } //!< <p/>
    /** @brief Return the square. */
    friend Ball<F, FE> sqr(const Ball<F, FE>& x) {
        return Operations<Ball<F, FE>>::_sqr(x);
    } //!< <p/>
    /** @brief Return the multiplicative inverse. */
    friend Ball<F, FE> rec(const Ball<F, FE>& x) {
        return Operations<Ball<F, FE>>::_rec(x);
    } //!< <p/>

    /** @brief Return the sum of the operands. */
    friend Ball<F, FE> add(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Operations<Ball<F, FE>>::_add(x1, x2);
    } //!< <p/>
    /** @brief Return the difference of the operands. */
    friend Ball<F, FE> sub(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Operations<Ball<F, FE>>::_sub(x1, x2);
    } //!< <p/>
    /** @brief Return the product of the operands. */
    friend Ball<F, FE> mul(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Operations<Ball<F, FE>>::_mul(x1, x2);
    } //!< <p/>
    /** @brief Return the quotient of the operands. */
    friend Ball<F, FE> div(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Operations<Ball<F, FE>>::_div(x1, x2);
    } //!< <p/>
    /** @brief Compute the product of the first two operands plus the third. */
    friend Ball<F, FE> fma(const Ball<F, FE>& x1, const Ball<F, FE>& x2, const Ball<F, FE>& x3) {
        return Operations<Ball<F, FE>>::_fma(x1, x2, x3);
    } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend Ball<F, FE> pow(const Ball<F, FE>& x, Int n) {
        return Operations<Ball<F, FE>>::_pow(x, n);
    } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend Ball<F, FE> pow(const Ball<F, FE>& x, Nat m) {
        return Operations<Ball<F, FE>>::_pow(x, m);
    } //!< <p/>
    //!@}

    //!@{
    //! \name Algebraic and transcendental operations
    /** @brief Return the square root. */
    friend Ball<F, FE> sqrt(const Ball<F, FE>& x) {
        return Ball<F, FE>(sqrt(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    /** @brief Return the exponential. */
    friend Ball<F, FE> exp(const Ball<F, FE>& x) {
        return Ball<F, FE>(exp(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    /** @brief Return the natural logarithm. */
    friend Ball<F, FE> log(const Ball<F, FE>& x) {
        return Ball<F, FE>(log(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    /** @brief Return the sine. */
    friend Ball<F, FE> sin(const Ball<F, FE>& x) {
        return Ball<F, FE>(sin(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    /** @brief Return the cosine. */
    friend Ball<F, FE> cos(const Ball<F, FE>& x) {
        return Ball<F, FE>(cos(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    /** @brief Return the tangent. */
    friend Ball<F, FE> tan(const Ball<F, FE>& x) {
        return Ball<F, FE>(tan(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    /** @brief Return the hyperbolic tangent. */
    friend Ball<F, FE> tanh(const Ball<F, FE>& x) {
        return Ball<F, FE>(tanh(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    /** @brief Return the inverse sine. */
    friend Ball<F, FE> asin(const Ball<F, FE>& x) {
        return Ball<F, FE>(asin(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    /** @brief Return the inverse cosine. */
    friend Ball<F, FE> acos(const Ball<F, FE>& x) {
        return Ball<F, FE>(acos(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    /** @brief Return the inverse tangent. */
    friend Ball<F, FE> atan(const Ball<F, FE>& x) {
        return Ball<F, FE>(atan(Bounds<F>(x)), x.error_precision());
    } //!< <p/>
    //!@}

    //!@{
    //! \name Lattice operations
    /** @brief Return the absolute value. */
    friend Ball<F, FE> abs(const Ball<F, FE>& x) {
        return Operations<Ball<F, FE>>::_abs(x);
    } //!< <p/>
    /** @brief Return the maximum of the operands. */
    friend Ball<F, FE> max(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Operations<Ball<F, FE>>::_max(x1, x2);
    } //!< <p/>
    /** @brief Return the minimum of the operands. */
    friend Ball<F, FE> min(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Operations<Ball<F, FE>>::_min(x1, x2);
    } //!< <p/>
    /** @brief Return a lower bound on the absolute value. */
    friend PositiveLowerBound<F> mig(const Ball<F, FE>& x) {
        return Operations<Ball<F, FE>>::_mig(x);
    } //!< <p/>
    /** @brief Return an upper bound on the absolute value. */
    friend PositiveUpperBound<F> mag(const Ball<F, FE>& x) {
        return Operations<Ball<F, FE>>::_mag(x);
    } //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operations
    /** @brief Return the sign of the value. */
    friend ValidatedKleenean sgn(const Ball<F, FE>& x) {
        if (x._v > x._e) {
            return true;
        } else if (x._v < -x._e) {
            return false;
        } else {
            return indeterminate;
        }
    } //!< <p/>
    //! \brief Equality comparison operator. Tests equality of represented real-point value.
    friend LogicalType<ValidatedTag> eq(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Operations<Ball<F, FE>>::_eq(x1, x2);
    } //!< <p/>
    //! \brief Strict less-than comparison operator. Tests equality of represented real-point value.
    friend LogicalType<ValidatedTag> lt(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Operations<Ball<F, FE>>::_lt(x1, x2);
    } //!< <p/>
      //!@}

    /*
    public:
        // Mixed Ball-ValidatedNumber operations
        template<AValidatedNumber Y> friend Ball<F,FE> operator+(Ball<F,FE> const& x1, Y const& y2) {
            if constexpr (AnExactDyadic<Y>) { return operator+(x1,F(y2,x1.precision())); }
            else { return operator+(x1,Ball<F,FE>(y2,x1.precision(),x1.error_precision())); } }
        template<AValidatedNumber Y> friend Ball<F,FE> operator-(Ball<F,FE> const& x1, Y const& y2) {
            if constexpr (AnExactDyadic<Y>) { return operator-(x1,F(y2,x1.precision())); }
            else { return operator-(x1,Ball<F,FE>(y2,x1.precision(),x1.error_precision())); } }
        template<AValidatedNumber Y> friend Ball<F,FE> operator*(Ball<F,FE> const& x1, Y const& y2) {
            if constexpr (AnExactDyadic<Y>) { return operator*(x1,F(y2,x1.precision())); }
            else { return operator*(x1,Ball<F,FE>(y2,x1.precision(),x1.error_precision())); } }
        template<AValidatedNumber Y> friend Ball<F,FE> operator/(Ball<F,FE> const& x1, Y const& y2) {
            if constexpr (AnExactDyadic<Y>) { return operator/(x1,F(y2,x1.precision())); }
            else { return operator/(x1,Ball<F,FE>(y2,x1.precision(),x1.error_precision())); } }
        template<AValidatedNumber Y> friend Ball<F,FE> operator+(Y const& y1, Ball<F,FE> const& x2) {
            if constexpr (AnExactDyadic<Y>) { return operator+(F(y1,x2.precision()),x2); }
            else { return operator+(Ball<F,FE>(y1,x2.precision(),x2.error_precision()),x2); } }
        template<AValidatedNumber Y> friend Ball<F,FE> operator-(Y const& y1, Ball<F,FE> const& x2) {
            if constexpr (AnExactDyadic<Y>) { return operator-(F(y1,x2.precision()),x2); }
            else { return operator-(Ball<F,FE>(y1,x2.precision(),x2.error_precision()),x2); } }
        template<AValidatedNumber Y> friend Ball<F,FE> operator*(Y const& y1, Ball<F,FE> const& x2) {
            if constexpr (AnExactDyadic<Y>) { return operator*(F(y1,x2.precision()),x2); }
            else { return operator*(Ball<F,FE>(y1,x2.precision(),x2.error_precision()),x2); } }
        template<AValidatedNumber Y> friend Ball<F,FE> operator/(Y const& y1, Ball<F,FE> const& x2) {
            if constexpr (AnExactDyadic<Y>) { return operator/(F(y1,x2.precision()),x2); }
            else { return operator/(Ball<F,FE>(y1,x2.precision(),x2.error_precision()),x2); } }
    */

  public:
    /** @brief Return the sum of the operands. */
    friend Ball<F, FE> add(const Ball<F, FE>& x1, const F& x2) { return add(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the difference of the operands. */
    friend Ball<F, FE> sub(const Ball<F, FE>& x1, const F& x2) { return sub(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the product of the operands. */
    friend Ball<F, FE> mul(const Ball<F, FE>& x1, const F& x2) { return mul(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the quotient of the operands. */
    friend Ball<F, FE> div(const Ball<F, FE>& x1, const F& x2) { return div(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the sum of the operands. */
    friend Ball<F, FE> add(const F& x1, const Ball<F, FE>& x2) { return add(Ball<F, FE>(x1, x2.error_precision()), x2); }
    /** @brief Return the difference of the operands. */
    friend Ball<F, FE> sub(const F& x1, const Ball<F, FE>& x2) { return sub(Ball<F, FE>(x1, x2.error_precision()), x2); }
    /** @brief Return the product of the operands. */
    friend Ball<F, FE> mul(const F& x1, const Ball<F, FE>& x2) { return mul(Ball<F, FE>(x1, x2.error_precision()), x2); }
    /** @brief Return the quotient of the operands. */
    friend Ball<F, FE> div(const F& x1, const Ball<F, FE>& x2) { return div(Ball<F, FE>(x1, x2.error_precision()), x2); }

    /** @brief Return the maximum of the operands. */
    friend Ball<F, FE> max(const Ball<F, FE>& x1, const F& x2) { return max(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the minimum of the operands. */
    friend Ball<F, FE> min(const Ball<F, FE>& x1, const F& x2) { return min(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the maximum of the operands. */
    friend Ball<F, FE> max(const F& x1, const Ball<F, FE>& x2) { return max(Ball<F, FE>(x1, x2.error_precision()), x2); }
    /** @brief Return the minimum of the operands. */
    friend Ball<F, FE> min(const F& x1, const Ball<F, FE>& x2) { return min(Ball<F, FE>(x1, x2.error_precision()), x2); }

    /** @brief Return the sum or unary positive value. */
    friend Ball<F, FE> operator+(const Ball<F, FE>& x1, const F& x2) { return operator+(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the difference or additive inverse. */
    friend Ball<F, FE> operator-(const Ball<F, FE>& x1, const F& x2) { return operator-(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the product. */
    friend Ball<F, FE> operator*(const Ball<F, FE>& x1, const F& x2) { return operator*(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the quotient. */
    friend Ball<F, FE> operator/(const Ball<F, FE>& x1, const F& x2) { return operator/(x1, Ball<F, FE>(x2, x1.error_precision())); }
    /** @brief Return the sum or unary positive value. */
    friend Ball<F, FE> operator+(const F& x1, const Ball<F, FE>& x2) { return operator+(Ball<F, FE>(x1, x2.error_precision()), x2); }
    /** @brief Return the difference or additive inverse. */
    friend Ball<F, FE> operator-(const F& x1, const Ball<F, FE>& x2) { return operator-(Ball<F, FE>(x1, x2.error_precision()), x2); }
    /** @brief Return the product. */
    friend Ball<F, FE> operator*(const F& x1, const Ball<F, FE>& x2) { return operator*(Ball<F, FE>(x1, x2.error_precision()), x2); }
    /** @brief Return the quotient. */
    friend Ball<F, FE> operator/(const F& x1, const Ball<F, FE>& x2) { return operator/(Ball<F, FE>(x1, x2.error_precision()), x2); }

    /** @brief Test equality of the operands. */
    friend ValidatedKleenean operator==(const Ball<F, FE>& x1, const F& x2) { return x1 == Ball<F, FE>(x2, x1.error_precision()); }
    /** @brief Test inequality of the operands. */
    friend ValidatedKleenean operator!=(const Ball<F, FE>& x1, const F& x2) { return x1 != Ball<F, FE>(x2, x1.error_precision()); }
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const Ball<F, FE>& x1, const F& x2) { return x1 < Ball<F, FE>(x2, x1.error_precision()); }
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(const Ball<F, FE>& x1, const F& x2) { return x1 > Ball<F, FE>(x2, x1.error_precision()); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedKleenean operator<=(const Ball<F, FE>& x1, const F& x2) { return x1 <= Ball<F, FE>(x2, x1.error_precision()); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedKleenean operator>=(const Ball<F, FE>& x1, const F& x2) { return x1 >= Ball<F, FE>(x2, x1.error_precision()); }
    /** @brief Test equality of the operands. */
    friend ValidatedKleenean operator==(const F& x1, const Ball<F, FE>& x2) { return Ball<F, FE>(x1, x2.error_precision()) == x2; }
    /** @brief Test inequality of the operands. */
    friend ValidatedKleenean operator!=(const F& x1, const Ball<F, FE>& x2) { return Ball<F, FE>(x1, x2.error_precision()) != x2; }
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const F& x1, const Ball<F, FE>& x2) { return Ball<F, FE>(x1, x2.error_precision()) < x2; }
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(const F& x1, const Ball<F, FE>& x2) { return Ball<F, FE>(x1, x2.error_precision()) > x2; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedKleenean operator<=(const F& x1, const Ball<F, FE>& x2) { return Ball<F, FE>(x1, x2.error_precision()) <= x2; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedKleenean operator>=(const F& x1, const Ball<F, FE>& x2) { return Ball<F, FE>(x1, x2.error_precision()) >= x2; }

  public:
    //!@{
    //! \name Rounding operations

    //! Round lower and upper bounds to nearest integer values.
    friend Ball<F, FE> round(const Ball<F, FE>& x) {
        return Ball<F, FE>(Bounds<F>(round(x.lower_raw()), round(x.upper_raw())), x.error_precision());
    }
    //! Round outward by 1 ulp. i.e. increase the error.
    friend Ball<F, FE> widen(const Ball<F, FE>& x) {
        const F m = F::min(x.precision());
        return Ball<F, FE>(Bounds<F>(sub(down, x.lower_raw(), m), add(up, x.upper_raw(), m)));
    }
    //! Round inward by 1 ulp. i.e. decrease the error.
    friend Ball<F, FE> narrow(const Ball<F, FE>& x) {
        const F m = F::min(x.precision());
        return Ball<F, FE>(Bounds<F>(add(up, x.lower_raw(), m), add(down, x.upper_raw(), m)));
    }
    //! Truncate to lower precision.
    friend Ball<F, FE> trunc(const Ball<F, FE>& x) {
        return Operations<Ball<F, FE>>::_trunc(x);
    }
    //! Truncate to lower precision.
    friend Ball<F, FE> trunc(const Ball<F, FE>& x, Nat n) {
        return Operations<Ball<F, FE>>::_trunc(x, n);
    }
    //! Round to the nearest integer to the midpoint.
    friend Integer cast_integer(const Ball<F, FE>& x) {
        return Operations<Ball<F, FE>>::_cast_integer(x);
    }
    //!@}

    //!@{
    //! \name Special value tests

    //! Tests whether \a x is NaN (not-a-number).
    friend Bool is_nan(const Ball<F, FE>& x) {
        return is_nan(x._v) || is_nan(x._e);
    }
    //! Tests whether \a x is a model of zero.
    friend auto is_zero(const Ball<F, FE>& x) -> LogicalType<ValidatedTag> {
        const F zero = nul(x.lower_raw());
        if (x.lower_raw() > zero || x.upper_raw() < zero) {
            return false;
        } else if (x.lower_raw() == zero && x.upper_raw() == zero) {
            return true;
        } else {
            return indeterminate;
        }
    }
    //! Tests whether \a x is a model of a positive number.
    friend auto is_positive(const Ball<F, FE>& x) -> LogicalType<ValidatedTag> {
        const F zero = nul(x.lower_raw());
        if (x.lower_raw() >= zero) {
            return true;
        } else if (x.upper_raw() < zero) {
            return false;
        } else {
            return indeterminate;
        }
    }
    //!@}

    //!@{
    //! \name Validated information tests and operations

    //! Tests is \a x1 and \a x2 have the same representation i.e. the same centre value and error bound.
    friend Bool same(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return x1._v == x2._v && x1._e == x2._e;
    }
    //! Tests is \a x1 is a valid approximation for \a x2 i.e. if \a x2 is within the error of the centre value of \a x1.
    friend Bool models(const Ball<F, FE>& x1, const F& x2) {
        return x1.lower_raw() <= x2 && x2 <= x1.upper_raw();
    }
    //! Tests is \a x1 and \a x2 are consistent with being a model of the same value i.e. they intersect.
    friend Bool consistent(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return (x1._v >= x2._v ? sub(up, x1._v, x2._v) : sub(up, x2._v, x1._v)) <= add(down, x1._e, x2._e);
    }
    //! Tests is \a x1 and \a x2 are inconsistent with being a model of the same value i.e. they are disjoint.
    friend Bool inconsistent(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return (x1._v >= x2._v ? sub(down, x1._v, x2._v) : sub(down, x2._v, x1._v)) > add(up, x1._e, x2._e);
    }
    //! Tests is \a x1 is a tighter representation than \a x2.
    friend Bool refines(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return add(up, (x1._v >= x2._v ? sub(up, x1._v, x2._v) : sub(up, x2._v, x1._v)), x1._e) <= x2._e;
    }
    //! The common refinement of \a x1 and \x2. Requires \a x1 and \a x2 to be consistent.
    friend Ball<F, FE> refinement(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Ball<F, FE>(refinement(Bounds<F>(x1), Bounds<F>(x2)), x1.error_precision());
    }
    /** @brief Return an enclosure containing both operands. */
    friend Ball<F, FE> coarsening(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Ball<F, FE>(coarsening(Bounds<F>(x1), Bounds<F>(x2)), x1.error_precision());
    }
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const Ball<F, FE>& x) {
        return Operations<Ball<F, FE>>::_write(os, x);
    } //!< Write to an output stream.
    friend InputStream& operator>>(InputStream& is, Ball<F, FE>& x) {
        return Operations<Ball<F, FE>>::_read(is, x);
    } //!< Read from an input stream.
      //!@}

  private:

  public:
    F _v;
    FE _e;
};

/** @brief Construct a \c Ball with the specified characteristics. */
template<class PR>
Ball(ValidatedNumber, PR) -> Ball<RawFloatType<PR>>;
/** @brief Construct a \c Ball with the specified characteristics. */
template<class PR, class PRE>
Ball(ValidatedNumber, PR, PRE) -> Ball<RawFloatType<PR>, RawFloatType<PRE>>;
/** @brief Construct a \c Ball with the specified characteristics. */
template<class F, class FE>
Ball(F, FE) -> Ball<F, FE>;

// Mixed Ball-value operations
/** @brief Return the sum of the operands. */
template<ARawFloat F, class PRE>
Ball<F, Float<PRE>> add(const F& x1, const F& x2, PRE pre) {
    typedef Float<PRE> FE;
    return Operations<Ball<F, FE>>::_add(x1, x2, pre);
}
/** @brief Return the difference of the operands. */
template<ARawFloat F, class PRE>
Ball<F, Float<PRE>> sub(const F& x1, const F& x2, PRE pre) {
    typedef Float<PRE> FE;
    return Operations<Ball<F, FE>>::_sub(x1, x2, pre);
}
/** @brief Return the product of the operands. */
template<ARawFloat F, class PRE>
Ball<F, Float<PRE>> mul(const F& x1, const F& x2, PRE pre) {
    typedef Float<PRE> FE;
    return Operations<Ball<F, FE>>::_mul(x1, x2, pre);
}
/** @brief Return the quotient of the operands. */
template<ARawFloat F, class PRE>
Ball<F, Float<PRE>> div(const F& x1, const F& x2, PRE pre) {
    typedef Float<PRE> FE;
    return Operations<Ball<F, FE>>::_div(x1, x2, pre);
}

/** @brief Return a factory for compatible objects. */
template<class F, class FE>
inline FloatBallFactory<PrecisionType<F>, PrecisionType<FE>> factory(const Ball<F, FE>& flt) {
    return FloatBallFactory<PrecisionType<F>, PrecisionType<FE>>(flt.precision(), flt.error_precision());
}
template<class PR, class PRE>
inline FloatBall<PR, PRE> FloatBallFactory<PR, PRE>::create(const Number<ValidatedTag>& y) { return FloatBall<PR, PRE>(y, this->_pr, this->_pre); }
template<class PR, class PRE>
inline FloatBall<PR, PRE> FloatBallFactory<PR, PRE>::create(const Number<EffectiveTag>& y) { return FloatBall<PR, PRE>(y, this->_pr, this->_pre); }
template<class PR, class PRE>
inline FloatBall<PR, PRE> FloatBallFactory<PR, PRE>::create(const Number<ExactTag>& y) { return FloatBall<PR, PRE>(y, this->_pr, this->_pre); }
template<class PR, class PRE>
inline FloatBall<PR, PRE> FloatBallFactory<PR, PRE>::create(const Real& y) { return FloatBall<PR, PRE>(y, this->_pr, this->_pre); }
template<class PR, class PRE>
inline FloatBall<PR, PRE> FloatBallFactory<PR, PRE>::create(const Rational& y) { return FloatBall<PR, PRE>(y, this->_pr, this->_pre); }
template<class PR, class PRE>
inline FloatBall<PR, PRE> FloatBallFactory<PR, PRE>::create(const Dyadic& y) { return FloatBall<PR, PRE>(y, this->_pr, this->_pre); }
template<class PR, class PRE>
inline FloatBall<PR, PRE> FloatBallFactory<PR, PRE>::create(const Integer& y) { return FloatBall<PR, PRE>(y, this->_pr, this->_pre); }

template<class PR, class PRE>
inline PositiveFloatBall<PR, PRE> FloatBallFactory<PR, PRE>::create(const PositiveValidatedNumber& y) { return PositiveFloatBall<PR, PRE>(y, this->_pr, this->_pre); }

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<class F, class FE>
class Positive<Ball<F, FE>> : public Ball<F, FE> {
    using PR = typename Ball<F, FE>::PrecisionType;
    using PRE = typename Ball<F, FE>::ErrorPrecisionType;

  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : Bounds<F>() {}
    /** @brief Construct a \c Positive from m, pre. */
    template<BuiltinUnsignedIntegral M>
    Positive(M m, PRE pre) : Ball<F, FE>(m, pre) {}
    /** @brief Construct a \c Positive from pr, pre. */
    explicit Positive(PR pr, PRE pre) : Ball<F, FE>(pr, pre) {}
    /** @brief Construct a \c Positive from x. */
    explicit Positive(const Ball<F, FE>& x) : Ball<F, FE>(x) {}
};

/** @brief Convert to the positive representation. */
template<class F, class FE>
inline PositiveBall<F, FE> cast_positive(const Ball<F, FE>& x) {
    return PositiveBall<F, FE>(x);
}

/** @brief Operations supported by  objects. */
template<class F, class FE>
struct Operations<Ball<F, FE>> {
    typedef typename FE::PrecisionType PRE;

    /** @brief Construct error from the supplied arguments. */
    static FE _make_error(const F& e) {
        static_assert(SameAs<F, FE> or DefaultConstructible<PRE>);
        if constexpr (SameAs<F, FE>) {
            return e;
        } else if constexpr (DefaultConstructible<PRE>) {
            return FE(e, up, PRE());
        }
    }

    /** @brief Truncate the value. */
    static Ball<F, FE> _trunc(const Ball<F, FE>& x);
    /** @brief Truncate the value. */
    static Ball<F, FE> _trunc(const Ball<F, FE>& x, Nat n);

    /** @brief Return the additive zero with the operand characteristics. */
    static Ball<F, FE> _nul(const Ball<F, FE>& x) {
        return Ball<F, FE>(nul(x._v), nul(x._e));
    }

    /** @brief Return the unary positive value. */
    static Ball<F, FE> _pos(const Ball<F, FE>& x) {
        return Ball<F, FE>(pos(x._v), x._e);
    }

    /** @brief Return the additive inverse. */
    static Ball<F, FE> _neg(const Ball<F, FE>& x) {
        return Ball<F, FE>(neg(x._v), x._e);
    }

    /** @brief Divide the value by two. */
    static Ball<F, FE> _hlf(const Ball<F, FE>& x) {
        return Ball<F, FE>(hlf(x._v), hlf(x._e));
    }

    /** @brief Return the square. */
    static Ball<F, FE> _sqr(const Ball<F, FE>& x) {
        Ball<F, FE> r = x * x;
        if (r._e > r._v) {
            r._e = hlf(add(up, r._e, _make_error(r._v)));
            r._v = F(Dyadic(r._e), upward, x.precision());
        }
        return r;
    }

    /** @brief Return the multiplicative inverse. */
    static Ball<F, FE> _rec(const Ball<F, FE>& x) {
        // Use this code to find value same as reciprocal value
        auto rv = rec(approx, x._v);
        auto ru = rec(up, sub(down, x._v, x._e));
        auto rl = rec(down, add(up, x._v, x._e));
        auto re = max(sub(up, ru, rv), sub(up, rv, rl));
        return Ball<F, FE>(rv, _make_error(re));
#ifdef ARIADNE_UNDEFINED
        // Use this code to get same result as interval computation
        auto ru = rec(up, sub(down, x._v, x._e));
        auto rl = rec(down, add(up, x._v, x._e));
        auto re = hlf(sub(up, ru, rl));
        auto rv = hlf(add(near, rl, ru));
        return Ball<F, FE>(rv, re);
#endif
    }

    /** @brief Return the sum of the operands. */
    static Ball<F, FE> _add(const Ball<F, FE>& x, const Ball<F, FE>& y) {
        auto rv = add(near, x._v, y._v);
        auto ru = add(up, x._v, y._v);
        auto rl = add(down, x._v, y._v);
        auto se = add(up, x._e, y._e);
        auto ae = hlf(FE(sub(up, ru, rl), up, se.precision()));
        auto re = add(up, ae, se);
        return Ball<F, FE>(rv, re);
    }

    /** @brief Return the difference of the operands. */
    static Ball<F, FE> _sub(const Ball<F, FE>& x, const Ball<F, FE>& y) {
        auto rv = sub(near, x._v, y._v);
        auto ru = sub(up, x._v, y._v);
        auto rl = sub(down, x._v, y._v);
        auto ae = _make_error(hlf(sub(up, ru, rl)));
        auto re = add(up, ae, add(up, x._e, y._e));
        return Ball<F, FE>(rv, re);
    }

    /** @brief Return the product of the operands. */
    static Ball<F, FE> _mul(const Ball<F, FE>& x, const Ball<F, FE>& y) {
        auto rv = mul(near, x._v, y._v);
        auto ru = mul(up, x._v, y._v);
        auto rl = mul(down, x._v, y._v);
        auto re0 = _make_error(hlf(sub(up, ru, rl)));
        auto re1 = add(up, re0, mul(up, x._e, y._e));
        auto re2 = add(up, mul(up, _make_error(abs(x._v)), y._e), mul(up, x._e, _make_error(abs(y._v))));
        auto re = add(up, re1, re2);
        return Ball<F, FE>(rv, re);
    }

    /** @brief Return the quotient of the operands. */
    static Ball<F, FE> _div(const Ball<F, FE>& x, const Ball<F, FE>& y) {
        return x * rec(y);
    }

    /** @brief Compute the product of the first two operands plus the third. */
    static Ball<F, FE> _fma(const Ball<F, FE>& x, const Ball<F, FE>& y, const Ball<F, FE>& z) {
        return _add(_mul(x, y), z);
    }

    /** @brief Return the sum of the operands. */
    static Ball<F, FE> _add(const F& x, const F& y, PRE pre) {
        auto rv = add(near, x, y);
        auto ru = add(up, x, y);
        auto rl = add(down, x, y);
        auto re = FE(hlf(sub(up, ru, rl)), up, pre);
        return Ball<F, FE>(rv, re);
    }

    /** @brief Return the difference of the operands. */
    static Ball<F, FE> _sub(const F& x, const F& y, PRE pre) {
        auto rv = sub(near, x, y);
        auto ru = sub(up, x, y);
        auto rl = sub(down, x, y);
        auto re = FE(hlf(sub(up, ru, rl)), up, pre);
        return Ball<F, FE>(rv, re);
    }

    /** @brief Return the product of the operands. */
    static Ball<F, FE> _mul(const F& x, const F& y, PRE pre) {
        auto rv = mul(near, x, y);
        auto ru = mul(up, x, y);
        auto rl = mul(down, x, y);
        auto re = FE(hlf(sub(up, ru, rl)), up, pre);
        return Ball<F, FE>(rv, re);
    }

    /** @brief Return the quotient of the operands. */
    static Ball<F, FE> _div(const F& x, const F& y, PRE pre) {
        auto rv = div(near, x, y);
        auto ru = div(up, x, y);
        auto rl = div(down, x, y);
        auto re = FE(hlf(sub(up, ru, rl)), up, pre);
        return Ball<F, FE>(rv, re);
    }

    /** @brief Raise the value to the specified power. */
    static Ball<F, FE> _pow(const Ball<F, FE>& x, Nat m) {
        return Ball<F, FE>(pow(Bounds<F>(x), m));
    }

    /** @brief Raise the value to the specified power. */
    static Ball<F, FE> _pow(const Ball<F, FE>& x, Int n) {
        return Ball<F, FE>(pow(Bounds<F>(x), n));
    }

    /** @brief Return the square root. */
    static Ball<F, FE> _sqrt(const Ball<F, FE>& x) {
        return Ball<F, FE>(sqrt(Bounds<F>(x)));
    }

    /** @brief Return the exponential. */
    static Ball<F, FE> _exp(const Ball<F, FE>& x) {
        return Ball<F, FE>(exp(Bounds<F>(x)));
    }

    /** @brief Return the natural logarithm. */
    static Ball<F, FE> _log(const Ball<F, FE>& x) {
        return Ball<F, FE>(log(Bounds<F>(x)));
    }

    /** @brief Return the sine. */
    static Ball<F, FE> _sin(const Ball<F, FE>& x) {
        return Ball<F, FE>(sin(Bounds<F>(x)));
    }

    /** @brief Return the cosine. */
    static Ball<F, FE> _cos(const Ball<F, FE>& x) {
        return Ball<F, FE>(cos(Bounds<F>(x)));
    }

    /** @brief Return the tangent. */
    static Ball<F, FE> _tan(const Ball<F, FE>& x) {
        return Ball<F, FE>(tan(Bounds<F>(x)));
    }

    /** @brief Return the inverse sine. */
    static Ball<F, FE> _asin(const Ball<F, FE>& x) {
        return Ball<F, FE>(asin(Bounds<F>(x)));
    }

    /** @brief Return the inverse cosine. */
    static Ball<F, FE> _acos(const Ball<F, FE>& x) {
        return Ball<F, FE>(acos(Bounds<F>(x)));
    }

    /** @brief Return the inverse tangent. */
    static Ball<F, FE> _atan(const Ball<F, FE>& x) {
        return Ball<F, FE>(atan(Bounds<F>(x)));
    }

    /** @brief Return the absolute value. */
    static Ball<F, FE> _abs(const Ball<F, FE>& x) {
        if (x._e < abs(x._v)) {
            return x;
        } else {
            auto rv = hlf(add(up, abs(x._v), x._e));
            return Ball<F, FE>(rv, _make_error(rv));
        }
    }

    /** @brief Return the maximum of the operands. */
    static Ball<F, FE> _max(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return hlf((x1 + x2) + abs(x1 - x2));
    }

    /** @brief Return the minimum of the operands. */
    static Ball<F, FE> _min(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return hlf((x1 + x2) - abs(x1 - x2));
    }

    /** @brief Return an upper bound on the absolute value. */
    static PositiveUpperBound<F> _mag(const Ball<F, FE>& x) {
        return PositiveUpperBound<F>(add(up, abs(x._v), x._e));
    }

    /** @brief Return a lower bound on the absolute value. */
    static PositiveLowerBound<F> _mig(const Ball<F, FE>& x) {
        return PositiveLowerBound<F>(max(F(0, x._v.precision()), sub(down, abs(x._v), x._e)));
    }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    static ValidatedKleenean _eq(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Bounds<F>(x1) == Bounds<F>(x2);
    }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    static ValidatedKleenean _lt(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Bounds<F>(x1) < Bounds<F>(x2);
    }

    /** @brief Test whether the operands have the same representation. */
    static Bool _same(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return x1._v == x2._v && x1._e == x2._e;
    }

    /** @brief Test whether the enclosure contains the supplied value. */
    static Bool _models(const Ball<F, FE>& x1, const F& x2) {
        return (x1._v >= x2 ? sub(up, x1._v, x2) : sub(up, x2, x1._v)) <= x1._e;
    }

    /** @brief Test whether the enclosures admit a common value. */
    static Bool _consistent(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return consistent(Bounds<F>(x1), Bounds<F>(x2));
    }

    /** @brief Test whether the enclosures cannot represent a common value. */
    static Bool _inconsistent(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return inconsistent(Bounds<F>(x1), Bounds<F>(x2));
    }

    /** @brief Test whether the first enclosure is contained in the second. */
    static Bool _refines(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return (x1._v >= x2._v ? sub(up, x1._v, x2._v) : sub(up, x2._v, x1._v)) <= sub(down, x2._e, x1._e);
    }

    /** @brief Return the common refinement of the enclosures. */
    static Ball<F, FE> _refinement(const Ball<F, FE>& x1, const Ball<F, FE>& x2) {
        return Ball<F, FE>(refinement(Bounds<F>(x1), Bounds<F>(x2)));
    }

    /** @brief Convert to an integer using the selected integer representation. */
    static Integer _cast_integer(const Ball<F, FE>& x);

    /** @brief Write the value to the output stream. */
    static OutputStream& _write(OutputStream& os, const Ball<F, FE>& x);

    /** @brief Read the value from the input stream. */
    static InputStream& _read(InputStream& is, Ball<F, FE>& x);
};

} // namespace Ariadne

#endif
