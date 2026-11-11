/***************************************************************************
 *            numeric/float_bounds.hpp
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

/** @file numeric/float_bounds.hpp
 *  @brief Floating-point bounds for real numbers.
 */

#ifndef ARIADNE_FLOAT_BOUNDS_HPP
#define ARIADNE_FLOAT_BOUNDS_HPP

#include "float.decl.hpp"
#include "float_factory.hpp"
#include "float_operations.hpp"
#include "float_traits.hpp"
#include "foundation/logical.decl.hpp"
#include "number.decl.hpp"
#include "positive.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

struct DefaultTag;

//! \ingroup NumericModule
//! \brief Validated bounds on a number with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that direct construction from a floating-point number is prohibited, since <c>%FloatDPBounds(3.3,dp)</c> would the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%FloatDPBounds(3.3_decimal,dp)</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c Bounds use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c ValidatedKleenean, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[\underline{x},\overline{x}]\leq [\underline{y},\overline{y}]\f$ returns \c True if \f$\overline{x}\leq \underline{y}\f$, since in this case \f$x\leq x\f$ whenever \f$x_1\in[\underline{x},\overline{x}]\f$ and \f$y\in[\underline{y},\overline{y}]\f$, \c False if \f$\underline{x}>\overline{y}\f$, since in this case we know \f$x>y\f$, and \c Indeterminate otherwise, since in this case we can find \f$x,y\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[\underline{x},\overline{x}]\f$==\f$[\underline{y},\overline{y}]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//! To test equality of representation, use \c same(x,y)
//!
//! To obtain the lower and upper bounds of the possible values, use \c x.lower() and \c x.upper().
//! To obtain a best estimate of the value, use \c x.value(), which has an error at most \a x.error().
//! If \f$v\f$ and \f$e\f$ are the returned value and error for the bounds \f$[l,u]\f$, then it is guaranteed that \f$v-e\leq l\f$ and \f$v+e\geq u\f$ in exact arithmetic.
//!
//! To test if the bounds contain a number , use \c models(Bounds<F>,F), and to test if bounds are inconsistent use \c inconsistent(x,y), and to test if \c x provides a better approximation, use \c refines(x,y).
//! \sa Real, FloatDP, Float, Ball, UpperBound, LowerBound, Approximation.
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne validated bounds can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert FloatBounds literals of the form \c {a,b} to an FloatBounds in functions.
//!
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   FloatDPBounds({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   FloatDPBounds({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   FloatDPBounds([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
template<class F>
class Bounds
    : public DefineFloatOperations<Bounds<F>> {
  protected:
    typedef ValidatedTag P;
    typedef typename F::RoundingModeType RND;
    typedef typename F::PrecisionType PR;

  public:
    //! <p/>
    typedef P Paradigm;
    //! <p/>
    typedef Bounds<F> NumericType;
    //! <p/>
    typedef Number<P> GenericType;
    //! <p/>
    typedef F RawType;
    //! <p/>
    typedef PR PrecisionType;
    //! <p/>
    typedef PR CharacteristicsType;

  public:
    //! Construct bounds of zero with precision \a pr.
    explicit Bounds(PrecisionType pr) : _l(0.0_x, pr), _u(0.0_x, pr) {}
    //! Construct bounds with value \a v.
    Bounds(const RawType& v) : _l(v), _u(v) {}
    //! Construct a lower bound of value \a l and an upper bound of value \a u.
    explicit Bounds(const RawType& l, const RawType& u) : _l(l), _u(u) {
        ARIADNE_PRECONDITION(_l.precision() == _u.precision());
    }
    //! Construct from a lower bound \a lower and an upper bound \a upper.
    Bounds(const LowerBound<F>& lower, const UpperBound<F>& upper);
    /** @brief Construct a \c Bounds from lower, upper. */
    Bounds(const LowerBound<F>& lower, const ValidatedUpperNumber& upper);
    /** @brief Construct a \c Bounds from lower, upper. */
    Bounds(const ValidatedLowerNumber& lower, const UpperBound<F>& upper);
    //! Construct from a lower bound \a lower and an upper bound \a upper, using precsion \a pr.
    Bounds(const ValidatedLowerNumber& lower, const ValidatedUpperNumber& upper, PR pr);
    /** @brief Construct a \c Bounds from n1, n2, pr. */
    template<BuiltinIntegral N1, BuiltinIntegral N2>
    Bounds(N1 n1, N2 n2, PR pr) : _l(n1, pr), _u(n2, pr) {}
    /** @brief Construct a \c Bounds from dl, du, pr. */
    Bounds(const ExactDouble& dl, const ExactDouble& du, PrecisionType pr) : _l(dl, pr), _u(du, pr) {}
    /** @brief Construct a \c Bounds from wl, wu, pr. */
    Bounds(const Dyadic& wl, const Dyadic& wu, PrecisionType pr) : _l(wl, down, pr), _u(wu, up, pr) {}
    /** @brief Construct a \c Bounds from dl, du, pr. */
    Bounds(const Decimal& dl, const Decimal& du, PrecisionType pr) : _l(dl, down, pr), _u(du, up, pr) {}
    /** @brief Construct a \c Bounds from ql, qu, pr. */
    Bounds(const Rational& ql, const Rational& qu, PrecisionType pr) : _l(ql, down, pr), _u(qu, up, pr) {}

    /** @brief Construct a \c Bounds from x, pr. */
    template<class FF>
        requires ConstructibleFrom<F, FF, RND, PR>
    Bounds(const Bounds<FF>& x, PR pr) : Bounds(F(x.lower_raw(), downward, pr), F(x.upper_raw(), upward, pr)) {}
    /** @brief Construct a \c Bounds from x, pr. */
    Bounds(const Bounds<F>& x, PR pr) : _l(x._l, down, pr), _u(x._u, up, pr) {}

    /** @brief Construct a \c Bounds from n, pr. */
    template<BuiltinIntegral N>
    Bounds(N n, PR pr) : _l(n, down, pr), _u(n, up, pr) {}
    /** @brief Construct a \c Bounds from d, pr. */
    Bounds(const ExactDouble& d, PR pr) : _l(d, pr), _u(d, pr) {}
    /** @brief Construct a \c Bounds from t, pr. */
    Bounds(const TwoExp& t, PR pr) : _l(t, pr), _u(t, pr) {}
    /** @brief Construct a \c Bounds from z, pr. */
    Bounds(const Integer& z, PR pr) : _l(z, down, pr), _u(z, up, pr) {}
    /** @brief Construct a \c Bounds from w, pr. */
    Bounds(const Dyadic& w, PR pr) : _l(w, down, pr), _u(w, up, pr) {}
    /** @brief Construct a \c Bounds from d, pr. */
    Bounds(const Decimal& d, PR pr) : _l(d, down, pr), _u(d, up, pr) {}
    /** @brief Construct a \c Bounds from q, pr. */
    Bounds(const Rational& q, PR pr) : _l(q, down, pr), _u(q, up, pr) {}
    /** @brief Construct a \c Bounds from x, pr. */
    Bounds(const Real& x, PR pr);
    //! Construct from generic validated bounds \a y, using precsion \a pr.
    Bounds(const ValidatedNumber& y, PR pr);

    //! Convert from a ball.
    template<class FE>
    Bounds(const Ball<F, FE>& x);

    /** @brief Assign the supplied value and return this object. */
    Bounds<F>& operator=(const F& x) { return *this = Bounds<F>(x); }
    //! Assign from generic validated bounds \a y, keeping the same precision.
    Bounds<F>& operator=(const ValidatedNumber& y) { return *this = Bounds<F>(y, this->precision()); }
    /** @brief Assign the supplied value and return this object. */
    template<BuiltinIntegral N>
    Bounds<F>& operator=(N n) { return *this = ValidatedNumber(n); }

    //! Downcast to generic validated bounds.
    operator ValidatedNumber() const;

    //! Create bounds from the generic validated number \a y with the same precision as \a this.
    Bounds<F> create(const ValidatedNumber& y) const { return Bounds<F>(y, this->precision()); }

    //! The lower bound.
    const LowerBound<F> lower() const;
    //! The upper bound.
    const UpperBound<F> upper() const;
    //! An approximation to the actual value.
    const F value() const;
    //! Bounds on the approximation given by value().
    const Error<F> error() const;

    /** @brief Return the stored value. */
    friend F value(const Bounds<F>& x) { return x.value(); }
    /** @brief Return the error bound. */
    friend Error<F> error(const Bounds<F>& x) { return x.error(); }

    /** @brief Return the underlying lower-endpoint representation. */
    const RawType& lower_raw() const { return _l; }
    /** @brief Return the underlying upper-endpoint representation. */
    const RawType& upper_raw() const { return _u; }
    /** @brief Return the underlying value representation. */
    const RawType value_raw() const { return hlf(add(near, _l, _u)); }
    /** @brief Return the underlying error representation. */
    const RawType error_raw() const {
        RawType v = value_raw();
        return max(sub(up, _u, v), sub(up, v, _l));
    }
    //! Approximate by a builtin double-precision value. DEPRECATED
    double get_d() const { return value_raw().get_d(); }

    //! The precision of the floating-point type used.
    PrecisionType precision() const { return _u.precision(); }
    //! The compuational characteristics needed to create the bounds; equivalent to the precision.
    CharacteristicsType characteristics() const { return _u.precision(); }
    //! Downcast to generic validated bounds.
    GenericType generic() const { return this->operator GenericType(); }

    //! Add \a e to upper bound, and subtract \a e from lower bound.
    Bounds<F> pm(const Error<F>& e) const;

    /** @brief Round the value to an integer. */
    friend Approximation<F> round(const Approximation<F>& x);
    /** @brief Return the midpoint. */
    friend F midpoint(const Bounds<F>& x);

  public:
#ifdef DOXYGEN
    //!@{
    //! \name Arithmetic operators
    /** @brief Return the sum or unary positive value. */
    friend Bounds<F> operator+(const Bounds<F>& x);                       //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend Bounds<F> operator-(const Bounds<F>& x);                       //!< <p/>
    /** @brief Return the sum or unary positive value. */
    friend Bounds<F> operator+(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend Bounds<F> operator-(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    /** @brief Return the product. */
    friend Bounds<F> operator*(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    /** @brief Return the quotient. */
    friend Bounds<F> operator/(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    /** @brief Add the operand in place. */
    friend Bounds<F>& operator+=(Bounds<F>& x1, const Bounds<F>& x2);     //!< <p/>
    /** @brief Subtract the operand in place. */
    friend Bounds<F>& operator-=(Bounds<F>& x1, const Bounds<F>& x2);     //!< <p/>
    /** @brief Multiply by the operand in place. */
    friend Bounds<F>& operator*=(Bounds<F>& x1, const Bounds<F>& x2);     //!< <p/>
    /** @brief Divide by the operand in place. */
    friend Bounds<F>& operator/=(Bounds<F>& x1, const Bounds<F>& x2);     //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operators
    /** @brief Test equality of the operands. */
    friend ValidatedKleenean operator==(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    /** @brief Test inequality of the operands. */
    friend ValidatedKleenean operator!=(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedKleenean operator<=(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedKleenean operator>=(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const Bounds<F>& x1, const Bounds<F>& x2);  //!< <p/>
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(const Bounds<F>& x1, const Bounds<F>& x2);  //!< <p/>
    //!@}
#endif // DOXYGEN

    //!@{
    //! \name Arithmetic operations
    /** @brief Return the additive zero with the operand characteristics. */
    friend Bounds<F> nul(const Bounds<F>& x) {
        return Bounds<F>(nul(x._l), nul(x._u));
    } //!< <p/>
    /** @brief Return the unary positive value. */
    friend Bounds<F> pos(const Bounds<F>& x) {
        return Bounds<F>(pos(x._l), pos(x._u));
    } //!< <p/>
    /** @brief Return the additive inverse. */
    friend Bounds<F> neg(const Bounds<F>& x) {
        return Bounds<F>(neg(x._u), neg(x._l));
    } //!< <p/>
    /** @brief Divide the value by two. */
    friend Bounds<F> hlf(const Bounds<F>& x) {
        return Bounds<F>(hlf(x._l), hlf(x._u));
    } //!< <p/>
    /** @brief Return the square. */
    friend Bounds<F> sqr(const Bounds<F>& x) {
        if (x._l > 0) {
            return Bounds<F>(mul(down, x._l, x._l), mul(up, x._u, x._u));
        } //!< <p/>
        else if (x._u < 0) {
            return Bounds<F>(mul(down, x._u, x._u), mul(up, x._l, x._l));
        } //!< <p/>
        else {
            return Bounds<F>(nul(x._l), max(mul(up, x._l, x._l), mul(up, x._u, x._u)));
        }
    } //!< <p/>
    /** @brief Return the multiplicative inverse. */
    friend Bounds<F> rec(const Bounds<F>& x) {
        if (x._l > 0 || x._u < 0) {
            return Bounds<F>(rec(down, x._u), rec(up, x._l));
        } //!< <p/>
        if (x._l == 0 && x._u == 0) {
            F one(1, x.precision());
            return Bounds<F>(one, -one);
        }
        // The reciprocal image of an interval crossing zero is disconnected;
        // its interval hull is the whole real line.
        F inf_ = F::inf(x.precision());
        return Bounds<F>(-inf_, +inf_);
    } //!< <p/>

    /** @brief Return the sum of the operands. */
    friend Bounds<F> add(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(add(down, x1._l, x2._l), add(up, x1._u, x2._u));
    } //!< <p/>
    /** @brief Return the difference of the operands. */
    friend Bounds<F> sub(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(sub(down, x1._l, x2._u), sub(up, x1._u, x2._l));
    } //!< <p/>
    /** @brief Return the product of the operands. */
    friend Bounds<F> mul(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Operations<Bounds<F>>::_mul(x1, x2);
    } //!< <p/>
    /** @brief Return the quotient of the operands. */
    friend Bounds<F> div(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Operations<Bounds<F>>::_div(x1, x2);
    } //!< <p/>
    /** @brief Compute the product of the first two operands plus the third. */
    friend Bounds<F> fma(const Bounds<F>& x1, const Bounds<F>& x2, const Bounds<F>& x3) {
        return Operations<Bounds<F>>::_fma(x1, x2, x3);
    } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend Bounds<F> pow(const Bounds<F>& x, Int n) {
        if (n < 0) {
            return pow(rec(x), Nat(-n));
        } else
            return pow(x, Nat(n));
    } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend Bounds<F> pow(const Bounds<F>& x, Nat m) {
        Bounds<F> y = (m % 2 == 0) ? abs(x) : x;
        Int n = static_cast<Int>(m);
        return Bounds<F>(pow(down, y._l, n), pow(up, y._u, n));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Algebraic and transcendental operations
    /** @brief Return the square root. */
    friend Bounds<F> sqrt(const Bounds<F>& x) {
        return Bounds<F>(sqrt(down, x.lower_raw()), sqrt(up, x.upper_raw()));
    } //!< <p/>
    /** @brief Return the exponential. */
    friend Bounds<F> exp(const Bounds<F>& x) {
        return Bounds<F>(exp(down, x.lower_raw()), exp(up, x.upper_raw()));
    } //!< <p/>
    /** @brief Return the natural logarithm. */
    friend Bounds<F> log(const Bounds<F>& x) {
        return Bounds<F>(log(down, x.lower_raw()), log(up, x.upper_raw()));
    } //!< <p/>
    /** @brief Return the sine. */
    friend Bounds<F> sin(const Bounds<F>& x) {
        return Operations<Bounds<F>>::_sin(x);
    } //!< <p/>
    /** @brief Return the cosine. */
    friend Bounds<F> cos(const Bounds<F>& x) {
        return Operations<Bounds<F>>::_cos(x);
    } //!< <p/>
    /** @brief Return the tangent. */
    friend Bounds<F> tan(const Bounds<F>& x) {
        return Operations<Bounds<F>>::_tan(x);
    } //!< <p/>
    //! \brief Monotone rigorous enclosure of the logistic sigmoid.
    friend Bounds<F> sigmoid(const Bounds<F>& x) {
        return Bounds<F>(sigmoid(down, x.lower_raw()), sigmoid(up, x.upper_raw()));
    }
    //! \brief Monotone rigorous enclosure of softplus.
    friend Bounds<F> softplus(const Bounds<F>& x) {
        return Bounds<F>(softplus(down, x.lower_raw()), softplus(up, x.upper_raw()));
    }
    /** @brief Return the hyperbolic tangent. */
    friend Bounds<F> tanh(const Bounds<F>& x) {
        return Bounds<F>(tanh(down, x.lower_raw()), tanh(up, x.upper_raw()));
    } //!< <p/>
    /** @brief Return the inverse sine. */
    friend Bounds<F> asin(const Bounds<F>& x) {
        F minus_one(-1, x.precision());
        F plus_one(1, x.precision());
        F lower = max(x.lower_raw(), minus_one);
        F upper = min(x.upper_raw(), plus_one);
        if (lower > upper) {
            return Bounds<F>(plus_one, minus_one);
        }
        return Bounds<F>(asin(down, lower), asin(up, upper));
    } //!< <p/>
    /** @brief Return the inverse cosine. */
    friend Bounds<F> acos(const Bounds<F>& x) {
        F minus_one(-1, x.precision());
        F plus_one(1, x.precision());
        F lower = max(x.lower_raw(), minus_one);
        F upper = min(x.upper_raw(), plus_one);
        if (lower > upper) {
            return Bounds<F>(plus_one, minus_one);
        }
        return Bounds<F>(acos(down, upper), acos(up, lower));
    } //!< <p/>
    /** @brief Return the inverse tangent. */
    friend Bounds<F> atan(const Bounds<F>& x) {
        return Bounds<F>(atan(down, x._l), atan(up, x._u));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Lattice operations
    /** @brief Return the absolute value. */
    friend Bounds<F> abs(const Bounds<F>& x) {
        if (x.lower_raw() >= 0) {
            return Bounds<F>(x.lower_raw(), x.upper_raw());
        } //!< <p/>
        else if (x.upper_raw() <= 0) {
            return Bounds<F>(neg(x.upper_raw()), neg(x.lower_raw()));
        } //!< <p/>
        else {
            return Bounds<F>(F(0, x.precision()), max(neg(x.lower_raw()), x.upper_raw()));
        }
    } //!< <p/>
    /** @brief Return the maximum of the operands. */
    friend Bounds<F> max(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(max(x1.lower_raw(), x2.lower_raw()), max(x1.upper_raw(), x2.upper_raw()));
    } //!< <p/>
    /** @brief Return the minimum of the operands. */
    friend Bounds<F> min(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(min(x1.lower_raw(), x2.lower_raw()), min(x1.upper_raw(), x2.upper_raw()));
    } //!< <p/>
    /** @brief Return a lower bound on the absolute value. */
    friend PositiveLowerBound<F> mig(const Bounds<F>& x) {
        return PositiveLowerBound<F>(max(F(0, x.precision()), max(x._l, neg(x._u))));
    } //!< <p/>
    /** @brief Return an upper bound on the absolute value. */
    friend PositiveUpperBound<F> mag(const Bounds<F>& x) {
        return PositiveUpperBound<F>(max(neg(x._l), x._u));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operations
    /** @brief Return the sign of the value. */
    friend ValidatedKleenean sgn(const Bounds<F>& x) {
        if (x._l > 0) {
            return true;
        } else if (x._u < 0) {
            return false;
        } else {
            return indeterminate;
        }
    } //!< <p/>
    //! \brief Equality comparison operator. Tests equality of represented real-point value.
    friend LogicalType<ValidatedTag> eq(const Bounds<F>& x1, const Bounds<F>& x2) {
        if (x1.upper_raw() < x2.lower_raw() || x1.lower_raw() > x2.upper_raw()) {
            return false;
        } //!< <p/>
        else if (x1.lower_raw() == x2.upper_raw() && x1.upper_raw() == x2.lower_raw()) {
            return true;
        } //!< <p/>
        else {
            return indeterminate;
        }
    } //!< <p/>
    //! \brief Strict less-than comparison operator. Tests equality of represented real-point value.
    friend LogicalType<ValidatedTag> lt(const Bounds<F>& x1, const Bounds<F>& x2) {
        if (x1.upper_raw() < x2.lower_raw()) {
            return true;
        } //!< <p/>
        else if (x1.lower_raw() >= x2.upper_raw()) {
            return false;
        } //!< <p/>
        else {
            return indeterminate;
        }
    } //!< <p/>
    //!@}

    //!@{
    //! \name Rounding operations

    //! Round bounds to nearest integer values.
    friend Bounds<F> round(const Bounds<F>& x) {
        return Bounds<F>(round(x.lower_raw()), round(x.upper_raw()));
    }
    //! Round outward by 1 ulp.
    friend Bounds<F> widen(const Bounds<F>& x) {
        const F m = F::min(x.precision());
        return Bounds<F>(sub(down, x._l, m), add(up, x._u, m));
    }
    //! Round inward by 1 ulp.
    friend Bounds<F> narrow(const Bounds<F>& x) {
        const F m = F::min(x.precision());
        return Bounds<F>(add(up, x._l, m), add(down, x._u, m));
    }
    //! Truncate to lower precision.
    friend Bounds<F> trunc(const Bounds<F>& x) {
        return Operations<Bounds<F>>::_trunc(x);
    }
    //! Truncate to lower precision.
    friend Bounds<F> trunc(const Bounds<F>& x, Nat n) {
        return Operations<Bounds<F>>::_trunc(x, n);
    }
    //! Round to the nearest integer to the midpoint.
    friend Integer cast_integer(const Bounds<F>& x) {
        return Operations<Bounds<F>>::_cast_integer(x);
    }
    //@}

    //!@{
    //! \name Special value tests

    //! Tests whether \a x is NaN (not-a-number).
    friend Bool is_nan(const Bounds<F>& x) {
        return is_nan(x._l) || is_nan(x._u);
    }
    //! Tests whether \a x is a model of zero.
    friend auto is_zero(const Bounds<F>& x) -> LogicalType<ValidatedTag> {
        const F zero = nul(x.lower_raw());
        if (x.lower_raw() > zero || x.upper_raw() < zero) {
            return false;
        } else if (x.lower_raw() == zero && x.upper_raw() == zero) {
            return true;
        } else {
            return indeterminate;
        }
    }
    //! Tests whether \a x is a model of a positive number. Returns \c true if the lower bound is positive, \c false if the upper bound is strictly negative, and \c indeterminate otherwise.
    friend auto is_positive(const Bounds<F>& x) -> LogicalType<ValidatedTag> {
        const F zero = nul(x.lower_raw());
        if (x.lower_raw() >= zero) {
            return true;
        } else if (x.upper_raw() < zero) {
            return false;
        } else {
            return indeterminate;
        }
    } //!< <p/>
    //!@}

    //!@{
    //! \name Validated information tests and operations

    //! Tests is \a x1 and \a x2 have the same representation i.e. the same lower and upper bounds.
    friend Bool same(const Bounds<F>& x1, const Bounds<F>& x2) {
        return x1._l == x2._l && x1._u == x2._u;
    }
    //! Tests is \a x1 is a valid approximation for \a x2 i.e. if \a x2 is within the lower and upper bounds of \a x1.
    friend Bool models(const Bounds<F>& x1, const F& x2) {
        return x1._l <= x2 && x1._u >= x2;
    }
    //! Tests is \a x1 and \a x2 are consistent with being a model of the same value i.e. they intersect.
    friend Bool consistent(const Bounds<F>& x1, const Bounds<F>& x2) {
        return x1._l <= x2._u && x1._u >= x2._l;
    }
    //! Tests is \a x1 and \a x2 are inconsistent with being a model of the same value i.e. they are disjoint.
    friend Bool inconsistent(const Bounds<F>& x1, const Bounds<F>& x2) {
        return x1._l > x2._u || x1._u < x2._l;
    }
    //! Tests is \a x1 is a tighter representation than \a x2.
    friend Bool refines(const Bounds<F>& x1, const Bounds<F>& x2) {
        return x1._l >= x2._l && x1._u <= x2._u;
    }
    //! The common refinement of \a x1 and \x2. Requires \a x1 and \a x2 to be consistent.
    friend Bounds<F> refinement(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(max(x1._l, x2._l), min(x1._u, x2._u));
    }
    //! Bounds for a value which is in either \a x1 or \a x2.
    friend Bounds<F> coarsening(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(min(x1._l, x2._l), max(x1._u, x2._u));
    }
    //!@}

    // FIXME: Added functionality
    /** @brief Test whether the operands have the same representation. */
    friend Bool same(const Bounds<F>& x, const Bounds<Dyadic>& w);
    /** @brief Test whether the operands have the same representation. */
    friend Bool same(const Bounds<F>& x, const Dyadic& w) { return x._l == w && x._u == w; }
    /** @brief Test whether the enclosure contains the supplied value. */
    friend Bool models(const Bounds<F>& x, const Dyadic& w) { return x._l <= w && w <= x._u; }
    /** @brief Round the value to an integer. */
    friend Bounds<F> round(const Bounds<F>&);

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const Bounds<F>& x) {
        return Operations<Bounds<F>>::_write(os, x);
    } //!< Write to an output stream.
    friend InputStream& operator>>(InputStream& is, Bounds<F>& x) {
        return Operations<Bounds<F>>::_read(is, x);
    } //!< Read from an input stream.
      //!@}

  public:
    // Value operations returning Bounds
    /** @brief Return the sum or unary positive value. */
    friend Bounds<F> operator+(const F& x1, const F& x2) { return add(x1, x2); }
    /** @brief Return the difference or additive inverse. */
    friend Bounds<F> operator-(const F& x1, const F& x2) { return sub(x1, x2); }
    /** @brief Return the product. */
    friend Bounds<F> operator*(const F& x1, const F& x2) { return mul(x1, x2); }
    /** @brief Return the quotient. */
    friend Bounds<F> operator/(const F& x1, const F& x2) { return div(x1, x2); }
    /** @brief Return the sum of the operands. */
    friend Bounds<F> add(const F& x1, const F& x2) { return Bounds<F>{add(down, x1, x2), add(up, x1, x2)}; }
    /** @brief Return the difference of the operands. */
    friend Bounds<F> sub(const F& x1, const F& x2) { return Bounds<F>{sub(down, x1, x2), sub(up, x1, x2)}; }
    /** @brief Return the product of the operands. */
    friend Bounds<F> mul(const F& x1, const F& x2) { return Bounds<F>{mul(down, x1, x2), mul(up, x1, x2)}; }
    /** @brief Return the quotient of the operands. */
    friend Bounds<F> div(const F& x1, const F& x2) { return Bounds<F>{div(down, x1, x2), div(up, x1, x2)}; }
    /** @brief Raise the value to the specified power. */
    friend Bounds<F> pow(const F& x, Nat m) { return Bounds<F>{pow(down, x, m), pow(up, x, m)}; }
    /** @brief Raise the value to the specified power. */
    friend Bounds<F> pow(const F& x, Int n) { return Bounds<F>{pow(down, x, n), pow(up, x, n)}; }
    /** @brief Return the square. */
    friend Bounds<F> sqr(const F& x) { return Bounds<F>{sqr(down, x), sqr(up, x)}; }
    /** @brief Return the multiplicative inverse. */
    friend Bounds<F> rec(const F& x) { return Bounds<F>{rec(down, x), rec(up, x)}; }
    /** @brief Return the square root. */
    friend Bounds<F> sqrt(const F& x) { return Bounds<F>{sqrt(down, x), sqrt(up, x)}; }
    /** @brief Return the exponential. */
    friend Bounds<F> exp(const F& x) { return Bounds<F>{exp(down, x), exp(up, x)}; }
    /** @brief Return the natural logarithm. */
    friend Bounds<F> log(const F& x) { return Bounds<F>{log(down, x), log(up, x)}; }
    /** @brief Return the sine. */
    friend Bounds<F> sin(const F& x) { return Bounds<F>{sin(down, x), sin(up, x)}; }
    /** @brief Return the cosine. */
    friend Bounds<F> cos(const F& x) { return Bounds<F>{cos(down, x), cos(up, x)}; }
    /** @brief Return the tangent. */
    friend Bounds<F> tan(const F& x) { return tan(Bounds<F>(x)); }
    /** @brief Return the inverse sine. */
    friend Bounds<F> asin(const F& x) { return Bounds<F>{asin(down, x), asin(up, x)}; }
    /** @brief Return the inverse cosine. */
    friend Bounds<F> acos(const F& x) { return Bounds<F>{acos(down, x), acos(up, x)}; }
    /** @brief Return the inverse tangent. */
    friend Bounds<F> atan(const F& x) { return Bounds<F>{atan(down, x), atan(up, x)}; }

  public:
    // Mixed Bounds-value operations
    /** @brief Return the sum of the operands. */
    friend Bounds<F> add(const Bounds<F>& x1, const F& x2) {
        return Bounds<F>(add(down, x1._l, x2), add(up, x1._u, x2));
    }
    /** @brief Return the difference of the operands. */
    friend Bounds<F> sub(const Bounds<F>& x1, const F& x2) {
        return Bounds<F>(sub(down, x1._l, x2), sub(up, x1._u, x2));
    }
    /** @brief Return the product of the operands. */
    friend Bounds<F> mul(const Bounds<F>& x1, const F& x2) {
        return mul(x1, Bounds<F>(x2));
    }
    /** @brief Return the quotient of the operands. */
    friend Bounds<F> div(const Bounds<F>& x1, const F& x2) {
        return div(x1, Bounds<F>(x2));
    }
    /** @brief Return the sum of the operands. */
    friend Bounds<F> add(const F& x1, const Bounds<F>& x2) {
        return Bounds<F>(add(down, x1, x2._l), add(up, x1, x2._u));
    }
    /** @brief Return the difference of the operands. */
    friend Bounds<F> sub(const F& x1, const Bounds<F>& x2) {
        return Bounds<F>(sub(down, x1, x2._u), sub(up, x1, x2._l));
    }
    /** @brief Return the product of the operands. */
    friend Bounds<F> mul(const F& x1, const Bounds<F>& x2) {
        return mul(Bounds<F>(x1), x2);
    }
    /** @brief Return the quotient of the operands. */
    friend Bounds<F> div(const F& x1, const Bounds<F>& x2) {
        return div(Bounds<F>(x1), x2);
    }

    /** @brief Return the maximum of the operands. */
    friend Bounds<F> max(const Bounds<F>& x1, const F& x2) {
        return Bounds<F>(max(down, x1._l, x2), max(up, x1._u, x2));
    }
    /** @brief Return the minimum of the operands. */
    friend Bounds<F> min(const Bounds<F>& x1, const F& x2) {
        return Bounds<F>(min(down, x1._l, x2), min(up, x1._u, x2));
    }
    /** @brief Return the maximum of the operands. */
    friend Bounds<F> max(const F& x1, const Bounds<F>& x2) {
        return Bounds<F>(max(down, x1, x2._l), max(up, x1, x2._u));
    }
    /** @brief Return the minimum of the operands. */
    friend Bounds<F> min(const F& x1, const Bounds<F>& x2) {
        return Bounds<F>(min(down, x1, x2._l), min(up, x1, x2._u));
    }

    /** @brief Return the sum or unary positive value. */
    friend Bounds<F> operator+(const Bounds<F>& x1, const F& x2) { return add(x1, x2); }
    /** @brief Return the difference or additive inverse. */
    friend Bounds<F> operator-(const Bounds<F>& x1, const F& x2) { return sub(x1, x2); }
    /** @brief Return the product. */
    friend Bounds<F> operator*(const Bounds<F>& x1, const F& x2) { return mul(x1, x2); }
    /** @brief Return the quotient. */
    friend Bounds<F> operator/(const Bounds<F>& x1, const F& x2) { return div(x1, x2); }
    /** @brief Return the sum or unary positive value. */
    friend Bounds<F> operator+(const F& x1, const Bounds<F>& x2) { return add(x1, x2); }
    /** @brief Return the difference or additive inverse. */
    friend Bounds<F> operator-(const F& x1, const Bounds<F>& x2) { return sub(x1, x2); }
    /** @brief Return the product. */
    friend Bounds<F> operator*(const F& x1, const Bounds<F>& x2) { return mul(x1, x2); }
    /** @brief Return the quotient. */
    friend Bounds<F> operator/(const F& x1, const Bounds<F>& x2) { return div(x1, x2); }

    /** @brief Test equality of the operands. */
    friend ValidatedKleenean operator==(const Bounds<F>& x1, const F& x2) { return x1 == Bounds<F>(x2); }
    /** @brief Test inequality of the operands. */
    friend ValidatedKleenean operator!=(const Bounds<F>& x1, const F& x2) { return x1 != Bounds<F>(x2); }
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const Bounds<F>& x1, const F& x2) { return x1 < Bounds<F>(x2); }
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(const Bounds<F>& x1, const F& x2) { return x1 > Bounds<F>(x2); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedKleenean operator<=(const Bounds<F>& x1, const F& x2) { return x1 <= Bounds<F>(x2); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedKleenean operator>=(const Bounds<F>& x1, const F& x2) { return x1 >= Bounds<F>(x2); }
    /** @brief Test equality of the operands. */
    friend ValidatedKleenean operator==(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) == x2; }
    /** @brief Test inequality of the operands. */
    friend ValidatedKleenean operator!=(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) != x2; }
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) < x2; }
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) > x2; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedKleenean operator<=(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) <= x2; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedKleenean operator>=(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) >= x2; }

    /*
      public:
        // Mixed Bounds-ValidatedNumber operations
        template<AValidatedNumber Y> friend Bounds<F> operator+(Bounds<F> const& x1, Y const& y2) {
            if constexpr (AnExactDyadic<Y>) { return operator+(x1,F(y2,x1.precision())); }
            else { return operator+(x1,Bounds<F>(y2,x1.precision())); } }
        template<AValidatedNumber Y> friend Bounds<F> operator-(Bounds<F> const& x1, Y const& y2) {
            if constexpr (AnExactDyadic<Y>) { return operator-(x1,F(y2,x1.precision())); }
            else { return operator-(x1,Bounds<F>(y2,x1.precision())); } }
        template<AValidatedNumber Y> friend Bounds<F> operator*(Bounds<F> const& x1, Y const& y2) {
            if constexpr (AnExactDyadic<Y>) { return operator*(x1,F(y2,x1.precision())); }
            else { return operator*(x1,Bounds<F>(y2,x1.precision())); } }
        template<AValidatedNumber Y> friend Bounds<F> operator/(Bounds<F> const& x1, Y const& y2) {
            if constexpr (AnExactDyadic<Y>) { return operator/(x1,F(y2,x1.precision())); }
            else { return operator/(x1,Bounds<F>(y2,x1.precision())); } }
        template<AValidatedNumber Y> friend Bounds<F> operator+(Y const& y1, Bounds<F> const& x2) {
            if constexpr (AnExactDyadic<Y>) { return operator+(F(y1,x2.precision()),x2); }
            else { return operator+(Bounds<F>(y1,x2.precision()),x2); } }
        template<AValidatedNumber Y> friend Bounds<F> operator-(Y const& y1, Bounds<F> const& x2) {
            if constexpr (AnExactDyadic<Y>) { return operator-(F(y1,x2.precision()),x2); }
            else { return operator-(Bounds<F>(y1,x2.precision()),x2); } }
        template<AValidatedNumber Y> friend Bounds<F> operator*(Y const& y1, Bounds<F> const& x2) {
            if constexpr (AnExactDyadic<Y>) { return operator*(F(y1,x2.precision()),x2); }
            else { return operator*(Bounds<F>(y1,x2.precision()),x2); } }
        template<AValidatedNumber Y> friend Bounds<F> operator/(Y const& y1, Bounds<F> const& x2) {
            if constexpr (AnExactDyadic<Y>) { return operator/(F(y1,x2.precision()),x2); }
            else { return operator/(Bounds<F>(y1,x2.precision()),x2); } }
    */

  public:
    inline static Nat output_places = 8;
    //! %Set the number of output places used to display the number. DEPRECATED
    static Void set_output_places(Nat p) { output_places = p; }

  private:

  public:
    RawType _l, _u;
};

template<class F>
template<class FE>
Bounds<F>::Bounds(const Ball<F, FE>& x) : Bounds(x.lower_raw(), x.upper_raw()) {}

/** @brief Construct bounds from the supplied arguments. */
template<class FE>
inline Bounds<FE> make_bounds(const Error<FE>& e) { return pm(e); }
/** @brief Return the midpoint. */
FloatDP midpoint(const Bounds<FloatDP>& x); // DEPRECATED

/** @brief Evaluate one divided by one plus the exponential of the negated argument. */
template<ARawFloat F>
Bounds<F> sigmoid(const F& x) {
    return sigmoid(Bounds<F>(x));
}
/** @brief Evaluate the logarithm of one plus the exponential. */
template<ARawFloat F>
Bounds<F> softplus(const F& x) {
    return softplus(Bounds<F>(x));
}
/** @brief Return the hyperbolic tangent. */
template<ARawFloat F>
Bounds<F> tanh(const F& x) {
    return tanh(Bounds<F>(x));
}

/** @brief Construct a \c Bounds with the specified characteristics. */
template<class PR>
Bounds(ValidatedNumber, PR) -> Bounds<RawFloatType<PR>>;
/** @brief Construct a \c Bounds with the specified characteristics. */
template<class PR>
Bounds(ValidatedLowerNumber, ValidatedUpperNumber, PR) -> Bounds<RawFloatType<PR>>;
/** @brief Construct a \c Bounds with the specified characteristics. */
template<class F>
Bounds(F, F) -> Bounds<F>;

/** @brief Return a factory for compatible objects. */
template<class F>
inline FloatFactory<PrecisionType<F>> factory(const Bounds<F>& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
/** @brief Construct a compatible object from the supplied value and characteristics. */
template<class PR>
inline FloatBounds<PR> FloatFactory<PR>::create(const Number<ValidatedTag>& y) { return FloatBounds<PR>(y, _pr); }
/** @brief Construct a compatible object from the supplied value and characteristics. */
template<class PR>
inline FloatBounds<PR> FloatFactory<PR>::create(const Number<EffectiveTag>& y) { return FloatBounds<PR>(y, _pr); }
/** @brief Construct a compatible object from the supplied value and characteristics. */
template<class PR>
inline FloatBounds<PR> FloatFactory<PR>::create(const Number<ExactTag>& y) { return FloatBounds<PR>(y, _pr); }
template<class PR>
inline FloatBounds<PR> FloatFactory<PR>::create(const Real& y) { return FloatBounds<PR>(y, _pr); }
template<class PR>
inline FloatBounds<PR> FloatFactory<PR>::create(const Rational& y) { return FloatBounds<PR>(y, _pr); }
template<class PR>
inline FloatBounds<PR> FloatFactory<PR>::create(const Dyadic& y) { return FloatBounds<PR>(y, _pr); }
template<class PR>
inline FloatBounds<PR> FloatFactory<PR>::create(const Integer& y) { return FloatBounds<PR>(y, _pr); }

template<class PR>
inline PositiveFloatBounds<PR> FloatFactory<PR>::create(const PositiveValidatedNumber& y) { return PositiveFloatBounds<PR>(y, _pr); }

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<class F>
class Positive<Bounds<F>> : public Bounds<F>,
                            public DeclarePositiveFloatOperations<PositiveBounds<F>> {
    using typename Bounds<F>::PR;

  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : Bounds<F>() {}
    /** @brief Construct a \c Positive from pr. */
    explicit Positive(PR const& pr) : Bounds<F>(pr) {}
    /** @brief Construct a \c Positive from m, pr. */
    template<BuiltinUnsignedIntegral M>
    Positive(M m, PR pr) : Bounds<F>(m, pr) {}
    /** @brief Construct a \c Positive from x. */
    explicit Positive(const F& x) : Bounds<F>(x) {}
    /** @brief Construct a \c Positive from l, u. */
    explicit Positive(const F& l, const F& u) : Bounds<F>(l, u) {}
    /** @brief Construct a \c Positive from x. */
    explicit Positive(const Bounds<F>& x) : Bounds<F>(x) {}
    /** @brief Construct a \c Positive from xl, xu. */
    Positive(const Positive<LowerBound<F>>& xl, const Positive<UpperBound<F>>& xu) : Bounds<F>(xl, xu) {}
    /** @brief Construct a \c Positive from y, pr. */
    Positive(const PositiveValidatedNumber& y, PR pr) : Bounds<F>(y, pr) {}

  public:
    /** @brief Return the stored value. */
    Positive<F> value() const { return cast_positive(this->Bounds<F>::value()); }
    /** @brief Return the lower endpoint. */
    Positive<LowerBound<F>> lower() const { return cast_positive(this->Bounds<F>::lower()); }
    /** @brief Return the upper endpoint. */
    Positive<UpperBound<F>> upper() const { return cast_positive(this->Bounds<F>::upper()); }

  public:
    /** @brief Return the additive zero with the operand characteristics. */
    friend PositiveBounds<F> nul(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(nul(x._l), nul(x._u));
    }
    /** @brief Divide the value by two. */
    friend PositiveBounds<F> hlf(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(hlf(x._l), hlf(x._u));
    }
    /** @brief Return the square. */
    friend PositiveBounds<F> sqr(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(sqr(down, x._l), sqr(up, x._u));
    }
    /** @brief Return the multiplicative inverse. */
    friend PositiveBounds<F> rec(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(rec(down, x._u), rec(up, x._l));
    }
    /** @brief Return the sum of the operands. */
    friend PositiveBounds<F> add(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(add(down, x1._l, x2._l), add(up, x1._u, x2._u));
    }
    /** @brief Return the product of the operands. */
    friend PositiveBounds<F> mul(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(mul(down, x1._l, x2._l), mul(up, x1._u, x2._u));
    }
    /** @brief Return the quotient of the operands. */
    friend PositiveBounds<F> div(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(div(down, x1._l, x2._u), div(up, x1._u, x2._l));
    }
    /** @brief Raise the value to the specified power. */
    friend PositiveBounds<F> pow(const PositiveBounds<F>& x, Nat m) {
        return PositiveBounds<F>(pow(down, x._l, static_cast<Int>(m)), pow(up, x._u, static_cast<Int>(m)));
    }
    /** @brief Return the square root. */
    friend PositiveBounds<F> sqrt(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(sqrt(down, x._l), sqrt(up, x._u));
    }
    /** @brief Return the inverse tangent. */
    friend PositiveBounds<F> atan(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(atan(down, x._l), atan(up, x._u));
    }
    /** @brief Return the maximum of the operands. */
    friend PositiveBounds<F> max(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(max(down, x1._l, x2._l), max(up, x1._u, x2._u));
    }
    /** @brief Return the maximum of the operands. */
    friend PositiveBounds<F> max(const PositiveBounds<F>& x1, const Bounds<F>& x2) {
        return PositiveBounds<F>(max(down, x1._l, x2._l), max(up, x1._u, x2._u));
    }
    /** @brief Return the maximum of the operands. */
    friend PositiveBounds<F> max(const Bounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(max(down, x1._l, x2._l), max(up, x1._u, x2._u));
    }
    /** @brief Return the minimum of the operands. */
    friend PositiveBounds<F> min(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(min(down, x1._l, x2._l), min(up, x2._u, x2._u));
    }
    /** @brief Return the absolute value. */
    friend PositiveBounds<F> abs(const PositiveBounds<F>& x) { return PositiveBounds<F>(abs(x._l), abs(x._u)); }

    /** @brief Convert to the unsigned representation. */
    friend const Bounds<F>& cast_unsigned(const PositiveBounds<F>& x) { return x; }
};

/** @brief Convert to the positive representation. */
template<class F>
inline PositiveBounds<F> cast_positive(const Bounds<F>& x) {
    return PositiveBounds<F>(x);
}

/** @brief Operations supported by  objects. */
template<class F>
class Operations<Bounds<F>> {
    typedef PrecisionType<F> PR;

  public:
    /** @brief Return the product of the operands. */
    static Bounds<F> _mul(const Bounds<F>& x1, const Bounds<F>& x2);
    /** @brief Return the quotient of the operands. */
    static Bounds<F> _div(const Bounds<F>& x1, const Bounds<F>& x2);

    /** @brief Compute the product of the first two operands plus the third. */
    static Bounds<F> _fma(const Bounds<F>& x1, const Bounds<F>& x2, const Bounds<F>& x3);

    /** @brief Return pi at the requested precision. */
    static Bounds<F> _pi(PR pr);
    /** @brief Return the sine. */
    static Bounds<F> _sin(const Bounds<F>& x);
    /** @brief Return the cosine. */
    static Bounds<F> _cos(const Bounds<F>& x);
    /** @brief Return the tangent. */
    static Bounds<F> _tan(const Bounds<F>& x);

    /** @brief Truncate the value. */
    static Bounds<F> _trunc(const Bounds<F>& x);
    /** @brief Truncate the value. */
    static Bounds<F> _trunc(const Bounds<F>& x, Nat n);

    /** @brief Convert to an integer using the selected integer representation. */
    static Integer _cast_integer(const Bounds<F>& x);

    // Mixed Bounded - Exact operations
    /** @brief Return the sum of the operands. */
    static Bounds<F> _add(const Bounds<F>& x1, const F& x2);
    /** @brief Return the sum of the operands. */
    static Bounds<F> _add(const F& x1, const Bounds<F>& x2);
    /** @brief Return the difference of the operands. */
    static Bounds<F> _sub(const Bounds<F>& x1, const F& x2);
    /** @brief Return the difference of the operands. */
    static Bounds<F> _sub(const F& x1, const Bounds<F>& x2);
    /** @brief Return the product of the operands. */
    static Bounds<F> _mul(const Bounds<F>& x1, const F& x2);
    /** @brief Return the product of the operands. */
    static Bounds<F> _mul(const F& x1, const Bounds<F>& x2);
    /** @brief Return the quotient of the operands. */
    static Bounds<F> _div(const Bounds<F>& x1, const F& x2);
    /** @brief Return the quotient of the operands. */
    static Bounds<F> _div(const F& x1, const Bounds<F>& x2);

    /** @brief Write the value to the output stream. */
    static OutputStream& _write(OutputStream& os, const Bounds<F>& x);
    /** @brief Read the value from the input stream. */
    static InputStream& _read(InputStream& is, Bounds<F>& x);
};

template<class F>
inline auto Operations<Bounds<F>>::_mul(const Bounds<F>& x1, const Bounds<F>& x2) -> Bounds<F> {
    const F& x1l = x1._l;
    const F& x1u = x1._u;
    const F& x2l = x2._l;
    const F& x2u = x2._u;
    PR pr(min(x1.precision(), x2.precision()));
    F rl(pr), ru(pr);
    if (x1l >= 0) {
        if (x2l >= 0) {
            rl = mul(down, x1l, x2l);
            ru = mul(up, x1u, x2u);
        } else if (x2u <= 0) {
            rl = mul(down, x1u, x2l);
            ru = mul(up, x1l, x2u);
        } else {
            rl = mul(down, x1u, x2l);
            ru = mul(up, x1u, x2u);
        }
    } else if (x1u <= 0) {
        if (x2l >= 0) {
            rl = mul(down, x1l, x2u);
            ru = mul(up, x1u, x2l);
        } else if (x2u <= 0) {
            rl = mul(down, x1u, x2u);
            ru = mul(up, x1l, x2l);
        } else {
            rl = mul(down, x1l, x2u);
            ru = mul(up, x1l, x2l);
        }
    } else {
        if (x2l >= 0) {
            rl = mul(down, x1l, x2u);
            ru = mul(up, x1u, x2u);
        } else if (x2u <= 0) {
            rl = mul(down, x1u, x2l);
            ru = mul(up, x1l, x2l);
        } else {
            rl = min(mul(down, x1u, x2l), mul(down, x1l, x2u));
            ru = max(mul(up, x1l, x2l), mul(up, x1u, x2u));
        }
    }
    return Bounds<F>(rl, ru);
}

template<class F>
inline auto Operations<Bounds<F>>::_div(const Bounds<F>& x1, const Bounds<F>& x2) -> Bounds<F> {
    const F& x1l = x1.lower_raw();
    const F& x1u = x1.upper_raw();
    const F& x2l = x2.lower_raw();
    const F& x2u = x2.upper_raw();
    PR pr(min(x1.precision(), x2.precision()));
    F rl(pr), ru(pr);

    // IMPORTANT: Need to be careful when one of the bounds is 0, since if x2l=-0.0 and x1u>0, then x2l>=0 but x1u/x2l=-inf
    if (x2l > 0) {
        if (x1l >= 0) {
            rl = div(down, x1l, x2u);
            ru = div(up, x1u, x2l);
        } else if (x1u <= 0) {
            rl = div(down, x1l, x2l);
            ru = div(up, x1u, x2u);
        } else {
            rl = div(down, x1l, x2l);
            ru = div(up, x1u, x2l);
        }
    } else if (x2u < 0) {
        if (x1l >= 0) {
            rl = div(down, x1u, x2u);
            ru = div(up, x1l, x2l);
        } else if (x1u <= 0) {
            rl = div(down, x1u, x2l);
            ru = div(up, x1l, x2u);
        } else {
            rl = div(down, x1u, x2u);
            ru = div(up, x1l, x2u);
        }
    } else {
        // ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds x1, FloatBounds x2)","x1="<<x1<<", x2="<<x2);
        rl = -F::inf(pr);
        ru = +F::inf(pr);
    }
    return Bounds<F>(rl, ru);
}

template<class F>
inline auto Operations<Bounds<F>>::_fma(const Bounds<F>& x1, const Bounds<F>& x2, const Bounds<F>& x3) -> Bounds<F> {
    return add(mul(x1, x2), x3);
}

} // namespace Ariadne

#endif
