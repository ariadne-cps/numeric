/***************************************************************************
 *            numeric/float_bounds.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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
    Bounds(const LowerBound<F>& lower, const ValidatedUpperNumber& upper);
    Bounds(const ValidatedLowerNumber& lower, const UpperBound<F>& upper);
    //! Construct from a lower bound \a lower and an upper bound \a upper, using precsion \a pr.
    Bounds(const ValidatedLowerNumber& lower, const ValidatedUpperNumber& upper, PR pr);
    template<BuiltinIntegral N1, BuiltinIntegral N2>
    Bounds(N1 n1, N2 n2, PR pr) : _l(n1, pr), _u(n2, pr) {}
    Bounds(const ExactDouble& dl, const ExactDouble& du, PrecisionType pr) : _l(dl, pr), _u(du, pr) {}
    Bounds(const Dyadic& wl, const Dyadic& wu, PrecisionType pr) : _l(wl, down, pr), _u(wu, up, pr) {}
    Bounds(const Decimal& dl, const Decimal& du, PrecisionType pr) : _l(dl, down, pr), _u(du, up, pr) {}
    Bounds(const Rational& ql, const Rational& qu, PrecisionType pr) : _l(ql, down, pr), _u(qu, up, pr) {}

    template<class FF>
        requires ConstructibleFrom<F, FF, RND, PR>
    Bounds(const Bounds<FF>& x, PR pr) : Bounds(F(x.lower_raw(), downward, pr), F(x.upper_raw(), upward, pr)) {}
    Bounds(const Bounds<F>& x, PR pr) : _l(x._l, down, pr), _u(x._u, up, pr) {}

    template<BuiltinIntegral N>
    Bounds(N n, PR pr) : _l(n, down, pr), _u(n, up, pr) {}
    Bounds(const ExactDouble& d, PR pr) : _l(d, pr), _u(d, pr) {}
    Bounds(const TwoExp& t, PR pr) : _l(t, pr), _u(t, pr) {}
    Bounds(const Integer& z, PR pr) : _l(z, down, pr), _u(z, up, pr) {}
    Bounds(const Dyadic& w, PR pr) : _l(w, down, pr), _u(w, up, pr) {}
    Bounds(const Decimal& d, PR pr) : _l(d, down, pr), _u(d, up, pr) {}
    Bounds(const Rational& q, PR pr) : _l(q, down, pr), _u(q, up, pr) {}
    Bounds(const Real& x, PR pr);
    //! Construct from generic validated bounds \a y, using precsion \a pr.
    Bounds(const ValidatedNumber& y, PR pr);

    //! Convert from a ball.
    template<class FE>
    Bounds(const Ball<F, FE>& x);

    Bounds<F>& operator=(const F& x) { return *this = Bounds<F>(x); }
    //! Assign from generic validated bounds \a y, keeping the same precision.
    Bounds<F>& operator=(const ValidatedNumber& y) { return *this = Bounds<F>(y, this->precision()); }
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

    friend F value(const Bounds<F>& x) { return x.value(); }
    friend Error<F> error(const Bounds<F>& x) { return x.error(); }

    const RawType& lower_raw() const { return _l; }
    const RawType& upper_raw() const { return _u; }
    const RawType value_raw() const { return hlf(add(near, _l, _u)); }
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

    friend Approximation<F> round(const Approximation<F>& x);
    friend F midpoint(const Bounds<F>& x);

  public:
#ifdef DOXYGEN
    //!@{
    //! \name Arithmetic operators
    friend Bounds<F> operator+(const Bounds<F>& x);                       //!< <p/>
    friend Bounds<F> operator-(const Bounds<F>& x);                       //!< <p/>
    friend Bounds<F> operator+(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    friend Bounds<F> operator-(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    friend Bounds<F> operator*(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    friend Bounds<F> operator/(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    friend Bounds<F>& operator+=(Bounds<F>& x1, const Bounds<F>& x2);     //!< <p/>
    friend Bounds<F>& operator-=(Bounds<F>& x1, const Bounds<F>& x2);     //!< <p/>
    friend Bounds<F>& operator*=(Bounds<F>& x1, const Bounds<F>& x2);     //!< <p/>
    friend Bounds<F>& operator/=(Bounds<F>& x1, const Bounds<F>& x2);     //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operators
    friend ValidatedKleenean operator==(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    friend ValidatedKleenean operator!=(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    friend ValidatedKleenean operator<=(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    friend ValidatedKleenean operator>=(const Bounds<F>& x1, const Bounds<F>& x2); //!< <p/>
    friend ValidatedKleenean operator<(const Bounds<F>& x1, const Bounds<F>& x2);  //!< <p/>
    friend ValidatedKleenean operator>(const Bounds<F>& x1, const Bounds<F>& x2);  //!< <p/>
    //!@}
#endif // DOXYGEN

    //!@{
    //! \name Arithmetic operations
    friend Bounds<F> nul(const Bounds<F>& x) {
        return Bounds<F>(nul(x._l), nul(x._u));
    } //!< <p/>
    friend Bounds<F> pos(const Bounds<F>& x) {
        return Bounds<F>(pos(x._l), pos(x._u));
    } //!< <p/>
    friend Bounds<F> neg(const Bounds<F>& x) {
        return Bounds<F>(neg(x._u), neg(x._l));
    } //!< <p/>
    friend Bounds<F> hlf(const Bounds<F>& x) {
        return Bounds<F>(hlf(x._l), hlf(x._u));
    } //!< <p/>
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

    friend Bounds<F> add(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(add(down, x1._l, x2._l), add(up, x1._u, x2._u));
    } //!< <p/>
    friend Bounds<F> sub(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(sub(down, x1._l, x2._u), sub(up, x1._u, x2._l));
    } //!< <p/>
    friend Bounds<F> mul(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Operations<Bounds<F>>::_mul(x1, x2);
    } //!< <p/>
    friend Bounds<F> div(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Operations<Bounds<F>>::_div(x1, x2);
    } //!< <p/>
    friend Bounds<F> fma(const Bounds<F>& x1, const Bounds<F>& x2, const Bounds<F>& x3) {
        return Operations<Bounds<F>>::_fma(x1, x2, x3);
    } //!< <p/>
    friend Bounds<F> pow(const Bounds<F>& x, Int n) {
        if (n < 0) {
            return pow(rec(x), Nat(-n));
        } else
            return pow(x, Nat(n));
    } //!< <p/>
    friend Bounds<F> pow(const Bounds<F>& x, Nat m) {
        Bounds<F> y = (m % 2 == 0) ? abs(x) : x;
        Int n = static_cast<Int>(m);
        return Bounds<F>(pow(down, y._l, n), pow(up, y._u, n));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Algebraic and transcendental operations
    friend Bounds<F> sqrt(const Bounds<F>& x) {
        return Bounds<F>(sqrt(down, x.lower_raw()), sqrt(up, x.upper_raw()));
    } //!< <p/>
    friend Bounds<F> exp(const Bounds<F>& x) {
        return Bounds<F>(exp(down, x.lower_raw()), exp(up, x.upper_raw()));
    } //!< <p/>
    friend Bounds<F> log(const Bounds<F>& x) {
        return Bounds<F>(log(down, x.lower_raw()), log(up, x.upper_raw()));
    } //!< <p/>
    friend Bounds<F> sin(const Bounds<F>& x) {
        return Operations<Bounds<F>>::_sin(x);
    } //!< <p/>
    friend Bounds<F> cos(const Bounds<F>& x) {
        return Operations<Bounds<F>>::_cos(x);
    } //!< <p/>
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
    friend Bounds<F> tanh(const Bounds<F>& x) {
        return Bounds<F>(tanh(down, x.lower_raw()), tanh(up, x.upper_raw()));
    } //!< <p/>
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
    friend Bounds<F> atan(const Bounds<F>& x) {
        return Bounds<F>(atan(down, x._l), atan(up, x._u));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Lattice operations
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
    friend Bounds<F> max(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(max(x1.lower_raw(), x2.lower_raw()), max(x1.upper_raw(), x2.upper_raw()));
    } //!< <p/>
    friend Bounds<F> min(const Bounds<F>& x1, const Bounds<F>& x2) {
        return Bounds<F>(min(x1.lower_raw(), x2.lower_raw()), min(x1.upper_raw(), x2.upper_raw()));
    } //!< <p/>
    friend PositiveLowerBound<F> mig(const Bounds<F>& x) {
        return PositiveLowerBound<F>(max(F(0, x.precision()), max(x._l, neg(x._u))));
    } //!< <p/>
    friend PositiveUpperBound<F> mag(const Bounds<F>& x) {
        return PositiveUpperBound<F>(max(neg(x._l), x._u));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operations
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
    friend Bool same(const Bounds<F>& x, const Bounds<Dyadic>& w);
    friend Bool same(const Bounds<F>& x, const Dyadic& w) { return x._l == w && x._u == w; }
    friend Bool models(const Bounds<F>& x, const Dyadic& w) { return x._l <= w && w <= x._u; }
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
    friend Bounds<F> operator+(const F& x1, const F& x2) { return add(x1, x2); }
    friend Bounds<F> operator-(const F& x1, const F& x2) { return sub(x1, x2); }
    friend Bounds<F> operator*(const F& x1, const F& x2) { return mul(x1, x2); }
    friend Bounds<F> operator/(const F& x1, const F& x2) { return div(x1, x2); }
    friend Bounds<F> add(const F& x1, const F& x2) { return Bounds<F>{add(down, x1, x2), add(up, x1, x2)}; }
    friend Bounds<F> sub(const F& x1, const F& x2) { return Bounds<F>{sub(down, x1, x2), sub(up, x1, x2)}; }
    friend Bounds<F> mul(const F& x1, const F& x2) { return Bounds<F>{mul(down, x1, x2), mul(up, x1, x2)}; }
    friend Bounds<F> div(const F& x1, const F& x2) { return Bounds<F>{div(down, x1, x2), div(up, x1, x2)}; }
    friend Bounds<F> pow(const F& x, Nat m) { return Bounds<F>{pow(down, x, m), pow(up, x, m)}; }
    friend Bounds<F> pow(const F& x, Int n) { return Bounds<F>{pow(down, x, n), pow(up, x, n)}; }
    friend Bounds<F> sqr(const F& x) { return Bounds<F>{sqr(down, x), sqr(up, x)}; }
    friend Bounds<F> rec(const F& x) { return Bounds<F>{rec(down, x), rec(up, x)}; }
    friend Bounds<F> sqrt(const F& x) { return Bounds<F>{sqrt(down, x), sqrt(up, x)}; }
    friend Bounds<F> exp(const F& x) { return Bounds<F>{exp(down, x), exp(up, x)}; }
    friend Bounds<F> log(const F& x) { return Bounds<F>{log(down, x), log(up, x)}; }
    friend Bounds<F> sin(const F& x) { return Bounds<F>{sin(down, x), sin(up, x)}; }
    friend Bounds<F> cos(const F& x) { return Bounds<F>{cos(down, x), cos(up, x)}; }
    friend Bounds<F> tan(const F& x) { return tan(Bounds<F>(x)); }
    friend Bounds<F> asin(const F& x) { return Bounds<F>{asin(down, x), asin(up, x)}; }
    friend Bounds<F> acos(const F& x) { return Bounds<F>{acos(down, x), acos(up, x)}; }
    friend Bounds<F> atan(const F& x) { return Bounds<F>{atan(down, x), atan(up, x)}; }

  public:
    // Mixed Bounds-value operations
    friend Bounds<F> add(const Bounds<F>& x1, const F& x2) {
        return Bounds<F>(add(down, x1._l, x2), add(up, x1._u, x2));
    }
    friend Bounds<F> sub(const Bounds<F>& x1, const F& x2) {
        return Bounds<F>(sub(down, x1._l, x2), sub(up, x1._u, x2));
    }
    friend Bounds<F> mul(const Bounds<F>& x1, const F& x2) {
        return mul(x1, Bounds<F>(x2));
    }
    friend Bounds<F> div(const Bounds<F>& x1, const F& x2) {
        return div(x1, Bounds<F>(x2));
    }
    friend Bounds<F> add(const F& x1, const Bounds<F>& x2) {
        return Bounds<F>(add(down, x1, x2._l), add(up, x1, x2._u));
    }
    friend Bounds<F> sub(const F& x1, const Bounds<F>& x2) {
        return Bounds<F>(sub(down, x1, x2._u), sub(up, x1, x2._l));
    }
    friend Bounds<F> mul(const F& x1, const Bounds<F>& x2) {
        return mul(Bounds<F>(x1), x2);
    }
    friend Bounds<F> div(const F& x1, const Bounds<F>& x2) {
        return div(Bounds<F>(x1), x2);
    }

    friend Bounds<F> max(const Bounds<F>& x1, const F& x2) {
        return Bounds<F>(max(down, x1._l, x2), max(up, x1._u, x2));
    }
    friend Bounds<F> min(const Bounds<F>& x1, const F& x2) {
        return Bounds<F>(min(down, x1._l, x2), min(up, x1._u, x2));
    }
    friend Bounds<F> max(const F& x1, const Bounds<F>& x2) {
        return Bounds<F>(max(down, x1, x2._l), max(up, x1, x2._u));
    }
    friend Bounds<F> min(const F& x1, const Bounds<F>& x2) {
        return Bounds<F>(min(down, x1, x2._l), min(up, x1, x2._u));
    }

    friend Bounds<F> operator+(const Bounds<F>& x1, const F& x2) { return add(x1, x2); }
    friend Bounds<F> operator-(const Bounds<F>& x1, const F& x2) { return sub(x1, x2); }
    friend Bounds<F> operator*(const Bounds<F>& x1, const F& x2) { return mul(x1, x2); }
    friend Bounds<F> operator/(const Bounds<F>& x1, const F& x2) { return div(x1, x2); }
    friend Bounds<F> operator+(const F& x1, const Bounds<F>& x2) { return add(x1, x2); }
    friend Bounds<F> operator-(const F& x1, const Bounds<F>& x2) { return sub(x1, x2); }
    friend Bounds<F> operator*(const F& x1, const Bounds<F>& x2) { return mul(x1, x2); }
    friend Bounds<F> operator/(const F& x1, const Bounds<F>& x2) { return div(x1, x2); }

    friend ValidatedKleenean operator==(const Bounds<F>& x1, const F& x2) { return x1 == Bounds<F>(x2); }
    friend ValidatedKleenean operator!=(const Bounds<F>& x1, const F& x2) { return x1 != Bounds<F>(x2); }
    friend ValidatedKleenean operator<(const Bounds<F>& x1, const F& x2) { return x1 < Bounds<F>(x2); }
    friend ValidatedKleenean operator>(const Bounds<F>& x1, const F& x2) { return x1 > Bounds<F>(x2); }
    friend ValidatedKleenean operator<=(const Bounds<F>& x1, const F& x2) { return x1 <= Bounds<F>(x2); }
    friend ValidatedKleenean operator>=(const Bounds<F>& x1, const F& x2) { return x1 >= Bounds<F>(x2); }
    friend ValidatedKleenean operator==(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) == x2; }
    friend ValidatedKleenean operator!=(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) != x2; }
    friend ValidatedKleenean operator<(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) < x2; }
    friend ValidatedKleenean operator>(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) > x2; }
    friend ValidatedKleenean operator<=(const F& x1, const Bounds<F>& x2) { return Bounds<F>(x1) <= x2; }
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

template<class FE>
inline Bounds<FE> make_bounds(const Error<FE>& e) { return pm(e); }
FloatDP midpoint(const Bounds<FloatDP>& x); // DEPRECATED

template<ARawFloat F>
Bounds<F> sigmoid(const F& x) {
    return sigmoid(Bounds<F>(x));
}
template<ARawFloat F>
Bounds<F> softplus(const F& x) {
    return softplus(Bounds<F>(x));
}
template<ARawFloat F>
Bounds<F> tanh(const F& x) {
    return tanh(Bounds<F>(x));
}

template<class PR>
Bounds(ValidatedNumber, PR) -> Bounds<RawFloatType<PR>>;
template<class PR>
Bounds(ValidatedLowerNumber, ValidatedUpperNumber, PR) -> Bounds<RawFloatType<PR>>;
template<class F>
Bounds(F, F) -> Bounds<F>;

template<class F>
inline FloatFactory<PrecisionType<F>> factory(const Bounds<F>& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR>
inline FloatBounds<PR> FloatFactory<PR>::create(const Number<ValidatedTag>& y) { return FloatBounds<PR>(y, _pr); }
template<class PR>
inline FloatBounds<PR> FloatFactory<PR>::create(const Number<EffectiveTag>& y) { return FloatBounds<PR>(y, _pr); }
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

template<class F>
class Positive<Bounds<F>> : public Bounds<F>,
                            public DeclarePositiveFloatOperations<PositiveBounds<F>> {
    using typename Bounds<F>::PR;

  public:
    Positive() : Bounds<F>() {}
    explicit Positive(PR const& pr) : Bounds<F>(pr) {}
    template<BuiltinUnsignedIntegral M>
    Positive(M m, PR pr) : Bounds<F>(m, pr) {}
    explicit Positive(const F& x) : Bounds<F>(x) {}
    explicit Positive(const F& l, const F& u) : Bounds<F>(l, u) {}
    explicit Positive(const Bounds<F>& x) : Bounds<F>(x) {}
    Positive(const Positive<LowerBound<F>>& xl, const Positive<UpperBound<F>>& xu) : Bounds<F>(xl, xu) {}
    Positive(const PositiveValidatedNumber& y, PR pr) : Bounds<F>(y, pr) {}

  public:
    Positive<F> value() const { return cast_positive(this->Bounds<F>::value()); }
    Positive<LowerBound<F>> lower() const { return cast_positive(this->Bounds<F>::lower()); }
    Positive<UpperBound<F>> upper() const { return cast_positive(this->Bounds<F>::upper()); }

  public:
    friend PositiveBounds<F> nul(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(nul(x._l), nul(x._u));
    }
    friend PositiveBounds<F> hlf(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(hlf(x._l), hlf(x._u));
    }
    friend PositiveBounds<F> sqr(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(sqr(down, x._l), sqr(up, x._u));
    }
    friend PositiveBounds<F> rec(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(rec(down, x._u), rec(up, x._l));
    }
    friend PositiveBounds<F> add(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(add(down, x1._l, x2._l), add(up, x1._u, x2._u));
    }
    friend PositiveBounds<F> mul(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(mul(down, x1._l, x2._l), mul(up, x1._u, x2._u));
    }
    friend PositiveBounds<F> div(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(div(down, x1._l, x2._u), div(up, x1._u, x2._l));
    }
    friend PositiveBounds<F> pow(const PositiveBounds<F>& x, Nat m) {
        return PositiveBounds<F>(pow(down, x._l, static_cast<Int>(m)), pow(up, x._u, static_cast<Int>(m)));
    }
    friend PositiveBounds<F> sqrt(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(sqrt(down, x._l), sqrt(up, x._u));
    }
    friend PositiveBounds<F> atan(const PositiveBounds<F>& x) {
        return PositiveBounds<F>(atan(down, x._l), atan(up, x._u));
    }
    friend PositiveBounds<F> max(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(max(down, x1._l, x2._l), max(up, x1._u, x2._u));
    }
    friend PositiveBounds<F> max(const PositiveBounds<F>& x1, const Bounds<F>& x2) {
        return PositiveBounds<F>(max(down, x1._l, x2._l), max(up, x1._u, x2._u));
    }
    friend PositiveBounds<F> max(const Bounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(max(down, x1._l, x2._l), max(up, x1._u, x2._u));
    }
    friend PositiveBounds<F> min(const PositiveBounds<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveBounds<F>(min(down, x1._l, x2._l), min(up, x2._u, x2._u));
    }
    friend PositiveBounds<F> abs(const PositiveBounds<F>& x) { return PositiveBounds<F>(abs(x._l), abs(x._u)); }

    friend const Bounds<F>& cast_unsigned(const PositiveBounds<F>& x) { return x; }
};

template<class F>
inline PositiveBounds<F> cast_positive(const Bounds<F>& x) {
    return PositiveBounds<F>(x);
}

template<class F>
class Operations<Bounds<F>> {
    typedef PrecisionType<F> PR;

  public:
    static Bounds<F> _mul(const Bounds<F>& x1, const Bounds<F>& x2);
    static Bounds<F> _div(const Bounds<F>& x1, const Bounds<F>& x2);

    static Bounds<F> _fma(const Bounds<F>& x1, const Bounds<F>& x2, const Bounds<F>& x3);

    static Bounds<F> _pi(PR pr);
    static Bounds<F> _sin(const Bounds<F>& x);
    static Bounds<F> _cos(const Bounds<F>& x);
    static Bounds<F> _tan(const Bounds<F>& x);

    static Bounds<F> _trunc(const Bounds<F>& x);
    static Bounds<F> _trunc(const Bounds<F>& x, Nat n);

    static Integer _cast_integer(const Bounds<F>& x);

    // Mixed Bounded - Exact operations
    static Bounds<F> _add(const Bounds<F>& x1, const F& x2);
    static Bounds<F> _add(const F& x1, const Bounds<F>& x2);
    static Bounds<F> _sub(const Bounds<F>& x1, const F& x2);
    static Bounds<F> _sub(const F& x1, const Bounds<F>& x2);
    static Bounds<F> _mul(const Bounds<F>& x1, const F& x2);
    static Bounds<F> _mul(const F& x1, const Bounds<F>& x2);
    static Bounds<F> _div(const Bounds<F>& x1, const F& x2);
    static Bounds<F> _div(const F& x1, const Bounds<F>& x2);

    static OutputStream& _write(OutputStream& os, const Bounds<F>& x);
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
