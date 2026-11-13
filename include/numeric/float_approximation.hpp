/***************************************************************************
 *            numeric/float_approximation.hpp
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

/** @file numeric/float_approximation.hpp
 *  @brief Floating-point approximations to real numbers.
 */

#ifndef ARIADNE_FLOAT_APPROXIMATION_HPP
#define ARIADNE_FLOAT_APPROXIMATION_HPP

#include "float.decl.hpp"
#include "float_factory.hpp"
#include "float_operations.hpp"
#include "float_traits.hpp"
#include "foundation/logical.decl.hpp"
#include "number.decl.hpp"
#include "positive.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating point number approximations to real numbers supporting approxiamate arithmetic.
//! \details
//! The \c Approximation<F> class represents approximate floating-point numbers.
//! Operations are performed approximately, with no guarantees on the output.
//! \sa Real, NaiveReal, FloatDP, Float, Ball, Bounds.
template<class F>
class Approximation
    : public DefineFloatOperations<Approximation<F>> {
  protected:
    typedef ApproximateTag P;
    typedef typename F::RoundingModeType RND;
    typedef typename F::PrecisionType PR;

  public:
    //! <p/>
    typedef ApproximateTag Paradigm;
    //! <p/>
    typedef Approximation<F> NumericType;
    //! <p/>
    typedef ApproximateNumber GenericType;
    //! <p/>
    typedef F RawType;
    //! <p/>
    typedef PR PrecisionType;
    //! <p/>
    typedef PR CharacteristicsType;

  public:
    //! <p/>
    /** @brief Construct a \c Approximation from pr. */
    explicit Approximation(PrecisionType pr) : _a(0.0_x, pr) {}
    //! <p/>
    /** @brief Construct a \c Approximation from a. */
    Approximation(const RawType& a) : _a(a) {}

    /** @brief Construct a \c Approximation from d, pr. */
    Approximation(double d, PR pr) : _a(cast_exact(d), near, pr) {}
    /** @brief Construct a \c Approximation from d, pr. */
    Approximation(ApproximateDouble d, PR pr) : _a(cast_exact(d), near, pr) {}
    /** @brief Construct a \c Approximation from d, pr. */
    Approximation(const ExactDouble& d, PR pr) : _a(d, pr) {}
    /** @brief Construct a \c Approximation from t, pr. */
    Approximation(const TwoExp& t, PR pr) : _a(t, pr) {}
    /** @brief Construct a \c Approximation from z, pr. */
    Approximation(const Integer& z, PR pr) : _a(z, near, pr) {}
    /** @brief Construct a \c Approximation from w, pr. */
    Approximation(const Dyadic& w, PR pr) : _a(w, near, pr) {}
    /** @brief Construct a \c Approximation from d, pr. */
    Approximation(const Decimal& d, PR pr) : _a(d, near, pr) {}
    /** @brief Construct a \c Approximation from q, pr. */
    Approximation(const Rational& q, PR pr) : _a(q, near, pr) {}
    /** @brief Construct a \c Approximation from r, pr. */
    Approximation(const Real& r, PR pr); // : _a(r.get(pr)) { }
    /** @brief Construct a \c Approximation from x, pr. */
    Approximation(const Approximation<F>& x, PR pr) : _a(x.raw(), near, pr) {}
    //! <p/>
    /** @brief Construct a \c Approximation from y, pr. */
    Approximation(const ApproximateNumber& y, PR pr); // : _a(y.get(ApproximateTag(),pr)) { }

    //! <p/>
    /** @brief Construct a \c Approximation from x, pr. */
    template<class FF>
        requires Constructible<F, FF, RND, PR>
    Approximation(const Approximation<FF>& x, PR pr) : _a(x.raw(), near, pr) {}

    //! <p/>
    /** @brief Construct a \c Approximation from x. */
    Approximation(const Error<F>& x); // FIXME: Remove
    //! <p/>
    /** @brief Construct a \c Approximation from x. */
    template<class FE>
    Approximation(const Ball<F, FE>& x);
    //! <p/>
    /** @brief Construct a \c Approximation from x. */
    Approximation(const Bounds<F>& x);
    //! <p/>
    /** @brief Construct a \c Approximation from x. */
    Approximation(const UpperBound<F>& x);
    //! <p/>
    /** @brief Construct a \c Approximation from x. */
    Approximation(const LowerBound<F>& x);

    //! <p/>
    /** @brief Assign the supplied value and return this object. */
    template<BuiltinIntegral N>
    Approximation<F>& operator=(N n) {
        this->_a = n;
        return *this;
    }
    //! <p/>
    /** @brief Assign the supplied value and return this object. */
    template<BuiltinFloatingPoint D>
    Approximation<F>& operator=(D x) {
        this->_a = ExactDouble(x);
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    Approximation<F>& operator=(const LowerBound<F>& x) { return *this = Approximation<F>(x); }
    /** @brief Assign the supplied value and return this object. */
    Approximation<F>& operator=(const UpperBound<F>& x) { return *this = Approximation<F>(x); }
    /** @brief Assign the supplied value and return this object. */
    Approximation<F>& operator=(const Bounds<F>& x) { return *this = Approximation<F>(x); }
    /** @brief Assign the supplied value and return this object. */
    Approximation<F>& operator=(const F& x) { return *this = Approximation<F>(x); }
    //! <p/>
    /** @brief Assign the supplied value and return this object. */
    Approximation<F>& operator=(const ApproximateNumber& y) { return *this = Approximation<F>(y, this->precision()); }
    //! <p/>
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    Approximation<F> create(const ApproximateNumber& y) const { return Approximation<F>(y, this->precision()); }

    //! <p/>
    /** @brief Convert to approximate number. */
    operator ApproximateNumber() const;

    //! <p/>
    /** @brief Return the working precision. */
    PrecisionType precision() const { return _a.precision(); }
    //! <p/>
    /** @brief Return the characteristics required to construct compatible values. */
    CharacteristicsType characteristics() const { return _a.precision(); }
    //! <p/>
    /** @brief Return the value through the generic numeric interface. */
    GenericType generic() const { return this->operator GenericType(); }
    //! <p/>
    /** @brief Convert to raw type. */
    explicit operator RawType() const { return this->_a; }
    //! <p/>
    /** @brief Convert to approximate double. */
    explicit operator ApproximateDouble() const { return this->_a.get_d(); }
    //! <p/>
    /** @brief Return the underlying numeric representation. */
    const RawType& raw() const { return this->_a; }
    //! <p/>
    /** @brief Return the underlying numeric representation. */
    RawType& raw() { return this->_a; }
    //! <p/> DEPRECATED
    double get_d() const { return this->_a.get_d(); }

  public:
#ifdef DOXYGEN
    //!@{
    //! \name Arithmetic operators
    /** @brief Return the sum or unary positive value. */
    friend Approximation<F> operator+(const Approximation<F>& x); //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend Approximation<F> operator-(const Approximation<F>& x); //!< <p/>
    /** @brief Return the sum or unary positive value. */
    friend Approximation<F> operator+(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend Approximation<F> operator-(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Return the product. */
    friend Approximation<F> operator*(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Return the quotient. */
    friend Approximation<F> operator/(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Add the operand in place. */
    friend Approximation<F>& operator+=(Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Subtract the operand in place. */
    friend Approximation<F>& operator-=(Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Multiply by the operand in place. */
    friend Approximation<F>& operator*=(Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Divide by the operand in place. */
    friend Approximation<F>& operator/=(Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operators
    /** @brief Test equality of the operands. */
    friend ApproximateKleenean operator==(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Test inequality of the operands. */
    friend ApproximateKleenean operator!=(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ApproximateKleenean operator<=(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ApproximateKleenean operator>=(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend ApproximateKleenean operator<(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is greater than the second. */
    friend ApproximateKleenean operator>(const Approximation<F>& x1, const Approximation<F>& x2); //!< <p/>
    //!@}
#endif // DOXYGEN

    //!@{
    //! \name Arithmetic operations
    /** @brief Return the additive zero with the operand characteristics. */
    friend Approximation<F> nul(const Approximation<F>& x) {
        return Approximation<F>(nul(x._a));
    } //!< <p/>
    /** @brief Return the unary positive value. */
    friend Approximation<F> pos(const Approximation<F>& x) {
        return Approximation<F>(pos(x._a));
    } //!< <p/>
    /** @brief Return the additive inverse. */
    friend Approximation<F> neg(const Approximation<F>& x) {
        return Approximation<F>(neg(x._a));
    } //!< <p/>
    /** @brief Divide the value by two. */
    friend Approximation<F> hlf(const Approximation<F>& x) {
        return Approximation<F>(hlf(x._a));
    } //!< <p/>
    /** @brief Return the square. */
    friend Approximation<F> sqr(const Approximation<F>& x) {
        return Approximation<F>(mul(near, x._a, x._a));
    } //!< <p/>
    /** @brief Return the multiplicative inverse. */
    friend Approximation<F> rec(const Approximation<F>& x) {
        return Approximation<F>(rec(near, x._a));
    } //!< <p/>

    /** @brief Return the sum of the operands. */
    friend Approximation<F> add(const Approximation<F>& x1, const Approximation<F>& x2) {
        return Approximation<F>(add(near, x1._a, x2._a));
    } //!< <p/>
    /** @brief Return the difference of the operands. */
    friend Approximation<F> sub(const Approximation<F>& x1, const Approximation<F>& x2) {
        return Approximation<F>(sub(near, x1._a, x2._a));
    } //!< <p/>
    /** @brief Return the product of the operands. */
    friend Approximation<F> mul(const Approximation<F>& x1, const Approximation<F>& x2) {
        return Approximation<F>(mul(near, x1._a, x2._a));
    } //!< <p/>
    /** @brief Return the quotient of the operands. */
    friend Approximation<F> div(const Approximation<F>& x1, const Approximation<F>& x2) {
        return Approximation<F>(div(near, x1._a, x2._a));
    } //!< <p/>
    /** @brief Compute the product of the first two operands plus the third. */
    friend Approximation<F> fma(const Approximation<F>& x1, const Approximation<F>& x2, const Approximation<F>& x3) {
        return Approximation<F>(fma(near, x1._a, x2._a, x3._a));
    } //!< <p/>

    /** @brief Raise the value to the specified power. */
    friend Approximation<F> pow(const Approximation<F>& x, Nat m) {
        return Approximation<F>(pow(approx, x._a, static_cast<Int>(m)));
    } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend Approximation<F> pow(const Approximation<F>& x, Int n) {
        return Approximation<F>(pow(approx, x._a, n));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Algebraic and transcendental operations
    /** @brief Return the square root. */
    friend Approximation<F> sqrt(const Approximation<F>& x) {
        return Approximation<F>(sqrt(approx, x._a));
    } //!< <p/>
    /** @brief Return the exponential. */
    friend Approximation<F> exp(const Approximation<F>& x) {
        return Approximation<F>(exp(approx, x._a));
    } //!< <p/>
    /** @brief Return the natural logarithm. */
    friend Approximation<F> log(const Approximation<F>& x) {
        return Approximation<F>(log(approx, x._a));
    } //!< <p/>
    /** @brief Return the sine. */
    friend Approximation<F> sin(const Approximation<F>& x) {
        return Approximation<F>(sin(approx, x._a));
    } //!< <p/>
    /** @brief Return the cosine. */
    friend Approximation<F> cos(const Approximation<F>& x) {
        return Approximation<F>(cos(approx, x._a));
    } //!< <p/>
    /** @brief Return the tangent. */
    friend Approximation<F> tan(const Approximation<F>& x) {
        return Approximation<F>(tan(approx, x._a));
    } //!< <p/>
    //! \brief Logistic sigmoid rounded approximately to the active precision.
    friend Approximation<F> sigmoid(const Approximation<F>& x) {
        return Approximation<F>(sigmoid(near, x._a));
    }
    //! \brief Softplus rounded approximately to the active precision.
    friend Approximation<F> softplus(const Approximation<F>& x) {
        return Approximation<F>(softplus(near, x._a));
    }
    /** @brief Return the hyperbolic tangent. */
    friend Approximation<F> tanh(const Approximation<F>& x) { return Approximation<F>(tanh(approx, x._a)); } //!< <p/>
    /** @brief Return the inverse sine. */
    friend Approximation<F> asin(const Approximation<F>& x) {
        return Approximation<F>(asin(approx, x._a));
    } //!< <p/>
    /** @brief Return the inverse cosine. */
    friend Approximation<F> acos(const Approximation<F>& x) {
        return Approximation<F>(acos(approx, x._a));
    } //!< <p/>
    /** @brief Return the inverse tangent. */
    friend Approximation<F> atan(const Approximation<F>& x) {
        return Approximation<F>(atan(approx, x._a));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Lattice operations
    /** @brief Return the absolute value. */
    friend Approximation<F> abs(const Approximation<F>& x) {
        return Approximation<F>(abs(x._a));
    } //!< <p/>
    /** @brief Return the maximum of the operands. */
    friend Approximation<F> max(const Approximation<F>& x, const Approximation<F>& y) {
        return Approximation<F>(max(x._a, y._a));
    } //!< <p/>
    /** @brief Return the minimum of the operands. */
    friend Approximation<F> min(const Approximation<F>& x, const Approximation<F>& y) {
        return Approximation<F>(min(x._a, y._a));
    } //!< <p/>
    /** @brief Return an upper bound on the absolute value. */
    friend PositiveApproximation<F> mag(const Approximation<F>& x) {
        return PositiveApproximation<F>(abs(x._a));
    } //!< <p/>
    /** @brief Return a lower bound on the absolute value. */
    friend PositiveApproximation<F> mig(const Approximation<F>& x) {
        return PositiveApproximation<F>(abs(x._a));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operations
    /** @brief Return the sign of the value. */
    friend ApproximateKleenean sgn(const Approximation<F>& x) {
        return ApproximateKleenean(x._a > 0 ? LogicalValue::LIKELY : x._a < 0 ? LogicalValue::UNLIKELY
                                                                              : LogicalValue::INDETERMINATE);
    } //!< <p/>
    /** @brief Test equality of the operands. */
    friend ApproximateKleenean eq(const Approximation<F>& x1, const Approximation<F>& x2) {
        return x1._a == x2._a;
    } //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend ApproximateKleenean lt(const Approximation<F>& x1, const Approximation<F>& x2) {
        return x1._a < x2._a;
    } //!< <p/>
    //@}

    //!@{
    //! \name Rounding operations
    /** @brief Return the greatest integer not exceeding the value. */
    friend Approximation<F> floor(const Approximation<F>& x) {
        return Approximation<F>(floor(x._a));
    } //!< <p/>
    /** @brief Return the least integer not less than the value. */
    friend Approximation<F> ceil(const Approximation<F>& x) {
        return Approximation<F>(ceil(x._a));
    } //!< <p/>
    /** @brief Round the value to an integer. */
    friend Approximation<F> round(const Approximation<F>& x) {
        return Approximation<F>(round(x._a));
    } //!< <p/>
    /** @brief Convert to an integer using the selected integer representation. */
    friend Integer cast_integer(const Approximation<F>& x) {
        return round(static_cast<Dyadic>(x._a));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Special value tests
    /** @brief Test whether the value is not a number. */
    friend Bool is_nan(const Approximation<F>& x) {
        return is_nan(x._a);
    } //!< <p/>
    //!@}

    //!@{
    //! \name Approximate information tests
    friend Bool same(const Approximation<F>& x1, const Approximation<F>& x2) {
        return x1._a == x2._a;
    } //!< Test if \a x1 and \a x2 are the same approximation. They might still possibly represent different numbers.
    //!}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const Approximation<F>& x) {
        os << "~";
        write(os, x.raw(), DecimalPrecision{Approximation<F>::output_places}, to_nearest);
        os << "~";
        return os;
    } //!< Write to an output stream.
    friend InputStream& operator>>(InputStream& is, Approximation<F>& x) {
        char c;
        is >> c;
        assert(c == '~');
        is >> x._a;
        is >> c;
        assert(c == '~');
        return is;
    } //!< Read from an input stream.
      //!@}

  public:
    inline static Nat output_places = 4;
    //! <p/>
    /** @brief Set the output places. */
    static Void set_output_places(Nat p) { output_places = p; }
    //! <p/>
    /** @brief Return an enclosure enlarged by the supplied symmetric error. */
    Approximation<F> pm(Approximation<F>) { return *this; }

  public:
    RawType _a;
};

template<class F>
template<class FE>
Approximation<F>::Approximation(const Ball<F, FE>& x) : Approximation(x.value_raw()) {}

/** @brief Construct a \c Approximation with the specified characteristics. */
template<class PR>
Approximation(ApproximateNumber, PR) -> Approximation<RawFloatType<PR>>;
/** @brief Construct a \c Approximation with the specified characteristics. */
template<class F>
Approximation(F) -> Approximation<F>;

/** @brief Return a factory for compatible objects. */
template<class F>
inline FloatFactory<PrecisionType<F>> factory(const Approximation<F>& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR>
inline FloatApproximation<PR> FloatFactory<PR>::create(const ApproximateNumber& y) { return FloatApproximation<PR>(y, _pr); }
template<class PR>
inline PositiveFloatApproximation<PR> FloatFactory<PR>::create(const PositiveApproximateNumber& y) { return PositiveFloatApproximation<PR>(y, _pr); }
template<class PR>
template<BuiltinFloatingPoint D>
inline FloatApproximation<PR> FloatFactory<PR>::create(const D& y) { return FloatApproximation<PR>(y, _pr); }

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<class F>
class Positive<Approximation<F>> : public Approximation<F>,
                                   public DeclarePositiveFloatOperations<PositiveApproximation<F>>,
                                   public DefineSemiFieldOperators<PositiveApproximation<F>>,
                                   public DefineConcreteGenericOperators<PositiveApproximation<F>> {
    using typename Approximation<F>::PR;

  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : Approximation<F>() {}
    /** @brief Construct a \c Positive from m. */
    template<BuiltinUnsignedIntegral M>
    Positive(M m) : Approximation<F>(m) {}
    /** @brief Construct a \c Positive from pr. */
    explicit Positive(PR const& pr) : Approximation<F>(pr) {}
    /** @brief Construct a \c Positive from x. */
    explicit Positive(const F& x) : Approximation<F>(x) {}
    /** @brief Construct a \c Positive from x. */
    explicit Positive(const Approximation<F>& x) : Approximation<F>(x) {}
    /** @brief Construct a \c Positive from y, pr. */
    Positive(const PositiveApproximateNumber& y, PR pr) : Approximation<F>(y, pr) {}
    /** @brief Copy-construct the object from the supplied instance. */
    Positive(const PositiveLowerBound<F>& x) : Approximation<F>(x) {}
    /** @brief Copy-construct the object from the supplied instance. */
    Positive(const PositiveUpperBound<F>& x) : Approximation<F>(x) {}
    /** @brief Copy-construct the object from the supplied instance. */
    Positive(const Positive<F>& x) : Approximation<F>(x) {}
    /** @brief Construct a \c Positive from x. */
    Positive(const Error<F>& x) : Approximation<F>(x) {}

  public:
    /** @brief Return the additive zero with the operand characteristics. */
    friend PositiveApproximation<F> nul(const PositiveApproximation<F>& x) { return PositiveApproximation<F>(nul(x.raw())); }
    /** @brief Divide the value by two. */
    friend PositiveApproximation<F> hlf(const PositiveApproximation<F>& x) { return PositiveApproximation<F>(hlf(x.raw())); }
    /** @brief Return the square. */
    friend PositiveApproximation<F> sqr(const PositiveApproximation<F>& x) { return PositiveApproximation<F>(sqr(near, x.raw())); }
    /** @brief Return the multiplicative inverse. */
    friend PositiveApproximation<F> rec(const PositiveApproximation<F>& x) { return PositiveApproximation<F>(rec(near, x.raw())); }
    /** @brief Return the sum of the operands. */
    friend PositiveApproximation<F> add(const PositiveApproximation<F>& x1, const PositiveApproximation<F>& x2) {
        return PositiveApproximation<F>(add(near, x1.raw(), x2.raw()));
    }
    /** @brief Return the product of the operands. */
    friend PositiveApproximation<F> mul(const PositiveApproximation<F>& x1, const PositiveApproximation<F>& x2) {
        return PositiveApproximation<F>(mul(near, x1.raw(), x2.raw()));
    }
    /** @brief Return the quotient of the operands. */
    friend PositiveApproximation<F> div(const PositiveApproximation<F>& x1, const PositiveApproximation<F>& x2) {
        return PositiveApproximation<F>(div(near, x1.raw(), x2.raw()));
    }
    /** @brief Raise the value to the specified power. */
    friend PositiveApproximation<F> pow(const PositiveApproximation<F>& x, Nat m) { return PositiveApproximation<F>(pow(near, x.raw(), static_cast<Int>(m))); }
    /** @brief Raise the value to the specified power. */
    friend PositiveApproximation<F> pow(const PositiveApproximation<F>& x, Int n) { return PositiveApproximation<F>(pow(near, x.raw(), n)); }
    /** @brief Return the square root. */
    friend PositiveApproximation<F> sqrt(const PositiveApproximation<F>& x) { return PositiveApproximation<F>(sqrt(near, x.raw())); }
    /** @brief Return the exponential. */
    friend PositiveApproximation<F> exp(const PositiveApproximation<F>& x) { return PositiveApproximation<F>(exp(near, x.raw())); }
    /** @brief Return the natural logarithm. */
    friend Approximation<F> log(const PositiveApproximation<F>& x) { return Approximation<F>(log(near, x.raw())); }
    /** @brief Return the inverse tangent. */
    friend PositiveApproximation<F> atan(const PositiveApproximation<F>& x) { return PositiveApproximation<F>(atan(near, x.raw())); }
    /** @brief Return the maximum of the operands. */
    friend PositiveApproximation<F> max(const PositiveApproximation<F>& x1, const PositiveApproximation<F>& x2) {
        return PositiveApproximation<F>(max(near, x1.raw(), x2.raw()));
    }
    /** @brief Return the maximum of the operands. */
    friend PositiveApproximation<F> max(const PositiveApproximation<F>& x1, const Approximation<F>& x2) {
        return PositiveApproximation<F>(max(near, x1.raw(), x2.raw()));
    }
    /** @brief Return the maximum of the operands. */
    friend PositiveApproximation<F> max(const Approximation<F>& x1, const PositiveApproximation<F>& x2) {
        return PositiveApproximation<F>(max(near, x1.raw(), x2.raw()));
    }
    /** @brief Return the minimum of the operands. */
    friend PositiveApproximation<F> min(const PositiveApproximation<F>& x1, const PositiveApproximation<F>& x2) {
        return PositiveApproximation<F>(min(near, x1.raw(), x2.raw()));
    }
    /** @brief Return the absolute value. */
    friend PositiveApproximation<F> abs(const PositiveApproximation<F>& x) { return PositiveApproximation<F>(abs(x)); }
};

/** @brief Convert to the positive representation. */
template<class F>
inline PositiveApproximation<F> cast_positive(const Approximation<F>& x) {
    return PositiveApproximation<F>(x);
}

} // namespace Ariadne

#endif
