/***************************************************************************
 *            numeric/complex.hpp
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

/** @file numeric/complex.hpp
 *  @brief Complex numbers.
 */

#ifndef ARIADNE_COMPLEX_HPP
#define ARIADNE_COMPLEX_HPP

#include <stdexcept>

#include "numeric/api.hpp"
#include "numeric/integer.hpp"
#include "numeric/real.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

template<class X>
class Positive;

template<>
class Positive<Real> : public Real {
  public:
    using Real::Real;
    Positive(const Real& r) : Real(r) {}
};

template<class X>
class Complex;
template<class X>
Complex(const X& re, const X& im) -> Complex<X>;

template<class X>
struct IsComplex : False {};
template<class X>
struct IsComplex<Complex<X>> : True {};

template<class X>
concept AComplex = IsComplex<X>::value;
template<class X>
concept AReal = not AComplex<X>;

struct PolarTag {};

Real sqrt(const Integer& z);

template<class X>
using ModulusType = decltype(sqrt(add(sqr(declval<X>()), sqr(declval<X>()))));
template<class X>
using ArgumentType = decltype(atan(div(declval<X>(), declval<X>())));

double atan(double);
template<class X>
inline auto atan2(const X& x, const X& y) -> decltype(atan(y / x)) {
    if (decide(x > 0)) {
        return atan(y / x);
    } else if (decide(y > 0)) {
        return pi / 2 - atan(x / y);
    } else if (decide(y < 0)) {
        return pi / (-2) - atan(x / y);
    } else {
        ARIADNE_THROW(std::runtime_error, "atan2(x,y)", "x=" << x << " and y=" << y << " could both be zero.");
    }
}

namespace Detail {
template<class X>
decltype(auto) complex_mig(const X& re, const X& im) {
    return cast_positive(sqrt(add(sqr(mig(re)), sqr(mig(im)))));
}
} // namespace Detail

class DefineComplexOperations {
    template<class X1, class X2>
    friend Complex<SumType<X1, X2>> add(const Complex<X1>& z1, const Complex<X2>& z2) {
        return Complex<SumType<X1, X2>>(z1._re + z2._re, z1._im + z2._im);
    } //!< \brief Sum \a z1+z2.
    template<class X1, class X2>
    friend Complex<DifferenceType<X1, X2>> sub(const Complex<X1>& z1, const Complex<X2>& z2) {
        return Complex<DifferenceType<X1, X2>>(z1._re - z2._re, z1._im - z2._im);
    } //!< \brief Difference \a z1-z2.
    template<class X1, class X2>
    friend Complex<ArithmeticType<X1, X2>> mul(const Complex<X1>& z1, const Complex<X2>& z2) {
        return Complex<ArithmeticType<X1, X2>>(z1._re * z2._re - z1._im * z2._im, z1._re * z2._im + z1._im * z2._re);
    } //!< \brief Product \a z1×z2.
    template<class X1, class X2>
    friend Complex<ArithmeticType<X1, X2>> div(const Complex<X1>& z1, const Complex<X2>& z2) {
        auto ns = add(sqr(z2._re), sqr(z2._im));
        return Complex<ArithmeticType<X1, X2>>((z1._re * z2._re + z1._im * z2._im) / ns, (z1._im * z2._re - z1._re * z2._im) / ns);
    } //!< \brief Quotient \a z1÷z2.

    template<class X1, class X2>
        requires AReal<X2>
    friend Complex<ProductType<X1, X2>> add(const Complex<X1>& z1, X2 const& x2) {
        return Complex<SumType<X1, X2>>(z1._re + x2, z1._im);
    }
    template<class X1, class X2>
        requires AReal<X1>
    friend Complex<ProductType<X1, X2>> add(X1 const& x1, const Complex<X2>& z2) {
        return Complex<SumType<X1, X2>>(x1 + z2._re, z2._im);
    }
    template<class X1, class X2>
        requires AReal<X2>
    friend Complex<ProductType<X1, X2>> sub(const Complex<X1>& z1, X2 const& x2) {
        return Complex<DifferenceType<X1, X2>>(z1._re - x2, z1._im);
    }
    template<class X1, class X2>
        requires AReal<X1>
    friend Complex<ProductType<X1, X2>> sub(X1 const& x1, const Complex<X2>& z2) {
        return Complex<DifferenceType<X1, X2>>(x1 - z2._re, -z2._im);
    }
    template<class X1, class X2>
        requires AReal<X2>
    friend Complex<ProductType<X1, X2>> mul(const Complex<X1>& z1, X2 const& x2) {
        return Complex<ProductType<X1, X2>>(z1._re * x2, z1._im * x2);
    }
    template<class X1, class X2>
        requires AReal<X1>
    friend Complex<ProductType<X1, X2>> mul(X1 const& x1, const Complex<X2>& z2) {
        return Complex<ProductType<X1, X2>>(x1 * z2._re, x1 * z2._im);
    }
    template<class X1, class X2>
        requires AReal<X1>
    friend Complex<QuotientType<X1, X2>> div(X1 const& x1, const Complex<X2>& z2) {
        auto ns = add(sqr(z2._re), sqr(z2._im));
        return Complex<ArithmeticType<X1, X2>>(x1 * z2._re / ns, -x1 * z2._im / ns);
    }
    template<class X1, class X2>
        requires AReal<X2>
    friend Complex<QuotientType<X1, X2>> div(const Complex<X1>& z1, X2 const& x2) {
        return Complex<ArithmeticType<X1, X2>>(z1._re / x2, z1._im / x2);
    }

    template<class X1, class X2>
    friend decltype(auto) operator+(const Complex<X1>& z1, const Complex<X2>& z2) { return add(z1, z2); }
    template<class X1, class X2>
    friend decltype(auto) operator-(const Complex<X1>& z1, const Complex<X2>& z2) { return sub(z1, z2); }
    template<class X1, class X2>
    friend decltype(auto) operator*(const Complex<X1>& z1, const Complex<X2>& z2) { return mul(z1, z2); }
    template<class X1, class X2>
    friend decltype(auto) operator/(const Complex<X1>& z1, const Complex<X2>& z2) { return div(z1, z2); }

    template<class X1, class X2>
    friend decltype(auto) operator+(const Complex<X1>& z1, X2 const& x2) { return add(z1, x2); }
    template<class X1, class X2>
    friend decltype(auto) operator+(X1 const& x1, const Complex<X2>& z2) { return add(x1, z2); }
    template<class X1, class X2>
    friend decltype(auto) operator-(const Complex<X1>& z1, X2 const& x2) { return sub(z1, x2); }
    template<class X1, class X2>
    friend decltype(auto) operator-(X1 const& x1, const Complex<X2>& z2) { return sub(x1, z2); }
    template<class X1, class X2>
    friend decltype(auto) operator*(const Complex<X1>& z1, X2 const& x2) { return mul(z1, x2); }
    template<class X1, class X2>
    friend decltype(auto) operator*(X1 const& x1, const Complex<X2>& z2) { return mul(x1, z2); }
    template<class X1, class X2>
    friend decltype(auto) operator/(const Complex<X1>& z1, X2 const& x2) { return div(z1, x2); }
    template<class X1, class X2>
    friend decltype(auto) operator/(X1 const& x1, const Complex<X2>& z2) { return div(x1, z2); }

    //!@{
    //! \name Comparison operations and operators.
    template<class X1, class X2>
    friend EqualityType<X1, X2> operator==(const Complex<X1>& z1, const Complex<X2>& z2) {
        return (z1._re == z2._re) && (z1._im == z2._im);
    } //!< Equality.
    template<class X1, class X2>
    friend InequalityType<X1, X2> operator!=(const Complex<X1>& z1, const Complex<X2>& z2) {
        return (z1._re != z2._re) || (z1._im != z2._im);
    } //!< Inequality.

    template<class X1, class X2>
        requires AReal<X2>
    friend EqualityType<X1, X2> operator==(const Complex<X1>& z1, X2 const& x2) {
        return (z1._re == x2) && (z1._im == 0);
    } //!< Equality.
    template<class X1, class X2>
        requires AReal<X2>
    friend EqualityType<X1, X2> operator==(X1 const& x1, const Complex<X2>& z2) {
        return (x1 == z2._re) && (0 == z2._im);
    } //!< Equality.

    template<class X1, class X2>
        requires AReal<X2>
    friend InequalityType<X1, X2> operator!=(const Complex<X1>& z1, X2 const& x2) {
        return (z1._re != x2) || (z1._im != 0);
    } //!< Equality.
    template<class X1, class X2>
        requires AReal<X2>
    friend InequalityType<X1, X2> operator!=(X1 const& x1, const Complex<X2>& z2) {
        return (x1 != z2._re) || (0 != z2._im);
    } //!< Equality.

    //!@}

    template<class X>
    friend Complex<X> nul(const Complex<X>& z) {
        return Complex<X>(nul(z._re), nul(z._im));
    } //!< Zero value \a 0.
    template<class X>
    friend Complex<X> pos(const Complex<X>& z) {
        return Complex<X>(pos(z._re), pos(z._im));
    } //!< Identity \a +z.
    template<class X>
    friend Complex<NegationType<X>> neg(const Complex<X>& z) {
        return Complex<NegationType<X>>(-z._re, -z._im);
    } //!< Negative \a -z.
    template<class X>
    friend Complex<X> hlf(const Complex<X>& z) {
        return Complex<X>(hlf(z._re), hlf(z._im));
    } //!< Half \a z÷2.
    template<class X>
    friend Complex<ArithmeticType<X>> sqr(const Complex<X>& z) {
        return Complex<ArithmeticType<X>>(sqr(z._re) - sqr(z._im), 2 * z._re * z._im);
    } //!< Square \a z<sup>2</sup>.
    template<class X>
    friend Complex<ArithmeticType<X>> rec(const Complex<X>& z) {
        auto ns = add(sqr(z._re), sqr(z._im));
        return Complex<ArithmeticType<X>>(z._re / ns, -z._im / ns);
    } //!< Reciprocal \a 1/z.
    template<class X>
    friend Complex<ArithmeticType<X>> pow(const Complex<X>& z, Int n); //!< \brief Power \a z<sup>n</sup>.
};

template<class X>
using TranscendentalType = decltype(sin(declval<X>()));

template<class X>
Complex<TranscendentalType<X>> make_complex_from_polar(const X& r, const X& th) {
    typedef decltype(exp(r) * cos(th)) R;
    return Complex<R>(exp(r) * cos(th), exp(r) * sin(th));
}

//! \ingroup NumericModule
//! \brief %Complex number type over real numbers \a X supporting basic arithmetic and exponential functions.
template<class X>
class Complex
    : DefineComplexOperations {
  private:

  public:
    X _re;
    X _im;
    friend class DefineComplexOperations;

  public:
    typedef typename X::Paradigm Paradigm;
    typedef X RealType;

  public:
    //!@{
    //! \name Constructors
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    Complex(PRS... prs) : _re(0u, prs...), _im(0u, prs...) {} //!< Default constructor yields the number 0.
    Complex(const X& x) : _re(x), _im(nul(x)) {}              //!< Construct the number \a x + <i>i</i> \a y.
    Complex(const X& x, const X& y) : _re(x), _im(y) {}       //!< Construct the number \a x + <i>i</i> \a y.

    // template<class POL> requires Same<POL,PolarTag> and Convertible<TranscendentalType<X>,X>
    //     Complex(X const& r, X const& th, POL) : _re(exp(r)*cos(th)), _im(exp(r)*sin(th)) { } //!< Construct the number \a exp(r)*(cos(th)+i*sin(th)).

    Complex(const X& r, const X& th, PolarTag) : _re(exp(r) * cos(th)), _im(exp(r) * sin(th)) {} //!< Construct the number \a exp(r)*(cos(th)+i*sin(th)).

    template<class Y>
        requires Convertible<Y, X>
    Complex(const Complex<Y>& z) : Complex(X(z._re), X(z._im)) {}
    template<class Y>
        requires Convertible<Y, X>
    Complex(const Y& y) : Complex(X(y)) {}
    template<class Y, class... PRS>
        requires Constructible<X, Y, PRS...>
    Complex(const Y& x, const Y& y, PRS... prs) : _re(x, prs...), _im(y, prs...) {}
    template<class Y, class... PRS>
        requires Constructible<X, Y, PRS...>
    Complex(const Complex<Y>& z, PRS... prs) : _re(z.real_part(), prs...), _im(z.imaginary_part(), prs...) {}
    template<class Y>
        requires Assignable<X, Y>
    Complex<X>& operator=(const Y& c) {
        _re = c;
        _im = 0;
        return *this;
    }
    //!@}

    //!@{
    //! \name Extract in Cartesian and polar coordinates
    const X& real_part() const { return _re; }      //!< The real part \a x of \a x+iy
    const X& imaginary_part() const { return _im; } //!< The imaginary part \a y of \a x+iy

    ModulusType<X> modulus() const { return cast_positive(sqrt(add(sqr(this->_re), sqr(this->_im)))); } //!< The modulus (absolute value) \f$r=\sqrt{x^2+y^2}\f$ of \f$x+iy\f$.
    ArgumentType<X> argument() const { return atan2(this->_re, this->_im); }                            //!< The argument (absolute value) \f$\theta=\atan(y/x)\f$ of \f$x+iy\f$.

    template<class... PRS>
    decltype(auto) get(PRS... prs) {
        typedef decltype(this->_re.get(prs...)) R;
        return Complex<R>(this->_re.get(prs...), this->_im.get(prs...));
    }
    //!@}

    //!@{
    //! \name Standard arithmetic operators
    friend Complex<X> operator+(const Complex<X>& z) { return pos(z); }                             //!< Unary plus.
    friend Complex<X> operator-(const Complex<X>& z) { return neg(z); }                             //!< Unary minus.
    friend Complex<X> operator+(const Complex<X>& z1, const Complex<X>& z2) { return add(z1, z2); } //!< Plus.
    friend Complex<X> operator-(const Complex<X>& z1, const Complex<X>& z2) { return sub(z1, z2); } //!< Minus.
    friend Complex<X> operator*(const Complex<X>& z1, const Complex<X>& z2) { return mul(z1, z2); } //!< Times.
    friend Complex<X> operator/(const Complex<X>& z1, const Complex<X>& z2) { return div(z1, z2); } //!< Divides.
    friend Complex<X>& operator+=(Complex<X>& z1, const Complex<X>& z2) {
        z1._re += z2._re;
        z1._im += z2._im;
        return z1;
    } //!< Inplace plus.
    friend Complex<X>& operator-=(Complex<X>& z1, const Complex<X>& z2) {
        z1._re -= z2._re;
        z1._im -= z2._im;
        return z1;
    } //!< Inplace minus.
    friend Complex<X>& operator*=(Complex<X>& z1, const Complex<X>& z2) { return z1 = mul(z1, z2); } //!< Inplace times.
    friend Complex<X>& operator/=(Complex<X>& z1, const Complex<X>& z2) { return z1 = div(z1, z2); } //!< Inplace divides.
    //!@}

    //!@{
    //! \name Named arithmetical functions

    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend Complex<X> add(const Complex<X>& z1, const Complex<X>& z2) {
        return Complex<X>(z1._re + z2._re, z1._im + z2._im);
    } //!< \brief Sum \a z1+z2.
    friend Complex<X> sub(const Complex<X>& z1, const Complex<X>& z2) {
        return Complex<X>(z1._re - z2._re, z1._im - z2._im);
    } //!< \brief Difference \a z1-z2.
    friend Complex<X> mul(const Complex<X>& z1, const Complex<X>& z2) {
        return Complex<X>(z1._re * z2._re - z1._im * z2._im, z1._re * z2._im + z1._im * z2._re);
    } //!< \brief Product \a z1×z2.
    friend Complex<X> div(const Complex<X>& z1, const Complex<X>& z2) {
        X ns = add(sqr(z2._re), sqr(z2._im));
        return Complex<X>((z1._re * z2._re + z1._im * z2._im) / ns, (z1._im * z2._re - z1._re * z2._im) / ns);
    } //!< \brief Quotient \a z1÷z2.

    friend Complex<X> sqrt(const Complex<X>& z) { //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
        auto r = sqrt(z.modulus());
        auto th = z.argument() / 2;
        return Complex(r * cos(th), r * sin(th));
    }
    friend Complex<X> exp(const Complex<X>& z) {
        X a = exp(z._re);
        return Complex<X>(a * cos(z._im), a * sin(z._im));
    } //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend Complex<X> log(const Complex<X>& z) {
        auto rs = sqr(z._re) + sqr(z._im);
        auto th = arg(z);
        return Complex<X>(log(rs) / 2, th);
    }
    //!< The natural logarithm of \a z. Currently requires \a Re(z) ≥ 0.
    //!@}

    //!@{
    //! \name Special complex number functions
    friend ModulusType<X> abs(const Complex<X>& z) { return z.modulus(); }
    friend decltype(auto) mag(const Complex<X>& z) { return cast_positive(sqrt(add(sqr(mag(z._re)), sqr(mag(z._im))))); } //!< Maximum possible absolute value \a |r|.
    friend decltype(auto) mig(const Complex<X>& z) { return Detail::complex_mig(z._re, z._im); }                          //!< Minimum possible absolute value \a |r|.
    friend ArgumentType<X> arg(const Complex<X>& z) { return z.argument(); }
    friend Complex<X> conj(const Complex<X>& z) { return Complex<X>(z._re, -z._im); }
    //!@}

    //!@{
    //! Operations based on the metric structure.
    friend Positive<X> dist(const Complex<X>& z1, const Complex<X>& z2) { return cast_positive(abs(sub(z1, z2))); }
    //< The distance |\a r <sub>1</sub>-\a r <sub>2</sub>| between \a r<sub>1</sub> and \a r<sub>2</sub>.
    //!@}

    //!@{
    //! \name Operations on the representation.
    friend Bool same(const Complex<X>&, const Complex<X>&); //!< Test equivalence of representation.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const Complex<X>& r) {
        return os << r._re << "+i*" << r._im;
    } //< Write to an output stream.
      //!@}

  private:
    static X _arg(const Complex<X>& z);
};

namespace Constants {
extern ARIADNE_NUMERIC_API const Complex<Integer> i;
}

} // namespace Ariadne

#endif
