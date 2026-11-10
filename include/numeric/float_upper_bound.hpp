/***************************************************************************
 *            numeric/float_upper_bound.hpp
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

/*! \file numeric/float_upper_bound.hpp
 *  \brief Floating-point upper bounds for real numbers.
 */

#ifndef ARIADNE_FLOAT_UPPER_BOUND_HPP
#define ARIADNE_FLOAT_UPPER_BOUND_HPP

#include "float.decl.hpp"
#include "float_factory.hpp"
#include "float_operations.hpp"
#include "float_traits.hpp"
#include "foundation/logical.decl.hpp"
#include "number.decl.hpp"
#include "positive.hpp"
#include "upper_number.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for real numbers.
//! \sa UpperReal, FloatDP, Bounds, LowerBound.
template<class F>
class UpperBound
    : public DefineDirectedFloatOperations<UpperBound<F>, LowerBound<F>>,
      public DefineFloatOperations<Approximation<F>> {
  protected:
    typedef UpperTag P;
    typedef typename F::RoundingModeType RND;
    typedef typename F::PrecisionType PR;

  public:
    //! <p/>
    typedef UpperTag Paradigm;
    //! <p/>
    typedef UpperBound<F> NumericType;
    //! <p/>
    typedef ValidatedUpperNumber GenericType;
    //! <p/>
    typedef F RawType;
    //! <p/>
    typedef PR PrecisionType;
    //! <p/>
    typedef PR CharacteristicsType;

  public:
    //! A upper bound of zero with precision \a pr.
    explicit UpperBound(PrecisionType pr) : _u(0.0_x, pr) {}
    //! A upper bound with value \a u.
    UpperBound(const RawType& u) : _u(u) {}

    template<BuiltinIntegral N>
    UpperBound(N n, PR pr) : UpperBound(ExactDouble(n), pr) {}
    UpperBound(const ExactDouble& d, PR pr) : _u(d, pr) {}
    UpperBound(const TwoExp& t, PR pr) : _u(t, pr) {}
    UpperBound(const Integer& z, PR pr) : _u(z, up, pr) {}
    UpperBound(const Dyadic& w, PR pr) : _u(w, up, pr) {}
    UpperBound(const Decimal& d, PR pr) : _u(d, up, pr) {}
    UpperBound(const Rational& q, PR pr) : _u(q, up, pr) {}
    UpperBound(const Real& r, PR pr);
    UpperBound(const F& x, PR pr);
    UpperBound(const Bounds<F>& x, PR pr);
    UpperBound(const UpperBound<F>& x, PR pr);
    //! A upper bound of type \p F from a generic upper bound \a y.
    UpperBound(const ValidatedUpperNumber& y, PR pr);
    template<class FF>
        requires Constructible<F, FF, RND, PR>
    UpperBound(const UpperBound<FF>& x, PR pr) : _u(x.raw(), up, pr) {}
    // Preserve the exact source value and round outward; do not bounce
    // through a generic ValidatedUpperNumber and recursively dispatch.
    template<class FF>
        requires(!Constructible<F, FF, RND, PR> && Constructible<Rational, FF> && Constructible<F, Rational, RND, PR>)
    UpperBound(const UpperBound<FF>& x, PR pr)
        : _u(Rational(x.raw()), up, pr) {}

    //! Convert from upper \em and lower bounds on a number.
    UpperBound(const Bounds<F>& x);
    template<class FE>
    UpperBound(const Ball<F, FE>& x);
    UpperBound(const Error<F>& x); // FIXME: Remove

    UpperBound<F>& operator=(const F& x) { return *this = UpperBound<F>(x); }
    UpperBound<F>& operator=(const Bounds<F>& x) { return *this = x.upper(); }
    //! Assign from the upper bound \a y, keeping the same properties.
    UpperBound<F>& operator=(const ValidatedUpperNumber& y) { return *this = UpperBound<F>(y, this->precision()); }
    //! Create a upper bound from the generic upper bound \a y with the same properties as \a this.
    UpperBound<F> create(const ValidatedUpperNumber& y) const { return UpperBound<F>(y, this->precision()); }
    //! Create a lower bound from the generic lower bound \a y with the same properties as \a this.
    LowerBound<F> create(const ValidatedLowerNumber& y) const { return LowerBound<F>(y, this->precision()); }

    //! Downcast to a generic upper bound.
    operator ValidatedUpperNumber() const;

    //! The precision of the floating-point type used.
    PrecisionType precision() const { return _u.precision(); }
    //! The compuational characteristics needed to create the upper bound; equivalent to the precision.
    CharacteristicsType characteristics() const { return _u.precision(); }
    //! Downcast to a generic upper bound.
    GenericType generic() const;
    //! The raw data used to represent the upper bound.
    const RawType& raw() const { return _u; }
    //! A mutable reference to the raw data used to represent the upper bound.
    RawType& raw() { return _u; }
    //! Over-approximate by a builtin double-precision value. DEPRECATED
    double get_d() const { return _u.get_d(); }

  public:
#ifdef DOXYGEN
    //!@{
    //! \name Arithmetic operators
    friend UpperBound<F> operator+(const UpperBound<F>& x);                           //!< <p/>
    friend UpperBound<F> operator-(const LowerBound<F>& x);                           //!< <p/>
    friend LowerBound<F> operator-(const UpperBound<F>& x);                           //!< <p/>
    friend UpperBound<F> operator+(const UpperBound<F>& x1, const UpperBound<F>& x2); //!< <p/>
    friend UpperBound<F> operator-(const UpperBound<F>& x1, const LowerBound<F>& x2); //!< <p/>
    friend LowerBound<F> operator-(const LowerBound<F>& x1, const UpperBound<F>& x2); //!< <p/>
    friend UpperBound<F>& operator+=(UpperBound<F>& x1, const UpperBound<F>& x2);     //!< <p/>
    friend UpperBound<F>& operator-=(UpperBound<F>& x1, const LowerBound<F>& x2);     //!< <p/>
    friend LowerBound<F>& operator-=(LowerBound<F>& x1, const UpperBound<F>& x2);     //!< <p/>

    friend Positive<UpperBound<F>> operator*(const Positive<UpperBound<F>>& x1, const Positive<UpperBound<F>>& x2); //!< <p/>
    friend Positive<UpperBound<F>> operator/(const Positive<UpperBound<F>>& x1, const Positive<LowerBound<F>>& x2); //!< <p/>
    friend Positive<LowerBound<F>> operator/(const Positive<LowerBound<F>>& x1, const Positive<UpperBound<F>>& x2); //!< <p/>
    friend Positive<UpperBound<F>>& operator+=(Positive<UpperBound<F>>& x1, const Positive<UpperBound<F>>& x2);     //!< <p/>
    friend Positive<UpperBound<F>>& operator*=(Positive<UpperBound<F>>& x1, const Positive<UpperBound<F>>& x2);     //!< <p/>
    friend Positive<UpperBound<F>>& operator/=(Positive<UpperBound<F>>& x1, const Positive<LowerBound<F>>& x2);     //!< <p/>
    friend Positive<LowerBound<F>>& operator/=(Positive<LowerBound<F>>& x1, const Positive<UpperBound<F>>& x2);     //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operators
    friend ValidatedLowerKleenean operator<=(const UpperBound<F>& x1, const LowerBound<F>& x2); //!< <p/>
    friend ValidatedUpperKleenean operator>=(const UpperBound<F>& x1, const LowerBound<F>& x2); //!< <p/>
    friend ValidatedLowerKleenean operator<(const UpperBound<F>& x1, const LowerBound<F>& x2);  //!< <p/>
    friend ValidatedUpperKleenean operator>(const UpperBound<F>& x1, const LowerBound<F>& x2);  //!< <p/>
    friend ValidatedUpperKleenean operator<=(const LowerBound<F>& x1, const UpperBound<F>& x2); //!< <p/>
    friend ValidatedLowerKleenean operator>=(const LowerBound<F>& x1, const UpperBound<F>& x2); //!< <p/>
    friend ValidatedUpperKleenean operator<(const LowerBound<F>& x1, const UpperBound<F>& x2);  //!< <p/>
    friend ValidatedLowerKleenean operator>(const LowerBound<F>& x1, const UpperBound<F>& x2);  //!< <p/>
    //!@}
#endif // DOXYGEN

    //!@{
    //! \name Monotone arithmetic operations
    friend UpperBound<F> nul(const UpperBound<F>& x) {
        return UpperBound<F>(nul(x._u));
    } //!< <p/>
    friend UpperBound<F> pos(const UpperBound<F>& x) {
        return UpperBound<F>(pos(x._u));
    } //!< <p/>
    friend LowerBound<F> neg(const UpperBound<F>& x) {
        return LowerBound<F>(neg(x._u));
    } //!< <p/>
    friend UpperBound<F> hlf(const UpperBound<F>& x) {
        return UpperBound<F>(hlf(x._u));
    } //!< <p/>

    friend UpperBound<F> add(const UpperBound<F>& x1, const UpperBound<F>& x2) {
        return UpperBound<F>(add(up, x1._u, x2._u));
    } //!< <p/>
    friend UpperBound<F> sub(const UpperBound<F>& x1, const LowerBound<F>& x2) {
        return UpperBound<F>(sub(up, x1._u, x2._l));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Monotone algebraic and transcendental operations
    friend Positive<UpperBound<F>> sqrt(const Positive<UpperBound<F>>& x); //!< <p/>
    friend Positive<UpperBound<F>> exp(const UpperBound<F>& x);            //!< <p/>
    friend UpperBound<F> log(const Positive<UpperBound<F>>& x);            //!< <p/>
    friend UpperBound<F> atan(const UpperBound<F>& x);                     //!< <p/>
    friend UpperBound<F> tanh(const UpperBound<F>& x) {
        return UpperBound<F>(tanh(Bounds<F>(x.raw())).upper_raw());
    }                                                                      //!< <p/>
    friend Positive<UpperBound<F>> atan(const Positive<UpperBound<F>>& x); //!< <p/>
    //!@}

    //!@{
    //! \name Lattice operations
    friend UpperBound<F> max(const UpperBound<F>& x1, const UpperBound<F>& x2) {
        return UpperBound<F>(max(x1._u, x2._u));
    } //!< <p/>
    friend UpperBound<F> min(const UpperBound<F>& x1, const UpperBound<F>& x2) {
        return UpperBound<F>(min(x1._u, x2._u));
    } //!< <p/>
    friend Approximation<F> abs(const UpperBound<F>& x) {
        return abs(Approximation<F>(x));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operations
    friend ValidatedNegatedSierpinskian eq(const UpperBound<F>& x1, const LowerBound<F>& x2) {
        if (x1._u < x2._l) {
            return false;
        } //!< <p/>
        else {
            return ValidatedNegatedSierpinskian(LogicalValue::INDETERMINATE);
        }
    } //!< <p/>
    friend ValidatedLowerKleenean lt(const UpperBound<F>& x1, const LowerBound<F>& x2) {
        if (x1._u < x2._l) {
            return true;
        } //!< <p/>
        else {
            return ValidatedLowerKleenean(LogicalValue::UNLIKELY);
        }
    } //!< <p/>
    //!@}

    //!@{
    //! \name Rounding operations
    friend Integer cast_integer(const UpperBound<F>& x) {
        return ceil(static_cast<Dyadic>(x._u));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Special value tests
    friend Bool is_nan(const UpperBound<F>& x) {
        return is_nan(x._u);
    } //!< <p/>
    //!@{

    //! \name Validated information tests and operations
    friend Bool same(const UpperBound<F>& x1, const UpperBound<F>& x2) {
        return x1._u == x2._u;
    } //!< <p/>
    friend Bool refines(const UpperBound<F>& x1, const UpperBound<F>& x2) {
        return x1._u <= x2._u;
    } //!< <p/>
    friend UpperBound<F> refinement(const UpperBound<F>& x1, const UpperBound<F>& x2) {
        return UpperBound<F>(min(x1._u, x2._u));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const UpperBound<F>& x) {
        return write(os, x.raw(), DecimalPrecision{Bounds<F>::output_places}, upward);
    } //!< Write to an output stream.
    friend InputStream& operator>>(InputStream&, UpperBound<F>&) {
        ARIADNE_NOT_IMPLEMENTED;
    } //!< Read from an input stream.
      //!@}

  public:
    friend LowerBound<F> neg(const UpperBound<F>& x);
    friend LowerBound<F> sub(const LowerBound<F>& x1, const UpperBound<F>& x2);

  public:
    friend UpperBound<F> operator*(const PositiveBounds<F>& x1, const UpperBound<F>& x2) {
        return UpperBound<F>(mul(up, x2.raw() >= 0 ? x1.upper().raw() : x1.lower().raw(), x2.raw()));
    }
    friend UpperBound<F> operator*(const UpperBound<F>& x1, const PositiveBounds<F>& x2) {
        return UpperBound<F>(mul(up, x1.raw(), x1.raw() >= 0 ? x2.upper().raw() : x2.lower().raw()));
    }
    friend UpperBound<F> operator/(const UpperBound<F>& x1, const PositiveBounds<F>& x2) {
        return UpperBound<F>(div(up, x1.raw(), x1.raw() >= 0 ? x2.lower().raw() : x2.upper().raw()));
    }
    // Needed to prevent ambiguity; useful as implementation is easier than PositiveBounds version.
    friend UpperBound<F> operator*(const Positive<F>& x1, const UpperBound<F>& x2) {
        return UpperBound<F>(mul(up, x1, x2.raw()));
    }
    friend UpperBound<F> operator*(const UpperBound<F>& x1, const Positive<F>& x2) {
        return UpperBound<F>(mul(up, x1.raw(), x2));
    }
    friend UpperBound<F> operator/(const UpperBound<F>& x1, const Positive<F>& x2) {
        return UpperBound<F>(div(up, x1.raw(), x2));
    }

  private:

  public:
    RawType _u;
};

template<class F>
template<class FE>
UpperBound<F>::UpperBound(const Ball<F, FE>& x) : UpperBound(x.upper_raw()) {}

template<class PR>
UpperBound(ValidatedUpperNumber, PR) -> UpperBound<RawFloatType<PR>>;
template<class F>
UpperBound(F) -> UpperBound<F>;

template<class F>
inline FloatFactory<PrecisionType<F>> factory(const UpperBound<F>& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR>
inline FloatUpperBound<PR> FloatFactory<PR>::create(const ValidatedUpperNumber& y) { return FloatUpperBound<PR>(y, _pr); }
template<class PR>
inline PositiveFloatUpperBound<PR> FloatFactory<PR>::create(const PositiveValidatedUpperNumber& y) { return PositiveFloatUpperBound<PR>(y, _pr); }

template<class F>
class Positive<UpperBound<F>> : public UpperBound<F>,
                                DefineConcreteGenericOperators<PositiveUpperBound<F>> {
    using typename UpperBound<F>::PR;

  public:
    Positive() : UpperBound<F>() {}
    explicit Positive(PR const& pr) : UpperBound<F>(pr) {}
    explicit Positive(const F& x) : UpperBound<F>(x) {}
    template<BuiltinUnsignedIntegral M>
    Positive(M m, PR pr) : UpperBound<F>(m, pr) {}
    template<BuiltinUnsignedIntegral M>
    Positive<F> create(M m) const { return Positive<F>(m, this->precision()); }
    explicit Positive(const UpperBound<F>& x) : UpperBound<F>(x) { ARIADNE_PRECONDITION_MSG(!(this->_u < 0), "x=" << x); }
    Positive(const PositiveValidatedUpperNumber& y, PR pr) : UpperBound<F>(y, pr) { ARIADNE_PRECONDITION_MSG(!(this->_u < 0), "y=" << y); }
    Positive(const Positive<F>& x) : UpperBound<F>(x) {}
    Positive(const PositiveBounds<F>& x) : UpperBound<F>(x) {}

  public:
    friend PositiveUpperBound<F> nul(const PositiveUpperBound<F>& x) {
        return PositiveUpperBound<F>(nul(x.raw()));
    }
    friend PositiveUpperBound<F> pos(const PositiveUpperBound<F>& x) {
        return PositiveUpperBound<F>(pos(x.raw()));
    }
    friend PositiveUpperBound<F> sqr(const PositiveUpperBound<F>& x) {
        return PositiveUpperBound<F>(sqr(up, x.raw()));
    }
    friend PositiveUpperBound<F> rec(const PositiveLowerBound<F>& x) {
        return PositiveUpperBound<F>(rec(up, x.raw()));
    }
    friend PositiveUpperBound<F> add(const PositiveUpperBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveUpperBound<F>(add(up, x1.raw(), x2.raw()));
    }
    friend PositiveUpperBound<F> mul(const PositiveUpperBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveUpperBound<F>(mul(up, x1.raw(), x2.raw()));
    }
    friend PositiveUpperBound<F> div(const PositiveUpperBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveUpperBound<F>(div(up, x1.raw(), x2.raw()));
    }
    friend PositiveUpperBound<F> max(const PositiveUpperBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveUpperBound<F>(max(x1.raw(), x2.raw()));
    }
    friend PositiveUpperBound<F> max(const PositiveUpperBound<F>& x1, const UpperBound<F>& x2) {
        return PositiveUpperBound<F>(max(x1.raw(), x2.raw()));
    }
    friend PositiveUpperBound<F> max(const UpperBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveUpperBound<F>(max(x1.raw(), x2.raw()));
    }
    friend PositiveUpperBound<F> min(const PositiveUpperBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveUpperBound<F>(min(x1.raw(), x2.raw()));
    }

    friend PositiveLowerBound<F> rec(const PositiveUpperBound<F>& x);
    friend PositiveLowerBound<F> div(const PositiveLowerBound<F>& x1, const PositiveUpperBound<F>& x2);

    friend PositiveUpperBound<F> pow(const PositiveUpperBound<F>& x, Nat m) {
        return PositiveUpperBound<F>(pow(up, x._u, static_cast<Int>(m)));
    }
    friend Approximation<F> pow(const UpperBound<F>& x, Int n) {
        return pow(Approximation<F>(x), n);
    }
    // FIXME: Implement pow for Natural/Integer
    friend PositiveUpperBound<F> pow(const PositiveUpperBound<F>& x, const Natural& m) {
        return PositiveUpperBound<F>(pow(up, x._u, static_cast<Int>(m.get_si())));
    }

    friend PositiveUpperBound<F> sqrt(const PositiveUpperBound<F>& x) {
        return PositiveUpperBound<F>(sqrt(up, x.raw()));
    }
    friend PositiveUpperBound<F> exp(const UpperBound<F>& x) {
        return PositiveUpperBound<F>(exp(up, x.raw()));
    }
    friend UpperBound<F> log(const PositiveUpperBound<F>& x) {
        return UpperBound<F>(log(up, x.raw()));
    }
    friend UpperBound<F> atan(const UpperBound<F>& x) {
        return UpperBound<F>(atan(up, x.raw()));
    }
    friend PositiveUpperBound<F> atan(const PositiveUpperBound<F>& x) {
        return PositiveUpperBound<F>(atan(up, x.raw()));
    }

  public:
    friend PositiveUpperBound<F> operator+(const PositiveUpperBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveUpperBound<F>(add(up, x1.raw(), x2.raw()));
    }
    friend PositiveUpperBound<F> operator*(const PositiveUpperBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveUpperBound<F>(mul(up, x1.raw(), x2.raw()));
    }
    friend PositiveUpperBound<F> operator/(const PositiveUpperBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveUpperBound<F>(div(up, x1.raw(), x2.raw()));
    }
    friend PositiveUpperBound<F>& operator+=(PositiveUpperBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return x1 = x1 + x2;
    }
    friend PositiveUpperBound<F>& operator*=(PositiveUpperBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return x1 = x1 * x2;
    }
    friend PositiveUpperBound<F>& operator/=(PositiveUpperBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return x1 = x1 / x2;
    }

  public:
    // Needed to prevent ambiguity
    friend PositiveUpperBound<F> operator*(const PositiveBounds<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveUpperBound<F>(mul(up, x1.upper().raw(), x2.raw()));
    }
    friend PositiveUpperBound<F> operator*(const PositiveUpperBound<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveUpperBound<F>(mul(up, x1.raw(), x2.upper().raw()));
    }
    friend PositiveUpperBound<F> operator/(const PositiveUpperBound<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveUpperBound<F>(div(up, x1.raw(), x2.lower().raw()));
    }
    friend PositiveUpperBound<F> operator*(const Positive<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveUpperBound<F>(mul(up, x1, x2.raw()));
    }
    friend PositiveUpperBound<F> operator*(const PositiveUpperBound<F>& x1, const Positive<F>& x2) {
        return PositiveUpperBound<F>(mul(up, x1.raw(), x2));
    }
    friend PositiveUpperBound<F> operator/(const PositiveUpperBound<F>& x1, const Positive<F>& x2) {
        return PositiveUpperBound<F>(div(up, x1.raw(), x2));
    }
    friend PositiveUpperBound<F> operator/(const Positive<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveUpperBound<F>(div(down, x1, x2.raw()));
    }
    friend PositiveUpperBound<F> operator/(const PositiveUpperBound<F>& x1, Nat m2) {
        return PositiveUpperBound<F>(div(up, x1.raw(), F(m2, x1.precision())));
    }
    friend PositiveLowerBound<F> operator/(const Positive<F>& x1, const PositiveUpperBound<F>& x2);
};

template<class F>
inline PositiveUpperBound<F> cast_positive(const UpperBound<F>& x) {
    return PositiveUpperBound<F>(x);
}

} // namespace Ariadne

#endif
