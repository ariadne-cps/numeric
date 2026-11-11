/***************************************************************************
 *            numeric/float_lower_bound.hpp
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

/** @file numeric/float_lower_bound.hpp
 *  @brief Floating-point lower bounds for real numbers.
 */

#ifndef FLOAT_LOWER_BOUND_H
#define FLOAT_LOWER_BOUND_H

#include "float.decl.hpp"
#include "float_factory.hpp"
#include "float_operations.hpp"
#include "float_traits.hpp"
#include "foundation/logical.decl.hpp"
#include "number.decl.hpp"
#include "positive.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating-point lower bounds for real numbers.
//! \sa LowerReal, FloatDP, Bounds, UpperBound.
template<class F>
class LowerBound
    : public DefineDirectedFloatOperations<LowerBound<F>, UpperBound<F>>,
      public DefineFloatOperations<Approximation<F>> {
  protected:
    typedef LowerTag P;
    typedef typename F::RoundingModeType RND;
    typedef typename F::PrecisionType PR;

  public:
    //! <p/>
    typedef LowerTag Paradigm;
    //! <p/>
    typedef LowerBound<F> NumericType;
    //! <p/>
    typedef ValidatedLowerNumber GenericType;
    //! <p/>
    typedef F RawType;
    //! <p/>
    typedef PR PrecisionType;
    //! <p/>
    typedef PR CharacteristicsType;

  public:
    //! A lower bound of zero with precision \a pr.
    explicit LowerBound(PrecisionType pr) : _l(0.0_x, pr) {}
    //! A lower bound with value \a l.
    LowerBound(const RawType& l) : _l(l) {}

    /** @brief Construct a \c LowerBound from n, pr. */
    template<BuiltinIntegral N>
    LowerBound(N n, PR pr) : LowerBound(ExactDouble(n), pr) {}
    /** @brief Construct a \c LowerBound from d, pr. */
    LowerBound(const ExactDouble& d, PR pr) : _l(d, pr) {}
    /** @brief Construct a \c LowerBound from t, pr. */
    LowerBound(const TwoExp& t, PR pr) : _l(t, pr) {}
    /** @brief Construct a \c LowerBound from z, pr. */
    LowerBound(const Integer& z, PR pr) : _l(z, down, pr) {}
    /** @brief Construct a \c LowerBound from w, pr. */
    LowerBound(const Dyadic& w, PR pr) : _l(w, down, pr) {}
    /** @brief Construct a \c LowerBound from d, pr. */
    LowerBound(const Decimal& d, PR pr) : _l(d, down, pr) {}
    /** @brief Construct a \c LowerBound from q, pr. */
    LowerBound(const Rational& q, PR pr) : _l(q, down, pr) {}
    /** @brief Construct a \c LowerBound from r, pr. */
    LowerBound(const Real& r, PR pr);
    /** @brief Construct a \c LowerBound from x, pr. */
    LowerBound(const F& x, PR pr);
    /** @brief Construct a \c LowerBound from x, pr. */
    LowerBound(const Bounds<F>& x, PR pr);
    /** @brief Construct a \c LowerBound from x, pr. */
    LowerBound(const LowerBound<F>& x, PR pr);
    //! A lower bound of type \p F from a generic lower bound \a y.
    LowerBound(const ValidatedLowerNumber& y, PR pr);
    /** @brief Construct a \c LowerBound from x, pr. */
    template<class FF>
        requires Constructible<F, FF, RND, PR>
    LowerBound(const LowerBound<FF>& x, PR pr) : _l(x.raw(), down, pr) {}
    // Exact cross-backend conversion without going through a generic
    // ValidatedLowerNumber, which would recursively redispatch to itself.
    /** @brief Construct a \c LowerBound from x, pr. */
    template<class FF>
        requires(!Constructible<F, FF, RND, PR> && Constructible<Rational, FF> && Constructible<F, Rational, RND, PR>)
    LowerBound(const LowerBound<FF>& x, PR pr)
        : _l(Rational(x.raw()), down, pr) {}

    //! Convert from lower \em and upper bounds on a number.
    LowerBound(const Bounds<F>& x);
    /** @brief Construct a \c LowerBound from x. */
    template<class FE>
    LowerBound(const Ball<F, FE>& x);

    /** @brief Assign the supplied value and return this object. */
    LowerBound<F>& operator=(const F& x) { return *this = LowerBound<F>(x); }
    /** @brief Assign the supplied value and return this object. */
    LowerBound<F>& operator=(const Bounds<F>& x) { return *this = x.lower(); }
    //! Assign from the lower bound \a y, keeping the same properties.
    LowerBound<F>& operator=(const ValidatedLowerNumber& y) { return *this = LowerBound<F>(y, this->precision()); }
    //! Create a lower bound from the generic lower bound \a y with the same properties as \a this.
    LowerBound<F> create(const ValidatedLowerNumber& y) const { return LowerBound<F>(y, this->precision()); }
    //! Create a upper bound from the generic upper bound \a y with the same properties as \a this.
    UpperBound<F> create(const ValidatedUpperNumber& y) const { return UpperBound<F>(y, this->precision()); }

    //! Downcast to a generic lower bound.
    operator ValidatedLowerNumber() const;

    //! The precision of the floating-point type used.
    PrecisionType precision() const { return _l.precision(); }
    //! The compuational characteristics needed to create the lower bound; equivalent to the precision.
    CharacteristicsType characteristics() const { return _l.precision(); }
    //! Downcast to a generic lower bound.
    GenericType generic() const;
    //! The raw data used to represent the lower bound.
    const RawType& raw() const { return _l; }
    //! A mutable reference to the raw data used to represent the lower bound.
    RawType& raw() { return _l; }
    //! Under-approximate by a builtin double-precision value. DEPRECATED \deprecated
    double get_d() const { return _l.get_d(); }

  public:
#ifdef DOXYGEN
    //!@{
    //! \name Arithmetic operators
    /** @brief Return the sum or unary positive value. */
    friend LowerBound<F> operator+(const LowerBound<F>& x);                           //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend LowerBound<F> operator-(const UpperBound<F>& x);                           //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend UpperBound<F> operator-(const LowerBound<F>& x);                           //!< <p/>
    /** @brief Return the sum or unary positive value. */
    friend LowerBound<F> operator+(const LowerBound<F>& x1, const LowerBound<F>& x2); //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend LowerBound<F> operator-(const LowerBound<F>& x1, const UpperBound<F>& x2); //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend UpperBound<F> operator-(const UpperBound<F>& x1, const LowerBound<F>& x2); //!< <p/>
    /** @brief Add the operand in place. */
    friend LowerBound<F>& operator+=(LowerBound<F>& x1, const LowerBound<F>& x2);     //!< <p/>
    /** @brief Subtract the operand in place. */
    friend LowerBound<F>& operator-=(LowerBound<F>& x1, const UpperBound<F>& x2);     //!< <p/>
    /** @brief Subtract the operand in place. */
    friend UpperBound<F>& operator-=(UpperBound<F>& x1, const LowerBound<F>& x2);     //!< <p/>

    /** @brief Return the sum or unary positive value. */
    friend Positive<LowerBound<F>> operator+(const Positive<LowerBound<F>>& x1, const Positive<LowerBound<F>>& x2); //!< <p/>
    /** @brief Return the product. */
    friend Positive<LowerBound<F>> operator*(const Positive<LowerBound<F>>& x1, const Positive<LowerBound<F>>& x2); //!< <p/>
    /** @brief Return the quotient. */
    friend Positive<LowerBound<F>> operator/(const Positive<LowerBound<F>>& x1, const Positive<UpperBound<F>>& x2); //!< <p/>
    /** @brief Return the quotient. */
    friend Positive<UpperBound<F>> operator/(const Positive<UpperBound<F>>& x1, const Positive<LowerBound<F>>& x2); //!< <p/>
    /** @brief Add the operand in place. */
    friend Positive<LowerBound<F>>& operator+=(Positive<LowerBound<F>>& x1, const Positive<LowerBound<F>>& x2);     //!< <p/>
    /** @brief Multiply by the operand in place. */
    friend Positive<LowerBound<F>>& operator*=(Positive<LowerBound<F>>& x1, const Positive<LowerBound<F>>& x2);     //!< <p/>
    /** @brief Divide by the operand in place. */
    friend Positive<LowerBound<F>>& operator/=(Positive<LowerBound<F>>& x1, const Positive<UpperBound<F>>& x2);     //!< <p/>
    /** @brief Divide by the operand in place. */
    friend Positive<UpperBound<F>>& operator/=(Positive<UpperBound<F>>& x1, const Positive<LowerBound<F>>& x2);     //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operators
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedUpperKleenean operator<=(const LowerBound<F>& x1, const UpperBound<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedLowerKleenean operator>=(const LowerBound<F>& x1, const UpperBound<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedUpperKleenean operator<(const LowerBound<F>& x1, const UpperBound<F>& x2);  //!< <p/>
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedLowerKleenean operator>(const LowerBound<F>& x1, const UpperBound<F>& x2);  //!< <p/>
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedLowerKleenean operator<=(const UpperBound<F>& x1, const LowerBound<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedUpperKleenean operator>=(const UpperBound<F>& x1, const LowerBound<F>& x2); //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedLowerKleenean operator<(const UpperBound<F>& x1, const LowerBound<F>& x2);  //!< <p/>
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedUpperKleenean operator>(const UpperBound<F>& x1, const LowerBound<F>& x2);  //!< <p/>
    //!@}
#endif // DOXYGEN

    //!@{
    //! \name Monotone arithmetic operations
    /** @brief Return the additive zero with the operand characteristics. */
    friend LowerBound<F> nul(const LowerBound<F>& x) {
        return LowerBound<F>(nul(x._l));
    } //!< <p/>
    /** @brief Return the unary positive value. */
    friend LowerBound<F> pos(const LowerBound<F>& x) {
        return LowerBound<F>(pos(x._l));
    } //!< <p/>
    /** @brief Return the additive inverse. */
    friend UpperBound<F> neg(const LowerBound<F>& x) {
        return UpperBound<F>(neg(x._l));
    } //!< <p/>
    /** @brief Divide the value by two. */
    friend LowerBound<F> hlf(const LowerBound<F>& x) {
        return LowerBound<F>(hlf(x._l));
    } //!< <p/>

    /** @brief Return the sum of the operands. */
    friend LowerBound<F> add(const LowerBound<F>& x1, const LowerBound<F>& x2) {
        return LowerBound<F>(add(down, x1._l, x2._l));
    } //!< <p/>
      //    friend Approximation<F> sub(LowerBound<F> const& x1, LowerBound<F> const& x2) {
      //        return UpperBound<F>(sub(near,x1._l,x2._l)); } //!< <p/>
    /** @brief Return the difference of the operands. */
    friend LowerBound<F> sub(const LowerBound<F>& x1, const UpperBound<F>& x2) {
        return LowerBound<F>(sub(down, x1._l, x2._u));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Monotone algebraic and transcendental operations
    /** @brief Return the square root. */
    friend LowerBound<F> sqrt(const LowerBound<F>& x) {
        return LowerBound<F>(sqrt(down, x.raw()));
    } //!< <p/>
    /** @brief Return the exponential. */
    friend LowerBound<F> exp(const LowerBound<F>& x) {
        return LowerBound<F>(exp(down, x.raw()));
    } //!< <p/>
    /** @brief Return the natural logarithm. */
    friend LowerBound<F> log(const LowerBound<F>& x) {
        return LowerBound<F>(log(down, x.raw()));
    } //!< <p/>
    /** @brief Return the inverse tangent. */
    friend LowerBound<F> atan(const LowerBound<F>& x) {
        return LowerBound<F>(atan(down, x.raw()));
    } //!< <p/>
    /** @brief Return the hyperbolic tangent. */
    friend LowerBound<F> tanh(const LowerBound<F>& x) {
        return LowerBound<F>(tanh(Bounds<F>(x.raw())).lower_raw());
    } //!< <p/>
    //!@}

    //!@{
    //! \name Lattice operations
    /** @brief Return the absolute value. */
    friend Approximation<F> abs(const LowerBound<F>& x) {
        return abs(Approximation<F>(x));
    } //!< <p/>
    /** @brief Return the maximum of the operands. */
    friend LowerBound<F> max(const LowerBound<F>& x1, const LowerBound<F>& x2) {
        return LowerBound<F>(max(x1._l, x2._l));
    } //!< <p/>
    /** @brief Return the minimum of the operands. */
    friend LowerBound<F> min(const LowerBound<F>& x1, const LowerBound<F>& x2) {
        return LowerBound<F>(min(x1._l, x2._l));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operations
    /** @brief Test equality of the operands. */
    friend ValidatedNegatedSierpinskian eq(const LowerBound<F>& x1, const UpperBound<F>& x2) {
        if (x1._l > x2._u) {
            return false;
        } //!< <p/>
        else {
            return ValidatedNegatedSierpinskian(LogicalValue::INDETERMINATE);
        }
    } //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedUpperKleenean lt(const LowerBound<F>& x1, const UpperBound<F>& x2) {
        if (x1._l >= x2._u) {
            return false;
        } //!< <p/>
        else {
            return ValidatedUpperKleenean(LogicalValue::LIKELY);
        }
    } //!< <p/>
    //!@}

    //!@{
    //! \name Rounding operations
    /** @brief Convert to an integer using the selected integer representation. */
    friend Integer cast_integer(const LowerBound<F>& x) {
        return floor(static_cast<Dyadic>(x._l));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Special value tests
    /** @brief Test whether the value is not a number. */
    friend Bool is_nan(const LowerBound<F>& x) {
        return is_nan(x._l);
    } //!< <p/>
    //!@}

    //!@{
    //! \name Validated information tests and operations
    /** @brief Test whether the operands have the same representation. */
    friend Bool same(const LowerBound<F>& x1, const LowerBound<F>& x2) {
        return x1._l == x2._l;
    } //!< <p/>
    /** @brief Test whether the first enclosure is contained in the second. */
    friend Bool refines(const LowerBound<F>& x1, const LowerBound<F>& x2) {
        return x1._l >= x2._l;
    } //!< <p/>
    /** @brief Return the common refinement of the enclosures. */
    friend LowerBound<F> refinement(const LowerBound<F>& x1, const LowerBound<F>& x2) {
        return LowerBound<F>(max(x1._l, x2._l));
    } //!< <p/>
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const LowerBound<F>& x) {
        return write(os, x.raw(), DecimalPrecision{Bounds<F>::output_places}, downward);
    } //!< Write to an output stream.
    friend InputStream& operator>>(InputStream&, LowerBound<F>&) {
        ARIADNE_NOT_IMPLEMENTED;
    } //!< Read from an input stream.
      //!@}

  public:
    /** @brief Return the additive inverse. */
    friend UpperBound<F> neg(const LowerBound<F>& x);
    /** @brief Return the difference of the operands. */
    friend UpperBound<F> sub(const UpperBound<F>& x1, const LowerBound<F>& x2);

  public:
    /** @brief Return the product. */
    friend LowerBound<F> operator*(const PositiveBounds<F>& x1, const LowerBound<F>& x2) {
        return LowerBound<F>(mul(down, x2.raw() >= 0 ? x1.lower().raw() : x1.upper().raw(), x2.raw()));
    }
    /** @brief Return the product. */
    friend LowerBound<F> operator*(const LowerBound<F>& x1, const PositiveBounds<F>& x2) {
        return LowerBound<F>(mul(down, x1.raw(), x1.raw() >= 0 ? x2.lower().raw() : x2.upper().raw()));
    }
    /** @brief Return the quotient. */
    friend LowerBound<F> operator/(const LowerBound<F>& x1, const PositiveBounds<F>& x2) {
        return LowerBound<F>(div(down, x1.raw(), x1.raw() >= 0 ? x2.upper().raw() : x2.lower().raw()));
    }
    // Needed to prevent ambiguity; useful as implementation is easier than PositiveBounds version.
    /** @brief Return the product. */
    friend LowerBound<F> operator*(const Positive<F>& x1, const LowerBound<F>& x2) {
        return LowerBound<F>(mul(down, x1, x2.raw()));
    }
    /** @brief Return the product. */
    friend LowerBound<F> operator*(const LowerBound<F>& x1, const Positive<F>& x2) {
        return LowerBound<F>(mul(down, x1.raw(), x2));
    }
    /** @brief Return the quotient. */
    friend LowerBound<F> operator/(const LowerBound<F>& x1, const Positive<F>& x2) {
        return LowerBound<F>(div(down, x1.raw(), x2));
    }

  private:

  public:
    RawType _l;
};

template<class F>
template<class FE>
LowerBound<F>::LowerBound(const Ball<F, FE>& x) : LowerBound(x.lower_raw()) {}

/** @brief Construct a \c LowerBound with the specified characteristics. */
template<class PR>
LowerBound(ValidatedLowerNumber, PR) -> LowerBound<RawFloatType<PR>>;
/** @brief Construct a \c LowerBound with the specified characteristics. */
template<class F>
LowerBound(F) -> LowerBound<F>;

/** @brief Return a factory for compatible objects. */
template<class F>
inline FloatFactory<PrecisionType<F>> factory(const LowerBound<F>& flt) { return FloatFactory<PrecisionType<F>>(flt.precision()); }
template<class PR>
inline FloatLowerBound<PR> FloatFactory<PR>::create(const ValidatedLowerNumber& y) { return FloatLowerBound<PR>(y, _pr); }
template<class PR>
inline PositiveFloatLowerBound<PR> FloatFactory<PR>::create(const PositiveValidatedLowerNumber& y) { return PositiveFloatLowerBound<PR>(y, _pr); }

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<class F>
class Positive<LowerBound<F>> : public LowerBound<F>,
                                DefineConcreteGenericOperators<PositiveLowerBound<F>> {
    using typename LowerBound<F>::PR;

  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : LowerBound<F>() {}
    /** @brief Construct a \c Positive from m, pr. */
    template<BuiltinUnsignedIntegral M>
    Positive(M m, PR const& pr) : LowerBound<F>(m, pr) {}
    /** @brief Construct a \c Positive from pr. */
    explicit Positive(PR const& pr) : LowerBound<F>(pr) {}
    /** @brief Construct a \c Positive from x. */
    explicit Positive(const F& x) : LowerBound<F>(x) {}
    /** @brief Construct a \c Positive from x. */
    explicit Positive(const LowerBound<F>& x) : LowerBound<F>(x) {}
    /** @brief Construct a \c Positive from y, pr. */
    Positive(const PositiveValidatedLowerNumber& y, PR pr) : LowerBound<F>(y, pr) {}
    /** @brief Copy-construct the object from the supplied instance. */
    Positive(const Positive<F>& x) : LowerBound<F>(x) {}
    /** @brief Copy-construct the object from the supplied instance. */
    Positive(const PositiveBounds<F>& x) : LowerBound<F>(x) {}

  public:
    /** @brief Return the additive zero with the operand characteristics. */
    friend PositiveLowerBound<F> nul(const PositiveLowerBound<F>& x) {
        return PositiveLowerBound<F>(nul(x.raw()));
    }
    /** @brief Return the unary positive value. */
    friend PositiveLowerBound<F> pos(const PositiveLowerBound<F>& x) {
        return PositiveLowerBound<F>(pos(x.raw()));
    }
    /** @brief Return the square. */
    friend PositiveLowerBound<F> sqr(const PositiveLowerBound<F>& x) {
        return PositiveLowerBound<F>(sqr(down, x.raw()));
    }
    /** @brief Return the multiplicative inverse. */
    friend PositiveLowerBound<F> rec(const PositiveUpperBound<F>& x) {
        return PositiveLowerBound<F>(rec(down, x.raw()));
    }
    /** @brief Return the sum of the operands. */
    friend PositiveLowerBound<F> add(const PositiveLowerBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveLowerBound<F>(add(down, x1.raw(), x2.raw()));
    }
    /** @brief Return the product of the operands. */
    friend PositiveLowerBound<F> mul(const PositiveLowerBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveLowerBound<F>(mul(down, x1.raw(), x2.raw()));
    }
    /** @brief Return the quotient of the operands. */
    friend PositiveLowerBound<F> div(const PositiveLowerBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveLowerBound<F>(div(down, x1.raw(), x2.raw()));
    }
    /** @brief Return the maximum of the operands. */
    friend PositiveLowerBound<F> max(const PositiveLowerBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveLowerBound<F>(max(x1.raw(), x2.raw()));
    }
    /** @brief Return the maximum of the operands. */
    friend PositiveLowerBound<F> max(const PositiveLowerBound<F>& x1, const LowerBound<F>& x2) {
        return PositiveLowerBound<F>(max(x1.raw(), x2.raw()));
    }
    /** @brief Return the maximum of the operands. */
    friend PositiveLowerBound<F> max(const LowerBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveLowerBound<F>(max(x1.raw(), x2.raw()));
    }
    /** @brief Return the minimum of the operands. */
    friend PositiveLowerBound<F> min(const PositiveLowerBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveLowerBound<F>(min(x1.raw(), x2.raw()));
    }

    /** @brief Return the multiplicative inverse. */
    friend PositiveUpperBound<F> rec(const PositiveLowerBound<F>& x);
    /** @brief Return the quotient of the operands. */
    friend PositiveUpperBound<F> div(const PositiveUpperBound<F>& x1, const PositiveLowerBound<F>& x2);

    /** @brief Raise the value to the specified power. */
    friend PositiveLowerBound<F> pow(const PositiveLowerBound<F>& x, Nat m) {
        return PositiveLowerBound<F>(pow(down, x.raw(), static_cast<int>(m)));
    }
    /** @brief Raise the value to the specified power. */
    friend Approximation<F> pow(const LowerBound<F>& x, Int n) {
        return pow(Approximation<F>(x), n);
    }

  public:
    /** @brief Return the sum or unary positive value. */
    friend PositiveLowerBound<F> operator+(const PositiveLowerBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveLowerBound<F>(add(down, x1.raw(), x2.raw()));
    }
    /** @brief Return the product. */
    friend PositiveLowerBound<F> operator*(const PositiveLowerBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveLowerBound<F>(mul(down, x1.raw(), x2.raw()));
    }
    /** @brief Return the quotient. */
    friend PositiveLowerBound<F> operator/(const PositiveLowerBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveLowerBound<F>(mul(down, x1.raw(), x2.raw()));
    }
    /** @brief Add the operand in place. */
    friend PositiveLowerBound<F>& operator+=(PositiveLowerBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return x1 = x1 + x2;
    }
    /** @brief Multiply by the operand in place. */
    friend PositiveLowerBound<F>& operator*=(PositiveLowerBound<F>& x1, const PositiveLowerBound<F>& x2) {
        return x1 = x1 * x2;
    }
    /** @brief Divide by the operand in place. */
    friend PositiveLowerBound<F>& operator/=(PositiveLowerBound<F>& x1, const PositiveUpperBound<F>& x2) {
        return x1 = x1 / x2;
    }

  public:
    // Needed to prevent ambiguity
    /** @brief Return the product. */
    friend PositiveLowerBound<F> operator*(const PositiveBounds<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveLowerBound<F>(mul(down, x1.lower().raw(), x2.raw()));
    }
    /** @brief Return the product. */
    friend PositiveLowerBound<F> operator*(const PositiveLowerBound<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveLowerBound<F>(mul(down, x1.raw(), x2.lower().raw()));
    }
    /** @brief Return the quotient. */
    friend PositiveLowerBound<F> operator/(const PositiveLowerBound<F>& x1, const PositiveBounds<F>& x2) {
        return PositiveLowerBound<F>(div(down, x1.raw(), x2.upper().raw()));
    }
    /** @brief Return the product. */
    friend PositiveLowerBound<F> operator*(const Positive<F>& x1, const PositiveLowerBound<F>& x2) {
        return PositiveLowerBound<F>(mul(down, x1, x2.raw()));
    }
    /** @brief Return the product. */
    friend PositiveLowerBound<F> operator*(const PositiveLowerBound<F>& x1, const Positive<F>& x2) {
        return PositiveLowerBound<F>(mul(down, x1.raw(), x2));
    }
    /** @brief Return the quotient. */
    friend PositiveLowerBound<F> operator/(const PositiveLowerBound<F>& x1, const Positive<F>& x2) {
        return PositiveLowerBound<F>(div(down, x1.raw(), x2));
    }
    /** @brief Return the quotient. */
    friend PositiveLowerBound<F> operator/(const Positive<F>& x1, const PositiveUpperBound<F>& x2) {
        return PositiveLowerBound<F>(div(down, x1, x2.raw()));
    }
    /** @brief Return the quotient. */
    friend PositiveLowerBound<F> operator/(const PositiveLowerBound<F>& x1, Nat m2) {
        return PositiveLowerBound<F>(div(down, x1.raw(), F(m2, x1.precision())));
    }
    /** @brief Return the quotient. */
    friend PositiveUpperBound<F> operator/(const Positive<F>& x1, const PositiveLowerBound<F>& x2);
};

/** @brief Convert to the positive representation. */
template<class F>
inline PositiveLowerBound<F> cast_positive(const LowerBound<F>& x) {
    return PositiveLowerBound<F>(x);
}

} // namespace Ariadne

#endif
