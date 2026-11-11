/***************************************************************************
 *            numeric/float_error.hpp
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

/** @file numeric/float_error.hpp
 *  @brief Floating-point error bounds for metric spaces.
 */

#ifndef ARIADNE_FLOAT_ERROR_HPP
#define ARIADNE_FLOAT_ERROR_HPP

#include "builtin.hpp"
#include "float.decl.hpp"
#include "float_approximation.hpp"
#include "float_factory.hpp"
#include "float_lower_bound.hpp"
#include "float_operations.hpp"
#include "float_traits.hpp"
#include "float_upper_bound.hpp"
#include "foundation/logical.decl.hpp"
#include "number.decl.hpp"
#include "positive.hpp"
#include "utility/macros.hpp"
#include "utility/metaprogramming.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Floating-point upper bounds for positive real numbers, suitable for use as an upper bound for an error in a metric space.
//! \sa FloatDP, Float, UpperBound.
template<class F>
class Error
    : public DefineDirectedFloatOperations<UpperBound<F>, LowerBound<F>>

    ,
      public DeclarePositiveDirectedNumericOperations<PositiveUpperBound<F>, PositiveLowerBound<F>>,
      public ProvideConcreteGenericDirectedSemiFieldOperations<PositiveUpperBound<F>, PositiveLowerBound<F>, Nat, Nat> {
    using PR = typename F::PrecisionType;

  private:

  public:
    F _e;

  public:
    typedef PositiveValidatedUpperNumber GenericType;
    typedef PR PrecisionType;
    typedef PR CharacteristicsType;

  public:
    /** @brief Construct a \c Error from x. */
    Error(const Positive<F>& x);
    /** @brief Construct a \c Error from x. */
    Error(const PositiveBounds<F>& x);
    //! Convert from a positive upper bound \a x on the value to be used to represent an error-bound.
    Error(const PositiveUpperBound<F>& x) : _e(x._u) {}

  public:
    /** @brief Construct a \c Error from pr. */
    explicit Error(PR const& pr) : _e(pr) {}
    //! Treat \a a as an upper-bound for an error.
    //! \pre Requires that \a e is positive (or zero).
    //! \todo Check error not being negative to allow for NaN as a valid input.
    explicit Error(const F& e) : _e(e) { ARIADNE_PRECONDITION_MSG(!(this->_e < 0), "e=" << *this); }
    //! Treat the natural number \a m as an upper-bound for an error, represented with precision \a pr.
    template<BuiltinUnsignedIntegral M>
    Error(M m, PR pr) : _e(m, pr) {}
    /** @brief Construct a \c Error from x. */
    explicit Error(const UpperBound<F>& x) : _e(x._u) {
        ARIADNE_PRECONDITION_MSG(!(this->_e < 0), "e=" << *this);
    }
    //! Treat \a y as an upper-bound for an error, represented with precision \a pr.
    //! \pre Requires that \a y is not definitely strictly negative.
    explicit Error(const ValidatedUpperNumber& y, PR pr) : Error(UpperBound<F>(y, pr)) {}
    /** @brief Construct a \c Error from d, pr. */
    explicit Error(const ExactDouble& d, PR pr) : Error(UpperBound<F>(d, pr)) {}
    /** @brief Construct a \c Error from t, pr. */
    explicit Error(const TwoExp& t, PR pr) : Error(UpperBound<F>(t, pr)) {}
    //! Assign from the natural number \a m, keeping the same precision.
    template<BuiltinUnsignedIntegral M>
    Error<F>& operator=(M m) {
        this->_e = F(m, this->precision());
        return *this;
    }
    //! Assign from the generic error bound \a y, keeping the same precision.
    Error<F>& operator=(ValidatedErrorNumber y);
    /** @brief Convert to validated error number. */
    operator ValidatedErrorNumber() const;
    /** @brief Return an enclosure enlarged by the supplied symmetric error. */
    friend Bounds<F> pm(const Error<F>& x) { return Bounds<F>(-x._e, +x._e); } //!< <p/>

  public:
    //! Downcast to a generic error bound.
    ValidatedErrorNumber generic() const;
    //! The precision of the floating-point type used.
    PrecisionType precision() const { return _e.precision(); }
    //! The compuational characteristics needed to create the error bound; equivalent to the precision.
    CharacteristicsType characteristics() const { return _e.precision(); }
    //! The raw data used to represent the error bound.
    const F& raw() const { return _e; }
    //! A mutable reference to the raw data used to represent the error bound.
    F& raw() { return _e; }

  public:
    //!@{
    //! \name Arithmetic operators
    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(const Error<F>& x1, const Error<F>& x2) { return Error<F>(add(up, x1._e, x2._e)); } //!< <p/>
    /** @brief Return the product. */
    friend Error<F> operator*(const Error<F>& x1, const Error<F>& x2) { return Error<F>(mul(up, x1._e, x2._e)); } //!< <p/>
    /** @brief Add the operand in place. */
    friend Error<F>& operator+=(Error<F>& x1, const Error<F>& x2) { return x1 = x1 + x2; }                        //!< <p/>
    /** @brief Multiply by the operand in place. */
    friend Error<F>& operator*=(Error<F>& x1, const Error<F>& x2) { return x1 = x1 * x2; }                        //!< <p/>

    /** @brief Return the sum or unary positive value. */
    friend UpperBound<F> operator+(const Error<F>& x) { return UpperBound<F>(+x._e); }                              //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend LowerBound<F> operator-(const Error<F>& x) { return LowerBound<F>(-x._e); }                              //!< <p/>
    /** @brief Return the sum or unary positive value. */
    friend UpperBound<F> operator+(const F& x1, const Error<F>& x2) { return UpperBound<F>(add(up, x1, x2._e)); }   //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend LowerBound<F> operator-(const F& x1, const Error<F>& x2) { return LowerBound<F>(sub(down, x1, x2._e)); } //!< <p/>
    //!@}

    //!@{
    //! \name Monotone arithmetic operations
    /** @brief Return the additive zero with the operand characteristics. */
    friend Error<F> nul(const Error<F>& x) { return Error<F>(0u, x.precision()); }                                                     //!< <p/>
    /** @brief Return the unary positive value. */
    friend UpperBound<F> pos(const Error<F>& x) { return UpperBound<F>(pos(x._e)); }                                                   //!< <p/>
    /** @brief Return the additive inverse. */
    friend LowerBound<F> neg(const Error<F>& x) { return LowerBound<F>(neg(x._e)); }                                                   //!< <p/>
    /** @brief Return the sum of the operands. */
    friend Error<F> add(const Error<F>& x1, const Error<F>& x2) { return Error<F>(add(up, x1._e, x2._e)); }                            //!< <p/>
    /** @brief Return the product of the operands. */
    friend Error<F> mul(const Error<F>& x1, const Error<F>& x2) { return Error<F>(mul(up, x1._e, x2._e)); }                            //!< <p/>
    /** @brief Return the square. */
    friend Error<F> sqr(const Error<F>& x) { return Error<F>(sqr(up, x._e)); }                                                         //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend Error<F> pow(const Error<F>& x, Nat m) { return Error<F>(pow(up, x._e, static_cast<Int>(m))); }                             //!< <p/>
    /** @brief Compute the product of the first two operands plus the third. */
    friend Error<F> fma(const Error<F>& x1, const Error<F>& x2, const Error<F>& x3) { return Error<F>(fma(up, x1._e, x2._e, x3._e)); } //!< <p/>
    //!@}

    //!@{
    //! \name Monotone algebraic and transcendental operations
    /** @brief Return the square root. */
    friend Error<F> sqrt(const Error<F>& x) { return Error<F>(sqrt(up, x._e)); }                                       //!< <p/>
    /** @brief Return the exponential. */
    friend Error<F> exp(const Error<F>& x) { return Error<F>(exp(up, x._e)); }                                         //!< <p/>
    /** @brief Return the natural logarithm. */
    friend UpperBound<F> log(const Error<F>& x) { return UpperBound<F>(log(up, x._e)); }                               //!< <p/>
    /** @brief Return the base-two logarithm. */
    friend UpperBound<F> log2(const Error<F>& x) { return log(x) / cast_positive(log(Bounds<F>(2u, x.precision()))); } //!< <p/>
    //!@}

    //!@{
    //! \name Lattice operations
    /** @brief Return the maximum of the operands. */
    friend Error<F> max(const Error<F>& x1, const Error<F>& x2) { return Error<F>(max(x1._e, x2._e)); } //!< <p/>
    /** @brief Return the minimum of the operands. */
    friend Error<F> min(const Error<F>& x1, const Error<F>& x2) { return Error<F>(min(x1._e, x2._e)); } //!< <p/>
    /** @brief Return the absolute value. */
    friend Error<F> abs(const Error<F>& x) { return x; }                                                //!< <p/>
    /** @brief Return an upper bound on the absolute value. */
    friend Error<F> mag(const Error<F>& x) { return x; }                                                //!< <p/>
                                                                                                        //!@}

    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(const Error<F>& x1, const Positive<F>& x2) { return Error<F>(add(up, x1._e, x2)); }
    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(const Positive<F>& x1, const Error<F>& x2) { return Error<F>(add(up, x1, x2._e)); }
    /** @brief Return the product. */
    friend Error<F> operator*(const Error<F>& x1, const Positive<F>& x2) { return Error<F>(mul(up, x1._e, x2)); }
    /** @brief Return the product. */
    friend Error<F> operator*(const Positive<F>& x1, const Error<F>& x2) { return Error<F>(mul(up, x1, x2._e)); }

    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(const Error<F>& x1, const PositiveBounds<F>& x2) { return Error<F>(add(up, x1._e, x2._u)); }
    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(const PositiveBounds<F>& x1, const Error<F>& x2) { return Error<F>(add(up, x1._u, x2._e)); }
    /** @brief Return the product. */
    friend Error<F> operator*(const Error<F>& x1, const PositiveBounds<F>& x2) { return Error<F>(mul(up, x1._e, x2._u)); }
    /** @brief Return the product. */
    friend Error<F> operator*(const PositiveBounds<F>& x1, const Error<F>& x2) { return Error<F>(mul(up, x1._u, x2._e)); }

    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(const Error<F>& x1, const PositiveUpperBound<F>& x2) { return Error<F>(add(up, x1._e, x2._u)); }
    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(const PositiveUpperBound<F>& x1, const Error<F>& x2) { return Error<F>(add(up, x1._u, x2._e)); }
    /** @brief Return the product. */
    friend Error<F> operator*(const Error<F>& x1, const PositiveUpperBound<F>& x2) { return Error<F>(mul(up, x1._e, x2._u)); }
    /** @brief Return the product. */
    friend Error<F> operator*(const PositiveUpperBound<F>& x1, const Error<F>& x2) { return Error<F>(mul(up, x1._u, x2._e)); }

    /** @brief Return the quotient. */
    friend Error<F> operator/(const Error<F>& x1, const PositiveLowerBound<F>& x2) { return Error<F>(div(up, x1._e, x2._l)); }

    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(const Error<F>& x1, Nat m2) { return Error<F>(add(up, x1._e, F(m2, up, x1.precision()))); }
    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(Nat m1, const Error<F>& x2) { return Error<F>(add(up, F(m1, up, x2.precision()), x2._e)); }
    /** @brief Return the product. */
    friend Error<F> operator*(const Error<F>& x1, Nat m2) { return Error<F>(mul(up, x1._e, F(m2, up, x1.precision()))); }
    /** @brief Return the product. */
    friend Error<F> operator*(Nat m1, const Error<F>& x2) { return Error<F>(mul(up, F(m1, up, x2.precision()), x2._e)); }
    /** @brief Return the quotient. */
    friend Error<F> operator/(const Error<F>& x1, Nat m2) { return Error<F>(div(up, x1._e, F(m2, down, x1.precision()))); }

    /** @brief Return the difference or additive inverse. */
    friend Approximation<F> operator-(const Error<F>& x1, const Error<F>& x2) { return Approximation<F>(sub(up, x1._e, x2._e)); }
    /** @brief Return the difference or additive inverse. */
    friend LowerBound<F> operator-(const PositiveBounds<F>& x1, const Error<F>& x2) { return LowerBound<F>(sub(down, x1._l, x2._e)); }
    /** @brief Return the difference or additive inverse. */
    friend UpperBound<F> operator-(const UpperBound<F>& x1, const LowerBound<F>& x2);
    /** @brief Return the difference or additive inverse. */
    friend LowerBound<F> operator-(const LowerBound<F>& x1, const UpperBound<F>& x2);

    /*
        friend PositiveUpperBound<F> operator+(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
        friend PositiveUpperBound<F> operator*(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
        friend PositiveUpperBound<F> add(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
        friend PositiveUpperBound<F> mul(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
        friend PositiveUpperBound<F> max(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
        friend PositiveUpperBound<F> min(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2);
    */
    /** @brief Add the operand in place. */
    friend Error<F>& operator+=(Error<F>& x1, const PositiveUpperBound<F>& x2) { return x1 = x1 + x2; }
    /** @brief Add the operand in place. */
    friend Error<F>& operator+=(Error<F>& x1, const PositiveBounds<F>& x2) { return x1 = x1 + x2; }
    /** @brief Add the operand in place. */
    friend Error<F>& operator+=(Error<F>& x1, const Positive<F>& x2) { return x1 = x1 + x2; }

    /** @brief Multiply by the operand in place. */
    friend Error<F>& operator*=(Error<F>& x1, const PositiveUpperBound<F>& x2) { return x1 = x1 * x2; }
    /** @brief Multiply by the operand in place. */
    friend Error<F>& operator*=(Error<F>& x1, const PositiveBounds<F>& x2) { return x1 = x1 * x2; }
    /** @brief Multiply by the operand in place. */
    friend Error<F>& operator*=(Error<F>& x1, const Positive<F>& x2) { return x1 = x1 * x2; }

  public:
    //!@{
    //! \name Validated information tests and operations
    /** @brief Test whether the operands have the same representation. */
    friend Bool same(const Error<F>& x1, const Error<F>& x2) { return x1._e == x2._e; }
    /** @brief Test whether the operands have the same representation. */
    friend Bool same(const Error<F>& x1, const Dyadic& x2) { return x1._e == x2; }
    //! <p/>
    /** @brief Test whether the first enclosure is contained in the second. */
    friend Bool refines(const Error<F>& x1, const Error<F>& x2) { return x1._e <= x2._e; }
    //! <p/>
    /** @brief Return the common refinement of the enclosures. */
    friend Error<F> refinement(const Error<F>& x1, const Error<F>& x2) { return Error<F>(min(x1._e, x2._e)); }
    //!@}

  public:
    /** @brief Convert to the positive representation. */
    friend const Error<F>& cast_positive(const Error<F>& x) { return x; }

  public:
    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const Error<F>& x) {
        return write(os, x.raw(), DecimalPrecision{Error<F>::output_places}, upward);
    } //!< Write to an output stream.
    friend InputStream& operator>>(InputStream& is, Error<F>& x) {
        F e(x.precision());
        is >> e;
        if (!is) {
            return is;
        }
        if (e < 0) {
            ARIADNE_THROW(std::runtime_error, "operator>>(InputStream&,Error<F>&)", "negative error bound " << e);
        }
        x.raw() = e;
        return is;
    } //!< Read from an input stream.
      //!@}

  public:
    inline static Nat output_places = 3;
    //! Set the number of decimal places used for the output. DEPRECATED
    static Void set_output_places(Nat p) { output_places = p; }
};

/** @brief Construct a \c Error with the specified characteristics. */
template<class PR>
Error(ValidatedUpperNumber, PR) -> Error<RawFloatType<PR>>;
/** @brief Construct a \c Error with the specified characteristics. */
template<class F>
Error(F) -> Error<F>;

} // namespace Ariadne

#endif
