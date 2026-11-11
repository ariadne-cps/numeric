/***************************************************************************
 *            numeric/rounded_float.hpp
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

/** @file numeric/rounded_float.hpp
 *  @brief Floating-point numbers supporting rounding-mode changes and rounded arithmetic.
 */

#ifndef ARIADNE_ROUNDED_FLOAT_HPP
#define ARIADNE_ROUNDED_FLOAT_HPP

#include <algorithm> // For std::max, std::min
#include <cmath>     // For std::floor std::ceil etc
#include <iosfwd>    // For std::floor std::ceil etc
#include <limits>    // For std::numeric_limits<double>

#include "numeric/builtin.hpp"
#include "numeric/declarations.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/operators.hpp"
#include "numeric/rounding.hpp"
#include "numeric/sign.hpp"

namespace Ariadne {

class Rational;
enum class Comparison : ComparableEnumerationType;

//! \ingroup NumericModule
//! \brief The precision of a FloatDP object. Since this is fixed, the class is only a tag; all objects are equal.
//! \relates FloatDP
class DoublePrecision;
using DP = DoublePrecision;

template<class FLT>
class Rounded;

using RoundedFloatDP = Rounded<FloatDP>;
using RoundedFloatMP = Rounded<FloatMP>;

/** @brief Numeric wrapper using the selected rounding policy. */
template<>
class Rounded<FloatDP> {
    FloatDP _flt;

  public:
    typedef RawTag Paradigm;
    typedef FloatDP FloatType;
    typedef Rounded<FloatDP> NumericType;
    typedef DoublePrecision PrecisionType;
    typedef DoublePrecision CharacteristicsType;
    typedef BuiltinRoundingModeType RoundingModeType;

  public:
    /** @brief Return the rounding mode. */
    static RoundingModeType get_rounding_mode() { return FloatDP::get_rounding_mode(); }
    /** @brief Set the rounding mode. */
    static Void set_rounding_mode(RoundingModeType rnd) { FloatDP::set_rounding_mode(rnd); }
    /** @brief Set the rounding downward. */
    static Void set_rounding_downward() { FloatDP::set_rounding_downward(); }
    /** @brief Set the rounding upward. */
    static Void set_rounding_upward() { FloatDP::set_rounding_upward(); }
    /** @brief Set the rounding to nearest. */
    static Void set_rounding_to_nearest() { FloatDP::set_rounding_to_nearest(); }
    /** @brief Set the rounding toward zero. */
    static Void set_rounding_toward_zero() { FloatDP::set_rounding_toward_zero(); }

  private:
    /** @brief Construct a \c Rounded from d. */
    explicit Rounded(double d) : _flt(d) {}

  public:
    /** @brief Construct a \c Rounded with the specified characteristics. */
    explicit Rounded(PrecisionType) : _flt() {}
    /** @brief Construct a \c Rounded from x. */
    Rounded(FloatType x) : _flt(x.dbl) {}

    /** @brief Construct a \c Rounded from n. */
    template<BuiltinIntegral N>
    Rounded(N n, PrecisionType) : _flt(n) {}
    /** @brief Construct a \c Rounded from d. */
    Rounded(ExactDouble d, PrecisionType) : _flt(d.get_d()) {}
    /** @brief Construct a \c Rounded from w, pr. */
    Rounded(const Dyadic& w, PrecisionType pr) : Rounded(FloatType(w, pr)) {}
    /** @brief Construct a \c Rounded from x. */
    Rounded(FloatType x, PrecisionType) : _flt(x.dbl) {}
    /** @brief Construct a \c Rounded from x. */
    Rounded(Rounded<FloatType> x, PrecisionType) : _flt(x._flt) {}
    /** @brief Construct a \c Rounded from y, pr. */
    template<class Y>
        requires Constructible<FloatType, Y, RoundingModeType, PrecisionType> and (not BuiltinIntegral<Y>)
    Rounded(const Y& y, PrecisionType pr) : Rounded(FloatType(y, FloatType::get_rounding_mode(), pr)) {}

    /** @brief Assign the supplied value and return this object. */
    template<BuiltinIntegral N>
    Rounded<FloatType>& operator=(N n) {
        this->_flt = n;
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    Rounded<FloatType>& operator=(FloatType x) {
        this->_flt = x;
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    Rounded<FloatType>& operator=(ExactDouble x) {
        this->_flt = x;
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    template<class Y>
        requires Constructible<FloatType, Y, RoundingModeType, PrecisionType> and (not BuiltinIntegral<Y>)
    Rounded<FloatType> operator=(const Y& y) { return (*this) = FloatType(y, FloatType::get_rounding_mode(), this->precision()); }

    // Approximation<FloatDP> contains a FloatDP member; it is not a FloatDP base subobject.
    // Reinterpreting the wrapper as FloatDP is not guaranteed by the C++ object model.
    /** @brief Construct a \c Rounded from x. */
    Rounded(const Approximation<FloatDP>& x);
    /** @brief Convert to approximation< float dp >. */
    operator Approximation<FloatDP>() const;

    /** @brief Return the underlying data pointer. */
    double data() const { return this->_flt.dbl; }
    /** @brief Return the underlying numeric representation. */
    const FloatType& raw() const { return this->_flt; }
    /** @brief Convert to float type. */
    explicit operator FloatType() const { return this->_flt; }
    /** @brief Return the working precision. */
    PrecisionType precision() const { return DoublePrecision(); }
    /** @brief Return the characteristics required to construct compatible values. */
    CharacteristicsType characteristics() const { return DoublePrecision(); }

    // Non-finiteness tests
    /** @brief Test whether the value is not a number. */
    friend Bool is_nan(Rounded<FloatType> x) { return std::isnan(x._flt.dbl); }

    // Integer conversions
    /** @brief Return the greatest integer not exceeding the value. */
    friend Rounded<FloatDP> floor(Rounded<FloatDP> x) { return Rounded<FloatDP>(std::floor(x._flt.dbl)); }
    /** @brief Return the least integer not less than the value. */
    friend Rounded<FloatDP> ceil(Rounded<FloatDP> x) { return Rounded<FloatDP>(std::ceil(x._flt.dbl)); }
    /** @brief Round the value to an integer. */
    friend Rounded<FloatDP> round(Rounded<FloatDP> x) { return Rounded<FloatDP>(std::round(x._flt.dbl)); }

    // Correctly rounded arithmetic
    /** @brief Return the additive zero with the operand characteristics. */
    friend Rounded<FloatDP> nul(Rounded<FloatDP> x) { return Rounded<FloatDP>(nul_rnd(x._flt.dbl)); }
    /** @brief Return the unary positive value. */
    friend Rounded<FloatDP> pos(Rounded<FloatDP> x) { return Rounded<FloatDP>(pos_rnd(x._flt.dbl)); }
    /** @brief Return the additive inverse. */
    friend Rounded<FloatDP> neg(Rounded<FloatDP> x) { return Rounded<FloatDP>(neg_rnd(x._flt.dbl)); }
    /** @brief Return the square. */
    friend Rounded<FloatDP> sqr(Rounded<FloatDP> x) { return Rounded<FloatDP>(sqr_rnd(x._flt.dbl)); }
    /** @brief Divide the value by two. */
    friend Rounded<FloatDP> hlf(Rounded<FloatDP> x) { return Rounded<FloatDP>(hlf_rnd(x._flt.dbl)); }
    /** @brief Return the multiplicative inverse. */
    friend Rounded<FloatDP> rec(Rounded<FloatDP> x) { return Rounded<FloatDP>(rec_rnd(x._flt.dbl)); }
    /** @brief Return the sum of the operands. */
    friend Rounded<FloatDP> add(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(add_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Return the difference of the operands. */
    friend Rounded<FloatDP> sub(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(sub_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Return the product of the operands. */
    friend Rounded<FloatDP> mul(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(mul_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Return the quotient of the operands. */
    friend Rounded<FloatDP> div(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(div_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Compute the product of the first two operands plus the third. */
    friend Rounded<FloatDP> fma(Rounded<FloatDP> x1, Rounded<FloatDP> x2, Rounded<FloatDP> x3) { return Rounded<FloatDP>(fma_rnd(x1._flt.dbl, x2._flt.dbl, x3._flt.dbl)); }
    /** @brief Raise the value to the specified power. */
    friend Rounded<FloatDP> pow(Rounded<FloatDP> x, Int n) { return Rounded<FloatDP>(pow_rnd(x._flt.dbl, n)); }
    /** @brief Return the square root. */
    friend Rounded<FloatDP> sqrt(Rounded<FloatDP> x) { return Rounded<FloatDP>(sqrt_rnd(x._flt.dbl)); }
    /** @brief Return the exponential. */
    friend Rounded<FloatDP> exp(Rounded<FloatDP> x) { return Rounded<FloatDP>(exp(FloatDP::get_rounding_mode(), x._flt)); }
    /** @brief Return the natural logarithm. */
    friend Rounded<FloatDP> log(Rounded<FloatDP> x) { return Rounded<FloatDP>(log(FloatDP::get_rounding_mode(), x._flt)); }
    /** @brief Return the sine. */
    friend Rounded<FloatDP> sin(Rounded<FloatDP> x) { return Rounded<FloatDP>(sin(FloatDP::get_rounding_mode(), x._flt)); }
    /** @brief Return the cosine. */
    friend Rounded<FloatDP> cos(Rounded<FloatDP> x) { return Rounded<FloatDP>(cos(FloatDP::get_rounding_mode(), x._flt)); }
    /** @brief Return the tangent. */
    friend Rounded<FloatDP> tan(Rounded<FloatDP> x) { return Rounded<FloatDP>(tan(FloatDP::get_rounding_mode(), x._flt)); }
    /** @brief Return the hyperbolic tangent. */
    friend Rounded<FloatDP> tanh(Rounded<FloatDP> x) { return Rounded<FloatDP>(tanh(FloatDP::get_rounding_mode(), x._flt)); }
    /** @brief Return the inverse sine. */
    friend Rounded<FloatDP> asin(Rounded<FloatDP> x) { return Rounded<FloatDP>(asin(FloatDP::get_rounding_mode(), x._flt)); }
    /** @brief Return the inverse cosine. */
    friend Rounded<FloatDP> acos(Rounded<FloatDP> x) { return Rounded<FloatDP>(acos(FloatDP::get_rounding_mode(), x._flt)); }
    /** @brief Return the inverse tangent. */
    friend Rounded<FloatDP> atan(Rounded<FloatDP> x) { return Rounded<FloatDP>(atan(FloatDP::get_rounding_mode(), x._flt)); }
    /** @brief Return pi at the requested precision. */
    static Rounded<FloatDP> pi(PrecisionType) { return Rounded<FloatDP>(pi_rnd()); }

    /** @brief Return the absolute value. */
    friend Rounded<FloatDP> abs(Rounded<FloatDP> x) { return Rounded<FloatDP>(abs_rnd(x._flt.dbl)); }
    /** @brief Return the maximum of the operands. */
    friend Rounded<FloatDP> max(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(max_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Return the minimum of the operands. */
    friend Rounded<FloatDP> min(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(min_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Return an upper bound on the absolute value. */
    friend Rounded<FloatDP> mag(Rounded<FloatDP> x) { return Rounded<FloatDP>(abs_rnd(x._flt.dbl)); }
    /** @brief Return a lower bound on the absolute value. */
    friend Rounded<FloatDP> mig(Rounded<FloatDP> x) { return Rounded<FloatDP>(abs_rnd(x._flt.dbl)); }

    /** @brief Return the midpoint value. */
    friend Rounded<FloatDP> med(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(hlf_rnd(add_rnd(x1._flt.dbl, x2._flt.dbl))); }
    /** @brief Return the radius. */
    friend Rounded<FloatDP> rad(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(hlf_rnd(sub_rnd(x2._flt.dbl, x1._flt.dbl))); }

    /** @brief Return the sum or unary positive value. */
    friend Rounded<FloatDP> operator+(Rounded<FloatDP> x) { return Rounded<FloatDP>(pos_rnd(x._flt.dbl)); }
    /** @brief Return the difference or additive inverse. */
    friend Rounded<FloatDP> operator-(Rounded<FloatDP> x) { return Rounded<FloatDP>(neg_rnd(x._flt.dbl)); }
    /** @brief Return the sum or unary positive value. */
    friend Rounded<FloatDP> operator+(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(add_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Return the difference or additive inverse. */
    friend Rounded<FloatDP> operator-(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(sub_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Return the product. */
    friend Rounded<FloatDP> operator*(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(mul_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Return the quotient. */
    friend Rounded<FloatDP> operator/(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(div_rnd(x1._flt.dbl, x2._flt.dbl)); }
    /** @brief Add the operand in place. */
    friend Rounded<FloatDP>& operator+=(Rounded<FloatDP>& x1, Rounded<FloatDP> x2) {
        x1._flt.dbl = add_rnd(x1._flt.dbl, x2._flt.dbl);
        return x1;
    }
    /** @brief Subtract the operand in place. */
    friend Rounded<FloatDP>& operator-=(Rounded<FloatDP>& x1, Rounded<FloatDP> x2) {
        x1._flt.dbl = sub_rnd(x1._flt.dbl, x2._flt.dbl);
        return x1;
    }
    /** @brief Multiply by the operand in place. */
    friend Rounded<FloatDP>& operator*=(Rounded<FloatDP>& x1, Rounded<FloatDP> x2) {
        x1._flt.dbl = mul_rnd(x1._flt.dbl, x2._flt.dbl);
        return x1;
    }
    /** @brief Divide by the operand in place. */
    friend Rounded<FloatDP>& operator/=(Rounded<FloatDP>& x1, Rounded<FloatDP> x2) {
        x1._flt.dbl = div_rnd(x1._flt.dbl, x2._flt.dbl);
        return x1;
    }

    /** @brief Return the sum or unary positive value. */
    friend Rounded<FloatDP> operator+(Rounded<FloatDP> x1, ExactDouble x2) { return Rounded<FloatDP>(add_rnd(x1._flt.dbl, x2.get_d())); }
    /** @brief Return the difference or additive inverse. */
    friend Rounded<FloatDP> operator-(Rounded<FloatDP> x1, ExactDouble x2) { return Rounded<FloatDP>(sub_rnd(x1._flt.dbl, x2.get_d())); }
    /** @brief Return the product. */
    friend Rounded<FloatDP> operator*(Rounded<FloatDP> x1, ExactDouble x2) { return Rounded<FloatDP>(mul_rnd(x1._flt.dbl, x2.get_d())); }
    /** @brief Return the quotient. */
    friend Rounded<FloatDP> operator/(Rounded<FloatDP> x1, ExactDouble x2) { return Rounded<FloatDP>(div_rnd(x1._flt.dbl, x2.get_d())); }
    /** @brief Return the sum or unary positive value. */
    friend Rounded<FloatDP> operator+(ExactDouble x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(add_rnd(x1.get_d(), x2._flt.dbl)); }
    /** @brief Return the difference or additive inverse. */
    friend Rounded<FloatDP> operator-(ExactDouble x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(sub_rnd(x1.get_d(), x2._flt.dbl)); }
    /** @brief Return the product. */
    friend Rounded<FloatDP> operator*(ExactDouble x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(mul_rnd(x1.get_d(), x2._flt.dbl)); }
    /** @brief Return the quotient. */
    friend Rounded<FloatDP> operator/(ExactDouble x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(div_rnd(x1.get_d(), x2._flt.dbl)); }
    /** @brief Add the operand in place. */
    friend Rounded<FloatDP>& operator+=(Rounded<FloatDP>& x1, ExactDouble x2) { return x1 = x1 + x2; }
    /** @brief Subtract the operand in place. */
    friend Rounded<FloatDP>& operator-=(Rounded<FloatDP>& x1, ExactDouble x2) { return x1 = x1 - x2; }
    /** @brief Multiply by the operand in place. */
    friend Rounded<FloatDP>& operator*=(Rounded<FloatDP>& x1, ExactDouble x2) { return x1 = x1 * x2; }
    /** @brief Divide by the operand in place. */
    friend Rounded<FloatDP>& operator/=(Rounded<FloatDP>& x1, ExactDouble x2) { return x1 = x1 / x2; }

    /** @brief Return the maximum of the operands. */
    friend Rounded<FloatDP> max(Rounded<FloatDP> x1, ExactDouble x2) { return Rounded<FloatDP>(max_rnd(x1._flt.dbl, x2.get_d())); }
    /** @brief Return the minimum of the operands. */
    friend Rounded<FloatDP> min(Rounded<FloatDP> x1, ExactDouble x2) { return Rounded<FloatDP>(min_rnd(x1._flt.dbl, x2.get_d())); }
    /** @brief Return the maximum of the operands. */
    friend Rounded<FloatDP> max(ExactDouble x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(max_rnd(x1.get_d(), x2._flt.dbl)); }
    /** @brief Return the minimum of the operands. */
    friend Rounded<FloatDP> min(ExactDouble x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(min_rnd(x1.get_d(), x2._flt.dbl)); }

    /** @brief Test equality of the operands. */
    friend ApproximateKleenean operator==(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1._flt.dbl == x2._flt.dbl; }
    /** @brief Test inequality of the operands. */
    friend ApproximateKleenean operator!=(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1._flt.dbl != x2._flt.dbl; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ApproximateKleenean operator<=(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1._flt.dbl <= x2._flt.dbl; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ApproximateKleenean operator>=(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1._flt.dbl >= x2._flt.dbl; }
    /** @brief Test whether the first operand is less than the second. */
    friend ApproximateKleenean operator<(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1._flt.dbl < x2._flt.dbl; }
    /** @brief Test whether the first operand is greater than the second. */
    friend ApproximateKleenean operator>(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1._flt.dbl > x2._flt.dbl; }

    /** @brief Test equality of the operands. */
    friend ApproximateKleenean operator==(Rounded<FloatDP> x1, ExactDouble x2) { return x1._flt.dbl == x2.get_d(); }
    /** @brief Test inequality of the operands. */
    friend ApproximateKleenean operator!=(Rounded<FloatDP> x1, ExactDouble x2) { return x1._flt.dbl != x2.get_d(); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ApproximateKleenean operator<=(Rounded<FloatDP> x1, ExactDouble x2) { return x1._flt.dbl <= x2.get_d(); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ApproximateKleenean operator>=(Rounded<FloatDP> x1, ExactDouble x2) { return x1._flt.dbl >= x2.get_d(); }
    /** @brief Test whether the first operand is less than the second. */
    friend ApproximateKleenean operator<(Rounded<FloatDP> x1, ExactDouble x2) { return x1._flt.dbl < x2.get_d(); }
    /** @brief Test whether the first operand is greater than the second. */
    friend ApproximateKleenean operator>(Rounded<FloatDP> x1, ExactDouble x2) { return x1._flt.dbl > x2.get_d(); }

    /** @brief Test equality of the operands. */
    friend ApproximateKleenean operator==(ExactDouble x1, Rounded<FloatDP> x2) { return x1.get_d() == x2._flt.dbl; }
    /** @brief Test inequality of the operands. */
    friend ApproximateKleenean operator!=(ExactDouble x1, Rounded<FloatDP> x2) { return x1.get_d() != x2._flt.dbl; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ApproximateKleenean operator<=(ExactDouble x1, Rounded<FloatDP> x2) { return x1.get_d() <= x2._flt.dbl; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ApproximateKleenean operator>=(ExactDouble x1, Rounded<FloatDP> x2) { return x1.get_d() >= x2._flt.dbl; }
    /** @brief Test whether the first operand is less than the second. */
    friend ApproximateKleenean operator<(ExactDouble x1, Rounded<FloatDP> x2) { return x1.get_d() < x2._flt.dbl; }
    /** @brief Test whether the first operand is greater than the second. */
    friend ApproximateKleenean operator>(ExactDouble x1, Rounded<FloatDP> x2) { return x1.get_d() > x2._flt.dbl; }

    /** @brief Test whether the operands have the same representation. */
    friend Bool same(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1._flt == x2._flt; }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Rounded<FloatDP>& x) { return os << x.data(); }
};

/** @brief Numeric wrapper using the selected rounding policy. */
template<class FLT>
class Rounded {
    static const CurrentRoundingMode _rnd;
    FLT _flt;

  public:
    typedef RawTag Paradigm;
    typedef FLT FloatType;
    typedef Rounded<FloatType> NumericType;
    typedef typename FloatType::PrecisionType PrecisionType;
    typedef typename FloatType::PrecisionType CharacteristicsType;
    typedef typename FloatType::RoundingModeType RoundingModeType;

  public:
    /** @brief Return the rounding mode. */
    static RoundingModeType get_rounding_mode() { return FloatType::get_rounding_mode(); }
    /** @brief Set the rounding mode. */
    static Void set_rounding_mode(RoundingModeType rnd) { FloatType::set_rounding_mode(rnd); }
    /** @brief Set the rounding downward. */
    static Void set_rounding_downward() { FloatType::set_rounding_downward(); }
    /** @brief Set the rounding upward. */
    static Void set_rounding_upward() { FloatType::set_rounding_upward(); }
    /** @brief Set the rounding to nearest. */
    static Void set_rounding_to_nearest() { FloatType::set_rounding_to_nearest(); }
    /** @brief Set the rounding toward zero. */
    static Void set_rounding_toward_zero() { FloatType::set_rounding_toward_zero(); }

  public:
    /** @brief Return the underlying data pointer. */
    FloatType data() const { return _flt; }
    /** @brief Return the working precision. */
    PrecisionType precision() const { return this->_flt.precision(); }
    /** @brief Return the characteristics required to construct compatible values. */
    CharacteristicsType characteristics() const { return this->_flt.precision(); }

    /** @brief Construct a \c Rounded from pr. */
    explicit Rounded(PrecisionType pr) : _flt(pr) {}
    /** @brief Construct a \c Rounded from x. */
    Rounded(FloatType x) : _flt(x) {}
    /** @brief Construct a \c Rounded from x. */
    Rounded(Rounded<FloatType> x, PrecisionType) : _flt(x._flt) {}

    /** @brief Construct a \c Rounded from y, pr. */
    template<class Y>
        requires Constructible<FloatType, Y, RoundingModeType, PrecisionType>
    Rounded(const Y& y, PrecisionType pr) : Rounded(FloatType(y, FloatType::get_rounding_mode(), pr)) {}

    /** @brief Assign the supplied value and return this object. */
    template<class Y>
        requires Constructible<FloatType, Y, RoundingModeType, PrecisionType>
    Rounded<FloatType> operator=(const Y& y) {
        this->_flt = FloatType(y, FloatType::get_rounding_mode(), this->precision());
        return *this;
    }

    /** @brief Convert to float type. */
    explicit operator FloatType() const { return FloatType(this->_flt); }
    /** @brief Return the underlying numeric representation. */
    const FloatType& raw() const { return this->_flt; }

    // Approximation<FloatType> contains a FloatType member; it is not a FloatType base subobject.
    // Reinterpreting the wrapper as FloatType is not guaranteed by the C++ object model.
    /** @brief Construct a \c Rounded from x. */
    Rounded(const Approximation<FloatType>& x) : Rounded(x.raw()) {}
    /** @brief Convert to approximation< float type >. */
    operator Approximation<FloatType>() const;

    // Non-finiteness tests
    /** @brief Test whether the value is not a number. */
    friend Bool is_nan(Rounded<FloatType> x) { return is_nan(x._flt); }

    // Integer conversions
    /** @brief Return the greatest integer not exceeding the value. */
    friend Rounded<FloatType> floor(Rounded<FloatType> x) { return Rounded<FloatType>(floor(x._flt)); }
    /** @brief Return the least integer not less than the value. */
    friend Rounded<FloatType> ceil(Rounded<FloatType> x) { return Rounded<FloatType>(ceil(x._flt)); }
    /** @brief Round the value to an integer. */
    friend Rounded<FloatType> round(Rounded<FloatType> x) { return Rounded<FloatType>(round(x._flt)); }

    // Correctly rounded arithmetic
    /** @brief Return the additive zero with the operand characteristics. */
    friend Rounded<FloatType> nul(Rounded<FloatType> x) { return Rounded<FloatType>(nul(x._flt)); }
    /** @brief Return the unary positive value. */
    friend Rounded<FloatType> pos(Rounded<FloatType> x) { return Rounded<FloatType>(pos(x._flt)); }
    /** @brief Return the additive inverse. */
    friend Rounded<FloatType> neg(Rounded<FloatType> x) { return Rounded<FloatType>(neg(x._flt)); }
    /** @brief Divide the value by two. */
    friend Rounded<FloatType> hlf(Rounded<FloatType> x) { return Rounded<FloatType>(hlf(x._flt)); }
    /** @brief Return the square. */
    friend Rounded<FloatType> sqr(Rounded<FloatType> x) { return Rounded<FloatType>(sqr(_rnd, x._flt)); }
    /** @brief Return the multiplicative inverse. */
    friend Rounded<FloatType> rec(Rounded<FloatType> x) { return Rounded<FloatType>(rec(_rnd, x._flt)); }
    /** @brief Return the sum of the operands. */
    friend Rounded<FloatType> add(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(add(_rnd, x1._flt, x2._flt)); }
    /** @brief Return the difference of the operands. */
    friend Rounded<FloatType> sub(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(sub(_rnd, x1._flt, x2._flt)); }
    /** @brief Return the product of the operands. */
    friend Rounded<FloatType> mul(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(mul(_rnd, x1._flt, x2._flt)); }
    /** @brief Return the quotient of the operands. */
    friend Rounded<FloatType> div(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(div(_rnd, x1._flt, x2._flt)); }
    /** @brief Compute the product of the first two operands plus the third. */
    friend Rounded<FloatType> fma(Rounded<FloatType> x1, Rounded<FloatType> x2, Rounded<FloatType> x3) { return Rounded<FloatType>(fma(_rnd, x1._flt, x2._flt, x3._flt)); }
    /** @brief Raise the value to the specified power. */
    friend Rounded<FloatType> pow(Rounded<FloatType> x, Int n) { return Rounded<FloatType>(pow(_rnd, x._flt, n)); }
    /** @brief Return the square root. */
    friend Rounded<FloatType> sqrt(Rounded<FloatType> x) { return Rounded<FloatType>(sqrt(_rnd, x._flt)); }
    /** @brief Return the exponential. */
    friend Rounded<FloatType> exp(Rounded<FloatType> x) { return Rounded<FloatType>(exp(_rnd, x._flt)); }
    /** @brief Return the natural logarithm. */
    friend Rounded<FloatType> log(Rounded<FloatType> x) { return Rounded<FloatType>(log(_rnd, x._flt)); }
    /** @brief Return the sine. */
    friend Rounded<FloatType> sin(Rounded<FloatType> x) { return Rounded<FloatType>(sin(_rnd, x._flt)); }
    /** @brief Return the cosine. */
    friend Rounded<FloatType> cos(Rounded<FloatType> x) { return Rounded<FloatType>(cos(_rnd, x._flt)); }
    /** @brief Return the tangent. */
    friend Rounded<FloatType> tan(Rounded<FloatType> x) { return Rounded<FloatType>(tan(_rnd, x._flt)); }
    /** @brief Return the hyperbolic tangent. */
    friend Rounded<FloatType> tanh(Rounded<FloatType> x) { return Rounded<FloatType>(tanh(_rnd, x._flt)); }
    /** @brief Return the inverse sine. */
    friend Rounded<FloatType> asin(Rounded<FloatType> x) { return Rounded<FloatType>(asin(_rnd, x._flt)); }
    /** @brief Return the inverse cosine. */
    friend Rounded<FloatType> acos(Rounded<FloatType> x) { return Rounded<FloatType>(acos(_rnd, x._flt)); }
    /** @brief Return the inverse tangent. */
    friend Rounded<FloatType> atan(Rounded<FloatType> x) { return Rounded<FloatType>(atan(_rnd, x._flt)); }
    /** @brief Return pi at the requested precision. */
    static Rounded<FloatType> pi(PrecisionType pr) { return Rounded<FloatType>(FloatType::pi(FloatType::get_rounding_mode(), pr)); }

    /** @brief Return the absolute value. */
    friend Rounded<FloatType> abs(Rounded<FloatType> x) { return Rounded<FloatType>(abs(x._flt)); }
    /** @brief Return the maximum of the operands. */
    friend Rounded<FloatType> max(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(max(x1._flt, x2._flt)); }
    /** @brief Return the minimum of the operands. */
    friend Rounded<FloatType> min(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(min(x1._flt, x2._flt)); }
    /** @brief Return an upper bound on the absolute value. */
    friend Rounded<FloatType> mag(Rounded<FloatType> x) { return Rounded<FloatType>(abs(x._flt)); }
    /** @brief Return a lower bound on the absolute value. */
    friend Rounded<FloatType> mig(Rounded<FloatType> x) { return Rounded<FloatType>(abs(x._flt)); }

    /** @brief Return the midpoint value. */
    friend Rounded<FloatType> med(Rounded<FloatType> x1, Rounded<FloatType> x2) {
        return Rounded<FloatType>(hlf(add(CurrentRoundingMode(), x1._flt, x2._flt)));
    }
    /** @brief Return the radius. */
    friend Rounded<FloatType> rad(Rounded<FloatType> x1, Rounded<FloatType> x2) {
        return Rounded<FloatType>(hlf(sub(CurrentRoundingMode(), x2._flt, x1._flt)));
    }

    /** @brief Return the sum or unary positive value. */
    friend Rounded<FloatType> operator+(Rounded<FloatType> x) { return Rounded<FloatType>(pos(x._flt)); }
    /** @brief Return the difference or additive inverse. */
    friend Rounded<FloatType> operator-(Rounded<FloatType> x) { return Rounded<FloatType>(neg(x._flt)); }
    /** @brief Return the sum or unary positive value. */
    friend Rounded<FloatType> operator+(Rounded<FloatType> x1, Rounded<FloatType> x2) {
        return Rounded<FloatType>(add(_rnd, x1._flt, x2._flt));
    }
    /** @brief Return the difference or additive inverse. */
    friend Rounded<FloatType> operator-(Rounded<FloatType> x1, Rounded<FloatType> x2) {
        return Rounded<FloatType>(sub(_rnd, x1._flt, x2._flt));
    }
    /** @brief Return the product. */
    friend Rounded<FloatType> operator*(Rounded<FloatType> x1, Rounded<FloatType> x2) {
        return Rounded<FloatType>(mul(_rnd, x1._flt, x2._flt));
    }
    /** @brief Return the quotient. */
    friend Rounded<FloatType> operator/(Rounded<FloatType> x1, Rounded<FloatType> x2) {
        return Rounded<FloatType>(div(_rnd, x1._flt, x2._flt));
    }
    /** @brief Add the operand in place. */
    friend Rounded<FloatType>& operator+=(Rounded<FloatType>& x1, Rounded<FloatType> x2) { return x1 = x1 + x2; }
    /** @brief Subtract the operand in place. */
    friend Rounded<FloatType>& operator-=(Rounded<FloatType>& x1, Rounded<FloatType> x2) { return x1 = x1 - x2; }
    /** @brief Multiply by the operand in place. */
    friend Rounded<FloatType>& operator*=(Rounded<FloatType>& x1, Rounded<FloatType> x2) { return x1 = x1 * x2; }
    /** @brief Divide by the operand in place. */
    friend Rounded<FloatType>& operator/=(Rounded<FloatType>& x1, Rounded<FloatType> x2) { return x1 = x1 / x2; }

    // Convert mixed ExactDouble operands through the selected backend and
    // its active rounding mode; never hard-wire the independent exact arithmetic implementation.
    /** @brief Return the sum or unary positive value. */
    friend Rounded<FloatType> operator+(Rounded<FloatType> x1, ExactDouble x2) { return x1 + Rounded<FloatType>(x2, x1.precision()); }
    /** @brief Return the difference or additive inverse. */
    friend Rounded<FloatType> operator-(Rounded<FloatType> x1, ExactDouble x2) { return x1 - Rounded<FloatType>(x2, x1.precision()); }
    /** @brief Return the product. */
    friend Rounded<FloatType> operator*(Rounded<FloatType> x1, ExactDouble x2) { return x1 * Rounded<FloatType>(x2, x1.precision()); }
    /** @brief Return the quotient. */
    friend Rounded<FloatType> operator/(Rounded<FloatType> x1, ExactDouble x2) { return x1 / Rounded<FloatType>(x2, x1.precision()); }
    /** @brief Return the sum or unary positive value. */
    friend Rounded<FloatType> operator+(ExactDouble x1, Rounded<FloatType> x2) { return Rounded<FloatType>(x1, x2.precision()) + x2; }
    /** @brief Return the difference or additive inverse. */
    friend Rounded<FloatType> operator-(ExactDouble x1, Rounded<FloatType> x2) { return Rounded<FloatType>(x1, x2.precision()) - x2; }
    /** @brief Return the product. */
    friend Rounded<FloatType> operator*(ExactDouble x1, Rounded<FloatType> x2) { return Rounded<FloatType>(x1, x2.precision()) * x2; }
    /** @brief Return the quotient. */
    friend Rounded<FloatType> operator/(ExactDouble x1, Rounded<FloatType> x2) { return Rounded<FloatType>(x1, x2.precision()) / x2; }
    /** @brief Add the operand in place. */
    friend Rounded<FloatType>& operator+=(Rounded<FloatType>& x1, ExactDouble x2) { return x1 = x1 + x2; }
    /** @brief Subtract the operand in place. */
    friend Rounded<FloatType>& operator-=(Rounded<FloatType>& x1, ExactDouble x2) { return x1 = x1 - x2; }
    /** @brief Multiply by the operand in place. */
    friend Rounded<FloatType>& operator*=(Rounded<FloatType>& x1, ExactDouble x2) { return x1 = x1 * x2; }
    /** @brief Divide by the operand in place. */
    friend Rounded<FloatType>& operator/=(Rounded<FloatType>& x1, ExactDouble x2) { return x1 = x1 / x2; }

    /** @brief Return the maximum of the operands. */
    friend Rounded<FloatType> max(Rounded<FloatType> x1, ExactDouble x2) { return Rounded<FloatType>(max(x1._flt, x2)); }
    /** @brief Return the minimum of the operands. */
    friend Rounded<FloatType> min(Rounded<FloatType> x1, ExactDouble x2) { return Rounded<FloatType>(min(x1._flt, x2)); }
    /** @brief Return the maximum of the operands. */
    friend Rounded<FloatType> max(ExactDouble x1, Rounded<FloatType> x2) { return Rounded<FloatType>(max(x1, x2._flt)); }
    /** @brief Return the minimum of the operands. */
    friend Rounded<FloatType> min(ExactDouble x1, Rounded<FloatType> x2) { return Rounded<FloatType>(min(x1, x2._flt)); }

    /** @brief Test equality of the operands. */
    friend ApproximateKleenean operator==(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt == x2._flt; }
    /** @brief Test inequality of the operands. */
    friend ApproximateKleenean operator!=(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt != x2._flt; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ApproximateKleenean operator<=(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt <= x2._flt; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ApproximateKleenean operator>=(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt >= x2._flt; }
    /** @brief Test whether the first operand is less than the second. */
    friend ApproximateKleenean operator<(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt < x2._flt; }
    /** @brief Test whether the first operand is greater than the second. */
    friend ApproximateKleenean operator>(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt > x2._flt; }

    /** @brief Test equality of the operands. */
    friend ApproximateKleenean operator==(Rounded<FloatType> x1, ExactDouble x2) { return x1._flt == x2; }
    /** @brief Test inequality of the operands. */
    friend ApproximateKleenean operator!=(Rounded<FloatType> x1, ExactDouble x2) { return x1._flt != x2; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ApproximateKleenean operator<=(Rounded<FloatType> x1, ExactDouble x2) { return x1._flt <= x2; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ApproximateKleenean operator>=(Rounded<FloatType> x1, ExactDouble x2) { return x1._flt >= x2; }
    /** @brief Test whether the first operand is less than the second. */
    friend ApproximateKleenean operator<(Rounded<FloatType> x1, ExactDouble x2) { return x1._flt < x2; }
    /** @brief Test whether the first operand is greater than the second. */
    friend ApproximateKleenean operator>(Rounded<FloatType> x1, ExactDouble x2) { return x1._flt > x2; }

    /** @brief Test equality of the operands. */
    friend ApproximateKleenean operator==(ExactDouble x1, Rounded<FloatType> x2) { return x1 == x2._flt; }
    /** @brief Test inequality of the operands. */
    friend ApproximateKleenean operator!=(ExactDouble x1, Rounded<FloatType> x2) { return x1 != x2._flt; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ApproximateKleenean operator<=(ExactDouble x1, Rounded<FloatType> x2) { return x1 <= x2._flt; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ApproximateKleenean operator>=(ExactDouble x1, Rounded<FloatType> x2) { return x1 >= x2._flt; }
    /** @brief Test whether the first operand is less than the second. */
    friend ApproximateKleenean operator<(ExactDouble x1, Rounded<FloatType> x2) { return x1 < x2._flt; }
    /** @brief Test whether the first operand is greater than the second. */
    friend ApproximateKleenean operator>(ExactDouble x1, Rounded<FloatType> x2) { return x1 > x2._flt; }

    /** @brief Test whether the operands have the same representation. */
    friend Bool same(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt == x2._flt; }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Rounded<FloatType>& x) { return os << x.data(); }
};

template<class FLT>
const CurrentRoundingMode Rounded<FLT>::_rnd = CurrentRoundingMode();

} // namespace Ariadne

#endif
