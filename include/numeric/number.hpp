/***************************************************************************
 *            numeric/number.hpp
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

/** @file numeric/number.hpp
 *  @brief Generic numbers
 */

#ifndef ARIADNE_NUMBER_HPP
#define ARIADNE_NUMBER_HPP

#include "arithmetic.hpp"
#include "dyadic.hpp"
#include "float.decl.hpp"
#include "foundation/logical.decl.hpp"
#include "foundation/logical.hpp"
#include "foundation/paradigm.hpp"
#include "integer.hpp"
#include "number.decl.hpp"
#include "number_interface.hpp"
#include "rational.hpp"
#include "real.hpp"
#include "utility/handle.hpp"
#include "utility/prototype.hpp"

namespace Ariadne {

/************ Number *********************************************************/

template<class X>
struct IsNumber;

class NumberInterface;

template<class P>
class Number;
/** @brief Compile-time predicate testing whether a type is number. */
template<class P>
struct IsNumber<Number<P>> : True {};

/** @brief Dispatch exception for the supported operand types. */
struct DispatchException : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

/** @brief Compile-time trait for is precision. */
template<class PR>
struct IsPrecisionTrait : False {};
/** @brief Compile-time trait for is precision. */
template<>
struct IsPrecisionTrait<DoublePrecision> : True {};
/** @brief Compile-time trait for is precision. */
template<>
struct IsPrecisionTrait<MultiplePrecision> : True {};

template<class PR>
concept IsPrecision = IsPrecisionTrait<PR>::value;

template<class P>
Positive<Number<P>> cast_positive(Number<P> y);

/** @brief Declare number operators for the supported operand types. */
class DeclareNumberOperators {
    /** @brief Return the sum or unary positive value. */
    friend ApproximateNumber operator+(const ApproximateNumber&, const ApproximateNumber&);
    /** @brief Return the difference or additive inverse. */
    friend ApproximateNumber operator-(const ApproximateNumber&, const ApproximateNumber&);
    /** @brief Return the product. */
    friend ApproximateNumber operator*(const ApproximateNumber&, const ApproximateNumber&);
    /** @brief Return the quotient. */
    friend ApproximateNumber operator/(const ApproximateNumber&, const ApproximateNumber&);

    /** @brief Return the sum or unary positive value. */
    friend ValidatedLowerNumber operator+(const ValidatedLowerNumber& y1, const ValidatedLowerNumber& y2);
    /** @brief Return the difference or additive inverse. */
    friend ValidatedLowerNumber operator-(const ValidatedLowerNumber& y1, const ValidatedUpperNumber& y2);
    /** @brief Return the sum or unary positive value. */
    friend ValidatedUpperNumber operator+(const ValidatedUpperNumber& y1, const ValidatedUpperNumber& y2);
    /** @brief Return the difference or additive inverse. */
    friend ValidatedUpperNumber operator-(const ValidatedUpperNumber& y1, const ValidatedLowerNumber& y2);

    /** @brief Return the sum or unary positive value. */
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator+(const N& y1, const D& d2) { return y1 + Number<ApproximateTag>(d2); }
    /** @brief Return the difference or additive inverse. */
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator-(const N& y1, const D& d2) { return y1 - Number<ApproximateTag>(d2); }
    /** @brief Return the product. */
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator*(const N& y1, const D& d2) { return y1 * Number<ApproximateTag>(d2); }
    /** @brief Return the quotient. */
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator/(const N& y1, const D& d2) { return y1 / Number<ApproximateTag>(d2); }
    /** @brief Return the sum or unary positive value. */
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator+(const D& d1, const N& y2) { return Number<ApproximateTag>(d1) + y2; }
    /** @brief Return the difference or additive inverse. */
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator-(const D& d1, const N& y2) { return Number<ApproximateTag>(d1) - y2; }
    /** @brief Return the product. */
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator*(const D& d1, const N& y2) { return Number<ApproximateTag>(d1) * y2; }
    /** @brief Return the quotient. */
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator/(const D& d1, const N& y2) { return Number<ApproximateTag>(d1) / y2; }

    /** @brief Return the sum or unary positive value. */
    template<ConcreteNumber R, class P>
    friend decltype(auto)
        operator+(const R& r1, const Number<P>& y2) { return Number<Paradigm<R>>(r1) + y2; }
    /** @brief Return the sum or unary positive value. */
    template<ConcreteNumber R, class P>
    friend decltype(auto)
        operator+(const Number<P>& y1, const R& r2) { return y1 + Number<Paradigm<R>>(r2); }
};

template<class X, class P>
concept ConvertibleBuiltinFloatingPointToNumber = Same<P, ApproximateTag> and BuiltinFloatingPoint<X>;
template<class X, class P>
concept ConvertibleViaRealToNumber = WeakerThan<P, ParadigmTag<X>> and Convertible<X, Real>;
template<class X, class P>
concept ConvertibleViaNumberToNumber = WeakerThan<P, ParadigmTag<X>> and (not Convertible<X, Real>) and Convertible<X, Number<ParadigmTag<X>>>;

//! \ingroup NumericModule
//! \brief Generic numbers with computational paradigm \a P,  which may be %EffectiveTag, %ValidatedTag, %UpperTag, %LowerTag or %ApproximateTag.
template<class P>
class Number
    : public Handle<const NumberInterface>,
      public DeclareNumberOperators {
    static_assert(IsParadigm<P>, "P must be a paradigm");
    static_assert(Same<P, ExactTag> or Same<P, EffectiveTag> or Same<P, ValidatedTag> or Same<P, ApproximateTag>);

    template<class PP>
    friend class Number;
    template<class PP>
    friend class LowerNumber;
    template<class PP>
    friend class UpperNumber;

    template<class PR>
    using ResultFloatType = FloatType<Weaker<P, ValidatedTag>, PR>;

    template<class X>
    static const bool IsGettableAs = ANumber<X> and WeakerThan<typename X::Paradigm, P> and (not Same<typename X::Paradigm, ExactTag>);
    /** @brief Convert the concrete value through the generic number interface. */
    template<class X>
    static Number<ParadigmTag<X>> _convert_via_number(const X& x) {
        using StrongNumber = Number<ParadigmTag<X>>;
        return x.operator StrongNumber();
    }

  public:
    typedef NumberInterface Interface;
    typedef P Paradigm;
    typedef Number<P> NumericType;

  public:
    using Handle<const Interface>::Handle;
    /** @brief Construct a \c Number from p. */
    explicit Number(NumberInterface* p) : Handle<const NumberInterface>(p) {}

  private:

  public:
    /** @brief Return the handle to the represented object. */
    Handle<const NumberInterface> handle() const { return *this; }
    /** @brief Copy-construct the object from the supplied instance. */
    explicit Number(Handle<const NumberInterface> h) : Handle<const NumberInterface>(h) {}

  public:
    /** @brief Construct a \c Number with the specified characteristics. */
    Number() : Number(Integer(0)) {}

    // Construct from a Number of a stronger paradigm
    /** @brief Copy-construct the object from the supplied instance. */
    template<StrongerThan<P> SP>
    Number(const Number<SP>& y) : Number(y.handle()) {}

    //! Construct from a builtin integer
    template<BuiltinIntegral N>
    Number(const N& n) : Number(Integer(n)) {}
    // Construct from a builtin floating-point number
    /** @brief Construct a \c Number from x. */
    template<ConvertibleBuiltinFloatingPointToNumber<P> X>
    Number(const X& x) : Number(Dyadic(x)) {}

    // Construct from a type which is convertible to Real.
    /** @brief Construct a \c Number from x. */
    template<ConvertibleViaRealToNumber<P> X>
    Number(const X& x) : Number(_convert_via_number(x)) {}

    // Construct from a type which is convertible to another Number type.
    // TODO: Decide conversion characteristics from concrete type to Number<P>
    /** @brief Construct a \c Number from x. */
    template<ConvertibleViaNumberToNumber<P> X>
    explicit Number(const X& x) : Number(_convert_via_number(x)) {}

    //! \brief Get the value of the number as a double-precision floating-point type
    ResultFloatType<DoublePrecision> get(const DoublePrecision& prec) const { return this->ref()._get(P(), prec); }
    //! \brief Get the value of the number as a multiple-precision floating-point type
    //! Native multiple precision, preserving the same paradigm and rounding contract.
    ResultFloatType<MultiplePrecision> get(const MultiplePrecision& prec) const { return this->ref()._get(P(), prec); }

    //! \brief Get the value of the number as a floating-point ball with the given precision and error precision.
    template<class PR, class PRE>
        requires IsPrecision<PR> and IsPrecision<PRE> and WeakerThan<ValidatedTag, P>
    FloatBall<PR, PRE> get(PR const& prec, PRE const& errprec) const { return this->ref()._get(P(), prec, errprec); }

    //! \brief Try to dynamic_cast the object to concrete type \a X.
    template<class X>
    X extract() const;

    /** @brief Return the sum or unary positive value. */
    friend Number<P> operator+(const Number<P>& y) { return pos(y); }                            //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend Number<P> operator-(const Number<P>& y) { return neg(y); }                            //!< <p/>
    /** @brief Return the sum or unary positive value. */
    friend Number<P> operator+(const Number<P>& y1, const Number<P>& y2) { return add(y1, y2); } //!< <p/>
    /** @brief Return the difference or additive inverse. */
    friend Number<P> operator-(const Number<P>& y1, const Number<P>& y2) { return sub(y1, y2); } //!< <p/>
    /** @brief Return the product. */
    friend Number<P> operator*(const Number<P>& y1, const Number<P>& y2) { return mul(y1, y2); } //!< <p/>
    /** @brief Return the quotient. */
    friend Number<P> operator/(const Number<P>& y1, const Number<P>& y2) { return div(y1, y2); } //!< <p/>
    /** @brief Add the operand in place. */
    friend Number<P>& operator+=(Number<P>& y1, const Number<P>& y2) { return y1 = y1 + y2; }    //!< <p/>
    /** @brief Subtract the operand in place. */
    friend Number<P>& operator-=(Number<P>& y1, const Number<P>& y2) { return y1 = y1 - y2; }    //!< <p/>
    /** @brief Multiply by the operand in place. */
    friend Number<P>& operator*=(Number<P>& y1, const Number<P>& y2) { return y1 = y1 * y2; }    //!< <p/>
    /** @brief Divide by the operand in place. */
    friend Number<P>& operator/=(Number<P>& y1, const Number<P>& y2) { return y1 = y1 / y2; }    //!< <p/>

    /** @brief Return the additive zero with the operand characteristics. */
    friend Number<P> nul(const Number<P>& y) { return _apply<P>(Nul(), y); }                            //!< <p/>
    /** @brief Return the unary positive value. */
    friend Number<P> pos(const Number<P>& y) { return _apply<P>(Pos(), y); }                            //!< <p/>
    /** @brief Return the additive inverse. */
    friend Number<P> neg(const Number<P>& y) { return _apply<P>(Neg(), y); }                            //!< <p/>
    /** @brief Return the square. */
    friend Number<P> sqr(const Number<P>& y) { return _apply<P>(Sqr(), y); }                            //!< <p/>
    /** @brief Divide the value by two. */
    friend Number<P> hlf(const Number<P>& y) { return _apply<P>(Hlf(), y); }                            //!< <p/>
    /** @brief Return the multiplicative inverse. */
    friend Number<P> rec(const Number<P>& y) { return _apply<P>(Rec(), y); }                            //!< <p/>
    /** @brief Return the sum of the operands. */
    friend Number<P> add(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Add(), y1, y2); } //!< <p/>
    /** @brief Return the difference of the operands. */
    friend Number<P> sub(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Sub(), y1, y2); } //!< <p/>
    /** @brief Return the product of the operands. */
    friend Number<P> mul(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Mul(), y1, y2); } //!< <p/>
    /** @brief Return the quotient of the operands. */
    friend Number<P> div(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Div(), y1, y2); } //!< <p/>
    /** @brief Return the square root. */
    friend Number<P> sqrt(const Number<P>& y) { return _apply<P>(Sqrt(), y); }                          //!< <p/>
    /** @brief Return the exponential. */
    friend Number<P> exp(const Number<P>& y) { return _apply<P>(Exp(), y); }                            //!< <p/>
    /** @brief Return the natural logarithm. */
    friend Number<P> log(const Number<P>& y) { return _apply<P>(Log(), y); }                            //!< <p/>
    /** @brief Return the sine. */
    friend Number<P> sin(const Number<P>& y) { return _apply<P>(Sin(), y); }                            //!< <p/>
    /** @brief Return the cosine. */
    friend Number<P> cos(const Number<P>& y) { return _apply<P>(Cos(), y); }                            //!< <p/>
    /** @brief Return the tangent. */
    friend Number<P> tan(const Number<P>& y) { return _apply<P>(Tan(), y); }                            //!< <p/>
    /** @brief Return the hyperbolic tangent. */
    friend Number<P> tanh(const Number<P>& y) { return _apply<P>(Tanh(), y); }                          //!< <p/>
    /** @brief Return the inverse sine. */
    friend Number<P> asin(const Number<P>& y) { return _apply<P>(Asin(), y); }                          //!< <p/>
    /** @brief Return the inverse cosine. */
    friend Number<P> acos(const Number<P>& y) { return _apply<P>(Acos(), y); }                          //!< <p/>
    /** @brief Return the inverse tangent. */
    friend Number<P> atan(const Number<P>& y) { return _apply<P>(Atan(), y); }                          //!< <p/>

    /** @brief Raise the value to the specified power. */
    friend Number<P> pow(const Number<P>& y, Nat m) { return _apply<P>(Pow(), y, m); } //!< <p/>
    /** @brief Raise the value to the specified power. */
    friend Number<P> pow(const Number<P>& y, Int n) { return _apply<P>(Pow(), y, n); } //!< <p/>

    /** @brief Return the absolute value. */
    friend Positive<Number<P>> abs(const Number<P>& y) { return Positive<Number<P>>(_apply<P>(Abs(), y)); } //!< <p/>
    /** @brief Return the maximum of the operands. */
    friend Number<P> max(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Max(), y1, y2); }     //!< <p/>
    /** @brief Return the minimum of the operands. */
    friend Number<P> min(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Min(), y1, y2); }     //!< <p/>

    /** @brief Test equality of the operands. */
    friend LogicalType<EqualityInformation<P>> operator==(const Number<P>& y1, const Number<P>& y2) {
        return Detail::logical_type_from_owned_pointer<EqualityInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Equal()), y2.ptr()));
    }
    /** @brief Test whether the first operand is less than the second. */
    friend LogicalType<LessThanInformation<P>> operator<(const Number<P>& y1, const Number<P>& y2) {
        return Detail::logical_type_from_owned_pointer<LessThanInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Less()), y2.ptr()));
    }
    /** @brief Test whether the first operand is greater than the second. */
    friend LogicalType<LessThanInformation<P>> operator>(const Number<P>& y1, const Number<P>& y2) { return (y2 < y1); }                         //!< <p/>
    /** @brief Test inequality of the operands. */
    friend LogicalNegationType<LogicalType<EqualityInformation<P>>> operator!=(const Number<P>& y1, const Number<P>& y2) { return !(y1 == y2); } //!< <p/>
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend LogicalType<LessThanInformation<P>> operator<=(const Number<P>& y1, const Number<P>& y2) { return !(y1 > y2); }                       //!< <p/>
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend LogicalType<LessThanInformation<P>> operator>=(const Number<P>& y1, const Number<P>& y2) { return !(y1 < y2); }                       //!< <p/>

    //! \brief The name of the dynamic type held by the %Number object.
    String class_name() const { return this->ref()._class_name(); }

    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const Number<P>& y) { return y.ref()._write(os); }

  private:
    /** @brief Apply the supplied operation to the operands. */
    template<class RP, class OP>
    static Number<RP> _apply(OP op, const Number<P>& y) { return Number<RP>(y.ref()._apply(op)); }
    /** @brief Apply the supplied operation to the operands. */
    template<class RP, class OP>
    static Number<RP> _apply(OP op, const Number<P>& y, Int n) { return Number<RP>(y.ref()._apply(op, n)); }
    /** @brief Apply the supplied operation to the operands. */
    template<class RP, class OP, class P1, class P2>
    static Number<RP> _apply(OP op, const Number<P1>& y1, const Number<P2>& y2) {
        return Number<RP>(y1.ref()._apply(op, &y2.ref()));
    }
};

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<class P>
class Positive<Number<P>> : public Number<P> {
    /** @brief Return the unsigned magnitude of the integer. */
    friend const Number<P>& unsign(const Positive<Number<P>>& y) { return y; }

  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : Number<P>() {}
    /** @brief Construct a \c Positive from y. */
    explicit Positive(const Number<P>& y) : Number<P>(y) {}
    /** @brief Construct a \c Positive from n. */
    template<BuiltinUnsignedIntegral N>
    Positive(N n) : Number<P>(n) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<class N>
        requires Convertible<N, ExactNumber>
    Positive(const Positive<N>& n) : Number<P>(ExactNumber(static_cast<const N&>(n))) {}
    /** @brief Construct a \c Positive from n. */
    template<class N>
        requires Constructible<ExactNumber, N> and (not BuiltinIntegral<N>)
    explicit Positive(const N& n) : Number<P>(ExactNumber(n)) {}
    /** @brief Convert to number< p >. */
    explicit operator Number<P>() const { return *this; }

    /** @brief Return the sum or unary positive value. */
    friend Positive<Number<P>> operator+(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    /** @brief Return the product. */
    friend Positive<Number<P>> operator*(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(mul(unsign(y1), unsign(y2)));
    }
    /** @brief Return the sum of the operands. */
    friend Positive<Number<P>> add(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    /** @brief Return the product of the operands. */
    friend Positive<Number<P>> mul(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(mul(unsign(y1), unsign(y2)));
    }
    /** @brief Return the maximum of the operands. */
    friend Positive<Number<P>> max(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(max(unsign(y1), unsign(y2)));
    }
    /** @brief Return the absolute value. */
    friend Positive<Number<P>> abs(const Positive<Number<P>>& y) { return y; }
};

/** @brief Convert to the generic representation. */
template<class X>
decltype(auto) cast_generic(const X& x) { return x.generic(); }

/** @brief Convert to the positive representation. */
template<class P>
Positive<Number<P>> cast_positive(Number<P> y) {
    return Positive<Number<P>>(y);
}
/** @brief Convert to the exact representation. */
template<class P>
Positive<ExactNumber> cast_exact(const Positive<Number<P>>& y) {
    return cast_positive(cast_exact(unsign(y)));
}
/** @brief Convert to the exact representation. */
template<class P>
Positive<ExactNumber> cast_exact(const Positive<UpperNumber<P>>& y) {
    return cast_positive(cast_exact(unsign(y)));
}
/** @brief Convert to the exact representation. */
template<class P>
Positive<ExactNumber> cast_exact(const Positive<LowerNumber<P>>& y) {
    return cast_positive(cast_exact(unsign(y)));
}

} // namespace Ariadne

#endif
