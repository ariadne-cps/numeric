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
template<class P>
struct IsNumber<Number<P>> : True {};

struct DispatchException : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

template<class PR>
struct IsPrecisionTrait : False {};
template<>
struct IsPrecisionTrait<DoublePrecision> : True {};
template<>
struct IsPrecisionTrait<MultiplePrecision> : True {};

template<class PR>
concept IsPrecision = IsPrecisionTrait<PR>::value;

template<class P>
Positive<Number<P>> cast_positive(Number<P> y);

class DeclareNumberOperators {
    friend ApproximateNumber operator+(const ApproximateNumber&, const ApproximateNumber&);
    friend ApproximateNumber operator-(const ApproximateNumber&, const ApproximateNumber&);
    friend ApproximateNumber operator*(const ApproximateNumber&, const ApproximateNumber&);
    friend ApproximateNumber operator/(const ApproximateNumber&, const ApproximateNumber&);

    friend ValidatedLowerNumber operator+(const ValidatedLowerNumber& y1, const ValidatedLowerNumber& y2);
    friend ValidatedLowerNumber operator-(const ValidatedLowerNumber& y1, const ValidatedUpperNumber& y2);
    friend ValidatedUpperNumber operator+(const ValidatedUpperNumber& y1, const ValidatedUpperNumber& y2);
    friend ValidatedUpperNumber operator-(const ValidatedUpperNumber& y1, const ValidatedLowerNumber& y2);

    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator+(const N& y1, const D& d2) { return y1 + Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator-(const N& y1, const D& d2) { return y1 - Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator*(const N& y1, const D& d2) { return y1 * Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator/(const N& y1, const D& d2) { return y1 / Number<ApproximateTag>(d2); }
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator+(const D& d1, const N& y2) { return Number<ApproximateTag>(d1) + y2; }
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator-(const D& d1, const N& y2) { return Number<ApproximateTag>(d1) - y2; }
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator*(const D& d1, const N& y2) { return Number<ApproximateTag>(d1) * y2; }
    template<GenericNumber N, BuiltinFloatingPoint D>
    friend decltype(auto)
        operator/(const D& d1, const N& y2) { return Number<ApproximateTag>(d1) / y2; }

    template<ConcreteNumber R, class P>
    friend decltype(auto)
        operator+(const R& r1, const Number<P>& y2) { return Number<Paradigm<R>>(r1) + y2; }
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
    explicit Number(NumberInterface* p) : Handle<const NumberInterface>(p) {}

  private:

  public:
    Handle<const NumberInterface> handle() const { return *this; }
    explicit Number(Handle<const NumberInterface> h) : Handle<const NumberInterface>(h) {}

  public:
    Number() : Number(Integer(0)) {}

    // Construct from a Number of a stronger paradigm
    template<StrongerThan<P> SP>
    Number(const Number<SP>& y) : Number(y.handle()) {}

    //! Construct from a builtin integer
    template<BuiltinIntegral N>
    Number(const N& n) : Number(Integer(n)) {}
    // Construct from a builtin floating-point number
    template<ConvertibleBuiltinFloatingPointToNumber<P> X>
    Number(const X& x) : Number(Dyadic(x)) {}

    // Construct from a type which is convertible to Real.
    template<ConvertibleViaRealToNumber<P> X>
    Number(const X& x) : Number(_convert_via_number(x)) {}

    // Construct from a type which is convertible to another Number type.
    // TODO: Decide conversion characteristics from concrete type to Number<P>
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

    friend Number<P> operator+(const Number<P>& y) { return pos(y); }                            //!< <p/>
    friend Number<P> operator-(const Number<P>& y) { return neg(y); }                            //!< <p/>
    friend Number<P> operator+(const Number<P>& y1, const Number<P>& y2) { return add(y1, y2); } //!< <p/>
    friend Number<P> operator-(const Number<P>& y1, const Number<P>& y2) { return sub(y1, y2); } //!< <p/>
    friend Number<P> operator*(const Number<P>& y1, const Number<P>& y2) { return mul(y1, y2); } //!< <p/>
    friend Number<P> operator/(const Number<P>& y1, const Number<P>& y2) { return div(y1, y2); } //!< <p/>
    friend Number<P>& operator+=(Number<P>& y1, const Number<P>& y2) { return y1 = y1 + y2; }    //!< <p/>
    friend Number<P>& operator-=(Number<P>& y1, const Number<P>& y2) { return y1 = y1 - y2; }    //!< <p/>
    friend Number<P>& operator*=(Number<P>& y1, const Number<P>& y2) { return y1 = y1 * y2; }    //!< <p/>
    friend Number<P>& operator/=(Number<P>& y1, const Number<P>& y2) { return y1 = y1 / y2; }    //!< <p/>

    friend Number<P> nul(const Number<P>& y) { return _apply<P>(Nul(), y); }                            //!< <p/>
    friend Number<P> pos(const Number<P>& y) { return _apply<P>(Pos(), y); }                            //!< <p/>
    friend Number<P> neg(const Number<P>& y) { return _apply<P>(Neg(), y); }                            //!< <p/>
    friend Number<P> sqr(const Number<P>& y) { return _apply<P>(Sqr(), y); }                            //!< <p/>
    friend Number<P> hlf(const Number<P>& y) { return _apply<P>(Hlf(), y); }                            //!< <p/>
    friend Number<P> rec(const Number<P>& y) { return _apply<P>(Rec(), y); }                            //!< <p/>
    friend Number<P> add(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Add(), y1, y2); } //!< <p/>
    friend Number<P> sub(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Sub(), y1, y2); } //!< <p/>
    friend Number<P> mul(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Mul(), y1, y2); } //!< <p/>
    friend Number<P> div(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Div(), y1, y2); } //!< <p/>
    friend Number<P> sqrt(const Number<P>& y) { return _apply<P>(Sqrt(), y); }                          //!< <p/>
    friend Number<P> exp(const Number<P>& y) { return _apply<P>(Exp(), y); }                            //!< <p/>
    friend Number<P> log(const Number<P>& y) { return _apply<P>(Log(), y); }                            //!< <p/>
    friend Number<P> sin(const Number<P>& y) { return _apply<P>(Sin(), y); }                            //!< <p/>
    friend Number<P> cos(const Number<P>& y) { return _apply<P>(Cos(), y); }                            //!< <p/>
    friend Number<P> tan(const Number<P>& y) { return _apply<P>(Tan(), y); }                            //!< <p/>
    friend Number<P> tanh(const Number<P>& y) { return _apply<P>(Tanh(), y); }                          //!< <p/>
    friend Number<P> asin(const Number<P>& y) { return _apply<P>(Asin(), y); }                          //!< <p/>
    friend Number<P> acos(const Number<P>& y) { return _apply<P>(Acos(), y); }                          //!< <p/>
    friend Number<P> atan(const Number<P>& y) { return _apply<P>(Atan(), y); }                          //!< <p/>

    friend Number<P> pow(const Number<P>& y, Nat m) { return _apply<P>(Pow(), y, m); } //!< <p/>
    friend Number<P> pow(const Number<P>& y, Int n) { return _apply<P>(Pow(), y, n); } //!< <p/>

    friend Positive<Number<P>> abs(const Number<P>& y) { return Positive<Number<P>>(_apply<P>(Abs(), y)); } //!< <p/>
    friend Number<P> max(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Max(), y1, y2); }     //!< <p/>
    friend Number<P> min(const Number<P>& y1, const Number<P>& y2) { return _apply<P>(Min(), y1, y2); }     //!< <p/>

    friend LogicalType<EqualityInformation<P>> operator==(const Number<P>& y1, const Number<P>& y2) {
        return Detail::logical_type_from_owned_pointer<EqualityInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Equal()), y2.ptr()));
    }
    friend LogicalType<LessThanInformation<P>> operator<(const Number<P>& y1, const Number<P>& y2) {
        return Detail::logical_type_from_owned_pointer<LessThanInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Less()), y2.ptr()));
    }
    friend LogicalType<LessThanInformation<P>> operator>(const Number<P>& y1, const Number<P>& y2) { return (y2 < y1); }                         //!< <p/>
    friend LogicalNegationType<LogicalType<EqualityInformation<P>>> operator!=(const Number<P>& y1, const Number<P>& y2) { return !(y1 == y2); } //!< <p/>
    friend LogicalType<LessThanInformation<P>> operator<=(const Number<P>& y1, const Number<P>& y2) { return !(y1 > y2); }                       //!< <p/>
    friend LogicalType<LessThanInformation<P>> operator>=(const Number<P>& y1, const Number<P>& y2) { return !(y1 < y2); }                       //!< <p/>

    //! \brief The name of the dynamic type held by the %Number object.
    String class_name() const { return this->ref()._class_name(); }

    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const Number<P>& y) { return y.ref()._write(os); }

  private:
    template<class RP, class OP>
    static Number<RP> _apply(OP op, const Number<P>& y) { return Number<RP>(y.ref()._apply(op)); }
    template<class RP, class OP>
    static Number<RP> _apply(OP op, const Number<P>& y, Int n) { return Number<RP>(y.ref()._apply(op, n)); }
    template<class RP, class OP, class P1, class P2>
    static Number<RP> _apply(OP op, const Number<P1>& y1, const Number<P2>& y2) {
        return Number<RP>(y1.ref()._apply(op, &y2.ref()));
    }
};

template<class P>
class Positive<Number<P>> : public Number<P> {
    friend const Number<P>& unsign(const Positive<Number<P>>& y) { return y; }

  public:
    Positive() : Number<P>() {}
    explicit Positive(const Number<P>& y) : Number<P>(y) {}
    template<BuiltinUnsignedIntegral N>
    Positive(N n) : Number<P>(n) {}
    template<class N>
        requires Convertible<N, ExactNumber>
    Positive(const Positive<N>& n) : Number<P>(ExactNumber(static_cast<const N&>(n))) {}
    template<class N>
        requires Constructible<ExactNumber, N> and (not BuiltinIntegral<N>)
    explicit Positive(const N& n) : Number<P>(ExactNumber(n)) {}
    explicit operator Number<P>() const { return *this; }

    friend Positive<Number<P>> operator+(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    friend Positive<Number<P>> operator*(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(mul(unsign(y1), unsign(y2)));
    }
    friend Positive<Number<P>> add(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    friend Positive<Number<P>> mul(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(mul(unsign(y1), unsign(y2)));
    }
    friend Positive<Number<P>> max(const Positive<Number<P>>& y1, const Positive<Number<P>>& y2) {
        return cast_positive(max(unsign(y1), unsign(y2)));
    }
    friend Positive<Number<P>> abs(const Positive<Number<P>>& y) { return y; }
};

template<class X>
decltype(auto) cast_generic(const X& x) { return x.generic(); }

template<class P>
Positive<Number<P>> cast_positive(Number<P> y) {
    return Positive<Number<P>>(y);
}
template<class P>
Positive<ExactNumber> cast_exact(const Positive<Number<P>>& y) {
    return cast_positive(cast_exact(unsign(y)));
}
template<class P>
Positive<ExactNumber> cast_exact(const Positive<UpperNumber<P>>& y) {
    return cast_positive(cast_exact(unsign(y)));
}
template<class P>
Positive<ExactNumber> cast_exact(const Positive<LowerNumber<P>>& y) {
    return cast_positive(cast_exact(unsign(y)));
}

} // namespace Ariadne

#endif
