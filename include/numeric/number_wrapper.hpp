/***************************************************************************
 *            numeric/number_wrapper.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/** @file numeric/number_wrapper.hpp
 *  @brief Declarations and utilities for number wrapper in Ariadne Numeric.
 */

#ifndef ARIADNE_NUMBER_WRAPPER_HPP
#define ARIADNE_NUMBER_WRAPPER_HPP

#include <cassert>

#include "builtin.hpp"
#include "float_approximation.hpp"
#include "float_ball.hpp"
#include "float_bounds.hpp"
#include "float_error.hpp"
#include "float_lower_bound.hpp"
#include "float_upper_bound.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "foundation/logical.hpp"
#include "foundation/paradigm.hpp"
#include "number.hpp"
#include "number_interface.hpp"
#include "numeric/operators.hpp"
#include "utility/dispatching.hpp"
#include "utility/module.hpp"

namespace Ariadne {

/************ Number *********************************************************/

// Declare fallbacks for Integer, Dyadic and Rational numbers
Real sqrt(const Real&);
Real exp(const Real&);
Real log(const Real&);
Real sin(const Real&);
Real cos(const Real&);
Real tan(const Real&);
Real tanh(const Real&);
Real asin(const Real&);
Real acos(const Real&);
Real atan(const Real&);

// Declare fallbacks for use by Lower/UpperBound
Approximation<FloatDP> div(const Approximation<FloatDP>& x1, const Approximation<FloatDP>& x2);
Approximation<FloatDP> sqr(const Approximation<FloatDP>& x);
Approximation<FloatDP> rec(const Approximation<FloatDP>& x);
Approximation<FloatDP> sin(const Approximation<FloatDP>& x);
Approximation<FloatDP> cos(const Approximation<FloatDP>& x);
Approximation<FloatDP> tan(const Approximation<FloatDP>& x);
Approximation<FloatDP> tanh(const Approximation<FloatDP>& x);
Approximation<FloatDP> asin(const Approximation<FloatDP>& x);
Approximation<FloatDP> acos(const Approximation<FloatDP>& x);
Approximation<FloatDP> pow(const Approximation<FloatDP>& x, const Integer& z);

Approximation<FloatMP> div(const Approximation<FloatMP>& x1, const Approximation<FloatMP>& x2);
Approximation<FloatMP> sqr(const Approximation<FloatMP>& x);
Approximation<FloatMP> rec(const Approximation<FloatMP>& x);
Approximation<FloatMP> sin(const Approximation<FloatMP>& x);
Approximation<FloatMP> cos(const Approximation<FloatMP>& x);
Approximation<FloatMP> tan(const Approximation<FloatMP>& x);
Approximation<FloatMP> tanh(const Approximation<FloatMP>& x);
Approximation<FloatMP> asin(const Approximation<FloatMP>& x);
Approximation<FloatMP> acos(const Approximation<FloatMP>& x);
Approximation<FloatMP> pow(const Approximation<FloatMP>& x, const Integer& z);

// Declare fallbacks for use by Dyadic/RationalBounds
template<class B>
concept AlgebraicBounds = SameAs<B, DyadicBounds> || SameAs<B, DecimalBounds> || SameAs<B, RationalBounds>;
template<AlgebraicBounds B>
B sqrt(const B& x) {
    ARIADNE_THROW(DispatchException, "sqrt(AlgebraicBounds)", "Unsupported operation sqrt(" << x << ")");
}
template<AlgebraicBounds B>
B exp(const B& x) {
    ARIADNE_THROW(DispatchException, "exp(AlgebraicBounds)", "Unsupported operation exp(" << x << ")");
}
template<AlgebraicBounds B>
B log(const B& x) {
    ARIADNE_THROW(DispatchException, "log(AlgebraicBounds)", "Unsupported operation log(" << x << ")");
}
template<AlgebraicBounds B>
B sin(const B& x) {
    ARIADNE_THROW(DispatchException, "sin(AlgebraicBounds)", "Unsupported operation sin(" << x << ")");
}
template<AlgebraicBounds B>
B cos(const B& x) {
    ARIADNE_THROW(DispatchException, "cos(AlgebraicBounds)", "Unsupported operation cos(" << x << ")");
}
template<AlgebraicBounds B>
B tan(const B& x) {
    ARIADNE_THROW(DispatchException, "tan(AlgebraicBounds)", "Unsupported operation tan(" << x << ")");
}
template<AlgebraicBounds B>
B tanh(const B& x) {
    ARIADNE_THROW(DispatchException, "tanh(AlgebraicBounds)", "Unsupported operation tanh(" << x << ")");
}
template<AlgebraicBounds B>
B asin(const B& x) {
    ARIADNE_THROW(DispatchException, "asin(AlgebraicBounds)", "Unsupported operation asin(" << x << ")");
}
template<AlgebraicBounds B>
B acos(const B& x) {
    ARIADNE_THROW(DispatchException, "acos(AlgebraicBounds)", "Unsupported operation acos(" << x << ")");
}
template<AlgebraicBounds B>
B atan(const B& x) {
    ARIADNE_THROW(DispatchException, "atan(AlgebraicBounds)", "Unsupported operation atan(" << x << ")");
}
template<AlgebraicBounds B>
ValidatedKleenean operator>(const B& y1, Int n2) { return y1 > B(n2); }

class NumberInterface;
template<class X>
class NumberMixin;
template<class X>
class NumberWrapper;

template<>
struct InterfaceTraits<NumberInterface> {
    template<class X>
    using MixinType = NumberMixin<X>;
    template<class X>
    using WrapperType = NumberWrapper<X>;
};

template<class X>
inline const X* extract(const NumberInterface* y) {
    return dynamic_cast<const NumberWrapper<X>*>(y);
}

inline OutputStream& operator<<(OutputStream& os, const NumberInterface& y) { return y._write(os); }

// FIXME: Should test for other potential infinities
inline Comparison cmp(const NumberInterface& y1, const NumberInterface& y2) {
    Comparison res;
    const FloatDP* x1 = extract<FloatDP>(&y1);
    const FloatDP* x2 = extract<FloatDP>(&y2);
    if (x1) {
        if (x2) {
            res = cmp(ExactDouble(x1->get_d()), ExactDouble(x2->get_d()));
        } else {
            res = cmp(ExactDouble(x1->get_d()), y2._get_q());
        }
    } else {
        if (x2) {
            res = cmp(y1._get_q(), ExactDouble(x2->get_d()));
        } else {
            res = cmp(y1._get_q(), y2._get_q());
        }
    }
    return res;
}

/*
template<class OP, class X1, class X2> NumberInterface* _concrete_apply(OP op, X1 const& x1, X2 const& x2) {
    if constexpr (CanApply<OP,X1,X2>::value) {
        return _make_number_wrapper(op(x1,x2));
    } else {
        String yc1=class_name<X1>(); String yc2=class_name<X2>();
        ARIADNE_THROW(DispatchException,op<<"(Number y1, Number y2) with y1="<<x1<<", y2="<<x2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
    }
}
*/

template<class R>
NumberInterface* _make_number_wrapper(const R& r) { return new NumberWrapper<R>(r); }
struct MakeNumberWrapper {
    template<class R>
    NumberInterface* operator()(const R& r) const { return new NumberWrapper<R>(r); }
};

template<class X>
inline X make_unsigned(X x) { return x; }
template<class X>
inline X make_unsigned(Positive<X> px) {
#ifndef __clang__
    return px;
#else
    return std::move(px);
#endif
}
inline Real make_unsigned(PositiveReal px) {
#ifndef __clang__
    return px;
#else
    return std::move(px);
#endif
}
inline Integer make_unsigned(Natural px) {
#ifndef __clang__
    return px;
#else
    return std::move(px);
#endif
}

// Instantiate a unary operation only when it is defined for the concrete
// numeric type.  Directed lower/upper bounds deliberately do not expose all
// elementary functions (e.g. squaring an arbitrary one-sided bound).
template<class R, class OP, class X>
inline R _concrete_apply_unary(OP op, const X& x) {
    static_assert(Same<R, NumberInterface*>);
    // These AlgebraicBounds overloads exist only to report unsupported
    // transcendental operations. They always throw: do not instantiate a
    // value-returning expression around them (MSVC ARM64 diagnoses it as
    // unreachable code under /WX).
    if constexpr (AlgebraicBounds<X> &&
                  (Same<OP, Sqrt> || Same<OP, Exp> || Same<OP, Log> ||
                   Same<OP, Sin> || Same<OP, Cos> || Same<OP, Tan> ||
                   Same<OP, Tanh> || Same<OP, Asin> || Same<OP, Acos> ||
                   Same<OP, Atan>)) {
        ARIADNE_THROW(DispatchException, "_concrete_apply(UnaryElementaryOperator)",
                      "Unsupported unary operation on " << class_name<X>());
    } else if constexpr (requires { op(x); }) {
        return _make_number_wrapper(make_unsigned(op(x)));
    } else {
        ARIADNE_THROW(DispatchException, "_concrete_apply(UnaryElementaryOperator)",
                      "Unsupported unary operation on " << class_name<X>());
    }
}

template<class R, class X>
inline R _concrete_apply(UnaryElementaryOperator op, const X& x) {
    static_assert(Same<R, NumberInterface*>);
    switch (op.code()) {
        case OperatorCode::NUL:
            return _concrete_apply_unary<R>(Nul(), x);
        case OperatorCode::POS:
            return _concrete_apply_unary<R>(Pos(), x);
        case OperatorCode::NEG:
            return _concrete_apply_unary<R>(Neg(), x);
        case OperatorCode::SQR:
            return _concrete_apply_unary<R>(Sqr(), x);
        case OperatorCode::HLF:
            return _concrete_apply_unary<R>(Hlf(), x);
        case OperatorCode::REC:
            return _concrete_apply_unary<R>(Rec(), x);
        case OperatorCode::SQRT:
            return _concrete_apply_unary<R>(Sqrt(), x);
        case OperatorCode::EXP:
            return _concrete_apply_unary<R>(Exp(), x);
        case OperatorCode::LOG:
            return _concrete_apply_unary<R>(Log(), x);
        case OperatorCode::SIN:
            return _concrete_apply_unary<R>(Sin(), x);
        case OperatorCode::COS:
            return _concrete_apply_unary<R>(Cos(), x);
        case OperatorCode::TAN:
            return _concrete_apply_unary<R>(Tan(), x);
        case OperatorCode::TANH:
            return _concrete_apply_unary<R>(Tanh(), x);
        case OperatorCode::ASIN:
            return _concrete_apply_unary<R>(Asin(), x);
        case OperatorCode::ACOS:
            return _concrete_apply_unary<R>(Acos(), x);
        case OperatorCode::ATAN:
            return _concrete_apply_unary<R>(Atan(), x);
        case OperatorCode::ABS:
            return _concrete_apply_unary<R>(Abs(), x);
        default:
            ARIADNE_THROW(DispatchException, "_concrete_apply(UnaryElementaryOperator,X const&)", "Unsupported unary operator code " << op.code());
    }
}

// Preserve native overload dispatch for every supported binary operator.
// A missing overload is an unsupported type/operator combination, not a
// reason to instantiate ill-formed code or to change the numerical backend.
template<class R, class X1, class X2>
inline R _concrete_apply(BinaryElementaryOperator op, X1 const& x1, X2 const& x2) {
    static_assert(Same<R, NumberInterface*>);
    return op.accept([&x1, &x2, &op](auto _op) -> R {
        if constexpr (requires { _op(x1, x2); }) {
            return _make_number_wrapper(_op(x1, x2));
        } else {
            ARIADNE_THROW(DispatchException, "_concrete_apply(BinaryElementaryOperator)",
                          "Unsupported " << op << " on " << class_name<X1>() << " and " << class_name<X2>());
        }
    });
}

template<class R, class X1, class X2>
inline R _concrete_apply_max_or_min(BinaryElementaryOperator op, X1 const& x1, X2 const& x2) {
    if (op.code() == BinaryElementaryOperator(Max()).code()) {
        return _make_number_wrapper(max(x1, x2));
    } else if (op.code() == BinaryElementaryOperator(Min()).code()) {
        return _make_number_wrapper(min(x1, x2));
    }
    String yc1 = class_name<X1>();
    String yc2 = class_name<X2>();
    ARIADNE_THROW(DispatchException, op << "(Number y1, Number y2) with y1=" << x1 << ", y2=" << x2, "No dispatch for " << op << "(" << yc1 << ", " << yc2 << ")");
}

// FIXME: Prefer symbolic dispatch
template<class R, ARawFloat F>
inline R _concrete_apply(BinaryElementaryOperator op, const F& x1, const F& x2) {
    return _concrete_apply_max_or_min<R>(op, x1, x2);
}
template<class R, ARawFloat F>
inline R _concrete_apply(BinaryElementaryOperator op, const F& x1, const Integer& z2) {
    return _concrete_apply_max_or_min<R>(op, x1, z2);
}
template<class R, ARawFloat F>
inline R _concrete_apply(BinaryElementaryOperator op, const Integer& z1, const F& x2) {
    return _concrete_apply_max_or_min<R>(op, z1, x2);
}
template<class R, ARawFloat F>
inline R _concrete_apply(BinaryElementaryOperator op, const F& x1, const Dyadic& w2) {
    return _concrete_apply_max_or_min<R>(op, x1, w2);
}
template<class R, ARawFloat F>
inline R _concrete_apply(BinaryElementaryOperator op, const Dyadic& w1, const F& x2) {
    return _concrete_apply_max_or_min<R>(op, w1, x2);
}
template<class R, ARawFloat F>
inline R _concrete_apply(BinaryElementaryOperator op, const F& x1, const Rational& q2) {
    String yc1 = class_name<F>();
    String yc2 = class_name<Rational>();
    ARIADNE_THROW(DispatchException, op << "(Number y1, Number y2) with y1=" << x1 << ", y2=" << q2, "No dispatch for " << op << "(" << yc1 << ", " << yc2 << ")");
}
template<class R, ARawFloat F>
inline R _concrete_apply(BinaryElementaryOperator op, const Rational& q1, const F& x2) {
    String yc1 = class_name<Rational>();
    String yc2 = class_name<F>();
    ARIADNE_THROW(DispatchException, op << "(Number y1, Number y2) with y1=" << q1 << ", y2=" << x2, "No dispatch for " << op << "(" << yc1 << ", " << yc2 << ")");
}

template<class R, class OP, class X1, class X2>
inline R _concrete_operator_apply(OP op, X1 const& x1, X2 const& x2) {
    auto res = op(x1, x2);
    if constexpr (Same<decltype(res), Bool>) {
        return new_logical_pointer_from_value(LogicalValue(res));
    } else if constexpr (Same<decltype(res.repr()), const LogicalValue&>) {
        return new_logical_pointer_from_value(res.repr());
    } else {
        static_assert(Same<decltype(res.repr()), const LogicalHandle&>);
        return res.repr().pointer()->_copy();
    }
}

template<class R, class X1, class X2>
inline R _concrete_apply(BinaryComparisonOperator op, X1 const& x1, X2 const& x2) {
    static_assert(Same<R, LogicalInterface*>);
    return op.accept([&x1, &x2](auto _op) { return _concrete_operator_apply<R>(_op, x1, x2); });
}

class AlgebraicNumberInterface;
class ValidatedAlgebraicNumberInterface;
template<class F>
class ConcreteNumberInterface;
template<class F, class FE>
class ConcreteBallInterface;

template<>
struct Managed<AlgebraicNumberInterface> {
    typedef Aware<Integer, Dyadic, Rational, Real> Types;
    // FIXME: typedef Aware<ExactDouble,Integer,Dyadic,Rational,Real> Types;
};
template<>
struct Managed<ValidatedAlgebraicNumberInterface> {
    typedef Aware<DyadicBounds, RationalBounds> Types;
    // FIXME: typedef Aware<ExactDouble,Integer,Dyadic,Rational,Real> Types;
};
template<class F>
struct Managed<ConcreteNumberInterface<F>> {
    typedef Aware<F, Ball<F>, Bounds<F>, UpperBound<F>, LowerBound<F>, Approximation<F>> Types;
};
template<class F, class FE>
struct Managed<ConcreteBallInterface<F, FE>> {
    typedef Aware<Ball<F, FE>> Types;
};

template<>
struct DispatcherTraits<ExactDouble> {
    typedef AlgebraicNumberInterface Interface;
};
template<>
struct DispatcherTraits<Integer> {
    typedef AlgebraicNumberInterface Interface;
};
template<>
struct DispatcherTraits<Dyadic> {
    typedef AlgebraicNumberInterface Interface;
};
template<>
struct DispatcherTraits<Rational> {
    typedef AlgebraicNumberInterface Interface;
};
template<>
struct DispatcherTraits<Real> {
    typedef AlgebraicNumberInterface Interface;
};

template<>
struct DispatcherTraits<DyadicBounds> {
    typedef ValidatedAlgebraicNumberInterface Interface;
};
template<>
struct DispatcherTraits<RationalBounds> {
    typedef ValidatedAlgebraicNumberInterface Interface;
};

template<>
struct DispatcherTraits<FloatDP> {
    typedef ConcreteNumberInterface<FloatDP> Interface;
};
template<>
struct DispatcherTraits<FloatMP> {
    typedef ConcreteNumberInterface<FloatMP> Interface;
};

template<class F>
struct DispatcherTraits<Ball<F>> {
    typedef ConcreteNumberInterface<F> Interface;
};
template<class F>
struct DispatcherTraits<Bounds<F>> {
    typedef ConcreteNumberInterface<F> Interface;
};
template<class F>
struct DispatcherTraits<UpperBound<F>> {
    typedef ConcreteNumberInterface<F> Interface;
};
template<class F>
struct DispatcherTraits<LowerBound<F>> {
    typedef ConcreteNumberInterface<F> Interface;
};
template<class F>
struct DispatcherTraits<Approximation<F>> {
    typedef ConcreteNumberInterface<F> Interface;
};

template<class F, class FE>
struct DispatcherTraits<Ball<F, FE>> {
    typedef ConcreteBallInterface<F, FE> Interface;
};

class AlgebraicNumberInterface
    : public virtual SelfOperableInterface<NumberInterface*, BinaryElementaryOperator, ManagedTypes<AlgebraicNumberInterface>> {
};

class ValidatedAlgebraicNumberInterface
    : public virtual SelfOperableInterface<NumberInterface*, BinaryElementaryOperator, ManagedTypes<ValidatedAlgebraicNumberInterface>> {
};

template<class F>
class ConcreteNumberInterface
    : public virtual SelfOperableInterface<NumberInterface*, BinaryElementaryOperator, ManagedTypes<ConcreteNumberInterface<F>>>,
      public virtual OperableInterface<NumberInterface*, BinaryElementaryOperator, ManagedTypes<AlgebraicNumberInterface>> {
};

template<class F, class FE>
class ConcreteBallInterface
    : public virtual SelfOperableInterface<NumberInterface*, BinaryElementaryOperator, ManagedTypes<ConcreteBallInterface<F, FE>>>,
      public virtual OperableInterface<NumberInterface*, BinaryElementaryOperator, ManagedTypes<AlgebraicNumberInterface>>
//    , public virtual OperableInterface<NumberInterface*,BinaryElementaryOperator,ManagedTypes<ConcreteNumberInterface<F>>>
{
};

template<class X, class DI = DispatcherInterface<X>>
class ElementaryBinaryNumberDispatcherMixin;

template<class Y>
class ElementaryBinaryNumberDispatcherMixin<Y, AlgebraicNumberInterface>
    : public SelfOperableMixin<Y, NumberInterface, NumberInterface*, BinaryElementaryOperator, ManagedTypes<AlgebraicNumberInterface>> {
};

template<class Y>
class ElementaryBinaryNumberDispatcherMixin<Y, ValidatedAlgebraicNumberInterface>
    : public SelfOperableMixin<Y, NumberInterface, NumberInterface*, BinaryElementaryOperator, ManagedTypes<ValidatedAlgebraicNumberInterface>> {
};

template<class Y, class F>
class ElementaryBinaryNumberDispatcherMixin<Y, ConcreteNumberInterface<F>>
    : public SelfOperableMixin<Y, NumberInterface, NumberInterface*, BinaryElementaryOperator, ManagedTypes<ConcreteNumberInterface<F>>>,
      public OperableMixin<Y, NumberInterface, NumberInterface*, BinaryElementaryOperator, ManagedTypes<AlgebraicNumberInterface>> {
};

template<class Y, class F, class FE>
class ElementaryBinaryNumberDispatcherMixin<Y, ConcreteBallInterface<F, FE>>
    : public SelfOperableMixin<Y, NumberInterface, NumberInterface*, BinaryElementaryOperator, ManagedTypes<ConcreteBallInterface<F, FE>>>,
      public OperableMixin<Y, NumberInterface, NumberInterface*, BinaryElementaryOperator, ManagedTypes<AlgebraicNumberInterface>>
//    , public OperableMixin<Y,NumberInterface,NumberInterface*,BinaryElementaryOperator,ManagedTypes<ConcreteNumberInterface<F>>>
{
};

template<class X, class DI = DispatcherInterface<X>>
class ComparisonBinaryNumberDispatcherMixin;

template<class Y>
class ComparisonBinaryNumberDispatcherMixin<Y, AlgebraicNumberInterface>
    : public SelfOperableMixin<Y, NumberInterface, LogicalInterface*, BinaryComparisonOperator, ManagedTypes<AlgebraicNumberInterface>> {
};

template<class Y>
class ComparisonBinaryNumberDispatcherMixin<Y, ValidatedAlgebraicNumberInterface>
    : public SelfOperableMixin<Y, NumberInterface, LogicalInterface*, BinaryComparisonOperator, ManagedTypes<ValidatedAlgebraicNumberInterface>> {
};

template<class Y, class F>
class ComparisonBinaryNumberDispatcherMixin<Y, ConcreteNumberInterface<F>>
    : public SelfOperableMixin<Y, NumberInterface, LogicalInterface*, BinaryComparisonOperator, ManagedTypes<ConcreteNumberInterface<F>>>,
      public OperableMixin<Y, NumberInterface, LogicalInterface*, BinaryComparisonOperator, ManagedTypes<AlgebraicNumberInterface>> {
};

template<class Y, class F, class FE>
class ComparisonBinaryNumberDispatcherMixin<Y, ConcreteBallInterface<F, FE>>
    : public SelfOperableMixin<Y, NumberInterface, LogicalInterface*, BinaryComparisonOperator, ManagedTypes<ConcreteBallInterface<F, FE>>>,
      public OperableMixin<Y, NumberInterface, LogicalInterface*, BinaryComparisonOperator, ManagedTypes<AlgebraicNumberInterface>> {
};

inline LogicalInterface* make_symbolic(BinaryComparisonOperator op, const NumberInterface* yp1, const NumberInterface* yp2) {
    String yc1 = yp1->_class_name();
    String yc2 = yp2->_class_name();
    ARIADNE_THROW(DispatchException, op << "(Number y1, Number y2) with y1=" << *yp1 << ", y2=" << *yp2, "No dispatch for " << op << "(" << yc1 << ", " << yc2 << ")");
}

template<class OP>
inline NumberInterface* make_symbolic(OP op, const NumberInterface* yp1, const NumberInterface* yp2) {
    //    Handle<NumberInterface> y1(const_cast<NumberInterface*>(yp1)->shared_from_this());
    //    Handle<NumberInterface> y2(const_cast<NumberInterface*>(yp2)->shared_from_this());
    String yc1 = yp1->_class_name();
    String yc2 = yp2->_class_name();
    ARIADNE_THROW(DispatchException, op << "(Number y1, Number y2) with y1=" << *yp1 << ", y2=" << *yp2, "No dispatch for " << op << "(" << yc1 << ", " << yc2 << ")");
}

template<class X>
struct ElementaryUnaryNumberOperationsMixin
    : public UnaryOperationMixin<X, NumberInterface, UnaryElementaryOperator, NumberInterface> {};
template<class X>
struct ElementaryBinaryNumberOperationsMixin
    : public BinaryOperationMixin<X, NumberInterface, BinaryElementaryOperator, NumberInterface> {};
template<class X, class N = Int>
struct ElementaryGradedNumberOperationsMixin
    : public GradedOperationMixin<X, NumberInterface, GradedElementaryOperator, NumberInterface, N> {};

template<class X>
struct ComparisonBinaryNumberOperationsMixin
    : public BinaryOperationMixin<X, LogicalInterface, BinaryComparisonOperator, NumberInterface> {};

inline OutputStream& operator<<(OutputStream& os, ParadigmCode cd) {
    switch (cd) {
        case ParadigmCode::EXACT:
            return os << "EXACT";
        case ParadigmCode::VALIDATED:
            return os << "VALIDATED";
        case ParadigmCode::APPROXIMATE:
            return os << "APPROXIMATE";
        default:
            return os << "UNKNOWN";
    }
}

template<class X>
class NumberGetterMixin : public virtual NumberInterface {
  public:
    //  operator X const& () const { return static_cast<NumberMixin<X>const&>(*this); }
    static const X& _cast(const NumberGetterMixin<X>& self) { return static_cast<const NumberWrapper<X>&>(static_cast<const NumberMixin<X>&>(self)); }
    static X& _cast(NumberGetterMixin<X>& self) { return static_cast<NumberWrapper<X>&>(static_cast<NumberMixin<X>&>(self)); }

    typedef Paradigm<X> P;
    friend class Number<P>;

    virtual NumberInterface* _copy() const override { return new NumberWrapper<X>(_cast(*this)); }
    virtual NumberInterface* _move() override { return new NumberWrapper<X>(std::move(_cast(*this))); }

    virtual LogicalInterface* _is_pos() const override {
        auto res = (_cast(*this) > 0);
        if constexpr (Same<decltype(res), bool>) {
            return new_logical_pointer_from_value(LogicalValue(res));
        } else if constexpr (Same<decltype(res.repr()), const LogicalValue&>) {
            return new_logical_pointer_from_value(res.repr());
        } else {
            static_assert(Same<decltype(res.repr()), const LogicalHandle&>);
            return res.repr().pointer()->_copy();
        }
    }

    virtual Rational _get_q() const override {
        return this->_get_as<Rational>();
    }

    virtual FloatDPBall _get(MetricTag, DoublePrecision pr, DoublePrecision) const override {
        return this->_get_as<FloatDPBall>(pr);
    }
    virtual FloatDPBounds _get(OrderTag, DoublePrecision pr) const override {
        return this->_get_as<FloatDPBounds>(pr);
    }
    virtual FloatDPUpperBound _get(UpperTag, DoublePrecision pr) const override {
        return this->_get_as<FloatDPUpperBound>(pr);
    }
    virtual FloatDPLowerBound _get(LowerTag, DoublePrecision pr) const override {
        return this->_get_as<FloatDPLowerBound>(pr);
    }
    virtual FloatDPApproximation _get(ApproximateTag, DoublePrecision pr) const override {
        return this->_get_as<FloatDPApproximation>(pr);
    }

    // Keep the same paradigm-specific return types and the supplied
    // Preserve the supplied centre/error precision in the native backend.
    virtual FloatMPDPBall _get(MetricTag, MultiplePrecision pr, DoublePrecision pre) const override {
        return this->_get_as<FloatMPDPBall>(pr, pre);
    }
    virtual FloatMPBall _get(MetricTag, MultiplePrecision pr, MultiplePrecision pre) const override {
        return this->_get_as<FloatMPBall>(pr, pre);
    }
    virtual FloatMPBounds _get(OrderTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPBounds>(pr);
    }
    virtual FloatMPUpperBound _get(UpperTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPUpperBound>(pr);
    }
    virtual FloatMPLowerBound _get(LowerTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPLowerBound>(pr);
    }
    virtual FloatMPApproximation _get(ApproximateTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPApproximation>(pr);
    }

    virtual ParadigmCode _paradigm() const override { return P::code(); }
    virtual String _class_name() const override { return class_name<X>(); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << _cast(*this); }

  private:
    template<class R>
    inline R _get_as() const {
        if constexpr (Constructible<R, X>) {
            return static_cast<R>(_cast(*this));
        } else {
            throw ParadigmError();
        }
    }
    template<class R, class PR>
    inline R _get_as(PR pr) const {
        if constexpr (Constructible<R, X, PR>) {
            return R(_cast(*this), pr);
        } else {
            (void)pr;
            throw ParadigmError();
        }
    }
    template<class R, class PR, class PRE>
    inline R _get_as(PR pr, PRE pre) const {
        if constexpr (Constructible<R, X, PR, PRE>) {
            return R(_cast(*this), pr, pre);
        } else {
            (void)pr;
            (void)pre;
            throw ParadigmError();
        }
    }
};

template<class X>
class NumberMixin
    : public ElementaryBinaryNumberDispatcherMixin<X>,
      public ComparisonBinaryNumberDispatcherMixin<X>,
      public NumberGetterMixin<X> {
    const X& _value() const { return NumberGetterMixin<X>::_cast(*this); }

  public:
    NumberInterface* _apply(UnaryElementaryOperator op) const override {
        return _concrete_apply<NumberInterface*>(op, _value());
    }
    NumberInterface* _apply(BinaryElementaryOperator op, const NumberInterface* other) const override {
        return Ariadne::_apply<NumberInterface*, X>(_value(), op, static_cast<const NumberInterface*>(this), other);
    }
    NumberInterface* _rapply(BinaryElementaryOperator op, const NumberInterface* other) const override {
        return Ariadne::_rapply<NumberInterface*, X>(_value(), op, static_cast<const NumberInterface*>(this), other);
    }
    NumberInterface* _apply(GradedElementaryOperator op, Int n) const override {
        return _make_number_wrapper(op(_value(), n));
    }
    LogicalInterface* _apply(BinaryComparisonOperator op, const NumberInterface* other) const override {
        return Ariadne::_apply<LogicalInterface*, X>(_value(), op, static_cast<const NumberInterface*>(this), other);
    }
    LogicalInterface* _rapply(BinaryComparisonOperator op, const NumberInterface* other) const override {
        return Ariadne::_rapply<LogicalInterface*, X>(_value(), op, static_cast<const NumberInterface*>(this), other);
    }

    operator const X&() const { return static_cast<const NumberWrapper<X>&>(*this); }
    operator X&() { return static_cast<NumberWrapper<X>&>(*this); }
};

template<class X>
class NumberWrapper
    : public X,
      public NumberMixin<X> {
    inline static const X& _cast(const NumberWrapper<X>& x) { return static_cast<const X&>(x); }
    static_assert(not Same<X, Handle<NumberInterface>>, "X must be a concrete number, not a handle");
    static_assert(not Same<X, Number<Paradigm<X>>>, "X must be a concrete number, not a generic number");

  public:
    NumberWrapper(const X& a) : X(a) {}
    NumberWrapper(X&& a) : X(std::forward<X>(a)) {}
};

} // namespace Ariadne

#endif /* ARIADNE_NUMBER_WRAPPER_HPP */
