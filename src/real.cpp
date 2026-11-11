/***************************************************************************
 *            real.cpp
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

#include <limits>
#include <stdexcept>

#include "accuracy.hpp"
#include "concepts.hpp"
#include "decimal.hpp"
#include "dyadic.hpp"
#include "float_bounds.hpp"
#include "floatmp.hpp"
#include "foundation/logical.hpp"
#include "integer.hpp"
#include "number_wrapper.hpp"
#include "numeric/operators.hpp"
#include "rational.hpp"
#include "real_interface.hpp"
#include "reals.hpp"
#include "sequence.hpp"
#include "twoexp.hpp"
#include "utility/module.hpp"

namespace Ariadne {

OutputStream& operator<<(OutputStream& os, const Bits& bits) {
    return os << static_cast<unsigned long int>(bits) << "_bits";
}

Accuracy::Accuracy(Bits precision)
    : _error(1, static_cast<Nat>(precision)) {
}

OutputStream& operator<<(OutputStream& os, const Accuracy& acc) {
    return os << "Accuracy(" << ScientificWriter()(acc.error()) << ")";
}

Bounds<FloatDP> Bounds<Dyadic>::get(DoublePrecision pr) const { return Bounds<FloatDP>(_l, _u, pr); }

template<>
Ball<FloatDP>::Ball(const Real& r, DoublePrecision pr) : Ball(r.compute_get(Effort(53), pr)) {}
template<>
Bounds<FloatDP>::Bounds(const Real& r, DoublePrecision pr) : Bounds<FloatDP>(r.compute_get(Effort(53), pr)) {}
template<>
UpperBound<FloatDP>::UpperBound(const Real& r, DoublePrecision pr) : UpperBound<FloatDP>(r.upper().compute_get(Effort(53), pr)) {}
template<>
LowerBound<FloatDP>::LowerBound(const Real& r, DoublePrecision pr) : LowerBound<FloatDP>(r.lower().compute_get(Effort(53), pr)) {}
template<>
Approximation<FloatDP>::Approximation(const Real& r, DoublePrecision pr) : Approximation<FloatDP>(r.compute_get(Effort(53), pr)) {}

template<class X>
struct ValidatedRealWrapper;

template<>
struct ValidatedRealWrapper<DyadicBounds> : public ValidatedRealInterface,
                                            public DyadicBounds {
    typedef DyadicBounds X;
    ValidatedRealWrapper(const X& x) : X(x) {}
    virtual DyadicBounds _get() const override final { return *this; }
    virtual FloatDPBounds _get(DoublePrecision pr) const override final { return FloatDPBounds(this->_get(), pr); }
    virtual FloatMPBounds _get(MultiplePrecision pr) const override final { return FloatMPBounds(this->_get(), pr); }
    virtual OutputStream& _write(OutputStream& os) const override final { return os << static_cast<const DyadicBounds&>(*this); }
};

struct RealBase : RealInterface {
    virtual ValidatedReal _compute(Accuracy acc) const override;
    virtual ValidatedReal _compute(Effort eff) const override;
};

ValidatedReal RealBase::_compute(Effort effort) const {
    return ValidatedReal(this->_compute_get(effort));
}

ValidatedReal RealBase::_compute(Accuracy accuracy) const {
    Nat effort = 1;
    Dyadic error_bound = accuracy.error();
    Dyadic error = 2 * error_bound;
    DyadicBounds res(0);
    while (!(error < error_bound)) {
        res = this->_compute_get(Effort(effort));
        error = hlf(res.upper_raw() - res.lower_raw());
        effort += 1;
    }
    return ValidatedReal(res);
}

MultiplePrecision precision(Nat pr) { return MultiplePrecision(pr); }

static const Nat PRECISION_PER_UNIT_EFFORT = 16;

struct RealExpressionBase : RealBase {
    virtual DyadicBounds _compute_get(Effort eff) const override {
        if (eff.work() * PRECISION_PER_UNIT_EFFORT <= 52) {
            return DyadicBounds(this->_compute_get(eff, double_precision));
        } else {
            return DyadicBounds(this->_compute_get(eff, precision(eff.work() * PRECISION_PER_UNIT_EFFORT)));
        }
    }
    virtual FloatDPBounds _compute_get(Effort eff, DoublePrecision pr) const override = 0;
    virtual FloatMPBounds _compute_get(Effort eff, MultiplePrecision pr) const override = 0;
};

template<class O, class A>
struct UnaryNumericExpression {
    O _op;
    A _arg;
    UnaryNumericExpression(O o, A a) : _op(o), _arg(a) {}
    OutputStream& _write_expression(OutputStream& os) const {
        return os << _op.code() << "(" << _arg << ")";
    }
};

template<class O, class A1, class A2>
struct BinaryNumericExpression {
    O _op;
    A1 _arg1;
    A2 _arg2;
    BinaryNumericExpression(O o, A1 a1, A2 a2) : _op(o), _arg1(a1), _arg2(a2) {}
    OutputStream& _write_expression(OutputStream& os) const {
        return os << _op.code() << "(" << _arg1 << "," << _arg2 << ")";
    }
};

template<class O, class A, class N>
struct GradedNumericExpression {
    O _op;
    A _arg;
    N _num;
    GradedNumericExpression(O o, A a, N n) : _op(o), _arg(a), _num(n) {}
    OutputStream& _write_expression(OutputStream& os) const {
        return os << _op.code() << "(" << _arg << "," << _num << ")";
    }
};

template<class O, class... AS>
struct RealWrapper;

template<class O, class A>
struct RealWrapper<O, A>
    : virtual RealExpressionBase,
      UnaryNumericExpression<O, A>,
      FloatDPBounds {
    using Expression = UnaryNumericExpression<O, A>;
    RealWrapper(O o, A a) : Expression(o, a), FloatDPBounds(this->_op(this->_arg.get(dp))) {}
    virtual FloatDPBounds _compute_get(Effort, DoublePrecision) const { return static_cast<FloatDPBounds>(*this); }
    virtual FloatMPBounds _compute_get(Effort, MultiplePrecision pr) const { return this->_op(this->_arg.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return this->_write_expression(os); }
};

template<class O, class A1, class A2>
struct RealWrapper<O, A1, A2>
    : virtual RealExpressionBase,
      BinaryNumericExpression<O, A1, A2>,
      FloatDPBounds {
    using Expression = BinaryNumericExpression<O, A1, A2>;
    RealWrapper(O o, A1 a1, A2 a2) : Expression(o, a1, a2), FloatDPBounds(this->_op(this->_arg1.get(dp), this->_arg2.get(dp))) {}
    virtual FloatDPBounds _compute_get(Effort, DoublePrecision) const { return static_cast<FloatDPBounds>(*this); }
    virtual FloatMPBounds _compute_get(Effort eff, MultiplePrecision pr) const { return this->_op(this->_arg1.compute_get(eff, pr), this->_arg2.compute_get(eff, pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return this->_write_expression(os); }
};

template<class A, class N>
struct RealWrapper<Pow, A, N>
    : virtual RealExpressionBase,
      GradedNumericExpression<Pow, A, N>,
      FloatDPBounds {
    using Expression = GradedNumericExpression<Pow, A, N>;
    RealWrapper(Pow o, A a, N n) : Expression(o, a, n), FloatDPBounds(this->_op(this->_arg.get(dp), n)) {}
    virtual FloatDPBounds _compute_get(Effort, DoublePrecision) const { return static_cast<FloatDPBounds>(*this); }
    virtual FloatMPBounds _compute_get(Effort eff, MultiplePrecision pr) const { return this->_op(this->_arg.compute_get(eff, pr), this->_num); }
    virtual OutputStream& _write(OutputStream& os) const { return this->_write_expression(os); }
};

template<class X>
struct RealWrapper<Cnst, X> : RealExpressionBase,
                              FloatDPBounds {
    X _c;

  public:
    RealWrapper(const X& x) : FloatDPBounds(x, dp), _c(x) {}
    virtual FloatDPBounds _compute_get(Effort, DoublePrecision) const { return static_cast<const FloatDPBounds&>(*this); }
    virtual FloatMPBounds _compute_get(Effort, MultiplePrecision pr) const { return FloatMPBounds(this->_c, pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<>
struct RealWrapper<Cnst, EffectiveNumber> : RealExpressionBase,
                                            FloatDPBounds {
    typedef EffectiveNumber X;
    X _c;

  public:
    RealWrapper(const X& x) : FloatDPBounds(x, dp), _c(x) {}
    virtual FloatDPBounds _compute_get(Effort, DoublePrecision) const { return static_cast<const FloatDPBounds&>(*this); }
    virtual FloatMPBounds _compute_get(Effort, MultiplePrecision pr) const { return FloatMPBounds(this->_c, pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_c; }
};

template<class Y>
struct RealLimit;

template<>
struct RealLimit<Real> : RealBase {
    Sequence<Real> _seq;

  public:
    RealLimit(const Sequence<Real>& seq) : _seq(seq) {}
    virtual DyadicBounds _compute_get(Effort eff) const {
        Nat n = eff.work() + 1u;
        return _seq[n].compute_get(eff).pm(Dyadic(TwoExp(-static_cast<Int>(n))));
    }
    virtual FloatDPBounds _compute_get(Effort eff, DoublePrecision pr) const {
        Nat n = eff.work() + 1u;
        return _seq[n].compute_get(eff, pr) + (FloatDPBounds(-(TwoExp(-static_cast<Int>(n))), +TwoExp(-static_cast<Int>(n)), pr));
    }
    virtual FloatMPBounds _compute_get(Effort eff, MultiplePrecision pr) const {
        Nat n = eff.work() + 1u;
        return _seq[n].compute_get(eff, pr) + FloatMPBounds(-(TwoExp(-static_cast<Int>(n))), +TwoExp(-static_cast<Int>(n)), pr);
    }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << "{" << _seq[0u] << ", " << _seq[1u] << ", " << _seq[2u] << ", ... }";
    }
};

template<>
struct RealLimit<Dyadic> : RealBase {
    Sequence<Dyadic> _seq;

  public:
    RealLimit(const Sequence<Dyadic>& seq) : _seq(seq) {}
    virtual DyadicBounds _compute_get(Effort eff) const {
        Nat n = eff.work() + 1u;
        return DyadicBounds(_seq[n] - (TwoExp(-static_cast<Int>(n))), _seq[n] + (TwoExp(-static_cast<Int>(n))));
    }
    virtual FloatDPBounds _compute_get(Effort eff, DoublePrecision pr) const {
        return this->_compute(eff).get(pr);
    }
    virtual FloatMPBounds _compute_get(Effort eff, MultiplePrecision pr) const {
        return this->_compute(eff).get(pr);
    }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << "{" << _seq[0u] << ", " << _seq[1u] << ", " << _seq[2u] << ", ... }";
    }
};

template<>
struct RealLimit<DyadicBounds> : RealBase {
    Sequence<DyadicBounds> _seq;

  public:
    RealLimit(const ConvergentSequence<DyadicBounds>& seq) : _seq(seq) {}
    virtual DyadicBounds _compute_get(Effort eff) const {
        return _seq[eff.work()];
    }
    virtual FloatDPBounds _compute_get(Effort eff, DoublePrecision pr) const {
        return _seq[eff.work()].get(pr);
    }
    virtual FloatMPBounds _compute_get(Effort eff, MultiplePrecision pr) const {
        return _seq[eff.work()].get(pr);
    }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << "{" << _seq[0u] << ", " << _seq[1u] << ", " << _seq[2u] << ", ... }";
    }
};

template<class O, class... A>
inline Real make_real(O o, A... a) {
    return Real(std::make_shared<RealWrapper<O, A...>>(o, a...));
}

Real::Real(SharedPointer<const Interface> p) : Handle<const Interface>(p) {}

Real::Real(const ConvergentSequence<DyadicBounds>& seq) : Real(std::make_shared<RealLimit<DyadicBounds>>(seq)) {}
Real::Real(const FastCauchySequence<Dyadic>& seq) : Real(std::make_shared<RealLimit<Dyadic>>(seq)) {}

UpperReal Real::upper() const { return UpperReal(this->_ptr); }
LowerReal Real::lower() const { return LowerReal(this->_ptr); }

double Real::get_d() const { return this->compute_get(Effort(53), dp).get_d(); }

/*
template<class PR> FloatBall<PR>::Float(Real const& x) : FloatBall<PR>(x.lower(),x.upper()) { }
template<class PR> FloatBounds<PR>::Float(Real const& x) : FloatBounds<PR>(x.lower(),x.upper()) { }
template<class PR> FloatUpperBound<PR>::Float(Real const& x) : FloatUpperBound<PR>(x.upper()) { }
template<class PR> FloatLowerBound<PR>::Float(Real const& x) : FloatLowerBound<PR>(x.lower()) { }
template<class PR> FloatApproximation<PR>::Float(Real const& x) : FloatApproximation<PR>(x.approx()) { }
*/

Real::Real(std::uint64_t m, Void*) : Real(Integer(m)) {}
Real::Real(std::int64_t n, Void*) : Real(Integer(n)) {}

Real::Real() : Real(Integer(0)) {}
NaiveReal::NaiveReal() : Handle<const Interface>(Real().managed_pointer()) {}
Real::Real(ExactDouble d) : Real(std::make_shared<RealWrapper<Cnst, ExactDouble>>(d)) {}
// Real::Real(ExactDouble d) : Real(Dyadic(d)) { }
Real::Real(const Integer& z) : Real(std::make_shared<RealWrapper<Cnst, Integer>>(z)) {}
Real::Real(const Dyadic& w) : Real(std::make_shared<RealWrapper<Cnst, Dyadic>>(w)) {}
Real::Real(const Decimal& d) : Real(std::make_shared<RealWrapper<Cnst, Decimal>>(d)) {}
Real::Real(const Rational& q) : Real(std::make_shared<RealWrapper<Cnst, Rational>>(q)) {}
Real::Real(EffectiveNumber q) : Real(std::make_shared<RealWrapper<Cnst, EffectiveNumber>>(q)) {}

Real add(const Real& x1, const Real& x2) { return make_real(Add(), x1, x2); }
Real sub(const Real& x1, const Real& x2) { return make_real(Sub(), x1, x2); }
Real mul(const Real& x1, const Real& x2) { return make_real(Mul(), x1, x2); }
Real div(const Real& x1, const Real& x2) { return make_real(Div(), x1, x2); }
Real pow(const Real& x1, Nat m2) { return make_real(Pow(), x1, Int(m2)); }
Real pow(const Real& x1, Int n2) { return make_real(Pow(), x1, n2); }
Real nul(const Real&) { return Real(0); }
Real pos(const Real& x) { return make_real(Pos(), x); }
Real neg(const Real& x) { return make_real(Neg(), x); }
Real hlf(const Real& x) { return make_real(Hlf(), x); }
Real sqr(const Real& x) { return make_real(Sqr(), x); }
Real rec(const Real& x) { return make_real(Rec(), x); }
Real sqrt(const Real& x) { return make_real(Sqrt(), x); }
Real exp(const Real& x) { return make_real(Exp(), x); }
Real log(const Real& x) { return make_real(Log(), x); }
Real sin(const Real& x) { return make_real(Sin(), x); }
Real cos(const Real& x) { return make_real(Cos(), x); }
Real tan(const Real& x) { return make_real(Tan(), x); }
Real tanh(const Real& x) { return make_real(Tanh(), x); }
Real asin(const Real& x) { return make_real(Asin(), x); }
Real acos(const Real& x) { return make_real(Acos(), x); }
Real atan(const Real& x) { return make_real(Atan(), x); }

PositiveReal abs(const Real& x) { return PositiveReal(make_real(Abs(), x)); }
Real max(const Real& x1, const Real& x2) { return make_real(Max(), x1, x2); }
Real min(const Real& x1, const Real& x2) { return make_real(Min(), x1, x2); }

Real limit(const FastCauchySequence<Real>& seq) { return Real(std::make_shared<RealLimit<Real>>(seq)); }
Real limit(const FastCauchySequence<Dyadic>& seq) { return Real(std::make_shared<RealLimit<Dyadic>>(seq)); }

Real choose(const Case<LowerKleenean, Real>& c1, const Case<LowerKleenean, Real>& c2) {
    if (choose(c1.condition(), c2.condition())) {
        return c1.term();
    } else {
        return c2.term();
    }
}

class WhenRealExpression : public RealBase {
    UpperKleenean _p1, _p2;
    Real _r1, _r2;

  public:
    WhenRealExpression(const UpperKleenean& p1, const Real& r1, const UpperKleenean& p2, const Real& r2) : _p1(p1), _p2(p2), _r1(r1), _r2(r2) {}
    WhenRealExpression(const Case<UpperKleenean, Real>& c1, const Case<UpperKleenean, Real>& c2)
        : WhenRealExpression(c1.condition(), c1.term(), c2.condition(), c2.term()) {}
    virtual DyadicBounds _compute_get(Effort eff) const;
    virtual FloatDPBounds _compute_get(Effort, DoublePrecision pr) const { return this->_compute(Effort(0u)).get(pr); }
    virtual FloatMPBounds _compute_get(Effort, MultiplePrecision pr) const { return this->_compute(Effort(static_cast<Nat>(pr.bits()))).get(pr); }
    friend OutputStream& operator<<(OutputStream& os, const WhenRealExpression& r) { return r._write(os); }

  public:
    virtual OutputStream& _write(OutputStream& os) const;
};
Real when(const Case<UpperKleenean, Real>& c1, const Case<UpperKleenean, Real>& c2) {
    return Real(std::make_shared<WhenRealExpression>(c1, c2));
}

DyadicBounds WhenRealExpression::_compute_get(Effort eff) const {
    ValidatedUpperKleenean cp1 = _p1.check(eff);
    ValidatedUpperKleenean cp2 = _p2.check(eff);

    ARIADNE_ASSERT_MSG(possibly(cp1) or possibly(cp2), "Unsatisfiable when-expression " << *this);
    if (not possibly(cp1)) {
        return _r2.compute_get(eff);
    }
    if (not possibly(cp2)) {
        return _r1.compute_get(eff);
    }

    ValidatedReal vr1 = _r1.compute(eff);
    ValidatedReal vr2 = _r2.compute(eff);

    DyadicBounds w1 = vr1.get();
    DyadicBounds w2 = vr2.get();

    return coarsening(w1, w2);
}

OutputStream& WhenRealExpression::_write(OutputStream& os) const {
    return os << "when(" << _p1 << " => " << _r1 << " & " << _p2 << " => " << _r2 << ")";
}

// May be true if r>a; may be false if r<b
// Equivalent to choose(r>a,r<b)
Boolean nondeterministic_greater(const Real& r, const Rational& a, const Rational& b) {
    ARIADNE_PRECONDITION(a < b);
    Effort eff(0);
refine:
    DyadicBounds x = r.compute_get(eff);
    if (x.lower_raw() > a) {
        return true;
    } else if (x.upper_raw() < b) {
        return false;
    }
    ++eff;
    goto refine;
}

PositiveUpperReal mag(const Real& r) { return abs(r); }
FloatDPError mag(const Real& r, DoublePrecision pr) { return mag(r.compute_get(Effort(53), pr)); }

OutputStream& operator<<(OutputStream& os, const Real& x) { return x._ptr->_write(os); }

Bool same(const Real& r1, const Real& r2) {
    // FIXME: Use symbolic approach
    DoublePrecision pr;
    FloatDPBounds x1(r1, pr);
    FloatDPBounds x2(r2, pr);
    return x1.lower_raw() == x2.upper_raw() && x1.upper_raw() == x2.lower_raw();
}

PositiveReal dist(const Real& r1, const Real& r2) { return abs(sub(r1, r2)); }

template<class O, class... ARGS>
struct LogicalWrapper;

template<class O>
struct LogicalWrapper<O, Real>
    : virtual LogicalInterface,
      UnaryNumericExpression<O, Real> {
    using Expression = UnaryNumericExpression<O, Real>;
    LogicalWrapper(O o, Real a) : Expression(o, a) {}
    virtual LogicalInterface* _copy() const;
    virtual LogicalValue _check(Effort e) const;
    virtual OutputStream& _write(OutputStream& os) const {
        return this->_write_expression(os);
    }
};

template<class O>
LogicalInterface* LogicalWrapper<O, Real>::_copy() const {
    return new LogicalWrapper<O, Real>(*this);
}

template<class O>
LogicalValue LogicalWrapper<O, Real>::_check(Effort e) const {
    return static_cast<LogicalValue>(this->_op(this->_arg.compute(e)));
    //    if(e==0u) { DoublePrecision p; return static_cast<LogicalValue>(this->_op(this->_arg.get(p))); }
    //    else { MultiplePrecision p(e*64); return static_cast<LogicalValue>(this->_op(this->_arg.get(p))); }
}

template<class O>
struct LogicalWrapper<O, Real, Real>
    : virtual LogicalInterface,
      BinaryNumericExpression<O, Real, Real> {
    using Expression = BinaryNumericExpression<O, Real, Real>;
    LogicalWrapper(O o, Real a1, Real a2) : Expression(o, a1, a2) {}
    virtual LogicalInterface* _copy() const;
    virtual LogicalValue _check(Effort e) const;
    virtual OutputStream& _write(OutputStream& os) const {
        return this->_write_expression(os);
    }
};

template<class O>
LogicalInterface* LogicalWrapper<O, Real, Real>::_copy() const {
    return new LogicalWrapper<O, Real, Real>(*this);
}
template<class O>
LogicalValue LogicalWrapper<O, Real, Real>::_check(Effort e) const {
    return static_cast<LogicalValue>(this->_op(this->_arg1.compute(e), this->_arg2.compute(e)));
    //    if(e==0u) { DoublePrecision p; return static_cast<LogicalValue>(this->_op(this->_arg1.get(p),this->_arg2.get(p))); }
    //    else { MultiplePrecision p(e*64); return static_cast<LogicalValue>(this->_op(this->_arg1.get(p),this->_arg2.get(p))); }
}

template<class R, class O, class... ARGS>
R make_logical(O op, ARGS... args) {
    return R(LogicalHandle(std::make_shared<LogicalWrapper<O, ARGS...>>(op, args...)));
}
template<class P, class O, class... ARGS>
LogicalType<P> make_logical_type(O op, ARGS... args) {
    return LogicalType<P>(std::make_shared<LogicalWrapper<O, ARGS...>>(op, args...));
}

NegatedSierpinskian eq(const Real& x1, const Real& x2) { return make_logical<NegatedSierpinskian>(Equal(), x1, x2); }
Kleenean lt(const Real& x1, const Real& x2) { return make_logical<Kleenean>(Less(), x1, x2); }

Falsifyable operator==(const Real& x1, const Real& x2) { return make_logical<NegatedSierpinskian>(Equal(), x1, x2); }
Verifyable operator!=(const Real& x1, const Real& x2) { return make_logical<Sierpinskian>(Unequal(), x1, x2); }
Quasidecidable operator<(const Real& x1, const Real& x2) { return make_logical<Kleenean>(Less(), x1, x2); }
Quasidecidable operator>(const Real& x1, const Real& x2) { return make_logical<Kleenean>(Gtr(), x1, x2); }
Quasidecidable operator<=(const Real& x1, const Real& x2) { return make_logical<Kleenean>(Leq(), x1, x2); }
Quasidecidable operator>=(const Real& x1, const Real& x2) { return make_logical<Kleenean>(Geq(), x1, x2); }

Kleenean sgn(const Real& x) { return make_logical<Kleenean>(Sgn(), x); }

Integer round(const Real& r) {
    DyadicBounds wb = r.compute(Accuracy(1_bits)).get();
    return round(hlf(wb.lower() + wb.upper()));
}

DyadicBounds Real::compute_get(Effort e) const { return this->_ptr->_compute(e); }
FloatDPBounds Real::get(DoublePrecision pr) const { return this->_ptr->_compute_get(Effort(0), pr); }
FloatMPBounds Real::get(MultiplePrecision pr) const {
    const std::size_t bits = pr.bits();
    if (bits > std::numeric_limits<Nat>::max()) {
        throw std::overflow_error("Real: precision exceeds the Effort range");
    }
    return this->compute_get(Effort(static_cast<Nat>(bits)), pr);
}
FloatDPBounds Real::compute_using(DoublePrecision pr) const { return this->_ptr->_compute_get(Effort(0), pr); }
FloatMPBounds Real::compute_using(MultiplePrecision pr) const {
    return this->get(pr);
}

template<>
String class_name<Real>() { return "Real"; }
template<>
String class_name<PositiveReal>() { return "PositiveReal"; }

const Real pi = 4 * atan(1); // Real(3.1415926535897930, 3.141592653589793238, 3.1415926535897936);
const Real infinity = Real(operator""_x(std::numeric_limits<double>::infinity()));

ValidatedReal Real::compute(Effort eff) const {
    return this->_ptr->_compute(eff);
}

ValidatedReal Real::compute(Accuracy accuracy) const {
    return this->_ptr->_compute(accuracy);
}

FloatDPBounds Real::compute_get(Effort eff, DoublePrecision pr) const {
    return this->_ptr->_compute_get(eff, pr);
}

FloatMPBounds Real::compute_get(Effort eff, MultiplePrecision pr) const {
    return this->_ptr->_compute_get(eff, pr);
}

ValidatedKleenean check_sgn(Real r, Effort eff) {
    auto x = r.compute_get(eff);
    if (definitely(x.lower_raw() > 0)) {
        return true;
    } else if (definitely(x.upper_raw() < 0)) {
        return false;
    } else {
        return indeterminate;
    }
}

LowerReal::LowerReal(SharedPointer<const Interface> p) : Handle<const Interface>(p) {
}

LowerReal::LowerReal(Real r) : LowerReal(r.managed_pointer()) {
}

ValidatedLowerReal LowerReal::compute(Effort eff) const {
    return ValidatedLowerReal(this->_ptr->_compute(eff));
}

DyadicLowerBound LowerReal::compute_get(Effort eff) const {
    return this->_ptr->_compute(eff).get();
}

FloatDPLowerBound LowerReal::compute_get(Effort eff, DoublePrecision pr) const {
    return this->_ptr->_compute(eff).get(pr);
}

FloatMPLowerBound LowerReal::compute_get(Effort eff, MultiplePrecision pr) const {
    return this->_ptr->_compute(eff).get(pr);
}

OutputStream& operator<<(OutputStream& os, const LowerReal& x) {
    return x.pointer()->_write(os);
}

UpperReal::UpperReal(SharedPointer<const Interface> p) : Handle<const Interface>(p) {
}

UpperReal::UpperReal(Real r) : UpperReal(r.managed_pointer()) {
}

ValidatedUpperReal UpperReal::compute(Effort eff) const {
    return this->_ptr->_compute(eff);
}

DyadicUpperBound UpperReal::compute_get(Effort eff) const {
    return this->_ptr->_compute(eff).get();
}

FloatDPUpperBound UpperReal::compute_get(Effort eff, DoublePrecision pr) const {
    return this->_ptr->_compute(eff).get(pr);
}

FloatMPUpperBound UpperReal::compute_get(Effort eff, MultiplePrecision pr) const {
    return this->_ptr->_compute(eff).get(pr);
}

OutputStream& operator<<(OutputStream& os, const UpperReal& x) {
    return x.pointer()->_write(os);
}

inline const Real& cast_real(const LowerReal& lr) { return reinterpret_cast<const Real&>(lr); }
inline const Real& cast_real(const UpperReal& ur) { return reinterpret_cast<const Real&>(ur); }
inline const Real& make_signed(const PositiveReal& pr) { return pr; }
inline const LowerReal& make_lower(const Real& r) { return reinterpret_cast<const LowerReal&>(r); }
inline const UpperReal& make_upper(const Real& r) { return reinterpret_cast<const UpperReal&>(r); }

PositiveReal cast_positive(const Real& pr) { return static_cast<const PositiveReal&>(pr); }

LowerReal max(const LowerReal& lr1, const LowerReal& lr2) { return make_lower(max(cast_real(lr1), cast_real(lr2))); }
LowerReal min(const LowerReal& lr1, const LowerReal& lr2) { return make_lower(min(cast_real(lr1), cast_real(lr2))); }
Real min(const LowerReal& lr1, const Real& r2) { return min(cast_real(lr1), r2); }
Real min(const Real& r1, const LowerReal& lr2) { return min(r1, cast_real(lr2)); }

UpperReal max(const UpperReal& ur1, const UpperReal& ur2) { return make_upper(max(cast_real(ur1), cast_real(ur2))); }
Real max(const UpperReal& ur1, const Real& r2) { return max(cast_real(ur1), r2); }
Real max(Real r1, const UpperReal& ur2) { return max(r1, cast_real(ur2)); }
UpperReal min(const UpperReal& ur1, const UpperReal& ur2) { return make_upper(min(cast_real(ur1), cast_real(ur2))); }

LowerReal neg(const UpperReal& ur) { return make_lower(neg(cast_real(ur))); }
UpperReal neg(const LowerReal& lr) { return make_upper(neg(cast_real(lr))); }
LowerReal add(const LowerReal& lr1, const LowerReal& lr2) { return make_lower(add(cast_real(lr1), cast_real(lr2))); }
UpperReal add(const UpperReal& ur1, const UpperReal& ur2) { return make_upper(add(cast_real(ur1), cast_real(ur2))); }

PositiveDyadicBounds PositiveReal::compute_get(Effort eff) const {
    return PositiveDyadicBounds(this->_ptr->_compute(eff).get());
}

PositiveDyadicLowerBound PositiveLowerReal::compute_get(Effort eff) const {
    return PositiveDyadicLowerBound(this->_ptr->_compute(eff).get());
}

PositiveDyadicUpperBound PositiveUpperReal::compute_get(Effort eff) const {
    return PositiveDyadicUpperBound(this->_ptr->_compute(eff).get());
}

PositiveReal max(const PositiveReal& pr1, const PositiveReal& pr2) { return cast_positive(max(make_signed(pr1), make_signed(pr2))); }
PositiveReal min(const PositiveReal& pr1, const PositiveReal& pr2) { return cast_positive(min(make_signed(pr1), make_signed(pr2))); }
PositiveReal add(const PositiveReal& pr1, const PositiveReal& pr2) { return cast_positive(add(make_signed(pr1), make_signed(pr2))); }
PositiveReal mul(const PositiveReal& pr1, const PositiveReal& pr2) { return cast_positive(mul(make_signed(pr1), make_signed(pr2))); }
PositiveReal div(const PositiveReal& pr1, const PositiveReal& pr2) { return cast_positive(div(make_signed(pr1), make_signed(pr2))); }
PositiveReal rec(const PositiveReal& pr) { return cast_positive(rec(make_signed(pr))); }

PositiveUpperReal rec(const PositiveLowerReal& plr) { return cast_positive(rec(cast_real(plr))); }
PositiveLowerReal rec(const PositiveUpperReal& pur) { return cast_positive(rec(cast_real(pur))); }
PositiveLowerReal add(const PositiveLowerReal& plr1, const PositiveLowerReal& plr2) { return cast_positive(add(cast_real(plr1), cast_real(plr2))); }
PositiveUpperReal add(const PositiveUpperReal& pur1, const PositiveUpperReal& pur2) { return cast_positive(add(cast_real(pur1), cast_real(pur2))); }
PositiveLowerReal mul(const PositiveLowerReal& plr1, const PositiveLowerReal& plr2) { return cast_positive(mul(cast_real(plr1), cast_real(plr2))); }
PositiveUpperReal mul(const PositiveUpperReal& pur1, const PositiveUpperReal& pur2) { return cast_positive(mul(cast_real(pur1), cast_real(pur2))); }
PositiveLowerReal div(const PositiveLowerReal& plr1, const PositiveUpperReal& pur2) { return cast_positive(div(cast_real(plr1), cast_real(pur2))); }
PositiveUpperReal div(const PositiveUpperReal& pur1, const PositiveLowerReal& plr2) { return cast_positive(div(cast_real(pur1), cast_real(plr2))); }

LowerReal mul(const LowerReal& lr1, const PositiveReal& pr2) { return mul(cast_real(lr1), make_signed(pr2)); }
UpperReal mul(const UpperReal& ur1, const PositiveReal& pr2) { return mul(cast_real(ur1), make_signed(pr2)); }
LowerReal mul(const PositiveReal& pr1, const LowerReal& lr2) { return mul(make_signed(pr1), cast_real(lr2)); }
UpperReal mul(const PositiveReal& pr1, const UpperReal& ur2) { return mul(make_signed(pr1), cast_real(ur2)); }
LowerReal div(const LowerReal& lr1, const PositiveReal& pr2) { return div(cast_real(lr1), make_signed(pr2)); }
UpperReal div(const UpperReal& ur1, const PositiveReal& pr2) { return div(cast_real(ur1), make_signed(pr2)); }
LowerReal div(const PositiveReal& pr1, const UpperReal& ur2) { return div(make_signed(pr1), cast_real(ur2)); }
UpperReal div(const PositiveReal& pr1, const LowerReal& lr2) { return div(make_signed(pr1), cast_real(lr2)); }

static_assert(Constructible<FloatDP, Dyadic, FloatDP::RoundingModeType, FloatDP::PrecisionType>);
static_assert(Constructible<FloatMP, Dyadic, FloatMP::RoundingModeType, FloatMP::PrecisionType>);

ValidatedReal::ValidatedReal(const DyadicBounds& y)
    : Handle<const Interface>(std::make_shared<ValidatedRealWrapper<DyadicBounds>>(y)) {}
ValidatedReal::operator DyadicBounds() const { return this->get(); }
DyadicBounds ValidatedReal::get() const { return this->_ptr->_get(); }
FloatDPBounds ValidatedReal::get(DoublePrecision pr) const { return this->_ptr->_get(pr); }
FloatMPBounds ValidatedReal::get(MultiplePrecision pr) const {
    return this->_ptr->_get(pr);
}

OutputStream& operator<<(OutputStream& os, const ValidatedReal& vr) { return vr._ptr->_write(os); }

ValidatedLowerReal::ValidatedLowerReal(const DyadicLowerBound& y)
    : Handle<const Interface>(std::make_shared<ValidatedRealWrapper<DyadicBounds>>(y.raw())) {}
ValidatedLowerReal::ValidatedLowerReal(const ValidatedReal& r) : Handle<const Interface>(r.managed_pointer()) {}
DyadicLowerBound ValidatedLowerReal::get() const { return this->_ptr->_get(); }
FloatDPLowerBound ValidatedLowerReal::get(DoublePrecision pr) const { return FloatDPLowerBound(this->get().raw(), pr); }
FloatMPLowerBound ValidatedLowerReal::get(MultiplePrecision pr) const {
    return FloatMPLowerBound(this->get().raw(), pr);
}
OutputStream& operator<<(OutputStream& os, const ValidatedLowerReal& lr) { return lr._ptr->_write(os); }

ValidatedUpperReal::ValidatedUpperReal(const DyadicUpperBound& y)
    : Handle<const Interface>(std::make_shared<ValidatedRealWrapper<DyadicBounds>>(y.raw())) {}
ValidatedUpperReal::ValidatedUpperReal(const ValidatedReal& r) : Handle<const Interface>(r.managed_pointer()) {}
DyadicUpperBound ValidatedUpperReal::get() const { return this->_ptr->_get(); }
FloatDPUpperBound ValidatedUpperReal::get(DoublePrecision pr) const { return FloatDPUpperBound(this->get().raw(), pr); }
FloatMPUpperBound ValidatedUpperReal::get(MultiplePrecision pr) const {
    return FloatMPUpperBound(this->get().raw(), pr);
}
OutputStream& operator<<(OutputStream& os, const ValidatedUpperReal& ur) { return ur._ptr->_write(os); }

ApproximateDouble::ApproximateDouble(const Real& r) : _d(FloatDPApproximation(r, dp).raw().get_d()) {}

static_assert(TranscendentalField<Real>);
static_assert(OrderedLattice<Real>);

} // namespace Ariadne
