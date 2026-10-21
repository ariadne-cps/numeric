/***************************************************************************
 * numeric/new_dyadic_extensions.hpp
 * Native Dyadic wrapper API, staged alongside the legacy backend for D48.
 ****************************************************************************/
#ifndef ARIADNE_NEW_DYADIC_EXTENSIONS_HPP
#define ARIADNE_NEW_DYADIC_EXTENSIONS_HPP

#include "numeric/new_dyadic.hpp"
#include "numeric/arithmetic.hpp"
#include "utility/handle.hpp"
#include "utility/writable.hpp"
#include "foundation/logical.hpp"

namespace Ariadne {

inline OutputStream& operator<<(OutputStream& os,NewDyadic const& value) {
    return os << value.literal();
}

using PositiveNewDyadic = Positive<NewDyadic>;

//! \relates NewDyadic
//! \name Type synonyms
//!@{
//! \ingroup NumericModule
using NewDyadicBounds = Bounds<NewDyadic>; //!< Alias for dyadic bounds on a number.
using NewDyadicBall = Ball<NewDyadic,NewDyadic>; //!< Alias for ball about a number with dyadic value and error.
using NewDyadicUpperBound = UpperBound<NewDyadic>; //!< Alias for dyadic upper bound for a number.
using NewDyadicLowerBound = LowerBound<NewDyadic>; //!< Alias for dyadic lower bound for a number.
using NewDyadicApproximation = Approximation<NewDyadic>; //!< Alias for dyadic approximateion to a number.

using PositiveNewDyadicApproximation = Positive<Approximation<NewDyadic>>; //!< <p/>
using PositiveNewDyadicLowerBound = Positive<LowerBound<NewDyadic>>; //!< <p/>
using PositiveNewDyadicUpperBound = Positive<UpperBound<NewDyadic>>; //!< <p/>
using PositiveNewDyadicBounds = Positive<Bounds<NewDyadic>>; //!< <p/>
//!@}

template<> class Bounds<NewDyadic> {
    NewDyadic _l, _u;
  public:
    typedef ValidatedTag Paradigm;
    Bounds(NewDyadic w) : _l(w), _u(w) { }
    template<ConvertibleTo<NewDyadic> X> Bounds(X const& x)
        : Bounds(NewDyadic(x)) { }
    Bounds(NewDyadic l, NewDyadic u) : _l(l), _u(u) { }
    template<class X> requires Constructible<NewDyadic,X> Bounds(Bounds<X> const& x)
        : Bounds(NewDyadic(x.lower_raw()),NewDyadic(x.upper_raw())) { }
    Bounds<NewDyadic> pm(NewDyadic e) { return Bounds<NewDyadic>(_l-e,_u+e); }
    NewDyadic lower() const { return _l; }
    NewDyadic upper() const { return _u; }
    NewDyadic lower_raw() const { return _l; }
    NewDyadic upper_raw() const { return _u; }
    Bounds<FloatDP> get(DoublePrecision pr) const;
    friend Bounds<NewDyadic> operator+(Bounds<NewDyadic> const& w) { return Bounds<NewDyadic>(+w._l,w._u); }
    friend Bounds<NewDyadic> operator-(Bounds<NewDyadic> const& w) { return Bounds<NewDyadic>(-w._u,-w._l); }
    friend Bounds<NewDyadic> operator+(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        return Bounds<NewDyadic>(w1._l+w2._l,w1._u+w2._u); }
    friend Bounds<NewDyadic> operator-(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        return Bounds<NewDyadic>(w1._l-w2._u,w1._u-w2._l); }
    friend Bounds<NewDyadic> operator*(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        // Exact multiplication enclosure: all four vertices, no rounding.
        NewDyadic p=w1._l*w2._l, q=w1._l*w2._u;
        NewDyadic r=w1._u*w2._l, s=w1._u*w2._u;
        return Bounds<NewDyadic>(min(min(p,q),min(r,s)),
                                 max(max(p,q),max(r,s)));
    }
    friend Bounds<Rational> operator/(Bounds<Rational> const& w1, Bounds<Rational> const& w2);
    friend Bounds<NewDyadic> add(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        return Bounds<NewDyadic>(w1._l+w2._l,w1._u+w2._u); }
    friend Bounds<NewDyadic> sub(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        return Bounds<NewDyadic>(w1._l-w2._u,w1._u-w2._l); }
    friend Bounds<NewDyadic> mul(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) { return w1*w2; }
    friend Bounds<Rational> div(Bounds<Rational> const& w1, Bounds<Rational> const& w2);
    friend Bounds<NewDyadic> pow(Bounds<NewDyadic> const& w, Nat m) {
        if (m==0u) { return Bounds<NewDyadic>(NewDyadic(1)); }
        if (m%2u!=0u) {
            return Bounds<NewDyadic>(pow(w._l,m),pow(w._u,m));
        }
        if (w._l>=0) {
            return Bounds<NewDyadic>(pow(w._l,m),pow(w._u,m));
        }
        if (w._u<=0) {
            return Bounds<NewDyadic>(pow(w._u,m),pow(w._l,m));
        }
        return Bounds<NewDyadic>(NewDyadic(0),
                                 max(pow(w._l,m),pow(w._u,m)));
    }
    friend Bounds<NewDyadic> pow(Bounds<NewDyadic> const& w, Int m) {
        if (m<0) { throw std::domain_error("negative NewDyadic interval power"); }
        return pow(w,static_cast<Nat>(m));
    }
    friend Bounds<NewDyadic> nul(Bounds<NewDyadic> const& w) {
        return Bounds<NewDyadic>(nul(w._l),nul(w._u)); }
    friend Bounds<NewDyadic> pos(Bounds<NewDyadic> const& w) {
        return Bounds<NewDyadic>(pos(w._l),pos(w._u)); }
    friend Bounds<NewDyadic> neg(Bounds<NewDyadic> const& w) {
        return Bounds<NewDyadic>(neg(w._u),neg(w._l)); }
    friend Bounds<NewDyadic> sqr(Bounds<NewDyadic> const& w) {
        if(w._l>0) { return Bounds<NewDyadic>(sqr(w._l),sqr(w._u)); }
        else if(w._u<0) { return Bounds<NewDyadic>(sqr(w._u),sqr(w._l)); }
        else { return Bounds<NewDyadic>(nul(w._l),max(sqr(w._l),sqr(w._u))); } }
    friend Bounds<NewDyadic> hlf(Bounds<NewDyadic> const& w) { return Bounds<NewDyadic>(hlf(w._l),hlf(w._u)); }
    friend Bounds<Rational> rec(Bounds<Rational> const& w);
    friend Bounds<NewDyadic> abs(Bounds<NewDyadic> const& w) { return Bounds<NewDyadic>(max(min(w._l,-w._u),0),max(-w._l,w._u)); }
    friend Bounds<NewDyadic> max(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) { return Bounds<NewDyadic>(max(w1._l,w2._l),max(w1._u,w2._u)); }
    friend Bounds<NewDyadic> min(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) { return Bounds<NewDyadic>(min(w1._l,w2._l),min(w1._u,w2._u)); }

    friend ValidatedKleenean sgn(Bounds<NewDyadic> const& w) {
        if (w._l>0) { return true; } else if (w._u<0) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator==(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        if (w1._l>=w2._u && w1._u<=w2._l) { return true; } else if (w1._u< w2._l || w1._l> w2._u) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator!=(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        if (w1._u< w2._l || w1._l> w2._u) { return true; } else if (w1._l>=w2._u && w1._u<=w2._l) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator<=(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        if (w1._u<=w2._l) { return true; } else if (w1._l> w2._u) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator>=(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        if (w1._l>=w2._u) { return true; } else if (w1._u< w2._l) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator< (Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        if (w1._u< w2._l) { return true; } else if (w1._l>=w2._u) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator> (Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        if (w1._l> w2._u) { return true; } else if (w1._u<=w2._l) { return false; } else { return indeterminate; } }
    friend Bool consistent(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        return w1._l<=w2._u && w1._u >= w2._l; }
    friend Boolean inconsistent(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        return w1._l>w2._u || w1._u < w2._l; }
    friend Boolean refines(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        return w1._l>=w2._l and w1._u<=w2._u; }
    friend Bounds<NewDyadic> refinement(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        return Bounds<NewDyadic>(max(w1._l,w2._l),min(w1._u,w2._u)); }
    friend Bounds<NewDyadic> coarsening(Bounds<NewDyadic> const& w1, Bounds<NewDyadic> const& w2) {
        return Bounds<NewDyadic>(min(w1._l,w2._l),max(w1._u,w2._u)); }
    friend OutputStream& operator<<(OutputStream& os, Bounds<NewDyadic> y) { return os << "[" << y._l << ":" << y._u << "]"; }
};


template<> class Ball<NewDyadic,NewDyadic> {
    NewDyadic _v, _e;
  public:
    Ball(NewDyadic w) : _v(w), _e(0) { }
    Ball(NewDyadic v, NewDyadic e) : _v(v), _e(e) { }
    template<class X, class XE> requires Constructible<NewDyadic,X> and Constructible<NewDyadic,XE>
    Ball(Ball<X,XE> const& x)
        : Ball(NewDyadic(x.value_raw()),NewDyadic(x.error_raw())) { }
    explicit Ball(Bounds<NewDyadic> const& w) : _v(hlf(w.lower()+w.upper())), _e(hlf(w.upper()-w.lower())) { }
    explicit operator Bounds<NewDyadic>() const { return Bounds<NewDyadic>(_v-_e,_v+_e); }
    NewDyadic value() const { return _v; }
    NewDyadic error() const { return _e; }
    NewDyadic value_raw() const { return _v; }
    NewDyadic error_raw() const { return _e; }
    friend NewDyadicBall operator+(NewDyadicBall const& w1, NewDyadicBall const& w2) { return NewDyadicBall(w1._v+w2._v,w1._e+w2._e); }
    friend NewDyadicBall operator-(NewDyadicBall const& w1, NewDyadicBall const& w2) { return NewDyadicBall(w1._v-w2._v,w1._e+w2._e); }
    friend NewDyadicBall operator*(NewDyadicBall const& w1, NewDyadicBall const& w2) { return NewDyadicBall(w1._v*w2._v,abs(w1._v)*w2._e+w1._e*abs(w2._v)+w1._e*w2._e); }
    friend NewDyadicBall abs(NewDyadicBall const& w) {
        if (abs(w._v)>=w._e) { return NewDyadicBall(abs(w._v),w._e); } else { NewDyadic av=hlf(max(w._e-w._v,w._v+w._e)); return NewDyadicBall(av,av); } }
    friend ValidatedKleenean operator<(NewDyadicBall const& w1, NewDyadicBall const& w2) {
        if (w1._v+w1._e<w2._v-w2._e) { return true; } else if (w1._v-w1._e >= w2._v+w2._e) { return false; } else { return indeterminate; } }
    friend Boolean inconsistent(NewDyadicBall const& w1, NewDyadicBall const& w2) { return abs(w1._v-w2._v) > w1._e+w2._e; }
    friend Boolean refines(NewDyadicBall const& w1, NewDyadicBall const& w2) { return abs(w1._v-w2._v)+w1._e <= w2._e; }
    friend NewDyadicBall refinement(NewDyadicBall const& w1, NewDyadicBall const& w2) { return NewDyadicBall(refinement(NewDyadicBounds(w1),NewDyadicBounds(w2))); }
    friend NewDyadicBall coarsening(NewDyadicBall const& w1, NewDyadicBall const& w2) { return NewDyadicBall(coarsening(NewDyadicBounds(w1),NewDyadicBounds(w2))); }
    friend OutputStream& operator<<(OutputStream& os, NewDyadicBall y) { return os << "[" << y._v << ":" << y._e << "]"; }
};

template<> class LowerBound<NewDyadic> {
    NewDyadic _l;
  public:
    LowerBound(NewDyadic l) : _l(l) { }
    LowerBound(Bounds<NewDyadic> lu) : _l(lu.lower_raw()) { }
    LowerBound(LowerBound<FloatDP> const& x);
    NewDyadic raw() const { return _l; }
    LowerBound<FloatDP> get(DoublePrecision pr) const;
    friend OutputStream& operator<<(OutputStream& os, LowerBound<NewDyadic> const& y) { return os << y._l << ":"; }
};

template<> class UpperBound<NewDyadic> {
    NewDyadic _u;
  public:
    explicit UpperBound(NewDyadic u) : _u(u) { }
    UpperBound(Bounds<NewDyadic> lu) : _u(lu.upper_raw()) { }
    UpperBound(UpperBound<FloatDP> const& x);
    NewDyadic raw() const { return _u; }
    UpperBound<FloatDP> get(DoublePrecision pr) const;
    friend OutputStream& operator<<(OutputStream& os, UpperBound<NewDyadic> const& y) { return os << ":" << y._u; }
};

template<> class Approximation<NewDyadic> {
    NewDyadic _a;
  public:
    Approximation(NewDyadic a) : _a(a) { }
    Approximation(Approximation<FloatDP> const& x);
    NewDyadic raw() const { return _a; }
    Approximation<FloatDP> get(DoublePrecision pr) const;
    friend OutputStream& operator<<(OutputStream& os, Approximation<NewDyadic> const& y) {
        return os << y._a;
    }
};

template<> class Positive<Bounds<NewDyadic>> : public Bounds<NewDyadic> { public: Positive(Bounds<NewDyadic> w) : Bounds<NewDyadic>(w) { }; };
template<> class Positive<LowerBound<NewDyadic>> : public LowerBound<NewDyadic> { public: Positive(LowerBound<NewDyadic> w) : LowerBound<NewDyadic>(w) { }; };
template<> class Positive<UpperBound<NewDyadic>> : public UpperBound<NewDyadic> { public: Positive(UpperBound<NewDyadic> w) : UpperBound<NewDyadic>(w) { } ; };
template<> class Positive<Approximation<NewDyadic>> : public Approximation<NewDyadic> { public: Positive(Approximation<NewDyadic> w) : Approximation<NewDyadic>(w) { } };


} // namespace Ariadne
#endif // ARIADNE_NEW_DYADIC_EXTENSIONS_HPP
