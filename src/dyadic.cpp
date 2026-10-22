/***************************************************************************
 * numeric/dyadic.cpp
 * Native exact dyadic: public interfaces and writers.
 ****************************************************************************/
#include "numeric/dyadic.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatap.hpp"
#include "numeric/rational.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/decimal.hpp"
#include "numeric/double.hpp"
#include "numeric/concepts.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_approximation.hpp"
#include "utility/string.hpp"
#include "utility/writable.hpp"
#include <limits>

namespace Ariadne {

const Dyadic infty=Dyadic::inf();
Rational::Rational(Dyadic const& x) { *this=x.to_rational(); }
Dyadic::Dyadic(FloatDP const& x) : Dyadic(Rational(x)) { }
Dyadic::Dyadic(FloatAP const& x) : Dyadic(Rational(x)) { }

Dyadic operator+(TwoExp y) { return Dyadic(y); }
Dyadic operator-(TwoExp y) { return -Dyadic(y); }
Dyadic operator*(Dyadic x,TwoExp y) { return x*Dyadic(y); }
Dyadic operator/(Dyadic x,TwoExp y) { return x*Dyadic(rec(y)); }
Dyadic operator/(Int n,Two) { return Dyadic(n,1u); }
Dyadic operator/(Integer n,Two) { return Dyadic(std::move(n),1u); }
Dyadic operator/(ExactDouble d,Two) { return hlf(Dyadic(d)); }
Dyadic operator/(Dyadic x,Two) { return hlf(x); }
OutputStream& operator<<(OutputStream& os,TwoExp x) { return os<<"2^"<<x.exponent(); }
Writer<Dyadic> Dyadic::_default_writer(new DecimalWriter());
Void Dyadic::set_default_writer(Writer<Dyadic> writer) { _default_writer=writer; }
Writer<Dyadic> Dyadic::default_writer() { return _default_writer; }
OutputStream& operator<<(OutputStream& os,Dyadic const& x) {
    return os<<Dyadic::default_writer()(x);
}
template<class X> inline OutputStream& write_infinite(OutputStream& os, X const& x) {
    if(is_nan(x)) {
        os << "NaN";
    } else {
        os << (sgn(x)==Sign::POSITIVE ? "" : "-") << "inf";
    }
    return os;
}

auto DecimalWriter::_write(OutputStream& os, Dyadic const& x) const -> OutputStream& {
    if(is_finite(x)) {
        Dyadic w=x;
        if(w<0) { os << "-"; w=-w; }
        Integer z=floor(w);
        os << z << ".";
        w-=z;
        while (w!=0) {
            w*=10;
            z=floor(w);
            w-=z;
            os << z;
        }
    } else {
        write_infinite(os,x);
    }
    return os;
}

auto ScientificWriter::_write(OutputStream& os, Dyadic const& x) const -> OutputStream& {
    if(is_finite(x)) {
        Dyadic w=x;
        if(w<0) { os << "-"; w=-w; }
        int e=0;
        while (w>=2) { w=hlf(w); e=e+1; }
        while (w<1) { w=w*2; e=e-1; }
        os << DecimalWriter()(w);
        os << "*2^" << e;
    } else {
        write_infinite(os,x);
    }
    return os;
}


auto FractionWriter::_write(OutputStream& os, Dyadic const& x) const -> OutputStream& {
    if(is_finite(x)) {
        const Rational q(x);
        os << q.numerator();
        Int exp = log2floor(q.denominator());
        if (exp!=0) { if(exp==1) { os << "/2"; } else { os << "/2^" << exp; } }
    } else {
        write_infinite(os,x);
    }
    return os;
}

auto RepresentationWriter<Dyadic>::_write(OutputStream& os, Dyadic const& x) const -> OutputStream& {
    return os << "Dyadic(" << x.mantissa() << "," << x.exponent() << "u)";
}


Dyadic make_dyadic(unsigned long long int n) { return Dyadic(Integer(n)); }
Dyadic hlf(Integer const& n) { return hlf(Dyadic(n)); }
// Exact dyadic endpoint values with directed rounding to each target.
Bounds<FloatAP> Bounds<Dyadic>::get(ArbitraryPrecision pr) const {
    return Bounds<FloatAP>(FloatAP(lower_raw(),FloatAP::ROUND_DOWNWARD,pr),
                           FloatAP(upper_raw(),FloatAP::ROUND_UPWARD,pr));
}
LowerBound<FloatAP> LowerBound<Dyadic>::get(ArbitraryPrecision pr) const {
    return LowerBound<FloatAP>(FloatAP(raw(),FloatAP::ROUND_DOWNWARD,pr));
}
UpperBound<FloatAP> UpperBound<Dyadic>::get(ArbitraryPrecision pr) const {
    return UpperBound<FloatAP>(FloatAP(raw(),FloatAP::ROUND_UPWARD,pr));
}
Approximation<FloatDP> Approximation<Dyadic>::get(DoublePrecision pr) const {
    return Approximation<FloatDP>(FloatDP(raw(),FloatDP::ROUND_TO_NEAREST,pr));
}
Approximation<FloatAP> Approximation<Dyadic>::get(ArbitraryPrecision pr) const {
    return Approximation<FloatAP>(FloatAP(raw(),FloatAP::ROUND_TO_NEAREST,pr));
}
template<> String class_name<Dyadic>() { return "Dyadic"; }
template<> String class_name<DyadicBounds>() { return "DyadicBounds"; }

static_assert(DyadicRing<Dyadic>);
static_assert(OrderedLattice<Dyadic>);
} // namespace Ariadne
