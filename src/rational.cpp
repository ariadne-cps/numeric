/***************************************************************************
 *            numeric/rational.cpp
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

/*! \file numeric/rational.cpp
 *  \brief
 */



#include "utility/stdlib.hpp"
#include "macro/general.hpp"
#include "utility/typedefs.hpp"

#include "rational.hpp"
#include "foundation/logical.hpp"
#include "number.hpp"
#include "builtin.hpp"
#include "integer.hpp"
#include "dyadic.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "sign.hpp"
#include "extended.hpp"

#include "concepts.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace Ariadne {

class InvalidRationalLiteralException {
  public:
    InvalidRationalLiteralException(StringType) { }
};

namespace {
Sign integer_sign(Integer const& n) { return sgn(n); }
Integer integer_add(Integer const& a,Integer const& b) { return add(a,b); }
Integer integer_sub(Integer const& a,Integer const& b) { return sub(a,b); }
Integer integer_mul(Integer const& a,Integer const& b) { return mul(a,b); }
Integer integer_neg(Integer const& n) { return neg(n); }
Comparison integer_cmp(Integer const& a,Integer const& b) { return cmp(a,b); }
BigUnsigned unsigned_mul(BigUnsigned const& a,BigUnsigned const& b) { return mul(a,b); }
// R3: compare bounded cross-products directly in fixed stack storage.
// The two product spans are independent of both inputs; no temporary
// heap-owned BigUnsigned is needed for magnitudes up to eight limbs.
int compare_cross_products(BigUnsigned const& an,BigUnsigned const& bd,
                           BigUnsigned const& bn,BigUnsigned const& ad) {
    constexpr std::size_t max_operand_limbs=8u;
    if (an.size()<=max_operand_limbs && bd.size()<=max_operand_limbs &&
            bn.size()<=max_operand_limbs && ad.size()<=max_operand_limbs) {
        constexpr std::size_t product_capacity=2u*max_operand_limbs;
        std::array<BigUnsigned::Limb,product_capacity> left{};
        std::array<BigUnsigned::Limb,product_capacity> right{};
        const std::size_t left_size=an.size()+bd.size();
        const std::size_t right_size=bn.size()+ad.size();
        BigUnsigned::mul_into(an.limbs(),bd.limbs(),
            BigUnsigned::MutableLimbView(left.data(),left_size));
        BigUnsigned::mul_into(bn.limbs(),ad.limbs(),
            BigUnsigned::MutableLimbView(right.data(),right_size));
        // Products can have a zero high limb. Compare all potentially
        // significant words and leave unwritten trailing words zero.
        for (std::size_t i=std::max(left_size,right_size);i!=0u;--i) {
            if (left[i-1u]!=right[i-1u]) {
                return left[i-1u]<right[i-1u] ? -1 : 1;
            }
        }
        return 0;
    }
    return compare(unsigned_mul(an,bd),unsigned_mul(bn,ad));
}
Bool unsigned_is_one(BigUnsigned const& n) { return n.size()==1u && n.limb(0u)==1u; }
Integer multiply_magnitude(Integer const& n,BigUnsigned const& m) {
    if (unsigned_is_one(m)) { return n; }
    return Integer::from_magnitude(mul(n.magnitude(),m),n.negative());
}
BigUnsigned exact_quotient(BigUnsigned const& n,BigUnsigned const& d) {
    return unsigned_is_one(d) ? n : quot(n,d);
}
Integer exact_quotient(Integer const& n,BigUnsigned const& d) {
    return Integer::from_magnitude(exact_quotient(n.magnitude(),d),n.negative());
}
}

template<> class ExtensionOperations<Rational> {
    friend class Rational;
    friend class ExtendedOperations<Rational>;

    static Bool is_nan(Rational const& q) { return q._kind==Rational::Kind::NAN_VALUE; }
    static Bool is_inf(Rational const& q) {
        return q._kind==Rational::Kind::POSITIVE_INFINITY || q._kind==Rational::Kind::NEGATIVE_INFINITY;
    }
    static Bool is_finite(Rational const& q) { return q._kind==Rational::Kind::FINITE; }
    static Bool is_zero(Rational const& q) { return is_finite(q) && q._numerator.magnitude().is_zero(); }
    static Sign sgn(Rational const& q) { return integer_sign(q._numerator); }

    static Void set_nan(Rational& q) { q=Rational::nan(); }
    static Void set_inf(Rational& q, Sign s) { q=Rational::inf(s); }
    static Void set_zero(Rational& q) { q=Rational(); }
};

template<> class FiniteOperations<Rational> {
    friend class ExtendedOperations<Rational>;

    static Rational reduced(Integer n,BigUnsigned d) {
        return Rational(std::move(n),std::move(d),Rational::Reduced{});
    }
    static Void set(Rational& r,Rational const& q) { r=q; }

    static Void add_sub(Rational& r,Rational const& a,Rational const& b,bool subtract) {
        if (a._numerator.magnitude().is_zero()) { r=subtract ? reduced(integer_neg(b._numerator),b._denominator) : b; return; }
        if (b._numerator.magnitude().is_zero()) { r=a; return; }
        const BigUnsigned g=gcd(a._denominator,b._denominator);
        const BigUnsigned ad=exact_quotient(a._denominator,g);
        const BigUnsigned bd=exact_quotient(b._denominator,g);
        const Integer an=multiply_magnitude(a._numerator,bd);
        const Integer bn=multiply_magnitude(b._numerator,ad);
        Integer n=subtract ? integer_sub(an,bn) : integer_add(an,bn);
        if (n.magnitude().is_zero()) { r=Rational(); return; }
        // Only a factor of the original denominator gcd can remain.
        const BigUnsigned h=gcd(n.magnitude(),g);
        r=reduced(exact_quotient(n,h),unsigned_mul(ad,exact_quotient(b._denominator,h)));
    }
    static Void add(Rational& r,Rational const& a,Rational const& b) { add_sub(r,a,b,false); }
    static Void sub(Rational& r,Rational const& a,Rational const& b) { add_sub(r,a,b,true); }
    static Void mul(Rational& r,Rational const& a,Rational const& b) {
        if (&a==&b) { sqr(r,a); return; }
        if (a._numerator.magnitude().is_zero() || b._numerator.magnitude().is_zero()) {
            r=Rational(); return;
        }
        const BigUnsigned g=gcd(a._numerator.magnitude(),b._denominator);
        const BigUnsigned h=gcd(b._numerator.magnitude(),a._denominator);
        r=reduced(integer_mul(exact_quotient(a._numerator,g),exact_quotient(b._numerator,h)),
                  unsigned_mul(exact_quotient(a._denominator,h),exact_quotient(b._denominator,g)));
    }
    static Void div(Rational& r,Rational const& a,Rational const& b) {
        if (b._numerator.magnitude().is_zero()) {
            r=Rational::inf(integer_sign(a._numerator)); return;
        }
        if (a._numerator.magnitude().is_zero()) { r=Rational(); return; }
        const BigUnsigned g=gcd(a._numerator.magnitude(),b._numerator.magnitude());
        const BigUnsigned h=gcd(a._denominator,b._denominator);
        Integer n=Integer::from_magnitude(
            unsigned_mul(exact_quotient(a._numerator.magnitude(),g),exact_quotient(b._denominator,h)),
            a._numerator.negative()!=b._numerator.negative());
        r=reduced(std::move(n),unsigned_mul(exact_quotient(a._denominator,h),
                                          exact_quotient(b._numerator.magnitude(),g)));
    }
    static Void neg(Rational& r,Rational const& q) { r=reduced(integer_neg(q._numerator),q._denominator); }
    static Void hlf(Rational& r,Rational const& q) {
        if ((q._numerator.magnitude().limb(0u)&1u)==0u) {
            r=reduced(Integer::from_magnitude(shift_right(q._numerator.magnitude(),1u),q._numerator.negative()),q._denominator);
        } else { r=reduced(q._numerator,shift_left(q._denominator,1u)); }
    }
    static Void sqr(Rational& r,Rational const& q) {
        r=reduced(Integer::from_magnitude(square(q._numerator.magnitude()),false),square(q._denominator));
    }
    static Void rec(Rational& r,Rational const& q) {
        if (q._numerator.magnitude().is_zero()) { r=Rational::inf(); return; }
        r=reduced(Integer::from_magnitude(q._denominator,q._numerator.negative()),q._numerator.magnitude());
    }
    static Void max(Rational& r,Rational const& a,Rational const& b) { r=cmp(a,b)>=Comparison::EQUAL ? a : b; }
    static Void min(Rational& r,Rational const& a,Rational const& b) { r=cmp(a,b)<=Comparison::EQUAL ? a : b; }
    static Void abs(Rational& r,Rational const& q) {
        r=reduced(Integer::from_magnitude(q._numerator.magnitude(),false),q._denominator);
    }
    static Comparison cmp(Rational const& a,Rational const& b) {
        if (a._numerator.negative()!=b._numerator.negative()) {
            return a._numerator.negative() ? Comparison::LESS : Comparison::GREATER;
        }
        if (compare(a._denominator,b._denominator)==0) { return integer_cmp(a._numerator,b._numerator); }
        int order=compare_cross_products(a._numerator.magnitude(),b._denominator,
                                         b._numerator.magnitude(),a._denominator);
        if (a._numerator.negative()) { order=-order; }
        return order<0 ? Comparison::LESS : order>0 ? Comparison::GREATER : Comparison::EQUAL;
    }
    static Boolean eq(Rational const& a,Rational const& b) {
        return a._numerator==b._numerator && compare(a._denominator,b._denominator)==0;
    }
};

Rational rec(Integer const& z) {
    return Rational(1,z);
}

Rational pow(Integer const& z, Int n) {
    const Nat magnitude=n<0 ? Nat(-(n+1))+Nat(1) : Nat(n);
    const Integer value=pow(z,magnitude);
    return n<0 ? Rational(1,value) : Rational(value);
}

Rational::Rational(Integer n,BigUnsigned d,Reduced)
    : _numerator(std::move(n)),_denominator(std::move(d))
{
    if (_numerator.magnitude().is_zero()) { _denominator=BigUnsigned(1u); }
}

Rational::Rational(Integer const& z) : _numerator(z) { }

Rational::Rational(Integer const& n,Integer const& d) {
    if (d.magnitude().is_zero()) { *this=inf(integer_sign(n)); return; }
    if (n.magnitude().is_zero()) { return; }
    const BigUnsigned g=gcd(n.magnitude(),d.magnitude());
    _numerator=Integer::from_magnitude(exact_quotient(n.magnitude(),g),n.negative()!=d.negative());
    _denominator=exact_quotient(d.magnitude(),g);
}

Rational::Rational(Int64 n) : Rational(Integer(n)) { }

Rational::Rational(ExactDouble const& x) {
    const double d=x.get_d();
    if (std::isnan(d)) { *this=nan(); return; }
    if (std::isinf(d)) { *this=inf(d<0.0 ? Sign::NEGATIVE : Sign::POSITIVE); return; }
    if (d==0.0) { return; }
    int exponent=0;
    const double fraction=std::frexp(std::abs(d),&exponent);
    std::uint64_t significand=static_cast<std::uint64_t>(std::ldexp(fraction,53));
    const int trailing=std::countr_zero(significand);
    significand>>=trailing;
    exponent=exponent-53+trailing;
    BigUnsigned magnitude(significand);
    if (exponent>=0) { magnitude=shift_left(magnitude,static_cast<std::size_t>(exponent)); }
    else { _denominator=shift_left(BigUnsigned(1u),static_cast<std::size_t>(-exponent)); }
    _numerator=Integer::from_magnitude(std::move(magnitude),d<0.0);
}

Rational::Rational(FloatDP const& x) : Rational(ExactDouble(x.get_d())) { }
Rational::Rational(FloatMP const& x) : Rational(x.operator Dyadic()) { }

Rational::Rational(String const& str) {
    // Decimal integer or numerator/denominator; permit ASCII whitespace as
    // the old GMP reader did, but reject extra separators and missing digits.
    std::string token;
    for (char c : str) {
        if (c!=' ' && c!='\t' && c!='\n' && c!='\r' && c!='\f' && c!='\v') { token.push_back(c); }
    }
    try {
        const auto slash=token.find('/');
        const Integer n(String(token.substr(0u,slash)));
        const Integer d=slash==std::string::npos ? Integer(1) : Integer(String(token.substr(slash+1u)));
        *this=Rational(n,d);
    } catch (std::invalid_argument const&) {
        ARIADNE_THROW(std::runtime_error,"Rational(string)",
                      "String \""<<str<<"\" does not have a valid rational number format.");
    }
}

Rational::Rational(Rational&& q) : Rational() {
    std::swap(_numerator,q._numerator);
    _denominator.swap(q._denominator);
    std::swap(_kind,q._kind);
}

Rational& Rational::operator=(Rational const& q) {
    if (this!=&q) {
        Rational copy(q);
        *this=std::move(copy);
    }
    return *this;
}

Rational& Rational::operator=(Rational&& q) {
    std::swap(_numerator,q._numerator);
    _denominator.swap(q._denominator);
    std::swap(_kind,q._kind);
    return *this;
}

double Rational::get_d() const {
    if (is_nan(*this)) { return std::numeric_limits<double>::quiet_NaN(); }
    const bool negative=_numerator.negative();
    const double sign=negative ? -1.0 : 1.0;
    if (is_inf(*this)) { return sign*std::numeric_limits<double>::infinity(); }
    const BigUnsigned& n=_numerator.magnitude();
    if (n.is_zero()) { return 0.0; }
    const std::size_t nb=n.bit_length(),db=_denominator.bit_length();
    if (nb>db && nb-db>1024u) { return sign*std::numeric_limits<double>::infinity(); }
    if (db>nb && db-nb>1075u) { return 0.0; }
    int exponent=nb>=db ? static_cast<int>(nb-db) : -static_cast<int>(db-nb);
    const int order=exponent>=0
        ? compare(n,shift_left(_denominator,static_cast<std::size_t>(exponent)))
        : compare(shift_left(n,static_cast<std::size_t>(-exponent)),_denominator);
    if (order<0) { --exponent; }
    if (exponent>1023) { return sign*std::numeric_limits<double>::infinity(); }
    if (exponent< -1074) { return 0.0; }
    // Truncate toward zero, including subnormals, independent of the active
    // floating rounding mode. The <=53-bit integer is exactly representable.
    const int scale=std::max(exponent-52,-1074);
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    // R7: certify the exact normal binary64 significand from the high
    // one or two limbs. n=N*B+nl, d=D*B+dl, 0<=nl,dl<B.
    // If N*2^shift=q*D+r, then r>=q proves the lower bound and
    // r+2^shift<D proves the strict upper bound for the full quotient.
    // q==2^52 also has a valid lower bound from the exact exponent.
    // Ambiguous or overflowing prefixes retain the complete division.
    if (scale<0 && scale>=-63 && exponent>=-1022 &&
            n.size()==_denominator.size() && n.size()>=2u) {
        __extension__ using Wide=unsigned __int128;
        const unsigned shift=static_cast<unsigned>(-scale);
        const Wide power=Wide{1}<<shift;
        auto certify=[&](Wide high_n,Wide high_d,
                         std::uint64_t& significand) -> bool {
            if (high_n>(~Wide{0}>>shift)) { return false; }
            const Wide scaled=high_n<<shift;
            const Wide candidate=scaled/high_d;
            const Wide remainder=scaled%high_d;
            constexpr Wide min_normal_significand=Wide{1}<<52u;
            constexpr Wide max_normal_significand=Wide{1}<<53u;
            if (candidate<min_normal_significand ||
                    candidate>=max_normal_significand) {
                return false;
            }
            const bool lower_proven=remainder>=candidate ||
                candidate==min_normal_significand;
            const bool upper_proven=high_d>power &&
                remainder<high_d-power;
            if (!lower_proven || !upper_proven) { return false; }
            significand=static_cast<std::uint64_t>(candidate);
            return true;
        };
        const std::size_t high=n.size()-1u;
        const Wide high_n=n.limb(high);
        const Wide high_d=_denominator.limb(high);
        std::uint64_t significand=0u;
        if (certify(high_n,high_d,significand)) {
            return std::copysign(
                std::ldexp(static_cast<double>(significand),scale),sign);
        }
        const Wide top_n=(high_n<<64u)|n.limb(high-1u);
        const Wide top_d=(high_d<<64u)|_denominator.limb(high-1u);
        if (certify(top_n,top_d,significand)) {
            return std::copysign(
                std::ldexp(static_cast<double>(significand),scale),sign);
        }
    }
#endif
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    // R5: short finite ratios need at most 128 bits of exact arithmetic.
    // Do not use floating approximation: the conversion must remain
    // identical under every active host rounding mode.
    if (n.size()==1u && _denominator.size()==1u &&
            scale<0 && scale>=-63) {
        __extension__ using Wide=unsigned __int128;
        const unsigned shift=static_cast<unsigned>(-scale);
        const Wide numerator=static_cast<Wide>(n.limb(0u))<<shift;
        const auto significand=static_cast<std::uint64_t>(
            numerator/_denominator.limb(0u));
        return std::copysign(std::ldexp(static_cast<double>(significand),scale),sign);
    }
#endif
    // R5: for similarly sized large operands, compute the quotient from
    // three high limbs and *prove* it equals the full exact quotient.
    // Writing n=N*B+nl, d=D*B+dl (0<=nl,dl<B), let
    // N*2^s=q*D+r. Then r>=q proves the lower bound; r+2^s<D
    // proves the upper bound. For normal values, the exact binary
    // exponent independently guarantees q>=2^52.
    // Ambiguous prefixes always use the original full-width division.
    if (scale<0 && scale>=-63 && exponent>=-1022 &&
            n.size()==_denominator.size() && n.size()>4u) {
        constexpr std::size_t retained_limbs=3u;
        const std::size_t drop=n.size()-retained_limbs;
        const unsigned shift=static_cast<unsigned>(-scale);
        const BigUnsigned high_n(n.limbs().subspan(drop));
        const BigUnsigned high_d(_denominator.limbs().subspan(drop));
        auto [candidate,remainder]=divrem(shift_left(high_n,shift),high_d);
        const bool lower_proven=compare(remainder,candidate)>=0 ||
            (candidate.size()==1u &&
             candidate.limb(0u)==(std::uint64_t{1}<<52u));
        if (candidate.size()==1u && lower_proven &&
                compare(add(remainder,shift_left(BigUnsigned(1u),shift)),
                        high_d)<0) {
            return std::copysign(
                std::ldexp(static_cast<double>(candidate.limb(0u)),scale),sign);
        }
    }
    const BigUnsigned digits=scale<0
        ? quot(shift_left(n,static_cast<std::size_t>(-scale)),_denominator)
        : quot(n,shift_left(_denominator,static_cast<std::size_t>(scale)));
    return std::copysign(std::ldexp(static_cast<double>(digits.limb(0u)),scale),sign);
}

Integer Rational::get_num() const { return _numerator; }
Integer Rational::get_den() const { return Integer::from_magnitude(_denominator,false); }
Integer Rational::numerator() const { return get_num(); }
Natural Rational::denominator() const { return Natural(get_den()); }

Rational operator/(Integer const& z1, Integer const& z2) {
    return Rational(z1,z2);
}


namespace { Rational const& unsign(PositiveRational const& q) { return q; } }

Rational max(Rational const& q1, Rational const& q2) {
    Rational r; ExtendedOperations<Rational>::max(r,q1,q2); return r;
}

PositiveRational max(Rational const& q1, PositiveRational const& q2) {
    return cast_positive(max(q1,unsign(q2)));
}
PositiveRational max(PositiveRational const& q1, Rational const& q2) {
    return cast_positive(max(unsign(q1),q2));
}
PositiveRational max(PositiveRational const& q1, PositiveRational const& q2) {
    return cast_positive(max(unsign(q1),unsign(q2)));
}

Rational min(Rational const& q1, Rational const& q2) {
    Rational r; ExtendedOperations<Rational>::min(r,q1,q2); return r;
}
PositiveRational min(PositiveRational const& q1, PositiveRational const& q2) {
    return cast_positive(min(unsign(q1),unsign(q2)));
}

Positive<Rational> abs(Rational const& q) {
    Positive<Rational> r; ExtendedOperations<Rational>::abs(r,q); return r;
}

Rational mag(Rational const& q) {
    Rational r; ExtendedOperations<Rational>::abs(r,q); return r;
}

Rational mig(Rational const& q) {
    Rational r; ExtendedOperations<Rational>::abs(r,q); return r;
}

Rational nul(Rational const&) {
    return Rational(0);
}

Rational pos(Rational const& q) {
    return q;
}

Rational neg(Rational const& q) {
    Rational r; ExtendedOperations<Rational>::neg(r,q); return r;
}

Positive<Rational> sqr(Rational const& q) {
    Positive<Rational> r; ExtendedOperations<Rational>::sqr(r,q); return r;
}

Rational hlf(Rational const& q) {
    Rational r; ExtendedOperations<Rational>::hlf(r,q); return r;
}

Rational rec(Rational const& q) {
    Rational r; ExtendedOperations<Rational>::rec(r,q); return r;
}

Rational add(Rational const& q1, Rational const& q2) {
    Rational r; ExtendedOperations<Rational>::add(r,q1,q2); return r;
}

Rational sub(Rational const& q1, Rational const& q2) {
    Rational r; ExtendedOperations<Rational>::sub(r,q1,q2); return r;
}

Rational mul(Rational const& q1, Rational const& q2) {
    Rational r; ExtendedOperations<Rational>::mul(r,q1,q2); return r;
}

Rational div(Rational const& q1, Rational const& q2) {
    Rational r; ExtendedOperations<Rational>::div(r,q1,q2); return r;
}

Rational div(Integer const& z1, Integer const& z2) {
    return Rational(z1,z2);
}

Rational fma(Rational const& a,Rational const& b,Rational const& c) {
    return add(mul(a,b),c);
}

Rational pow(Rational const& q, Nat m) {
    if (m==0u) { return Rational(1); }
    if (is_nan(q)) { return Rational::nan(); }
    if (is_inf(q)) {
        return Rational::inf(q._numerator.negative() && (m%2u)!=0u ? Sign::NEGATIVE : Sign::POSITIVE);
    }
    return Rational(Integer::from_magnitude(Ariadne::pow(q._numerator.magnitude(),m),
                                            q._numerator.negative() && (m%2u)!=0u),
                    Ariadne::pow(q._denominator,m),Rational::Reduced{});
}

Rational pow(Rational const& q, Int n) {
    const Nat magnitude=n<0 ? Nat(-(n+1))+Nat(1) : Nat(n);
    return pow(n<0 ? rec(q) : q,magnitude);
}

Boolean eq(Rational const& q1, Rational const& q2) {
    return ExtendedOperations<Rational>::eq(q1,q2);
}

Boolean lt(Rational const& q1, Rational const& q2) {
    return cmp(q1,q2)<Comparison::EQUAL;
}

Rational Rational::inf(Sign sgn) {
    Rational q(sgn==Sign::ZERO ? 0 : sgn==Sign::POSITIVE ? +1 : -1);
    q._denominator=BigUnsigned();
    q._kind=sgn==Sign::ZERO ? Kind::NAN_VALUE : sgn==Sign::POSITIVE ? Kind::POSITIVE_INFINITY : Kind::NEGATIVE_INFINITY;
    return q;
}

Rational Rational::inf() {
    return Rational::inf(Sign::POSITIVE);
}

Rational Rational::nan() {
    return Rational::inf(Sign::ZERO);
}

Bool is_nan(Rational const& q) {
    return ExtendedOperations<Rational>::is_nan(q);
}

Bool is_inf(Rational const& q) {
    return ExtendedOperations<Rational>::is_inf(q);
}

Bool is_finite(Rational const& q) {
    return ExtendedOperations<Rational>::is_finite(q);
}

Bool is_zero(Rational const& q) {
    return ExtendedOperations<Rational>::is_zero(q);
}

Sign sgn(Rational const& q) {
    return ExtendedOperations<Rational>::sgn(q);
}

Integer round(Rational const& q) {
    ARIADNE_PRECONDITION(is_finite(q));
    auto division=divrem(q._numerator.magnitude(),q._denominator);
    if (compare(shift_left(division.second,1u),q._denominator)>=0) {
        division.first=add(division.first,BigUnsigned(1u));
    }
    return Integer::from_magnitude(std::move(division.first),q._numerator.negative());
}

Integer ceil(Rational const& q) {
    ARIADNE_PRECONDITION(is_finite(q));
    // R4: for |n|<d, ceil is 1 iff the numerator is positive.
    // Canonical zero remains 0; other fractions retain exact division.
    if (compare(q._numerator.magnitude(),q._denominator)<0) {
        return !q._numerator.negative() && !q._numerator.magnitude().is_zero()
            ? Integer(1) : Integer(0);
    }
    auto division=divrem(q._numerator.magnitude(),q._denominator);
    if (!q._numerator.negative() && !division.second.is_zero()) {
        division.first=add(division.first,BigUnsigned(1u));
    }
    return Integer::from_magnitude(std::move(division.first),q._numerator.negative());
}

Integer floor(Rational const& q) {
    ARIADNE_PRECONDITION(is_finite(q));
    // R4: proper fractions have a known quotient and sign-correct floor.
    // Avoid quotient/remainder materialization, including for large limbs.
    if (compare(q._numerator.magnitude(),q._denominator)<0) {
        return q._numerator.negative() ? Integer(-1) : Integer(0);
    }
    auto division=divrem(q._numerator.magnitude(),q._denominator);
    if (q._numerator.negative() && !division.second.is_zero()) {
        division.first=add(division.first,BigUnsigned(1u));
    }
    return Integer::from_magnitude(std::move(division.first),q._numerator.negative());
}


Comparison cmp(Rational const& q1, Rational const& q2) {
    return ExtendedOperations<Rational>::cmp(q1,q2);
}

Comparison cmp(Rational const& q1, ExactDouble const& x2) {
    double d2=x2.get_d();
    ARIADNE_PRECONDITION(not std::isnan(d2));
    return cmp(q1,Rational(x2));
}

Comparison cmp(ExactDouble const& x1, Rational const& q2) {
    return Comparison(-(int)cmp(q2,x1));
}

Rational operator""_q(unsigned long long int n) {
    return Rational(operator""_z(n));
}

Rational operator""_q(long double x) {
    static const uint32_t max_cf_coef = std::numeric_limits<uint32_t>::max();
    static const std::size_t N=11;

    // See if the long double value is an exact single-precision value
    volatile float sx=static_cast<float>(x);
    if(static_cast<long double>(sx)==x) {
        return Rational(ExactDouble(sx));
    }
    // Compute the continued fraction expansion of x, storing coefficients in cf
    // Stop if a coefficient is larger than max_cf_coef, since then the result
    //   is accurate to within 2^{-32}
    long double t=x;
    Int s=1;
    if(x<0) { s=-1; t=-x; }
    long double cf[N+1];
    cf[N]=0;
    std::size_t i;
    for(i=0; i!=N; ++i) {
        cf[i]=std::floor(t);
        t=1/(t-cf[i]);
        if(t>max_cf_coef) { break; }
    }
    if(i==N)
        ARIADNE_THROW(InvalidRationalLiteralException,"Rational operator""_q(long double)",
                      "x="<<x<<" is not a sufficiently close approximation to a simple rational number.");
    // Compute the result from the continued fraction coefficients
    Rational q = Rational(ExactDouble(cf[i]));
    while(i!=0) {
        --i;
        q=Rational(ExactDouble(cf[i]))+rec(q);
    }
    if(s==-1) { q=-q; }
    double xd=static_cast<double>(x);
    Rational xq=Rational(ExactDouble(xd));
    //volatile double qd=q.get_d();
    double ae=std::abs((q-xq).get_d());
    double re=ae/std::max(1.0,std::abs(xd));
    if(re>std::numeric_limits<double>::epsilon())
        ARIADNE_THROW(InvalidRationalLiteralException,"Rational operator""_q(long double)",
                      "Rational approximation q="<<q<<" to x="<<x<<"="<<xd<<" has error "<<ae<<" and relative error "<<re<<" while is larger than machine epsilon");
    return q;
}

InputStream& operator>>(InputStream& is, Rational& q) {
    String token;
    is >> token;
    if(is) { q=Rational(token); }
    return is;
}

OutputStream& operator<<(OutputStream& os, Rational const& q1) {
    Integer num=q1.get_num();
    Integer den=q1.get_den();
    os << num;
    if(den!=Integer(1)) { os << "/" << den; }
    return os;
}

struct RoundExact { };
inline RoundExact opposite(RoundExact) { return RoundExact(); }
template<class Y> inline decltype(auto) sqr(RoundExact, Y const& y) { return sqr(y); }
template<class Y> inline decltype(auto) rec(RoundExact, Y const& y) { return rec(y); }
template<class Y1, class Y2> inline decltype(auto) mul(RoundExact, Y1 const& y1, Y2 const& y2) { return y1*y2; }
template<class Y1, class Y2> inline decltype(auto) div(RoundExact, Y1 const& y1, Y2 const& y2) { return y1/y2; }

template<class RNDUP, class Y> inline auto _mul(RNDUP rnd_up, Bounds<Y> const& y1, Bounds<Y> const& y2) -> Bounds<Y>
{
    auto rnd_down=opposite(rnd_up);
    const Y& y1l=y1.lower_raw(); const Y& y1u=y1.upper_raw();
    const Y& y2l=y2.lower_raw(); const Y& y2u=y2.upper_raw();
    if(y1l>=0) {
        if(y2l>=0) {
            return Bounds<Y>(mul(rnd_down,y1l,y2l),mul(rnd_up,y1u,y2u));
        } else if(y2u<=0) {
            return Bounds<Y>(mul(rnd_down,y1u,y2l),mul(rnd_up,y1l,y2u));
        } else {
            return Bounds<Y>(mul(rnd_down,y1u,y2l),mul(rnd_up,y1u,y2u));
        }
    }
    else if(y1u<=0) {
        if(y2l>=0) {
            return Bounds<Y>(mul(rnd_down,y1l,y2u),mul(rnd_up,y1u,y2l));
        } else if(y2u<=0) {
            return Bounds<Y>(mul(rnd_down,y1u,y2u),mul(rnd_up,y1l,y2l));
        } else {
            return Bounds<Y>(mul(rnd_down,y1l,y2u),mul(rnd_up,y1l,y2l));
        }
    } else {
        if(y2l>=0) {
            return Bounds<Y>(mul(rnd_down,y1l,y2u),mul(rnd_up,y1u,y2u));
        } else if(y2u<=0) {
            return Bounds<Y>(mul(rnd_down,y1u,y2l),mul(rnd_up,y1l,y2l));
        } else {
            return Bounds<Y>(min(mul(rnd_down,y1u,y2l),mul(rnd_down,y1l,y2u)),max(mul(rnd_up,y1l,y2l),mul(rnd_up,y1u,y2u)));
        }
    }
}

template<class RNDUP, class Y> inline auto _div(RNDUP rnd_up, Bounds<Y> const& y1, Bounds<Y> const& y2) -> Bounds<Y>
{
    auto rnd_down=opposite(rnd_up);
    const Y& y1l=y1.lower_raw(); const Y& y1u=y1.upper_raw();
    const Y& y2l=y2.lower_raw(); const Y& y2u=y2.upper_raw();
    // IMPORTANT: Need to be careful when one of the bounds is 0, since if y2l=-0.0 and y1u>0, then y2l>=0 but y1u/y2l=-inf
    if(y2l>0) {
        if(y1l>=0) {
            return Bounds<Y>(div(rnd_down,y1l,y2u),div(rnd_up,y1u,y2l));
        } else if(y1u<=0) {
            return Bounds<Y>(div(rnd_down,y1l,y2l),div(rnd_up,y1u,y2u));
        } else {
            return Bounds<Y>(div(rnd_down,y1l,y2l),div(rnd_up,y1u,y2l));
        }
    }
    else if(y2u<0) {
        if(y1l>=0) {
            return Bounds<Y>(div(rnd_down,y1u,y2u),div(rnd_up,y1l,y2l));
        } else if(y1u<=0) {
            return Bounds<Y>(div(rnd_down,y1u,y2l),div(rnd_up,y1l,y2u));
        } else {
            return Bounds<Y>(div(rnd_down,y1u,y2u),div(rnd_up,y1l,y2u));
        }
    }
    else {
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds y1, FloatBounds y2)","y1="<<y1<<", y2="<<y2);
        if constexpr(HasPrecisionType<Y>) {
            auto pr=max(y1.precision(),y2.precision());
            return Bounds<Y>(-Y::inf(pr),+Y::inf(pr));
        } else {
            return Bounds<Y>(-Y::inf(),+Y::inf());
        }
    }
}

template<class RNDUP, class Y> Bounds<Y> _sqr(RNDUP rnd_up, Bounds<Y> const& y) {
    auto rnd_down=opposite(rnd_up);
    const Y& yl=y.lower_raw(); const Y& yu=y.upper_raw();
    if(yl>0) {
        return Bounds<Y>(sqr(rnd_down,yl),sqr(rnd_up,yu));
    } else if(yu<0) {
        return Bounds<Y>(sqr(rnd_down,yu),sqr(rnd_up,yl));
    } else {
        return Bounds<Y>(nul(yl),max(sqr(rnd_up,yl),sqr(rnd_up,yu)));
    }
}

template<class RNDUP, class Y> Bounds<Y> _rec(RNDUP rnd_up, Bounds<Y> const& y) {
    auto rnd_down=opposite(rnd_up);
    const Y& yl=y.lower_raw(); const Y& yu=y.upper_raw();
    if(yl>0 || yu<0) {
        return Bounds<Y>(rec(rnd_down,yu),rec(rnd_up,yl));
    } else {
        if constexpr(HasPrecisionType<Y>) {
            auto pr=y.precision();
            return Bounds<Y>(-Y::inf(pr),+Y::inf(pr));
        } else {
            return Bounds<Y>(-Y::inf(),+Y::inf());
        }
    }
}

RationalBounds operator*(RationalBounds const& q1, RationalBounds const& q2) {
    return _mul(RoundExact(),q1,q2); }
RationalBounds operator/(RationalBounds const& q1, RationalBounds const& q2) {
    return _div(RoundExact(),q1,q2); }


RationalBounds nul(RationalBounds const& q) {
    return RationalBounds(nul(q._l),nul(q._u)); }
RationalBounds pos(RationalBounds const& q) {
    return RationalBounds(pos(q._l),pos(q._u)); }
RationalBounds neg(RationalBounds const& q) {
    return RationalBounds(neg(q._u),neg(q._l)); }
RationalBounds hlf(RationalBounds const& q) {
    return RationalBounds(hlf(q._l),hlf(q._u)); }
RationalBounds sqr(RationalBounds const& q) {
    return _sqr(RoundExact(),q); }
RationalBounds rec(RationalBounds const& q) {
    return _rec(RoundExact(),q); }

RationalBounds add(RationalBounds const& q1, RationalBounds const& q2) {
    return RationalBounds(q1._l+q2._l,q1._u+q2._u); }
RationalBounds sub(RationalBounds const& q1, RationalBounds const& q2) {
    return RationalBounds(q1._l-q2._u,q1._u-q2._l); }
RationalBounds mul(RationalBounds const& q1, RationalBounds const& q2) {
    return _mul(RoundExact(),q1,q2); }
RationalBounds div(RationalBounds const& q1, RationalBounds const& q2) {
    return _div(RoundExact(),q1,q2); }
RationalBounds pow(RationalBounds const& q, Nat m) {
    RationalBounds r = (m%2==0) ? abs(q) : q;
    Int n=static_cast<Int>(m);
    return RationalBounds(pow(r._l,n),pow(r._u,n)); }
RationalBounds pow(RationalBounds const& q, Int n) {
    if(n<0) { return pow(rec(q),Nat(-n)); } else return pow(q,Nat(n));}

static_assert(not Convertible<Int,PositiveRational>);
RationalBounds abs(RationalBounds const& q) {
    if(q._l>=0) { return q; }
    if(q._u<=0) { return RationalBounds(-q._u,-q._l); }
    return RationalBounds(0,max(-q._l,q._u)); }
RationalBounds max(RationalBounds const& q1, RationalBounds const& q2) {
    return RationalBounds(max(q1._l,q2._l),max(q1._u,q2._u)); }
RationalBounds min(RationalBounds const& q1, RationalBounds const& q2) {
    return RationalBounds(min(q1._l,q2._l),min(q1._u,q2._u)); }

ValidatedKleenean operator==(RationalBounds const& q1, RationalBounds const& q2) {
    if (q1._l >= q2._u && q1._u <= q2._l) { return true; } else if (q1._l > q2._u || q1._u < q2._l) { return false; } else { return indeterminate; } }
ValidatedKleenean operator!=(RationalBounds const& q1, RationalBounds const& q2) {
    if (q1._l > q2._u || q1._u < q2._l) { return true; } else if (q1._l >= q2._u && q1._u <= q2._l) { return false; } else { return indeterminate; } }
ValidatedKleenean operator<=(RationalBounds const& q1, RationalBounds const& q2) {
    if (q1._u <= q2._l) { return true; } else if (q1._l > q2._u) { return false; } else { return indeterminate; } }
ValidatedKleenean operator>=(RationalBounds const& q1, RationalBounds const& q2) {
    if (q1._l >= q2._u) { return true; } else if (q1._u < q2._l) { return false; } else { return indeterminate; } }
ValidatedKleenean operator< (RationalBounds const& q1, RationalBounds const& q2) {
    if (q1._u < q2._l) { return true; } else if (q1._l >= q2._u) { return false; } else { return indeterminate; } }
ValidatedKleenean operator> (RationalBounds const& q1, RationalBounds const& q2) {
    if (q1._l > q2._u) { return true; } else if (q1._u <= q2._l) { return false; } else { return indeterminate; } }

template<> String class_name<Rational>() { return "Rational"; }
template<> String class_name<RationalBounds>() { return "RationalBounds"; }


// TODO: Move to dyadic.cpp
DyadicBounds operator*(DyadicBounds const& w1, DyadicBounds const& w2) {
    return _mul(RoundExact(),w1,w2); }

DyadicBounds mul(DyadicBounds const& w1, DyadicBounds const& w2) {
    return _mul(RoundExact(),w1,w2); }

    DyadicBounds pow(DyadicBounds const& w, Int n) {
    ARIADNE_PRECONDITION(n>=0);
    return pow(w,static_cast<Nat>(n)); }

static_assert(Field<Rational>);
static_assert(Lattice<Rational>);

} // namespace Ariadne
