/***************************************************************************
 *            numeric/float64.cpp
 *
 *  Copyright  2008-22  Pieter Collins
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

#include "utility/standard.hpp"

#include <iostream>
#include <iomanip>
#include <cassert>
#include <bit>
#include <cstdint>
#include <cmath>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <utility>



#include "macro/general.hpp"
#include "numeric/builtin.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/rounding.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/rounded_float.hpp"
#include "numeric/float_approximation.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"

#include "numeric/concepts.hpp"

namespace Ariadne {

/************  Publicly-accessible rounding-mode changing *******************/

static const double _pi_down=3.1415926535897931;
static const double _pi_near=3.1415926535897931;
static const double _pi_up  =3.1415926535897936;

int abslog10floor(double x);


static_assert(not GenericNumber<FloatDP>);

namespace {
// Only run this path under ambient RN.  Conservative operand gates
// exclude any nonzero residual smaller than 2^-904, including the
// otherwise dangerous subnormal/underflow corner cases.
constexpr double certified_floor=0x1p-400;
double adjusted(double n,double residual,FloatDP::RoundingModeType mode) {
    if ((mode==FloatDP::ROUND_UPWARD && residual>0.0) ||
        (mode==FloatDP::ROUND_DOWNWARD && residual<0.0)) {
        // D10: all callers certify finite normal binary64 results.
        // Step the IEEE-754 encoding by one, reversing for negatives.
        // This also preserves normal/subnormal and overflow boundaries.
        auto bits=std::bit_cast<std::uint64_t>(n);
        if ((n>0.0)==(mode==FloatDP::ROUND_UPWARD)) { ++bits; }
        else { --bits; }
        return std::bit_cast<double>(bits);
    }
    return n;
}
double two_sum_error(double a,double b,double s) {
    const double t=s-a;
    return (a-(s-t))+(b-t);
}
bool exact_add(FloatDP::RoundingModeType mode,double a,double b,double& result) {
    const double s=a+b;
    if (!std::isnormal(s)||!std::isfinite(a)||!std::isfinite(b))return false;
    result=adjusted(s,two_sum_error(a,b,s),mode);
    return true;
}
bool exact_mul(FloatDP::RoundingModeType mode,double a,double b,double& result) {
    const double p=a*b;
    if (!std::isnormal(p)||std::fabs(a)<certified_floor
        ||std::fabs(b)<certified_floor)return false;
    result=adjusted(p,std::fma(a,b,-p),mode);
    return true;
}
bool exact_div(FloatDP::RoundingModeType mode,double a,double b,double& result) {
    const double q=a/b;
    if (!std::isnormal(q)||std::fabs(a)<certified_floor
        ||std::fabs(b)<certified_floor||std::fabs(q)<certified_floor)return false;
    const double e=std::fma(-q,b,a); // sign(a/b-q) = sign(e/b)
    result=adjusted(q,std::signbit(b)?-e:e,mode);
    return true;
}
bool exact_sqrt(FloatDP::RoundingModeType mode,double a,double& result) {
    const double r=std::sqrt(a);
    if (!std::isnormal(r)||r<certified_floor)return false;
    result=adjusted(r,std::fma(-r,r,a),mode);
    return true;
}
bool exact_fma(FloatDP::RoundingModeType mode,double a,double b,
               double c,double& result) {
    const double r=std::fma(a,b,c);
    if (!std::isnormal(r)||std::fabs(a)<certified_floor
        ||std::fabs(b)<certified_floor||std::fabs(c)<certified_floor
        ||std::fabs(r)<certified_floor)return false;
    // D8: Sterbenz guarantees c-r is exact if c and r have the same
    // sign and are within a factor of two.  One further correctly
    // rounded FMA then gives the sign of the exact a*b+c-r residual.
    // The existing 2^-400 gate keeps every nonzero residual normal
    // or representable subnormal, so zero cannot hide underflow.
    if ((r>0.0 && c>0.0)||(r<0.0 && c<0.0)) {
        const double ar=std::fabs(r);
        const double ac=std::fabs(c);
        if (ar>=ac*0.5 && ac>=ar*0.5) {
            result=adjusted(r,std::fma(a,b,c-r),mode);
            return true;
        }
    }
    // Retain D7's independently certified expansion when the
    // Sterbenz fast path cannot be used.
    const double p=a*b;
    if (!std::isnormal(p)||std::fabs(p)<certified_floor)return false;
    const double s=p+c;
    if (s!=r)return false;
    const double e=std::fma(a,b,-p)+two_sum_error(p,c,s);
    result=adjusted(r,e,mode);
    return true;
}
bool exact_pow(FloatDP::RoundingModeType mode,double x,Nat n,double& result) {
    if (!(x>0.0)||!std::isfinite(x)||x<certified_floor||n>16u)return false;
    double p=x,r=1.0;
    bool first_factor=true;
    while(n!=0u) {
        if ((n&1u)!=0u) {
            // Multiplying by the initial exact 1.0 cannot round.
            if (first_factor) { r=p; first_factor=false; }
            else if (!exact_mul(mode,r,p,r)) { return false; }
        }
        n/=2u;
        // The last square is discarded by the original algorithm.
        if (n!=0u && !exact_mul(mode,p,p,p)) { return false; }
    }

    result=r;
    return true;
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
// D12: keep the complete D11 FMA3 kernels, but expose the instruction
// through a local GNU/Clang x86-64 inline-assembly primitive, rather
// than an out-of-line function with a different target attribute.
// The caller selects this path only after checking CPU FMA3 support.
// On Windows and other architectures the original portable kernel stays.
__attribute__((always_inline)) inline double fma3_inline(double a,
                                                           double b,double c) {
    __asm__ volatile("vfmadd231sd %2, %1, %0"
        : "+x"(c) : "x"(a), "x"(b));
    return c;
}
__attribute__((always_inline)) inline
bool exact_mul_fma3(FloatDP::RoundingModeType mode,
                    double a,double b,double& result) {
    const double p=a*b;
    if (!std::isnormal(p)||std::fabs(a)<certified_floor
        ||std::fabs(b)<certified_floor)return false;
    result=adjusted(p,fma3_inline(a,b,-p),mode);
    return true;
}
__attribute__((always_inline)) inline
bool exact_div_fma3(FloatDP::RoundingModeType mode,
                    double a,double b,double& result) {
    const double q=a/b;
    if (!std::isnormal(q)||std::fabs(a)<certified_floor
        ||std::fabs(b)<certified_floor||std::fabs(q)<certified_floor)return false;
    const double e=fma3_inline(-q,b,a);
    result=adjusted(q,std::signbit(b)?-e:e,mode);
    return true;
}
__attribute__((always_inline)) inline
bool exact_fma_fma3(FloatDP::RoundingModeType mode,
                    double a,double b,double c,double& result) {
    const double r=fma3_inline(a,b,c);
    if (!std::isnormal(r)||std::fabs(a)<certified_floor
        ||std::fabs(b)<certified_floor||std::fabs(c)<certified_floor
        ||std::fabs(r)<certified_floor)return false;
    if ((r>0.0 && c>0.0)||(r<0.0 && c<0.0)) {
        const double ar=std::fabs(r);
        const double ac=std::fabs(c);
        if (ar>=ac*0.5 && ac>=ar*0.5) {
            result=adjusted(r,fma3_inline(a,b,c-r),mode);
            return true;
        }
    }
    const double p=a*b;
    if (!std::isnormal(p)||std::fabs(p)<certified_floor)return false;
    const double s=p+c;
    if (s!=r)return false;
    const double e=fma3_inline(a,b,-p)+two_sum_error(p,c,s);
    result=adjusted(r,e,mode);
    return true;
}
__attribute__((always_inline)) inline
bool exact_pow_fma3(FloatDP::RoundingModeType mode,
                    double x,Nat n,double& result) {
    if (!(x>0.0)||!std::isfinite(x)||x<certified_floor||n>16u)return false;
#if defined(__GNUC__) && !defined(__clang__)
    // D13: GCC does not fully specialize the general repeated-squaring
    // loop for degree seven. Use the identical four directed products,
    // in their original order, with the same certification/fallback gates.
    // Clang already optimizes this loop effectively; leave it unchanged.
    if (n==7u) {
        double square,cube,fourth,seventh;
        if (!exact_mul_fma3(mode,x,x,square)) return false;
        if (!exact_mul_fma3(mode,x,square,cube)) return false;
        if (!exact_mul_fma3(mode,square,square,fourth)) return false;
        if (!exact_mul_fma3(mode,cube,fourth,seventh)) return false;
        result=seventh;
        return true;
    }
#endif
    double p=x,r=1.0;
    bool first_factor=true;
    while(n!=0u) {
        if ((n&1u)!=0u) {
            if (first_factor) { r=p; first_factor=false; }
            else if (!exact_mul_fma3(mode,r,p,r)) { return false; }
        }
        n/=2u;
        if (n!=0u && !exact_mul_fma3(mode,p,p,p)) { return false; }
    }
    result=r;
    return true;
}
bool cpu_has_fma3() {
#if defined(__FMA__)
    return true;
#else
    static const bool available=__builtin_cpu_supports("fma");
    return available;
#endif
}
#endif

bool exact_mul_dispatch(FloatDP::RoundingModeType mode,
                        double a,double b,double& result) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
    if (cpu_has_fma3()) return exact_mul_fma3(mode,a,b,result);
#endif
    return exact_mul(mode,a,b,result);
}
bool exact_div_dispatch(FloatDP::RoundingModeType mode,
                        double a,double b,double& result) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
    if (cpu_has_fma3()) return exact_div_fma3(mode,a,b,result);
#endif
    return exact_div(mode,a,b,result);
}
bool exact_fma_dispatch(FloatDP::RoundingModeType mode,
                        double a,double b,double c,double& result) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
    if (cpu_has_fma3()) return exact_fma_fma3(mode,a,b,c,result);
#endif
    return exact_fma(mode,a,b,c,result);
}
bool exact_pow_dispatch(FloatDP::RoundingModeType mode,
                        double a,Nat n,double& result) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
    if (cpu_has_fma3()) return exact_pow_fma3(mode,a,n,result);
#endif
    return exact_pow(mode,a,n,result);
}

} // namespace

FloatDP FloatDP::_certified_add(RoundingModeType mode,FloatDP a,FloatDP b) {
    double r;
    if(get_rounding_mode()==ROUND_TO_NEAREST && exact_add(mode,a.dbl,b.dbl,r))return FloatDP(r);
    return _apply_rnd(&add_rnd,mode,a,b);
}
FloatDP FloatDP::_certified_sub(RoundingModeType mode,FloatDP a,FloatDP b) {
    double r;
    if(get_rounding_mode()==ROUND_TO_NEAREST && exact_add(mode,a.dbl,-b.dbl,r))return FloatDP(r);
    return _apply_rnd(&sub_rnd,mode,a,b);
}
FloatDP FloatDP::_certified_mul(RoundingModeType mode,FloatDP a,FloatDP b) {
    double r;
    if(get_rounding_mode()==ROUND_TO_NEAREST && exact_mul_dispatch(mode,a.dbl,b.dbl,r))return FloatDP(r);
    return _apply_rnd(&mul_rnd,mode,a,b);
}
FloatDP FloatDP::_certified_div(RoundingModeType mode,FloatDP a,FloatDP b) {
    double r;
    if(get_rounding_mode()==ROUND_TO_NEAREST && exact_div_dispatch(mode,a.dbl,b.dbl,r))return FloatDP(r);
    return _apply_rnd(&div_rnd,mode,a,b);
}
FloatDP FloatDP::_certified_fma(RoundingModeType mode,FloatDP a,FloatDP b,FloatDP c) {
    double r;
    if(get_rounding_mode()==ROUND_TO_NEAREST && exact_fma_dispatch(mode,a.dbl,b.dbl,c.dbl,r))return FloatDP(r);
    return _apply_rnd(&fma_rnd,mode,a,b,c);
}
FloatDP FloatDP::_certified_sqrt(RoundingModeType mode,FloatDP a) {
    double r;
    if(get_rounding_mode()==ROUND_TO_NEAREST && exact_sqrt(mode,a.dbl,r))return FloatDP(r);
    return _apply_rnd(&sqrt_rnd,mode,a);
}
FloatDP FloatDP::_certified_pow(RoundingModeType mode,FloatDP a,Nat n) {
    double r;
    if(get_rounding_mode()==ROUND_TO_NEAREST && exact_pow_dispatch(mode,a.dbl,n,r))return FloatDP(r);
    return _apply_rnd((double(*)(double,unsigned int))&pow_rnd,mode,a,n);
}

Rounded<FloatDP>::Rounded(Approximation<FloatDP> const& x)
    : Rounded(x.raw())
{
}

FloatDP::Float(ExactDouble const& d, PrecisionType)
    : FloatDP(d.get_d())
{
}

FloatDP::Float(TwoExp const& t, PrecisionType)
    : FloatDP(std::ldexp(1.0,t.exponent()))
{
    ARIADNE_PRECONDITION(Dyadic(*this)==Dyadic(t));
}

FloatDP::Float(Dyadic const& w, PrecisionType)
    : FloatDP(w.get_d())
{
    ARIADNE_PRECONDITION(Dyadic(*this)==w || is_nan(w));
}

FloatDP::Float(String const& str, PrecisionType pr)
    : FloatDP(Dyadic(str),pr)
{
}

FloatDP::Float(Integer const& x, RoundingModeType rnd, PrecisionType pr)
    : FloatDP(Dyadic(x),rnd,pr)
{
}

namespace {

double rounded_dyadic_value(Dyadic const& w, BuiltinRoundingModeType rnd) {
    double d=w.get_d();
    if (is_finite(w)) {
        auto old_rnd=FloatDP::get_rounding_mode();
        if(rnd==FloatDP::ROUND_UPWARD) {
            FloatDP::set_rounding_upward();
            while (Dyadic(d)<w) { d=d+std::numeric_limits<double>::min(); }
            FloatDP::set_rounding_mode(old_rnd);
        }
        if(rnd==FloatDP::ROUND_DOWNWARD) {
            FloatDP::set_rounding_downward();
            while (Dyadic(d)>w) { d=d-std::numeric_limits<double>::min(); }
            FloatDP::set_rounding_mode(old_rnd);
        }
    }
    return d;
}

inline Rational cast_rational(double d) { return Rational(ExactDouble(d)); }

// D9: Convert small finite fractions without constructing a second Rational,
// comparing large cross-products or modifying the processor rounding mode.
// A quotient of two nonzero uint64 values is always a normal binary64 value.
bool rounded_short_rational_value(Rational const& q,
                                  BuiltinRoundingModeType rnd,double& result) {
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    if (!is_finite(q)) { return false; }
    const Integer& numerator=q.numerator_reference();
    const BigUnsigned& magnitude=numerator.magnitude();
    const BigUnsigned& denominator=q.denominator_reference();
    if (magnitude.size()!=1u || denominator.size()!=1u) { return false; }
    const std::uint64_t n=magnitude.limb(0u);
    const std::uint64_t d=denominator.limb(0u);
    __extension__ using Wide=unsigned __int128;
    const int bits_n=std::bit_width(n);
    const int bits_d=std::bit_width(d);
    const int candidate=bits_n-bits_d;
    const bool below=candidate>=0
        ? Wide(n)<(Wide(d)<<static_cast<unsigned>(candidate))
        : (Wide(n)<<static_cast<unsigned>(-candidate))<Wide(d);
    // For nonzero 64-bit operands the exact exponent is in [-64,63].
    const int exponent=candidate-static_cast<int>(below);
    const int scale=exponent-52;
    Wide dividend=n;
    Wide divisor=d;
    if (scale<=0) { dividend<<=static_cast<unsigned>(-scale); }
    else { divisor<<=static_cast<unsigned>(scale); }
    // Both scaled magnitudes fit in 128 bits. The quotient is exactly
    // in [2^52,2^53), hence is representable as a double without rounding.
    std::uint64_t significand=0u;
    Wide remainder=0u;
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    if (scale<=0 && (dividend>>64u)<Wide(d)) {
        // One native 128/64 DIVQ, also returning the exact remainder.
        std::uint64_t low=static_cast<std::uint64_t>(dividend);
        std::uint64_t high=static_cast<std::uint64_t>(dividend>>64u);
        __asm__ volatile("divq %[divisor]"
            : "+a"(low), "+d"(high)
            : [divisor] "r"(d)
            : "cc");
        significand=low;
        remainder=high;
    } else
#endif
    {
        significand=static_cast<std::uint64_t>(dividend/divisor);
        remainder=dividend%divisor;
    }
    const bool negative=numerator.negative();
    const bool away=(rnd==FloatDP::ROUND_UPWARD && !negative) ||
                    (rnd==FloatDP::ROUND_DOWNWARD && negative);
    // Preserve the existing nearest/toward-zero truncation contract;
    // the remainder exactly distinguishes an inexact directed result.
    if (away && remainder!=0u) { ++significand; }
    result=std::copysign(
        std::ldexp(static_cast<double>(significand),scale),
        negative ? -1.0 : 1.0);
    return true;
#else
    (void)q; (void)rnd; (void)result;
    return false;
#endif
}

double rounded_rational_value(Rational const& q, BuiltinRoundingModeType rnd) {
    double fast;
    if (rounded_short_rational_value(q,rnd,fast)) { return fast; }
    double d=q.get_d();
    if (is_finite(q)) {
        auto old_rnd=FloatDP::get_rounding_mode();
        if(rnd==FloatDP::ROUND_UPWARD) {
            FloatDP::set_rounding_upward();
            while (cast_rational(d)<q) { d=d+std::numeric_limits<double>::min(); }
            FloatDP::set_rounding_mode(old_rnd);
        }
        if(rnd==FloatDP::ROUND_DOWNWARD) {
            FloatDP::set_rounding_downward();
            while (cast_rational(d)>q) { d=d-std::numeric_limits<double>::min(); }
            FloatDP::set_rounding_mode(old_rnd);
        }
    }
    return d;
}

} // namespace

FloatDP::Float(Dyadic const& w, RoundingModeType rnd, PrecisionType)
    : FloatDP(rounded_dyadic_value(w,rnd))
{
}

FloatDP::Float(Decimal const& dec, RoundingModeType rnd, PrecisionType pr)
    : FloatDP(Rational(dec),rnd,pr)
{
}

FloatDP::Float(Rational const& q, RoundingModeType rnd, PrecisionType)
    : FloatDP(rounded_rational_value(q,rnd))
{
}

// Convert the exact FloatMP value directly at double precision with the
// requested rounding. Preserve the non-finite cases and signed zero, which
// are not faithfully represented by ordinary finite rational arithmetic.
FloatDP::Float(FloatMP const& value, RoundingModeType rnd, PrecisionType pr)
    : FloatDP(value.is_nan() ? FloatDP::nan(pr)
              : value.is_inf() ? FloatDP::inf(value.negative() ? Sign::NEGATIVE : Sign::POSITIVE,pr)
              : value.is_zero() ? FloatDP(ExactDouble(value.negative() ? -0.0 : 0.0),pr)
              : FloatDP(static_cast<Rational>(value),rnd,pr))
{
    // Reject unsupported rounding modes for every kind of FloatMP input,
    // including exact values, infinities, NaNs, and signed zero.
    ARIADNE_PRECONDITION(rnd==ROUND_TO_NEAREST || rnd==ROUND_DOWNWARD
                      || rnd==ROUND_UPWARD || rnd==ROUND_TOWARD_ZERO);
}

namespace {

enum class FloatDPTranscendental {
    sine, cosine, tangent, arctangent, arcsine, arccosine,
    exponential, logarithm, hyperbolic_tangent
};

// FloatDP's certified near-half kernel is deliberately independent of
// FloatMP.  Coefficients are outward dyadic enclosures at scale 2^-85,
// derived from the repository's independently generated 128-bit exact
// transcendental oracles.  The few derivative bounds that are not direct
// function values are obtained by exact integer identities:
//   tan' = 1+tan^2, atan' = 4/5,
//   asin' = -acos' = 2/sqrt(3), tanh' = 1-tanh^2.
// Runtime evaluation therefore needs only the shared unsigned-limb
// multiplier plus binary64 bit manipulation; no FloatMP object, precision
// tag, adaptive series, or multiprecision transcendental cache is touched.
struct FloatDPFixed85 {
    std::uint64_t low;
    std::uint64_t high;
};

struct FloatDPCenter {
    FloatDPFixed85 lower;
    FloatDPFixed85 upper;
    FloatDPFixed85 slope_lower;
    FloatDPFixed85 slope_upper;
    unsigned error;
    bool slope_negative;
    bool result_negative;
};

constexpr unsigned FLOATDP_FIXED_BITS=85u;

constexpr FloatDPFixed85 fixed85(std::uint64_t low,std::uint64_t high) {
    return {low,high};
}

template<FloatDPTranscendental Function>
constexpr FloatDPCenter floatdp_center(bool negative) {
    if constexpr(Function==FloatDPTranscendental::sine) {
        return {
            fixed85(0x3a2582f7f43b25e1ULL,0x00000000000f5774ULL),
            fixed85(0x3a2582f7f43b25e2ULL,0x00000000000f5774ULL),
            fixed85(0x065b7d4f9db7bbb3ULL,0x00000000001c1528ULL),
            fixed85(0x065b7d4f9db7bbb4ULL,0x00000000001c1528ULL),
            1u,false,negative
        };
    }
    else if constexpr(Function==FloatDPTranscendental::cosine) {
        return {
            fixed85(0x065b7d4f9db7bbb3ULL,0x00000000001c1528ULL),
            fixed85(0x065b7d4f9db7bbb4ULL,0x00000000001c1528ULL),
            fixed85(0x3a2582f7f43b25e1ULL,0x00000000000f5774ULL),
            fixed85(0x3a2582f7f43b25e2ULL,0x00000000000f5774ULL),
            1u,true,false
        };
    }
    else if constexpr(Function==FloatDPTranscendental::tangent) {
        return {
            fixed85(0x5bf3474a43179648ULL,0x0000000000117b4fULL),
            fixed85(0x5bf3474a43179649ULL,0x0000000000117b4fULL),
            fixed85(0x7c8a28fc590a79f8ULL,0x0000000000298cdfULL),
            fixed85(0x7c8a28fc590a79faULL,0x0000000000298cdfULL),
            6u,false,negative
        };
    }
    else if constexpr(Function==FloatDPTranscendental::arctangent) {
        return {
            fixed85(0x82b0dda7b456fe44ULL,0x00000000000ed633ULL),
            fixed85(0x82b0dda7b456fe45ULL,0x00000000000ed633ULL),
            fixed85(0x9999999999999999ULL,0x0000000000199999ULL),
            fixed85(0x999999999999999aULL,0x0000000000199999ULL),
            1u,false,negative
        };
    }
    else if constexpr(Function==FloatDPTranscendental::arcsine) {
        return {
            fixed85(0x382d73658465bb32ULL,0x000000000010c152ULL),
            fixed85(0x382d73658465bb33ULL,0x000000000010c152ULL),
            fixed85(0x8b2066389a431f03ULL,0x000000000024f34eULL),
            fixed85(0x8b2066389a431f04ULL,0x000000000024f34eULL),
            2u,false,negative
        };
    }
    else if constexpr(Function==FloatDPTranscendental::arccosine) {
        if(negative) {
            return {
                fixed85(0xe0b5cd961196eccbULL,0x0000000000430548ULL),
                fixed85(0xe0b5cd961196ecccULL,0x0000000000430548ULL),
                fixed85(0x8b2066389a431f03ULL,0x000000000024f34eULL),
                fixed85(0x8b2066389a431f04ULL,0x000000000024f34eULL),
                2u,false,false
            };
        }
        return {
            fixed85(0x705ae6cb08cb7665ULL,0x00000000002182a4ULL),
            fixed85(0x705ae6cb08cb7666ULL,0x00000000002182a4ULL),
            fixed85(0x8b2066389a431f03ULL,0x000000000024f34eULL),
            fixed85(0x8b2066389a431f04ULL,0x000000000024f34eULL),
            2u,true,false
        };
    }
    else if constexpr(Function==FloatDPTranscendental::hyperbolic_tangent) {
        return {
            fixed85(0xebab4579b2894771ULL,0x00000000000ec9a9ULL),
            fixed85(0xebab4579b2894772ULL,0x00000000000ec9a9ULL),
            fixed85(0x6fa343940dfe7c96ULL,0x0000000000192a94ULL),
            fixed85(0x6fa343940dfe7c98ULL,0x0000000000192a94ULL),
            1u,false,negative
        };
    }
    else if constexpr(Function==FloatDPTranscendental::exponential) {
        if(negative) {
            return {
                fixed85(0xfc6f9609fe7aceb4ULL,0x00000000001368b2ULL),
                fixed85(0xfc6f9609fe7aceb5ULL,0x00000000001368b2ULL),
                fixed85(0xfc6f9609fe7aceb4ULL,0x00000000001368b2ULL),
                fixed85(0xfc6f9609fe7aceb5ULL,0x00000000001368b2ULL),
                2u,true,false
            };
        }
        return {
            fixed85(0x1c3c0d3792e5bfdfULL,0x000000000034c253ULL),
            fixed85(0x1c3c0d3792e5bfe0ULL,0x000000000034c253ULL),
            fixed85(0x1c3c0d3792e5bfdfULL,0x000000000034c253ULL),
            fixed85(0x1c3c0d3792e5bfe0ULL,0x000000000034c253ULL),
            2u,false,false
        };
    }
    else if constexpr(Function==FloatDPTranscendental::logarithm) {
        return {
            fixed85(0xfefa39ef35793c76ULL,0x0000000000162e42ULL),
            fixed85(0xfefa39ef35793c77ULL,0x0000000000162e42ULL),
            fixed85(0x0000000000000000ULL,0x0000000000400000ULL),
            fixed85(0x0000000000000000ULL,0x0000000000400000ULL),
            8u,true,true
        };
    }
    else {
        static_assert(Function==FloatDPTranscendental::logarithm,
                      "Unknown FloatDP transcendental");
    }
}

// Exhaust the complete near-half input lattice at compile time. This checks
// the stored certificates, rather than relying on a sample of test inputs:
// 2048 binary64 values below 1/2 and 1025 at or above it. The model uses
// 32-bit partial products, independently of the runtime limb multiplier.
consteval std::uint64_t floatdp_certificate_correction(
    FloatDPFixed85 slope,std::uint64_t offset,unsigned bits,bool round_upward) {
    const std::uint64_t bottom=(slope.low&0xffffffffu)*offset;
    const std::uint64_t top=(slope.low>>32u)*offset+(bottom>>32u);
    const std::uint64_t low=(top<<32u)|(bottom&0xffffffffu);
    const std::uint64_t high=slope.high*offset+(top>>32u);
    const std::uint64_t remainder=low&((std::uint64_t{1}<<bits)-1u);
    return (low>>bits)+(high<<(64u-bits))+(round_upward && remainder!=0u);
}

consteval FloatDPFixed85 floatdp_certificate_endpoint(
    FloatDPFixed85 center,std::uint64_t correction,bool subtract) {
    const std::uint64_t low=subtract ? center.low-correction : center.low+correction;
    const std::uint64_t high=subtract ? center.high-(center.low<correction)
                                      : center.high+(low<center.low);
    return {low,high};
}

consteval bool floatdp_verify_certificate_slice(
    FloatDPCenter center,unsigned bits,unsigned first,unsigned last) {
    // Leave at least a whole limb of margin at the endpoint width limits.
    // All displacements are at most 2^44+8, and all products fit 98 bits.
    if(center.lower.high<=(std::uint64_t{1}<<19u) ||
       center.upper.high>=(std::uint64_t{1}<<23u)-1u ||
       center.slope_upper.high>=(std::uint64_t{1}<<23u) ||
       center.slope_lower.high>=(std::uint64_t{1}<<23u) || center.error>8u)
        return false;
    for(unsigned offset=first;offset<=last;++offset) {
        const auto cl=floatdp_certificate_correction(center.slope_lower,offset,bits,false);
        const auto cu=floatdp_certificate_correction(center.slope_upper,offset,bits,true);
        const bool subtract=(bits==54u)!=center.slope_negative;
        // Applying the error separately keeps this model entirely unsigned.
        const auto lower=floatdp_certificate_endpoint(
            floatdp_certificate_endpoint(center.lower,subtract ? cu : cl,subtract),
            center.error,true);
        const auto upper=floatdp_certificate_endpoint(
            floatdp_certificate_endpoint(center.upper,subtract ? cl : cu,subtract),
            center.error,false);
        const unsigned width=128u-static_cast<unsigned>(std::countl_zero(lower.high));
        const unsigned upper_width=128u-static_cast<unsigned>(std::countl_zero(upper.high));
        if(width<84u || width>87u || upper_width!=width) return false;
        const unsigned shift=width-53u;
        const auto retained=[shift](FloatDPFixed85 value) {
            return (value.low>>shift)|(value.high<<(64u-shift));
        };
        const auto l=retained(lower),u=retained(upper);
        if(l!=u || l==(std::uint64_t{1}<<53u)-1u) return false;
        const std::uint64_t mask=(std::uint64_t{1}<<shift)-1u;
        const std::uint64_t half=std::uint64_t{1}<<(shift-1u);
        const auto lr=lower.low&mask,ur=upper.low&mask;
        // Both endpoints are inexact, never ties, with the same nearest
        // result. Equal truncations also prove agreement in directed modes.
        if(lr==0u || ur==0u || lr==half || ur==half || (lr>half)!=(ur>half))
            return false;
    }
    return true;
}

// Separate small constant expressions respect default compiler evaluation
// limits; they require no nonstandard increase of the constexpr step budget.
template<FloatDPCenter Center,unsigned Bits,unsigned Slice>
struct FloatDPCertificateSlice {
    static constexpr unsigned first=Bits==54u ? 1u+Slice*256u : Slice*256u;
    static constexpr unsigned last=Bits==54u ? (Slice+1u)*256u
                                             : (Slice==4u ? 1024u : (Slice+1u)*256u-1u);
    static constexpr bool valid=floatdp_verify_certificate_slice(Center,Bits,first,last);
    static_assert(valid,"FloatDP near-half certificate no longer permits specialized rounding");
};

template<FloatDPTranscendental Function,bool Negative,std::size_t... I>
consteval bool floatdp_verify_certificate(std::index_sequence<I...>) {
    constexpr auto center=floatdp_center<Function>(Negative);
    return (FloatDPCertificateSlice<center,54u,I>::valid && ...) &&
           FloatDPCertificateSlice<center,53u,0u>::valid &&
           FloatDPCertificateSlice<center,53u,1u>::valid &&
           FloatDPCertificateSlice<center,53u,2u>::valid &&
           FloatDPCertificateSlice<center,53u,3u>::valid &&
           FloatDPCertificateSlice<center,53u,4u>::valid;
}

template<FloatDPTranscendental Function>
consteval bool floatdp_verify_certificates() {
    return floatdp_verify_certificate<Function,false>(std::make_index_sequence<8u>{}) &&
           floatdp_verify_certificate<Function,true>(std::make_index_sequence<8u>{});
}

static_assert(floatdp_verify_certificates<FloatDPTranscendental::sine>());
static_assert(floatdp_verify_certificates<FloatDPTranscendental::cosine>());
static_assert(floatdp_verify_certificates<FloatDPTranscendental::tangent>());
static_assert(floatdp_verify_certificates<FloatDPTranscendental::arctangent>());
static_assert(floatdp_verify_certificates<FloatDPTranscendental::arcsine>());
static_assert(floatdp_verify_certificates<FloatDPTranscendental::arccosine>());
static_assert(floatdp_verify_certificates<FloatDPTranscendental::exponential>());
static_assert(floatdp_verify_certificates<FloatDPTranscendental::logarithm>());
static_assert(floatdp_verify_certificates<FloatDPTranscendental::hyperbolic_tangent>());

std::uint64_t floatdp_correction(FloatDPFixed85 slope,std::uint64_t offset,
                                 unsigned denominator_bits,bool round_upward) {
    // F1: a 2x1-limb product with portable native uint64_t arithmetic.
    // The exhaustive certificate bounds slope<2^87 and offset<=2^11;
    // 32-bit partial products and carry additions cannot overflow.
    // No general multiprecision multiplication or 128-bit extension is needed.
    const std::uint64_t bottom=(slope.low&0xffffffffu)*offset;
    const std::uint64_t top=(slope.low>>32u)*offset+(bottom>>32u);
    const std::uint64_t low=(top<<32u)|(bottom&0xffffffffu);
    const std::uint64_t high=slope.high*offset+(top>>32u);

    // In the certified window denominator_bits is exactly 53 or 54.
    const std::uint64_t remainder_mask=(std::uint64_t{1}<<denominator_bits)-1u;
    const std::uint64_t retained=(low>>denominator_bits)|
        (high<<(64u-denominator_bits));
    return retained+static_cast<std::uint64_t>(
        round_upward && (low&remainder_mask)!=0u);
}

constexpr FloatDPFixed85 floatdp_endpoint(FloatDPFixed85 center,
                                         std::int64_t displacement) {
    // F2: the certified endpoint is exactly two uint64_t words.
    // Apply the one-word signed correction with a native carry or borrow.
    // The existing consteval certificate guarantees no 128-bit overflow
    // or underflow for any argument accepted by this short kernel.
    const std::uint64_t absolute=displacement<0
        ? std::uint64_t{0}-static_cast<std::uint64_t>(displacement)
        : static_cast<std::uint64_t>(displacement);
    if(displacement<0) {
        const std::uint64_t low=center.low-absolute;
        return {low,center.high-static_cast<std::uint64_t>(
            center.low<absolute)};
    }
    const std::uint64_t low=center.low+absolute;
    return {low,center.high+static_cast<std::uint64_t>(
        low<center.low)};
}

// Explicit exact word-boundary checks complement the source-keyed
// binary64 transcendental oracles and exhaustive center certificates.
static_assert(floatdp_endpoint({0u,3u},-1).low==
              std::numeric_limits<std::uint64_t>::max());
static_assert(floatdp_endpoint({0u,3u},-1).high==2u);
static_assert(floatdp_endpoint({std::numeric_limits<std::uint64_t>::max(),3u},
                                1).low==0u);
static_assert(floatdp_endpoint({std::numeric_limits<std::uint64_t>::max(),3u},
                                1).high==4u);
static_assert(floatdp_endpoint({1u,3u},-1).low==0u);

// F3: build an IEEE-754 binary64 normal number from the rigorously
// certified 53-bit normalized significand and unbiased exponent.
// Fixed85 certificates confine the output exponent to [-2,+1], exclude
// rounding carries, and guarantee every endpoint is a normal result.
static_assert(std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::radix==2 &&
              std::numeric_limits<double>::digits==53 &&
              sizeof(double)==sizeof(std::uint64_t));
constexpr std::uint64_t floatdp_normal_bits(bool negative,int highest,
                                           std::uint64_t significand) {
    return (static_cast<std::uint64_t>(negative)<<63u) |
           (static_cast<std::uint64_t>(highest+1023)<<52u) |
           (significand&((std::uint64_t{1}<<52u)-1u));
}
static_assert(floatdp_normal_bits(false,0,std::uint64_t{1}<<52u)==
              0x3ff0000000000000ULL);
static_assert(floatdp_normal_bits(true,-1,std::uint64_t{1}<<52u)==
              0xbfe0000000000000ULL);
static_assert(floatdp_normal_bits(false,-2,std::uint64_t{3}<<51u)==
              0x3fd8000000000000ULL);
static_assert(floatdp_normal_bits(false,1,std::uint64_t{5}<<50u)==
              0x4004000000000000ULL);

// The exhaustive compile-time certificate proves that both endpoints are
// inexact, never halfway, never carry into a new exponent, and round alike
// in every supported mode. This private rounder uses exactly that domain.
bool floatdp_round_fixed(FloatDPFixed85 magnitude,
                         bool negative,FloatDP::RoundingModeType mode,double& output) {
    const std::uint64_t low=magnitude.low;
    const std::uint64_t high=magnitude.high;
    const unsigned width=128u-static_cast<unsigned>(std::countl_zero(high));
    const unsigned shift=width-53u;
    std::uint64_t retained=(low>>shift)|(high<<(64u-shift));
    bool increment=false;
    if(mode==FloatDP::ROUND_TO_NEAREST) {
        const std::uint64_t halfway=std::uint64_t{1}<<(shift-1u);
        increment=(low&halfway)!=0u;
    } else if(mode==FloatDP::ROUND_UPWARD) {
        increment=!negative;
    } else if(mode==FloatDP::ROUND_DOWNWARD) {
        increment=negative;
    } else if(mode!=FloatDP::ROUND_TOWARD_ZERO) {
        return false;
    }
    if(increment) ++retained;
    const int highest=static_cast<int>(width)-1-static_cast<int>(FLOATDP_FIXED_BITS);
    // The former ldexp step only scaled an exact 53-bit integer by a
    // power of two. Emit the same bits directly without floating-point
    // arithmetic, regardless of the ambient hardware rounding mode.
    output=std::bit_cast<double>(floatdp_normal_bits(negative,highest,retained));
    return true;
}

template<FloatDPTranscendental Function>
bool floatdp_certified_transcendental(
                                      FloatDP::RoundingModeType mode,
                                      double input,double& output) {
    // F6: every near-half operation, including exp_nearest, is covered
    // by the existing exhaustive consteval nearest-even certificate.
    // No special nearest-mode exclusion remains; inputs outside the
    // certified window still use their original fallback.

    const std::uint64_t encoded=std::bit_cast<std::uint64_t>(input);
    const bool negative=(encoded>>63u)!=0u;
    const std::uint64_t exponent=(encoded>>52u)&0x7ffu;
    if(exponent!=1021u && exponent!=1022u) return false;
    if constexpr(Function==FloatDPTranscendental::logarithm) {
        if(negative) return false;
    }

    const unsigned denominator_bits=static_cast<unsigned>(1023u+52u-exponent);
    const std::uint64_t half=std::uint64_t{1}<<(denominator_bits-1u);
    const std::uint64_t magnitude=
        (encoded&((std::uint64_t{1}<<52u)-1u))|(std::uint64_t{1}<<52u);
    const bool offset_negative=magnitude<half;
    const std::uint64_t offset=offset_negative ? half-magnitude : magnitude-half;

    // At scale B=85 the quadratic remainder certificates used by all
    // centers are bounded by their small integer error constants whenever
    // |delta| <= 2^-43.
    constexpr unsigned cutoff=(FLOATDP_FIXED_BITS+1u)/2u;
    // denominator_bits is 53 or 54, both greater than cutoff=43.
    if(offset>(std::uint64_t{1}<<(denominator_bits-cutoff))) {
        return false;
    }

    const FloatDPCenter center=floatdp_center<Function>(negative);
    const bool subtract=offset_negative!=center.slope_negative;
    const std::uint64_t correction=floatdp_correction(
        subtract ? center.slope_upper : center.slope_lower,
        offset,denominator_bits,subtract);
    const std::int64_t lower_offset=
        (subtract ? -static_cast<std::int64_t>(correction)
                  :  static_cast<std::int64_t>(correction))
        -static_cast<std::int64_t>(center.error);
    const auto lower_magnitude=floatdp_endpoint(center.lower,lower_offset);
    return floatdp_round_fixed(lower_magnitude,center.result_negative,mode,output);
}

} // namespace

FloatDP exp(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental<FloatDPTranscendental::exponential>(mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&exp_rnd,mode,x);
}

FloatDP log(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental<FloatDPTranscendental::logarithm>(mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&log_rnd,mode,x);
}

FloatDP sin(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental<FloatDPTranscendental::sine>(mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&sin_rnd,mode,x);
}

FloatDP cos(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental<FloatDPTranscendental::cosine>(mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&cos_rnd,mode,x);
}

FloatDP tan(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental<FloatDPTranscendental::tangent>(mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&tan_rnd,mode,x);
}

FloatDP atan(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental<FloatDPTranscendental::arctangent>(mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&atan_rnd,mode,x);
}

FloatDP asin(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental<FloatDPTranscendental::arcsine>(mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&asin_rnd,mode,x);
}

FloatDP acos(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental<FloatDPTranscendental::arccosine>(mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&acos_rnd,mode,x);
}

FloatDP tanh(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental<FloatDPTranscendental::hyperbolic_tangent>(
            mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&tanh_rnd,mode,x);
}
FloatDP::Float(FloatDP const& x, RoundingModeType, PrecisionType)
    : FloatDP(x)
{
}

FloatDP& FloatDP::operator=(ExactDouble const& x) {
    this->dbl=x.get_d();
    return *this;
}

FloatDP& FloatDP::operator=(Dyadic const& w) {
    if (is_finite(w)) {
        this->dbl=w.get_d();
        ARIADNE_PRECONDITION(*this==w);
    } else if (is_nan(w)) {
        *this=FloatDP::nan(this->precision());
    } else {
        *this=FloatDP::inf(sgn(w),this->precision());
    }
    return *this;
}

FloatDP::operator Dyadic () const {
    return Dyadic(this->dbl);
}

FloatDP::operator Rational () const {
    return Rational(ExactDouble(this->dbl));
}

FloatDP pow_rnd(FloatDP x, Int n)
{
    return FloatDP(pow_rnd(x.dbl,n));
}

FloatDP sqrt_rnd(FloatDP x)
{
    return FloatDP(sqrt_rnd(x.dbl));
}

FloatDP exp_rnd(FloatDP x)
{
    return FloatDP(exp_rnd(x.dbl));
}

FloatDP log_rnd(FloatDP x)
{
    return FloatDP(log_rnd(x.dbl));
}

FloatDP sin_rnd(FloatDP x)
{
    return FloatDP(sin_rnd(x.dbl));
}

FloatDP cos_rnd(FloatDP x)
{
    return FloatDP(cos_rnd(x.dbl));
}

FloatDP tan_rnd(FloatDP x)
{
    return FloatDP(tan_rnd(x.dbl));
}

FloatDP tanh_rnd(FloatDP x)
{
    return FloatDP(tanh_rnd(x.dbl));
}

FloatDP atan_rnd(FloatDP x)
{
    return FloatDP(atan_rnd(x.dbl));
}

FloatDP FloatDP::pi(BuiltinRoundingModeType rnd, DoublePrecision) {
    ARIADNE_PRECONDITION(rnd==ROUND_TO_NEAREST || rnd==ROUND_UPWARD || rnd==ROUND_DOWNWARD);
    if(rnd==ROUND_UPWARD) { return FloatDP(_pi_up); }
    if(rnd==ROUND_DOWNWARD) { return FloatDP(_pi_down); }
    return FloatDP(_pi_near);
}

FloatDP::RoundingModeType FloatDP::get_rounding_mode() { return Ariadne::get_builtin_rounding_mode(); }
Void FloatDP::set_rounding_mode(RoundingModeType rnd) {
    ARIADNE_PRECONDITION(rnd==ROUND_TO_NEAREST || rnd==ROUND_DOWNWARD ||
                         rnd==ROUND_UPWARD || rnd==ROUND_TOWARD_ZERO);
    Ariadne::set_builtin_rounding_mode(rnd);
}
Void FloatDP::set_rounding_downward() { Ariadne::set_builtin_rounding_downward(); }
Void FloatDP::set_rounding_upward() { Ariadne::set_builtin_rounding_upward(); }
Void FloatDP::set_rounding_to_nearest() { Ariadne::set_builtin_rounding_to_nearest(); }
Void FloatDP::set_rounding_toward_zero() { Ariadne::set_builtin_rounding_toward_zero(); }

FloatDP::PrecisionType FloatDP::get_default_precision() { return FloatDP::PrecisionType(); }
FloatDP::PrecisionType FloatDP::precision() const { return FloatDP::PrecisionType(); }
Void FloatDP::set_precision(FloatDP::PrecisionType) { }

FloatDP FloatDP::min(PrecisionType) { return FloatDP(std::numeric_limits<double>::min()); }
FloatDP FloatDP::max(PrecisionType) { return FloatDP(std::numeric_limits<double>::max()); }
FloatDP FloatDP::eps(PrecisionType) { return FloatDP(std::numeric_limits<double>::epsilon()); }

FloatDP FloatDP::inf(Sign sgn, PrecisionType) {
    switch (sgn) {
    case Sign::POSITIVE: return FloatDP(std::numeric_limits<double>::infinity());
    case Sign::NEGATIVE: return FloatDP(-std::numeric_limits<double>::infinity());
    default: return FloatDP(std::numeric_limits<double>::quiet_NaN());
    }
}
FloatDP FloatDP::inf(PrecisionType pr) { return FloatDP::inf(Sign::POSITIVE,pr); }
FloatDP FloatDP::nan(PrecisionType pr) { return FloatDP::inf(Sign::ZERO,pr); }

Integer cast_integer(FloatDP x) {
    Dyadic w(x);
    Integer z=round(w);
    ARIADNE_PRECONDITION_MSG(z==w,"Cannot cast non-integral value "<<z<<" to an Integer");
    return z;
}

template<class R, class A> R integer_cast(A const&);
template<> Nat integer_cast<Nat,FloatDP>(FloatDP const& x) { return static_cast<Nat>(x.dbl); }
template<> Int integer_cast<Int,FloatDP>(FloatDP const& x) { return static_cast<Int>(x.dbl); }


Comparison cmp(FloatDP x1, Rational const& q2) {
    if(std::isfinite(x1.get_d())) { return cmp(Rational(x1),q2); }
    else { return x1.get_d()>0.0 ? Comparison::GREATER : Comparison::LESS; }
}
Comparison cmp(Rational const& q1, FloatDP x2) {
    if(std::isfinite(x2.get_d())) { return cmp(q1,Rational(x2)); }
    else { return x2.get_d()>0.0 ? Comparison::LESS : Comparison::GREATER; }
}

/*
OutputStream& write(OutputStream& os, FloatDP const& x, Nat bits, BuiltinRoundingModeType rnd) {
    Nat dgts = std::ceil(bits*std::log(2))/std::log(10);
    FloatDP::RoundingModeType old_rnd=FloatDP::get_rounding_mode();
    FloatDP::set_rounding_mode(rnd);
    os << std::setprecision(dgts) << x.get_d();
    FloatDP::set_rounding_mode(old_rnd);
    return os;
}
*/

// Preserve exact binary64 inputs, including signed zero and IEEE
// special values. Decimal rounding uses the native multiple-precision
// formatter without foreign floating-point temporaries.
namespace {
FloatMP::RoundingModeType native_decimal_mode(BuiltinRoundingModeType mode) {
    ARIADNE_PRECONDITION(mode==ROUND_TO_NEAREST ||
                         mode==ROUND_DOWNWARD || mode==ROUND_UPWARD);
    return mode==ROUND_TO_NEAREST ? FloatMP::ROUND_TO_NEAREST :
           mode==ROUND_DOWNWARD ? FloatMP::ROUND_DOWNWARD :
                                  FloatMP::ROUND_UPWARD;
}
FloatMP native_binary64(FloatDP const& x) {
    constexpr unsigned long bits=53u;
    const double value=x.get_d();
    if(std::isnan(value)) return FloatMP::nan(bits);
    if(std::isinf(value)) return FloatMP::infinity(bits,std::signbit(value));
    if(value==0.0) return FloatMP::zero(bits,std::signbit(value));
    return FloatMP(ExactDouble(value),MultiplePrecision(std::size_t{bits}));
}
} // namespace

String FloatDP::literal() const {
    StringStream ss;
    write(ss,*this,DecimalPrecision(53),near);
    return ss.str();
}

String FloatDP::literal(RoundingModeType rnd) const {
    StringStream ss;
    write(ss,*this,DecimalPrecision(18),rnd);
    return ss.str();
}

OutputStream& write(OutputStream& os, FloatDP const& x, DecimalPlaces plcs, BuiltinRoundingModeType rnd) {
    return os << print(native_binary64(x),plcs,native_decimal_mode(rnd));
}
OutputStream& write(OutputStream& os, FloatDP const& x, DecimalPrecision figs, BuiltinRoundingModeType rnd) {
    return os << print(native_binary64(x),figs,native_decimal_mode(rnd));
}

OutputStream& repr(OutputStream& os, FloatDP const& x) {
    return write(os,x,DecimalPrecision(17),near);
}
OutputStream& repr(OutputStream& os, FloatDP const& x, BuiltinRoundingModeType rnd) {
    static const DecimalPrecision dgts(18);
    return write(os,x,dgts,rnd);
}

OutputStream& operator<<(OutputStream& os, FloatDP const& x) {
    return repr(os,x);
}

InputStream& operator>>(InputStream& is, FloatDP& x) {
    is >> std::ws;
    if (!is) { return is; }

    std::string str;
    int c=is.peek();
    if (c=='+' || c=='-') { str.push_back(static_cast<char>(is.get())); c=is.peek(); }

    Bool have_digits=false;
    while (c>='0' && c<='9') {
        have_digits=true;
        str.push_back(static_cast<char>(is.get()));
        c=is.peek();
    }
    if (c=='.') {
        str.push_back(static_cast<char>(is.get()));
        c=is.peek();
        while (c>='0' && c<='9') {
            have_digits=true;
            str.push_back(static_cast<char>(is.get()));
            c=is.peek();
        }
    }

    if (!have_digits) { is.setstate(std::ios::failbit); return is; }

    if (c=='e' || c=='E') {
        str.push_back(static_cast<char>(is.get()));
        c=is.peek();
        if (c=='+' || c=='-') { str.push_back(static_cast<char>(is.get())); c=is.peek(); }
        Bool have_exponent_digits=false;
        while (c>='0' && c<='9') {
            have_exponent_digits=true;
            str.push_back(static_cast<char>(is.get()));
            c=is.peek();
        }
        if (!have_exponent_digits) { is.setstate(std::ios::failbit); return is; }
    }

    std::istringstream parser(str);
    parser.imbue(std::locale::classic());

    double r=0.0;

    FloatDP::RoundingModeType old_rnd=FloatDP::get_rounding_mode();
    FloatDP::set_rounding_to_nearest();
    Bool parse_failed=not static_cast<Bool>(parser >> r);
    FloatDP::set_rounding_mode(old_rnd);

    if (parse_failed) {
        is.setstate(std::ios::failbit);
        return is;
    }
    x.dbl=r;
    return is;
}

Nat FloatDP::output_places=16;
Void FloatDP::set_output_places(Nat pl) { output_places=pl; }

Float32::operator FloatDP() const { return FloatDP(cast_exact((double)this->flt),dp); }

template<> String class_name<double>() { return "double"; }

template<> String class_name<ApproximateDouble>() { return "ApproximateDouble"; }
template<> String class_name<ExactDouble>() { return "ExactDouble"; }


template<> String class_name<DoublePrecision>() { return "DoublePrecision"; }
template<> String class_name<FloatDP>() { return "FloatDP"; }

template<> String class_name<Rounded<FloatDP>>() { return "Rounded<FloatDP>"; }

Positive<FloatDP> abs(FloatDP x) { return Positive<FloatDP>(FloatDP(std::fabs(x.dbl))); }
Positive<UpperBound<FloatDP>> mag(FloatDP x) { return Positive<UpperBound<FloatDP>>(abs(x)); }
Positive<LowerBound<FloatDP>> mig(FloatDP x) { return Positive<LowerBound<FloatDP>>(abs(x)); }

static_assert(RoundedTranscendentalField<FloatDP>);
static_assert(OrderedLattice<FloatDP>);

} // namespace Ariadne
