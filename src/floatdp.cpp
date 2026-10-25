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
#include <array>
#include <bit>
#include <cstdint>
#include <cmath>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>



#include "macro/general.hpp"
#include "numeric/builtin.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/rounding.hpp"
#include "numeric/big_unsigned.hpp"
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

double rounded_rational_value(Rational const& q, BuiltinRoundingModeType rnd) {
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
// FloatMP.  Coefficients are outward dyadic enclosures at scale 2^-96,
// derived from the repository's independently generated 128-bit exact
// transcendental oracles.  The few derivative bounds that are not direct
// function values are obtained by exact integer identities:
//   tan' = 1+tan^2, atan' = 4/5,
//   asin' = -acos' = 2/sqrt(3), tanh' = 1-tanh^2.
// Runtime evaluation therefore needs only the shared unsigned-limb
// multiplier plus binary64 bit manipulation; no FloatMP object, precision
// tag, adaptive series, or multiprecision transcendental cache is touched.
struct FloatDPFixed96 {
    std::uint64_t low;
    std::uint64_t high;
};

struct FloatDPCenter {
    FloatDPFixed96 lower;
    FloatDPFixed96 upper;
    FloatDPFixed96 slope_lower;
    FloatDPFixed96 slope_upper;
    unsigned error;
    bool slope_negative;
    bool result_negative;
};

constexpr unsigned FLOATDP_FIXED_BITS=96u;

constexpr FloatDPFixed96 fixed96(std::uint64_t low,std::uint64_t high) {
    return {low,high};
}

FloatDPCenter floatdp_center(FloatDPTranscendental function,bool negative) {
    switch(function) {
        case FloatDPTranscendental::sine:
            return {
                fixed96(0x2c17bfa1d92f0d93ULL,0x000000007abba1d1ULL),
                fixed96(0x2c17bfa1d92f0d94ULL,0x000000007abba1d1ULL),
                fixed96(0xdbea7cedbddd9da2ULL,0x00000000e0a94032ULL),
                fixed96(0xdbea7cedbddd9da3ULL,0x00000000e0a94032ULL),
                1u,false,negative
            };
        case FloatDPTranscendental::cosine:
            return {
                fixed96(0xdbea7cedbddd9da2ULL,0x00000000e0a94032ULL),
                fixed96(0xdbea7cedbddd9da3ULL,0x00000000e0a94032ULL),
                fixed96(0x2c17bfa1d92f0d93ULL,0x000000007abba1d1ULL),
                fixed96(0x2c17bfa1d92f0d94ULL,0x000000007abba1d1ULL),
                1u,true,false
            };
        case FloatDPTranscendental::tangent:
            return {
                fixed96(0x9a3a5218bcb2403cULL,0x000000008bda7adfULL),
                fixed96(0x9a3a5218bcb2403dULL,0x000000008bda7adfULL),
                fixed96(0x5147e2c853cfc458ULL,0x000000014c66fbe4ULL),
                fixed96(0x5147e2c853cfc45aULL,0x000000014c66fbe4ULL),
                6u,false,negative
            };
        case FloatDPTranscendental::arctangent:
            return {
                fixed96(0x86ed3da2b7f222f6ULL,0x0000000076b19c15ULL),
                fixed96(0x86ed3da2b7f222f7ULL,0x0000000076b19c15ULL),
                fixed96(0xccccccccccccccccULL,0x00000000ccccccccULL),
                fixed96(0xcccccccccccccccdULL,0x00000000ccccccccULL),
                1u,false,negative
            };
        case FloatDPTranscendental::arcsine:
            return {
                fixed96(0x6b9b2c232dd99707ULL,0x00000000860a91c1ULL),
                fixed96(0x6b9b2c232dd99708ULL,0x00000000860a91c1ULL),
                fixed96(0x0331c4d218f81e4aULL,0x00000001279a7459ULL),
                fixed96(0x0331c4d218f81e4bULL,0x00000001279a7459ULL),
                2u,false,negative
            };
        case FloatDPTranscendental::arccosine:
            if(negative) {
                return {
                    fixed96(0xae6cb08cb7665c1eULL,0x00000002182a4705ULL),
                    fixed96(0xae6cb08cb7665c1fULL,0x00000002182a4705ULL),
                    fixed96(0x0331c4d218f81e4aULL,0x00000001279a7459ULL),
                    fixed96(0x0331c4d218f81e4bULL,0x00000001279a7459ULL),
                    2u,false,false
                };
            }
            return {
                fixed96(0xd73658465bb32e0fULL,0x000000010c152382ULL),
                fixed96(0xd73658465bb32e10ULL,0x000000010c152382ULL),
                fixed96(0x0331c4d218f81e4aULL,0x00000001279a7459ULL),
                fixed96(0x0331c4d218f81e4bULL,0x00000001279a7459ULL),
                2u,true,false
            };
        case FloatDPTranscendental::hyperbolic_tangent:
            return {
                fixed96(0x5a2bcd944a3b8871ULL,0x00000000764d4f5dULL),
                fixed96(0x5a2bcd944a3b8872ULL,0x00000000764d4f5dULL),
                fixed96(0x1a1ca06ff3e4becdULL,0x00000000c954a37dULL),
                fixed96(0x1a1ca06ff3e4becfULL,0x00000000c954a37dULL),
                1u,false,negative
            };
        case FloatDPTranscendental::exponential:
            if(negative) {
                return {
                    fixed96(0x7cb04ff3d675a355ULL,0x000000009b4597e3ULL),
                    fixed96(0x7cb04ff3d675a356ULL,0x000000009b4597e3ULL),
                    fixed96(0x7cb04ff3d675a355ULL,0x000000009b4597e3ULL),
                    fixed96(0x7cb04ff3d675a356ULL,0x000000009b4597e3ULL),
                    2u,true,false
                };
            }
            return {
                fixed96(0xe069bc972dfefab6ULL,0x00000001a61298e1ULL),
                fixed96(0xe069bc972dfefab7ULL,0x00000001a61298e1ULL),
                fixed96(0xe069bc972dfefab6ULL,0x00000001a61298e1ULL),
                fixed96(0xe069bc972dfefab7ULL,0x00000001a61298e1ULL),
                2u,false,false
            };
        case FloatDPTranscendental::logarithm:
            return {
                fixed96(0xd1cf79abc9e3b398ULL,0x00000000b17217f7ULL),
                fixed96(0xd1cf79abc9e3b399ULL,0x00000000b17217f7ULL),
                fixed96(0x0000000000000000ULL,0x0000000200000000ULL),
                fixed96(0x0000000000000000ULL,0x0000000200000000ULL),
                8u,true,true
            };
        default:
            throw std::logic_error("Unknown FloatDP transcendental");
    }
}

std::uint64_t floatdp_correction(FloatDPFixed96 slope,std::uint64_t offset,
                                 unsigned denominator_bits,bool round_upward) {
    using Limb=BigUnsigned::Limb;
    const std::array<Limb,2u> slope_words={slope.low,slope.high};
    const Limb offset_word=offset;
    std::array<Limb,3u> product={0u,0u,0u};
    BigUnsigned::mul_into(
        BigUnsigned::LimbView(slope_words.data(),slope_words.size()),
        BigUnsigned::LimbView(&offset_word,1u),
        BigUnsigned::MutableLimbView(product.data(),product.size()));

    ARIADNE_ASSERT(denominator_bits>0u && denominator_bits<64u);
    std::uint64_t retained=product[0u]>>denominator_bits;
    retained|=product[1u]<<(64u-denominator_bits);
    if(round_upward) {
        const std::uint64_t mask=(std::uint64_t{1}<<denominator_bits)-1u;
        if((product[0u]&mask)!=0u) ++retained;
    }
    ARIADNE_ASSERT(product[2u]==0u);
    return retained;
}

std::array<BigUnsigned::Limb,2u>
floatdp_endpoint(FloatDPFixed96 center,std::int64_t displacement) {
    using Limb=BigUnsigned::Limb;
    const std::array<Limb,2u> input={center.low,center.high};
    std::array<Limb,2u> output={0u,0u};
    const std::uint64_t absolute=displacement<0
        ? std::uint64_t{0}-static_cast<std::uint64_t>(displacement)
        : static_cast<std::uint64_t>(displacement);
    const Limb increment=absolute;
    if(displacement<0) {
        BigUnsigned::sub_into(
            BigUnsigned::LimbView(input.data(),input.size()),
            BigUnsigned::LimbView(&increment,1u),
            BigUnsigned::MutableLimbView(output.data(),output.size()));
    } else {
        BigUnsigned::add_into(
            BigUnsigned::LimbView(input.data(),input.size()),
            BigUnsigned::LimbView(&increment,1u),
            BigUnsigned::MutableLimbView(output.data(),output.size()));
    }
    return output;
}

bool floatdp_round_fixed(std::array<BigUnsigned::Limb,2u> const& magnitude,
                         bool negative,FloatDP::RoundingModeType mode,double& output) {
    const std::uint64_t low=magnitude[0u];
    const std::uint64_t high=magnitude[1u];
    if(high==0u) return false;

    const unsigned high_width=64u-static_cast<unsigned>(std::countl_zero(high));
    const unsigned width=64u+high_width;
    if(width<=53u || width>117u) return false;
    const unsigned shift=width-53u;
    if(shift==0u || shift>=64u) return false;

    std::uint64_t retained=(low>>shift)|(high<<(64u-shift));
    const std::uint64_t discarded_mask=(std::uint64_t{1}<<shift)-1u;
    const bool inexact=(low&discarded_mask)!=0u;
    bool increment=false;

    if(mode==FloatDP::ROUND_TO_NEAREST) {
        const std::uint64_t halfway=std::uint64_t{1}<<(shift-1u);
        const bool above_half=(low&(halfway-1u))!=0u;
        increment=(low&halfway)!=0u && (above_half || (retained&1u)!=0u);
    } else if(mode==FloatDP::ROUND_UPWARD) {
        increment=!negative && inexact;
    } else if(mode==FloatDP::ROUND_DOWNWARD) {
        increment=negative && inexact;
    } else if(mode!=FloatDP::ROUND_TOWARD_ZERO) {
        return false;
    }

    int highest=static_cast<int>(width)-1-static_cast<int>(FLOATDP_FIXED_BITS);
    if(increment) {
        ++retained;
        if(retained==(std::uint64_t{1}<<53u)) {
            retained>>=1u;
            ++highest;
        }
    }

    double result=std::ldexp(static_cast<double>(retained),highest-52);
    output=negative ? -result : result;
    return true;
}

bool floatdp_certified_transcendental(FloatDPTranscendental function,
                                      FloatDP::RoundingModeType mode,
                                      double input,double& output) {
    // Preserve the established inexpensive nearest-mode implementations
    // where the certified near-half path was never selected.
    if(mode==FloatDP::ROUND_TO_NEAREST &&
       function!=FloatDPTranscendental::sine &&
       function!=FloatDPTranscendental::cosine &&
       function!=FloatDPTranscendental::arcsine &&
       function!=FloatDPTranscendental::arccosine &&
       function!=FloatDPTranscendental::hyperbolic_tangent) {
        return false;
    }

    const std::uint64_t encoded=std::bit_cast<std::uint64_t>(input);
    const bool negative=(encoded>>63u)!=0u;
    const std::uint64_t exponent=(encoded>>52u)&0x7ffu;
    if(exponent!=1021u && exponent!=1022u) return false;
    if(function==FloatDPTranscendental::logarithm && negative) return false;

    const unsigned denominator_bits=static_cast<unsigned>(1023u+52u-exponent);
    const std::uint64_t half=std::uint64_t{1}<<(denominator_bits-1u);
    const std::uint64_t magnitude=
        (encoded&((std::uint64_t{1}<<52u)-1u))|(std::uint64_t{1}<<52u);
    const bool offset_negative=magnitude<half;
    const std::uint64_t offset=offset_negative ? half-magnitude : magnitude-half;

    // At scale B=96 the quadratic remainder certificates used by all
    // centers are bounded by their small integer error constants whenever
    // |delta| <= 2^-48.
    constexpr unsigned cutoff=(FLOATDP_FIXED_BITS+1u)/2u;
    if(offset!=0u &&
       (cutoff>denominator_bits ||
        offset>(std::uint64_t{1}<<(denominator_bits-cutoff)))) {
        return false;
    }

    const FloatDPCenter center=floatdp_center(function,negative);
    const std::uint64_t correction_lower=floatdp_correction(
        center.slope_lower,offset,denominator_bits,false);
    const std::uint64_t correction_upper=floatdp_correction(
        center.slope_upper,offset,denominator_bits,true);
    const bool subtract=offset_negative!=center.slope_negative;

    const std::int64_t lower_offset=
        (subtract ? -static_cast<std::int64_t>(correction_upper)
                  :  static_cast<std::int64_t>(correction_lower))
        -static_cast<std::int64_t>(center.error);
    const std::int64_t upper_offset=
        (subtract ? -static_cast<std::int64_t>(correction_lower)
                  :  static_cast<std::int64_t>(correction_upper))
        +static_cast<std::int64_t>(center.error);

    const auto lower_magnitude=floatdp_endpoint(center.lower,lower_offset);
    const auto upper_magnitude=floatdp_endpoint(center.upper,upper_offset);
    double lower_result=0.0,upper_result=0.0;
    if(!floatdp_round_fixed(lower_magnitude,center.result_negative,mode,lower_result) ||
       !floatdp_round_fixed(upper_magnitude,center.result_negative,mode,upper_result)) {
        return false;
    }
    if(std::bit_cast<std::uint64_t>(lower_result)!=
       std::bit_cast<std::uint64_t>(upper_result)) {
        return false;
    }
    output=lower_result;
    return true;
}

} // namespace

FloatDP exp(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental(FloatDPTranscendental::exponential,mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&exp_rnd,mode,x);
}

FloatDP log(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental(FloatDPTranscendental::logarithm,mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&log_rnd,mode,x);
}

FloatDP sin(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental(FloatDPTranscendental::sine,mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&sin_rnd,mode,x);
}

FloatDP cos(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental(FloatDPTranscendental::cosine,mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&cos_rnd,mode,x);
}

FloatDP tan(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental(FloatDPTranscendental::tangent,mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&tan_rnd,mode,x);
}

FloatDP atan(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental(FloatDPTranscendental::arctangent,mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&atan_rnd,mode,x);
}

FloatDP asin(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental(FloatDPTranscendental::arcsine,mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&asin_rnd,mode,x);
}

FloatDP acos(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental(FloatDPTranscendental::arccosine,mode,x.dbl,output)) return FloatDP(output);
    return FloatDP::_apply_rnd(&acos_rnd,mode,x);
}

FloatDP tanh(FloatDP::RoundingModeType mode,FloatDP x) {
    double output;
    if (floatdp_certified_transcendental(FloatDPTranscendental::hyperbolic_tangent,
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
Void FloatDP::set_rounding_mode(RoundingModeType rnd) { Ariadne::set_builtin_rounding_mode(rnd); }
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
    ARIADNE_ASSERT_MSG(z==w,"Cannot cast non-integral value "<<z<<" to an Integer");
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
