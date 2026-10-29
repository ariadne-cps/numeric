/***************************************************************************
 *            numeric/floatmp.cpp
 *
 *  Copyright  2026  Luca Geretti
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

/*! \file numeric/floatmp.cpp
 *  \brief Native multiple-precision floating-point operations.
 */

#include "numeric/floatmp.hpp"
#include "numeric/detail/checked_arithmetic.hpp"
#include "numeric/detail/fixed_point.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include "numeric/dyadic.hpp"

#include "numeric/double.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/rounding.hpp"
#include "numeric/positive.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/decimal.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <climits>
#include <type_traits>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <istream>
#include <limits>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Ariadne {

// All unsigned arithmetic in FloatMP uses the repository's optimized
// 64-bit-limb BigUnsigned; the nonnegative significand is owned directly.
namespace {
bool floatmp_test_bit(BigUnsigned::LimbView n,std::size_t bit) noexcept {
    // round_short/round_scaled_integer retain at least two bits, so their
    // guard/parity bit is inside the canonical input. get_d_general only
    // reaches this helper for discarded<=width, including the smallest tie.
    return ((n[bit/64u]>>(bit%64u))&1u)!=0u;
}
bool floatmp_test_bit(BigUnsigned const& n,std::size_t bit) noexcept {
    return floatmp_test_bit(n.limbs(),bit);
}
bool floatmp_any_bits_below(BigUnsigned::LimbView digits,std::size_t limit) noexcept {
    return BigUnsigned::any_bits_below(digits,limit);
}
bool floatmp_any_bits_below(BigUnsigned const& n,std::size_t limit) noexcept {
    return floatmp_any_bits_below(n.limbs(),limit);
}
std::size_t floatmp_popcount(BigUnsigned const& n) noexcept {
    std::size_t count=0u;
    for(auto word:n.limbs()) count+=static_cast<std::size_t>(std::popcount(word));
    return count;
}
// Binary-ratio conversion and decimal formatting use the same exact
// nearest/even decision: compare twice the remainder with its divisor.
bool floatmp_half_even_increment(BigUnsigned const& remainder,
                                 BigUnsigned const& divisor,bool odd) {
    const int order=compare(shift_left(remainder,1u),divisor);
    return order>0 || (order==0 && odd);
}
BigUnsigned add_big_unsigned(BigUnsigned const&,BigUnsigned const&);
} // namespace
// Round a borrowed, canonical nonzero magnitude to at most 64 bits.
// All bit inspection precedes the write, so the view may alias result.
void FloatMP::round_short(FloatMP& result,BigUnsigned::LimbView magnitude,
                          Exponent scale,Rounding mode) {
    const auto precision=result._precision;
    const auto width=floatmp_exponent_from_bit_length(
        (magnitude.size()-1u)*64u+64u-static_cast<unsigned>(std::countl_zero(magnitude.back())));
    const auto shift=checked_sub(width,static_cast<Exponent>(precision));
    result._exponent=checked_add(scale,shift);
    // The retained significand fits a limb even when the exact input
    // (e.g. a product or power) spans many limbs. Read its leading
    // window and guard/sticky bits without constructing shifted owners.
    std::uint64_t retained;
    bool increment=false;
    if (shift <= 0) {
        retained=magnitude[0u] << static_cast<unsigned>(-shift);
    } else {
        const auto dropped=static_cast<std::size_t>(shift);
        const auto word=dropped/64u;
        const auto offset=static_cast<unsigned>(dropped%64u);
        retained=magnitude[word]>>offset;
        if (offset!=0u && word+1u<magnitude.size()) retained|=magnitude[word+1u]<<(64u-offset);
        switch(mode.value()) {
            case Rounding::ToNearest:
                increment=floatmp_test_bit(magnitude,dropped-1u) &&
                    ((retained&1u)!=0u ||
                     floatmp_any_bits_below(magnitude,dropped-1u));
                break;
            case Rounding::Downward:
                increment=result._negative && floatmp_any_bits_below(magnitude,dropped);
                break;
            case Rounding::Upward:
                increment=!result._negative && floatmp_any_bits_below(magnitude,dropped);
                break;
            default: // TowardZero: the fourth value of the closed rounding type.
                break;
        }
    }
    if (increment) ++retained;
    // A 64-bit carry wraps to zero; narrower carries remain in the limb.
    if (retained==0u || (precision<64u && retained==(std::uint64_t{1}<<precision))) {
        retained=std::uint64_t{1}<<(precision-1u);
        if (result._exponent==std::numeric_limits<Exponent>::max()) {
            result=infinity(precision,result._negative);
            return;
        }
        result._exponent=checked_add(result._exponent,1);
    }
    result._significand=BigUnsigned(retained);
}

// Round a positive native magnitude in place, preserving tie-to-even.
void FloatMP::round_scaled_integer(FloatMP& result, Exponent scale, Rounding mode) {
    const auto precision = result._precision;
    if (precision<=64u) {
        round_short(result,result._significand.limbs(),scale,mode);
        return;
    }
    const Exponent shift = checked_sub(bit_length(result._significand),
                                       static_cast<Exponent>(precision));
    result._exponent = checked_add(scale, shift);
    if (shift == 0) return;
    if (shift < 0) {
        result._significand=shift_left(result._significand,
            static_cast<std::size_t>(shift_amount(checked_sub(0,shift))));
        return;
    }
    const std::size_t dropped=static_cast<std::size_t>(shift_amount(shift));
    const bool inexact=floatmp_any_bits_below(result._significand,dropped);
    bool increment = false;
    if (inexact) {
        switch (mode.value()) {
            case Rounding::ToNearest:
                increment=floatmp_test_bit(result._significand,dropped-1u) &&
                    (floatmp_test_bit(result._significand,dropped) ||
                     floatmp_any_bits_below(result._significand,dropped-1u));
                break;
            case Rounding::Downward:
                increment = result._negative;
                break;
            case Rounding::Upward:
                increment = !result._negative;
                break;
            default: // TowardZero: the fourth value of the closed rounding type.
                break;
        }
    }
    result._significand=shift_right(result._significand,dropped);
    if (increment) result._significand=add_big_unsigned(
        result._significand,BigUnsigned(1u));
    if (bit_length(result._significand)>static_cast<Exponent>(precision)) {
        result._significand=shift_right(result._significand,1u);
        // Crossing the maximal finite exponent by a rounded carry yields infinity.
        if (result._exponent == std::numeric_limits<Exponent>::max()) {
            result = infinity(precision, result._negative);
            return;
        }
        result._exponent = checked_add(result._exponent, 1);
    }
}

namespace {
// BigUnsigned's hidden-friend add relies on ADL, which is suppressed
// by FloatMP::add during lookup inside FloatMP member functions.
BigUnsigned add_big_unsigned(BigUnsigned const& a, BigUnsigned const& b) {
    return add(a,b);
}
// ADL for BigUnsigned's hidden friends is suppressed inside FloatMP::mul.
BigUnsigned mul_big_unsigned(BigUnsigned const& a, BigUnsigned const& b) {
    return mul(a,b);
}
BigUnsigned sub_big_unsigned(BigUnsigned const& a, BigUnsigned const& b) {
    return sub(a,b);
}
} // namespace

// D64: exact division and remainder now use the shared native
// BigUnsigned kernels, including for the old <=64-bit fast path.
// GMP is temporarily confined to the boundary wrapper and FloatMP
// significand storage until the remaining storage migration.
FloatMP FloatMP::from_native_ratio(Integer const& numerator,
                                  Integer const& denominator,
                                  Exponent scale, unsigned long precision,
                                  Rounding mode, bool negative_zero,
                                  bool invert_sign) {
    return from_unsigned_ratio(numerator.magnitude(),denominator.magnitude(),
        scale,precision,mode,
        (numerator.negative()!=denominator.negative())!=invert_sign,negative_zero);
}

FloatMP FloatMP::from_unsigned_ratio(BigUnsigned const& n,BigUnsigned const& d,
                                    Exponent scale,unsigned long precision,
                                    Rounding mode,bool negative,bool negative_zero) {
    if (d.is_zero()) {
        throw std::domain_error("FloatMP zero denominator");
    }
    FloatMP result(precision);
    if (n.is_zero()) {
        result.set_zero(negative_zero);
        return result;
    }
    result._negative=negative;

    // Power-of-two denominators are exact binary scalings. Preserve
    // the fast path used by the old mpq_t constructor without creating a
    // GMP rational or executing a general multi-limb division.
    bool denominator_is_power_of_two = true;
    bool saw_nonzero_limb = false;
    for (auto limb : d.limbs()) {
        if (limb == 0u) continue;
        if (saw_nonzero_limb || !std::has_single_bit(limb)) {
            denominator_is_power_of_two = false;
            break;
        }
        saw_nonzero_limb = true;
    }
    if (denominator_is_power_of_two) {
        const Exponent denominator_shift =
            floatmp_exponent_from_bit_length(d.bit_length()-1u);
        result._significand=n;
        round_scaled_integer(result,checked_sub(scale,denominator_shift),mode);
        return result;
    }

    if (precision<=64u && n.size()==1u && d.size()==1u) {
        // Normalize a short ratio using scalar metadata, then use the
        // shared exact 128/64 limb divider. The scaled numerator needs at
        // most p+bit_length(d) bits; the quotient always fits one limb.
        const std::uint64_t numerator=n.limb(0u),denominator=d.limb(0u);
        Exponent k=static_cast<Exponent>(std::countl_zero(denominator))-
                   static_cast<Exponent>(std::countl_zero(numerator));
        if (k>=0 ? numerator<(denominator<<static_cast<unsigned>(k))
                 : (numerator<<static_cast<unsigned>(-k))<denominator) --k;
        const Exponent shift=static_cast<Exponent>(precision-1u)-k;
        std::array<std::uint64_t,2u> dividend{},quotient{};
        std::uint64_t divisor=denominator;
        if (shift>=64) {
            dividend[1u]=numerator<<static_cast<unsigned>(shift-64);
        } else if (shift>=0) {
            const auto s=static_cast<unsigned>(shift);
            dividend[0u]=numerator<<s;
            if (s!=0u) dividend[1u]=numerator>>(64u-s);
        } else {
            dividend[0u]=numerator;
            divisor<<=static_cast<unsigned>(-shift);
        }
        const std::uint64_t remainder=BigUnsigned::divrem_limb_into(
            BigUnsigned::LimbView(dividend.data(),dividend[1u]!=0u ? 2u : 1u),
            divisor,quotient);
        std::uint64_t retained=quotient[0u];
        bool increment=false;
        if (remainder!=0u) {
            switch(mode.value()) {
                case Rounding::ToNearest:
                    increment=remainder>divisor-remainder ||
                        (remainder==divisor-remainder && (retained&1u)!=0u);
                    break;
                case Rounding::Downward: increment=negative; break;
                case Rounding::Upward: increment=!negative; break;
                default: // TowardZero: the fourth value of the closed rounding type.
                break;
            }
        }
        result._exponent=checked_add(scale,k-static_cast<Exponent>(precision-1u));
        if (increment) ++retained;
        // The normalized quotient is positive. At p=64 it cannot reach
        // UINT64_MAX: for k<0, its gap below 2^64 is >=2^64/d>1.
        // For k>=0, q=UINT64_MAX would force S=2^(63-k)>=d;
        // the gap is then >=S/d>=1, with equality only for the excluded
        // power-of-two denominator d=S. Increment therefore cannot wrap.
        if (precision<64u && retained==(std::uint64_t{1}<<precision)) {
            retained=std::uint64_t{1}<<(precision-1u);
            if (result._exponent==std::numeric_limits<Exponent>::max()) return infinity(precision,negative);
            result._exponent=checked_add(result._exponent,1);
        }
        result._significand=BigUnsigned(retained);
        return result;
    }

    Exponent k=checked_sub(floatmp_exponent_from_bit_length(n.bit_length()),
                           floatmp_exponent_from_bit_length(d.bit_length()));
    if (k>=0) {
        if (compare(n,shift_left(d,shift_amount(k)))<0) {
            k=checked_sub(k,1);
        }
    } else {
        if (compare(shift_left(n,shift_amount(checked_sub(0,k))),d)<0) {
            k=checked_sub(k,1);
        }
    }

    Exponent const shift=checked_sub(static_cast<Exponent>(precision-1u),k);
    BigUnsigned dividend;
    BigUnsigned divisor;
    if (shift>=0) {
        dividend=shift_left(n,shift_amount(shift));
        divisor=d;
    } else {
        dividend=n;
        divisor=shift_left(d,shift_amount(checked_sub(0,shift)));
    }
    auto [quotient,remainder]=divrem(dividend,divisor);
    bool increment=false;
    if (!remainder.is_zero()) {
        switch (mode.value()) {
            case Rounding::ToNearest:
                increment=floatmp_half_even_increment(remainder,divisor,(quotient.limb(0u)&1u)!=0u);
                break;
            case Rounding::Downward: increment=result._negative; break;
            case Rounding::Upward: increment=!result._negative; break;
            default: // TowardZero: the fourth value of the closed rounding type.
                break;
        }
    }
    if (increment) quotient=add_big_unsigned(quotient,BigUnsigned(1u));
    result._exponent=checked_add(
        scale,checked_sub(k,static_cast<Exponent>(precision-1u)));
    if (quotient.bit_length()>precision) {
        quotient=shift_right(quotient,1u);
        if (result._exponent==std::numeric_limits<Exponent>::max()) {
            return infinity(precision,result._negative);
        }
        result._exponent=checked_add(result._exponent,1);
    }
    result._significand=std::move(quotient);
    return result;
}

FloatMP FloatMP::negate(FloatMP x) {
    if (!x.is_nan()) x._negative = !x._negative;
    return x;
}

FloatMP FloatMP::add(FloatMP const& a, FloatMP const& b,
unsigned long p, Rounding mode) {
    return add_sub(a, b, p, mode, false);
}

FloatMP FloatMP::sub(FloatMP const& a, FloatMP const& b,
unsigned long p, Rounding mode) {
    return add_sub(a, b, p, mode, true);
}

FloatMP FloatMP::add_sub(FloatMP const& a, FloatMP const& b,
unsigned long p, Rounding mode, bool subtract) {
    const bool b_negative = b._negative != subtract;
    if (a.is_nan() || b.is_nan()) return nan(p);
    if (a.is_inf() || b.is_inf()) {
        if (a.is_inf() && b.is_inf() && a._negative != b_negative) return nan(p);
        return infinity(p, a.is_inf() ? a._negative : b_negative);
    }
    if (a.is_zero() && b.is_zero()) {
        return zero(p, a._negative == b_negative ? a._negative : mode == Rounding::Downward);
    }
    if (a.is_zero()) return FloatMP(subtract ? negate(b) : b, mode, PrecisionType(p));
    if (b.is_zero()) return FloatMP(a, mode, PrecisionType(p));

    // Align the exact integers for ordinary gaps; for immense gaps use
    // a rigorously small, signed sticky term before a single rounding.
    // Unlike the old approximation shortcut, this preserves directed
    // subtraction and ties, including when the tiny operand opposes the
    // large operand.
    const bool a_high = a._exponent >= b._exponent;
    FloatMP const& high = a_high ? a : b;
    FloatMP const& low = a_high ? b : a;
    const bool high_negative = a_high ? a._negative : b_negative;
    const bool low_negative = a_high ? b_negative : a._negative;
    // Unsigned subtraction handles even INT64_MAX - INT64_MIN exactly.
    const std::uint64_t gap = static_cast<std::uint64_t>(high._exponent) -
                              static_cast<std::uint64_t>(low._exponent);
    // This conservative gap bound places the smaller operand below all
    // noncoincident rounding boundaries even when p > high.precision.
    const std::uint64_t high_bits = static_cast<std::uint64_t>(high._precision);
    const std::uint64_t low_bits = static_cast<std::uint64_t>(low._precision);
    const bool far_apart = gap > high_bits && gap - high_bits > low_bits &&
        gap - high_bits - low_bits > static_cast<std::uint64_t>(p) + 2u;
    if (far_apart) {
        // Preserve the certified signed-sticky approximation using native
        // Integer for signed arithmetic and BigUnsigned for exact shifts.
        const Exponent padding=checked_add(
            checked_add(static_cast<Exponent>(high._precision),
                        static_cast<Exponent>(p)),2);
        const auto shift=static_cast<std::size_t>(shift_amount(padding));
        Integer exact=Integer::from_magnitude(
            shift_left(high._significand,shift),high_negative);
        exact=low_negative ? exact-Integer(1) : exact+Integer(1);
        const Integer denominator=Integer::from_magnitude(
            shift_left(BigUnsigned(1u),shift),false);
        return from_native_ratio(exact,denominator,high._exponent,p,mode);
    }

    // For short results, shared shifted-combine kernels avoid an owning
    // shifted operand, including product/addend alignment in a fused op.
    if (p<=64u && gap<64u) {
        FloatMP result(p);
        result._negative=high_negative;
        const auto shift=static_cast<unsigned>(gap);
        if (high_negative==low_negative) {
            result._significand=BigUnsigned::add_shifted(
                high._significand.limbs(),shift,low._significand.limbs(),0u);
        } else {
            // A larger exponent need not mean a larger magnitude when
            // input precisions differ. Compare without a shifted owner.
            const auto high_width=static_cast<std::uint64_t>(high._significand.bit_length())+gap;
            const auto low_width=static_cast<std::uint64_t>(low._significand.bit_length());
            bool high_larger=high_width>low_width;
            if (high_width==low_width) {
                // Equal widths are uncommon for mixed precisions. The
                // shared exact comparison retains cancellation semantics.
                high_larger=compare(shift_left(high._significand,shift),low._significand)>=0;
            }
            result._significand=high_larger
                ? BigUnsigned::sub_shifted(high._significand.limbs(),shift,low._significand.limbs(),0u)
                : BigUnsigned::sub_shifted(low._significand.limbs(),0u,high._significand.limbs(),shift);
            if (!high_larger) result._negative=!result._negative;
        }
        if (result._significand.is_zero()) result.set_zero(mode==Rounding::Downward);
        else round_scaled_integer(result,low._exponent,mode);
        return result;
    }

    // For ordinary larger gaps, align and combine positive magnitudes natively.
    // Preserve the independent sign metadata and a single final rounding.
    FloatMP result(p);
    result._negative=high_negative;
    BigUnsigned aligned=shift_left(high._significand,
        static_cast<std::size_t>(
            shift_amount(checked_sub(high._exponent,low._exponent))));
    if(a._negative==b_negative) {
        aligned=add_big_unsigned(aligned,low._significand);
    } else {
        if(compare(aligned,low._significand)<0) {
            aligned=sub_big_unsigned(low._significand,aligned);
            result._negative=!result._negative;
        } else {
            aligned=sub_big_unsigned(aligned,low._significand);
        }
    }
    if(aligned.is_zero()) {
        result.set_zero(mode==Rounding::Downward);
        return result;
    }
    result._significand=std::move(aligned);
    round_scaled_integer(result,low._exponent,mode);
    return result;
}

FloatMP FloatMP::mul(FloatMP const& a, FloatMP const& b,
unsigned long p, Rounding mode) {
    if (a.is_nan() || b.is_nan() || (a.is_inf() && b.is_zero()) || (a.is_zero() && b.is_inf())) return nan(p);
    const bool neg = a._negative != b._negative;
    if (a.is_inf() || b.is_inf()) return infinity(p, neg);
    if (a.is_zero() || b.is_zero()) return zero(p, neg);
    FloatMP result(p);
    result._negative = neg;
    if (p<=64u && a._significand.size()<=2u && b._significand.size()<=2u) {
        std::array<BigUnsigned::Limb,4u> product;
        const auto count=a._significand.size()+b._significand.size();
        BigUnsigned::mul_into(a._significand.limbs(),b._significand.limbs(),
            BigUnsigned::MutableLimbView(product.data(),count));
        auto used=count;
        while (product[used-1u]==0u) --used;
        round_short(result,BigUnsigned::LimbView(product.data(),used),
                    checked_add(a._exponent,b._exponent),mode);
        return result;
    }
    result._significand=mul_big_unsigned(a._significand,b._significand);
    round_scaled_integer(result,checked_add(a._exponent,b._exponent),mode);
    return result;
}

FloatMP FloatMP::div(FloatMP const& a, FloatMP const& b,
unsigned long p, Rounding mode) {
    if (a.is_nan() || b.is_nan() || (a.is_zero() && b.is_zero()) || (a.is_inf() && b.is_inf())) return nan(p);
    const bool neg = a._negative != b._negative;
    if (a.is_inf() || b.is_zero()) return infinity(p, neg);
    if (a.is_zero() || b.is_inf()) return zero(p, neg);
    // Borrow the unsigned magnitudes; no signed Integer copies are needed.
    return from_unsigned_ratio(a._significand,b._significand,
                             checked_sub(a._exponent,b._exponent),
                             p,mode,neg);
}

FloatMP FloatMP::fma(FloatMP const& a, FloatMP const& b,
FloatMP const& c, unsigned long p, Rounding mode) {
    if (a.is_nan() || b.is_nan() || c.is_nan() || (a.is_zero() && b.is_inf()) || (a.is_inf() && b.is_zero())) return nan(p);
    const bool neg = a._negative != b._negative;
    if (a.is_inf() || b.is_inf()) {
        if (c.is_inf() && c._negative != neg) return nan(p);
        return infinity(p, neg);
    }
    if (c.is_inf()) return infinity(p, c._negative);
    if (a.is_zero() || b.is_zero()) return add(zero(p, neg), c, p, mode);
    // Short fused a*b+c: align the exact 128-bit product and the exact
    // one-limb addend directly, then round ONCE. This avoids creating a
    // full-precision temporary FloatMP and re-entering add_sub's generic
    // special-value and exponent dispatch. All arithmetic remains exact.
    if (p<=64u && a._significand.size()==1u &&
        b._significand.size()==1u && c._significand.size()==1u) {
        const Exponent product_scale=checked_add(a._exponent,b._exponent);
        const bool product_high=product_scale>=c._exponent;
        const std::uint64_t gap=product_high
            ? static_cast<std::uint64_t>(product_scale)-
              static_cast<std::uint64_t>(c._exponent)
            : static_cast<std::uint64_t>(c._exponent)-
              static_cast<std::uint64_t>(product_scale);
        if (gap<64u) {
            std::array<BigUnsigned::Limb,2u> product{};
            BigUnsigned::mul_into(a._significand.limbs(),
                b._significand.limbs(),product);
            const auto used=product[1u]==0u ? 1u : 2u;
            const BigUnsigned::LimbView product_words(product.data(),used);
            const auto addend_words=c._significand.limbs();
            const unsigned product_shift=product_high
                ? static_cast<unsigned>(gap) : 0u;
            const unsigned addend_shift=product_high
                ? 0u : static_cast<unsigned>(gap);
            FloatMP result(p);
            result._negative=neg;
            if (neg==c._negative) {
                result._significand=BigUnsigned::add_shifted(
                    product_words,product_shift,addend_words,addend_shift);
            } else {
                // Different signs can cancel. Decide the exact magnitude
                // ordering before subtracting; equal leading widths require
                // a full limb comparison, never a rounded approximation.
                const std::size_t product_width=
                    (used-1u)*64u+64u-static_cast<unsigned>(
                        std::countl_zero(product[used-1u]))+product_shift;
                const std::size_t addend_width=
                    c._significand.bit_length()+addend_shift;
                const bool product_larger=product_width>addend_width ||
                    (product_width==addend_width &&
                     compare(shift_left(BigUnsigned(product_words),product_shift),
                             shift_left(c._significand,addend_shift))>=0);
                result._significand=product_larger
                    ? BigUnsigned::sub_shifted(product_words,product_shift,
                                               addend_words,addend_shift)
                    : BigUnsigned::sub_shifted(addend_words,addend_shift,
                                               product_words,product_shift);
                if (!product_larger) result._negative=!result._negative;
            }
            if (result._significand.is_zero()) {
                result.set_zero(mode==Rounding::Downward);
            } else {
                round_scaled_integer(result,
                    product_high ? c._exponent : product_scale,mode);
            }
            return result;
        }
    }
    // Wide operands and large gaps retain the proven general single-round
    // fused path, including signed-sticky alignment and mixed precisions.
    BigUnsigned product=mul_big_unsigned(a._significand,b._significand);
    FloatMP exact(static_cast<unsigned long>(product.bit_length()));
    exact._significand=std::move(product);
    exact._negative = neg;
    exact._exponent = checked_add(a._exponent, b._exponent);
    return add(exact, c, p, mode);
}

// Integer owns a native signed magnitude; no mpz_t import is required,
// including for INT64_MIN.
FloatMP::Float(std::int64_t n, PrecisionType pr)
    : Float(Integer(n), ROUND_TO_NEAREST, pr) { }

namespace {
thread_local FloatMP::RoundingModeType floatmp_current_rounding =
    FloatMP::ROUND_TO_NEAREST;
thread_local std::size_t floatmp_default_precision = 53;
}

FloatMP::RoundingModeType FloatMP::get_rounding_mode() noexcept {
    return floatmp_current_rounding;
}

FloatMP::PrecisionType FloatMP::get_default_precision() {
    return PrecisionType(floatmp_default_precision);
}

void FloatMP::set_default_precision(PrecisionType precision) {
    const FloatMP validated(precision);
    floatmp_default_precision = validated.precision().bits();
}

void FloatMP::set_rounding_mode(RoundingModeType mode) {
    // Every defined Rounding value denotes one of the four supported modes.
    floatmp_current_rounding=mode;
}

FloatMP::Rounding::Rounding(RoundToNearest const&) noexcept : _value(ToNearest) { }
FloatMP::Rounding::Rounding(RoundDownward const&) noexcept : _value(Downward) { }
FloatMP::Rounding::Rounding(RoundUpward const&) noexcept : _value(Upward) { }
FloatMP::Rounding::Rounding(RoundTowardZero const&) noexcept : _value(TowardZero) { }
FloatMP::Rounding::Rounding(RoundApproximately const&) noexcept : _value(ToNearest) { }
FloatMP::Rounding::Rounding(CurrentRoundingMode const&) noexcept
    : _value(FloatMP::get_rounding_mode().value()) { }

FloatMP::Rounding::Rounding(Ariadne::Rounding const& mode) noexcept
    : _value(static_cast<BuiltinRoundingModeType>(mode)==Ariadne::ROUND_DOWNWARD ? Downward
             : static_cast<BuiltinRoundingModeType>(mode)==Ariadne::ROUND_UPWARD ? Upward
             : static_cast<BuiltinRoundingModeType>(mode)==Ariadne::ROUND_TOWARD_ZERO ? TowardZero
             : ToNearest) { }

// Forward current-rounding operations to the existing native implementation.
FloatMP sqr(CurrentRoundingMode, FloatMP const& x) { return sqr(FloatMP::get_rounding_mode(),x); }
FloatMP rec(CurrentRoundingMode, FloatMP const& x) { return rec(FloatMP::get_rounding_mode(),x); }
FloatMP add(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2) { return add(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP sub(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2) { return sub(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP mul(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2) { return mul(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP div(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2) { return div(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP fma(CurrentRoundingMode, FloatMP const& x1, FloatMP const& x2, FloatMP const& x3) { return fma(FloatMP::get_rounding_mode(),x1,x2,x3); }
FloatMP pow(CurrentRoundingMode, FloatMP const& x, Nat m) { return pow(FloatMP::get_rounding_mode(),x,m); }
FloatMP pow(CurrentRoundingMode, FloatMP const& x, Int n) { return pow(FloatMP::get_rounding_mode(),x,n); }
FloatMP sqrt(CurrentRoundingMode, FloatMP const& x) { return sqrt(FloatMP::get_rounding_mode(),x); }
FloatMP exp(CurrentRoundingMode, FloatMP const& x) { return exp(FloatMP::get_rounding_mode(),x); }
FloatMP log(CurrentRoundingMode, FloatMP const& x) { return log(FloatMP::get_rounding_mode(),x); }
FloatMP sin(CurrentRoundingMode, FloatMP const& x) { return sin(FloatMP::get_rounding_mode(),x); }
FloatMP cos(CurrentRoundingMode, FloatMP const& x) { return cos(FloatMP::get_rounding_mode(),x); }
FloatMP tan(CurrentRoundingMode, FloatMP const& x) { return tan(FloatMP::get_rounding_mode(),x); }
FloatMP tanh(CurrentRoundingMode, FloatMP const& x) { return tanh(FloatMP::get_rounding_mode(),x); }
FloatMP asin(CurrentRoundingMode, FloatMP const& x) { return asin(FloatMP::get_rounding_mode(),x); }
FloatMP acos(CurrentRoundingMode, FloatMP const& x) { return acos(FloatMP::get_rounding_mode(),x); }
FloatMP atan(CurrentRoundingMode, FloatMP const& x) { return atan(FloatMP::get_rounding_mode(),x); }
FloatMP FloatMP::pi(CurrentRoundingMode, MultiplePrecision pr) { return pi(FloatMP::get_rounding_mode(),pr); }
FloatMP add(CurrentRoundingMode, FloatMP const& x1, FloatDP const& x2) { return add(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP sub(CurrentRoundingMode, FloatMP const& x1, FloatDP const& x2) { return sub(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP mul(CurrentRoundingMode, FloatMP const& x1, FloatDP const& x2) { return mul(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP div(CurrentRoundingMode, FloatMP const& x1, FloatDP const& x2) { return div(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP add(CurrentRoundingMode, FloatDP const& x1, FloatMP const& x2) { return add(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP sub(CurrentRoundingMode, FloatDP const& x1, FloatMP const& x2) { return sub(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP mul(CurrentRoundingMode, FloatDP const& x1, FloatMP const& x2) { return mul(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP div(CurrentRoundingMode, FloatDP const& x1, FloatMP const& x2) { return div(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP add(CurrentRoundingMode, FloatMP const& x1, ExactDouble x2) { return add(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP sub(CurrentRoundingMode, FloatMP const& x1, ExactDouble x2) { return sub(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP mul(CurrentRoundingMode, FloatMP const& x1, ExactDouble x2) { return mul(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP div(CurrentRoundingMode, FloatMP const& x1, ExactDouble x2) { return div(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP add(CurrentRoundingMode, ExactDouble x1, FloatMP const& x2) { return add(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP sub(CurrentRoundingMode, ExactDouble x1, FloatMP const& x2) { return sub(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP mul(CurrentRoundingMode, ExactDouble x1, FloatMP const& x2) { return mul(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP div(CurrentRoundingMode, ExactDouble x1, FloatMP const& x2) { return div(FloatMP::get_rounding_mode(),x1,x2); }

FloatMP::Float(FloatDP const& value, PrecisionType precision)
    : Float(ExactDouble(value.get_d()), ROUND_TO_NEAREST, precision) { }

FloatMP::Float(FloatDP const& value, RoundingModeType mode, PrecisionType precision)
    : Float(ExactDouble(value.get_d()), mode, precision) { }

FloatMP::Float(FloatMP const& value, PrecisionType precision)
    : Float(value, ROUND_TO_NEAREST, precision) { }

FloatMP::Float(FloatMP const& value, RoundingModeType mode, PrecisionType precision)
    : Float(precision) {
    if (value._precision == _precision) {
        *this = value;
        return;
    }
    if (value.is_nan()) {
        *this = nan(precision);
    } else if (value.is_inf()) {
        *this = infinity(static_cast<unsigned long>(precision.bits()), value.negative());
    } else if (value.is_zero()) {
        *this = zero(static_cast<unsigned long>(precision.bits()), value.negative());
    } else {
        // Rescale the existing binary significand without materializing
        // significand * 2^exponent as an enormous Rational. This is
        // essential at the native exponent limits.
        _significand=value._significand;
        _negative=value._negative;
        round_scaled_integer(*this,value._exponent,mode);
    }
}

FloatMP& FloatMP::operator=(Integer const& n) {
    return *this = FloatMP(n, get_rounding_mode(), precision());
}

FloatMP& FloatMP::operator=(Dyadic const& n) {
    // Preserve the exact-assignment contract of FloatMP. Constructing the
    // replacement first also leaves *this unchanged if it is inexact.
    return *this = FloatMP(n, precision());
}

FloatMP& FloatMP::operator=(TwoExp const& n) {
    return *this = FloatMP(n, precision());
}

FloatMP& FloatMP::operator=(ExactDouble const& n) {
    return *this = FloatMP(n, get_rounding_mode(), precision());
}

FloatMP& FloatMP::assign_unsigned(std::uint64_t n) {
    return *this = FloatMP(Integer(n), get_rounding_mode(), precision());
}

FloatMP::operator Dyadic() const {
    if (is_nan()) {
        return Dyadic::nan();
    }
    if (is_inf()) {
        return Dyadic::inf(negative() ? Sign::NEGATIVE : Sign::POSITIVE);
    }
    // Export the exact native binary significand/exponent.
    BigUnsigned magnitude=_significand;
    if (magnitude.is_zero()) { return Dyadic(0); }
    if (_exponent>=0) {
        return Dyadic(Integer::from_magnitude(
            shift_left(magnitude,shift_amount(_exponent)),_negative));
    }
    const auto shift=shift_amount(checked_sub(0,_exponent));
    return Dyadic(Integer::from_magnitude(std::move(magnitude),_negative),shift);
}

Integer cast_integer(FloatMP const& value) {
    if (!value.is_finite()) {
        throw std::domain_error("FloatMP cannot cast a nonfinite value to Integer");
    }
    const Rational exact = static_cast<Rational>(value);
    if (exact.denominator()!=1) {
        throw std::domain_error("FloatMP is not an integer");
    }
    return exact.numerator();
}

Comparison cmp(FloatMP const& a, Rational const& b) {
    if (a.is_nan() || is_nan(b)) {
        return Comparison::INCOMPARABLE;
    }
    return cmp(static_cast<Rational>(a), b);
}

Comparison cmp(Rational const& a, FloatMP const& b) {
    if (b.is_nan() || is_nan(a)) {
        return Comparison::INCOMPARABLE;
    }
    return cmp(a, static_cast<Rational>(b));
}

FloatMP::Float(ExactDouble const& value, PrecisionType precision)
    : Float(value, ROUND_TO_NEAREST, precision) { }

FloatMP::Float(Dyadic const& value, PrecisionType precision)
    : Float(Rational(value), ROUND_TO_NEAREST, precision) {
    // Unlike the explicitly rounded constructor, this conversion must be
    // exact. In particular, a Dyadic with more significant bits than the
    // requested precision cannot be silently rounded.
    if (!this->is_nan() && cmp(*this, Rational(value)) != Comparison::EQUAL) {
        throw std::runtime_error("FloatMP: Dyadic is not exactly representable at the requested precision");
    }
}

FloatMP::Float(TwoExp const& value, PrecisionType precision)
    : Float(Dyadic(value), precision) { }

FloatMP::Float(std::string const& value, PrecisionType precision)
    : Float(Dyadic(value), precision) { }

FloatMP::Float(ExactDouble const& value, RoundingModeType mode,
               PrecisionType precision)
    : Float(Rational(value), mode, precision) {
    // Rational has one canonical zero. Preserve the source binary sign
    // before it is lost at this otherwise exact conversion boundary.
    if (value.get_d()==0.0) set_zero(std::signbit(value.get_d()));
}

FloatMP::Float(Dyadic const& value, RoundingModeType mode,
               PrecisionType precision)
    : Float(Rational(value), mode, precision) { }

FloatMP::Float(Decimal const& value, RoundingModeType mode,
               PrecisionType precision)
    : Float(Rational(value), mode, precision) { }

void FloatMP::set_precision(PrecisionType precision) {
    if (precision.bits() == _precision) {
        return;
    }
    if (is_nan()) {
        *this = nan(precision);
    } else if (is_inf()) {
        *this = infinity(static_cast<unsigned long>(precision.bits()), _negative);
    } else if (is_zero()) {
        *this = zero(static_cast<unsigned long>(precision.bits()), _negative);
    } else {
        *this = FloatMP(static_cast<Rational>(*this), ROUND_TO_NEAREST, precision);
    }
}


Nat FloatMP::output_places=20;
void FloatMP::set_output_places(Nat places) { output_places=places; }

namespace {
// Format an exact binary rational to fixed decimal places using integer
// division and a single correctly directed decimal rounding step.
std::string decimal_fixed(FloatMP const& x, unsigned int places,
                          FloatMP::RoundingModeType mode) {
    if (x.is_nan()) return "nan";
    if (x.is_inf()) return x.negative() ? "-inf" : "inf";
    // Preserve the exact dyadic value rather than constructing a
    // canonical Rational and a temporary GMP rational. Decimal digits
    // are obtained by one native unsigned division and a single exact
    // directed rounding decision.
    BigUnsigned numerator=x.significand();
    BigUnsigned denominator(1u);
    const auto binary_exponent=x.exponent();
    if(binary_exponent>0) {
        numerator=shift_left(numerator,
            static_cast<std::size_t>(binary_exponent));
    } else if(binary_exponent<0) {
        // Avoid negating INT64_MIN; preserve the checked shift domain.
        const auto distance=static_cast<std::uint64_t>(
            -(binary_exponent+1))+1u;
        const auto address_distance=Detail::checked_unsigned_cast<std::size_t>(distance,
            "FloatMP decimal exponent shift out of range");
        denominator=shift_left(denominator,address_distance);
    }
    const bool neg=x.negative();
    const BigUnsigned scale=pow(BigUnsigned(10u),
        static_cast<std::uint64_t>(places));
    auto [quotient,remainder]=divrem(mul(numerator,scale),denominator);
    if (!remainder.is_zero()) {
        bool increase=false;
        switch(mode.value()) {
            case FloatMP::Rounding::ToNearest:
                increase=floatmp_half_even_increment(remainder,denominator,
                    (quotient.limb(0u)&BigUnsigned::Limb{1})!=0u);
                break;
            case FloatMP::Rounding::Downward: increase=neg; break;
            case FloatMP::Rounding::Upward: increase=!neg; break;
            default: // TowardZero: the fourth value of the closed rounding type.
                break;
        }
        if(increase) quotient=add(quotient,BigUnsigned(1u));
    }
    std::string digits=to_decimal(quotient);
    if(places>0) {
        if(digits.size()<=places) digits.insert(0,places+1-digits.size(),'0');
        digits.insert(digits.size()-places,1,'.');
    } else digits.push_back('.');
    if(neg) digits.insert(digits.begin(),'-');
    return digits;
}
}

std::string print(FloatMP const& x, DecimalPlaces places, FloatMP::RoundingModeType mode) {
    return decimal_fixed(x,static_cast<unsigned int>(places),mode);
}

std::string print(FloatMP const& x, DecimalPrecision figures, FloatMP::RoundingModeType mode) {
    if(!x.is_finite() || x.is_zero()) return decimal_fixed(x,0,mode);
    double absolute;
    {
        // A directed binary64 conversion can turn a tiny source into the
        // smallest subnormal, giving it the wrong decimal scale. Use a
        // nearest conversion for this estimate and restore the caller's mode.
        RoundingSession<FloatMP> session;
        session.set(FloatMP::ROUND_TO_NEAREST);
        absolute=std::fabs(x.get_d());
    }
    // DecimalPrecision is a count of significant figures rather than places.
    // For values outside the range of double, derive the base-10 exponent
    // from the full integer significand and binary exponent.
    long long magnitude=0;
    if(std::isfinite(absolute) && absolute>0) {
        magnitude=static_cast<long long>(std::floor(std::log10(absolute)));
    } else {
        const double estimated=(static_cast<double>(x.significand().bit_length())-1.0
                               +static_cast<double>(x.exponent()))*0.3010299956639812;
        // Both the significand width and exponent fit Exponent. Their
        // sum times log10(2) lies strictly inside the signed 64-bit range.
        // An int is insufficient for valid extreme binary exponents.
        magnitude=static_cast<long long>(std::floor(estimated));
    }
    const long long places=static_cast<long long>(static_cast<unsigned int>(figures))-1-magnitude;
    const auto decimal_places=places>0 ? Detail::checked_unsigned_cast<unsigned int>(
        static_cast<unsigned long long>(places),"FloatMP decimal places out of range") : 0u;
    return decimal_fixed(x,decimal_places,mode);
}

OutputStream& write(OutputStream& out, FloatMP const& x, DecimalPlaces places, FloatMP::RoundingModeType mode) {
    return out<<print(x,places,mode);
}

OutputStream& write(OutputStream& out, FloatMP const& x, DecimalPrecision figures, FloatMP::RoundingModeType mode) {
    return out<<print(x,figures,mode);
}

OutputStream& repr(OutputStream& out, FloatMP const& x, FloatMP::RoundingModeType mode) {
    return write(out,x,DecimalPrecision(static_cast<Nat>(x.precision().bits()/3+4)),mode);
}

OutputStream& repr(OutputStream& out, FloatMP const& x) {
    return repr(out,x,FloatMP::get_rounding_mode());
}

std::string FloatMP::literal() const {
    return literal(get_rounding_mode());
}

std::string FloatMP::literal(RoundingModeType mode) const {
    std::ostringstream out;
    repr(out,*this,mode);
    return out.str();
}

OutputStream& operator<<(OutputStream& os, FloatMP const& value) {
    return os << static_cast<Rational>(value);
}

InputStream& operator>>(InputStream& is, FloatMP& value) {
    // Follow the established token boundaries, leaving punctuation for the caller.
    std::string token;
    int c = is.get();
    while (c == ' ' || c == '\t' || c == '\n') {
        c = is.get();
    }
    while ((c >= '0' && c <= '9') || c == '+' || c == '-' ||
           c == '.' || c == 'e' || c == 'E') {
        token.push_back(static_cast<char>(c));
        c = is.get();
    }
    const bool at_eof = c == std::char_traits<char>::eof();
    if (!at_eof) {
        is.putback(static_cast<char>(c));
    }

    try {
        const std::size_t exponent_start = token.find_first_of("eE");
        const String mantissa = token.substr(0, exponent_start);
        bool has_digit = false;
        for (char digit : mantissa) {
            // Tokenization permits only digits, signs and '.' before e/E.
            if (digit >= '0') {
                has_digit = true;
                break;
            }
        }
        if (!has_digit) {
            throw std::invalid_argument("FloatMP invalid decimal literal");
        }
        Rational exact = static_cast<Rational>(Decimal(mantissa));
        if (exponent_start != std::string::npos) {
            const std::string suffix = token.substr(exponent_start + 1);
            std::size_t parsed = 0;
            const int exponent = std::stoi(suffix, &parsed);
            if (parsed != suffix.size()) {
                throw std::invalid_argument("FloatMP invalid decimal exponent");
            }
            exact *= pow(Rational(10), static_cast<Int>(exponent));
        }
        // Input is rounded once from its exact decimal value, independently
        // of the mutable arithmetic rounding mode.
        FloatMP parsed(exact, FloatMP::ROUND_TO_NEAREST, value.precision());
        value = std::move(parsed);
    } catch (const std::exception&) {
        is.setstate(std::ios::failbit);
        return is;
    }
    if (at_eof) {
        is.clear(is.rdstate() & ~std::ios::failbit);
    }
    return is;
}

Positive<FloatMP> abs(FloatMP const& x) {
    return cast_positive(x.negative() ? neg(x) : x);
}

PositiveUpperBound<FloatMP> mag(FloatMP const& x) {
    return PositiveUpperBound<FloatMP>(abs(x));
}

PositiveLowerBound<FloatMP> mig(FloatMP const& x) {
    return PositiveLowerBound<FloatMP>(abs(x));
}

// The immediate neighbors are computed on the p-bit integer significand;
// the exponent only changes on a carry or at a power-of-two boundary.
FloatMP next(RoundUpward, FloatMP const& x) {
    if (x.is_nan() || (x.is_inf() && !x.negative())) return x;
    if (x.is_inf()) return neg(FloatMP::max(x.precision()));
    if (x.is_zero()) return FloatMP::min(x.precision());
    FloatMP r(x);
    const bool all_ones = floatmp_popcount(r._significand) == r._precision;
    const bool first = floatmp_popcount(r._significand) == 1;
    if (r._negative) {
        if (first) {
            if (r._exponent == std::numeric_limits<FloatMP::Exponent>::min()) {
                return FloatMP::zero(r._precision, true);
            }
            r._significand=sub_big_unsigned(
                shift_left(BigUnsigned(1u),r._precision),BigUnsigned(1u));
            --r._exponent;
        } else {
            r._significand=sub_big_unsigned(r._significand,BigUnsigned(1u));
        }
    } else if (all_ones) {
        if (r._exponent == std::numeric_limits<FloatMP::Exponent>::max()) {
            return FloatMP::inf(r.precision());
        }
        r._significand=shift_left(BigUnsigned(1u),r._precision-1u);
        ++r._exponent;
    } else {
        r._significand=add_big_unsigned(r._significand,BigUnsigned(1u));
    }
    return r;
}

FloatMP next(RoundDownward, FloatMP const& x) {
    if (x.is_nan() || (x.is_inf() && x.negative())) return x;
    if (x.is_inf()) return FloatMP::max(x.precision());
    if (x.is_zero()) return neg(FloatMP::min(x.precision()));
    FloatMP r(x);
    const bool all_ones = floatmp_popcount(r._significand) == r._precision;
    const bool first = floatmp_popcount(r._significand) == 1;
    if (r._negative) {
        if (all_ones) {
            if (r._exponent == std::numeric_limits<FloatMP::Exponent>::max()) {
                return FloatMP::inf(Sign::NEGATIVE, r.precision());
            }
            r._significand=shift_left(BigUnsigned(1u),r._precision-1u);
            ++r._exponent;
        } else {
            r._significand=add_big_unsigned(r._significand,BigUnsigned(1u));
        }
    } else if (first) {
        if (r._exponent == std::numeric_limits<FloatMP::Exponent>::min()) {
            return FloatMP::zero(r._precision);
        }
        r._significand=sub_big_unsigned(
                shift_left(BigUnsigned(1u),r._precision),BigUnsigned(1u));
        --r._exponent;
    } else {
        r._significand=sub_big_unsigned(r._significand,BigUnsigned(1u));
    }
    return r;
}

// Integral rounding works directly on the BigUnsigned significand.
// In particular, signed minimum exponents never trigger huge 2^k allocations.
namespace {


bool has_bit(BigUnsigned const& magnitude,std::size_t index) {
    return (magnitude.limb(index/64u) &
            (std::uint64_t{1}<<(index%64u)))!=0u;
}

FloatMP signed_integral_float(BigUnsigned magnitude,bool negative,
                              MultiplePrecision precision) {
    // floor/ceil return subunit cases before calling this helper. For round,
    // a shift equal to the width retains zero but its leading guard bit
    // increments that zero to one. Every caller therefore supplies |n|>=1.
    return FloatMP(Integer::from_magnitude(std::move(magnitude),negative),
                   FloatMP::ROUND_TO_NEAREST,precision);
}

} // namespace

FloatMP floor(FloatMP const& x) {
    if(!x.is_finite() || x.is_zero() || x._exponent>=0) return x;
    const std::uint64_t shift=std::uint64_t{0}-
        static_cast<std::uint64_t>(x._exponent);
    const BigUnsigned& magnitude=x._significand;
    if(shift>=magnitude.bit_length()) {
        return FloatMP(x._negative ? std::int64_t{-1} : std::int64_t{0},x.precision());
    }
    const std::size_t dropped=static_cast<std::size_t>(shift);
    BigUnsigned integral=shift_right(magnitude,dropped);
    if(x._negative && floatmp_any_bits_below(magnitude,dropped)) {
        integral=add_big_unsigned(integral,BigUnsigned(1u));
    }
    return signed_integral_float(std::move(integral),x._negative,x.precision());
}

FloatMP ceil(FloatMP const& x) {
    if(!x.is_finite() || x.is_zero() || x._exponent>=0) return x;
    const std::uint64_t shift=std::uint64_t{0}-
        static_cast<std::uint64_t>(x._exponent);
    const BigUnsigned& magnitude=x._significand;
    if(shift>=magnitude.bit_length()) {
        return FloatMP(x._negative ? std::int64_t{0} : std::int64_t{1},x.precision());
    }
    const std::size_t dropped=static_cast<std::size_t>(shift);
    BigUnsigned integral=shift_right(magnitude,dropped);
    if(!x._negative && floatmp_any_bits_below(magnitude,dropped)) {
        integral=add_big_unsigned(integral,BigUnsigned(1u));
    }
    return signed_integral_float(std::move(integral),x._negative,x.precision());
}

FloatMP round(FloatMP const& x) {
    if(!x.is_finite() || x.is_zero() || x._exponent>=0) return x;
    const std::uint64_t shift=std::uint64_t{0}-
        static_cast<std::uint64_t>(x._exponent);
    const BigUnsigned& magnitude=x._significand;
    const std::size_t bits=magnitude.bit_length();
    if(shift>bits) {
        return FloatMP::zero(x._precision,x._negative);
    }
    const std::size_t dropped=static_cast<std::size_t>(shift);
    BigUnsigned integral=shift_right(magnitude,dropped);
    // Preserve the historical nearest-integer/ties-away-from-zero contract.
    if(has_bit(magnitude,dropped-1u)) {
        integral=add_big_unsigned(integral,BigUnsigned(1u));
    }
    return signed_integral_float(std::move(integral),x._negative,x.precision());
}

FloatMP abs(FloatMP::RoundingModeType, FloatMP const& x) {
    return x.negative() ? neg(x) : x;
}

FloatMP mag(FloatMP::RoundingModeType mode, FloatMP const& x) {
    return abs(mode, x);
}

FloatMP& iadd(FloatMP::RoundingModeType mode, FloatMP& a, FloatMP const& b) {
    return a = FloatMP::add(a, b, static_cast<unsigned long>(a.precision().bits()), mode);
}

FloatMP& isub(FloatMP::RoundingModeType mode, FloatMP& a, FloatMP const& b) {
    return a = FloatMP::sub(a, b, static_cast<unsigned long>(a.precision().bits()), mode);
}

FloatMP& imul(FloatMP::RoundingModeType mode, FloatMP& a, FloatMP const& b) {
    return a = FloatMP::mul(a, b, static_cast<unsigned long>(a.precision().bits()), mode);
}

FloatMP& idiv(FloatMP::RoundingModeType mode, FloatMP& a, FloatMP const& b) {
    return a = FloatMP::div(a, b, static_cast<unsigned long>(a.precision().bits()), mode);
}

FloatMP add(FloatMP::RoundingModeType mode, FloatMP const& a, FloatDP const& b) {
    return FloatMP(Rational(a) + Rational(b), mode, a.precision());
}

FloatMP add(FloatMP::RoundingModeType mode, FloatDP const& a, FloatMP const& b) {
    return FloatMP(Rational(a) + Rational(b), mode, b.precision());
}

FloatMP sub(FloatMP::RoundingModeType mode, FloatMP const& a, FloatDP const& b) {
    return FloatMP(Rational(a) - Rational(b), mode, a.precision());
}

FloatMP sub(FloatMP::RoundingModeType mode, FloatDP const& a, FloatMP const& b) {
    return FloatMP(Rational(a) - Rational(b), mode, b.precision());
}

FloatMP mul(FloatMP::RoundingModeType mode, FloatMP const& a, FloatDP const& b) {
    return FloatMP(Rational(a) * Rational(b), mode, a.precision());
}

FloatMP mul(FloatMP::RoundingModeType mode, FloatDP const& a, FloatMP const& b) {
    return FloatMP(Rational(a) * Rational(b), mode, b.precision());
}

FloatMP div(FloatMP::RoundingModeType mode, FloatMP const& a, FloatDP const& b) {
    return FloatMP(Rational(a) / Rational(b), mode, a.precision());
}

FloatMP div(FloatMP::RoundingModeType mode, FloatDP const& a, FloatMP const& b) {
    return FloatMP(Rational(a) / Rational(b), mode, b.precision());
}

FloatMP add(FloatMP::RoundingModeType mode, FloatMP const& a, ExactDouble const& b) {
    return FloatMP(Rational(a) + Rational(b), mode, a.precision());
}

FloatMP add(FloatMP::RoundingModeType mode, ExactDouble const& a, FloatMP const& b) {
    return FloatMP(Rational(a) + Rational(b), mode, b.precision());
}

FloatMP sub(FloatMP::RoundingModeType mode, FloatMP const& a, ExactDouble const& b) {
    return FloatMP(Rational(a) - Rational(b), mode, a.precision());
}

FloatMP sub(FloatMP::RoundingModeType mode, ExactDouble const& a, FloatMP const& b) {
    return FloatMP(Rational(a) - Rational(b), mode, b.precision());
}

FloatMP mul(FloatMP::RoundingModeType mode, FloatMP const& a, ExactDouble const& b) {
    return FloatMP(Rational(a) * Rational(b), mode, a.precision());
}

FloatMP mul(FloatMP::RoundingModeType mode, ExactDouble const& a, FloatMP const& b) {
    return FloatMP(Rational(a) * Rational(b), mode, b.precision());
}

FloatMP div(FloatMP::RoundingModeType mode, FloatMP const& a, ExactDouble const& b) {
    return FloatMP(Rational(a) / Rational(b), mode, a.precision());
}

FloatMP div(FloatMP::RoundingModeType mode, ExactDouble const& a, FloatMP const& b) {
    return FloatMP(Rational(a) / Rational(b), mode, b.precision());
}


FloatMP operator*(TwoExp const& t, FloatMP const& x) {
    return shft(x, t.exponent());
}

FloatMP operator*(FloatMP const& x, TwoExp const& t) {
    return shft(x, t.exponent());
}

FloatMP operator/(FloatMP const& x, TwoExp const& t) {
    return shft(x, -t.exponent());
}

FloatMP mul(FloatMP const& x, TwoExp const& t) {
    return shft(x, t.exponent());
}

FloatMP div(FloatMP const& x, TwoExp const& t) {
    return shft(x, -t.exponent());
}

FloatMP& operator*=(FloatMP& x, TwoExp const& t) {
    return x = x * t;
}

FloatMP& operator/=(FloatMP& x, TwoExp const& t) {
    return x = x / t;
}

// Native fused multiply-add: one outward rounding per endpoint.
Bounds<FloatMP> fma(FloatMP const& a, FloatMP const& b, FloatMP const& c) {
    return Bounds<FloatMP>(fma(FloatMP::ROUND_DOWNWARD,a,b,c),
                           fma(FloatMP::ROUND_UPWARD,a,b,c));
}

// Other outward-rounded FloatMP operations are defined by Bounds<F> in
// float_bounds.hpp; do not define competing overloads here.
FloatMP pow(FloatMP::RoundingModeType mode, FloatMP const& x, Int exponent) {
    if (exponent >= 0) {
        return pow(mode, x, static_cast<Nat>(exponent));
    }
    const auto bits = static_cast<unsigned long>(x.precision().bits());
    const bool odd = (exponent % 2) != 0;
    if (x.is_nan()) return FloatMP::nan(bits);
    if (x.is_inf()) return FloatMP::zero(bits, x.negative() && odd);
    if (x.is_zero()) return FloatMP::infinity(bits, x.negative() && odd);
    // The reciprocal of a rounded power need not be correctly rounded:
    // compute the exact rational power, then round once.
    const Rational exact = pow(static_cast<Rational>(x), exponent);
    return FloatMP(exact, mode, x.precision());
}

Comparison cmp(FloatMP const& a, FloatDP const& b) {
    return cmp(a, Rational(b));
}

Comparison cmp(FloatDP const& a, FloatMP const& b) {
    return cmp(Rational(a), b);
}

bool operator==(FloatMP const& a, FloatDP const& b) { return cmp(a,b)==Comparison::EQUAL; }
bool operator!=(FloatMP const& a, FloatDP const& b) { return cmp(a,b)!=Comparison::EQUAL; }
bool operator< (FloatMP const& a, FloatDP const& b) { return cmp(a,b)==Comparison::LESS; }
bool operator> (FloatMP const& a, FloatDP const& b) { return cmp(a,b)==Comparison::GREATER; }
bool operator<=(FloatMP const& a, FloatDP const& b) {
    auto c=cmp(a,b); return c==Comparison::LESS || c==Comparison::EQUAL;
}

bool operator>=(FloatMP const& a, FloatDP const& b) {
    auto c=cmp(a,b); return c==Comparison::GREATER || c==Comparison::EQUAL;
}

bool operator==(FloatDP const& a, FloatMP const& b) { return cmp(a,b)==Comparison::EQUAL; }
bool operator!=(FloatDP const& a, FloatMP const& b) { return cmp(a,b)!=Comparison::EQUAL; }
bool operator< (FloatDP const& a, FloatMP const& b) { return cmp(a,b)==Comparison::LESS; }
bool operator> (FloatDP const& a, FloatMP const& b) { return cmp(a,b)==Comparison::GREATER; }
bool operator<=(FloatDP const& a, FloatMP const& b) {
    auto c=cmp(a,b); return c==Comparison::LESS || c==Comparison::EQUAL;
}

bool operator>=(FloatDP const& a, FloatMP const& b) {
    auto c=cmp(a,b); return c==Comparison::GREATER || c==Comparison::EQUAL;
}

FloatMP shft(FloatMP const& value, Int amount) {
    FloatMP result(value);
    if (result.is_finite() && !result.is_zero()) {
        result._exponent = FloatMP::checked_add(
            result._exponent, static_cast<FloatMP::Exponent>(amount));
    }
    return result;
}

FloatMP::Float(Integer const& n, RoundingModeType mode, PrecisionType pr)
: Float(from_native_ratio(n,Integer(1u),0,
                          static_cast<unsigned long>(pr.bits()),mode)) {
}

namespace {
// Rational declares is_nan/is_inf as friend functions: ADL finds them,
// whereas qualified lookup (Ariadne::is_nan) does not.
bool rational_is_nan(Rational const& q) {
    return is_nan(q);
}

bool rational_is_inf(Rational const& q) {
    return is_inf(q);
}
}

FloatMP FloatMP::from_ratio(Integer const& numerator,
                            Integer const& denominator,
                            RoundingModeType mode,PrecisionType precision) {
    return from_native_ratio(numerator,denominator,0,
                             static_cast<unsigned long>(precision.bits()),mode);
}

FloatMP::Float(Rational const& q, RoundingModeType mode, PrecisionType pr)
: Float(pr) {
    if (rational_is_nan(q)) {
        *this=nan(_precision);
    }
    else if (rational_is_inf(q)) {
        *this=infinity(_precision,q.numerator().negative());
    } else {
        // Rational already owns canonical Integer/BigUnsigned values.
        // Avoid exporting them to GMP merely to import them back.
        *this=from_native_ratio(q.numerator(),q.denominator(),0,_precision,mode);
    }
}

FloatMP::operator Rational() const {
    if (is_nan()) return Rational::nan();
    if (is_inf()) return Rational::inf(negative() ? Sign::NEGATIVE : Sign::POSITIVE);
    BigUnsigned n=_significand;
    if (_exponent>0) { n=shift_left(n,shift_amount(_exponent)); }
    Integer numerator=Integer::from_magnitude(std::move(n),_negative);
    if (_exponent<0) return Rational::from_dyadic(std::move(numerator),
        shift_amount(checked_sub(0,_exponent)));
    return Rational(std::move(numerator));
}

double FloatMP::get_d_general() const {
    // FloatMP stores an exact binary significand and an exponent. Converting
    // through Rational builds/canonicalizes a potentially enormous 2^k
    // denominator and also silently truncates rather than respecting the
    // active rounding mode. Extract the leading bits directly.
    if (is_nan()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (is_inf()) {
        const double infinity = std::numeric_limits<double>::infinity();
        return _negative ? -infinity : infinity;
    }
    if (is_zero()) {
        return _negative ? -0.0 : 0.0;
    }

    constexpr Exponent digits = std::numeric_limits<double>::digits;
    constexpr Exponent largest = std::numeric_limits<double>::max_exponent - 1;
    constexpr Exponent min_normal = std::numeric_limits<double>::min_exponent - 1;
    constexpr Exponent min_subnormal = min_normal - digits + 1;
    const Rounding::Value mode = get_rounding_mode().value();
    const auto signed_value = [this](double value) {
        return _negative ? -value : value;
    };
    const auto rounds_away = [this,mode]() {
        return mode == Rounding::ToNearest ||
               (mode == Rounding::Upward && !_negative) ||
               (mode == Rounding::Downward && _negative);
    };

    // The finite significand is positive and has at most _precision bits,
    // which itself fits Exponent by construction.
    const Exponent top=static_cast<Exponent>(_significand.bit_length()-1u);
    if (_exponent > 0 && top > std::numeric_limits<Exponent>::max() - _exponent) {
        return signed_value(rounds_away() ?
                            std::numeric_limits<double>::infinity() :
                            std::numeric_limits<double>::max());
    }
    const Exponent highest = _exponent + top;
    if (highest > largest) {
        return signed_value(rounds_away() ?
                            std::numeric_limits<double>::infinity() :
                            std::numeric_limits<double>::max());
    }
    if (highest < min_subnormal - 1) {
        const bool directed_away = (mode == Rounding::Upward && !_negative) ||
                                   (mode == Rounding::Downward && _negative);
        return signed_value(directed_away ?
                            std::numeric_limits<double>::denorm_min() : 0.0);
    }

    // Extract the 53 leading bits exactly, truncating any low bits.
    const std::size_t width=_significand.bit_length();
    const BigUnsigned top_bits=width>static_cast<std::size_t>(digits)
        ? shift_right(_significand,width-static_cast<std::size_t>(digits))
        : shift_left(_significand,static_cast<std::size_t>(digits)-width);
    const std::uint64_t leading=top_bits.limb(0u);
    const Exponent subnormal_shift = highest < min_normal
                                    ? min_normal - highest : 0;
    std::uint64_t mantissa = leading >> static_cast<unsigned>(subnormal_shift);

    // The original significand bits below this cut decide every rounding
    // mode, including IEEE halfway ties and the transition to subnormals.
    const Exponent discarded = top - (digits - 1) + subnormal_shift;
    if (discarded > 0) {
        const std::size_t dropped=static_cast<std::size_t>(discarded);
        bool increment = false;
        switch (mode) {
            case Rounding::ToNearest:
                increment=Detail::nearest_even_binary_increment(_significand.limbs(),dropped,
                    (mantissa&1u)!=0u);
                break;
            case Rounding::Downward:
                increment=_negative &&
                    floatmp_any_bits_below(_significand,dropped);
                break;
            case Rounding::Upward:
                increment=!_negative &&
                    floatmp_any_bits_below(_significand,dropped);
                break;
            default: // TowardZero: the fourth value of the closed rounding type.
                break;
        }
        if (increment) {
            ++mantissa;
        }
    }

    if (highest >= min_normal) {
        // A carry may move the leading bit into the next exponent, or
        // overflow the largest finite double after rounding.
        if (mantissa == (std::uint64_t{1} << digits)) {
            mantissa >>= 1u;
            if (highest == largest) {
                // A carry requires increment: nearest rounding, or a directed
                // mode away from zero. All three overflow to signed infinity.
                return signed_value(std::numeric_limits<double>::infinity());
            }
            return signed_value(std::ldexp(
                static_cast<double>(mantissa),
                static_cast<int>(highest - digits + 2)));
        }
        return signed_value(std::ldexp(
            static_cast<double>(mantissa),
            static_cast<int>(highest - digits + 1)));
    }
    // Integer multiples of the least subnormal are exactly representable.
    // mantissa == 2^52 also correctly produces the minimum normal double.
    return signed_value(std::ldexp(
        static_cast<double>(mantissa), static_cast<int>(min_subnormal)));
}

Comparison cmp(FloatMP const& a, FloatMP const& b) {
    if (a.is_nan() || b.is_nan()) {
        return Comparison::INCOMPARABLE;
    }
    if (a.is_inf() || b.is_inf()) {
        if (a.is_inf() && b.is_inf()) {
            return a.negative() == b.negative() ? Comparison::EQUAL :
                   (a.negative() ? Comparison::LESS : Comparison::GREATER);
        }
        return a.is_inf() ? (a.negative() ? Comparison::LESS : Comparison::GREATER) :
               (b.negative() ? Comparison::GREATER : Comparison::LESS);
    }

    // Signed zeros compare equal, independent of their stored sign.
    if (a.is_zero() || b.is_zero()) {
        if (a.is_zero() && b.is_zero()) {
            return Comparison::EQUAL;
        }
        if (a.is_zero()) {
            return b.negative() ? Comparison::GREATER : Comparison::LESS;
        }
        return a.negative() ? Comparison::LESS : Comparison::GREATER;
    }

    if (a._negative != b._negative) {
        return a._negative ? Comparison::LESS : Comparison::GREATER;
    }

    // Nonzero normalized operands of equal precision have the same
    // significand bit-length. Their exponents alone decide magnitude
    // whenever they differ, regardless of the significand values.
    // This avoids GMP size/shift work on the common equal-precision case.
    if (a._precision == b._precision && a._exponent != b._exponent) {
        const bool greater = (a._exponent > b._exponent) != a._negative;
        return greater ? Comparison::GREATER : Comparison::LESS;
    }

    // Compare significand * 2^exponent directly. A different most
    // significant binary position decides the order without allocating
    // any GMP temporaries or constructing Rational values.
    int magnitude = 0;
    if (a._exponent == b._exponent) {
        magnitude = compare(a._significand,b._significand);
    } else {
        const std::size_t a_bits=a._significand.bit_length();
        const std::size_t b_bits=b._significand.bit_length();

        if (a._exponent > b._exponent) {
            // Unsigned subtraction also handles a signed exponent gap
            // larger than INT64_MAX without undefined behavior.
            const std::uint64_t gap =
                static_cast<std::uint64_t>(a._exponent) -
                static_cast<std::uint64_t>(b._exponent);

            if (a_bits >= b_bits || gap > b_bits - a_bits) {
                magnitude = 1;
            } else if (gap < b_bits - a_bits) {
                magnitude = -1;
            } else {
                // The leading positions coincide. Shift the smaller
                // representation right; the discarded bits break ties.
                const auto shift=static_cast<std::size_t>(gap);
                const BigUnsigned truncated=shift_right(b._significand,shift);
                magnitude=compare(a._significand,truncated);
                if(magnitude==0 && floatmp_any_bits_below(b._significand,shift)) {
                    magnitude = -1;
                }
            }
        } else {
            const std::uint64_t gap =
                static_cast<std::uint64_t>(b._exponent) -
                static_cast<std::uint64_t>(a._exponent);

            if (b_bits >= a_bits || gap > a_bits - b_bits) {
                magnitude = -1;
            } else if (gap < a_bits - b_bits) {
                magnitude = 1;
            } else {
                const auto shift=static_cast<std::size_t>(gap);
                const BigUnsigned truncated=shift_right(a._significand,shift);
                magnitude=compare(truncated,b._significand);
                if(magnitude==0 && floatmp_any_bits_below(a._significand,shift)) {
                    magnitude = 1;
                }
            }
        }
    }

    if (a._negative) {
        magnitude = -magnitude;
    }
    return magnitude < 0 ? Comparison::LESS :
           (magnitude > 0 ? Comparison::GREATER : Comparison::EQUAL);
}

// The normalized significand has precisely p bits. After scaling it
// by p-1 or p binary positions, its integer square root also has p
// bits. The integer remainder determines the exact rounding direction.
// This avoids adaptive approximations and any artificial iteration cap.
FloatMP sqrt(FloatMP::RoundingModeType mode, FloatMP const& value) {
    const unsigned long bits = static_cast<unsigned long>(value.precision().bits());
    if (value.is_nan() || (value.negative() && !value.is_zero())) {
        return FloatMP::nan(bits);
    }
    if (value.is_inf() || value.is_zero()) {
        return value;
    }

    const unsigned long offset = bits - 1;
    const bool odd_exponent = value._exponent % 2 != 0;
    const bool odd_magnitude = (offset % 2 != 0) != odd_exponent;
    const unsigned long shift = offset + (odd_magnitude ? 1ul : 0ul);

    // Compute floor((exponent + p - 1)/2) without overflowing int64_t
    // and with mathematical floor rather than truncation toward zero.
    const FloatMP::Exponent half_exponent =
        value._exponent / 2 - ((value._exponent < 0 && odd_exponent) ? 1 : 0);
    FloatMP::Exponent result_exponent = FloatMP::checked_add(
        half_exponent, static_cast<FloatMP::Exponent>(offset / 2));
    if (odd_exponent && offset % 2 != 0) {
        result_exponent = FloatMP::checked_add(result_exponent, 1);
    }
    result_exponent = FloatMP::checked_sub(
        result_exponent, static_cast<FloatMP::Exponent>(offset));

    // Compute the exact integer square root and remainder with native
    // unsigned limbs. GMP is not part of the root/rounding calculation.
    // Only the final owned significand export remains transitional.
    const BigUnsigned& magnitude=value._significand;
    const BigUnsigned scaled=shift_left(magnitude,static_cast<std::size_t>(
        FloatMP::shift_amount(static_cast<FloatMP::Exponent>(shift))));
    auto [root,remainder]=sqrtrem(scaled);

    bool increment = false;
    if (!remainder.is_zero()) {
        switch (mode.value()) {
            case FloatMP::Rounding::ToNearest:
                // scaled=root*root+remainder. The exact nearest criterion
                // is remainder>root; a halfway tie is impossible.
                increment = compare(remainder,root)>0;
                break;
            case FloatMP::Rounding::Upward:
                increment = true;
                break;
            case FloatMP::Rounding::Downward:
            default: // TowardZero: the fourth value of the closed rounding type.
                break;
        }
    }
    if (increment) {
        root=add_big_unsigned(root,BigUnsigned(1u));
    }
    if (root.bit_length()>bits) {
        root=shift_right(root,1u);
        result_exponent = FloatMP::checked_add(result_exponent, 1);
    }

    FloatMP result(bits);
    result._significand=std::move(root);
    result._exponent = result_exponent;
    return result;
}

FloatMP pow(FloatMP::RoundingModeType mode, FloatMP const& value, Nat exponent) {
    const unsigned long bits = static_cast<unsigned long>(value.precision().bits());
    const unsigned long degree = static_cast<unsigned long>(exponent);
    if (degree == 0) {
        return FloatMP(std::int64_t{1}, value.precision());
    }
    if (value.is_nan()) {
        return FloatMP::nan(bits);
    }
    if (value.is_inf()) {
        return FloatMP::infinity(bits, value.negative() && (degree % 2 != 0));
    }
    if (value.is_zero()) {
        return FloatMP::zero(bits, value.negative() && (degree % 2 != 0));
    }

    // The operand is exactly sign * significand * 2^exponent.
    // A nonnegative integer power therefore has an exact integer
    // significand and a binary scale. Round ONCE, after powering:
    // multiplying already-rounded Floats would introduce double rounding.
    if (degree == 1ul) {
        return value;
    }

    // Calculate exponent*degree with checked binary multiplication.
    // Exception: the unnormalized product may exceed Exponent even when
    // normalization cancels it (e.g. very large powers of exact one).
    // Preserve the old general rational path in that rare situation.
    FloatMP::Exponent product_scale = 0;
    FloatMP::Exponent doubled_scale = value._exponent;
    unsigned long remaining_degree = degree;
    bool scale_overflow = false;
    try {
        while (remaining_degree != 0ul) {
            if ((remaining_degree & 1ul) != 0ul) {
                product_scale = FloatMP::checked_add(
                    product_scale, doubled_scale);
            }
            remaining_degree >>= 1u;
            if (remaining_degree != 0ul) {
                doubled_scale = FloatMP::checked_add(
                    doubled_scale, doubled_scale);
            }
        }
    } catch (std::overflow_error const&) {
        scale_overflow = true;
    }

    if (!scale_overflow) {
        FloatMP result(bits);
        result._negative = value._negative && ((degree & 1ul) != 0ul);
        if (bits<=64u && degree<=16ul) {
            // Near either bounding power of two, the scaled exact power
            // is 2^q +/- n*d + R, with 0<R<1/2 for d!=0 under this
            // bound. Its integer part and inexactness are therefore known
            // without constructing the full p*n-bit integer power.
            const auto word=value._significand.limb(0u);
            const auto lower=std::uint64_t{1}<<(bits-1u);
            const auto upper=bits==64u ? std::uint64_t{0} : std::uint64_t{1}<<bits;
            const auto below=word-lower,above=upper-word;
            const bool near_upper=above<below;
            const auto distance=near_upper ? above : below;
            const auto q=near_upper ? bits : bits-1u;
            if (distance==0u || (q>=2u &&
                distance<=(std::uint64_t{1}<<((q-2u)/2u))/degree)) {
                const auto correction=distance*degree;
                auto retained=near_upper ? upper-correction : lower+correction;
                if (distance!=0u) {
                    switch(mode.value()) {
                        case FloatMP::Rounding::ToNearest: break;
                        case FloatMP::Rounding::Downward:
                            if (result._negative) ++retained;
                            break;
                        case FloatMP::Rounding::Upward:
                            if (!result._negative) ++retained;
                            break;
                        default: // TowardZero.
                            break;
                    }
                }
                result._exponent=FloatMP::checked_add(product_scale,
                    static_cast<FloatMP::Exponent>(q*(degree-1ul)));
                result._significand=BigUnsigned(retained);
                return result;
            }
            // Exact short floating-point powers need at most 1024 bits.
            // Keep intermediate products in two disjoint stack buffers;
            // all multiplication still uses the shared limb substrate.
            // Only the final significand owns storage and is rounded.
            constexpr std::size_t capacity=17u;
            std::array<BigUnsigned::Limb,capacity> first,second;
            auto current=value._significand.limbs();
            bool use_first=true;
            const auto multiply=[&](BigUnsigned::LimbView right) {
                auto& buffer=use_first ? first : second;
                const auto count=current.size()+right.size();
                BigUnsigned::mul_into(current,right,
                    BigUnsigned::MutableLimbView(buffer.data(),count));
                auto used=count;
                // Both factors are nonzero and the full product fits the
                // supplied span. Normalization stops at its leading limb.
                while (buffer[used-1u]==0u) --used;
                current=BigUnsigned::LimbView(buffer.data(),used);
                use_first=!use_first;
            };
            for (auto bit=std::bit_floor(degree)>>1u;bit!=0ul;bit>>=1u) {
                multiply(current);
                if ((degree&bit)!=0ul) multiply(value._significand.limbs());
            }
            FloatMP::round_short(result,current,product_scale,mode);
            return result;
        }
        // Exact power on the already-owned unsigned magnitude.
        result._significand=pow(value._significand,
            static_cast<std::uint64_t>(degree));
        FloatMP::round_scaled_integer(result, product_scale, mode);
        return result;
    }

    // Rare normalization-cancellation/scale-overflow fallback.
    // Keep an exact rational result and one final rounding, without
    // roundtrips through mpq_t/mpz_t. Rational::pow is native.
    const Rational exact=static_cast<Rational>(value);
    return FloatMP(pow(exact,exponent),mode,value.precision());
}

template<> String class_name<MultiplePrecision>() { return "MultiplePrecision"; }
template<> String class_name<FloatMP>() { return "FloatMP"; }
template<> String class_name<Rounded<FloatMP>>() { return "Rounded<FloatMP>"; }


// Exact rational intervals and certified transcendental helpers live in
// one implementation namespace, shared by exp/log and trigonometry.
namespace Detail {

class Fraction {
    // Exact, canonical native Rational. Its arithmetic already performs
    // denominator-GCD factoring and cross-cancellation; reimplementing
    // those algorithms here would duplicate the optimized native core.
    Rational _value;

    explicit Fraction(Rational value) : _value(std::move(value)) { }

public:
    Fraction()=default;
    explicit Fraction(long integer) : _value(Integer(integer)) { }
    explicit Fraction(Integer const& integer) : _value(integer) { }

    // Public conversion owns the exact native exponent/shift checks too.
    // Callers dispatch IEEE special values before constructing a Fraction.
    explicit Fraction(FloatMP const& binary) : _value(binary.operator Rational()) { }

    Fraction(BigUnsigned const& magnitude,unsigned long denominator_bits)
        : _value(Rational::from_dyadic(
              Integer::from_magnitude(magnitude,false),
              static_cast<std::size_t>(denominator_bits))) { }

    Fraction(BigUnsigned const& magnitude,unsigned long denominator_bits,
             bool negative)
        : _value(Rational::from_dyadic(
              Integer::from_magnitude(magnitude,negative),
              static_cast<std::size_t>(denominator_bits))) { }

    Fraction(Fraction const&)=default;
    Fraction(Fraction&&)=default;
    Fraction& operator=(Fraction const&)=default;
    Fraction& operator=(Fraction&&)=default;

    Integer const& numerator() const noexcept {
        return _value.numerator_reference();
    }
    BigUnsigned const& denominator() const noexcept {
        return _value.denominator_reference();
    }

    std::optional<std::size_t> dyadic_denominator_bits() const noexcept {
        auto digits=denominator().limbs();
        bool found=false;
        std::size_t exponent=0u;
        for (std::size_t i=0u;i<digits.size();++i) {
            if (digits[i]==0u) continue;
            if (found || !std::has_single_bit(digits[i])) {
                return std::nullopt;
            }
            found=true;
            exponent=64u*i+
                static_cast<std::size_t>(std::countr_zero(digits[i]));
        }
        // Fraction owns a finite Rational with a positive denominator.
        // Reaching the end therefore found exactly one nonzero limb.
        return std::optional<std::size_t>(exponent);
    }
    FloatMP rounded(unsigned long bits,FloatMP::RoundingModeType mode) const {
        return FloatMP::from_ratio(numerator(),
            Integer::from_magnitude(denominator(),false),
            mode,FloatMP::PrecisionType(bits));
    }
    void halve() {
        _value=hlf(_value);
    }
    int sign() const noexcept {
        return numerator().magnitude().is_zero() ? 0 :
            numerator().negative() ? -1 : 1;
    }
    friend int compare(Fraction const& a,Fraction const& b) {
        const Comparison value=cmp(a._value,b._value);
        return value==Comparison::LESS ? -1 :
               value==Comparison::GREATER ? 1 : 0;
    }
    // Exact dyadic arithmetic bypasses the general Rational GCD:
    // align signed numerators for +/- or multiply numerators for *.
    // All non-dyadic operands retain the native Rational fallback.
    static Fraction dyadic_sum(Fraction const& a,Fraction const& b,
                               std::size_t a_bits,std::size_t b_bits,
                               bool subtract_rhs) {
        const std::size_t bits=std::max(a_bits,b_bits);
        const auto scaled=[](Integer const& n,std::size_t shift) {
            return shift==0u ? n : Integer::from_magnitude(
                shift_left(n.magnitude(),shift),n.negative());
        };
        const Integer left=scaled(a.numerator(),bits-a_bits);
        const Integer right=scaled(b.numerator(),bits-b_bits);
        return Fraction(Rational::from_dyadic(
            subtract_rhs ? sub(left,right) : add(left,right),bits));
    }
    friend Fraction operator+(Fraction const& a,Fraction const& b) {
        if (a.sign()==0 || b.sign()==0) return Fraction(a._value+b._value);
        const auto a_bits=a.dyadic_denominator_bits();
        const auto b_bits=b.dyadic_denominator_bits();
        if (a_bits && b_bits) {
            return dyadic_sum(a,b,*a_bits,*b_bits,false);
        }
        return Fraction(a._value+b._value);
    }
    friend Fraction operator-(Fraction const& a,Fraction const& b) {
        if (a.sign()==0 || b.sign()==0) return Fraction(a._value-b._value);
        const auto a_bits=a.dyadic_denominator_bits();
        const auto b_bits=b.dyadic_denominator_bits();
        if (a_bits && b_bits) {
            return dyadic_sum(a,b,*a_bits,*b_bits,true);
        }
        return Fraction(a._value-b._value);
    }
    friend Fraction operator*(Fraction const& a,Fraction const& b) {
        if (a.sign()==0 || b.sign()==0) return Fraction(a._value*b._value);
        const auto a_bits=a.dyadic_denominator_bits();
        const auto b_bits=b.dyadic_denominator_bits();
        if (a_bits && b_bits) {
            // An unrepresentable binary scale is a resource error. Check
            // scalar metadata before multiplying or allocating a denominator.
            const auto bits=Detail::checked_unsigned_add(*a_bits,*b_bits,
                "FloatMP dyadic product scale overflow");
            return Fraction(Rational::from_dyadic(
                mul(a.numerator(),b.numerator()),bits));
        }
        return Fraction(a._value*b._value);
    }
    friend Fraction operator/(Fraction const& a,Fraction const& b) {
        // Interval division excludes zero before reaching this operator.
        // Direct reductions divide by nonzero constants, positive x+1/x+2,
        // positive exponential bounds, or 1+sqrt(1-x*x), all by construction.
        return Fraction(a._value/b._value);
    }
};

struct Interval {
    Fraction lower;
    Fraction upper;
};

Fraction absolute(Fraction const& x) {
    return x.sign() < 0 ? Fraction{} - x : x;
}

Fraction rational_argument(FloatMP const& x) {
    return Fraction(x);
}

// Canonical scratch words and exact comparisons with a power-of-two
// stopping threshold are common to all fixed-point series.
void normalize_fixed_words(std::vector<BigUnsigned::Limb>& digits) {
    while (!digits.empty() && digits.back()==0u) digits.pop_back();
}

bool within_fixed_threshold(std::vector<BigUnsigned::Limb> const& digits,
                            unsigned long guard) {
    const auto word=static_cast<std::size_t>(guard/64ul);
    if (digits.size()<=word) return true;
    if (digits.size()>word+1u) return false;
    const auto leading=BigUnsigned::Limb{1}<<static_cast<unsigned>(guard%64ul);
    if (digits.back()!=leading) return digits.back()<leading;
    return !BigUnsigned::any_bits_below(digits,static_cast<std::size_t>(guard));
}

// Accumulate a fixed-point series with shared exact limb kernels. The
// caller proves sum>=addend for subtraction (alternating decreasing terms).
// Inputs and scratch are disjoint; swapping retains both buffers' capacities.
void accumulate_fixed_terms(std::vector<BigUnsigned::Limb>& sum,
                            BigUnsigned::LimbView addend,bool subtract,
                            std::vector<BigUnsigned::Limb>& scratch) {
    const std::size_t count=std::max(sum.size(),addend.size());
    scratch.resize(count+1u);
    if (subtract) {
        (void)BigUnsigned::sub_into(sum,addend,scratch);
    } else {
        scratch[count]=BigUnsigned::add_into(sum,addend,
            BigUnsigned::MutableLimbView(scratch.data(),count));
    }
    normalize_fixed_words(scratch);
    sum.swap(scratch);
}

// Shared floor/ceil(a*b/2^shift), with disjoint input/output/scratch.
// Zero coefficients, exact products and rounding carries follow the same
// native kernels for exp, log, sine/cosine and atan. All buffers retain capacity.
void scaled_fixed_product(BigUnsigned::LimbView left,BigUnsigned::LimbView right,
                          std::size_t shift,bool round_upward,
                          std::vector<BigUnsigned::Limb>& output,
                          std::vector<BigUnsigned::Limb>& product,
                          std::vector<BigUnsigned::Limb>& increment_scratch) {
    if (left.empty() || right.empty()) {
        output.clear();
        return;
    }
    // Existing 8-byte limb arrays bound the combined count below SIZE_MAX.
    product.resize(left.size()+right.size());
    BigUnsigned::mul_into(left,right,product);
    normalize_fixed_words(product);
    const bool discarded=round_upward && BigUnsigned::any_bits_below(product,shift);
    output.resize(BigUnsigned::shift_right_output_limbs(product.size(),shift));
    output.resize(BigUnsigned::shift_right_into(product,shift,output));
    if (discarded) {
        const BigUnsigned::Limb one=1u;
        accumulate_fixed_terms(output,BigUnsigned::LimbView(&one,1u),false,increment_scratch);
    }
}

// Evaluate exp on |x| <= 1/2 with a single downward-rounded
// fixed-point Taylor recurrence. Let Q=2^B, a=floor(Q*|x|),
// t_0=Q, t_n=floor(t_(n-1)*a/(Q*n)), and T_n=Q*|x|^n/n!.
// As |x| <= 1/2, T_n - t_n < 6 for every n >= 1:
// the recurrence error satisfies e_n < e_(n-1)/(2*n)+2/n+1.
// After summing through n, the omitted exponential tail is at
// most T_n <= t_n+6, so the complete result is enclosed by
// [sum(t_i), sum(t_i)+t_n+6*(n+1)] / Q. This avoids a second
// rounded recurrence without relying on unproved negligible errors.
// For negative x take the reciprocal, reversing endpoints.
Interval exp_series(Fraction const& x, unsigned long work) {
    // Only exp_interval calls this kernel: exp/tanh dispatch zero before
    // exact argument reduction, and all cached centers are nonzero.
    const Fraction magnitude = absolute(x);
    const unsigned long max_terms = work / 2ul + 64ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(max_terms));
    const unsigned long fixed_bits=Detail::checked_unsigned_add(work,guard,
        "FloatMP exponential work precision overflow");
    const auto fixed_shift_size=Detail::checked_unsigned_cast<std::size_t>(fixed_bits,
        "FloatMP exponential fixed-point shift overflow");

    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const auto fixed_shift=fixed_shift_size;
    const std::size_t whole=fixed_shift/64u;
    const unsigned partial=static_cast<unsigned>(fixed_shift%64u);
    const Integer& native_numerator=magnitude.numerator();
    const BigUnsigned& native_denominator=magnitude.denominator();
    const BigUnsigned argument=quot(
        shift_left(native_numerator.magnitude(),fixed_shift),
        native_denominator);

    // The five buffers retain their capacities across all Taylor terms.
    // In contrast to chained owning arithmetic, no temporary BigUnsigned
    // or GMP integer is constructed for multiplication, shifting, division,
    // or accumulation in the recurrence.
    Words term(whole+1u,0u);
    term[whole]=Limb{1}<<partial;
    Words sum_lower=term;
    Words product,shifted,next_term,next_sum;
    const auto below_threshold=[&](Words const& digits) {
        return within_fixed_threshold(digits,guard);
    };

    // t_n<=Q/(2^n*n!)<=Q/2^(2n-1), since n!>=2^(n-1).
    // Thus by n=floor(work/2)+1 the positive threshold 2^guard is
    // reached, strictly before max_terms. The counter cannot wrap.
    for (unsigned long n = 1ul;; ++n) {
        // t_n=floor(floor(t_(n-1)*a/2^B)/n). Reuse the same
        // multiplication, shifting and short-division buffers.
        // t_0=2^B is nonzero; an empty later term is below the positive
        // threshold and returns at the end of its own iteration.
        if (argument.is_zero()) {
            product.clear();
        } else {
            // Existing 8-byte limb arrays bound the combined count below SIZE_MAX.
            product.resize(term.size()+argument.size());
            BigUnsigned::mul_into(term,argument.limbs(),product);
            normalize_fixed_words(product);
        }

        shifted.resize(BigUnsigned::shift_right_output_limbs(product.size(),fixed_shift));
        shifted.resize(BigUnsigned::shift_right_into(product,fixed_shift,shifted));
        const BigUnsigned::LimbView right_shifted=shifted;
        next_term.resize(right_shifted.size());
        (void)BigUnsigned::divrem_limb_into(
            right_shifted,static_cast<Limb>(n),next_term);
        normalize_fixed_words(next_term);
        term.swap(next_term);

        const std::size_t count=std::max(sum_lower.size(),term.size());
        next_sum.resize(count+1u);
        const Limb carry=BigUnsigned::add_into(
            sum_lower,term,BigUnsigned::MutableLimbView(next_sum.data(),count));
        next_sum[count]=carry;
        normalize_fixed_words(next_sum);
        sum_lower.swap(next_sum);

        if (below_threshold(term)) {
            // The same proven enclosure: sum + term + 6*(n+1).
            const BigUnsigned lower(BigUnsigned::LimbView(
                sum_lower.data(),sum_lower.size()));
            const BigUnsigned last_term(BigUnsigned::LimbView(
                term.data(),term.size()));
            const BigUnsigned tail=mul(
                add(BigUnsigned(static_cast<Limb>(n)),BigUnsigned(1u)),
                BigUnsigned(6u));
            const BigUnsigned upper=add(add(lower,last_term),tail);
            const Interval positive{Fraction(lower,fixed_bits),
                                    Fraction(upper,fixed_bits)};
            if (x.sign() > 0) {
                return positive;
            }
            return {Fraction(1L) / positive.upper,
                    Fraction(1L) / positive.lower};
        }
    }
}

Interval exp_interval(Fraction x, unsigned long work) {
    // Preserve argument reduction, Taylor enclosure, and the exact
    // outward-directed repeated-squaring reconstruction.
    const unsigned long reduction_bits = std::min(
        16ul, std::max(6ul, work / 32ul));
    const Fraction reduction_limit =
        Fraction(1L) / Fraction(static_cast<long>(1ul << reduction_bits));

    unsigned halvings = 0;
    while (compare(absolute(x), reduction_limit) > 0) {
        if (++halvings > 2048u) {
            throw std::overflow_error("FloatMP exp argument range");
        }
        x.halve();
    }
    constexpr char reconstruction_error[]="FloatMP exponential reconstruction precision overflow";
    const auto fixed_bits=Detail::checked_unsigned_add(
        Detail::checked_unsigned_add(work,static_cast<unsigned long>(halvings),reconstruction_error),
        64ul,reconstruction_error);
    const unsigned long reduced_work=fixed_bits-56ul;
    Interval interval = exp_series(x, reduced_work);
    if (halvings == 0u) {
        return interval;
    }

    const auto fixed_shift_size=Detail::checked_unsigned_cast<std::size_t>(fixed_bits,
        "FloatMP exponential reconstruction shift overflow");
    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const std::size_t shift=fixed_shift_size;

    // Quantize the two initial certified rational bounds. The
    // lower is rounded down and the upper rounded up at the same
    // dyadic scale. Avoid canonical Rational/GCD round trips.
    const auto to_scaled=[&](Fraction const& bound,bool round_upward) {
        const Integer& numerator=bound.numerator();
        const BigUnsigned& denominator=bound.denominator();
        auto [value,remainder]=divrem(
            shift_left(numerator.magnitude(),shift),
            denominator);
        if (round_upward && !remainder.is_zero()) {
            value=add(value,BigUnsigned(1u));
        }
        const auto digits=value.limbs();
        return Words(digits.begin(),digits.end());
    };
    Words lower=to_scaled(interval.lower,false);
    Words upper=to_scaled(interval.upper,true);
    Words product,next_lower,next_upper;
    Words increment_scratch;
    const auto square_outward=[&](Words const& value,Words& output,bool round_upward) {
        scaled_fixed_product(value,value,shift,round_upward,output,product,increment_scratch);
    };

    // Square only positive endpoint magnitudes, preserving the
    // certified enclosure through monotonicity of squaring.
    for (unsigned i=0u;i<halvings;++i) {
        square_outward(lower,next_lower,false);
        lower.swap(next_lower);
        square_outward(upper,next_upper,true);
        upper.swap(next_upper);
    }
    return {
        Fraction(BigUnsigned(BigUnsigned::LimbView(
                     lower.data(),lower.size())),fixed_bits),
        Fraction(BigUnsigned(BigUnsigned::LimbView(
                     upper.data(),upper.size())),fixed_bits)
    };
}

// For |z| <= 1/3, log((1+z)/(1-z)) = 2 sum z^(2n+1)/(2n+1).
// Instead of growing exact rational denominators at every term, evaluate
// the series in fixed-point arithmetic with outward-rounded GMP integers.
// The first N terms are enclosed term by term. Since |z| <= 1/3,
// the remaining tail is at most (9/4)*3^(-2N-1), strictly below
// 2^(-work) for N = floor(work/3)+4.
Interval log_series(Fraction const& z, unsigned long work) {
    if (z.sign() == 0) {
        return {Fraction{}, Fraction{}};
    }

    const unsigned long terms = work / 3ul + 4ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(terms));
    const unsigned long fixed_bits=Detail::checked_unsigned_add(work,guard,
        "FloatMP logarithm work precision overflow");
    const auto fixed_shift_size=Detail::checked_unsigned_cast<std::size_t>(fixed_bits,
        "FloatMP logarithm fixed-point shift overflow");

    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const std::size_t shift=fixed_shift_size;
    const Fraction magnitude=absolute(z);
    const Integer& num=magnitude.numerator();
    const BigUnsigned& den=magnitude.denominator();
    // log_interval supplies z=(u-2^q)/(u+2^q), with u odd,
    // q>=1 and x!=1. Numerator/denominator are odd and coprime,
    // with denominator>1; log2_interval supplies 1/3. Thus nonzero
    // z is never dyadic, and division by its denominator is inexact.
    const BigUnsigned z_lower=quot(shift_left(num.magnitude(),shift),den);
    const BigUnsigned z_upper=add(z_lower,BigUnsigned(1u));

    // Match GMP's exact floor/ceil endpoints. Only the initial fraction
    // import and final Fraction construction use the legacy GMP carrier.
    // Every fixed-point product, power, integer division and sum below
    // operates on reusable native caller-owned limb buffers.
    Words power_lo(z_lower.limbs().begin(),z_lower.limbs().end());
    Words power_hi(z_upper.limbs().begin(),z_upper.limbs().end());
    Words square_lo,square_hi,sum_lo,sum_hi;
    Words product,next_power,quotient,tail;
    Words increment_scratch;
    const auto increment=[&](Words& digits) {
        const BigUnsigned::Limb one=1u;
        accumulate_fixed_terms(digits,BigUnsigned::LimbView(&one,1u),false,increment_scratch);
    };

    const auto scaled_product=[&](Words const& left,Words const& right,
                                   bool outward,Words& output) {
        scaled_fixed_product(left,right,shift,outward,output,product,increment_scratch);
    };

    Words next_sum;
    const auto accumulate=[&](Words& sum,Words const& addend) {
        accumulate_fixed_terms(sum,addend,false,next_sum);
    };
    const auto within_threshold=[&](Words const& digits) {
        return within_fixed_threshold(digits,guard);
    };

    scaled_product(power_lo,power_lo,false,square_lo);
    scaled_product(power_hi,power_hi,true,square_hi);

    // Q>=2^32 and y_hi/Q<=1/3+1/Q imply square_hi/Q<1/8.
    // Each upper power update is <=power/8+1, so after m updates
    // power<=3*Q/(8*8^m)+8/7. At m=terms=floor(work/3)+4,
    // 3*power<2^guard. Therefore the threshold terminates the loop
    // by terms iterations, without a separate counter-limit fallback.
    for (unsigned long n=0u;;++n) {
        // n < work/3+4 and work+guard fits, with guard>=32; 2*n+1 fits.
        const Limb divisor=static_cast<Limb>(2ul*n+1ul);
        quotient.resize(power_lo.size());
        (void)BigUnsigned::divrem_limb_into(power_lo,divisor,quotient);
        normalize_fixed_words(quotient);
        accumulate(sum_lo,quotient);

        quotient.resize(power_hi.size());
        const Limb upper_remainder=
            BigUnsigned::divrem_limb_into(power_hi,divisor,quotient);
        normalize_fixed_words(quotient);
        if (upper_remainder!=0u) increment(quotient);
        accumulate(sum_hi,quotient);

        scaled_product(power_lo,square_lo,false,next_power);
        power_lo.swap(next_power);
        scaled_product(power_hi,square_hi,true,next_power);
        power_hi.swap(next_power);

        // Nonzero z has a positive upward-rounded initial power and
        // square. Upward rounding preserves positivity at every update.
        const Limb three=3u;
        tail.resize(power_hi.size()+1u);
        BigUnsigned::mul_into(power_hi,BigUnsigned::LimbView(&three,1u),tail);
        normalize_fixed_words(tail);
        if (within_threshold(tail)) {
            break;
        }
    }

    const BigUnsigned lower=shift_left(
        BigUnsigned(BigUnsigned::LimbView(sum_lo.data(),sum_lo.size())),1u);
    const BigUnsigned upper=add(
        shift_left(BigUnsigned(BigUnsigned::LimbView(
            sum_hi.data(),sum_hi.size())),1u),
        BigUnsigned(BigUnsigned::LimbView(tail.data(),tail.size())));
    const Interval positive{Fraction(lower,fixed_bits),
                            Fraction(upper,fixed_bits)};
    return z.sign()<0
        ? Interval{Fraction{}-positive.upper,Fraction{}-positive.lower}
        : positive;
}

Interval log2_interval(unsigned long work) {
    // An enclosure computed at higher working precision remains valid
    // for all lower-precision requests. This cache contains an interval,
    // not a rounded scalar, and is independent of the rounding mode.
    struct Cache {
        unsigned long work = 0;
        std::optional<Interval> value;
    };
    static thread_local Cache cache;
    if (cache.value && cache.work >= work) {
        return *cache.value;
    }

    const Interval result = log_series(Fraction(1L) / Fraction(3L), work);
    cache.value = result;
    cache.work = work;
    return result;
}

Interval log_interval(FloatMP const& binary, unsigned long work) {
    // binary=s*2^e, 2^(p-1)<=s<2^p. Normalize only the p-bit
    // significand: materializing 2^|e| is unnecessary even at native
    // exponent limits. The exact power e+p-1 may exceed Exponent.
    const auto precision=static_cast<unsigned long>(binary.precision().bits());
    Fraction x(binary.significand(),precision-1ul);
    Integer power=add(Integer(binary.exponent()),Integer(precision-1ul));
    if (compare(x,Fraction(3L)/Fraction(2L))>0) {
        x.halve();
        power=add(power,Integer(1));
    }
    // x is now in [3/4,3/2], hence |(x-1)/(x+1)|<=1/5.
    const Fraction z = (x - Fraction(1L)) / (x + Fraction(1L));
    Interval result = log_series(z, work);
    if (!power.magnitude().is_zero()) {
        const Interval ln2 = log2_interval(work);
        const Fraction count(Integer::from_magnitude(power.magnitude(),false));
        if (!power.negative()) {
            result = {result.lower + count * ln2.lower,
                      result.upper + count * ln2.upper};
        } else {
            result = {result.lower - count * ln2.upper,
                      result.upper - count * ln2.lower};
        }
    }
    return result;
}

// Certified half-centered exponential shortcut. The sign-dependent
// constants are cached as rigorous outward DYADIC enclosures, not
// approximations rounded in a caller-selected mode. The general
// exp_interval remains the fallback everywhere outside the window.
// Compare |numerator|*2^binary_shift <= denominator exactly.
// The common case is decided from base-2 bit lengths without creating
// even a native multiprecision temporary. When their bit lengths match,
// compare the actual native limbs; this is the only ambiguous case.
// GMP is touched only to read the not-yet-migrated Fraction carrier.
bool within_dyadic_window(Fraction const& value,unsigned long binary_shift) {
    if (value.sign()==0) return true;
    const BigUnsigned& numerator=value.numerator().magnitude();
    const BigUnsigned& denominator=value.denominator();
    const std::size_t numerator_bits=numerator.bit_length();
    const std::size_t denominator_bits=denominator.bit_length();
    using BitCount=std::common_type_t<std::size_t,unsigned long>;
    const int order=Detail::compare_unsigned_sum(
        static_cast<BitCount>(numerator_bits),static_cast<BitCount>(binary_shift),
        static_cast<BitCount>(denominator_bits));
    if (order!=0) return order<0;
    // Equality proves the shift fits size_t, even on a narrower ABI.
    return compare(
        shift_left(numerator,static_cast<std::size_t>(binary_shift)),
        denominator)<=0;
}

bool exp_inside_half_window(Fraction const& delta, unsigned long work) {
    // certified starts with at least four guard bits and only increases
    // them. Public precision>=2 therefore gives work>=6. Thus |delta|<=2^(-ceil(work/2)) keeps the entire
    // segment around +/-1/2 inside [-3/4,3/4].
    const unsigned long cutoff_bits=work/2ul+work%2ul;
    if constexpr(Detail::unsigned_narrowing<std::size_t,unsigned long>) {
        if(!Detail::unsigned_fits<std::size_t>(cutoff_bits)) return false;
    }
    return within_dyadic_window(delta,cutoff_bits);
}

// Convert a signed certified rational endpoint to an exact outward
// dyadic quotient. Positive floor/ceil use quotient/remainder directly;
// for negative endpoints the direction of magnitude rounding reverses.
// This is a temporary boundary import until Fraction itself is native.
// Import one GMP-backed exact endpoint, but keep the integer quotient
// and its sign native. This avoids importing the result back out of a
// Fraction when a cached centered correction needs integer bounds.
Integer outward_dyadic_integer(Fraction const& exact,
                               unsigned long bits,
                               bool round_toward_positive) {
    const auto shift=Detail::checked_unsigned_cast<std::size_t>(bits,
        "FloatMP dyadic endpoint shift overflow");
    const Integer& numerator=exact.numerator();
    const BigUnsigned& denominator=exact.denominator();
    const bool negative=numerator.negative();
    const bool round_magnitude_up=negative
        ? !round_toward_positive : round_toward_positive;
    BigUnsigned magnitude=Detail::scaled_unsigned_quotient(
        numerator.magnitude(),denominator,shift,round_magnitude_up);
    return Integer::from_magnitude(std::move(magnitude),negative);
}

Fraction outward_dyadic_endpoint(Fraction const& exact,
                                 unsigned long bits,
                                 bool round_toward_positive) {
    const Integer result=outward_dyadic_integer(
        exact,bits,round_toward_positive);
    return Fraction(result.magnitude(),bits,result.negative());
}

Interval exp_outward_dyadic(Interval const& bounds,
                            unsigned long denominator_bits) {
    return {
        outward_dyadic_endpoint(bounds.lower,denominator_bits,false),
        outward_dyadic_endpoint(bounds.upper,denominator_bits,true)
    };
}

struct ExponentialHalfCenters {
    Interval positive;
    Interval negative;
};

ExponentialHalfCenters const& exponential_half_centers(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<ExponentialHalfCenters> enclosure;
    };
    static thread_local Cache cache;
    if (!cache.enclosure || cache.work < work) {
        const auto checked_work=Detail::checked_unsigned_add(work,40ul,
            "FloatMP exponential half precision overflow");
        const Fraction half = Fraction(1L) / Fraction(2L);
        const Interval positive = exp_interval(half, checked_work);
        // Both bounds are strictly positive, so reciprocity reverses
        // their order and yields an exact interval for exp(-1/2).
        const Interval negative{
            Fraction(1L) / positive.upper,
            Fraction(1L) / positive.lower
        };
        ExponentialHalfCenters fresh{
            exp_outward_dyadic(positive, work + 24ul),
            exp_outward_dyadic(negative, work + 24ul)
        };
        cache.enclosure = std::move(fresh);
        cache.work = work;
    }
    return *cache.enclosure;
}

Interval exp_half_centered_or_general(Fraction const& input,
                                     unsigned long work) {
    const Fraction half = Fraction(1L) / Fraction(2L);
    const Fraction center = input.sign() < 0 ? Fraction{} - half : half;
    const Fraction delta = input - center;
    if (!exp_inside_half_window(delta, work)) {
        return exp_interval(input, work);
    }

    // exp(center+delta)=exp(center)*(1+delta)+R.
    // On the segment [-3/4,3/4], |exp''(t)| < e < 3,
    // hence Taylor gives |R| < 3*delta^2/2 <= 2*delta^2.
    // 1+delta >= 3/4 is positive, so endpoint multiplication
    // preserves order. All operations below use exact rationals.
    const auto& cache = exponential_half_centers(work);
    const Interval& base = input.sign() < 0
        ? cache.negative : cache.positive;
    const Fraction multiplier = Fraction(1L) + delta;
    const Fraction remainder = Fraction(2L) * delta * delta;
    return {base.lower * multiplier - remainder,
            base.upper * multiplier + remainder};
}

// Certified log close to 1/2, analogous to the centered exp route.
// For x=1/2+delta with |delta|<=2^-ceil(work/2), work>=4,
// x lies in [1/4,3/4] and |log''(x)|=1/x^2<=16.
// Consequently log(x)=-log(2)+2*delta+R, |R|<=8*delta^2.
// Cache an outward interval of -log(2); the exact rational slope
// costs no extra certified constant. Other arguments use log_interval.
Interval log_half_centered_interval(Fraction const& delta,
                                    unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<Interval> negative_log_two;
    };
    static thread_local Cache cache;
    if (!cache.negative_log_two || cache.work < work) {
        const auto checked_work=Detail::checked_unsigned_add(work,40ul,
            "FloatMP centered log precision overflow");
        const Interval certified_log_two = log2_interval(checked_work);
        const Interval negative_center{
            Fraction{} - certified_log_two.upper,
            Fraction{} - certified_log_two.lower
        };
        // The helper rounds interval endpoints outward by exact GMP
        // floor/ceiling, including for negative values.
        cache.negative_log_two = exp_outward_dyadic(
            negative_center, work + 24ul);
        cache.work = work;
    }

    const Fraction linear_term = Fraction(2L) * delta;
    const Fraction maximum_remainder = Fraction(8L) * delta * delta;
    return {cache.negative_log_two->lower + linear_term - maximum_remainder,
            cache.negative_log_two->upper + linear_term + maximum_remainder};
}

template<class Evaluate>
FloatMP certified(Evaluate evaluate, FloatMP::RoundingModeType mode,
                  FloatMP::PrecisionType precision) {
    const unsigned long bits = static_cast<unsigned long>(precision.bits());
    unsigned long extra = 32;
    for (;;) {
        constexpr char precision_error[]="FloatMP certified precision overflow";
        const auto guarded_work=Detail::checked_unsigned_add(
            Detail::checked_unsigned_add(bits,extra,precision_error),16ul,precision_error);
        const auto work=guarded_work-16ul;
        const Interval interval = evaluate(work);
        const FloatMP lower = interval.lower.rounded(bits,mode);
        const FloatMP upper = interval.upper.rounded(bits,mode);
        if (lower == upper) {
            return lower;
        }
        extra=Detail::checked_unsigned_add(extra,extra,
            "FloatMP certified refinement overflow");
    }
}


// Arithmetic on these fractions is exact. The approximation error is
// represented only by the endpoints of Interval, never by floating point.

Interval point(Fraction const& value) {
    return {value, value};
}

Interval add(Interval const& a, Interval const& b) {
    return {a.lower + b.lower, a.upper + b.upper};
}

Interval subtract(Interval const& a, Interval const& b) {
    return {a.lower - b.upper, a.upper - b.lower};
}

Interval negative(Interval const& value) {
    return {Fraction{} - value.upper, Fraction{} - value.lower};
}

Interval multiply(Interval const& a, Interval const& b) {
    const Fraction candidates[4] = {
        a.lower * b.lower, a.lower * b.upper,
        a.upper * b.lower, a.upper * b.upper
    };
    Fraction lower(candidates[0]);
    Fraction upper(candidates[0]);
    for (unsigned i = 1; i < 4; ++i) {
        if (compare(candidates[i], lower) < 0) {
            lower = candidates[i];
        }
        if (compare(candidates[i], upper) > 0) {
            upper = candidates[i];
        }
    }
    return {std::move(lower), std::move(upper)};
}

// A finite dyadic argument cannot equal a trigonometric pole. A certified
// denominator containing zero requests more precision, rather than a domain error.
struct TrigonometricRefinementRequired { };

Interval divide(Interval const& numerator, Interval const& denominator) {
    if (denominator.lower.sign() <= 0 && denominator.upper.sign() >= 0) {
        throw TrigonometricRefinementRequired{};
    }

    // Division by a negative interval is the negation of division by
    // its positive reflection. Reuse one monotonic endpoint selection
    // for both signs, including numerators that straddle zero.
    if (denominator.upper.sign()<0) {
        return negative(divide(numerator,negative(denominator)));
    }
    if (numerator.lower.sign()>=0) {
        return {numerator.lower/denominator.upper,numerator.upper/denominator.lower};
    }
    if (numerator.upper.sign()<=0) {
        return {numerator.lower/denominator.lower,numerator.upper/denominator.upper};
    }
    return {numerator.lower/denominator.lower,numerator.upper/denominator.lower};
}

// For |x| <= 1, sine and cosine have alternating, decreasing Taylor
// terms. Instead of two directed fixed-point recurrences, calculate one
// downward recurrence and add a proven bound for every discarded bit.
// The cosine needed by tangent is still derived from sine by a directed
// integer square root, so only a single Taylor series is evaluated.
enum class TrigOutputs { sine, cosine, both };

std::pair<Interval, Interval> small_sine_cosine(Fraction const& x,
                                               unsigned long work,
                                               TrigOutputs outputs) {
    // Public sin/cos/tan dispatch zero before exact reduction. The other
    // caller computes the nonzero half-angle center, so x is nonzero.
    const unsigned long max_terms = work / 2ul + 48ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(max_terms));
    const unsigned long fixed_bits=Detail::checked_unsigned_add(work,guard,
        "FloatMP sine/cosine work precision overflow");
    const auto fixed_shift_size=Detail::checked_unsigned_cast<std::size_t>(fixed_bits,
        "FloatMP sine/cosine fixed-point shift overflow");

    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const std::size_t shift=fixed_shift_size;
    const Fraction magnitude=absolute(x);
    const Integer& numerator=magnitude.numerator();
    const BigUnsigned& denominator=magnitude.denominator();
    const BigUnsigned scale=shift_left(BigUnsigned(1u),shift);
    const BigUnsigned argument=quot(
        shift_left(numerator.magnitude(),shift),denominator);

    Words square,term,sum,product,next_term,quotient;
    const auto within_threshold=[&](Words const& digits) {
        return within_fixed_threshold(digits,guard);
    };

    // Before subtracting t_(2k+1), sum is nonnegative term pairs plus
    // t_(2k)>=t_(2k+1). Therefore every exact subtraction is nonnegative.
    Words next_sum;
    const auto accumulate=[&](Words const& addend,bool subtract) {
        accumulate_fixed_terms(sum,addend,subtract,next_sum);
    };

    // floor(a*b/2^B), with disjoint scratch and output limb storage.
    // Reuse product and result capacity across Taylor terms.
    Words increment_scratch;
    const auto scaled_product=[&](Words const& left,Words const& right,
                                   Words& output) {
        scaled_fixed_product(left,right,shift,false,output,product,increment_scratch);
    };
    Words argument_words(argument.limbs().begin(),argument.limbs().end());
    scaled_product(argument_words,argument_words,square);

    // When both functions are requested, derive cosine from the sine
    // enclosure; retain its native endpoints to avoid GMP round-trips.
    BigUnsigned sine_lower_nonnegative,sine_upper_nonnegative;
    const auto series=[&](bool sine) -> Interval {
        const auto start=sine ? argument.limbs() : scale.limbs();
        term.assign(start.begin(),start.end());
        sum.clear();

        // The first omitted term after iteration n is at most
        // Q/(2n+2)!<=Q/2^(2n+1). It therefore reaches 2^guard
        // before max_terms=work/2+48 iterations, without counter wrap.
        for (unsigned long n=0u;;++n) {
            accumulate(term,n%2ul!=0ul);

            // At n=0 both factors fit. Reaching n>0 means the previous
            // checked product, at least (2n-1)*(2n), fit in unsigned long.
            // Since unsigned long has at least 32 bits, this also bounds
            // the current linear factors by 2n+3<=ULONG_MAX.
            const unsigned long twice_n=2ul*n;
            const unsigned long a=sine ? twice_n+2ul : twice_n+1ul;
            const unsigned long b=sine ? twice_n+3ul : twice_n+2ul;
            const Limb divisor=static_cast<Limb>(Detail::checked_unsigned_multiply(
                a,b,"FloatMP sine/cosine factorial overflow"));

            scaled_product(term,square,next_term);
            quotient.resize(next_term.size());
            (void)BigUnsigned::divrem_limb_into(next_term,divisor,quotient);
            normalize_fixed_words(quotient);
            term.swap(quotient);

            if (within_threshold(term)) {
                // Identical certificate: 8*(n+1) accumulated units,
                // followed by the first omitted term plus 8 units.
                const BigUnsigned result(BigUnsigned::LimbView(
                    sum.data(),sum.size()));
                const BigUnsigned error=mul(
                    BigUnsigned(static_cast<Limb>(n+1ul)),BigUnsigned(8u));
                const BigUnsigned tail=add(
                    BigUnsigned(BigUnsigned::LimbView(
                        term.data(),term.size())),BigUnsigned(8u));
                const BigUnsigned lower_error=n%2ul==0ul
                    ? add(error,tail) : error;
                const BigUnsigned upper_error=n%2ul==0ul
                    ? error : add(error,tail);
                const bool lower_negative=compare(result,lower_error)<0;
                const BigUnsigned lower=lower_negative
                    ? sub(lower_error,result) : sub(result,lower_error);
                const BigUnsigned upper=add(result,upper_error);
                if (sine && outputs==TrigOutputs::both) {
                    sine_lower_nonnegative=lower_negative
                        ? BigUnsigned{} : lower;
                    sine_upper_nonnegative=upper;
                }
                return {Fraction(lower,fixed_bits,lower_negative),
                        Fraction(upper,fixed_bits)};
            }
        }
    };

    const bool want_sine=outputs!=TrigOutputs::cosine;
    Interval sine=want_sine ? series(true) : point(Fraction{});
    Interval cosine=point(Fraction{});

    if (outputs==TrigOutputs::cosine) {
        cosine=series(false);
    } else if (outputs==TrigOutputs::both) {
        // The fixed-point certificate gives U-sin(|x|)<3*2^-work:
        // at termination t_next<=2^guard and n<=floor(work/2), so
        // the two accumulated error allowances and tails total less
        // than 3*2^guard. Every caller supplies work>=6. Also
        // sin(|x|)<=sin(1)<=1-1/6+1/120=101/120, and 101/120+3*2^-6<8/9.
        // Thus 0<=L<=sin(|x|)<=U<8/9; clamping U or checking L>U
        // would only duplicate invariants of the enclosure construction.
        const BigUnsigned& lower=sine_lower_nonnegative;
        const BigUnsigned& upper=sine_upper_nonnegative;

        // sqrt(Q^2-s^2) is 2-Lipschitz throughout [L,U], since
        // U<8Q/9. One floor root plus twice the interval width and
        // one rounding unit therefore encloses both cosine endpoints.
        const BigUnsigned squared_scale=mul(scale,scale);
        const BigUnsigned cosine_lower=isqrt(sub(squared_scale,mul(upper,upper)));
        const BigUnsigned width=sub(upper,lower);
        const BigUnsigned cosine_upper=add(cosine_lower,
            add(shift_left(width,1u),BigUnsigned(1u)));
        cosine={Fraction(cosine_lower,fixed_bits),
                Fraction(cosine_upper,fixed_bits)};
    }

    if (want_sine && x.sign()<0) {
        sine=negative(sine);
    }
    return {std::move(sine),std::move(cosine)};
}

// Multiply an interval by a signed exact scalar without evaluating four
// products. This is used in the centered trigonometric enclosure.
Interval multiply_by_fraction(Interval const& input, Fraction const& scalar) {
    if (scalar.sign() >= 0) {
        return {input.lower * scalar, input.upper * scalar};
    }
    return {input.upper * scalar, input.lower * scalar};
}

// Cache certified enclosures of sin(1/2) and cos(1/2), rather than
// rounded scalars. The same (or tighter) intervals are valid for every
// later request of no greater working precision on this thread.
std::pair<Interval, Interval> const& half_trigonometric_intervals(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<std::pair<Interval, Interval>> bounds;
    };
    static thread_local Cache cache;
    if (!cache.bounds || cache.work < work) {
        cache.bounds = small_sine_cosine(
            Fraction(1L) / Fraction(2L), work, TrigOutputs::both);
        cache.work = work;
    }
    return *cache.bounds;
}

// The cached half-angle sine/cosine endpoints are certified rationals.
// Convert each one OUTWARD to a shared dyadic scale once per thread and
// work-precision upgrade, then perform the local Taylor enclosure using
// integer significands rather than generic rational arithmetic.
std::pair<Interval, Interval> centered_trigonometric_dyadic(
    Fraction const& argument, Fraction const& offset,
    unsigned long work, TrigOutputs output) {
    struct Cache {
        unsigned long work = 0;
        unsigned long bits = 0;
        Integer sine_lower, sine_upper, cosine_lower, cosine_upper;
    };
    static thread_local Cache fixed_center;

    if (fixed_center.work < work) {
        const auto checked_work=Detail::checked_unsigned_add(work,32ul,
            "FloatMP centered fixed-point precision overflow");
        const unsigned long target_bits = work + 16ul;
        const auto& [sine, cosine] =
            half_trigonometric_intervals(checked_work);
        // Keep signed native coefficients for the complete cached work
        // precision. Their construction is one-time, on cache upgrade.
        fixed_center.sine_lower =
            outward_dyadic_integer(sine.lower,target_bits,false);
        fixed_center.sine_upper =
            outward_dyadic_integer(sine.upper,target_bits,true);
        fixed_center.cosine_lower =
            outward_dyadic_integer(cosine.lower,target_bits,false);
        fixed_center.cosine_upper =
            outward_dyadic_integer(cosine.upper,target_bits,true);
        fixed_center.bits = target_bits;
        fixed_center.work = work;
    }

    // The caller established a canonical dyadic offset n/2^k,
    // with |offset| <= 2^-ceil(work/2). The carrier is native.
    const std::size_t denominator_shift =
        offset.dyadic_denominator_bits().value();
    const std::size_t shift=denominator_shift;
    const Integer& delta_numerator=offset.numerator();

    // The original Taylor certificate is |R|<=delta^2/2,
    // hence 2^(bits-work-1) integer units. The same bound is used
    // for both the sine and cosine first-order corrections.
    const Integer error=Integer::from_magnitude(
        shift_left(BigUnsigned(1u),
            static_cast<std::size_t>(fixed_center.bits-work-1ul)),false);

    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    Words product,scaled;
    Words increment_scratch;

    // floor or ceil of coefficient*n / 2^k. There are only two
    // products per requested function. Both share the same scratch;
    // no GMP product/division or canonical Rational is constructed.
    const auto multiply_outward=[&](Integer const& coefficient,
                                    bool round_toward_positive) -> Integer {
        const bool negative=coefficient.negative()!=delta_numerator.negative();
        const bool round_magnitude_up=negative ? !round_toward_positive : round_toward_positive;
        scaled_fixed_product(coefficient.magnitude().limbs(),delta_numerator.magnitude().limbs(),
            shift,round_magnitude_up,scaled,product,increment_scratch);
        return Integer::from_magnitude(BigUnsigned(
            BigUnsigned::LimbView(scaled.data(),scaled.size())),negative);
    };

    const auto certified_endpoint=[&](Integer const& exact) {
        return Fraction(exact.magnitude(),fixed_center.bits,exact.negative());
    };
    Interval sine_result=point(Fraction{});
    Interval cosine_result=point(Fraction{});
    if (output!=TrigOutputs::cosine) {
        // sin(1/2+d) = S + C*d + R, |R| <= d^2/2.
        const bool delta_negative=delta_numerator.negative();
        const Integer correction_lower=multiply_outward(
            delta_negative ? fixed_center.cosine_upper
                           : fixed_center.cosine_lower,false);
        const Integer correction_upper=multiply_outward(
            delta_negative ? fixed_center.cosine_lower
                           : fixed_center.cosine_upper,true);
        const Integer lower=fixed_center.sine_lower+correction_lower-error;
        const Integer upper=fixed_center.sine_upper+correction_upper+error;
        sine_result={certified_endpoint(lower),certified_endpoint(upper)};
        if (argument.sign()<0) sine_result=negative(sine_result);
    }
    if (output!=TrigOutputs::sine) {
        // cos(1/2+d) = C - S*d + R, |R| <= d^2/2.
        const Integer correction_lower=multiply_outward(
            delta_numerator.negative() ? fixed_center.sine_upper
                                       : fixed_center.sine_lower,false);
        const Integer correction_upper=multiply_outward(
            delta_numerator.negative() ? fixed_center.sine_lower
                                       : fixed_center.sine_upper,true);
        const Integer lower=fixed_center.cosine_lower-correction_upper-error;
        const Integer upper=fixed_center.cosine_upper-correction_lower+error;
        cosine_result={certified_endpoint(lower),certified_endpoint(upper)};
    }
    return {std::move(sine_result),std::move(cosine_result)};
}

// For finite delta, test |delta| <= 2^(-ceil(work/2)) using exact
// integer arithmetic. Then delta^2/2 <= 2^(-work-1) (work even)
// or even smaller (work odd), so a first-order Taylor enclosure
// has a certified remainder below one bit at the requested work.
bool sufficiently_near_half(Fraction const& delta, unsigned long work) {
    const unsigned long half_work=work/2ul+work%2ul;
    if constexpr(Detail::unsigned_narrowing<std::size_t,unsigned long>) {
        if(!Detail::unsigned_fits<std::size_t>(half_work)) return false;
    }
    return within_dyadic_window(delta,half_work);
}

// Evaluate atan(x) for |x| <= 1/2 as an alternating Taylor series.
// All intermediates are GMP fixed-point integers: rounding each positive
// power and each term outward yields an enclosure of the partial sum.
// The alternating remainder has the sign of the first omitted term and
// magnitude bounded by that term. No GMP rational is allocated per term.
Interval small_atan(Fraction const& x, unsigned long work) {
    if (x.sign() == 0) {
        return point(Fraction{});
    }

    const unsigned long max_terms = work / 2ul + 4ul;
    const unsigned long guard =
        32ul + static_cast<unsigned long>(std::bit_width(max_terms));
    const unsigned long fixed_bits=Detail::checked_unsigned_add(work,guard,
        "FloatMP arctangent work precision overflow");
    const auto fixed_shift_size=Detail::checked_unsigned_cast<std::size_t>(fixed_bits,
        "FloatMP arctangent fixed-point shift overflow");

    using Limb=BigUnsigned::Limb;
    using Words=std::vector<Limb>;
    const std::size_t shift=fixed_shift_size;
    const Fraction magnitude=absolute(x);
    const Integer numerator=magnitude.numerator();
    const BigUnsigned& denominator=magnitude.denominator();
    const BigUnsigned argument=quot(
        shift_left(numerator.magnitude(),shift),denominator);

    Words power(argument.limbs().begin(),argument.limbs().end());
    Words sum,square,term,product,next_power;

    // All Taylor powers are nonnegative. Output and product retain their
    // capacity across terms; no owning BigUnsigned is made in the loop.
    Words increment_scratch;
    const auto scaled_product=[&](Words const& left,Words const& right,
                                   Words& output) {
        scaled_fixed_product(left,right,shift,false,output,product,increment_scratch);
    };

    // Powers decrease for |x|<=1/2 and divisors increase. Before an
    // odd subtraction, sum contains nonnegative pairs plus the preceding
    // term, which is at least the current term. No borrow is possible.
    Words next_sum;
    const auto accumulate=[&](Words const& addend,bool subtract) {
        accumulate_fixed_terms(sum,addend,subtract,next_sum);
    };
    const auto under_limit=[&](Words const& digits) {
        return within_fixed_threshold(digits,guard);
    };

    // Q=2^B, a=floor(Q*|x|), b=floor(a*a/Q), then
    // p_(n+1)=floor(p_n*b/Q). Keep the original 4-unit power
    // and 5-unit per-term error certificates.
    scaled_product(power,power,square);
    // After iteration n, power<=Q/2^(2n+3), since |x|<=1/2.
    // It reaches 2^guard before max_terms=work/2+4 iterations,
    // so neither a counter-limit exit nor counter wrap is possible.
    for (unsigned long n=0u;;++n) {
        // n < work/2+4, so 2*n+1 <= work+7 < work+guard.
        const Limb divisor=static_cast<Limb>(2ul*n+1ul);
        term.resize(power.size());
        (void)BigUnsigned::divrem_limb_into(power,divisor,term);
        normalize_fixed_words(term);
        accumulate(term,n%2ul!=0ul);

        scaled_product(power,square,next_power);
        power.swap(next_power);

        if (under_limit(power)) {
            // The first omitted alternating term is <= power+4.
            // n+1 included terms introduce at most 5*(n+1) units.
            const BigUnsigned result(BigUnsigned::LimbView(
                sum.data(),sum.size()));
            const BigUnsigned error=mul(
                BigUnsigned(static_cast<Limb>(n+1ul)),BigUnsigned(5u));
            const BigUnsigned tail=add(
                BigUnsigned(BigUnsigned::LimbView(
                    power.data(),power.size())),BigUnsigned(4u));
            const BigUnsigned lower_error=n%2ul==0ul
                ? add(error,tail) : error;
            const BigUnsigned upper_error=n%2ul==0ul
                ? error : add(error,tail);

            const bool lower_negative=compare(result,lower_error)<0;
            const BigUnsigned lower=lower_negative
                ? sub(lower_error,result) : sub(result,lower_error);
            const BigUnsigned upper=add(result,upper_error);
            const Interval positive{
                Fraction(lower,fixed_bits,lower_negative),
                Fraction(upper,fixed_bits)
            };
            return x.sign()<0 ? negative(positive) : positive;
        }
    }
}

Interval pi_interval(unsigned long work) {
    // Keep a single certified enclosure per thread, at the greatest
    // requested work precision. A tighter enclosure is always valid
    // for a smaller work request, independently of rounding mode.
    struct Cache {
        unsigned long work = 0;
        std::optional<Interval> value;
    };
    static thread_local Cache cache;
    if (cache.value && cache.work >= work) {
        return *cache.value;
    }

    // Machin: pi = 16 atan(1/5) - 4 atan(1/239).
    const Interval a = small_atan(Fraction(1L) / Fraction(5L), work);
    const Interval b = small_atan(Fraction(1L) / Fraction(239L), work);
    const Interval result = subtract(multiply(point(Fraction(16L)), a),
                                     multiply(point(Fraction(4L)), b));
    cache.value = result;
    cache.work = work;
    return result;
}

// An interval for atan(1/2) computed at greater working precision
// remains certified for lower precision requests. Cache the enclosure,
// never a rounded scalar, so all four rounding modes share the result.
Interval atan_half_interval(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<Interval> value;
    };
    static thread_local Cache cache;
    if (cache.value && cache.work >= work) {
        return *cache.value;
    }
    const Interval result = small_atan(Fraction(1L) / Fraction(2L), work);
    cache.value = result;
    cache.work = work;
    return result;
}

Interval atan_interval(Fraction const& x, unsigned long work) {
    // The centered branch has additional certified-constant and rational
    // reduction costs. At moderate precision the direct Taylor kernel is
    // faster for small arguments (notably asin's argument near 2-sqrt(3)).
    // At high precision the saved Taylor iterations justify that overhead.
    const Fraction direct_limit = work <= 512ul
        ? Fraction(1L) / Fraction(3L) : Fraction(1L) / Fraction(4L);
    if (compare(absolute(x), direct_limit) < 0) {
        return small_atan(x, work);
    }
    if (x.sign() < 0) {
        return negative(atan_interval(Fraction{} - x, work));
    }

    const Fraction half_center_limit = Fraction(5L) / Fraction(7L);
    if (compare(x, half_center_limit) <= 0) {
        // atan(x) = atan(1/2) + atan((2*x-1)/(x+2)).
        // On [1/4,5/7] the second argument lies in [-2/9,3/19].
        // For x near 1/2 it is tiny, so one Taylor iteration often
        // replaces dozens or hundreds. Both intervals are certified.
        const Fraction reduced = (Fraction(2L) * x - Fraction(1L)) /
                                 (x + Fraction(2L));
        return add(atan_half_interval(work), small_atan(reduced, work));
    }

    const Interval quarter_pi = divide(pi_interval(work), point(Fraction(4L)));
    if (compare(x, Fraction(2L)) <= 0) {
        // atan(x) = pi/4 + atan((x-1)/(x+1)), |argument| <= 1/3.
        const Fraction reduced = (x - Fraction(1L)) / (x + Fraction(1L));
        return add(quarter_pi, small_atan(reduced, work));
    }
    // atan(x) = pi/2 - atan(1/x) for x > 2.
    return subtract(multiply(point(Fraction(2L)), quarter_pi),
                    small_atan(Fraction(1L) / x, work));
}

std::pair<Interval, Interval> sine_cosine_intervals(Fraction x,
                                                   unsigned long work,
                                                   TrigOutputs outputs) {
    unsigned halvings = 0;
    while (compare(absolute(x), Fraction(1L)) > 0) {
        if (++halvings > 4096u) {
            throw std::overflow_error("FloatMP trigonometric argument too large");
        }
        x.halve();
    }

    // For an argument very close to 1/2, use a cached certified
    // sin(1/2), cos(1/2) and the Taylor remainder bound
    // |f(1/2+d)-f(1/2)-f'(1/2)*d| <= d*d/2 for f=sin,cos.
    // This is valid for EVERY real argument, not merely binary inputs.
    // If d is not tiny relative to work, use the original Taylor kernel.
    // The *dyadic fixed-point* centered path is cheap even at 53-bit
    // target precision (85-bit initial work), once its rigorous cached
    // center is initialized. Keep the original exact neighborhood check
    // and Taylor fallback. Tangent needs BOTH intervals and retains its
    // separate certified path: the cached path is for standalone
    // sine/cosine only.
    if (halvings == 0u && work >= 64ul &&
        outputs != TrigOutputs::both) {
        const Fraction half = Fraction(1L) / Fraction(2L);
        const Fraction delta = absolute(x) - half;
        if (sufficiently_near_half(delta, work)) {
            // The sole entry is evaluate(sin/cos/tan, rational_argument(x)):
            // FloatMP inputs are dyadic, and exact halving, abs and subtraction
            // of 1/2 preserve that property. No general-rational path is needed.
            return centered_trigonometric_dyadic(x,delta,work,outputs);
        }
    }

    // Double-angle reconstruction needs BOTH values, regardless of
    // which function the caller eventually requests.
    const TrigOutputs kernel_outputs = halvings == 0u ? outputs : TrigOutputs::both;
    // Certification itself decides when the enclosure is tight enough;
    // one bit per doubling compensates argument-reduction amplification.
    const auto kernel_work=Detail::checked_unsigned_add(work,static_cast<unsigned long>(halvings),
        "FloatMP trigonometric reconstruction precision overflow");
    auto results = small_sine_cosine(x,kernel_work,kernel_outputs);
    for (unsigned n = 0; n < halvings; ++n) {
        const Interval sine = multiply(point(Fraction(2L)),
                                       multiply(results.first, results.second));
        const Interval cosine = subtract(multiply(results.second, results.second),
                                          multiply(results.first, results.first));
        results = {sine, cosine};
    }
    return results;
}

Interval square_root_interval(Fraction const& x, unsigned long work) {
    // The two callers supply 3 or 1-x*x with |x|<1. Public inverse
    // trigonometric dispatch handles endpoints and domain errors first.
    // The argument is an exact positive fraction. Its integer quotient
    // and exact square-root remainder are native; the final
    // interval endpoints stay in the native Fraction carrier.
    constexpr char shift_error[]="FloatMP radical precision too large";
    const auto work_shift=Detail::checked_unsigned_cast<std::size_t>(work,shift_error);
    const auto squared_shift=Detail::checked_unsigned_multiply(work_shift,std::size_t{2u},shift_error);
    const Integer& numerator=x.numerator();
    const BigUnsigned& denominator=x.denominator();
    auto [quotient,division_remainder]=divrem(
        shift_left(numerator.magnitude(),squared_shift),
        denominator);
    auto [root,square_remainder]=sqrtrem(quotient);
    const Fraction lower(root,work);
    // Equality to the original rational requires BOTH exact integer
    // division and an exact square. No extra root^2 * denominator
    // multiply or canonical Rational/GCD construction is necessary.
    if (division_remainder.is_zero() && square_remainder.is_zero()) {
        return point(lower);
    }
    return {lower,Fraction(add(root,BigUnsigned(1u)),work)};
}

Interval asin_interval(Fraction const& x, unsigned long work) {
    // asin(x) = 2 atan(x / (1 + sqrt(1-x*x))).
    // Unlike atan(x/sqrt(1-x*x)), the denominator remains >= 1.
    const Interval radical = square_root_interval(
        Fraction(1L) - x * x, work);
    const Interval denominator = add(point(Fraction(1L)), radical);
    const Interval argument = divide(point(x), denominator);

    // atan is monotone and globally 1-Lipschitz:
    // 0 <= atan(b)-atan(a) <= b-a for a<=b, because
    // atan'(t)=1/(1+t*t) is between zero and one.
    // Evaluating one endpoint, plus the exact input interval width,
    // avoids a second certified arctangent evaluation. Choose the
    // endpoint nearer to zero to reduce its Taylor-series work.
    const Fraction width = argument.upper - argument.lower;
    if (compare(absolute(argument.lower), absolute(argument.upper)) <= 0) {
        const Interval left = atan_interval(argument.lower, work);
        return {Fraction(2L) * left.lower,
                Fraction(2L) * (left.upper + width)};
    }
    const Interval right = atan_interval(argument.upper, work);
    return {Fraction(2L) * (right.lower - width),
            Fraction(2L) * right.upper};
}

// For y=|x|, u=exp(-2*y) lies in [0,1] and
// tanh(y)=(1-u)/(1+u), a decreasing function of u. Apply this
// transform to exact outward exponential endpoints before final rounding.
Interval tanh_interval(Fraction const& x,unsigned long work) {
    // Public tanh dispatches zero first; the cached center is exactly 1/2.
    const Fraction one(1L);
    const Interval exponential=exp_interval(Fraction(-2L)*absolute(x),work);
    const Interval positive{
        (one-exponential.upper)/(one+exponential.upper),
        (one-exponential.lower)/(one+exponential.lower)
    };
    return x.sign()<0 ? negative(positive) : positive;
}

// The same first-order tanh enclosure also serves multiple precision.
// Cache certified dyadic center/derivative bounds rather than repeatedly
// evaluating exp(-2*x) for arguments close to +/-1/2.
Interval tanh_centered_or_general(Fraction const& input,unsigned long work) {
    const Fraction delta=absolute(input)-Fraction(1L)/Fraction(2L);
    if (!exp_inside_half_window(delta,work)) return tanh_interval(input,work);
    struct Cache {
        unsigned long work=0;
        Interval center,slope;
    };
    static thread_local Cache cache;
    if (cache.work<work) {
        const auto checked_work=Detail::checked_unsigned_add(work,40ul,
            "FloatMP tanh center precision overflow");
        const auto center=tanh_interval(Fraction(1L)/Fraction(2L),checked_work);
        const Fraction one(1L);
        const Interval slope{one-center.upper*center.upper,
                             one-center.lower*center.lower};
        cache.center=exp_outward_dyadic(center,work+24ul);
        cache.slope=exp_outward_dyadic(slope,work+24ul);
        cache.work=work;
    }
    auto result=add(cache.center,multiply_by_fraction(cache.slope,delta));
    const Fraction error=delta*delta; // |tanh''|<=2, hence |R|<=delta^2
    result.lower=result.lower-error;
    result.upper=result.upper+error;
    return input.sign()<0 ? negative(result) : result;
}

using Function=FloatTranscendental;

// Replace exact division in repeated near-half tangent evaluations with
// dyadic outward endpoints. The initial general division is computed only
// on a cache miss or when the requested work precision increases.
Interval outward_dyadic_interval(Interval const& input, unsigned long bits) {
    return {
        outward_dyadic_endpoint(input.lower,bits,false),
        outward_dyadic_endpoint(input.upper,bits,true)
    };
}

// Certified linearization of atan around 1/2, for signed inputs.
// Cache an OUTWARD interval of atan(1/2), never a rounded scalar.
// For |delta|<=2^-ceil(work/2) and work>=4, y=|x|=1/2+delta
// lies inside [1/4,3/4]. There
// |atan''(y)|=2*y/(1+y*y)^2 <= 3/2, hence |R|<=delta^2.
// The derivative 4/5 is an EXACT rational, not a second cache.
Interval centered_atan_half_interval(Fraction const& input,
                                      Fraction const& delta,
                                      unsigned long work) {
    // This helper is reached only from certified(arctangent, x), where
    // x is a p-bit FloatMP and work >= p+32. Near 1/2, x (and therefore
    // delta) is dyadic with at most p+1 denominator bits. Hence the exact
    // integer correction below always has a nonnegative shift.
    const auto checked_work=Detail::checked_unsigned_add(work,40ul,
        "FloatMP centered atan precision overflow");
    const unsigned long bits = work + 24ul;
    // Preserve the rational computation on a target whose size_t cannot
    // represent the working binary scale. On ordinary 32/64-bit ABIs,
    // size_t is at least as wide as unsigned long.
#if SIZE_MAX < ULONG_MAX
    if (!Detail::unsigned_fits<std::size_t>(bits+2ul)) {
        struct Cache {
            unsigned long work = 0;
            std::optional<Interval> center;
        };
        static thread_local Cache cached_rational;
        if (!cached_rational.center || cached_rational.work < work) {
            cached_rational.center = outward_dyadic_interval(
                atan_half_interval(checked_work), work + 24ul);
            cached_rational.work = work;
        }
        const Fraction linear = Fraction(4L) * delta / Fraction(5L);
        const Fraction remainder = delta * delta;
        const Interval positive{
            cached_rational.center->lower + linear - remainder,
            cached_rational.center->upper + linear + remainder
        };
        return input.sign() < 0 ? negative(positive) : positive;
    }
#endif
    const std::size_t delta_shift = delta.dyadic_denominator_bits().value();

    struct Cache {
        unsigned long work = 0;
        unsigned long bits = 0;
        Integer lower;
        Integer upper;
    };
    static thread_local Cache centered;
    if (centered.work < work) {
        const Interval certified_center=atan_half_interval(checked_work);
        centered.lower=outward_dyadic_integer(
            certified_center.lower,bits,false);
        centered.upper=outward_dyadic_integer(
            certified_center.upper,bits,true);
        centered.bits=bits;
        centered.work=work;
    }

    // delta=n/2^k: (4/5)*delta*2^B=n*2^(B+2-k)/5.
    // Quotient and remainder are computed once on a native magnitude.
    // For negative n, floor/ceil exchange the upward magnitude
    // correction; exact divisibility must not change either bound.
    const Integer& delta_numerator=delta.numerator();
    const std::size_t shift=static_cast<std::size_t>(
        centered.bits+2ul-static_cast<unsigned long>(delta_shift));
    auto [quotient,remainder_by_five]=divrem_limb(
        shift_left(delta_numerator.magnitude(),shift),5u);
    const bool delta_negative=delta_numerator.negative();
    const bool nonexact=remainder_by_five!=0u;
    const Integer linear_lower=Integer::from_magnitude(
        delta_negative && nonexact
            ? add(quotient,BigUnsigned(1u)) : quotient,
        delta_negative);
    const Integer linear_upper=Integer::from_magnitude(
        !delta_negative && nonexact
            ? add(quotient,BigUnsigned(1u)) : quotient,
        delta_negative);

    // |atan''(y)|<=3/2 for y in [1/4,3/4], so the same
    // first-order Taylor remainder is at most 2^-work.
    const Integer error=Integer::from_magnitude(
        shift_left(BigUnsigned(1u),
                   static_cast<std::size_t>(centered.bits-work)),false);
    const Integer lower=centered.lower+linear_lower-error;
    const Integer upper=centered.upper+linear_upper+error;
    const Interval positive{
        Fraction(lower.magnitude(),centered.bits,lower.negative()),
        Fraction(upper.magnitude(),centered.bits,upper.negative())
    };
    return input.sign()<0 ? negative(positive) : positive;
}

struct TangentHalfEnclosures {
    Interval tangent;
    Interval derivative;
};

// Compute an interval for tan(1/2) from the existing certified
// sine/cosine kernel and round its endpoints OUTWARD to dyadics.
// Because tan(1/2)>0, the derivative 1+tan(1/2)^2 is monotone
// in the cached interval. Quantize the derivative outward too.
TangentHalfEnclosures const& tangent_half_enclosures(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<TangentHalfEnclosures> enclosures;
    };
    static thread_local Cache cache;
    if (!cache.enclosures || cache.work < work) {
        const auto checked_work=Detail::checked_unsigned_add(work,48ul,
            "FloatMP tangent cache work precision overflow");
        const unsigned long target_bits = work + 32ul;
        const auto& [sine, cosine] = half_trigonometric_intervals(checked_work);
        const Interval tangent_ratio_bounds = divide(sine, cosine);
        const Interval cached_tangent_bounds = outward_dyadic_interval(tangent_ratio_bounds, target_bits);
        // sin(1/2)>=23/48 and cos(1/2)<=1. The preceding kernel's
        // error is <3*2^-checked_work, checked_work>=54; outward
        // quantization loses <2^-target_bits, target_bits>=38. Its
        // tangent lower endpoint is therefore strictly positive.
        const Interval raw_derivative{
            Fraction(1L) + cached_tangent_bounds.lower * cached_tangent_bounds.lower,
            Fraction(1L) + cached_tangent_bounds.upper * cached_tangent_bounds.upper
        };
        TangentHalfEnclosures fresh_enclosures{
            cached_tangent_bounds,
            outward_dyadic_interval(raw_derivative, target_bits)
        };
        cache.enclosures = std::move(fresh_enclosures);
        cache.work = work;
    }
    return *cache.enclosures;
}

// On |x-1/2| <= 2^-ceil(work/2), Taylor's theorem gives
// tan(x) = T + (1+T^2)*delta + R, |R|<=6*delta^2.
// Indeed x stays in [1/4,3/4]; here |sin(x)|<=3/4 and
// cos(x)>=1-x*x/2 >= 23/32 > 1/2, hence |tan(x)|<=3/2,
// sec(x)^2<=4, |tan''(x)|=|2*tan(x)*sec(x)^2|<=12.
// Interval arithmetic with cached dyadic T and T' therefore
// gives certified bounds for any sign and rounding mode.
Interval centered_tangent_interval(Fraction const& x,
                                   Fraction const& delta,
                                   unsigned long work) {
    const auto& cache = tangent_half_enclosures(work);
    Interval result = add(cache.tangent,
                          multiply_by_fraction(cache.derivative, delta));
    const Fraction remainder = Fraction(6L) * delta * delta;
    result.lower = result.lower - remainder;
    result.upper = result.upper + remainder;
    return x.sign() < 0 ? negative(result) : result;
}

// For y=|x| close to 1/2, the inverse-trigonometric functions admit
// a certified local enclosure: asin(1/2+d)=pi/6 + (2/sqrt(3))*d + R,
// |R|<=2*d*d. This avoids the radical, rational quotient and atan
// sequence of the general inverse-trigonometric implementation.
// Cache OUTWARD enclosures, never rounded center values.
struct InverseHalfEnclosures {
    Interval arcsine_half;
    Interval arccosine_positive_half;
    Interval arccosine_negative_half;
    Interval slope_half;
};

InverseHalfEnclosures const& inverse_half_enclosures(unsigned long work) {
    struct Cache {
        unsigned long work = 0;
        std::optional<InverseHalfEnclosures> bounds;
    };
    static thread_local Cache center_cache;
    if (!center_cache.bounds || center_cache.work < work) {
        const auto checked_work=Detail::checked_unsigned_add(work,40ul,
            "FloatMP inverse-half precision overflow");
        const unsigned long dyadic_bits = work + 24ul;
        const Interval certified_pi = pi_interval(checked_work);
        const Interval asin_half_bounds =
            divide(certified_pi, point(Fraction(6L)));
        const Interval acos_positive_bounds =
            divide(certified_pi, point(Fraction(3L)));
        const Interval acos_negative_bounds =
            multiply(point(Fraction(2L)), acos_positive_bounds);
        const Interval certified_sqrt_three =
            square_root_interval(Fraction(3L), checked_work);
        const Interval slope_bounds =
            divide(point(Fraction(2L)), certified_sqrt_three);
        center_cache.bounds = InverseHalfEnclosures{
            outward_dyadic_interval(asin_half_bounds, dyadic_bits),
            outward_dyadic_interval(acos_positive_bounds, dyadic_bits),
            outward_dyadic_interval(acos_negative_bounds, dyadic_bits),
            outward_dyadic_interval(slope_bounds, dyadic_bits)
        };
        center_cache.work = work;
    }
    return *center_cache.bounds;
}

// When |d|<=2^-ceil(work/2) and work>=4, y=1/2+d is in
// [1/4,3/4]. There asin''(y)=y/(1-y*y)^(3/2) <24/7<4:
// y<=3/4, 1-y*y>=7/16 and sqrt(1-y*y)>1/2.
// Taylor's theorem gives |R| <= 2*d*d, for BOTH asin and
// acos. Use oddness of asin and acos(-y)=pi-acos(y).
Interval centered_inverse_half_interval(Fraction const& input,
                                        Fraction const& offset,
                                        unsigned long work,
                                        bool evaluate_arccosine) {
    const auto& cached = inverse_half_enclosures(work);
    const Interval linear_term_interval =
        multiply_by_fraction(cached.slope_half, offset);

    Interval enclosure;
    if (!evaluate_arccosine) {
        enclosure = add(cached.arcsine_half, linear_term_interval);
    } else if (input.sign() >= 0) {
        enclosure = subtract(cached.arccosine_positive_half,
                             linear_term_interval);
    } else {
        enclosure = add(cached.arccosine_negative_half,
                        linear_term_interval);
    }

    const Fraction maximum_remainder = Fraction(2L) * offset * offset;
    enclosure.lower = enclosure.lower - maximum_remainder;
    enclosure.upper = enclosure.upper + maximum_remainder;
    return (!evaluate_arccosine && input.sign() < 0)
        ? negative(enclosure) : enclosure;
}

// Short centered evaluations use a common fixed binary scale and the same
// rigorous Taylor bounds as the general near-half paths. No Rational is
// constructed per call, and cached data are outward enclosures of constants,
// shared by every rounding mode, never cached rounded function results.
struct FloatMPShortTranscendentals {
    struct Center {
        unsigned long bits=0;
        BigUnsigned lower,upper,slope_lower,slope_upper;
        unsigned error=0;
        bool slope_negative=false;
        bool result_negative=false;
    };

    static void coefficients(Center& target,Interval const& value,
                             Interval const& slope,unsigned long bits) {
        target.lower=outward_dyadic_integer(value.lower,bits,false).magnitude();
        target.upper=outward_dyadic_integer(value.upper,bits,true).magnitude();
        target.slope_lower=outward_dyadic_integer(slope.lower,bits,false).magnitude();
        target.slope_upper=outward_dyadic_integer(slope.upper,bits,true).magnitude();
    }

    template<Function OperationKind>
    static Center make_center(bool negative,unsigned long bits) {
        static_assert(OperationKind==Function::exponential ||
                      OperationKind==Function::logarithm ||
                      OperationKind==Function::hyperbolic_tangent ||
                      OperationKind==Function::sine || OperationKind==Function::cosine ||
                      OperationKind==Function::tangent || OperationKind==Function::arctangent ||
                      OperationKind==Function::arcsine || OperationKind==Function::arccosine);
        (void)negative; // The logarithm center is independent of this argument.
        Center result;
        if constexpr(OperationKind==Function::exponential) {
            const auto& centers=exponential_half_centers(bits);
            const auto& center=negative ? centers.negative : centers.positive;
            coefficients(result,center,center,bits);
            result.error=2u;
            result.slope_negative=negative;
        } else if constexpr(OperationKind==Function::hyperbolic_tangent) {
            const auto center=tanh_interval(Fraction(1L)/Fraction(2L),bits+40ul);
            const Fraction one(1L);
            const Interval slope{one-center.upper*center.upper,
                                 one-center.lower*center.lower};
            coefficients(result,center,slope,bits);
            // |tanh''|=2*|tanh|*(1-tanh^2)<=2 globally.
            result.error=1u;
            result.result_negative=negative;
        } else if constexpr(OperationKind==Function::logarithm) {
            coefficients(result,log2_interval(bits+32ul),point(Fraction(2L)),bits);
            result.error=8u;
            result.slope_negative=true;
            result.result_negative=true;
        } else if constexpr(OperationKind==Function::sine || OperationKind==Function::cosine) {
            const auto& [sine,cosine]=half_trigonometric_intervals(bits+32ul);
            if constexpr(OperationKind==Function::sine) {
                coefficients(result,sine,cosine,bits);
                result.result_negative=negative;
            } else {
                coefficients(result,cosine,sine,bits);
                result.slope_negative=true;
            }
            result.error=1u;
        } else if constexpr(OperationKind==Function::tangent) {
            const auto& center=tangent_half_enclosures(bits);
            coefficients(result,center.tangent,center.derivative,bits);
            result.error=6u;
            result.result_negative=negative;
        } else if constexpr(OperationKind==Function::arctangent) {
            coefficients(result,atan_half_interval(bits+32ul),
                divide(point(Fraction(4L)),point(Fraction(5L))),bits);
            result.error=1u;
            result.result_negative=negative;
        } else {
            const auto& center=inverse_half_enclosures(bits);
            if constexpr(OperationKind==Function::arcsine) {
                coefficients(result,center.arcsine_half,center.slope_half,bits);
                result.result_negative=negative;
            } else {
                const auto& value=negative ? center.arccosine_negative_half
                                          : center.arccosine_positive_half;
                coefficients(result,value,center.slope_half,bits);
                result.slope_negative=!negative;
            }
            result.error=2u;
        }
        result.bits=bits;
        return result;
    }

    // Each statically selected operation owns its two signed caches.
    // Statically selected center retrieval avoids an indirect callback.
    template<Function OperationKind>
    static Center& get_center(bool negative,unsigned long bits) {
        static thread_local std::array<Center,2u> centers;
        auto& center=centers[static_cast<unsigned>(negative)];
        if (center.bits!=bits) center=make_center<OperationKind>(negative,bits);
        return center;
    }

    // Return floor(slope*n/2^shift), and optionally its ceiling. Products
    // use the existing shared native multiplier, including portable builds.
    static std::uint64_t correction(BigUnsigned const& slope,std::uint64_t n,
                                     unsigned shift,bool round_upward) {
        if (n==0u || slope.is_zero()) return 0u;
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
        // The p<=64 centered kernel uses at most two slope limbs, one
        // offset limb and a shift in [1,64]. Use two exact 64x64
        // products and their carry instead of general limb dispatch.
        // A third product limb takes the exact original fallback.
        if (shift>0u && shift<=64u && slope.size()<=2u) {
            __extension__ using Wide=unsigned __int128;
            const Wide first=static_cast<Wide>(slope.limb(0u))*n;
            const Wide second=static_cast<Wide>(
                slope.size()==2u ? slope.limb(1u) : 0u)*n+(first>>64u);
            if ((second>>64u)==0u) {
                const auto low=static_cast<std::uint64_t>(first);
                const auto high=static_cast<std::uint64_t>(second);
                const auto retained=shift==64u ? high
                    : (low>>shift)|(high<<(64u-shift));
                const bool discarded=shift==64u ? low!=0u
                    : (low&((std::uint64_t{1}<<shift)-1u))!=0u;
                return retained+static_cast<std::uint64_t>(round_upward && discarded);
            }
        }
#endif
        // Keep portable exact multiplication and remainder handling
        // for MSVC and all out-of-contract inputs.
        std::array<BigUnsigned::Limb,3u> product;
        BigUnsigned::mul_into(slope.limbs(),BigUnsigned::LimbView(&n,1u),product);
        const auto word=shift/64u;
        const auto offset=shift%64u;
        std::uint64_t retained=product[word]>>offset;
        if (offset!=0u) retained|=product[word+1u]<<(64u-offset);
        if (round_upward && floatmp_any_bits_below(BigUnsigned::LimbView(product),shift)) ++retained;
        return retained;
    }

    // Both short and wide kernels certify exactly the same postcondition:
    // every enclosed value rounds to one finite FloatMP. Use the shared
    // exact comparison instead of duplicating representation-specific tests.
    static std::optional<FloatMP> certified_endpoints(FloatMP lower,FloatMP const& upper) {
        // Successful finite certificates can compare their normalized
        // representation without entering general numeric ordering.
        // Retain general equality for signed zeros and exceptional cases.
        if (lower._kind==FloatMP::Kind::Finite &&
            upper._kind==FloatMP::Kind::Finite &&
            lower._negative==upper._negative &&
            lower._exponent==upper._exponent &&
            (lower._precision<=64u && upper._precision<=64u
                ? lower._significand.limb(0u)==upper._significand.limb(0u)
                : compare(lower._significand,upper._significand)==0)) {
            return lower;
        }
        if (lower==upper) return lower;
        return std::nullopt;
    }

    template<Function OperationKind>
    static std::optional<FloatMP> evaluate(FloatMP::Rounding mode,
                                           FloatMP const& x) {
        const auto precision=x._precision;
        // The supported short window has an exact one-word signed offset
        // from |x|=1/2. Other arguments use the unchanged adaptive path.
        if (precision>64u || x._exponent < -64 || x._exponent >= 0) return std::nullopt;
        const unsigned denominator_bits=static_cast<unsigned>(-x._exponent);
        const auto half=std::uint64_t{1}<<(denominator_bits-1u);
        const auto magnitude=x._significand.limb(0u);
        const bool offset_negative=magnitude<half;
        const auto offset=offset_negative ? half-magnitude : magnitude-half;
        const unsigned long bits=precision+32ul;
        const auto cutoff=static_cast<unsigned>((bits+1ul)/2ul);
        if (offset!=0u && (cutoff>denominator_bits ||
            offset>(std::uint64_t{1}<<(denominator_bits-cutoff)))) return std::nullopt;

        const auto& center=get_center<OperationKind>(x._negative,bits);

        const auto lo=correction(center.slope_lower,offset,denominator_bits,false);
        const auto hi=correction(center.slope_upper,offset,denominator_bits,true);
        const bool subtract=offset_negative!=center.slope_negative;
        // Since |delta|<=2^-ceil(B/2), C*delta^2*2^B<=C.
        // The correction fits <=52 bits at B<=96, hence these signed
        // offsets cannot overflow. Centers stay positive in this window.
        const std::int64_t lower_offset=(subtract ? -static_cast<std::int64_t>(hi)
            : static_cast<std::int64_t>(lo))-static_cast<std::int64_t>(center.error);
        const std::int64_t upper_offset=(subtract ? -static_cast<std::int64_t>(lo)
            : static_cast<std::int64_t>(hi))+static_cast<std::int64_t>(center.error);
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
        // D148: with p in [53,64], cached positive endpoints occupy at
        // most 128 bits. Certify the two rounded scalar words before
        // creating an owned FloatMP; the portable original path remains
        // below for wider or otherwise unsupported magnitudes.
        if (precision>=53u && center.lower.size()>0u &&
            center.upper.size()>0u && center.lower.size()<=2u &&
            center.upper.size()<=2u) {
            __extension__ using Wide=unsigned __int128;
            const auto native=[&](BigUnsigned const& value) -> Wide {
                const Wide high=value.size()==2u ? value.limb(1u) : 0u;
                return (high<<64u)|value.limb(0u);
            };
            Wide low_value=native(center.lower);
            Wide high_value=native(center.upper);
            const auto adjust=[](Wide& value,std::int64_t displacement) {
                const Wide magnitude=displacement<0
                    ? static_cast<std::uint64_t>(-(displacement+1))+1u
                    : static_cast<std::uint64_t>(displacement);
                if (displacement<0) {
                    if (value<magnitude) return false;
                    value-=magnitude;
                } else {
                    if (value>(~Wide{0})-magnitude) return false;
                    value+=magnitude;
                }
                return true;
            };
            if (adjust(low_value,lower_offset) &&
                adjust(high_value,upper_offset)) {
                const auto width=[](Wide value) {
                    const auto high=static_cast<std::uint64_t>(value>>64u);
                    if (high!=0u) return 128u-static_cast<unsigned>(std::countl_zero(high));
                    return 64u-static_cast<unsigned>(
                        std::countl_zero(static_cast<std::uint64_t>(value)));
                };
                const unsigned lower_width=width(low_value);
                const unsigned upper_width=width(high_value);
                // Short-centered constants are strictly nonzero and
                // have >p bits. Keep round_short as the exact fallback
                // if a future center violates this proven local window.
                if (lower_width>precision && upper_width>precision &&
                    lower_width-precision<=64u && upper_width-precision<=64u) {
                    struct Rounded {
                        std::uint64_t significand;
                        FloatMP::Exponent exponent;
                    };
                    const auto round=[&](Wide value,unsigned shift) -> Rounded {
                        std::uint64_t kept=static_cast<std::uint64_t>(value>>shift);
                        const Wide discarded=value&((Wide{1}<<shift)-1u);
                        if (discarded!=0u) {
                            bool increment=false;
                            switch (mode.value()) {
                                case FloatMP::Rounding::ToNearest: {
                                    const Wide halfway=Wide{1}<<(shift-1u);
                                    increment=discarded>halfway ||
                                        (discarded==halfway && (kept&1u)!=0u);
                                    break;
                                }
                                case FloatMP::Rounding::Downward:
                                    increment=center.result_negative;
                                    break;
                                case FloatMP::Rounding::Upward:
                                    increment=!center.result_negative;
                                    break;
                                default: // TowardZero.
                                    break;
                            }
                            if (increment) ++kept;
                        }
                        FloatMP::Exponent exponent=
                            static_cast<FloatMP::Exponent>(shift)-
                            static_cast<FloatMP::Exponent>(bits);
                        if (kept==0u || (precision<64u &&
                            kept==(std::uint64_t{1}<<precision))) {
                            kept=std::uint64_t{1}<<(precision-1u);
                            ++exponent;
                        }
                        return {kept,exponent};
                    };
                    const auto low=round(low_value,lower_width-
                        static_cast<unsigned>(precision));
                    const auto high=round(high_value,upper_width-
                        static_cast<unsigned>(precision));
                    if (low.exponent!=high.exponent ||
                        low.significand!=high.significand) return std::nullopt;
                    FloatMP result(precision);
                    result._negative=center.result_negative;
                    result._exponent=low.exponent;
                    result._significand=BigUnsigned(low.significand);
                    return result;
                }
            }
        }
#endif
        const auto endpoint=[&](BigUnsigned const& constant,std::int64_t displacement) {
            const std::uint64_t absolute=static_cast<std::uint64_t>(
                displacement<0 ? -displacement : displacement);
            const BigUnsigned::LimbView increment(&absolute,1u);
            std::array<BigUnsigned::Limb,2u> output;
            if (displacement<0) BigUnsigned::sub_into(constant.limbs(),increment,output);
            else BigUnsigned::add_into(constant.limbs(),increment,output);
            FloatMP result(precision);
            result._negative=center.result_negative;
            FloatMP::round_short(result,BigUnsigned::LimbView(output.data(),
                output[1u]!=0u ? 2u : 1u),-static_cast<FloatMP::Exponent>(bits),mode);
            return result;
        };
        FloatMP lower=endpoint(center.lower,lower_offset);
        FloatMP upper=endpoint(center.upper,upper_offset);
        // Equal rounded magnitude endpoints certify both positive and
        // negative enclosures. An ambiguous result uses the existing Ziv
        // refinement, including all difficult halfway/directed cases.
        return certified_endpoints(std::move(lower),upper);
    }
    static std::optional<FloatMP> evaluate_tanh_wide(FloatMP::Rounding mode,
                                                     FloatMP const& x) {
        const auto precision=x._precision;
        if (precision<=64u) return std::nullopt;
        // Oversized working scales cannot be serviced by the general
        // certifier either. Validate scalar metadata before any allocation.
        const auto guarded_precision=Detail::checked_unsigned_add(precision,256ul,
            "FloatMP tanh working precision overflow");
        (void)floatmp_exponent_from_bit_length(
            Detail::checked_unsigned_cast<std::size_t>(guarded_precision,
                "FloatMP tanh working precision out of range"));
        const auto p=static_cast<FloatMP::Exponent>(precision);
        if (x._exponent!=-p && x._exponent!=-p-1) return std::nullopt;
        const bool offset_negative=x._exponent==-p-1;
        const auto denominator_bits=static_cast<std::size_t>(-x._exponent);
        const unsigned long bits=precision+32ul;
        BigUnsigned offset;
        const auto digits=x._significand.limbs();
        const auto leading_bit=std::uint64_t{1}<<((precision-1ul)%64ul);
        if (!offset_negative && digits.back()==leading_bit) {
            // Removing the leading half gives the exact offset. Trimming
            // zero high limbs avoids copying an entire high-precision input.
            offset=BigUnsigned(digits.first(digits.size()-1u));
        } else {
            const auto half=shift_left(BigUnsigned(1u),denominator_bits-1u);
            offset=offset_negative ? sub_big_unsigned(half,x._significand)
                                   : sub_big_unsigned(x._significand,half);
        }
        const auto maximum_offset_bit=denominator_bits-(bits+1ul)/2ul;
        const auto width=offset.bit_length();
        if (width>maximum_offset_bit+1u ||
            (width==maximum_offset_bit+1u &&
             floatmp_any_bits_below(offset,maximum_offset_bit))) return std::nullopt;

        static thread_local Center center;
        if (center.bits!=bits) {
            center=make_center<Function::hyperbolic_tangent>(false,bits);
        }
        // One slope product suffices: its truncation costs <1 unit, the
        // slope enclosure costs <=2*|delta|<1, and Taylor's remainder
        // costs <=1. Thus +/-3 units encloses every accumulated error.
        // At the center, W=bits+40. The exp series has width
        // <2^(-W-r-7); r squarings plus outward quantization give
        // <2^(-W-6). The tanh transform has slope<2, and forming
        // 1-tanh^2 costs another factor<1. Thus the scaled slope
        // enclosure has width<2^-45. Outward integer endpoints
        // differ by at most two units, as required by this certificate.
        const auto correction=shift_right(mul_big_unsigned(center.slope_lower,offset),
                                           denominator_bits);
        const auto enlarged=add_big_unsigned(correction,BigUnsigned(3u));
        const bool correction_at_least_three=compare(correction,BigUnsigned(3u))>=0;
        const auto reduced=correction_at_least_three
            ? sub_big_unsigned(correction,BigUnsigned(3u))
            : sub_big_unsigned(BigUnsigned(3u),correction);
        const auto endpoint=[&](BigUnsigned const& constant,bool upper) {
            const bool enlarge=upper!=offset_negative;
            const auto& displacement=enlarge ? enlarged : reduced;
            bool subtract=offset_negative;
            if (!enlarge && !correction_at_least_three) subtract=!subtract;
            FloatMP result(precision);
            result._negative=x._negative;
            result._significand=subtract ? sub_big_unsigned(constant,displacement)
                                        : add_big_unsigned(constant,displacement);
            FloatMP::round_scaled_integer(result,-static_cast<FloatMP::Exponent>(bits),mode);
            return result;
        };
        FloatMP lower=endpoint(center.lower,false);
        FloatMP upper=endpoint(center.upper,true);
        // 8/3<e<11/4 gives 5/11<tanh(1/2)<7/15. Here p>=65,
        // |delta|<=2^-49, and endpoint/rounding errors are <2^-47.
        // Both rounded magnitudes therefore stay strictly in (1/4,1/2)
        // and share exponent -p-1, even when their significands differ.
        return certified_endpoints(std::move(lower),upper);
    }

};


// The callers select their operation statically. The adaptive certification
// loop accepts its evaluator directly, keeping one shared refinement policy.
template<Function OperationKind>
Interval evaluate(Fraction const& x, unsigned long work) {
    (void)x; // Pi evaluators have no operand.
    static_assert(OperationKind!=Function::exponential && OperationKind!=Function::logarithm);
    // certified supplies work=precision+extra, extra>=4.
    if constexpr(OperationKind==Function::hyperbolic_tangent) {
        return tanh_centered_or_general(x,work);
    } else if constexpr(OperationKind==Function::sine) {
        return sine_cosine_intervals(x,work,TrigOutputs::sine).first;
    } else if constexpr(OperationKind==Function::cosine) {
        return sine_cosine_intervals(x,work,TrigOutputs::cosine).second;
    } else if constexpr(OperationKind==Function::tangent) {
        const Fraction delta=absolute(x)-Fraction(1L)/Fraction(2L);
        if (sufficiently_near_half(delta,work)) {
            return centered_tangent_interval(x,delta,work);
        }
        const auto [sine,cosine]=sine_cosine_intervals(x,work,TrigOutputs::both);
        return divide(sine,cosine);
    } else if constexpr(OperationKind==Function::arctangent) {
        const Fraction delta=absolute(x)-Fraction(1L)/Fraction(2L);
        if (sufficiently_near_half(delta,work)) {
            return centered_atan_half_interval(x,delta,work);
        }
        return atan_interval(x,work);
    } else if constexpr(OperationKind==Function::arcsine) {
        const Fraction offset=absolute(x)-Fraction(1L)/Fraction(2L);
        if (sufficiently_near_half(offset,work)) {
            return centered_inverse_half_interval(x,offset,work,false);
        }
        return asin_interval(x,work);
    } else if constexpr(OperationKind==Function::arccosine) {
        const Fraction offset=absolute(x)-Fraction(1L)/Fraction(2L);
        if (sufficiently_near_half(offset,work)) {
            return centered_inverse_half_interval(x,offset,work,true);
        }
        return subtract(divide(pi_interval(work),point(Fraction(2L))),asin_interval(x,work));
    } else if constexpr(OperationKind==Function::pi) {
        return pi_interval(work);
    } else if constexpr(OperationKind==Function::half_pi) {
        return divide(pi_interval(work),point(Fraction(2L)));
    } else {
        static_assert(OperationKind==Function::negative_half_pi);
        return negative(divide(pi_interval(work),point(Fraction(2L))));
    }
}

using IntervalEvaluator=Interval(*)(Fraction const&,unsigned long);

// A result is returned only when all real values in its certified
// enclosure round to exactly the same FloatMP. This includes directed
// rounding and correctly handles the nearest-even halfway boundaries.
FloatMP certified(IntervalEvaluator evaluator, Fraction const& x,
                  FloatMP::RoundingModeType mode, FloatMP::PrecisionType precision,
                  unsigned long extra=32u) {
    const unsigned long bits = static_cast<unsigned long>(precision.bits());
    for (;;) {
        constexpr char precision_error[]="FloatMP trigonometric precision too large";
        const auto guarded_work=Detail::checked_unsigned_add(
            Detail::checked_unsigned_add(bits,extra,precision_error),32ul,precision_error);
        const auto work=guarded_work-32ul;

        Interval bounds;
        try {
            bounds=evaluator(x,work);
        } catch(TrigonometricRefinementRequired const&) {
            extra=Detail::checked_unsigned_add(extra,extra,
                "FloatMP trigonometric pole refinement overflow");
            continue;
        }
        const FloatMP low = bounds.lower.rounded(bits,mode);
        const FloatMP high = bounds.upper.rounded(bits,mode);
        if (low == high) {
            return low;
        }

        extra=Detail::checked_unsigned_add(extra,extra,
            "FloatMP trigonometric refinement overflow");
    }
}

// The six public unary wrappers fix their operation at compile time.
// Pi constants and the separately handled exp/log/tanh APIs never enter here.
template<Function OperationKind>
FloatMP unary(FloatMP::RoundingModeType mode,FloatMP const& x) {
    static_assert(OperationKind==Function::sine || OperationKind==Function::cosine ||
                  OperationKind==Function::tangent || OperationKind==Function::arctangent ||
                  OperationKind==Function::arcsine || OperationKind==Function::arccosine);
    const unsigned long bits=static_cast<unsigned long>(x.precision().bits());
    if(x.is_nan()) return FloatMP::nan(bits);
    if(x.is_inf()) {
        if constexpr(OperationKind==Function::arctangent) {
            return certified(x.negative() ? evaluate<Function::negative_half_pi> : evaluate<Function::half_pi>,
                             Fraction{},mode,x.precision());
        } else {
            return FloatMP::nan(bits);
        }
    }
    if(x.is_zero()) {
        if constexpr(OperationKind==Function::cosine) {
            return FloatMP(std::int64_t{1},x.precision());
        } else if constexpr(OperationKind==Function::arccosine) {
            return certified(evaluate<OperationKind>,Fraction{},mode,x.precision());
        } else {
            return FloatMP::zero(bits,x.negative());
        }
    }
    if(auto short_result=FloatMPShortTranscendentals::evaluate<OperationKind>(mode,x))
        return std::move(*short_result);
    const Fraction argument=rational_argument(x);
    if constexpr(OperationKind==Function::arcsine || OperationKind==Function::arccosine) {
        if(compare(absolute(argument),Fraction(1L))>0) return FloatMP::nan(bits);
        if(compare(argument,Fraction(1L))==0) {
            if constexpr(OperationKind==Function::arccosine) {
                return FloatMP::zero(bits);
            } else {
                return certified(evaluate<Function::half_pi>,Fraction{},mode,x.precision());
            }
        }
        if(compare(argument,Fraction(-1L))==0) {
            if constexpr(OperationKind==Function::arcsine) {
                return certified(evaluate<Function::negative_half_pi>,Fraction{},mode,x.precision());
            } else {
                return certified(evaluate<Function::pi>,Fraction{},mode,x.precision());
            }
        }
    }
    if constexpr(OperationKind==Function::tangent) {
        // Tangent often needs only a few guard bits. Ambiguous bounds,
        // including a denominator spanning a pole, refine automatically.
        return certified(evaluate<OperationKind>,argument,mode,x.precision(),4ul);
    } else {
        return certified(evaluate<OperationKind>,argument,mode,x.precision());
    }
}

} // namespace Detail

FloatMP exp(FloatMP::RoundingModeType mode, FloatMP const& x) {
    using namespace Detail;
    const unsigned long bits = static_cast<unsigned long>(x.precision().bits());
    if (x.is_nan()) {
        return FloatMP::nan(bits);
    }
    if (x.is_inf()) {
        return x.negative() ? FloatMP::zero(bits) : FloatMP::infinity(bits);
    }
    if (x.is_zero()) {
        return FloatMP(std::int64_t{1}, x.precision());
    }
    if (auto short_result=FloatMPShortTranscendentals::evaluate<Function::exponential>(mode,x)) return std::move(*short_result);
    const Fraction input = rational_argument(x);
    return certified([&](unsigned long work) {
        return exp_half_centered_or_general(input, work);
    }, mode, x.precision());
}

FloatMP log(FloatMP::RoundingModeType mode, FloatMP const& x) {
    using namespace Detail;
    const unsigned long bits = static_cast<unsigned long>(x.precision().bits());
    if (x.is_nan() || (x.negative() && !x.is_zero())) {
        return FloatMP::nan(bits);
    }
    if (x.is_zero()) {
        return FloatMP::infinity(bits, true);
    }
    if (x.is_inf()) {
        return FloatMP::infinity(bits);
    }
    if (auto short_result=FloatMPShortTranscendentals::evaluate<Function::logarithm>(mode,x)) return std::move(*short_result);
    if (x.exponent()==-static_cast<FloatMP::Exponent>(bits-1ul) &&
        floatmp_popcount(x.significand())==1u) {
        return FloatMP::zero(bits);
    }
    return certified([&](unsigned long work) {
        // Only [1/4,1) can meet the centered shortcut's window. Its
        // dyadic denominator then has at most p+1 bits, independently
        // of the full native exponent range handled by log_interval.
        const auto half_exponent=-static_cast<FloatMP::Exponent>(bits);
        if (x.exponent()==half_exponent || x.exponent()==half_exponent-1) {
            const Fraction delta=rational_argument(x)-Fraction(1L)/Fraction(2L);
            if (exp_inside_half_window(delta,work)) {
                return log_half_centered_interval(delta,work);
            }
        }
        return log_interval(x,work);
    }, mode, x.precision());
}

FloatMP tanh(FloatMP::RoundingModeType mode,FloatMP const& x) {
    using namespace Detail;
    const auto bits=x._precision;
    if (x.is_nan()) return FloatMP::nan(bits);
    if (x.is_zero()) return x;
    if (x.is_inf()) return FloatMP(x._negative ? -1 : 1,x.precision());
    if (auto result=FloatMPShortTranscendentals::evaluate<Function::hyperbolic_tangent>(mode,x)) return std::move(*result);
    if (auto result=FloatMPShortTranscendentals::evaluate_tanh_wide(mode,x))
        return std::move(*result);

    const bool directed_toward_zero=mode==FloatMP::Rounding::TowardZero ||
        (mode==FloatMP::Rounding::Downward && !x._negative) ||
        (mode==FloatMP::Rounding::Upward && x._negative);
    // e>2 implies 1-tanh(y)<2*exp(-2*y)<2^(1-2*y).
    // At y>=ceil((p+2)/2), this is below half an ulp beneath 1.
    const FloatMP saturation(Integer((static_cast<std::uint64_t>(bits)+3u)/2u),x.precision());
    if (cmp(abs(x),saturation)!=Comparison::LESS) {
        FloatMP result(x._negative ? -1 : 1,x.precision());
        if (directed_toward_zero) {
            result._significand=sub_big_unsigned(shift_left(BigUnsigned(1u),bits),BigUnsigned(1u));
            result._exponent=-static_cast<FloatMP::Exponent>(bits);
        }
        return result;
    }
    const auto highest=FloatMP::checked_add(x._exponent,
        static_cast<FloatMP::Exponent>(bits-1u));
    // 0<y-tanh(y)<y^3. This conservative small-input window
    // places the difference below half the spacing beneath y, including
    // binade boundaries. The exact answer is strictly closer to zero.
    if (highest<=-static_cast<FloatMP::Exponent>(bits/2u+4u)) {
        if (!directed_toward_zero) return x;
        FloatMP result=x;
        if (floatmp_popcount(x._significand)==1u) {
            result._exponent=FloatMP::checked_sub(x._exponent,1);
            result._significand=sub_big_unsigned(shift_left(BigUnsigned(1u),bits),BigUnsigned(1u));
        } else result._significand=sub_big_unsigned(x._significand,BigUnsigned(1u));
        return result;
    }
    return certified(evaluate<Function::hyperbolic_tangent>,rational_argument(x),mode,x.precision());
}

FloatMP sin(FloatMP::RoundingModeType mode, FloatMP const& x) {
    using namespace Detail;
    return unary<Function::sine>(mode,x);
}

FloatMP cos(FloatMP::RoundingModeType mode, FloatMP const& x) {
    using namespace Detail;
    return unary<Function::cosine>(mode,x);
}

FloatMP tan(FloatMP::RoundingModeType mode, FloatMP const& x) {
    using namespace Detail;
    return unary<Function::tangent>(mode,x);
}

FloatMP atan(FloatMP::RoundingModeType mode, FloatMP const& x) {
    using namespace Detail;
    return unary<Function::arctangent>(mode,x);
}

FloatMP asin(FloatMP::RoundingModeType mode, FloatMP const& x) {
    using namespace Detail;
    return unary<Function::arcsine>(mode,x);
}

FloatMP acos(FloatMP::RoundingModeType mode, FloatMP const& x) {
    using namespace Detail;
    return unary<Function::arccosine>(mode,x);
}

FloatMP Float<MultiplePrecision>::pi(RoundingModeType mode, PrecisionType precision) {
    using namespace Detail;
    static_assert(Rounding::ToNearest==0 && Rounding::Downward==1 &&
                  Rounding::Upward==2 && Rounding::TowardZero==3);
    const auto index=static_cast<unsigned>(mode.value());

    // Unlike the certified interval cache, the rounded result is specific
    // to BOTH target precision and rounding mode. Keep four values for
    // the most recently used precision per thread; no global lock, no
    // dependence on a mutable default rounding mode.
    struct Cache {
        unsigned long bits = 0;
        std::optional<FloatMP> rounded[4];
    };
    static thread_local Cache cache;
    const unsigned long bits = static_cast<unsigned long>(precision.bits());
    if (cache.bits != bits) {
        for (auto& entry : cache.rounded) {
            entry.reset();
        }
        cache.bits = bits;
    }
    auto& entry = cache.rounded[index];
    if (!entry) {
        entry.emplace(certified(evaluate<Function::pi>, Fraction{}, mode, precision));
    }
    return *entry;
}

} // namespace Ariadne
