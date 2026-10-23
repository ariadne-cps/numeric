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
#include "floatap.hpp"
#include "sign.hpp"
#include "extended.hpp"

#include "concepts.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <numeric>
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
Integer integer_neg(Integer const& n) { return neg(n); }
Comparison integer_cmp(Integer const& a,Integer const& b) { return cmp(a,b); }
BigUnsigned unsigned_mul(BigUnsigned const& a,BigUnsigned const& b) { return mul(a,b); }
// R3: compare bounded cross-products directly in fixed stack storage.
// The two product spans are independent of both inputs; no temporary
// heap-owned BigUnsigned is needed for magnitudes up to eight limbs.
int compare_cross_products(BigUnsigned const& an,BigUnsigned const& bd,
                           BigUnsigned const& bn,BigUnsigned const& ad) {
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    // R8: four single-limb magnitudes need just two exact 64x64 products.
    // Avoid calling the generic fused multiply kernel and scanning buffers.
    if (an.size()==1u && bd.size()==1u &&
            bn.size()==1u && ad.size()==1u) {
        __extension__ using Wide=unsigned __int128;
        const Wide left=static_cast<Wide>(an.limb(0u))*bd.limb(0u);
        const Wide right=static_cast<Wide>(bn.limb(0u))*ad.limb(0u);
        return left<right ? -1 : left>right ? 1 : 0;
    }
#endif
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    // R17: certify ordering using only the top two limbs of each
    // magnitude. For x=U*B^k+x0 and y=V*B^m+y0, with
    // 0<=x0<B^k and 0<=y0<B^m:
    //   U*V <= x*y/B^(k+m) < (U+1)*(V+1).
    // Compare the *exact* 2x2-limb lower products against the
    // opposing upper product, formed by adding U+V+1 with carries.
    // A certificate is possible only with matching total limb
    // counts (the same power-of-B scale). Ambiguous cases fall back.
    if (an.size()>=2u && bd.size()>=2u &&
            bn.size()>=2u && ad.size()>=2u &&
            an.size()+bd.size()==bn.size()+ad.size() &&
            (an.size()>2u || bd.size()>2u ||
             bn.size()>2u || ad.size()>2u)) {
        using Limb=BigUnsigned::Limb;
        __extension__ using Wide=unsigned __int128;
        std::array<Limb,4u> left{},right{};
        BigUnsigned::mul_into(an.limbs().last(2u),bd.limbs().last(2u),
            BigUnsigned::MutableLimbView(left.data(),left.size()));
        BigUnsigned::mul_into(bn.limbs().last(2u),ad.limbs().last(2u),
            BigUnsigned::MutableLimbView(right.data(),right.size()));
        const auto compare_words=[](std::array<Limb,4u> const& x,
                                    std::array<Limb,4u> const& y) {
            for (std::size_t i=4u;i!=0u;--i) {
                if (x[i-1u]!=y[i-1u]) {
                    return x[i-1u]<y[i-1u] ? -1 : 1;
                }
            }
            return 0;
        };
        // The prefix upper bound is UV+U+V+1. The carry can
        // overflow 256 bits only for both all-ones prefixes;
        // in that case decline to certify from that bound.
        const auto upper_prefix=[](std::array<Limb,4u> product,
                                   BigUnsigned const& x,
                                   BigUnsigned const& y,
                                   bool& fits) {
            Wide carry=1u;
            for (std::size_t i=0u;i<4u;++i) {
                Wide sum=static_cast<Wide>(product[i])+carry;
                if (i<2u) {
                    sum+=x.limb(x.size()-2u+i);
                    sum+=y.limb(y.size()-2u+i);
                }
                product[i]=static_cast<Limb>(sum);
                carry=sum>>64u;
            }
            fits=carry==0u;
            return product;
        };
        bool fits=false;
        const auto upper_left=upper_prefix(left,an,bd,fits);
        if (fits && compare_words(upper_left,right)<=0) {
            return -1;
        }
        const auto upper_right=upper_prefix(right,bn,ad,fits);
        if (fits && compare_words(upper_right,left)<=0) {
            return 1;
        }
    }
#endif
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
// R28: exact comparison against a logically shifted magnitude.
// Compare bit lengths first; read the shifted limbs without allocating.
int compare_shifted_unsigned(BigUnsigned const& a,BigUnsigned const& b,
                             std::size_t shift) {
    using Limb=BigUnsigned::Limb;
    const std::size_t abits=a.bit_length(),bbits=b.bit_length();
    if (shift>std::numeric_limits<std::size_t>::max()-bbits) { return -1; }
    const std::size_t scaled_bits=bbits+shift;
    if (abits!=scaled_bits) { return abits<scaled_bits ? -1 : 1; }
    const auto left=a.limbs(),right=b.limbs();
    const std::size_t whole=shift/64u;
    const unsigned part=static_cast<unsigned>(shift%64u);
    for (std::size_t i=left.size();i!=0u;--i) {
        const std::size_t index=i-1u;
        Limb shifted=0u;
        if (index>=whole) {
            const std::size_t j=index-whole;
            if (j<right.size()) {
                shifted=part==0u ? right[j] : right[j]<<part;
            }
            if (part!=0u && j!=0u && j-1u<right.size()) {
                shifted|=right[j-1u]>>(64u-part);
            }
        }
        if (left[index]!=shifted) { return left[index]<shifted ? -1 : 1; }
    }
    return 0;
}
// R31: compute gcd(p*B^k+p0, q*B^k+q0) without a multi-limb
// Euclidean sequence when the middle radix-B digits are zero.
// This is a general, exact arithmetic property; dense values use
// the original gcd. The helper serves constructor, add/sub and mul.
BigUnsigned gcd_sparse_scalar(BigUnsigned const& sparse,
                              BigUnsigned::Limb gap) {
    using Limb=BigUnsigned::Limb;
    if (gap==0u) { return sparse; }
    if (gap==1u) { return BigUnsigned(Limb{1}); }
    if (gap>65535u) { return gcd(sparse,BigUnsigned(gap)); }

    // All partial products and their sum fit uint64 for gap<=65535.
    const Limb radix=(std::numeric_limits<Limb>::max()%gap+1u)%gap;
    Limb power=1u,base=radix;
    std::size_t exponent=sparse.size()-1u;
    while (exponent!=0u) {
        if ((exponent&1u)!=0u) { power=power*base%gap; }
        exponent>>=1u;
        if (exponent!=0u) { base=base*base%gap; }
    }
    const Limb remainder=((sparse.limb(sparse.size()-1u)%gap)*power%gap+
                           sparse.limb(0u)%gap)%gap;
    return BigUnsigned(std::gcd(gap,remainder));
}

BigUnsigned gcd_sparse_radix(BigUnsigned const& n,BigUnsigned const& d) {
    using Limb=BigUnsigned::Limb;
    if (n.size()!=d.size() || n.size()<2u || n.size()>128u) {
        return gcd(n,d);
    }
    // A dense pair normally fails at its first middle limb.
    if (n.size()>2u && (n.limb(1u)!=0u || d.limb(1u)!=0u)) {
        return gcd(n,d);
    }
    for (std::size_t i=2u;i+1u<n.size();++i) {
        if (n.limb(i)!=0u || d.limb(i)!=0u) { return gcd(n,d); }
    }

    // R30's exact high-coefficient quotient, retained first.
    // If larger=t*smaller+gap in the low limb only, Euclid gives
    // gcd(larger,smaller)=gcd(smaller,gap).
    const bool n_high_first=n.limb(n.size()-1u)>=d.limb(d.size()-1u);
    const BigUnsigned& larger=n_high_first ? n : d;
    const BigUnsigned& smaller=n_high_first ? d : n;
    const Limb high=larger.limb(larger.size()-1u);
    const Limb low_high=smaller.limb(smaller.size()-1u);
    if (high%low_high==0u) {
        const Limb factor=high/low_high;
        const Limb small_low=smaller.limb(0u);
        if (small_low==0u ||
                factor<=std::numeric_limits<Limb>::max()/small_low) {
            const Limb scaled=factor*small_low;
            const Limb low=larger.limb(0u);
            const Limb gap=low>=scaled ? low-scaled : scaled-low;
            return gcd_sparse_scalar(smaller,gap);
        }
    }

    // R31: for general high coefficients p,q, the determinant
    // delta=q*n0-p*d0 equals q*n-p*d, since their B^k terms cancel.
    // Every common divisor of n,d divides delta. Conversely, when
    // gcd(p,delta)=1, any divisor of n,delta divides p*d and hence
    // divides d; gcd(n,delta)=gcd(n,d). Symmetrically use d if
    // gcd(q,delta)=1. Otherwise retain the original full gcd.
    const Limb p=n.limb(n.size()-1u),q=d.limb(d.size()-1u);
    const Limb n0=n.limb(0u),d0=d.limb(0u);
    if ((n0!=0u && q>std::numeric_limits<Limb>::max()/n0) ||
            (d0!=0u && p>std::numeric_limits<Limb>::max()/d0)) {
        return gcd(n,d);
    }
    const Limb qn=q*n0,pd=p*d0;
    const Limb delta=qn>=pd ? qn-pd : pd-qn;
    if (delta!=0u) {
        if (std::gcd(p,delta)==1u) { return gcd_sparse_scalar(n,delta); }
        if (std::gcd(q,delta)==1u) { return gcd_sparse_scalar(d,delta); }
    }
    return gcd(n,d);
}
// R20: fuse bounded signed cross products without owning intermediates.
Integer signed_cross_sum(Integer const& p,BigUnsigned const& x,
                         Integer const& q,BigUnsigned const& y,bool subtract) {
    constexpr std::size_t max_limbs=8u;
    auto const& a=p.magnitude();
    auto const& b=q.magnitude();
    // Preserve the faster scalar path for four one-limb factors,
    // as well as the existing unit-factor and long-operand paths.
    const bool all_scalar=a.size()==1u && x.size()==1u &&
                          b.size()==1u && y.size()==1u;
    if (all_scalar || a.size()>max_limbs || x.size()>max_limbs ||
            b.size()>max_limbs || y.size()>max_limbs ||
            unsigned_is_one(x) || unsigned_is_one(y)) {
        const Integer lhs=multiply_magnitude(p,x);
        const Integer rhs=multiply_magnitude(q,y);
        return subtract ? integer_sub(lhs,rhs) : integer_add(lhs,rhs);
    }
    using Limb=BigUnsigned::Limb;
    std::array<Limb,17u> out;
    const std::size_t width=std::max(a.size()+x.size(),b.size()+y.size());
    const auto output=BigUnsigned::MutableLimbView(out.data(),width+1u);
    BigUnsigned::mul_into(a.limbs(),x.limbs(),output);
    const bool opposite=p.negative()!=(q.negative()!=subtract);
    const Limb borrow=BigUnsigned::accumulate_product_into(b.limbs(),y.limbs(),output,opposite);
    if (opposite && borrow!=0u) {
        // The guard limb makes a negative difference an unambiguous
        // two's-complement value. Restore its unsigned magnitude.
        Limb carry=1u;
        for (auto& word:output) {
            const Limb inverted=~word;
            word=inverted+carry;
            carry=static_cast<Limb>(word<inverted);
        }
    }
    return Integer::from_magnitude(BigUnsigned(output),p.negative()!=(opposite && borrow!=0u));
}
BigUnsigned exact_quotient(BigUnsigned const& n,BigUnsigned const& d) {
    if (unsigned_is_one(d)) { return n; }
    return n.size()<=1u ? quot(n,d) : quot_exact(n,d);
}
// R19a: add_sub still needs the signed overload for the final
// numerator reduction when the common factor h is greater than one.
Integer exact_quotient(Integer const& n,BigUnsigned const& d) {
    return Integer::from_magnitude(exact_quotient(n.magnitude(),d),n.negative());
}
// R19: construct a product from exact cross-cancelled magnitudes,
// borrowing each original operand whenever its divisor is one.
// quot(x,1) would materialize a full BigUnsigned copy unnecessarily;
// make the decision independently for each factor.
BigUnsigned multiply_cancelled(BigUnsigned const& lhs,BigUnsigned const& lhs_divisor,
                               BigUnsigned const& rhs,BigUnsigned const& rhs_divisor) {
    const bool left_unit=unsigned_is_one(lhs_divisor);
    const bool right_unit=unsigned_is_one(rhs_divisor);
    if (left_unit && right_unit) { return unsigned_mul(lhs,rhs); }
    constexpr std::size_t capacity=8u;
    if (lhs.size()<=capacity && rhs.size()<=capacity &&
            (lhs.size()>1u || rhs.size()>1u) &&
            lhs_divisor.size()==1u && rhs_divisor.size()==1u) {
        using Limb=BigUnsigned::Limb;
        std::array<Limb,capacity> left,right;
        std::array<Limb,2u*capacity> product;
        const auto cancelled=[](BigUnsigned const& n,BigUnsigned const& d,
                                bool unit,BigUnsigned::MutableLimbView output) {
            if (unit) { return n.limbs(); }
            BigUnsigned::divexact_limb_into(n.limbs(),d.limb(0u),output);
            std::size_t count=n.size();
            while (count!=0u && output[count-1u]==0u) { --count; }
            return BigUnsigned::LimbView(output.data(),count);
        };
        const auto a=cancelled(lhs,lhs_divisor,left_unit,left);
        const auto b=cancelled(rhs,rhs_divisor,right_unit,right);
        const std::size_t count=a.size()+b.size();
        BigUnsigned::mul_into(a,b,BigUnsigned::MutableLimbView(product.data(),count));
        return BigUnsigned(BigUnsigned::LimbView(product.data(),count));
    }
    if (unsigned_is_one(lhs_divisor)) {
        return unsigned_mul(lhs,exact_quotient(rhs,rhs_divisor));
    }
    if (unsigned_is_one(rhs_divisor)) { return unsigned_mul(exact_quotient(lhs,lhs_divisor),rhs); }
    return unsigned_mul(exact_quotient(lhs,lhs_divisor),exact_quotient(rhs,rhs_divisor));
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

    // Both owning components are complete before this call. Add/sub
    // has proved the nonzero result canonical, even when r aliases an
    // input. Avoid constructing and swapping a second Rational.
    static Void set_sum(Rational& r,Integer&& n,BigUnsigned&& d) noexcept {
        r._numerator=std::move(n);
        r._denominator=std::move(d);
        r._kind=Rational::Kind::FINITE;
    }

    static Void add_sub(Rational& r,Rational const& a,Rational const& b,bool subtract) {
        if (a._numerator.magnitude().is_zero()) { r=subtract ? reduced(integer_neg(b._numerator),b._denominator) : b; return; }
        if (b._numerator.magnitude().is_zero()) { r=a; return; }
        const BigUnsigned g=gcd_sparse_radix(a._denominator,b._denominator);
        if (unsigned_is_one(g)) {
            // R16: coprime canonical denominators cannot introduce a
            // common factor in the cross-sum numerator. Their product
            // is already the reduced denominator. Borrow both inputs
            // instead of copying them through exact_quotient(x,1),
            // and move the owning result numerator into the rational.
            Integer n=signed_cross_sum(a._numerator,b._denominator,
                                       b._numerator,a._denominator,subtract);
            if (n.magnitude().is_zero()) { r=Rational(); return; }
            set_sum(r,std::move(n),unsigned_mul(a._denominator,b._denominator));
            return;
        }
        const BigUnsigned ad=exact_quotient(a._denominator,g);
        const BigUnsigned bd=exact_quotient(b._denominator,g);
        Integer n=signed_cross_sum(a._numerator,bd,b._numerator,ad,subtract);
        if (n.magnitude().is_zero()) { r=Rational(); return; }
        // Only a factor of the original denominator gcd can remain.
        const BigUnsigned h=gcd(n.magnitude(),g);
        if (unsigned_is_one(h)) {
            // The numerator is already canonical; reuse its ownership
            // and the unreduced input denominator without copies.
            set_sum(r,std::move(n),unsigned_mul(ad,b._denominator));
            return;
        }
        set_sum(r,exact_quotient(n,h),unsigned_mul(ad,exact_quotient(b._denominator,h)));
    }
    static Void add(Rational& r,Rational const& a,Rational const& b) { add_sub(r,a,b,false); }
    static Void sub(Rational& r,Rational const& a,Rational const& b) { add_sub(r,a,b,true); }
    static Void mul(Rational& r,Rational const& a,Rational const& b) {
        if (&a==&b) { sqr(r,a); return; }
        if (a._numerator.magnitude().is_zero() || b._numerator.magnitude().is_zero()) {
            r=Rational(); return;
        }
        const BigUnsigned g=gcd_sparse_radix(a._numerator.magnitude(),b._denominator);
        const BigUnsigned h=gcd_sparse_radix(b._numerator.magnitude(),a._denominator);
        BigUnsigned numerator=multiply_cancelled(
            a._numerator.magnitude(),g,b._numerator.magnitude(),h);
        BigUnsigned denominator=multiply_cancelled(
            a._denominator,h,b._denominator,g);
        // Like add/sub, commit both fully evaluated owned components
        // directly: the intermediate Rational and its swap are unnecessary.
        set_sum(r,Integer::from_magnitude(
            std::move(numerator),a._numerator.negative()!=b._numerator.negative()),
            std::move(denominator));
    }
    static Void div(Rational& r,Rational const& a,Rational const& b) {
        if (b._numerator.magnitude().is_zero()) {
            r=Rational::inf(integer_sign(a._numerator)); return;
        }
        if (a._numerator.magnitude().is_zero()) { r=Rational(); return; }
        // R32: the same exact sparse-radix GCD used by mul also
        // applies to div's cross-cancellation pair. The borrowed
        // arguments remain live until both results are computed.
        const BigUnsigned g=gcd_sparse_radix(
            a._numerator.magnitude(),b._numerator.magnitude());
        const BigUnsigned h=gcd_sparse_radix(
            a._denominator,b._denominator);
        BigUnsigned numerator=multiply_cancelled(
            a._numerator.magnitude(),g,b._denominator,h);
        BigUnsigned denominator=multiply_cancelled(
            a._denominator,h,b._numerator.magnitude(),g);
        // Avoid an intermediate Rational and its move/swap, as in R31 mul.
        set_sum(r,Integer::from_magnitude(
            std::move(numerator),a._numerator.negative()!=b._numerator.negative()),
            std::move(denominator));
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

Rational Rational::from_dyadic(Integer numerator,std::size_t bits) {
    if (numerator.magnitude().is_zero()) return Rational{};
    const BigUnsigned& magnitude=numerator.magnitude();
    std::size_t trailing=0u;
    for (auto limb:magnitude.limbs()) {
        if (limb!=0u) {
            trailing+=static_cast<std::size_t>(std::countr_zero(limb));
            break;
        }
        trailing+=64u;
    }
    const std::size_t remove=std::min(trailing,bits);
    if (remove!=0u) {
        numerator=Integer::from_magnitude(
            shift_right(magnitude,remove),numerator.negative());
    }
    return Rational(std::move(numerator),
        shift_left(BigUnsigned(1u),bits-remove),Reduced{});
}

Rational::Rational(Integer const& n,Integer const& d) {
    if (d.magnitude().is_zero()) { *this=inf(integer_sign(n)); return; }
    if (n.magnitude().is_zero()) { return; }
    const BigUnsigned g=gcd_sparse_radix(n.magnitude(),d.magnitude());
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
Rational::Rational(FloatMP const& x) : Rational(x.operator Rational()) { }
Rational::Rational(FloatAP const& x) : Rational(x.operator Rational()) { }

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
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) && defined(__SIZEOF_INT128__) && \
    (defined(__GNUC__) || defined(__clang__))
    // R12: for two single-limb magnitudes, compute the exact binary
    // exponent without allocating/normalizing shifted BigUnsigned
    // values. The ratio always lies in the normal binary64 range.
    // Use the baseline x86-64 DIVQ instruction for the 53-bit quotient
    // rather than the compiler's general unsigned 128-bit division.
    if (n.size()==1u && _denominator.size()==1u) {
        const std::uint64_t numerator=n.limb(0u);
        const std::uint64_t denominator=_denominator.limb(0u);
        const int candidate=static_cast<int>(std::bit_width(numerator))-
                            static_cast<int>(std::bit_width(denominator));
        int exponent=candidate;
        // Shifts fit in 64 bits: the shifted operand has the same
        // bit length as the unshifted operand being compared against.
        const bool below=candidate>=0
            ? numerator<(denominator<<static_cast<unsigned>(candidate))
            : (numerator<<static_cast<unsigned>(-candidate))<denominator;
        if (below) { --exponent; }
        const int scale=exponent-52;
        if (scale>0) {
            const auto significand=numerator/
                (denominator<<static_cast<unsigned>(scale));
            return std::copysign(
                std::ldexp(static_cast<double>(significand),scale),sign);
        }
        __extension__ using Wide=unsigned __int128;
        const Wide scaled=static_cast<Wide>(numerator)<<
                          static_cast<unsigned>(-scale);
        std::uint64_t significand=static_cast<std::uint64_t>(scaled);
        std::uint64_t remainder=static_cast<std::uint64_t>(scaled>>64u);
        if (remainder<denominator) {
            __asm__ volatile("divq %[divisor]"
                : "+a"(significand), "+d"(remainder)
                : [divisor] "r"(denominator)
                : "cc");
            return std::copysign(
                std::ldexp(static_cast<double>(significand),scale),sign);
        }
        // Defensive portable division for any unexpected bound.
        const auto digits=static_cast<std::uint64_t>(scaled/denominator);
        return std::copysign(std::ldexp(static_cast<double>(digits),scale),sign);
    }
#endif
    const std::size_t nb=n.bit_length(),db=_denominator.bit_length();
    if (nb>db && nb-db>1024u) { return sign*std::numeric_limits<double>::infinity(); }
    if (db>nb && db-nb>1075u) { return 0.0; }
    int exponent=nb>=db ? static_cast<int>(nb-db) : -static_cast<int>(db-nb);
    const int order=exponent>=0
        ? compare_shifted_unsigned(n,_denominator,
                                   static_cast<std::size_t>(exponent))
        : -compare_shifted_unsigned(_denominator,n,
                                    static_cast<std::size_t>(-exponent));
    if (order<0) { --exponent; }
    if (exponent>1023) { return sign*std::numeric_limits<double>::infinity(); }
    if (exponent< -1074) { return 0.0; }
    // Truncate toward zero, including subnormals, independent of the active
    // floating rounding mode. The <=53-bit integer is exactly representable.
    const int scale=std::max(exponent-52,-1074);
#if defined(ARIADNE_BIG_UNSIGNED_ARM64) && defined(__SIZEOF_INT128__) && \
    (defined(__GNUC__) || defined(__clang__))
    // R13: for medium equal-width normal ratios, estimate the 53-bit
    // quotient from the highest limb; accept it ONLY after verifying
    // q*d <= n*2^shift < (q+1)*d with full exact integer arithmetic.
    // A rounding-mode-dependent initial estimate is harmless because
    // the certificate, not the estimate, controls the returned value.
    if (exponent>=-1022 && scale<=0 && scale>=-63 &&
            n.size()==_denominator.size() &&
            n.size()>4u && n.size()<=8u) {
        using Limb=BigUnsigned::Limb;
        constexpr std::size_t capacity=10u;
        const std::size_t count=n.size();
        const unsigned shift=static_cast<unsigned>(-scale);
        const double approximate=std::ldexp(
            static_cast<double>(n.limb(count-1u))/
            static_cast<double>(_denominator.limb(count-1u)),
            static_cast<int>(shift));
        constexpr std::uint64_t smallest=std::uint64_t{1}<<52u;
        constexpr std::uint64_t largest=(std::uint64_t{1}<<53u)-1u;
        // Avoid converting a floating value outside the uint64 domain.
        // Clamping also handles an estimate rounding to 2^53 at a
        // representational boundary. Exact tests decide acceptance.
        if (std::isfinite(approximate) && approximate>=0.0) {
            const std::uint64_t rounded_significand=static_cast<std::uint64_t>(
                std::clamp(approximate,static_cast<double>(smallest),
                          static_cast<double>(largest)));
            const std::uint64_t guessed=std::min(rounded_significand,largest);
            std::array<Limb,capacity> scaled{};
            std::array<Limb,capacity> product{};
            const unsigned offset=shift/64u;
            const unsigned low_shift=shift%64u;
            for (std::size_t i=0u;i<count;++i) {
                scaled[i+offset]|=n.limb(i)<<low_shift;
                if (low_shift!=0u) {
                    scaled[i+offset+1u]|=n.limb(i)>>(64u-low_shift);
                }
            }
            const auto compared=[&](std::uint64_t coefficient) {
                product.fill(0u);
                const std::array<Limb,1u> multiplier={coefficient};
                BigUnsigned::mul_into(_denominator.limbs(),multiplier,
                    BigUnsigned::MutableLimbView(product.data(),count+1u));
                for (std::size_t i=capacity;i!=0u;--i) {
                    if (scaled[i-1u]!=product[i-1u]) {
                        return scaled[i-1u]<product[i-1u] ? -1 : 1;
                    }
                }
                return 0;
            };
            if (compared(guessed)>=0 && compared(guessed+1u)<0) {
                return std::copysign(
                    std::ldexp(static_cast<double>(guessed),scale),sign);
            }
        }
    }
#endif
// R8: the native 128/128 prefix division significantly regressed ARM64
// while accelerating x86-64. Keep the exact R7 certificate on x86-64
// and use the earlier R5/fallback path on ARM64 and other targets.
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) && defined(__SIZEOF_INT128__) && \
    (defined(__GNUC__) || defined(__clang__))
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
        // R29: when both two-limb prefixes have a zero low word,
        // the exact 128/128 prefix division simplifies to 64/64.
        // N=high_n*B and D=high_d*B; with s=-scale,
        // q=floor((high_n*2^s)/high_d), r=(high_n*2^s)%high_d.
        // The full-prefix remainder is r*B. Since q<2^53<B,
        // r>0 proves its lower bound; q==2^52 is also proven
        // by the exponent. The upper bound is strict because
        // (high_d-r)*B>=B>2^s. Otherwise use R7's normal path.
        if (n.limb(high-1u)==0u && _denominator.limb(high-1u)==0u &&
                high_n<=static_cast<Wide>(
                    std::numeric_limits<std::uint64_t>::max()>>shift)) {
            const std::uint64_t scaled=
                static_cast<std::uint64_t>(high_n)<<shift;
            const std::uint64_t leading_d=static_cast<std::uint64_t>(high_d);
            const std::uint64_t candidate=scaled/leading_d;
            const std::uint64_t remainder=scaled%leading_d;
            constexpr std::uint64_t min_significand=std::uint64_t{1}<<52u;
            constexpr std::uint64_t max_significand=std::uint64_t{1}<<53u;
            if (candidate>=min_significand &&
                    candidate<max_significand &&
                    (remainder!=0u || candidate==min_significand)) {
                return std::copysign(
                    std::ldexp(static_cast<double>(candidate),scale),sign);
            }
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


static_assert(Field<Rational>);
static_assert(Lattice<Rational>);

} // namespace Ariadne
