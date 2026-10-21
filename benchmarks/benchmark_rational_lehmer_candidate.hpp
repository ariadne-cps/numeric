/***************************************************************************
 *            benchmark_rational_lehmer_candidate.hpp
 *
 *  Copyright 2026 Ariadne contributors
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

/*! \file benchmark_rational_lehmer_candidate.hpp
 *  \brief R25 nonproduction bounded Lehmer GCD diagnostic.
 *
 *  Never included in production sources. A guessed high-limb quotient is
 *  only used to build a determinant +/-1 integer matrix; that matrix
 *  is applied exactly to entire limb arrays. Invalid, nondecreasing or
 *  oversized results fall back to the verified full remainder kernel.
 */

#ifndef ARIADNE_BENCHMARK_RATIONAL_LEHMER_CANDIDATE_HPP
#define ARIADNE_BENCHMARK_RATIONAL_LEHMER_CANDIDATE_HPP

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <utility>

namespace AriadneBenchmark::LehmerDiagnostic {

constexpr std::size_t capacity=8u;
using Limb=Ariadne::BigUnsigned::Limb;
using Words=std::array<Limb,capacity>;

#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
__extension__ using SignedWide=__int128;
__extension__ using UnsignedWide=unsigned __int128;

inline std::size_t trim(Words const& words,std::size_t count) {
    while (count!=0u && words[count-1u]==0u) { --count; }
    return count;
}

inline int compare(Words const& a,std::size_t na,Words const& b,std::size_t nb) {
    if (na!=nb) { return na<nb ? -1 : 1; }
    for (std::size_t i=na;i!=0u;--i) {
        if (a[i-1u]!=b[i-1u]) { return a[i-1u]<b[i-1u] ? -1 : 1; }
    }
    return 0;
}

// Capture up to 127 leading bits, aligned for equally long operands.
// Discarding one low prefix bit guarantees a positive signed-128 value.
inline SignedWide prefix(Words const& words,std::size_t count) {
    return (static_cast<SignedWide>(words[count-1u])<<63u) |
           static_cast<SignedWide>(words[count-2u]>>1u);
}

// Exact row output from signed coefficients in base 2^64. Each coefficient
// has magnitude <=2^19, so all per-limb products and signed carries fit
// SignedWide. Conversion to Limb is the defined residue modulo 2^64.
inline bool apply_row(
        Words const& a,std::size_t na,Words const& b,std::size_t nb,
        SignedWide m,SignedWide n,Words& out,std::size_t& count)
{
    const SignedWide base=SignedWide{1}<<64u;
    SignedWide carry=0;
    const std::size_t width=std::max(na,nb);
    for (std::size_t i=0u;i<width;++i) {
        const SignedWide ai=i<na ? static_cast<SignedWide>(a[i]) : SignedWide{0};
        const SignedWide bi=i<nb ? static_cast<SignedWide>(b[i]) : SignedWide{0};
        const SignedWide sum=m*ai+n*bi+carry;
        const Limb low=static_cast<Limb>(sum);
        out[i]=low;
        carry=(sum-static_cast<SignedWide>(low))/base;
    }
    if (carry<0 || carry>static_cast<SignedWide>(~Limb{0})) { return false; }
    count=width;
    if (carry!=0) {
        if (count==capacity) { return false; }
        out[count++]=static_cast<Limb>(carry);
    }
    count=trim(out,count);
    return true;
}

// The matrix is a product of elementary Euclidean matrices
// [[0,1],[1,-q]], hence has determinant +/-1 and preserves gcd.
// Each transform is checked on entire original operands before use.
// A strict reduction of the maximum prevents nontermination even if
// truncated prefix quotients are approximate.
inline bool lehmer_batch(Words& a,std::size_t& na,Words& b,std::size_t& nb) {
    if (na!=nb || nb<2u) { return false; }
    SignedWide x=prefix(a,na),y=prefix(b,nb);
    SignedWide aa=1,bb=0,cc=0,dd=1;
    unsigned steps=0u;
    constexpr SignedWide limit=SignedWide{1}<<19u;
    while (steps<20u) {
        const SignedWide den1=y+cc,den2=y+dd;
        if (den1<=0 || den2<=0) { break; }
        const SignedWide num1=x+aa,num2=x+bb;
        if (num1<=0 || num2<=0) { break; }
        const SignedWide q1=num1/den1,q2=num2/den2;
        if (q1!=q2 || q1<1 || q1>10000) { break; }
        const SignedWide nexta=cc,nextb=dd;
        const SignedWide nextc=aa-q1*cc,nextd=bb-q1*dd;
        if (nexta>limit || nexta<-limit ||
                nextb>limit || nextb<-limit ||
                nextc>limit || nextc<-limit ||
                nextd>limit || nextd<-limit ||
                x-q1*y<0) {
            break;
        }
        const SignedWide next_y=x-q1*y;
        x=y; y=next_y;
        aa=nexta;bb=nextb;cc=nextc;dd=nextd;
        ++steps;
    }
    if (steps==0u) { return false; }
    Words next_a{},next_b{};
    std::size_t next_na=0u,next_nb=0u;
    if (!apply_row(a,na,b,nb,aa,bb,next_a,next_na) ||
            !apply_row(a,na,b,nb,cc,dd,next_b,next_nb)) {
        return false;
    }
    if (compare(next_a,next_na,next_b,next_nb)<0) {
        std::swap(next_a,next_b);
        std::swap(next_na,next_nb);
    }
    if (compare(next_a,next_na,a,na)>=0) { return false; }
    a=next_a; b=next_b;
    na=next_na;nb=next_nb;
    return true;
}

inline Limb scalar_remainder(Words const& a,std::size_t na,Limb divisor) {
    UnsignedWide result=0u;
    for (std::size_t i=na;i!=0u;--i) {
        result=((result<<64u)|a[i-1u])%divisor;
    }
    return static_cast<Limb>(result);
}

// Exact fallback using the existing public, alias-free division kernel.
// It produces both outputs, whereas production gcd selects remainder only;
// this cost is intentionally reflected in experimental measurements.
inline void full_remainder(
        Words const& a,std::size_t na,Words const& b,std::size_t nb,
        Words& result,std::size_t& nr)
{
    std::array<Limb,capacity> quotient;
    std::array<Limb,2u*capacity+1u> scratch;
    Ariadne::BigUnsigned::divrem_into(
        Ariadne::BigUnsigned::LimbView(a.data(),na),
        Ariadne::BigUnsigned::LimbView(b.data(),nb),
        Ariadne::BigUnsigned::MutableLimbView(quotient.data(),na),
        Ariadne::BigUnsigned::MutableLimbView(result.data(),nb),
        Ariadne::BigUnsigned::MutableLimbView(
            scratch.data(),Ariadne::BigUnsigned::divrem_scratch_limbs(na,nb)));
    nr=trim(result,nb);
}
#endif

inline Ariadne::BigUnsigned gcd(
        Ariadne::BigUnsigned const& lhs,Ariadne::BigUnsigned const& rhs)
{
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    if (lhs.size()>capacity || rhs.size()>capacity) {
        return Ariadne::gcd(lhs,rhs);
    }
    Words a{},b{};
    std::copy(lhs.limbs().begin(),lhs.limbs().end(),a.begin());
    std::copy(rhs.limbs().begin(),rhs.limbs().end(),b.begin());
    std::size_t na=lhs.size(),nb=rhs.size();
    if (compare(a,na,b,nb)<0) {
        std::swap(a,b);
        std::swap(na,nb);
    }
    while (nb!=0u) {
        if (nb==1u) {
            const Limb common=std::gcd(b[0u],scalar_remainder(a,na,b[0u]));
            return Ariadne::BigUnsigned(common);
        }
        if (lehmer_batch(a,na,b,nb)) { continue; }
        Words remainder{};
        std::size_t nr=0u;
        full_remainder(a,na,b,nb,remainder,nr);
        a=b;na=nb;
        b=remainder;nb=nr;
    }
    return Ariadne::BigUnsigned(Ariadne::BigUnsigned::LimbView(a.data(),na));
#else
    return Ariadne::gcd(lhs,rhs);
#endif
}

constexpr bool available =
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    true;
#else
    false;
#endif

} // namespace AriadneBenchmark::LehmerDiagnostic

#endif // ARIADNE_BENCHMARK_RATIONAL_LEHMER_CANDIDATE_HPP
