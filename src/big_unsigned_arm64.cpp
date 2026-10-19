/***************************************************************************
 *            numeric/big_unsigned_arm64.cpp
 *
 *  Copyright  2026  Ariadne contributors
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

/*! \file numeric/big_unsigned_arm64.cpp
 *  \brief Selected native carry-chain schoolbook for AArch64 GNU/Clang.
 */

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <cstdint>

#if !(defined(__aarch64__) || defined(_M_ARM64))
#error "The ARM64 multiplication kernel requires an AArch64 target"
#endif

#if !(defined(__GNUC__) || defined(__clang__))
#error "The ARM64 multiplication kernel requires GNU or Clang __int128"
#endif

namespace Ariadne::BigUnsignedDetail {

using Limb = BigUnsigned::Limb;

// Native ARMv8-A row-schoolbook with explicit two-carry handling.
// Unlike the shared __int128 fused implementation, inline assembly fixes
// the sequence: MUL, UMULH, two ADDS/ADC pairs per product. This path
// demonstrated integrated speedups on Apple Silicon at 1024-4096 bits.
// The intermediate is bounded by (B-1)^2 + 2(B-1) = B^2-1, B=2^64.
// Inputs/output cannot overlap. Output is fully overwritten, no allocation.
// B38 skips exact-zero outer limbs, retaining the dense carry-chain
// sequence unchanged for every nonzero row (including truncated products).
void mul_low_into_carry_arm64(BigUnsigned::LimbView a,
                              BigUnsigned::LimbView b,
                              BigUnsigned::MutableLimbView output) noexcept
{
    std::fill(output.begin(),output.end(),Limb{0});
    for (std::size_t i=0u; i<a.size() && i<output.size(); ++i) {
        const Limb ai=a[i];
        // B38: a zero row contributes no product or carry. Output was
        // cleared before the loop, and this row's final carry position
        // has not been visited by any earlier row. Skip every inner
        // multiplication without changing truncated/full semantics.
        if (ai==0u) {
            continue;
        }
        Limb carry=0u;
        const std::size_t terms=std::min(b.size(),output.size()-i);
        for (std::size_t j=0u; j<terms; ++j) {
            Limb lo,hi;
            const Limb bj=b[j];
            const Limb previous=output[i+j];
            __asm__ volatile(
                "mul %x[lo], %x[ai], %x[bj]\n\t"
                "umulh %x[hi], %x[ai], %x[bj]\n\t"
                "adds %x[lo], %x[lo], %x[previous]\n\t"
                "adc %x[hi], %x[hi], xzr\n\t"
                "adds %x[lo], %x[lo], %x[carry]\n\t"
                "adc %x[hi], %x[hi], xzr\n\t"
                : [lo] "=&r" (lo), [hi] "=&r" (hi)
                : [ai] "r" (ai), [bj] "r" (bj),
                  [previous] "r" (previous), [carry] "r" (carry)
                : "cc");
            output[i+j]=lo;
            carry=hi;
        }
        if (i+terms<output.size()) {
            output[i+terms]=carry;
        }
    }
}


// B37: symmetric AArch64 full square. Off-diagonal products a[i]*a[j]
// with i<j are computed once using the same native carry-chain
// instruction sequence as the selected generic ARM64 multiplier.
// After doubling the off-diagonal sum, add the n diagonal products.
// This uses n*(n+1)/2 wide products instead of n*n. The two linear
// carry passes are exact in radix B=2^64; the final carry is zero
// since a^2 < B^(2*n). No heap allocation, new ISA extension,
// global state or compiler-specific optimization setting is needed.
// Contract: n>=8; output.size()==2*n; all views are disjoint.
void square_into_symmetric_arm64(BigUnsigned::LimbView a,
                                  BigUnsigned::MutableLimbView output) noexcept
{
    __extension__ using Wide=unsigned __int128;
    const std::size_t n=a.size();
    std::fill(output.begin(),output.end(),Limb{0});

    // Accumulate the strictly upper triangle, one product per pair.
    // output[i+n] lies past all columns previously touched by row i.
    for (std::size_t i=0u;i+1u<n;++i) {
        const Limb ai=a[i];
        // B38: zero off-diagonal rows leave an untouched carry slot
        // output[i+n]. The later doubling and diagonal accumulation
        // still execute for every limb, including these zero rows.
        if (ai==0u) {
            continue;
        }
        Limb carry=0u;
        for (std::size_t j=i+1u;j<n;++j) {
            const std::size_t k=i+j;
            const Limb bj=a[j];
            const Limb previous=output[k];
            Limb lo,hi;
            __asm__ volatile(
                "mul %x[lo], %x[ai], %x[bj]\n\t"
                "umulh %x[hi], %x[ai], %x[bj]\n\t"
                "adds %x[lo], %x[lo], %x[previous]\n\t"
                "adc %x[hi], %x[hi], xzr\n\t"
                "adds %x[lo], %x[lo], %x[carry]\n\t"
                "adc %x[hi], %x[hi], xzr\n\t"
                : [lo] "=&r" (lo), [hi] "=&r" (hi)
                : [ai] "r" (ai), [bj] "r" (bj),
                  [previous] "r" (previous), [carry] "r" (carry)
                : "cc");
            output[k]=lo;
            carry=hi;
        }
        output[i+n]=carry;
    }

    // Add the symmetric lower triangle by doubling cross terms.
    Limb carry=0u;
    for (std::size_t k=0u;k<output.size();++k) {
        const Limb word=output[k];
        output[k]=(word<<1u)|carry;
        carry=word>>63u;
    }

    // Diagonal products occupy disjoint pairs of positions before
    // carry propagation. Wide additions preserve all carries.
    carry=0u;
    for (std::size_t i=0u;i<n;++i) {
        const Wide diagonal=static_cast<Wide>(a[i])*a[i];
        const std::size_t k=2u*i;
        const Wide low=static_cast<Wide>(output[k])+
                       static_cast<Limb>(diagonal)+carry;
        output[k]=static_cast<Limb>(low);
        carry=static_cast<Limb>(low>>64u);
        const Wide high=static_cast<Wide>(output[k+1u])+
                        static_cast<Limb>(diagonal>>64u)+carry;
        output[k+1u]=static_cast<Limb>(high);
        carry=static_cast<Limb>(high>>64u);
    }
}

} // namespace Ariadne::BigUnsignedDetail
