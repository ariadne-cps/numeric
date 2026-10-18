/***************************************************************************
 *            numeric/big_unsigned_x86-64.cpp
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

/*! \file numeric/big_unsigned_x86-64.cpp
 *  \brief Optional x86-64 baseline multiply kernel; no BMI2/ADX required.
 */

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <utility>

#if defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
#elif !((defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__))
#error "The BigUnsigned x86-64 kernel requires a supported x86-64 compiler"
#endif

namespace Ariadne::BigUnsignedDetail {

namespace {

using Limb = BigUnsigned::Limb;

// Every supported x86-64 CPU has a 64x64->128 unsigned product.
// MSVC offers _umul128; GCC/Clang lower unsigned __int128 to the
// baseline instruction. No optional CPU extension or runtime CPUID.
// All shifts are on unsigned values; no signed-overflow assumptions.
std::pair<Limb,Limb> multiply_wide_x86_64(Limb a,Limb b) noexcept
{
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned __int64 high = 0u;
    const unsigned __int64 low = _umul128(a,b,&high);
    return {low,high};
#else
    __extension__ using Wide = unsigned __int128;
    const Wide product = static_cast<Wide>(a)*static_cast<Wide>(b);
    return {static_cast<Limb>(product),
            static_cast<Limb>(product>>64u)};
#endif
}

} // namespace

// Same radix-B schoolbook recurrence as the portable reference, B=2^64.
// Every inner step satisfies product + existing_limb + carry < B^2,
// so high_product + two overflow bits is representable in one limb.
// Source spans and output must not overlap; all output words are set.
void mul_low_into_x86_64(BigUnsigned::LimbView a,
                         BigUnsigned::LimbView b,
                         BigUnsigned::MutableLimbView output) noexcept
{
    std::fill(output.begin(),output.end(),Limb{0});

    for (std::size_t i=0u; i<a.size() && i<output.size(); ++i) {
        Limb carry = 0u;
        const std::size_t terms = std::min(b.size(),output.size()-i);

        for (std::size_t j=0u; j<terms; ++j) {
            const auto [low,high] = multiply_wide_x86_64(a[i],b[j]);
            const Limb sum = output[i+j]+low;
            const Limb overflow1 = static_cast<Limb>(sum<low);
            const Limb value = sum+carry;
            const Limb overflow2 = static_cast<Limb>(value<sum);

            output[i+j] = value;
            carry = high+overflow1+overflow2;
        }

        if (i+terms < output.size()) {
            output[i+terms] = carry;
        }
    }
}

// Experimental column-wise (Comba) schoolbook multiplication in radix 2^64.
// The 192-bit lo:hi:top accumulator collects the 128-bit partial products
// for one output column plus the propagated high part of the preceding one.
// This structure writes each output digit exactly once. The low-word and
// high-word additions are unsigned, and their overflows extend 'top'.
// In an actual representable operand view, top is bounded by the number of
// products in a column; even in the all-ones case it fits a limb.
// No allocation, global state, BMI2, ADX or alignment beyond uint64_t.
// Input and output spans must be disjoint, including partial overlap.
void mul_low_into_comba_x86_64(BigUnsigned::LimbView a,
                               BigUnsigned::LimbView b,
                               BigUnsigned::MutableLimbView output) noexcept
{
    std::fill(output.begin(),output.end(),Limb{0});
    if (a.empty() || b.empty()) {
        return;
    }

    // The highest nonzero product column is at a.size()+b.size()-2.
    // Avoid overflow in the size sum even for a truncated operation.
    const std::size_t columns =
        a.size()>std::numeric_limits<std::size_t>::max()-b.size()
        ? output.size()
        : std::min(output.size(),a.size()+b.size()-1u);

    Limb lo=0u,hi=0u,top=0u;
    for (std::size_t k=0u; k<columns; ++k) {
        const std::size_t first=k>=b.size() ? k-b.size()+1u : 0u;
        const std::size_t last=std::min(k,a.size()-1u);
        for (std::size_t i=first; i<=last; ++i) {
            const auto [product_lo,product_hi]=
                multiply_wide_x86_64(a[i],b[k-i]);

            const Limb next_lo=lo+product_lo;
            const Limb carry_lo=static_cast<Limb>(next_lo<lo);
            const Limb next_hi=hi+product_hi;
            const Limb carry_hi=static_cast<Limb>(next_hi<hi);
            const Limb high_with_carry=next_hi+carry_lo;
            const Limb carry_hi_extra=static_cast<Limb>(
                high_with_carry<next_hi);

            lo=next_lo;
            hi=high_with_carry;
            top+=carry_hi+carry_hi_extra;
        }

        output[k]=lo;
        lo=hi;
        hi=top;
        top=0u;
    }
    if (columns<output.size()) {
        output[columns]=lo;
    }
}



// Experimental GNU/Clang x86-64 row-schoolbook with explicit carry flags.
// Equivalent to the selected shared wide-accumulator fused recurrence:
// total = a[i]*b[j] + output[i+j] + carry, total < 2^128.
// MULQ writes RDX:RAX; ADDQ/ADCQ propagate each of the two 64-bit
// additions into the high half. +&a and =&d forbid early-clobber
// overlap with inputs: incorrect allocation could corrupt a multiplier.
// No BMI2/ADX, extra scratch, allocations, or runtime CPU detection.
// Output and input spans must not alias; all output words are overwritten.
// On MSVC (no GNU inline asm), the existing x86-64 B3 path is reused.
void mul_low_into_carry_x86_64(BigUnsigned::LimbView a,
                               BigUnsigned::LimbView b,
                               BigUnsigned::MutableLimbView output) noexcept
{
#if (defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__) && !defined(_MSC_VER)
    std::fill(output.begin(),output.end(),Limb{0});
    for (std::size_t i=0u; i<a.size() && i<output.size(); ++i) {
        const Limb ai=a[i];
        Limb carry=0u;
        const std::size_t terms=std::min(b.size(),output.size()-i);
        for (std::size_t j=0u; j<terms; ++j) {
            Limb lo=ai,hi;
            const Limb previous=output[i+j];
            const Limb bj=b[j];
            __asm__ volatile(
                "mulq %[bj]\n\t"
                "addq %[previous], %%rax\n\t"
                "adcq $0, %%rdx\n\t"
                "addq %[carry], %%rax\n\t"
                "adcq $0, %%rdx\n\t"
                : "+&a" (lo), "=&d" (hi)
                : [bj] "r" (bj), [previous] "r" (previous),
                  [carry] "r" (carry)
                : "cc");
            output[i+j]=lo;
            carry=hi;
        }
        if (i+terms<output.size()) {
            output[i+terms]=carry;
        }
    }
#else
    mul_low_into_x86_64(a,b,output);
#endif
}

} // namespace Ariadne::BigUnsignedDetail
