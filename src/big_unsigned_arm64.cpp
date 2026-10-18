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
 *  \brief Selected native fused schoolbook for AArch64 GCC/Clang/AppleClang.
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

// Experimental ARMv8-A row-schoolbook with explicit two-carry handling.
// Unlike the selected __int128 implementation, inline assembly fixes the
// sequence: MUL, UMULH, two ADDS/ADC pairs per product. Performance is
// unknown until integrated macOS ARM64 measurements.
// The intermediate is bounded by (B-1)^2 + 2(B-1) = B^2-1, B=2^64.
// Inputs/output cannot overlap. Output is fully overwritten, no allocation.
void mul_low_into_carry_arm64(BigUnsigned::LimbView a,
                              BigUnsigned::LimbView b,
                              BigUnsigned::MutableLimbView output) noexcept
{
    std::fill(output.begin(),output.end(),Limb{0});
    for (std::size_t i=0u; i<a.size() && i<output.size(); ++i) {
        const Limb ai=a[i];
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

} // namespace Ariadne::BigUnsignedDetail
