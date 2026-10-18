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

} // namespace Ariadne::BigUnsignedDetail
