/***************************************************************************
 *            numeric/big_unsigned_fallback.cpp
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

/*! \file numeric/big_unsigned_fallback.cpp
 *  \brief Portable C++20 limb multiplication and reference backend.
 */

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <cstdint>
#include <span>
#include <utility>

namespace Ariadne {

namespace {

// This pure C++20 reference is available on all supported targets, even
// when an optimized architecture backend is selected for production.
// Keeping it compiled makes portable-vs-native tests and A/B benchmarks
// independent of ISA support.
// Portable 64x64 -> 128 multiply using only standard 64-bit arithmetic.
// Represent operands as a1*2^32 + a0 and b1*2^32 + b0.
std::pair<std::uint64_t,std::uint64_t> multiply_wide(
        std::uint64_t a,std::uint64_t b) noexcept
{
    constexpr std::uint64_t mask = 0xffffffffull;
    const std::uint64_t a0 = a & mask;
    const std::uint64_t a1 = a >> 32u;
    const std::uint64_t b0 = b & mask;
    const std::uint64_t b1 = b >> 32u;

    const std::uint64_t p00 = a0*b0;
    const std::uint64_t p01 = a0*b1;
    const std::uint64_t p10 = a1*b0;
    const std::uint64_t p11 = a1*b1;

    // The middle accumulator is at most 3*(2^32-1), so 64 bits suffice.
    const std::uint64_t middle = (p00 >> 32u) + (p01 & mask) + (p10 & mask);
    const std::uint64_t low = (p00 & mask) | (middle << 32u);
    const std::uint64_t high = p11 + (p01 >> 32u) +
                               (p10 >> 32u) + (middle >> 32u);
    return {low,high};
}

// Compute the product modulo 2^(64*out.size()), storing only low limbs.
// The destination must be distinct from both source spans; kernel zeroes it.
// At each step the carry is below 2^64: (B-1)^2 + 2(B-1) = B^2-1.
void multiply_limbs(std::span<std::uint64_t const> a,
                    std::span<std::uint64_t const> b,
                    std::span<std::uint64_t> out) noexcept
{
    std::fill(out.begin(),out.end(),std::uint64_t{0});

    for (std::size_t i=0u; i<a.size() && i<out.size(); ++i) {
        std::uint64_t carry = 0u;
        const std::size_t terms = std::min(b.size(),out.size()-i);

        for (std::size_t j=0u; j<terms; ++j) {
            const auto [low,high] = multiply_wide(a[i],b[j]);
            const std::uint64_t sum = out[i+j]+low;
            const std::uint64_t overflow1 = static_cast<std::uint64_t>(sum<low);
            const std::uint64_t result = sum+carry;
            const std::uint64_t overflow2 = static_cast<std::uint64_t>(result<sum);

            out[i+j] = result;
            carry = high+overflow1+overflow2;
        }

        if (i+terms < out.size()) {
            // No prior row can have written beyond index i+terms-1.
            out[i+terms] = carry;
        }
    }
}

} // namespace

void BigUnsigned::mul_low_into_portable(
        LimbView a,LimbView b,MutableLimbView output)
{
    multiply_limbs(a,b,output);
}

} // namespace Ariadne
