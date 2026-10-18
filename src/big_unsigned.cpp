/***************************************************************************
 *            numeric/big_unsigned.cpp
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

/*! \file numeric/big_unsigned.cpp
 *  \brief Portable representation and foundational unsigned limb operations.
 */

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Ariadne {

#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
namespace BigUnsignedDetail {
void mul_low_into_x86_64(BigUnsigned::LimbView a,
                         BigUnsigned::LimbView b,
                         BigUnsigned::MutableLimbView output) noexcept;
void mul_low_into_comba_x86_64(BigUnsigned::LimbView a,
                               BigUnsigned::LimbView b,
                               BigUnsigned::MutableLimbView output) noexcept;
void mul_low_into_fused_x86_64(BigUnsigned::LimbView a,
                               BigUnsigned::LimbView b,
                               BigUnsigned::MutableLimbView output) noexcept;
} // namespace BigUnsignedDetail
#endif

#if defined(ARIADNE_BIG_UNSIGNED_ARM64)
namespace BigUnsignedDetail {
void mul_low_into_fused_arm64(BigUnsigned::LimbView a,
                              BigUnsigned::LimbView b,
                              BigUnsigned::MutableLimbView output) noexcept;
} // namespace BigUnsignedDetail
#endif

BigUnsigned::MutableLimbView BigUnsigned::Scratch::acquire(std::size_t count) {
    if (count <= inline_capacity) {
        std::fill_n(_inline_limbs.data(),count,Limb{0});
        return {_inline_limbs.data(),count};
    }

    // A new request reuses existing vector capacity whenever possible.
    _heap.resize(count);
    std::fill(_heap.begin(),_heap.end(),Limb{0});
    return {_heap.data(),count};
}

std::size_t BigUnsigned::Scratch::capacity() const noexcept {
    return std::max(inline_capacity,_heap.capacity());
}

BigUnsigned::BigUnsigned(Limb value) noexcept {
    _inline_limbs[0] = value;
    _size = value != 0u ? 1u : 0u;
}

BigUnsigned::BigUnsigned(std::span<Limb const> digits) {
    set_limbs(digits);
}

BigUnsigned::BigUnsigned(BigUnsigned&& other) noexcept
    : _inline_limbs(other._inline_limbs),
      _heap(std::move(other._heap)),
      _size(other._size)
{
    other._size = 0u;
}

BigUnsigned& BigUnsigned::operator=(BigUnsigned&& other) noexcept {
    if (this != &other) {
        _inline_limbs = other._inline_limbs;
        _heap = std::move(other._heap);
        _size = other._size;

        other._size = 0u;
    }

    return *this;
}

BigUnsigned::Limb* BigUnsigned::mutable_data() noexcept {
    return _size <= inline_capacity ? _inline_limbs.data() : _heap.data();
}

void BigUnsigned::resize_for_result(std::size_t limbs) {
    if (limbs > inline_capacity) {
        _heap.resize(limbs,0u);
    } else {
        _inline_limbs.fill(0u);
    }

    _size = limbs;
}

void BigUnsigned::normalize() noexcept {
    auto* digits = mutable_data();

    while (_size != 0u && digits[_size-1u] == 0u) {
        --_size;
    }

    if (_size <= inline_capacity && !_heap.empty()) {
        // Read the significant digits before releasing heap storage.
        std::copy_n(_heap.data(),_size,_inline_limbs.data());
        _heap.clear();
    }
}

BigUnsigned::LimbView BigUnsigned::limbs() const noexcept {
    return _size <= inline_capacity
        ? std::span<Limb const>(_inline_limbs.data(),_size)
        : std::span<Limb const>(_heap.data(),_size);
}

BigUnsigned::Limb BigUnsigned::limb(std::size_t index) const noexcept {
    return index < _size ? limbs()[index] : 0u;
}

void BigUnsigned::assign(LimbView digits) {
    set_limbs(digits);
}

void BigUnsigned::set_limbs(std::span<Limb const> digits) {
    // The empty sequence is the unique representation of zero.
    while (!digits.empty() && digits.back() == 0u) {
        digits = digits.first(digits.size()-1u);
    }

    if (digits.size() <= inline_capacity) {
        // Copy before mutating state: digits may refer to our own storage.
        std::array<Limb,inline_capacity> replacement{};
        std::copy(digits.begin(),digits.end(),replacement.begin());

        _inline_limbs = replacement;
        _size = digits.size();
        _heap.clear();
    } else {
        // Construct first so allocation failures leave this value intact.
        std::vector<Limb> replacement(digits.begin(),digits.end());

        _heap.swap(replacement);
        _size = digits.size();
    }
}

std::size_t BigUnsigned::bit_length() const noexcept {
    if (_size == 0u) {
        return 0u;
    }

    return 64u * (_size-1u) +
        64u - static_cast<unsigned>(std::countl_zero(limb(_size-1u)));
}

void BigUnsigned::swap(BigUnsigned& other) noexcept {
    using std::swap;

    swap(_inline_limbs,other._inline_limbs);
    _heap.swap(other._heap);
    swap(_size,other._size);
}

int compare(BigUnsigned const& a, BigUnsigned const& b) noexcept {
    if (a._size != b._size) {
        return a._size < b._size ? -1 : 1;
    }

    for (std::size_t i=a._size; i!=0u; --i) {
        if (a.limb(i-1u) != b.limb(i-1u)) {
            return a.limb(i-1u) < b.limb(i-1u) ? -1 : 1;
        }
    }

    return 0;
}

BigUnsigned::Limb BigUnsigned::add_into(
        LimbView a,LimbView b,MutableLimbView output)
{
    const std::size_t count = std::max(a.size(),b.size());
    if (output.size() < count) {
        throw std::length_error("BigUnsigned add_into output too small");
    }

    Limb carry = 0u;
    for (std::size_t i=0u; i<count; ++i) {
        const Limb left = i<a.size() ? a[i] : 0u;
        const Limb right = i<b.size() ? b[i] : 0u;
        const Limb sum = left+right;
        const Limb value = sum+carry;

        carry = static_cast<Limb>(sum<left || value<sum);
        output[i] = value;
    }

    std::fill(output.subspan(count).begin(),output.end(),Limb{0});
    return carry;
}

BigUnsigned::Limb BigUnsigned::sub_into(
        LimbView a,LimbView b,MutableLimbView output)
{
    const std::size_t count = std::max(a.size(),b.size());
    if (output.size() < count) {
        throw std::length_error("BigUnsigned sub_into output too small");
    }

    Limb borrow = 0u;
    for (std::size_t i=0u; i<count; ++i) {
        const Limb left = i<a.size() ? a[i] : 0u;
        const Limb right = i<b.size() ? b[i] : 0u;
        const Limb difference = left-right;
        const Limb value = difference-borrow;

        borrow = static_cast<Limb>(left<right || difference<borrow);
        output[i] = value;
    }

    std::fill(output.subspan(count).begin(),output.end(),Limb{0});
    return borrow;
}

BigUnsigned add(BigUnsigned const& a, BigUnsigned const& b) {
    const std::size_t count = std::max(a._size,b._size);
    if (count == std::numeric_limits<std::size_t>::max()) {
        throw std::length_error("BigUnsigned addition size overflow");
    }

    BigUnsigned result;
    result.resize_for_result(count+1u);
    auto out = BigUnsigned::MutableLimbView(result.mutable_data(),count+1u);
    out[count] = BigUnsigned::add_into(a.limbs(),b.limbs(),out.first(count));
    result.normalize();
    return result;
}

BigUnsigned sub(BigUnsigned const& a, BigUnsigned const& b) {
    if (compare(a,b)<0) {
        throw std::underflow_error("BigUnsigned negative result");
    }

    BigUnsigned result;
    result.resize_for_result(a._size);
    const auto out = BigUnsigned::MutableLimbView(result.mutable_data(),a._size);
    (void)BigUnsigned::sub_into(a.limbs(),b.limbs(),out);
    result.normalize();
    return result;
}

BigUnsigned shift_left(BigUnsigned const& a, std::size_t bits) {
    if (a.is_zero()) {
        return BigUnsigned{};
    }

    const std::size_t words = bits/64u;
    const unsigned offset = static_cast<unsigned>(bits%64u);

    if (words >= std::numeric_limits<std::size_t>::max()-a._size-1u) {
        throw std::length_error("BigUnsigned left shift size overflow");
    }

    BigUnsigned result;
    result.resize_for_result(a._size+words+1u);
    auto* digits = result.mutable_data();

    BigUnsigned::Limb carry = 0u;
    for (std::size_t i=0u; i<a._size; ++i) {
        const auto value = a.limb(i);

        digits[words+i] = (value<<offset)|carry;
        carry = offset == 0u ? 0u : (value>>(64u-offset));
    }

    digits[words+a._size] = carry;
    result.normalize();

    return result;
}

BigUnsigned shift_right(BigUnsigned const& a, std::size_t bits) {
    const std::size_t words = bits/64u;

    if (words >= a._size) {
        return BigUnsigned{};
    }

    const unsigned offset = static_cast<unsigned>(bits%64u);
    BigUnsigned result;
    result.resize_for_result(a._size-words);
    auto* digits = result.mutable_data();

    for (std::size_t i=words; i<a._size; ++i) {
        const auto low = a.limb(i);
        const auto high = i+1u < a._size ? a.limb(i+1u) : 0u;

        digits[i-words] = offset == 0u ? low :
            ((low>>offset)|(high<<(64u-offset)));
    }

    result.normalize();

    return result;
}

// Portable 64x64 -> 128 multiply using only standard 64-bit arithmetic.
// Represent operands as a1*2^32 + a0 and b1*2^32 + b0.
namespace {

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

bool BigUnsigned::has_native_multiply() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) || defined(ARIADNE_BIG_UNSIGNED_ARM64)
    return true;
#else
    return false;
#endif
}

bool BigUnsigned::has_x86_64_multiply() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    return true;
#else
    return false;
#endif
}

bool BigUnsigned::selected_multiply_uses_fused() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) || defined(ARIADNE_BIG_UNSIGNED_ARM64_FUSED)
    return true;
#else
    return false;
#endif
}

void BigUnsigned::mul_low_into_comba(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    BigUnsignedDetail::mul_low_into_comba_x86_64(a,b,output);
#else
    mul_low_into_portable(a,b,output);
#endif
}

void BigUnsigned::mul_low_into_fused(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    BigUnsignedDetail::mul_low_into_fused_x86_64(a,b,output);
#elif defined(ARIADNE_BIG_UNSIGNED_ARM64)
    BigUnsignedDetail::mul_low_into_fused_arm64(a,b,output);
#else
    mul_low_into_portable(a,b,output);
#endif
}

void BigUnsigned::mul_low_into_schoolbook(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    // Reproduce the B3 selected behavior, including its one-limb
    // portable fallback, to preserve a strict A/B control.
    if (a.size() >= 2u && b.size() >= 2u) {
        BigUnsignedDetail::mul_low_into_x86_64(a,b,output);
        return;
    }
#endif
    mul_low_into_portable(a,b,output);
}

void BigUnsigned::mul_low_into(LimbView a,LimbView b,MutableLimbView output) {
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED)
    // GCC x86-64's selected fused kernel has measured integrated gains.
    BigUnsignedDetail::mul_low_into_fused_x86_64(a,b,output);
#elif defined(ARIADNE_BIG_UNSIGNED_ARM64_FUSED)
    // ARM64 selection is explicit opt-in until native macOS benchmarks.
    BigUnsignedDetail::mul_low_into_fused_arm64(a,b,output);
#else
    mul_low_into_schoolbook(a,b,output);
#endif
}

namespace {

void check_full_product_size(
        BigUnsigned::LimbView a,BigUnsigned::LimbView b,
        BigUnsigned::MutableLimbView output)
{
    if (a.empty() || b.empty()) {
        return;
    }
    if (a.size() > std::numeric_limits<std::size_t>::max()-b.size()
            || output.size() < a.size()+b.size()) {
        throw std::length_error("BigUnsigned mul_into output too small");
    }
}

} // namespace

void BigUnsigned::mul_into(LimbView a,LimbView b,MutableLimbView output) {
    check_full_product_size(a,b,output);
    mul_low_into(a,b,output);
}

void BigUnsigned::mul_into_portable(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_portable(a,b,output);
}

void BigUnsigned::mul_into_comba(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_comba(a,b,output);
}

void BigUnsigned::mul_into_fused(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_fused(a,b,output);
}

void BigUnsigned::mul_into_schoolbook(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_schoolbook(a,b,output);
}

BigUnsigned mul_low(BigUnsigned const& a,BigUnsigned const& b,
                    std::size_t limb_count)
{
    if (a.is_zero() || b.is_zero() || limb_count==0u) {
        return BigUnsigned{};
    }

    if (a._size > std::numeric_limits<std::size_t>::max()-b._size) {
        throw std::length_error("BigUnsigned product size overflow");
    }

    const std::size_t result_size = std::min(limb_count,a._size+b._size);
    BigUnsigned result;
    result.resize_for_result(result_size);
    BigUnsigned::mul_low_into(
        a.limbs(),b.limbs(),
        BigUnsigned::MutableLimbView(result.mutable_data(),result_size));
    result.normalize();
    return result;
}

BigUnsigned mul(BigUnsigned const& a,BigUnsigned const& b)
{
    if (a.is_zero() || b.is_zero()) {
        return BigUnsigned{};
    }

    if (a._size > std::numeric_limits<std::size_t>::max()-b._size) {
        throw std::length_error("BigUnsigned product size overflow");
    }

    BigUnsigned result;
    result.resize_for_result(a._size+b._size);
    BigUnsigned::mul_into(
        a.limbs(),b.limbs(),
        BigUnsigned::MutableLimbView(result.mutable_data(),a._size+b._size));
    result.normalize();
    return result;
}

BigUnsigned square(BigUnsigned const& a)
{
    return mul(a,a);
}

} // namespace Ariadne
