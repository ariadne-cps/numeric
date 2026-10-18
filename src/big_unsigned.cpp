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

std::span<BigUnsigned::Limb const> BigUnsigned::limbs() const noexcept {
    return _size <= inline_capacity
        ? std::span<Limb const>(_inline_limbs.data(),_size)
        : std::span<Limb const>(_heap.data(),_size);
}

BigUnsigned::Limb BigUnsigned::limb(std::size_t index) const noexcept {
    return index < _size ? limbs()[index] : 0u;
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

BigUnsigned add(BigUnsigned const& a, BigUnsigned const& b) {
    const std::size_t size = std::max(a._size,b._size);

    if (size == std::numeric_limits<std::size_t>::max()) {
        throw std::length_error("BigUnsigned addition size overflow");
    }

    BigUnsigned result;
    result.resize_for_result(size+1u);
    auto* digits = result.mutable_data();

    BigUnsigned::Limb carry = 0u;
    for (std::size_t i=0u; i<size; ++i) {
        const auto left = a.limb(i);
        const auto right = b.limb(i);
        const auto sum = left+right;
        const auto value = sum+carry;

        carry = static_cast<BigUnsigned::Limb>(sum<left || value<sum);
        digits[i] = value;
    }

    digits[size] = carry;
    result.normalize();

    return result;
}

BigUnsigned sub(BigUnsigned const& a, BigUnsigned const& b) {
    if (compare(a,b) < 0) {
        throw std::underflow_error("BigUnsigned negative result");
    }

    BigUnsigned result;
    result.resize_for_result(a._size);
    auto* digits = result.mutable_data();

    BigUnsigned::Limb borrow = 0u;
    for (std::size_t i=0u; i<a._size; ++i) {
        const auto left = a.limb(i);
        const auto right = b.limb(i);
        const auto difference = left-right;
        const auto value = difference-borrow;

        borrow = static_cast<BigUnsigned::Limb>(left<right || difference<borrow);
        digits[i] = value;
    }

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

} // namespace Ariadne
