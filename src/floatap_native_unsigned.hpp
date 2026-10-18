/***************************************************************************
 *            numeric/floatap_native_unsigned.hpp
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

/*! \file numeric/floatap_native_unsigned.hpp
 *  \brief GMP-free unsigned arithmetic for the eventual FloatAP significand.
 */

#ifndef ARIADNE_NUMERIC_FLOATAP_NATIVE_UNSIGNED_HPP
#define ARIADNE_NUMERIC_FLOATAP_NATIVE_UNSIGNED_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Ariadne::FloatAPNative {

// This type is private to FloatAP: it is not a proposed Integer/Rational
// facade. Little-endian 32-bit limbs keep all products within uint64_t
// on GCC, Clang and MSVC without compiler-specific wide integers.
//
// The first four limbs are inline. A dynamic buffer is used only when
// the normalized magnitude exceeds the inline capacity.
class Unsigned {
  public:
    using Limb = std::uint32_t;
    static constexpr std::size_t inline_capacity = 4u;
    static constexpr std::uint64_t base = std::uint64_t{1} << 32;

  private:
    std::array<Limb,inline_capacity> _inline{};
    std::vector<Limb> _large;
    std::size_t _size = 0u;

    Limb* mutable_data() noexcept {
        return _size <= inline_capacity ? _inline.data() : _large.data();
    }

  public:
    Unsigned() = default;

    explicit Unsigned(std::uint64_t value) {
        _inline[0] = static_cast<Limb>(value);
        _inline[1] = static_cast<Limb>(value >> 32);
        _size = _inline[1] != 0u ? 2u : (_inline[0] != 0u ? 1u : 0u);
    }

    std::size_t size() const noexcept { return _size; }
    bool is_zero() const noexcept { return _size == 0u; }

    Limb const* data() const noexcept {
        return _size <= inline_capacity ? _inline.data() : _large.data();
    }

    Limb word(std::size_t index) const noexcept {
        return index < _size ? data()[index] : 0u;
    }

    // Resize is also used for scratch buffers. The caller must trim()
    // any newly written leading zero limbs before comparison/return.
    void resize(std::size_t size) {
        if (size <= inline_capacity) {
            if (_size > inline_capacity) {
                std::copy_n(_large.data(),size,_inline.data());
                _large.clear();
            } else if (size > _size) {
                std::fill(_inline.begin()+static_cast<std::ptrdiff_t>(_size),
                          _inline.begin()+static_cast<std::ptrdiff_t>(size),0u);
            }
        } else if (_size <= inline_capacity) {
            _large.assign(size,0u);
            std::copy_n(_inline.data(),_size,_large.data());
        } else {
            _large.resize(size,0u);
        }
        _size = size;
    }

    void trim() noexcept {
        Limb const* const words = data();
        std::size_t trimmed = _size;
        while (trimmed != 0u && words[trimmed-1u] == 0u) { --trimmed; }
        if (_size > inline_capacity && trimmed <= inline_capacity) {
            std::copy_n(_large.data(),trimmed,_inline.data());
            _large.clear();
        }
        _size = trimmed;
    }

    void set_word(std::size_t index,Limb value) {
        if (index >= _size) {
            if (value == 0u) { return; }
            if (index == std::numeric_limits<std::size_t>::max()) {
                throw std::overflow_error("FloatAP native limb index overflow");
            }
            resize(index+1u);
        }
        mutable_data()[index] = value;
        trim();
    }

    void set_bit(std::size_t bit) {
        const std::size_t index=bit/32u;
        if (index >= _size) { resize(index+1u); }
        mutable_data()[index] |= Limb{1} << static_cast<unsigned>(bit%32u);
    }

    bool bit(std::size_t bit) const noexcept {
        const std::size_t index=bit/32u;
        return index < _size &&
            ((word(index) >> static_cast<unsigned>(bit%32u)) & 1u) != 0u;
    }

    std::size_t bit_length() const noexcept {
        return _size == 0u ? 0u :
            32u*(_size-1u) + 32u -
            static_cast<std::size_t>(std::countl_zero(word(_size-1u)));
    }

    friend int compare(Unsigned const& a,Unsigned const& b) noexcept {
        if (a._size != b._size) { return a._size < b._size ? -1 : 1; }
        for (std::size_t i=a._size; i!=0u; --i) {
            if (a.word(i-1u) != b.word(i-1u)) {
                return a.word(i-1u) < b.word(i-1u) ? -1 : 1;
            }
        }
        return 0;
    }

    friend bool operator==(Unsigned const& a,Unsigned const& b) noexcept {
        return compare(a,b)==0;
    }

    friend Unsigned add(Unsigned const& a,Unsigned const& b) {
        const std::size_t size=std::max(a.size(),b.size());
        if (size==std::numeric_limits<std::size_t>::max()) {
            throw std::overflow_error("FloatAP native addition size overflow");
        }
        Unsigned result;
        result.resize(size+1u);
        std::uint64_t carry=0u;
        for (std::size_t i=0u; i<size; ++i) {
            const std::uint64_t sum=std::uint64_t{a.word(i)}+
                                    b.word(i)+carry;
            result.mutable_data()[i]=static_cast<Limb>(sum);
            carry=sum >> 32;
        }
        result.mutable_data()[size]=static_cast<Limb>(carry);
        result.trim();
        return result;
    }

    // Unsigned subtraction requires a >= b. It never wraps at zero.
    friend Unsigned sub(Unsigned const& a,Unsigned const& b) {
        if (compare(a,b)<0) {
            throw std::domain_error("FloatAP native unsigned underflow");
        }
        Unsigned result;
        result.resize(a.size());
        std::uint64_t borrow=0u;
        for (std::size_t i=0u; i<a.size(); ++i) {
            const std::uint64_t subtraction=std::uint64_t{b.word(i)}+borrow;
            const std::uint64_t ai=a.word(i);
            result.mutable_data()[i]=static_cast<Limb>(ai-subtraction);
            borrow=ai<subtraction ? 1u : 0u;
        }
        result.trim();
        return result;
    }

    friend Unsigned shift_left(Unsigned const& value,std::size_t bits) {
        if (value.is_zero()) { return {}; }
        const std::size_t words=bits/32u;
        if (words >= std::numeric_limits<std::size_t>::max()-value.size()-1u) {
            throw std::overflow_error("FloatAP native shift size overflow");
        }
        Unsigned result;
        result.resize(value.size()+words+1u);
        const unsigned small=static_cast<unsigned>(bits%32u);
        std::uint64_t carry=0u;
        for (std::size_t i=0u; i<value.size(); ++i) {
            const std::uint64_t next=(std::uint64_t{value.word(i)}<<small)|carry;
            result.mutable_data()[i+words]=static_cast<Limb>(next);
            carry=next >> 32;
        }
        result.mutable_data()[value.size()+words]=static_cast<Limb>(carry);
        result.trim();
        return result;
    }

    friend Unsigned shift_right(Unsigned const& value,std::size_t bits) {
        const std::size_t words=bits/32u;
        if (value.size()<=words) { return {}; }
        const unsigned small=static_cast<unsigned>(bits%32u);
        Unsigned result;
        result.resize(value.size()-words);
        for (std::size_t i=0u; i<result.size(); ++i) {
            const std::uint32_t low=value.word(i+words);
            result.mutable_data()[i]=small == 0u ? low :
                static_cast<Limb>((low >> small) |
                    (std::uint64_t{value.word(i+words+1u)} << (32u-small)));
        }
        result.trim();
        return result;
    }

    friend Unsigned mul(Unsigned const& a,Unsigned const& b) {
        if (a.is_zero() || b.is_zero()) { return {}; }
        if (a.size()>std::numeric_limits<std::size_t>::max()-b.size()) {
            throw std::overflow_error("FloatAP native product size overflow");
        }
        Unsigned result;
        result.resize(a.size()+b.size());
        for (std::size_t i=0u; i<a.size(); ++i) {
            std::uint64_t carry=0u;
            for (std::size_t j=0u; j<b.size(); ++j) {
                const std::uint64_t product=std::uint64_t{a.word(i)}*b.word(j)
                    +result.mutable_data()[i+j]+carry;
                result.mutable_data()[i+j]=static_cast<Limb>(product);
                carry=product >> 32;
            }
            result.mutable_data()[i+b.size()]=static_cast<Limb>(carry);
        }
        result.trim();
        return result;
    }

    // Exact base-2^32 Knuth division. Output satisfies a = q*b+r,
    // with 0 <= r < b. Invalid zero divisors throw; all inputs stay
    // unchanged. Only standard 64-bit intermediate products are needed.
    friend std::pair<Unsigned,Unsigned> divrem(
            Unsigned const& a,Unsigned const& b)
    {
        if (b.is_zero()) {
            throw std::domain_error("FloatAP native divide by zero");
        }
        if (compare(a,b)<0) { return {{},a}; }
        if (b.size()==1u) {
            Unsigned q;
            q.resize(a.size());
            std::uint64_t r=0u;
            for (std::size_t i=a.size(); i!=0u; --i) {
                const std::uint64_t part=(r<<32)|a.word(i-1u);
                q.mutable_data()[i-1u]=static_cast<Limb>(part/b.word(0u));
                r=part%b.word(0u);
            }
            q.trim();
            return {std::move(q),Unsigned(r)};
        }
        const unsigned normalize=static_cast<unsigned>(
            std::countl_zero(b.word(b.size()-1u)));
        Unsigned v=shift_left(b,normalize);
        Unsigned u=shift_left(a,normalize);
        const std::size_t n=v.size();
        const std::size_t m=u.size()-n;
        if (u.size()==std::numeric_limits<std::size_t>::max()) {
            throw std::overflow_error("FloatAP native division size overflow");
        }
        u.resize(u.size()+1u);
        Unsigned q;
        q.resize(m+1u);
        Limb* const u_words=u.mutable_data();
        Limb* const q_words=q.mutable_data();
        for (std::size_t offset=m+1u; offset!=0u; --offset) {
            const std::size_t j=offset-1u;
            const std::uint64_t top=u_words[j+n];
            const std::uint64_t next=u_words[j+n-1u];
            std::uint64_t estimate,estimate_remainder;
            if (top>=v.word(n-1u)) {
                estimate=base-1u;
                estimate_remainder=next+v.word(n-1u);
            } else {
                const std::uint64_t numerator=(top<<32)|next;
                estimate=numerator/v.word(n-1u);
                estimate_remainder=numerator%v.word(n-1u);
            }
            while (estimate_remainder<base &&
                   estimate*v.word(n-2u)>
                   estimate_remainder*base+u_words[j+n-2u]) {
                --estimate;
                estimate_remainder+=v.word(n-1u);
            }
            std::uint64_t carry=0u,borrow=0u;
            for (std::size_t i=0u; i<n; ++i) {
                const std::uint64_t product=estimate*v.word(i)+carry;
                carry=product>>32;
                const std::uint64_t subtract=(product&(base-1u))+borrow;
                borrow=u_words[j+i]<subtract ? 1u : 0u;
                u_words[j+i]=static_cast<Limb>(
                    std::uint64_t{u_words[j+i]}-subtract);
            }
            const std::uint64_t last=carry+borrow;
            const bool overestimate=u_words[j+n]<last;
            u_words[j+n]=static_cast<Limb>(
                std::uint64_t{u_words[j+n]}-last);
            if (overestimate) {
                --estimate;
                std::uint64_t addition_carry=0u;
                for (std::size_t i=0u; i<n; ++i) {
                    const std::uint64_t sum=std::uint64_t{u_words[j+i]}+
                                            v.word(i)+addition_carry;
                    u_words[j+i]=static_cast<Limb>(sum);
                    addition_carry=sum>>32;
                }
                u_words[j+n]+=static_cast<Limb>(addition_carry);
            }
            q_words[j]=static_cast<Limb>(estimate);
        }
        q.trim();
        Unsigned normalized_remainder;
        normalized_remainder.resize(n);
        for (std::size_t i=0u; i<n; ++i) {
            normalized_remainder.mutable_data()[i]=u_words[i];
        }
        normalized_remainder.trim();
        return {std::move(q),shift_right(normalized_remainder,normalize)};
    }
};

} // namespace Ariadne::FloatAPNative

#endif // ARIADNE_NUMERIC_FLOATAP_NATIVE_UNSIGNED_HPP
