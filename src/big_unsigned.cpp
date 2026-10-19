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
 *  \brief Shared unsigned arithmetic, fused multiplication, and Karatsuba.
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
void mul_low_into_carry_x86_64(BigUnsigned::LimbView a,
                               BigUnsigned::LimbView b,
                               BigUnsigned::MutableLimbView output) noexcept;
} // namespace BigUnsignedDetail
#endif

#if defined(ARIADNE_BIG_UNSIGNED_ARM64)
namespace BigUnsignedDetail {
void mul_low_into_carry_arm64(BigUnsigned::LimbView a,
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

// Baseline C++20 division for an unsigned (high:low)/divisor digit.
// Precondition: 0 <= high < divisor and divisor != 0.
// Every iteration appends one bit to the remainder. A 65th bit is
// represented by overflow; in that case subtraction is unconditional.
// Since the pre-subtraction remainder is < 2*divisor, one subtraction
// is sufficient and the post-subtraction remainder is < divisor.
// The quotient fits one limb. No ISA-dependent code is used.
namespace {

std::pair<std::uint64_t,std::uint64_t> divide_wide_portable(
        std::uint64_t high,std::uint64_t low,std::uint64_t divisor) noexcept
{
    std::uint64_t remainder=high;
    std::uint64_t quotient=0u;
    for (unsigned bit=64u;bit!=0u;--bit) {
        const std::uint64_t overflow=remainder>>63u;
        remainder=(remainder<<1u)|((low>>(bit-1u))&std::uint64_t{1});
        if (overflow!=0u || remainder>=divisor) {
            remainder-=divisor;
            quotient|=std::uint64_t{1}<<(bit-1u);
        }
    }
    return {quotient,remainder};
}

} // namespace

BigUnsigned::Limb BigUnsigned::divrem_limb_into(
        LimbView a,Limb divisor,MutableLimbView quotient)
{
    if (divisor==0u) {
        throw std::domain_error("BigUnsigned division by zero");
    }
    if (quotient.size()<a.size()) {
        throw std::length_error("BigUnsigned quotient output too small");
    }

    Limb remainder=0u;
    for (std::size_t i=a.size();i!=0u;--i) {
        const auto [digit,next]=divide_wide_portable(
            remainder,a[i-1u],divisor);
        quotient[i-1u]=digit;
        remainder=next;
    }
    std::fill(quotient.subspan(a.size()).begin(),quotient.end(),Limb{0});
    return remainder;
}

std::pair<BigUnsigned,BigUnsigned::Limb> divrem_limb(
        BigUnsigned const& a,BigUnsigned::Limb divisor)
{
    BigUnsigned::Scratch scratch;
    auto digits=scratch.acquire(a.size());
    const auto remainder=BigUnsigned::divrem_limb_into(
        a.limbs(),divisor,digits);
    return {BigUnsigned(digits),remainder};
}


std::size_t BigUnsigned::divrem_scratch_limbs(
        std::size_t dividend_limbs,std::size_t divisor_limbs)
{
    if (divisor_limbs<2u || dividend_limbs<divisor_limbs) {
        return 0u;
    }
    const auto max=std::numeric_limits<std::size_t>::max();
    if (divisor_limbs>(max-2u)/2u ||
            dividend_limbs>max-(2u*divisor_limbs+2u)) {
        throw std::length_error("BigUnsigned division scratch size overflow");
    }
    // Normalized dividend (n+1), divisor (m), digit product (m+1).
    return dividend_limbs+2u*divisor_limbs+2u;
}

namespace {

// Compare equal-width, little-endian unsigned sequences.
bool greater_fixed(BigUnsigned::LimbView a,BigUnsigned::LimbView b) noexcept
{
    for (std::size_t i=a.size();i!=0u;--i) {
        if (a[i-1u]!=b[i-1u]) {
            return a[i-1u]>b[i-1u];
        }
    }
    return false;
}

// Subtract a smaller disjoint unsigned view from a mutable destination.
void subtract_in_place(BigUnsigned::MutableLimbView destination,
                       BigUnsigned::LimbView source) noexcept
{
    BigUnsigned::Limb borrow=0u;
    for (std::size_t i=0u;i<destination.size();++i) {
        const auto left=destination[i];
        const auto right=i<source.size()?source[i]:BigUnsigned::Limb{0};
        const auto difference=left-right;
        destination[i]=difference-borrow;
        borrow=static_cast<BigUnsigned::Limb>(
            left<right || difference<borrow);
    }
}

} // namespace

void BigUnsigned::divrem_into(
        LimbView a,LimbView b,MutableLimbView quotient,
        MutableLimbView remainder,MutableLimbView scratch)
{
    if (b.empty()) {
        throw std::domain_error("BigUnsigned division by zero");
    }
    if (quotient.size()<a.size() ||
            remainder.size()<std::min(a.size(),b.size())) {
        throw std::length_error("BigUnsigned division output too small");
    }
    const auto required=divrem_scratch_limbs(a.size(),b.size());
    if (scratch.size()<required) {
        throw std::length_error("BigUnsigned division scratch too small");
    }

    std::fill(quotient.begin(),quotient.end(),Limb{0});
    std::fill(remainder.begin(),remainder.end(),Limb{0});
    if (a.size()<b.size()) {
        std::copy(a.begin(),a.end(),remainder.begin());
        return;
    }
    if (b.size()==1u) {
        remainder[0]=divrem_limb_into(a,b[0],quotient);
        return;
    }

    const std::size_t n=a.size();
    const std::size_t m=b.size();
    auto u=scratch.first(n+1u);
    auto v=scratch.subspan(n+1u,m);
    auto product=scratch.subspan(n+1u+m,m+1u);

    // Normalize the divisor's highest bit; keep the dividend carry.
    const unsigned shift=static_cast<unsigned>(std::countl_zero(b.back()));
    Limb carry=0u;
    for (std::size_t i=0u;i<n;++i) {
        const auto digit=a[i];
        u[i]=(digit<<shift)|carry;
        carry=shift==0u?0u:digit>>(64u-shift);
    }
    u[n]=carry;
    carry=0u;
    for (std::size_t i=0u;i<m;++i) {
        const auto digit=b[i];
        v[i]=(digit<<shift)|carry;
        carry=shift==0u?0u:digit>>(64u-shift);
    }

    // Radix-2^64 long division: the normalized top-two-limb estimate
    // bounds the quotient digit from above. Correct its full product
    // downward before subtracting; at most two corrections are needed.
    for (std::size_t counter=n-m+1u;counter!=0u;--counter) {
        const std::size_t j=counter-1u;
        const Limb high=u[j+m];
        const Limb low=u[j+m-1u];
        Limb digit=high>=v[m-1u] ?
            std::numeric_limits<Limb>::max() :
            divide_wide_portable(high,low,v[m-1u]).first;
        mul_into(v,LimbView(&digit,1u),product);
        const auto window=LimbView(u.data()+j,m+1u);
        while (greater_fixed(product,window)) {
            --digit;
            subtract_in_place(product,v);
        }
        subtract_in_place(u.subspan(j,m+1u),product);
        quotient[j]=digit;
    }

    // Unnormalize the remainder; it is strictly smaller than b.
    for (std::size_t i=0u;i<m;++i) {
        const Limb high=i+1u<m?u[i+1u]:Limb{0};
        remainder[i]=shift==0u?u[i]:
            (u[i]>>shift)|(high<<(64u-shift));
    }
}

std::pair<BigUnsigned,BigUnsigned> divrem(
        BigUnsigned const& a,BigUnsigned const& b)
{
    if (b.is_zero()) {
        throw std::domain_error("BigUnsigned division by zero");
    }
    BigUnsigned quotient;
    BigUnsigned remainder;
    quotient.resize_for_result(a.size());
    const auto rem_size=std::min(a.size(),b.size());
    remainder.resize_for_result(rem_size);
    BigUnsigned::Scratch workspace;
    const auto scratch=workspace.acquire(
        BigUnsigned::divrem_scratch_limbs(a.size(),b.size()));
    BigUnsigned::divrem_into(
        a.limbs(),b.limbs(),
        BigUnsigned::MutableLimbView(quotient.mutable_data(),a.size()),
        BigUnsigned::MutableLimbView(remainder.mutable_data(),rem_size),
        scratch);
    quotient.normalize();
    remainder.normalize();
    return {std::move(quotient),std::move(remainder)};
}

namespace {

// Compiler-supported wide accumulation, shared by GCC/Clang x86-64 and ARM64.
// The intermediate is bounded by (B-1)^2+2*(B-1)=B^2-1, B=2^64.
// Preserve the historical MSVC x64 native B3 fallback. On other toolchains
// without unsigned __int128, retain portable multiplication as a fallback.
// Inputs and output must be disjoint; all output limbs are overwritten.
void multiply_limbs_fused(BigUnsigned::LimbView a,
                          BigUnsigned::LimbView b,
                          BigUnsigned::MutableLimbView output) noexcept
{
#if defined(_MSC_VER) && defined(_M_X64) && defined(ARIADNE_BIG_UNSIGNED_X86_64)
    BigUnsignedDetail::mul_low_into_x86_64(a,b,output);
#elif defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    __extension__ using Wide = unsigned __int128;
    std::fill(output.begin(),output.end(),BigUnsigned::Limb{0});
    for (std::size_t i=0u; i<a.size() && i<output.size(); ++i) {
        BigUnsigned::Limb carry=0u;
        const std::size_t terms=std::min(b.size(),output.size()-i);
        for (std::size_t j=0u; j<terms; ++j) {
            const Wide total=static_cast<Wide>(a[i])*b[j]
                             +output[i+j]+carry;
            output[i+j]=static_cast<BigUnsigned::Limb>(total);
            carry=static_cast<BigUnsigned::Limb>(total>>64u);
        }
        if (i+terms<output.size()) {
            output[i+terms]=carry;
        }
    }
#else
    BigUnsigned::mul_low_into_portable(a,b,output);
#endif
}

} // namespace

bool BigUnsigned::has_native_multiply() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) || defined(ARIADNE_BIG_UNSIGNED_ARM64)
    return true;
#else
    return false;
#endif
}

bool BigUnsigned::has_arm64_multiply() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_ARM64)
    return true;
#else
    return false;
#endif
}

bool BigUnsigned::has_x86_64_carry_multiply() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) && defined(__x86_64__) \
    && (defined(__GNUC__) || defined(__clang__)) && !defined(_MSC_VER)
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
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED)
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
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) || defined(ARIADNE_BIG_UNSIGNED_ARM64)
    multiply_limbs_fused(a,b,output);
#else
    mul_low_into_portable(a,b,output);
#endif
}

void BigUnsigned::mul_low_into_arm64_carry(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_ARM64)
    BigUnsignedDetail::mul_low_into_carry_arm64(a,b,output);
#else
    mul_low_into_portable(a,b,output);
#endif
}

void BigUnsigned::mul_low_into_x86_64_carry(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    BigUnsignedDetail::mul_low_into_carry_x86_64(a,b,output);
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
#if defined(ARIADNE_BIG_UNSIGNED_ARM64)
    // Apple Silicon measurements and native GMP differential tests validate
    // the explicit carry-chain path; fused remains available for comparison.
    BigUnsignedDetail::mul_low_into_carry_arm64(a,b,output);
#elif defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED)
    // GCC x86-64 continues to select the shared unsigned __int128 kernel.
    multiply_limbs_fused(a,b,output);
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

void BigUnsigned::mul_into_arm64_carry(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_arm64_carry(a,b,output);
}

void BigUnsigned::mul_into_x86_64_carry(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_x86_64_carry(a,b,output);
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


// Portable, allocation-free Karatsuba implementation (formerly separate TU).
// Kept in the common arithmetic source alongside selected limb kernels.
namespace {

using Limb = BigUnsigned::Limb;
using Read = BigUnsigned::LimbView;
using Write = BigUnsigned::MutableLimbView;

// This is a recursion cutoff, not a proven performance crossover.
// Small cases use the current architecture-selected schoolbook kernel.
constexpr std::size_t karatsuba_default_cutoff=8u;
constexpr std::size_t karatsuba_minimum_cutoff=8u;

std::size_t checked_sum(std::size_t a,std::size_t b)
{
    if (a>std::numeric_limits<std::size_t>::max()-b) {
        throw std::length_error("BigUnsigned Karatsuba size overflow");
    }
    return a+b;
}

bool should_split(std::size_t a,std::size_t b,
                  std::size_t cutoff) noexcept
{
    const std::size_t half=std::max(a,b)/2u;
    return std::min(a,b)>cutoff
        && a>half && b>half;
}

// Exact recursion-frame bound: sumA, sumB, cross product, plus the
// maximum of the three sequential recursive calls' workspaces.
// No nested frame overlaps the live parent sums or cross product.
std::size_t karatsuba_workspace(
        std::size_t a,std::size_t b,std::size_t cutoff)
{
    if (!should_split(a,b,cutoff)) {
        return 0u;
    }

    const std::size_t half=std::max(a,b)/2u;
    const std::size_t a0=std::min(a,half);
    const std::size_t b0=std::min(b,half);
    const std::size_t a1=a-a0;
    const std::size_t b1=b-b0;
    const std::size_t sa=checked_sum(std::max(a0,a1),1u);
    const std::size_t sb=checked_sum(std::max(b0,b1),1u);
    const std::size_t pair=checked_sum(sa,sb);
    const std::size_t frame=checked_sum(pair,pair);

    const std::size_t nested=std::max({
        karatsuba_workspace(a0,b0,cutoff),
        karatsuba_workspace(a1,b1,cutoff),
        karatsuba_workspace(sa,sb,cutoff)
    });
    return checked_sum(frame,nested);
}

// An extra high zero limb makes addition exact even when both halves
// are saturated. Arguments must not overlap the destination.
void sum_halves(Read a,Read b,Write output) noexcept
{
    Limb carry=0u;
    for (std::size_t i=0u;i<output.size();++i) {
        const Limb left=i<a.size() ? a[i] : Limb{0};
        const Limb right=i<b.size() ? b[i] : Limb{0};
        const Limb sum=left+right;
        const Limb value=sum+carry;
        carry=static_cast<Limb>(sum<left || value<sum);
        output[i]=value;
    }
}

// The caller has proved lhs >= rhs as unsigned integers. Subtracting
// from lhs itself avoids an extra allocated intermediate buffer.
void subtract_in_place(Write lhs,Read rhs)
{
    Limb borrow=0u;
    for (std::size_t i=0u;i<lhs.size();++i) {
        const Limb right=i<rhs.size() ? rhs[i] : Limb{0};
        const Limb left=lhs[i];
        const Limb difference=left-right;
        const Limb value=difference-borrow;
        borrow=static_cast<Limb>(left<right || difference<borrow);
        lhs[i]=value;
    }
    if (borrow!=0u) {
        throw std::logic_error("BigUnsigned Karatsuba negative cross term");
    }
}

// Accumulate the shifted cross term, touching only limbs affected by
// term or a propagated carry. In the no-carry case, the remaining output
// limbs are already correct and need not be read or rewritten.
// A widened sum-of-halves temporary can extend beyond the destination:
// every discarded word must still be checked for zero.
void add_shifted(Write output,Read term)
{
    const std::size_t common=std::min(output.size(),term.size());
    Limb carry=0u;
    std::size_t i=0u;
    for (;i<common;++i) {
        const Limb left=output[i];
        const Limb right=term[i];
        const Limb sum=left+right;
        const Limb value=sum+carry;
        carry=static_cast<Limb>(sum<left || value<sum);
        output[i]=value;
    }
    while (carry!=0u && i<output.size()) {
        const Limb left=output[i];
        output[i]=left+carry;
        carry=static_cast<Limb>(output[i]<left);
        ++i;
    }
    if (carry!=0u) {
        throw std::logic_error("BigUnsigned Karatsuba output overflow");
    }
    for (std::size_t j=output.size();j<term.size();++j) {
        if (term[j]!=0u) {
            throw std::logic_error("BigUnsigned Karatsuba nonzero high term");
        }
    }
}

// output.size()==a.size()+b.size(); scratch size is precomputed.
// z0 and z2 occupy disjoint output slices, so no extra result buffer
// is required for either of them. Only sums and cross term use scratch.
void karatsuba_recursive(
        Read a,Read b,Write output,Write scratch,
        std::size_t cutoff)
{
    if (!should_split(a.size(),b.size(),cutoff)) {
        BigUnsigned::mul_low_into(a,b,output);
        return;
    }

    // should_split() ensures both operands extend beyond half. Thus
    // a0.size()==b0.size()==half, z0 occupies [0,2*half) and z2
    // occupies [2*half,a.size()+b.size()). Their recursive base cases
    // overwrite every limb; clearing the entire parent is redundant.
    const std::size_t half=std::max(a.size(),b.size())/2u;
    const Read a0=a.first(std::min(a.size(),half));
    const Read b0=b.first(std::min(b.size(),half));
    const Read a1=a.subspan(a0.size());
    const Read b1=b.subspan(b0.size());

    karatsuba_recursive(
        a0,b0,output.first(a0.size()+b0.size()),scratch,cutoff);
    karatsuba_recursive(
        a1,b1,output.subspan(2u*half,a1.size()+b1.size()),
        scratch,cutoff);

    const std::size_t sa=std::max(a0.size(),a1.size())+1u;
    const std::size_t sb=std::max(b0.size(),b1.size())+1u;
    const std::size_t pair=sa+sb;

    Write sum_a=scratch.first(sa);
    Write sum_b=scratch.subspan(sa,sb);
    Write cross=scratch.subspan(pair,pair);
    Write nested=scratch.subspan(2u*pair);

    sum_halves(a0,a1,sum_a);
    sum_halves(b0,b1,sum_b);
    karatsuba_recursive(sum_a,sum_b,cross,nested,cutoff);

    // (a0+a1)(b0+b1)-a0*b0-a1*b1 = a0*b1+a1*b0.
    subtract_in_place(cross,output.first(a0.size()+b0.size()));
    subtract_in_place(
        cross,output.subspan(2u*half,a1.size()+b1.size()));
    add_shifted(output.subspan(half),cross);
}

} // namespace

std::size_t BigUnsigned::karatsuba_scratch_limbs(
        std::size_t a_limbs,std::size_t b_limbs)
{
    return karatsuba_scratch_limbs(
        a_limbs,b_limbs,karatsuba_default_cutoff);
}

std::size_t BigUnsigned::karatsuba_scratch_limbs(
        std::size_t a_limbs,std::size_t b_limbs,std::size_t cutoff)
{
    if (cutoff<karatsuba_minimum_cutoff) {
        throw std::invalid_argument("BigUnsigned Karatsuba cutoff must be >= 8");
    }
    (void)checked_sum(a_limbs,b_limbs);
    return karatsuba_workspace(a_limbs,b_limbs,cutoff);
}

void BigUnsigned::mul_into_karatsuba(
        LimbView a,LimbView b,
        MutableLimbView output,MutableLimbView scratch)
{
    mul_into_karatsuba(
        a,b,output,scratch,karatsuba_default_cutoff);
}

void BigUnsigned::mul_into_karatsuba(
        LimbView a,LimbView b,
        MutableLimbView output,MutableLimbView scratch,
        std::size_t cutoff)
{
    if (cutoff<karatsuba_minimum_cutoff) {
        throw std::invalid_argument("BigUnsigned Karatsuba cutoff must be >= 8");
    }
    if (a.empty() || b.empty()) {
        std::fill(output.begin(),output.end(),Limb{0});
        return;
    }

    const std::size_t count=checked_sum(a.size(),b.size());
    if (output.size()<count) {
        throw std::length_error("BigUnsigned Karatsuba output too small");
    }
    const std::size_t required=karatsuba_workspace(a.size(),b.size(),cutoff);
    if (scratch.size()<required) {
        throw std::length_error("BigUnsigned Karatsuba scratch too small");
    }

    // Recursion writes all count result limbs, including its base cases.
    // Only extra destination words need explicit clearing.
    std::fill(output.subspan(count).begin(),output.end(),Limb{0});
    karatsuba_recursive(
        a,b,output.first(count),scratch.first(required),cutoff);
}

} // namespace Ariadne
