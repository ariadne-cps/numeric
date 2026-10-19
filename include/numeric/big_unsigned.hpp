/***************************************************************************
 *            numeric/big_unsigned.hpp
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

/*! \file numeric/big_unsigned.hpp
 *  \brief Shared GMP-independent unsigned multiprecision arithmetic.
 */

#ifndef ARIADNE_BIG_UNSIGNED_HPP
#define ARIADNE_BIG_UNSIGNED_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Canonical nonnegative integer with little-endian 64-bit limbs.
//! \sa Integer, Rational, FloatAP
class BigUnsigned
{
  public:
    using Limb = std::uint64_t;
    static constexpr std::size_t inline_capacity = 4u;
    using LimbView = std::span<Limb const>;
    using MutableLimbView = std::span<Limb>;

    //! Reusable, independently owned, zeroed temporary limb storage.
    //! Every acquire invalidates the previous view; not thread-shared.
    class Scratch {
        static constexpr std::size_t inline_capacity = 16u;
        std::array<Limb,inline_capacity> _inline_limbs{};
        std::vector<Limb> _heap;

      public:
        MutableLimbView acquire(std::size_t count);
        std::size_t capacity() const noexcept;
    };

  private:
    std::array<Limb,inline_capacity> _inline_limbs{};
    std::vector<Limb> _heap;
    std::size_t _size = 0u;

    Limb* mutable_data() noexcept;
    void resize_for_result(std::size_t limbs);
    void normalize() noexcept;
    void set_limbs(std::span<Limb const> digits);

  public:
    BigUnsigned() = default;
    explicit BigUnsigned(Limb value) noexcept;
    explicit BigUnsigned(std::span<Limb const> digits);
    BigUnsigned(BigUnsigned const&) = default;
    BigUnsigned(BigUnsigned&& other) noexcept;
    BigUnsigned& operator=(BigUnsigned const&) = default;
    BigUnsigned& operator=(BigUnsigned&& other) noexcept;
    ~BigUnsigned() = default;

    std::size_t size() const noexcept { return _size; }
    bool is_zero() const noexcept { return _size == 0u; }
    LimbView limbs() const noexcept;
    //! Replace this magnitude with the canonical value of a limb view.
    //! Self-referential views are supported, including subspans.
    void assign(LimbView digits);

    //! Allocation-free exact limb kernels; inputs and output must be disjoint.
    //! Result is written to the first max(a.size(),b.size()) limbs;
    //! higher destination limbs are zero; the return value is carry (0 or 1).
    static Limb add_into(LimbView a,LimbView b,MutableLimbView output);
    //! Unsigned subtraction modulo 2^(64*output.size()) when output
    //! is large enough for all input limbs. Returns borrow (0 or 1).
    //! Inputs and output must be disjoint.
    static Limb sub_into(LimbView a,LimbView b,MutableLimbView output);
    //! Divide a nonnegative limb view by a nonzero 64-bit limb.
    //! The quotient uses caller-owned storage of at least a.size()
    //! limbs; any extra output limbs are zeroed. Returns the exact
    //! remainder (< divisor). Source and destination must not overlap.
    //! Throws on a zero divisor or insufficient destination size.
    //! This operation has no ISA-specific implementation.
    static Limb divrem_limb_into(
        LimbView a,Limb divisor,MutableLimbView quotient);
    //! Exact scratch for normalized source views of the supplied lengths.
    //! No scratch is needed for one-limb divisors or shorter dividends.
    //! Throws length_error on size overflow.
    static std::size_t divrem_scratch_limbs(
        std::size_t dividend_limbs,std::size_t divisor_limbs);
    //! Exact unsigned quotient and remainder in caller-owned storage.
    //! Input views must be canonical: no high zero limbs, except empty zero.
    //! Quotient requires a.size() limbs; remainder min(a.size(),b.size()).
    //! Extra destination limbs are set to zero. All five views must be
    //! pairwise disjoint; no allocation or mutation before size checks.
    //! Zero divisor raises domain_error, undersized storage length_error.
    static void divrem_into(
        LimbView a,LimbView b,MutableLimbView quotient,
        MutableLimbView remainder,MutableLimbView scratch);
    //! Product modulo 2^(64*output.size()); output is fully overwritten.
    //! Inputs and output must be disjoint. No intermediate allocation.
    //! Uses the selected architecture kernel when it is available.
    static void mul_low_into(LimbView a,LimbView b,MutableLimbView output);
    //! Identical modular product using only portable C++20 arithmetic.
    //! Implemented in big_unsigned_fallback.cpp, built on every target
    //! for independent differential tests and fair benchmarks.
    static void mul_low_into_portable(LimbView a,LimbView b,MutableLimbView output);
    //! Exact full product in caller storage; output must have at least
    //! a.size()+b.size() limbs for nonempty operands. Extra limbs are zeroed.
    //! Both input views must be disjoint from output. No allocations.
    static void mul_into(LimbView a,LimbView b,MutableLimbView output);
    //! Full product using the explicitly selected portable reference.
    static void mul_into_portable(LimbView a,LimbView b,MutableLimbView output);
    //! Experimental column-wise Comba kernel for benchmarking on x86-64.
    //! On other architectures, uses the portable schoolbook reference.
    //! Result is the product modulo 2^(64*output.size()) and completely
    //! overwrites output. Output must be disjoint from both inputs.
    //! Does not alter the architecture-selected production path.
    static void mul_low_into_comba(LimbView a,LimbView b,MutableLimbView output);
    //! Exact full Comba product into caller storage. Checks that output
    //! has a.size()+b.size() limbs for nonempty inputs; does not allocate.
    //! Overlap with either input is forbidden; extra limbs are zeroed.
    static void mul_into_comba(LimbView a,LimbView b,MutableLimbView output);
    //! Shared fused 128-bit accumulation schoolbook for GCC/Clang x86-64
    //! and ARM64. MSVC x86-64 retains the native B3 schoolbook fallback;
    //! other architectures retain the portable C++20 reference.
    //! Output is disjoint from inputs and fully overwritten; computes the
    //! exact product modulo 2^(64*output.size()) without allocation.
    static void mul_low_into_fused(LimbView a,LimbView b,MutableLimbView output);
    //! Explicit full fused product; same size and aliasing contracts
    //! as mul_into. Does not affect the production kernel selector.
    static void mul_into_fused(LimbView a,LimbView b,MutableLimbView output);
    //! Explicit selected ARM64 MUL/UMULH/ADDS/ADC carry-chain
    //! schoolbook. Elsewhere uses the portable reference; this entry
    //! point remains available for benchmarking and regression.
    static void mul_low_into_arm64_carry(
        LimbView a,LimbView b,MutableLimbView output);
    //! Exact full-product ARM64 carry-chain implementation.
    static void mul_into_arm64_carry(
        LimbView a,LimbView b,MutableLimbView output);
    //! Experimental x86-64 MULQ/ADDQ/ADCQ carry-chain row schoolbook.
    //! Native GNU/Clang assembly, MSVC B3 fallback, portable elsewhere.
    //! Fully overwrites a disjoint output, without scratch allocation.
    //! Does not change the production dispatch.
    static void mul_low_into_x86_64_carry(
        LimbView a,LimbView b,MutableLimbView output);
    //! Exact full-product counterpart with mul_into's output size contract.
    static void mul_into_x86_64_carry(
        LimbView a,LimbView b,MutableLimbView output);
    //! Force the B3 x86-64 row-schoolbook reference chosen before B6.
    //! On unsupported architectures, uses the C++20 portable reference.
    //! Product is exact modulo 2^(64*output.size()); no overlap allowed.
    static void mul_low_into_schoolbook(
        LimbView a,LimbView b,MutableLimbView output);
    //! Force B3 full-product selection using caller-owned output storage.
    static void mul_into_schoolbook(
        LimbView a,LimbView b,MutableLimbView output);
    //! B8 experimental Karatsuba multiplication (full product only).
    //! Caller supplies output of at least a.size()+b.size() limbs and
    //! disjoint scratch obtained from karatsuba_scratch_limbs(). Both
    //! output and scratch must be disjoint from each input and each other.
    //! Output is fully overwritten, with extra limbs zeroed. Neither
    //! this kernel nor its recursive calls allocate. Does NOT change the
    //! selected production multiply. The base case uses selected schoolbook.
    static void mul_into_karatsuba(
        LimbView a,LimbView b,MutableLimbView output,
        MutableLimbView scratch);
    //! Required B8 scratch size in 64-bit limbs for given input lengths.
    //! Returns zero for base cases; may throw on size arithmetic overflow.
    static std::size_t karatsuba_scratch_limbs(
        std::size_t a_limbs,std::size_t b_limbs);
    //! Explicit B8.1 variant with a caller-selected recursion cutoff
    //! (minimum 8 limbs); below cutoff uses current selected schoolbook.
    //! Scratch must be computed using the SAME cutoff. No global state.
    static void mul_into_karatsuba(
        LimbView a,LimbView b,MutableLimbView output,
        MutableLimbView scratch,std::size_t cutoff);
    //! Scratch size for a caller-selected Karatsuba recursion cutoff.
    static std::size_t karatsuba_scratch_limbs(
        std::size_t a_limbs,std::size_t b_limbs,std::size_t cutoff);
    //! Whether a supported x86-64 or ARM64 native kernel was built.
    //! Native ARM64 selects the specialized carry-chain kernel by default.
    static bool has_native_multiply() noexcept;
    //! Whether the ARM64 native fused/carry-chain kernels were built.
    static bool has_arm64_multiply() noexcept;
    //! Whether GNU/Clang x86-64 inline-assembly carry experiment exists.
    static bool has_x86_64_carry_multiply() noexcept;
    //! Whether the x86-64 native schoolbook/Comba kernels were built.
    static bool has_x86_64_multiply() noexcept;
    //! Whether production selects the common unsigned __int128 fused kernel.
    //! True only for GCC x86-64 with its independent selector enabled;
    //! native ARM64 instead selects the explicit carry-chain kernel.
    static bool selected_multiply_uses_fused() noexcept;
    Limb limb(std::size_t index) const noexcept;
    std::size_t bit_length() const noexcept;
    void swap(BigUnsigned& other) noexcept;

    //! Compare exact unsigned values: -1, 0 or +1.
    friend int compare(BigUnsigned const& a, BigUnsigned const& b) noexcept;

    //! Exact addition; inputs are never modified.
    friend BigUnsigned add(BigUnsigned const& a, BigUnsigned const& b);
    //! Exact subtraction; throws std::underflow_error if a < b.
    friend BigUnsigned sub(BigUnsigned const& a, BigUnsigned const& b);
    //! Divide by one limb; return normalized quotient and exact remainder.
    //! Throws std::domain_error for a zero divisor.
    friend std::pair<BigUnsigned,Limb> divrem_limb(
        BigUnsigned const& a,Limb divisor);
    //! Euclidean unsigned quotient and remainder: a = b*q+r, 0 <= r < b.
    //! Returns canonical magnitudes, throws domain_error for b=0.
    friend std::pair<BigUnsigned,BigUnsigned> divrem(
        BigUnsigned const& a,BigUnsigned const& b);
    //! Exact product of two unsigned magnitudes.
    friend BigUnsigned mul(BigUnsigned const& a, BigUnsigned const& b);

    //! Product modulo 2^(64*limb_count); zero limbs returns zero.
    friend BigUnsigned mul_low(BigUnsigned const& a, BigUnsigned const& b,
                               std::size_t limb_count);

    //! Exact square (currently uses the full multiplication kernel).
    friend BigUnsigned square(BigUnsigned const& a);

    //! Exact multiplication by 2^bits (no fixed precision limit).
    friend BigUnsigned shift_left(BigUnsigned const& a, std::size_t bits);
    //! Floor division by 2^bits.
    friend BigUnsigned shift_right(BigUnsigned const& a, std::size_t bits);
};

} // namespace Ariadne

#endif
