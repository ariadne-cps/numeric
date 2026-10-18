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
    //! Product modulo 2^(64*output.size()); output is fully overwritten.
    //! Inputs and output must be disjoint. No intermediate allocation.
    static void mul_low_into(LimbView a,LimbView b,MutableLimbView output);
    Limb limb(std::size_t index) const noexcept;
    std::size_t bit_length() const noexcept;
    void swap(BigUnsigned& other) noexcept;

    //! Compare exact unsigned values: -1, 0 or +1.
    friend int compare(BigUnsigned const& a, BigUnsigned const& b) noexcept;

    //! Exact addition; inputs are never modified.
    friend BigUnsigned add(BigUnsigned const& a, BigUnsigned const& b);
    //! Exact subtraction; throws std::underflow_error if a < b.
    friend BigUnsigned sub(BigUnsigned const& a, BigUnsigned const& b);
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
