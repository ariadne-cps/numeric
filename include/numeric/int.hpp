/***************************************************************************
 *            numeric/int.hpp
 *
 *  Copyright  2013-22  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @file numeric/int.hpp
 *  @brief Fixed-width integral classes
 */

#ifndef ARIADNE_INT_HPP
#define ARIADNE_INT_HPP

#include <limits>

#include "numeric/concepts.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

/************  Ints ********************************************************/

/** @brief Checked unsigned 32-bit integer value. */
class Nat32 {
    uint32_t _m;
    /** @brief Convert to an unsigned integer after checking the representable range. */
    static uint32_t _checked_unsigned(uint64_t m) {
        ARIADNE_PRECONDITION(m <= std::numeric_limits<uint32_t>::max());
        return static_cast<uint32_t>(m);
    }
    /** @brief Convert to a signed integer after checking the representable range. */
    static uint32_t _checked_signed(int64_t n) {
        ARIADNE_PRECONDITION(static_cast<uint64_t>(n) <= std::numeric_limits<uint32_t>::max());
        return static_cast<uint32_t>(n);
    }

  public:
    /** @brief Construct a \c Nat32 with the specified characteristics. */
    Nat32() : _m(0u) {}
    /** @brief Construct a \c Nat32 from m. */
    template<BuiltinUnsignedIntegral M>
    Nat32(M m) : _m(_checked_unsigned(static_cast<uint64_t>(m))) {}
    /** @brief Construct a \c Nat32 from n. */
    template<BuiltinSignedIntegral N>
    Nat32(N n) : _m(_checked_signed(static_cast<int64_t>(n))) {}
    /** @brief Convert to a built-in unsigned long integer. */
    uint32_t get_ui() const { return _m; }
};

/** @brief Checked unsigned 64-bit integer value. */
class Nat64 {
    uint64_t _m;
    /** @brief Convert to a signed integer after checking the representable range. */
    static uint64_t _checked_signed(int64_t n) {
        ARIADNE_PRECONDITION(n >= 0);
        return static_cast<uint64_t>(n);
    }

  public:
    /** @brief Construct a \c Nat64 with the specified characteristics. */
    Nat64() : _m(0u) {}
    /** @brief Construct a \c Nat64 from m. */
    template<BuiltinUnsignedIntegral M>
    Nat64(M m) : _m(static_cast<uint64_t>(m)) {}
    /** @brief Construct a \c Nat64 from n. */
    template<BuiltinSignedIntegral N>
    Nat64(N n) : _m(_checked_signed(static_cast<int64_t>(n))) {}
    /** @brief Convert to a built-in unsigned long integer. */
    uint64_t get_ui() const { return _m; }
};

/** @brief Checked signed 32-bit integer value. */
class Int32 {
    int32_t _n;
    /** @brief Convert to an unsigned integer after checking the representable range. */
    static int32_t _checked_unsigned(uint64_t m) {
        ARIADNE_PRECONDITION(m <= static_cast<uint64_t>(std::numeric_limits<int32_t>::max()));
        return static_cast<int32_t>(m);
    }
    /** @brief Convert to a signed integer after checking the representable range. */
    static int32_t _checked_signed(int64_t n) {
        ARIADNE_PRECONDITION(n >= std::numeric_limits<int32_t>::min());
        ARIADNE_PRECONDITION(n <= std::numeric_limits<int32_t>::max());
        return static_cast<int32_t>(n);
    }

  public:
    /** @brief Construct a \c Int32 with the specified characteristics. */
    Int32() : _n(0) {}
    /** @brief Construct a \c Int32 from m. */
    template<BuiltinUnsignedIntegral M>
    Int32(M m) : _n(_checked_unsigned(static_cast<uint64_t>(m))) {}
    /** @brief Construct a \c Int32 from n. */
    template<BuiltinSignedIntegral N>
    Int32(N n) : _n(_checked_signed(static_cast<int64_t>(n))) {}
    /** @brief Convert to a built-in signed long integer. */
    int32_t get_si() const { return _n; }
};

/** @brief Checked signed 64-bit integer value. */
class Int64 {
    int64_t _n;
    /** @brief Convert to an unsigned integer after checking the representable range. */
    static int64_t _checked_unsigned(uint64_t m) {
        ARIADNE_PRECONDITION(m <= static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));
        return static_cast<int64_t>(m);
    }

  public:
    /** @brief Construct a \c Int64 with the specified characteristics. */
    Int64() : _n(0) {}
    /** @brief Construct a \c Int64 from m. */
    template<BuiltinUnsignedIntegral M>
    Int64(M m) : _n(_checked_unsigned(static_cast<uint64_t>(m))) {}
    /** @brief Construct a \c Int64 from n. */
    template<BuiltinSignedIntegral N>
    Int64(N n) : _n(static_cast<int64_t>(n)) {}
    /** @brief Convert to a built-in signed long integer. */
    int64_t get_si() const { return _n; }
};

/** @brief Return the factorial. */
uint32_t fac(uint8_t n);
/** @brief Return the factorial. */
uint16_t fac(uint16_t n);
/** @brief Return the factorial. */
uint32_t fac(uint32_t n);
/** @brief Return the factorial. */
uint64_t fac(uint64_t n);
/** @brief Return the binomial coefficient. */
uint32_t bin(uint8_t n, uint8_t k);
/** @brief Return the binomial coefficient. */
uint16_t bin(uint16_t n, uint16_t k);
/** @brief Return the binomial coefficient. */
uint32_t bin(uint32_t n, uint32_t k);
/** @brief Return the binomial coefficient. */
uint64_t bin(uint64_t n, uint64_t k);

} // namespace Ariadne

#endif
