/***************************************************************************
 *            numeric/signed_integer.hpp
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

/*! \file numeric/signed_integer.hpp
 *  \brief GMP-independent signed magnitude for native Integer migration.
 */

#ifndef ARIADNE_SIGNED_INTEGER_HPP
#define ARIADNE_SIGNED_INTEGER_HPP

#include "numeric/big_unsigned.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Canonical signed magnitude based on BigUnsigned.
//! \sa Integer, BigUnsigned
//! Zero always has a nonnegative sign. No GMP storage or arithmetic.
class SignedInteger
{
    BigUnsigned _magnitude;
    bool _negative = false;

    SignedInteger(BigUnsigned magnitude,bool negative) noexcept;

  public:
    SignedInteger() = default;
    explicit SignedInteger(std::int64_t value) noexcept;
    explicit SignedInteger(std::uint64_t value) noexcept;
    //! Accept optional '+' or '-' followed by nonempty ASCII decimal digits.
    //! Whitespace and other characters raise std::invalid_argument.
    explicit SignedInteger(std::string_view decimal);

    static SignedInteger from_magnitude(BigUnsigned magnitude,bool negative)
        noexcept;

    BigUnsigned const& magnitude() const noexcept { return _magnitude; }
    bool is_zero() const noexcept { return _magnitude.is_zero(); }
    bool is_negative() const noexcept { return _negative; }
    //! Returns -1, 0 or 1.
    int sign() const noexcept;
    std::string decimal() const;

    //! Checked conversions. Out-of-range and negative-to-unsigned throw.
    std::uint64_t get_u64() const;
    std::int64_t get_i64() const;

    friend int compare(SignedInteger const& a,SignedInteger const& b) noexcept;
    friend SignedInteger neg(SignedInteger const& a);
    friend SignedInteger add(SignedInteger const& a,SignedInteger const& b);
    friend SignedInteger sub(SignedInteger const& a,SignedInteger const& b);
    friend SignedInteger mul(SignedInteger const& a,SignedInteger const& b);
    friend SignedInteger square(SignedInteger const& a);
    friend SignedInteger pow(SignedInteger const& a,std::uint64_t exponent);
    //! Division truncates toward zero. Remainder has the dividend's sign.
    //! For b!=0, a=b*q+r and abs(r)<abs(b). Throws on division by zero.
    friend std::pair<SignedInteger,SignedInteger> divrem(
        SignedInteger const& a,SignedInteger const& b);
    friend SignedInteger shift_left(SignedInteger const& a,std::size_t bits);
    //! Signed right shift truncates toward zero (not floor division).
    friend SignedInteger shift_right(SignedInteger const& a,std::size_t bits);
};

} // namespace Ariadne

#endif
