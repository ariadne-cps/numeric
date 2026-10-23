/***************************************************************************
 *            integer_gmp_bridge.hpp
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

#ifndef ARIADNE_INTEGER_GMP_BRIDGE_HPP
#define ARIADNE_INTEGER_GMP_BRIDGE_HPP

// Temporary implementation-boundary conversions. Integer itself stores
// only a native BigUnsigned magnitude and sign, never mpz_t.
#include "numeric/integer.hpp"

#include <gmp.h>

#include <array>
#include <cstdint>
#include <vector>

namespace Ariadne {

inline Integer integer_from_gmp(mpz_srcptr source)
{
    if (mpz_sgn(source)==0) {
        return Integer{};
    }
    // Four limbs fit inside BigUnsigned itself. Use the same inline
    // capacity for the GMP export buffer, avoiding a heap allocation
    // on the overwhelmingly common small-precision boundary path.
    // mpz_sizeinbase(x,2) is exact for nonzero x.
    const std::size_t bit_count=mpz_sizeinbase(source,2u);
    const std::size_t capacity=1u+(bit_count-1u)/64u;
    std::array<std::uint64_t,BigUnsigned::inline_capacity> inline_words{};
    std::vector<std::uint64_t> large_words;
    std::uint64_t* words=inline_words.data();
    if (capacity>inline_words.size()) {
        large_words.resize(capacity);
        words=large_words.data();
    }
    std::size_t used=0u;
    mpz_export(words,&used,-1,sizeof(std::uint64_t),0,0,source);
    return Integer::from_magnitude(
        BigUnsigned(std::span(words,used)),mpz_sgn(source)<0);
}

inline void integer_to_gmp(mpz_ptr destination,Integer const& source)
{
    const auto limbs=source.magnitude().limbs();
    mpz_import(destination,limbs.size(),-1,
               sizeof(BigUnsigned::Limb),0,0,limbs.data());
    if (source.negative()) {
        mpz_neg(destination,destination);
    }
}

// Used only at interfaces where another numeric type still stores GMP.
// This is not a second Integer implementation or permanent storage.
class GmpIntegerBoundary
{
    mpz_t _value;

  public:
    GmpIntegerBoundary() { mpz_init(_value); }
    explicit GmpIntegerBoundary(Integer const& source) : GmpIntegerBoundary() {
        integer_to_gmp(_value,source);
    }
    GmpIntegerBoundary(GmpIntegerBoundary const&) = delete;
    GmpIntegerBoundary& operator=(GmpIntegerBoundary const&) = delete;
    ~GmpIntegerBoundary() { mpz_clear(_value); }
    mpz_ptr mutable_value() noexcept { return _value; }
    mpz_srcptr value() const noexcept { return _value; }
};

} // namespace Ariadne

#endif
