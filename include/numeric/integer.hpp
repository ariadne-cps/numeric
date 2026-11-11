/***************************************************************************
 *            numeric/integer.hpp
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

/** @file numeric/integer.hpp
 *  @brief Declarations and utilities for integer in Ariadne Numeric.
 */

#ifndef ARIADNE_INTEGER_HPP
#define ARIADNE_INTEGER_HPP

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "foundation/logical.hpp"
#include "numeric/arithmetic.hpp"
#include "numeric/big_unsigned.hpp"
#include "numeric/int.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/sign.hpp"
#include "utility/macros.hpp"
#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

struct ExactTag;
class Integer;
template<>
struct IsNumber<Integer> : True {};

template<>
struct IsNumber<Nat> : True {};
template<>
struct IsNumber<Int> : True {};

//! \ingroup NumericModule
//! \brief Arbitrarily-sized integers.
//! \sa Natural, Dyadic, Decimal, Rational, Real
class Integer {
  private:
    BigUnsigned _magnitude;
    bool _negative = false;

    // Normalizes negative zero; used by the signed arithmetic functions.
    Integer(BigUnsigned&& magnitude, bool negative) noexcept
        : _magnitude(std::move(magnitude)),
          _negative(negative && !_magnitude.is_zero()) {
    }

  public:
    //! Construct a canonical signed magnitude without using GMP.
    static Integer from_magnitude(BigUnsigned magnitude, bool negative) noexcept {
        return Integer(std::move(magnitude), negative);
    }
    const BigUnsigned& magnitude() const noexcept { return _magnitude; }
    bool negative() const noexcept { return _negative; }

    typedef ExactTag Paradigm; //!< <p/>

  public:
    ~Integer() = default; //!< Destructor
    Integer();            //!< Default constructor yielding \c 0
    template<BuiltinUnsignedIntegral M>
    Integer(M m); //!< Convert from a builtin unsigned (positive) integral type
    template<BuiltinSignedIntegral N>
    Integer(N n);                      //!< Convert from a builtin signed (positive or negative) integral type
    explicit Integer(const String&);   //!< Construct from a string literal
    Integer(const Integer&) = default; //!< Copy constructor
    Integer(Integer&& z) noexcept      //!< Move constructor
        : _magnitude(std::move(z._magnitude)), _negative(z._negative) {
        z._negative = false;
    }
    Integer& operator=(const Integer&) = default; //!< %Assignment
    Integer& operator=(Integer&& z) noexcept {    //!< Move assignment
        if (this != &z) {
            _magnitude = std::move(z._magnitude);
            _negative = z._negative;
            z._negative = false;
        }
        return *this;
    }
    operator ExactNumber() const; //!< Convert to a generic exact number
    String literal() const;       //!< A string literal

    friend Integer operator""_z(unsigned long long int n); //!< An integer literal

    //!@{
    //! \name Arithmetic operators
    friend Integer& operator++(Integer& z);                                                //!< <p/>
    friend Integer& operator--(Integer& z);                                                //!< <p/>
    friend Integer operator+(const Integer& z) { return pos(z); }                          //!< <p/>
    friend Integer operator-(const Integer& z) { return neg(z); }                          //!< <p/>
    friend Integer operator+(const Integer& z1, const Integer& z2) { return add(z1, z2); } //!< <p/>
    friend Integer operator-(const Integer& z1, const Integer& z2) { return sub(z1, z2); } //!< <p/>
    friend Integer operator*(const Integer& z1, const Integer& z2) { return mul(z1, z2); } //!< <p/>
    friend Rational operator/(const Integer& z1, const Integer& z2);                       //!< <p/>
    friend Integer operator%(const Integer& z1, const Integer& z2);                        //!< <p/>
    friend Integer& operator+=(Integer& z1, const Integer& z2);                            //!< <p/>
    friend Integer& operator-=(Integer& z1, const Integer& z2);                            //!< <p/>
    friend Integer& operator*=(Integer& z1, const Integer& z2);                            //!< <p/>
    //!@}

    //!@{
    //! \name Comparison operators
    friend Boolean operator==(const Integer& z1, const Integer& z2) { return eq(z1, z2); }  //!< <p/>
    friend Boolean operator!=(const Integer& z1, const Integer& z2) { return !eq(z1, z2); } //!< <p/>
    friend Boolean operator<=(const Integer& z1, const Integer& z2) { return !lt(z2, z1); } //!< <p/>
    friend Boolean operator>=(const Integer& z1, const Integer& z2) { return !lt(z1, z2); } //!< <p/>
    friend Boolean operator<(const Integer& z1, const Integer& z2) { return lt(z1, z2); }   //!< <p/>
    friend Boolean operator>(const Integer& z1, const Integer& z2) { return lt(z2, z1); }   //!< <p/>
    //!@}

    //!@{
    //! \name Arithmetic operations
    friend Int log2floor(const Natural& n); //!< \f$\lfloor\log_{2}n\rfloor\f$. Returns \c -1 if \a n is \c 0.

    friend Integer nul(const Integer& z);        //!< Zero \a 0.
    friend Integer pos(const Integer& z);        //!< Identity \a +z.
    friend Integer neg(const Integer& z);        //!< Negative \a -z.
    friend Natural sqr(const Integer& z);        //!< Square \a z<sup>2</sup>.
    friend Dyadic hlf(const Integer& z);         //!< Half \a z÷2.
    friend Dyadic shft(const Integer& z, Int n); //!< Bit-shift \a z×2<sup>n</sup>.
    friend Rational rec(const Integer& z);       //!< Reciprocal \a 1/z.

    friend Integer add(const Integer& z1, const Integer& z2);                    //!< \brief Sum \a z1+z2.
    friend Integer sub(const Integer& z1, const Integer& z2);                    //!< \brief Difference \a z1-z2.
    friend Integer mul(const Integer& z1, const Integer& z2);                    //!< \brief Product \a z1×z2.
    friend Rational div(const Integer& z1, const Integer& z2);                   //!< \brief Quotient \a z1÷z2.
    friend Integer fma(const Integer& z1, const Integer& z2, const Integer& z3); //!< \brief Fused multiply-and-add \a z1×z2+z3.
    friend Integer pow(const Integer& z, Nat m);                                 //!< \brief Power \a z<sup>m</sup>.
    //! Exact floor square root of a nonnegative Integer; throws domain_error
    //! for negative values. Zero maps to canonical nonnegative zero.
    friend Integer isqrt(const Integer& z);
    //! Exact root and remainder: z=root*root+remainder, 0<=remainder<=2*root.
    //! Requires z>=0; both returned Integer values are nonnegative.
    friend std::pair<Integer, Integer> sqrtrem(const Integer& z);
    friend Rational pow(const Integer& z, Int n); //!< \brief Power \a z<sup>n</sup>.

    friend Integer quot(const Integer& z1, const Integer& z2); //!< <p/>
    friend Integer rem(const Integer& z1, const Integer& z2);  //!< <p/>
    //!@}

    //!@{
    //! \name Lattice operations
    friend Natural abs(const Integer& z);                     //!< Absolute value \a |z|.
    friend Integer min(const Integer& z1, const Integer& z2); //!< Minimum \a z1∧z2.
    friend Integer max(const Integer& z1, const Integer& z2); //!< Maximum \a z1∨z2.
    //!@}

    //!@{
    //! \name Comparison operations
    friend Sign sgn(const Integer& z);                           //!< The sign of \a z.
    friend Comparison cmp(const Integer& w1, const Integer& w2); //!< Compares which of \a w1 and \a w2 is larger.
    friend Boolean eq(const Integer& w1, const Integer& w2);     //!< Tests if \a w1 is equal to \a w2.
    friend Boolean lt(const Integer& w1, const Integer& w2);     //!< Tests if \a w1 is less than \a w2.
    //!@}

    //!@{
    //! \name Special value tests
    friend Bool is_nan(const Integer& z);    //!< Tests whether \a z is NaN (not-a-number).
    friend Bool is_inf(const Integer& z);    //!< Tests whether \a z is ±∞.
    friend Bool is_finite(const Integer& z); //!< Tests whether \a z is finite.
    friend Bool is_zero(const Integer& z);   //!< Tests whether \a z is zero.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const Integer& z); //!< <p/>
                                                                         //!@}

  public:
    template<BuiltinIntegral N>
    N get() const;           //!< Convert to a builtin integral value of type \a N. Throws a \c std::runtime_error if the value does not fit in type \a N.
    long int get_si() const; //!< Checked conversion to host long; throws if out of range.

  private:

  public:
    Integer(Nat32 m);
    Integer(Int32 n);
    Integer(Nat64 m);
    Integer(Int64 n);
};

//! Signed wrapper for the shared GMP-free magnitude square root.
Integer isqrt(const Integer& z);
std::pair<Integer, Integer> sqrtrem(const Integer& z);

template<BuiltinUnsignedIntegral M>
inline Integer::Integer(M m) : Integer(Nat64(m)) {}
template<BuiltinSignedIntegral N>
inline Integer::Integer(N n) : Integer(Int64(n)) {}
Integer operator""_z(unsigned long long int n);
Integer operator""_z(const char* str, std::size_t);

template<BuiltinIntegral N>
inline N Integer::get() const {
    if (_magnitude.size() > 1u) {
        throw std::out_of_range("Integer does not fit destination type");
    }
    const std::uint64_t value = _magnitude.limb(0u);
    if (_negative) {
        constexpr std::uint64_t min_magnitude =
            static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) + 1u;
        if (value > min_magnitude) {
            throw std::out_of_range("Integer does not fit destination type");
        }
        if constexpr (BuiltinUnsignedIntegral<N>) {
            throw std::out_of_range("Integer does not fit destination type");
        } else {
            const std::int64_t signed_value = value == min_magnitude
                                                ? std::numeric_limits<std::int64_t>::min()
                                                : -static_cast<std::int64_t>(value);
            if constexpr (std::numeric_limits<N>::digits < 63) {
                if (!std::in_range<N>(signed_value)) {
                    throw std::out_of_range("Integer does not fit destination type");
                }
            }
            return static_cast<N>(signed_value);
        }
    }
    if constexpr (std::numeric_limits<N>::digits < 64) {
        if (!std::in_range<N>(value)) {
            throw std::out_of_range("Integer does not fit destination type");
        }
    }
    return static_cast<N>(value);
}

template<class R, class A>
R integer_cast(const A& a) {
    return cast_integer(a).template get<R>();
}

template<>
class Positive<Integer> : public Integer {
  public:
    Positive() : Integer() {}
    template<BuiltinUnsignedIntegral M>
    Positive(M m) : Integer(m) {}
    Positive(int n) = delete;
    explicit Positive(const Integer& z) : Integer(z) { ARIADNE_PRECONDITION(z >= 0); }
};

//! \brief A positive integer.
class Natural : public Positive<Integer> {
  public:
    Natural() : Positive<Integer>() {}
    template<BuiltinUnsignedIntegral M>
    Natural(M m) : Positive<Integer>(m) {}
    Natural(int n) = delete;
    explicit Natural(const Integer& z) : Positive<Integer>(z) {}
    friend Natural& operator++(Natural& n) {
        ++static_cast<Integer&>(n);
        return n;
    }
    friend Natural& operator+=(Natural& n1, const Natural& n2) {
        static_cast<Integer&>(n1) += n2;
        return n1;
    }
    friend Natural operator+(const Natural& n1, const Natural& n2) { return Natural(static_cast<const Integer&>(n1) + static_cast<const Integer&>(n2)); }
    friend Natural operator*(const Natural& n1, const Natural& n2) { return Natural(static_cast<const Integer&>(n1) * static_cast<const Integer&>(n2)); }
};

inline Natural cast_positive(const Integer& z) { return Natural(z); }

} // namespace Ariadne

#endif
