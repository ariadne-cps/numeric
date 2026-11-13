/***************************************************************************
 *            numeric/decimal.hpp
 *
 *  Copyright  2014-20  Pieter Collins
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

/** @file numeric/decimal.hpp
 *  @brief ExactTag decimal numbers, useful for user input.
 */

#ifndef ARIADNE_DECIMAL_HPP
#define ARIADNE_DECIMAL_HPP

#include <iostream>
#include <string>

#include "numeric/integer.hpp"
#include "numeric/number.decl.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief A decimal number.
//! \sa Integer, Dyadic, Rational, Real
class Decimal {
  public:
    static const Integer _ten;
    Integer _p;
    Nat _q;

  public:
    typedef ExactTag Paradigm;
    //! \brief Default constructor creates the number 0 (zero).
    Decimal();
    //! \brief Construct the number p/10^q.
    explicit Decimal(Integer p, Nat q);
    /** @brief Disallow this construction to preserve the ownership contract. */
    explicit Decimal(Integer p, Int q) = delete;
    //! \brief Construct from a builtin integral type.
    template<BuiltinIntegral N>
    Decimal(N n);
    //! \brief Construct from a double-precision floating-point number representation.
    explicit Decimal(double d);
    //! \brief Construct from a string representation.
    explicit Decimal(const String&);
    //! \brief Convert from a dyadic.
    Decimal(const Dyadic&);
    //! \brief Convert to a rational number.
    operator Rational() const;
    //! \brief Convert to an generic number.
    operator ExactNumber() const;
    //! \brief A string literal
    String literal() const;

    //! \name Arithmetic operators
    //!@{
    friend Decimal operator+(const Decimal& d);                      //!< Unary plus \a +d.
    friend Decimal operator-(const Decimal& d);                      //!< Unary minus \a -d.
    friend Decimal operator+(const Decimal& d1, const Decimal& d2);  //!< Plus \a d1+d2.
    friend Decimal operator-(const Decimal& d1, const Decimal& d2);  //!< Minus \a d1-d2.
    friend Decimal operator*(const Decimal& d1, const Decimal& d2);  //!< Times \a d1*d2.
    friend Rational operator/(const Decimal& d1, const Decimal& d2); //!< Divides \a d1/d2, yielding a rational.
    friend Decimal& operator+=(Decimal& d1, const Decimal& d2);      //!< Inplace plus \a d1:=d1+d2.
    //!@}

    //!@{
    //! \name Comparison operators
    /** @brief Test equality of the operands. */
    friend Boolean operator==(const Decimal& d1, const Decimal& d2) { return cmp(d1, d2) == Comparison::EQUAL; } //!< <p/>
    /** @brief Test inequality of the operands. */
    friend Boolean operator!=(const Decimal& d1, const Decimal& d2) { return cmp(d1, d2) != Comparison::EQUAL; } //!< <p/>
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend Boolean operator<=(const Decimal& d1, const Decimal& d2) { return cmp(d1, d2) != Comparison::GREATER; } //!< <p/>
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend Boolean operator>=(const Decimal& d1, const Decimal& d2) { return cmp(d1, d2) != Comparison::LESS; } //!< <p/>
    /** @brief Test whether the first operand is less than the second. */
    friend Boolean operator<(const Decimal& d1, const Decimal& d2) { return cmp(d1, d2) == Comparison::LESS; } //!< <p/>
    /** @brief Test whether the first operand is greater than the second. */
    friend Boolean operator>(const Decimal& d1, const Decimal& d2) { return cmp(d1, d2) == Comparison::GREATER; } //!< <p/>
    //!@}

    //!@{
    //! \name Arithmetic operations
    friend Decimal nul(const Decimal& d); //!< Zero \a 0.
    friend Decimal sqr(const Decimal& d); //!< Square \a d<sup>2</sup>.
    friend Decimal hlf(const Decimal& d); //!< Half \a d/2.
    //!@}

    //!@{
    //! \name Lattice operations
    friend Decimal abs(const Decimal& d);                     //!< Absolute value \a |d|.
    friend Decimal min(const Decimal& d1, const Decimal& d2); //!< Minimum \a d1∧d2.
    friend Decimal max(const Decimal& d1, const Decimal& d2); //!< Maximum \a d1∨d2.
    //!@}

    //!@{
    //! \name Comparison operations
    friend Comparison cmp(const Decimal& d1, const Decimal& d2); //!< Comparison of two decimal values.
    //!@}

    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const Decimal& d);
    //! \brief Construct from an integer literal.
    friend Decimal operator""_decimal(unsigned long long int n);
    //! \brief Construct from a floating-point literal.
    friend Decimal operator""_decimal(long double dbl);
    //! \brief Construct from a string literal.
    friend Decimal operator""_decimal(const char* str, std::size_t);

    //! \brief Shorthand for operator""_decimal.
    friend Decimal operator""_dec(long double dbl);

    //! \brief Alternative for operator""_decimal for use in Python interface.
    friend Decimal dec_(long double dbl);

    /** @brief Put the value in canonical form. */
    void canonicalize();
};
/** @brief Construct the typed value from the literal. */
Decimal operator""_dec(unsigned long long int n);
/** @brief Construct the typed value from the literal. */
Decimal operator""_dec(long double dbl);
/** @brief Construct the typed value from the literal. */
Decimal operator""_dec(const char* str, std::size_t);
/** @brief Construct the typed value from the literal. */
Decimal operator""_decimal(unsigned long long int n);
/** @brief Construct the typed value from the literal. */
Decimal operator""_decimal(long double dbl);
/** @brief Construct the typed value from the literal. */
Decimal operator""_decimal(const char* str, std::size_t);

template<BuiltinIntegral N>
inline Decimal::Decimal(N n) : Decimal(n, 0u) {}

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<>
class Positive<Decimal> : public Decimal {
  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : Decimal() {}
    /** @brief Construct a \c Positive from m. */
    template<BuiltinUnsignedIntegral M>
    Positive(M m) : Decimal(m) {}
    /** @brief Disallow this construction to preserve the ownership contract. */
    Positive(int n) = delete;
    /** @brief Construct a \c Positive from w. */
    explicit Positive(const Decimal& w) : Decimal(w) { ARIADNE_ASSERT(w >= 0); }
};
/** @brief Convert to the positive representation. */
inline Positive<Decimal> cast_positive(const Decimal& w) { return Positive<Decimal>(w); }
using PositiveDecimal = Positive<Decimal>;

/** @brief Validated numeric enclosure with lower and upper endpoints. */
template<>
class Bounds<Decimal> {
    Decimal _l, _u;

  public:
    /** @brief Construct a \c Bounds from w. */
    Bounds(Decimal w) : _l(w), _u(w) {}
    /** @brief Construct a \c Bounds from x. */
    template<ConvertibleTo<Decimal> X>
    Bounds(const X& x)
        : Bounds(Decimal(x)) {}
    /** @brief Construct a \c Bounds from l, u. */
    Bounds(Decimal l, Decimal u) : _l(l), _u(u) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<class X>
        requires Constructible<Decimal, X>
    Bounds(const Bounds<X>& x)
        : Bounds(Decimal(x.lower_raw()), Decimal(x.upper_raw())) {}
    /** @brief Convert to validated number. */
    operator ValidatedNumber() const;
    /** @brief Return the lower endpoint. */
    Decimal lower() const { return _l; }
    /** @brief Return the upper endpoint. */
    Decimal upper() const { return _u; }
    /** @brief Return the underlying lower-endpoint representation. */
    Decimal lower_raw() const { return _l; }
    /** @brief Return the underlying upper-endpoint representation. */
    Decimal upper_raw() const { return _u; }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, Bounds<Decimal> y) { return os << "[" << y._l << ":" << y._u << "]"; }
};

} // namespace Ariadne

#endif
