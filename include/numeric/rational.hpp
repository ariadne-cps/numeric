/***************************************************************************
 *            numeric/rational.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/** @file numeric/rational.hpp
 *  @brief Declarations and utilities for rational in Ariadne Numeric.
 */

#ifndef ARIADNE_RATIONAL_HPP
#define ARIADNE_RATIONAL_HPP

#include "numeric/arithmetic.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/integer.hpp"
#include "numeric/number.decl.hpp"
#include "utility/metaprogramming.hpp"
#include "utility/string.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

class Int64;
template<class X>
class FiniteOperations;
template<class X>
class ExtensionOperations;

enum class Comparison : ComparableEnumerationType;

/************ Rational *******************************************************/

//! \ingroup NumericModule
//! \brief %Rational numbers.
//! Finite values are reduced n/d with d>0 and canonical zero 0/1.
//! NaN and signed infinities have an explicit kind, independent of GMP.
//! \sa Integer, Dyadic, Decimal, Real
class Rational {
  private:
    enum class Kind : unsigned char { FINITE,
                                      POSITIVE_INFINITY,
                                      NEGATIVE_INFINITY,
                                      NAN_VALUE };
    Integer _numerator;
    BigUnsigned _denominator{1u};
    Kind _kind = Kind::FINITE;

    // Finite values are always reduced, with positive denominator and 0/1.
    // Only arithmetic which has proved coprimality may use this constructor.
    /** @brief Helper selecting a reduced numeric representation. */
    struct Reduced {};
    /** @brief Construct a \c Rational from numerator, denominator. */
    Rational(Integer numerator, BigUnsigned denominator, Reduced);
    friend class FiniteOperations<Rational>;
    friend class ExtensionOperations<Rational>;

  public:
    typedef ExactTag Paradigm;
    typedef Rational NumericType;

  public:
    /** @brief Destroy the object and release its owned resources. */
    ~Rational() = default;
    /** @brief Construct a \c Rational with the specified characteristics. */
    Rational() = default;
    /** @brief Construct a \c Rational with the specified characteristics. */
    Rational(const Integer&, const Integer&);
    /** @brief Construct a \c Rational from n. */
    template<BuiltinIntegral N>
    Rational(N n);
    /** @brief Construct a \c Rational with the specified characteristics. */
    Rational(Int64);
    /** @brief Construct a \c Rational with the specified characteristics. */
    explicit Rational(const FloatDP&);
    /** @brief Construct a \c Rational with the specified characteristics. */
    explicit Rational(const FloatMP&);
    /** @brief Construct a \c Rational with the specified characteristics. */
    Rational(const ExactDouble&);
    /** @brief Construct a \c Rational with the specified characteristics. */
    Rational(const Integer&);
    //! Construct an exact n/2^bits without a generic multiprecision GCD.
    static Rational from_dyadic(Integer numerator, std::size_t bits);
    //! Borrow canonical finite components without allocating or copying.
    const Integer& numerator_reference() const noexcept { return _numerator; }
    /** @brief Return a reference to the denominator. */
    const BigUnsigned& denominator_reference() const noexcept { return _denominator; }
    /** @brief Construct a \c Rational with the specified characteristics. */
    Rational(const Dyadic&);
    /** @brief Construct a \c Rational with the specified characteristics. */
    explicit Rational(const String&);
    /** @brief Copy-construct the object from the supplied instance. */
    Rational(const Rational&) = default;
    /** @brief Move-construct the object from the supplied instance. */
    Rational(Rational&&);
    /** @brief Assign the supplied value and return this object. */
    Rational& operator=(const Rational&);
    /** @brief Assign the supplied value and return this object. */
    Rational& operator=(Rational&&);
    /** @brief Convert to exact number. */
    operator ExactNumber() const;

    /** @brief Return an infinite value. */
    static Rational inf(Sign sgn);
    /** @brief Return an infinite value. */
    static Rational inf();
    /** @brief Return a not-a-number value. */
    static Rational nan();

    /** @brief Return the numerator. */
    Integer get_num() const;
    /** @brief Return the denominator. */
    Integer get_den() const;
    /** @brief Return the numerator. */
    Integer numerator() const;
    /** @brief Return the denominator. */
    Natural denominator() const;

    //! \name Arithmetic operators
    //!@{
    /** @brief Return the sum or unary positive value. */
    friend Rational operator+(const Rational& q) { return pos(q); }
    /** @brief Return the difference or additive inverse. */
    friend Rational operator-(const Rational& q) { return neg(q); }
    /** @brief Return the sum or unary positive value. */
    friend Rational operator+(const Rational& q1, const Rational& q2) { return add(q1, q2); }
    /** @brief Return the difference or additive inverse. */
    friend Rational operator-(const Rational& q1, const Rational& q2) { return sub(q1, q2); }
    /** @brief Return the product. */
    friend Rational operator*(const Rational& q1, const Rational& q2) { return mul(q1, q2); }
    /** @brief Return the quotient. */
    friend Rational operator/(const Rational& q1, const Rational& q2) { return div(q1, q2); }
    /** @brief Add the operand in place. */
    friend Rational& operator+=(Rational& q1, const Rational& q2) { return q1 = add(q1, q2); }
    /** @brief Subtract the operand in place. */
    friend Rational& operator-=(Rational& q1, const Rational& q2) { return q1 = sub(q1, q2); }
    /** @brief Multiply by the operand in place. */
    friend Rational& operator*=(Rational& q1, const Rational& q2) { return q1 = mul(q1, q2); }
    /** @brief Divide by the operand in place. */
    friend Rational& operator/=(Rational& q1, const Rational& q2) { return q1 = div(q1, q2); }
    //@}

    /** @brief Return the quotient. */
    friend Rational operator/(const Integer& z1, const Integer& z2);

    //! \name Comparison operators
    //!@{
    /** @brief Test equality of the operands. */
    friend Boolean operator==(const Rational& q1, const Rational& q2) { return eq(q1, q2); }
    /** @brief Test inequality of the operands. */
    friend Boolean operator!=(const Rational& q1, const Rational& q2) { return !eq(q1, q2); }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend Boolean operator<=(const Rational& q1, const Rational& q2) { return !lt(q2, q1); }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend Boolean operator>=(const Rational& q1, const Rational& q2) { return !lt(q1, q2); }
    /** @brief Test whether the first operand is less than the second. */
    friend Boolean operator<(const Rational& q1, const Rational& q2) { return lt(q1, q2); }
    /** @brief Test whether the first operand is greater than the second. */
    friend Boolean operator>(const Rational& q1, const Rational& q2) { return lt(q2, q1); }
    //@}

    //! \name Named arithmetical functions
    //!@{
    friend Rational nul(const Rational& q);                                          //!< Zero \a 0.
    friend Rational pos(const Rational& q);                                          //!< Identity \a +q.
    friend Rational neg(const Rational& q);                                          //!< Negative \a -q.
    friend Rational hlf(const Rational& q);                                          //!< Half \a q÷2.
    friend Positive<Rational> sqr(const Rational& q);                                //!< Square \a q<sup>2</sup>.
    friend Rational rec(const Rational& q);                                          //!< Reciprocal \a 1/q.
    friend Rational add(const Rational& q1, const Rational& q2);                     //!< \brief Add \a q1+q2.
    friend Rational sub(const Rational& q1, const Rational& q2);                     //!< \brief Subtract \a q1-q2.
    friend Rational mul(const Rational& q1, const Rational& q2);                     //!< \brief Multiply \a q1×q2.
    friend Rational div(const Rational& q1, const Rational& q2);                     //!< \brief Divide \a q1÷q2.
    friend Rational fma(const Rational& q1, const Rational& q2, const Rational& q3); //!< \brief Fused multiply-and-add \a q1×q2+q3.
    friend Rational pow(const Rational& q, Int n);                                   //!< \brief Power \a q<sup>n</sup>.
    /** @brief Raise the value to the specified power. */
    friend Rational pow(const Rational& q, Nat n);
    //!@}

    /** @brief Return the square root. */
    friend Real sqrt(const Real&);
    /** @brief Return the exponential. */
    friend Real exp(const Real&);
    /** @brief Return the natural logarithm. */
    friend Real log(const Real&);
    /** @brief Return the sine. */
    friend Real sin(const Real&);
    /** @brief Return the cosine. */
    friend Real cos(const Real&);
    /** @brief Return the tangent. */
    friend Real tan(const Real&);
    /** @brief Return the inverse sine. */
    friend Real asin(const Real&);
    /** @brief Return the inverse cosine. */
    friend Real acos(const Real&);
    /** @brief Return the inverse tangent. */
    friend Real atan(const Real&);

    //! \name Lattice operations
    //!@{
    friend Positive<Rational> abs(const Rational& q);            //!< Absolute value \a |q|.
    friend Rational min(const Rational& q1, const Rational& q2); //!< Minimum \a q1∧q2.
    /** @brief Return the minimum of the operands. */
    friend Positive<Rational> min(const Positive<Rational>& q1, const Positive<Rational>& q2);
    friend Rational max(const Rational& q1, const Rational& q2); //!< Maximum \a q1∨q2.
    /** @brief Return the maximum of the operands. */
    friend Positive<Rational> max(const Rational& q1, const Positive<Rational>& q2);
    /** @brief Return the maximum of the operands. */
    friend Positive<Rational> max(const Positive<Rational>& q1, const Rational& q2);
    /** @brief Return the maximum of the operands. */
    friend Positive<Rational> max(const Positive<Rational>& q1, const Positive<Rational>& q2);

    /** @brief Return an upper bound on the absolute value. */
    friend Rational mag(const Rational& q);
    /** @brief Return a lower bound on the absolute value. */
    friend Rational mig(const Rational& q);
    //!@}

    //!@{
    //! \name Rounding operations
    friend Integer floor(const Rational& q); //!< Round \a q down to the nearest lower integer.
    friend Integer round(const Rational& q); //!< Round \a q to the nearest integer. %Rounding of halves is implementation-dependent.
    friend Integer ceil(const Rational& q);  //!< Round \a q up to the nearest higher integer.
    //!@}

    //!@{
    //! \name Comparison operations
    friend Sign sgn(const Rational& q);                            //!< The sign of \a q.
    friend Comparison cmp(const Rational& q1, const Rational& q2); //!< Compares which of \a q1 and \a q2 is larger.
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const Rational& q, Int n) { return cmp(q, Rational(n)); }
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(Int n, const Rational& q) { return cmp(Rational(n), q); }
    friend Boolean eq(const Rational& q1, const Rational& q2); //!< Tests if \a q1 is equal to \a q2.
    friend Boolean lt(const Rational& q1, const Rational& q2); //!< Tests if \a q1 is less than \a q2.
    //!@}

    //!@{
    //! \name Special value tests
    friend Bool is_nan(const Rational& q);    //!< Tests whether \a q is NaN (not-a-number).
    friend Bool is_inf(const Rational& q);    //!< Tests whether \a q is ±∞.
    friend Bool is_finite(const Rational& q); //!< Tests whether \a q is finite.
    friend Bool is_zero(const Rational& q);   //!< Tests whether \a q is zero.
    //!@}

    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const Rational& q1, const ExactDouble& d2);
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const ExactDouble& d1, const Rational& q2);

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const Rational& q); //!< Write to an output stream.
    friend InputStream& operator>>(InputStream& os, Rational& q);         //!< Read from an input stream.
    /** @brief Construct the typed value from the literal. */
    friend Rational operator""_q(long double x);
    //! \brief Alternative for operator""_q for use in Python interface.
    friend Rational q_(long double x);
    //!@}

  public:
    //! Binary64 conversion truncating toward zero; underflow yields +0.
    double get_d() const;

  private:
    friend class Dyadic;
};
/** @brief Compile-time predicate testing whether a type is number. */
template<>
struct IsNumber<Rational> : True {};
/** @brief Construct the typed value from the literal. */
Rational operator""_q(unsigned long long int n);
/** @brief Construct the typed value from the literal. */
Rational operator""_q(long double x);

template<BuiltinIntegral N>
inline Rational::Rational(N n) : Rational(Integer(n)) {}

/** @brief Read the value from the input stream. */
InputStream& operator>>(InputStream& is, Rational& q1);

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<>
class Positive<Rational> : public Rational {
  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : Rational() {}
    /** @brief Construct a \c Positive from m. */
    template<BuiltinUnsignedIntegral M>
    Positive(M m) : Rational(m) {}
    /** @brief Construct a \c Positive from q. */
    explicit Positive(const Rational& q) : Rational(q) { ARIADNE_ASSERT(q >= 0); }
};
/** @brief Convert to the positive representation. */
inline Positive<Rational> cast_positive(const Rational& q) { return Positive<Rational>(q); }

using PositiveRational = Positive<Rational>;

using RationalBounds = Bounds<Rational>; //!< Alias for rational bounds on a number. //!< \ingroup NumericModule

/** @brief Validated numeric enclosure with lower and upper endpoints. */
template<>
class Bounds<Rational> {
    Rational _l, _u;

  public:
    typedef ValidatedTag Paradigm;
    /** @brief Construct a \c Bounds from q. */
    Bounds(Rational q) : _l(q), _u(q) {}
    /** @brief Construct a \c Bounds from l, u. */
    Bounds(Rational l, Rational u) : _l(l), _u(u) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<class X>
        requires Constructible<Rational, X>
    Bounds(const Bounds<X>& x)
        : Bounds(Rational(x.lower_raw()), Rational(x.upper_raw())) {}
    /** @brief Convert to validated number. */
    operator ValidatedNumber() const;
    /** @brief Return an enclosure enlarged by the supplied symmetric error. */
    Bounds<Rational> pm(Rational e) { return RationalBounds(_l - e, _u + e); }
    /** @brief Return the lower endpoint. */
    Rational lower() const { return _l; }
    /** @brief Return the upper endpoint. */
    Rational upper() const { return _u; }
    /** @brief Return the underlying lower-endpoint representation. */
    Rational lower_raw() const { return _l; }
    /** @brief Return the underlying upper-endpoint representation. */
    Rational upper_raw() const { return _u; }
    /** @brief Return the sum or unary positive value. */
    friend RationalBounds operator+(const RationalBounds& q) { return RationalBounds(+q._l, +q._u); }
    /** @brief Return the difference or additive inverse. */
    friend RationalBounds operator-(const RationalBounds& q) { return RationalBounds(-q._u, -q._l); }
    /** @brief Return the sum or unary positive value. */
    friend RationalBounds operator+(const RationalBounds& q1, const RationalBounds& q2) { return RationalBounds(q1._l + q2._l, q1._u + q2._u); }
    /** @brief Return the difference or additive inverse. */
    friend RationalBounds operator-(const RationalBounds& q1, const RationalBounds& q2) { return RationalBounds(q1._l - q2._u, q1._u - q2._l); }
    /** @brief Return the product. */
    friend RationalBounds operator*(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Return the quotient. */
    friend RationalBounds operator/(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Return the sum of the operands. */
    friend RationalBounds add(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Return the difference of the operands. */
    friend RationalBounds sub(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Return the product of the operands. */
    friend RationalBounds mul(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Return the quotient of the operands. */
    friend RationalBounds div(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Raise the value to the specified power. */
    friend RationalBounds pow(const RationalBounds& q, Nat m);
    /** @brief Raise the value to the specified power. */
    friend RationalBounds pow(const RationalBounds& q, Int n);
    /** @brief Return the additive zero with the operand characteristics. */
    friend RationalBounds nul(const RationalBounds& q);
    /** @brief Return the unary positive value. */
    friend RationalBounds pos(const RationalBounds& q);
    /** @brief Return the additive inverse. */
    friend RationalBounds neg(const RationalBounds& q);
    /** @brief Divide the value by two. */
    friend RationalBounds hlf(const RationalBounds& q);
    /** @brief Return the square. */
    friend RationalBounds sqr(const RationalBounds& q);
    /** @brief Return the multiplicative inverse. */
    friend RationalBounds rec(const RationalBounds& q);
    /** @brief Return the absolute value. */
    friend RationalBounds abs(const RationalBounds& q);
    /** @brief Return the maximum of the operands. */
    friend RationalBounds max(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Return the minimum of the operands. */
    friend RationalBounds min(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Test equality of the operands. */
    friend ValidatedKleenean operator==(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Test inequality of the operands. */
    friend ValidatedKleenean operator!=(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedKleenean operator<=(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedKleenean operator>=(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const RationalBounds& q1, const RationalBounds& q2);
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(const RationalBounds& q1, const RationalBounds& q2);

    /** @brief Test whether the enclosures cannot represent a common value. */
    friend Boolean inconsistent(const RationalBounds& q1, const RationalBounds& q2) { return q1._l > q2._u || q1._u < q2._l; }
    /** @brief Test whether the first enclosure is contained in the second. */
    friend Boolean refines(const RationalBounds& q1, const RationalBounds& q2) { return q1._l >= q2._l and q1._u <= q2._u; }
    /** @brief Return the common refinement of the enclosures. */
    friend RationalBounds refinement(const RationalBounds& q1, const RationalBounds& q2) { return RationalBounds(max(q1._l, q2._l), min(q1._u, q2._u)); }
    /** @brief Return an enclosure containing both operands. */
    friend RationalBounds coarsening(const RationalBounds& q1, const RationalBounds& q2) { return RationalBounds(min(q1._l, q2._l), max(q1._u, q2._u)); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, RationalBounds y) { return os << "[" << y._l << ":" << y._u << "]"; }
};

} // namespace Ariadne

#endif
