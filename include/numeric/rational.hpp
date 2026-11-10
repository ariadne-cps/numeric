/***************************************************************************
 *            numeric/rational.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/rational.hpp
 *  \brief
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
    struct Reduced {};
    Rational(Integer numerator, BigUnsigned denominator, Reduced);
    friend class FiniteOperations<Rational>;
    friend class ExtensionOperations<Rational>;

  public:
    typedef ExactTag Paradigm;
    typedef Rational NumericType;

  public:
    ~Rational() = default;
    Rational() = default;
    Rational(const Integer&, const Integer&);
    template<BuiltinIntegral N>
    Rational(N n);
    Rational(Int64);
    explicit Rational(const FloatDP&);
    explicit Rational(const FloatMP&);
    Rational(const ExactDouble&);
    Rational(const Integer&);
    //! Construct an exact n/2^bits without a generic multiprecision GCD.
    static Rational from_dyadic(Integer numerator, std::size_t bits);
    //! Borrow canonical finite components without allocating or copying.
    const Integer& numerator_reference() const noexcept { return _numerator; }
    const BigUnsigned& denominator_reference() const noexcept { return _denominator; }
    Rational(const Dyadic&);
    explicit Rational(const String&);
    Rational(const Rational&) = default;
    Rational(Rational&&);
    Rational& operator=(const Rational&);
    Rational& operator=(Rational&&);
    operator ExactNumber() const;

    static Rational inf(Sign sgn);
    static Rational inf();
    static Rational nan();

    Integer get_num() const;
    Integer get_den() const;
    Integer numerator() const;
    Natural denominator() const;

    //! \name Arithmetic operators
    //!@{
    friend Rational operator+(const Rational& q) { return pos(q); }
    friend Rational operator-(const Rational& q) { return neg(q); }
    friend Rational operator+(const Rational& q1, const Rational& q2) { return add(q1, q2); }
    friend Rational operator-(const Rational& q1, const Rational& q2) { return sub(q1, q2); }
    friend Rational operator*(const Rational& q1, const Rational& q2) { return mul(q1, q2); }
    friend Rational operator/(const Rational& q1, const Rational& q2) { return div(q1, q2); }
    friend Rational& operator+=(Rational& q1, const Rational& q2) { return q1 = add(q1, q2); }
    friend Rational& operator-=(Rational& q1, const Rational& q2) { return q1 = sub(q1, q2); }
    friend Rational& operator*=(Rational& q1, const Rational& q2) { return q1 = mul(q1, q2); }
    friend Rational& operator/=(Rational& q1, const Rational& q2) { return q1 = div(q1, q2); }
    //@}

    friend Rational operator/(const Integer& z1, const Integer& z2);

    //! \name Comparison operators
    //!@{
    friend Boolean operator==(const Rational& q1, const Rational& q2) { return eq(q1, q2); }
    friend Boolean operator!=(const Rational& q1, const Rational& q2) { return !eq(q1, q2); }
    friend Boolean operator<=(const Rational& q1, const Rational& q2) { return !lt(q2, q1); }
    friend Boolean operator>=(const Rational& q1, const Rational& q2) { return !lt(q1, q2); }
    friend Boolean operator<(const Rational& q1, const Rational& q2) { return lt(q1, q2); }
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
    friend Rational pow(const Rational& q, Nat n);
    //!@}

    friend Real sqrt(const Real&);
    friend Real exp(const Real&);
    friend Real log(const Real&);
    friend Real sin(const Real&);
    friend Real cos(const Real&);
    friend Real tan(const Real&);
    friend Real asin(const Real&);
    friend Real acos(const Real&);
    friend Real atan(const Real&);

    //! \name Lattice operations
    //!@{
    friend Positive<Rational> abs(const Rational& q);            //!< Absolute value \a |q|.
    friend Rational min(const Rational& q1, const Rational& q2); //!< Minimum \a q1∧q2.
    friend Positive<Rational> min(const Positive<Rational>& q1, const Positive<Rational>& q2);
    friend Rational max(const Rational& q1, const Rational& q2); //!< Maximum \a q1∨q2.
    friend Positive<Rational> max(const Rational& q1, const Positive<Rational>& q2);
    friend Positive<Rational> max(const Positive<Rational>& q1, const Rational& q2);
    friend Positive<Rational> max(const Positive<Rational>& q1, const Positive<Rational>& q2);

    friend Rational mag(const Rational& q);
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
    friend Comparison cmp(const Rational& q, Int n) { return cmp(q, Rational(n)); }
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

    friend Comparison cmp(const Rational& q1, const ExactDouble& d2);
    friend Comparison cmp(const ExactDouble& d1, const Rational& q2);

    //!@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, const Rational& q); //!< Write to an output stream.
    friend InputStream& operator>>(InputStream& os, Rational& q);         //!< Read from an input stream.
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
template<>
struct IsNumber<Rational> : True {};
Rational operator""_q(unsigned long long int n);
Rational operator""_q(long double x);

template<BuiltinIntegral N>
inline Rational::Rational(N n) : Rational(Integer(n)) {}

InputStream& operator>>(InputStream& is, Rational& q1);

template<>
class Positive<Rational> : public Rational {
  public:
    Positive() : Rational() {}
    template<BuiltinUnsignedIntegral M>
    Positive(M m) : Rational(m) {}
    explicit Positive(const Rational& q) : Rational(q) { ARIADNE_ASSERT(q >= 0); }
};
inline Positive<Rational> cast_positive(const Rational& q) { return Positive<Rational>(q); }

using PositiveRational = Positive<Rational>;

using RationalBounds = Bounds<Rational>; //!< Alias for rational bounds on a number. //!< \ingroup NumericModule

template<>
class Bounds<Rational> {
    Rational _l, _u;

  public:
    typedef ValidatedTag Paradigm;
    Bounds(Rational q) : _l(q), _u(q) {}
    Bounds(Rational l, Rational u) : _l(l), _u(u) {}
    template<class X>
        requires Constructible<Rational, X>
    Bounds(const Bounds<X>& x)
        : Bounds(Rational(x.lower_raw()), Rational(x.upper_raw())) {}
    operator ValidatedNumber() const;
    Bounds<Rational> pm(Rational e) { return RationalBounds(_l - e, _u + e); }
    Rational lower() const { return _l; }
    Rational upper() const { return _u; }
    Rational lower_raw() const { return _l; }
    Rational upper_raw() const { return _u; }
    friend RationalBounds operator+(const RationalBounds& q) { return RationalBounds(+q._l, +q._u); }
    friend RationalBounds operator-(const RationalBounds& q) { return RationalBounds(-q._u, -q._l); }
    friend RationalBounds operator+(const RationalBounds& q1, const RationalBounds& q2) { return RationalBounds(q1._l + q2._l, q1._u + q2._u); }
    friend RationalBounds operator-(const RationalBounds& q1, const RationalBounds& q2) { return RationalBounds(q1._l - q2._u, q1._u - q2._l); }
    friend RationalBounds operator*(const RationalBounds& q1, const RationalBounds& q2);
    friend RationalBounds operator/(const RationalBounds& q1, const RationalBounds& q2);
    friend RationalBounds add(const RationalBounds& q1, const RationalBounds& q2);
    friend RationalBounds sub(const RationalBounds& q1, const RationalBounds& q2);
    friend RationalBounds mul(const RationalBounds& q1, const RationalBounds& q2);
    friend RationalBounds div(const RationalBounds& q1, const RationalBounds& q2);
    friend RationalBounds pow(const RationalBounds& q, Nat m);
    friend RationalBounds pow(const RationalBounds& q, Int n);
    friend RationalBounds nul(const RationalBounds& q);
    friend RationalBounds pos(const RationalBounds& q);
    friend RationalBounds neg(const RationalBounds& q);
    friend RationalBounds hlf(const RationalBounds& q);
    friend RationalBounds sqr(const RationalBounds& q);
    friend RationalBounds rec(const RationalBounds& q);
    friend RationalBounds abs(const RationalBounds& q);
    friend RationalBounds max(const RationalBounds& q1, const RationalBounds& q2);
    friend RationalBounds min(const RationalBounds& q1, const RationalBounds& q2);
    friend ValidatedKleenean operator==(const RationalBounds& q1, const RationalBounds& q2);
    friend ValidatedKleenean operator!=(const RationalBounds& q1, const RationalBounds& q2);
    friend ValidatedKleenean operator<=(const RationalBounds& q1, const RationalBounds& q2);
    friend ValidatedKleenean operator>=(const RationalBounds& q1, const RationalBounds& q2);
    friend ValidatedKleenean operator<(const RationalBounds& q1, const RationalBounds& q2);
    friend ValidatedKleenean operator>(const RationalBounds& q1, const RationalBounds& q2);

    friend Boolean inconsistent(const RationalBounds& q1, const RationalBounds& q2) { return q1._l > q2._u || q1._u < q2._l; }
    friend Boolean refines(const RationalBounds& q1, const RationalBounds& q2) { return q1._l >= q2._l and q1._u <= q2._u; }
    friend RationalBounds refinement(const RationalBounds& q1, const RationalBounds& q2) { return RationalBounds(max(q1._l, q2._l), min(q1._u, q2._u)); }
    friend RationalBounds coarsening(const RationalBounds& q1, const RationalBounds& q2) { return RationalBounds(min(q1._l, q2._l), max(q1._u, q2._u)); }
    friend OutputStream& operator<<(OutputStream& os, RationalBounds y) { return os << "[" << y._l << ":" << y._u << "]"; }
};

} // namespace Ariadne

#endif
