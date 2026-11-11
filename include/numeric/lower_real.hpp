/***************************************************************************
 *            numeric/lower_real.hpp
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

/** @file numeric/lower_real.hpp
 *  @brief Real numbers, including lower, upper and naive types, and validated versions.
 */

#ifndef ARIADNE_LOWER_REAL_HPP
#define ARIADNE_LOWER_REAL_HPP

#include "foundation/logical.decl.hpp"
#include "foundation/paradigm.hpp"
#include "numeric/approximate_real.hpp"
#include "numeric/arithmetic.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/real.hpp"
#include "numeric/sequence.hpp"
#include "numeric/validated_real.hpp"
#include "utility/handle.hpp"
#include "utility/pointer.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

class Accuracy;

class Real;
class LowerReal;
class UpperReal;
class NaiveReal;
class PositiveReal;
class PositiveLowerReal;
class PositiveUpperReal;
class PositiveNaiveReal;

class ValidatedReal;
class ValidatedUpperReal;
class ValidatedLowerReal;
class ApproximateReal;

class RealInterface;
using LowerRealInterface = RealInterface;
using UpperRealInterface = RealInterface;
using NaiveRealInterface = RealInterface;

class ValidatedRealInterface;

//! \ingroup NumericModule
//! \brief Lower real number type \f$\R_<\f$.
//! \details An effectively computable <em>lower real</em> is a real number for which it is possible to compute arbitrarily accurate (dyadic) rational lower bounds;
//! equivalently, a bounded increasing sequence of (dyadic) rationals converging to the number.
//! However, the convergence rate is not known, and there may be arbitrarily large upward jumps in the sequence.
//!
//! Multiplication and division are unsupported, since given \f$x \geq \underline{x}\f$ and \f$y\geq\underline{y}\f$,
//! we cannot deduce \f$x \times y\geq\underline{x}\times\underline{y}\f$ if the bounds are negative.
//! However, multiplication and reciprocation of positive lower reals <em>is</em> supported.
//!
//! %Positive lower real numbers are a natural type for representing the measure of an open set,
//! or the separation between two points in a metric space.
//! \sa Real, UpperReal, NaiveReal
class LowerReal
    : public Handle<const LowerRealInterface> {
  public:
    typedef EffectiveTag Paradigm;
    typedef LowerReal NumericType;

  public:
    /** @brief Construct a \c LowerReal with the specified characteristics. */
    LowerReal(Real);
    /** @brief Construct a \c LowerReal with the specified characteristics. */
    explicit LowerReal(SharedPointer<const Interface>);

  public:
    //!@{
    //! \name Computation of rigorous validated lower bounds on the number
    /** @brief Compute the value at the requested accuracy or effort. */
    ValidatedLowerReal compute(Effort eff) const;
    /** @brief Compute the value and extract the requested numeric representation. */
    DyadicLowerBound compute_get(Effort eff) const;
    /** @brief Compute the value and extract the requested numeric representation. */
    FloatDPLowerBound compute_get(Effort eff, DoublePrecision pr) const;
    /** @brief Compute the value and extract the requested numeric representation. */
    FloatMPLowerBound compute_get(Effort eff, MultiplePrecision pr) const;
    //!@}

  public:
    //!@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(const LowerReal &r) = delete;      //!< \em No absolute value operator!
    friend LowerReal min(const LowerReal &r1, const LowerReal &r2); //!< The mimimum of \a r1 and \a r2.
    friend LowerReal max(const LowerReal &r1, const LowerReal &r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //!@{
    //! \name Named arithmetical functions
    friend LowerReal pos(const LowerReal &r);                       //!< Identity \a +r.
    friend UpperReal neg(const LowerReal &r);                       //!< Negative \a -r.
    friend LowerReal neg(const UpperReal &r);                       //!< Negative \a -r.
    friend LowerReal hlf(const LowerReal &r);                       //!< Half \a r÷2.
    friend LowerReal add(const LowerReal &r1, const LowerReal &r2); //!< \brief Sum \a r1+r2.
    friend LowerReal sub(const LowerReal &r1, const UpperReal &r2); //!< \brief Difference \a r1-r2.
    friend UpperReal sub(const UpperReal &r1, const LowerReal &r2); //!< \brief Difference \a r1-r2.

    friend NaiveReal mul(const LowerReal &r1, const LowerReal &r2) = delete; //!< \brief \em No multiplication operator, since non-monotone!
    friend NaiveReal div(const LowerReal &r1, const UpperReal &r2) = delete; //!< \brief \em No division operator, since non-monotone!
    friend NaiveReal rec(const LowerReal &r) = delete;                       //!< \brief \em No reciprocal operation, since non-monotone!
    //!@}

    //!@{
    //! \name Named arithmetical functions on positive numbers
    friend LowerReal mul(const LowerReal &r1, const PositiveReal &r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend LowerReal mul(const PositiveReal &r1, const LowerReal &r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend LowerReal div(const LowerReal &r1, const PositiveReal &r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend LowerReal div(const PositiveReal &r1, const UpperReal &r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend UpperReal div(const PositiveReal &r1, const LowerReal &r2); //!< Division \a r1÷r2 preserves monotonicity.

    friend PositiveLowerReal add(const PositiveLowerReal &r1, const PositiveLowerReal &r2); //!< Addition \a r1+r2 preserves positivity.
    friend PositiveLowerReal mul(const PositiveLowerReal &r1, const PositiveLowerReal &r2); //!< Multiplication \a r1×r2 preserves positivity.
    friend PositiveLowerReal div(const PositiveLowerReal &r1, const PositiveUpperReal &r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveUpperReal div(const PositiveUpperReal &r1, const PositiveLowerReal &r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveLowerReal rec(const PositiveUpperReal &r);                               //!< Reciprocal \a 1/r preserves positivity.
    friend PositiveUpperReal rec(const PositiveLowerReal &r);                               //!< Reciprocal \a 1/r preserves positivity.
    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend PositiveLowerReal sqrt(const PositiveLowerReal &r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend PositiveLowerReal exp(const LowerReal &r);          //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend LowerReal log(const PositiveLowerReal &r);          //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend LowerReal atan(const LowerReal &r);                 //!< The arc-tangent of \a r.
    friend PositiveLowerReal atan(const PositiveLowerReal &r); //!< Arc-tangent preserves positivity

    friend NaiveReal sin(const LowerReal &r); //!< No sine operation, since non-monotone!
    friend NaiveReal cos(const LowerReal &r); //!< No cosine operation, since non-monotone!
    friend NaiveReal tan(const LowerReal &r); //!< No tangent operation, since non-monotone!
    //!@}

    //!@{
    //! \name Limit operations.
    friend LowerReal limit(const IncreasingSequence<Dyadic> &qlbs);
    //!< Create a real number from an increasing sequence of dyadic lower bounds.
    //!@}

    //!@{
    //! \name Input/output operations
    /** @brief Write the value to the output stream. */
    friend OutputStream &operator<<(OutputStream &os, const LowerReal &r); //< Write to an output stream.
    //!@}

    //!@{
    //! \name Mixed operations
    /** @brief Return the minimum of the operands. */
    friend Real min(const LowerReal &lr1, const Real &r2);
    /** @brief Return the minimum of the operands. */
    friend Real min(const Real &r1, const LowerReal &lr2);
    //!@}
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class PositiveLowerReal : public LowerReal {
  public:
    /** @brief Construct a \c PositiveLowerReal from r. */
    PositiveLowerReal(PositiveReal r) : LowerReal(r) {}
    /** @brief Construct a \c PositiveLowerReal from r. */
    explicit PositiveLowerReal(LowerReal r) : LowerReal(r) {}
    /** @brief Compute the value and extract the requested numeric representation. */
    PositiveLowerBound<Dyadic> compute_get(Effort) const;

  public:
    /** @brief Return the multiplicative inverse. */
    PositiveLowerReal rec(const PositiveUpperReal &);
    /** @brief Return the multiplicative inverse. */
    PositiveUpperReal rec(const PositiveLowerReal &);
    /** @brief Return the sum of the operands. */
    PositiveLowerReal add(const PositiveLowerReal &, const PositiveLowerReal &);
    /** @brief Return the product of the operands. */
    PositiveLowerReal mul(const PositiveLowerReal &, const PositiveLowerReal &);
    /** @brief Return the quotient of the operands. */
    PositiveLowerReal div(const PositiveLowerReal &, const PositiveUpperReal &);
    /** @brief Return the quotient of the operands. */
    PositiveUpperReal div(const PositiveUpperReal &, const PositiveLowerReal &);
};

/** @brief Return the multiplicative inverse. */
PositiveUpperReal rec(const PositiveLowerReal &plr);
/** @brief Return the sum of the operands. */
PositiveLowerReal add(const PositiveLowerReal &plr1, const PositiveLowerReal &plr2);
/** @brief Return the product of the operands. */
PositiveLowerReal mul(const PositiveLowerReal &plr1, const PositiveLowerReal &plr2);
/** @brief Return the quotient of the operands. */
PositiveLowerReal div(const PositiveLowerReal &plr1, const PositiveUpperReal &pur2);

/** @brief Validated lower information about a real number. */
class ValidatedLowerReal
    : public Handle<const ValidatedRealInterface> {
  public:
    /** @brief Construct a \c ValidatedLowerReal with the specified characteristics. */
    ValidatedLowerReal(const DyadicLowerBound &);
    /** @brief Construct a \c ValidatedLowerReal with the specified characteristics. */
    ValidatedLowerReal(const ValidatedReal &);
    /** @brief Return the requested element or stored value. */
    DyadicLowerBound get() const;
    /** @brief Return the requested element or stored value. */
    FloatDPLowerBound get(DoublePrecision) const;
    /** @brief Return the requested element or stored value. */
    FloatMPLowerBound get(MultiplePrecision) const;
    /** @brief Write the value to the output stream. */
    friend OutputStream &operator<<(OutputStream &, const ValidatedLowerReal &);
};

} // namespace Ariadne

#endif
