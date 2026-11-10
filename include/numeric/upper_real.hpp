/***************************************************************************
 *            numeric/upper_real.hpp
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

/*! \file numeric/upper_real.hpp
 *  \brief Real numbers, including lower, upper and naive types, and validated versions.
 */

#ifndef ARIADNE_UPPER_REAL_HPP
#define ARIADNE_UPPER_REAL_HPP

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
//! \brief Upper real number type \f$\R_>\f$.
//! \details An effectively computable <em>upper real</em> is a real number for which it is possible to compute arbitrarily accurate (dyadic) rational upper bounds;
//! equivalently, a bounded decreasing sequence of (dyadic) rationals converging to the number.
//! However, the convergence rate is not known, and there may be arbitrarily large downward jumps in the sequence.
//!
//! Multiplication and division are unsupported, since given \f$x \leq \overline{x}\f$ and \f$y\leq\overline{y}\f$,
//! we cannot deduce \f$x \times y\leq\overline{x}\times\overline{y}\f$.
//! However, multiplication and reciprocation of positive upper reals <em>is</em> supported.
//!
//! %Positive upper real numbers are a natural type for representing upper bounds of a distance or norm,
//! including the radius of a ball, or the measure of an closed set.
//! \sa Real, LowerReal, NaiveReal
class UpperReal
    : public Handle<const UpperRealInterface> {
  public:
    typedef EffectiveTag Paradigm;
    typedef UpperReal NumericType;

  public:
    UpperReal(Real);
    explicit UpperReal(SharedPointer<const Interface>);

  public:
    //!@{
    //! \name Computation of rigorous validated upper bounds on the number
    ValidatedUpperReal compute(Effort eff) const;
    DyadicUpperBound compute_get(Effort eff) const;
    FloatDPUpperBound compute_get(Effort eff, DoublePrecision pr) const;
    FloatMPUpperBound compute_get(Effort eff, MultiplePrecision pr) const;
    //!@}

  public:
    //!@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(const UpperReal &r) = delete;      //!< \em No absolute value operator!
    friend UpperReal min(const UpperReal &r1, const UpperReal &r2); //!< The mimimum of \a r1 and \a r2.
    friend UpperReal max(const UpperReal &r1, const UpperReal &r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //!@{
    //! \name Named arithmetical functions
    friend UpperReal pos(const UpperReal &r);                       //!< Identity \a +r.
    friend LowerReal neg(const UpperReal &r);                       //!< Negative \a -r.
    friend UpperReal neg(const LowerReal &r);                       //!< Negative \a -r.
    friend UpperReal hlf(const UpperReal &r);                       //!< Half \a r÷2.
    friend UpperReal add(const UpperReal &r1, const UpperReal &r2); //!< \brief Sum \a r1+r2.
    friend UpperReal sub(const UpperReal &r1, const LowerReal &r2); //!< \brief Difference \a r1-r2.
    friend LowerReal sub(const LowerReal &r1, const UpperReal &r2); //!< \brief Difference \a r1-r2.

    friend NaiveReal mul(const UpperReal &r1, const UpperReal &r2) = delete; //!< \brief \em No multiplication operator, since non-monotone!
    friend NaiveReal mul(const UpperReal &r1, const LowerReal &r2) = delete; //!< \brief \em No division operator, since non-monotone!
    friend NaiveReal rec(const UpperReal &r) = delete;                       //!< \brief \em No reciprocal operation, since non-monotone!
    //!@}

    //!@{
    //! \name Named arithmetical functions on positive numbers
    friend UpperReal mul(const UpperReal &r1, const PositiveReal &r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend UpperReal mul(const PositiveReal &r1, const UpperReal &r2); //!< Multiplication \a r1×r2 preserves monotonicity.
    friend UpperReal div(const UpperReal &r1, const PositiveReal &r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend UpperReal div(const PositiveReal &r1, const LowerReal &r2); //!< Division \a r1÷r2 preserves monotonicity.
    friend LowerReal div(const PositiveReal &r1, const UpperReal &r2); //!< Division \a r1÷r2 preserves monotonicity.

    friend PositiveUpperReal add(const PositiveUpperReal &r1, const PositiveUpperReal &r2); //!< Addition \a r1+r2 preserves positivity.
    friend PositiveUpperReal mul(const PositiveUpperReal &r1, const PositiveUpperReal &r2); //!< Multiplication \a r1×r2 preserves positivity.
    friend PositiveUpperReal div(const PositiveUpperReal &r1, const PositiveLowerReal &r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveLowerReal div(const PositiveLowerReal &r1, const PositiveUpperReal &r2); //!< Division \a r1÷r2 preserves positivity.
    friend PositiveUpperReal rec(const PositiveLowerReal &r);                               //!< Reciprocal \a 1/r preserves positivity.
    friend PositiveLowerReal rec(const PositiveUpperReal &r);                               //!< Reciprocal \a 1/r preserves positivity.
    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend PositiveUpperReal sqrt(const PositiveUpperReal &r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend PositiveUpperReal exp(const UpperReal &r);          //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend UpperReal log(const PositiveUpperReal &r);          //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend UpperReal atan(const UpperReal &r);                 //!< The arc-tangent of \a r.
    friend PositiveUpperReal atan(const PositiveUpperReal &r); //!< Arc-tangent preserves positivity
    friend NaiveReal sin(const UpperReal &r);                  //!< No sine operation, since non-monotone!
    friend NaiveReal cos(const UpperReal &r);                  //!< No cosine operation, since non-monotone!
    friend NaiveReal tan(const UpperReal &r);                  //!< No tangent operation, since non-monotone!
    //!@}

    //!@{
    //! \name Limit operations.
    friend UpperReal limit(const DecreasingSequence<Dyadic> &qubs);
    //!< Create a real number from an decreasing sequence of dyadic upper bounds.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream &operator<<(OutputStream &os, const UpperReal &r); //< Write to an output stream.
    //!@}

    //!@{
    //! \name Mixed operations
    friend Real max(const UpperReal &ur1, const Real &r2);
    friend Real max(Real r1, const UpperReal &ur2);
    //!@}
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable upper real numbers defined by conversion to concrete floats.
class PositiveUpperReal : public UpperReal {
  public:
    PositiveUpperReal(PositiveReal r) : UpperReal(r) {}
    explicit PositiveUpperReal(UpperReal r) : UpperReal(r) {}
    PositiveUpperBound<Dyadic> compute_get(Effort) const;

  public:
    PositiveUpperReal rec(const PositiveLowerReal &);
    PositiveLowerReal rec(const PositiveUpperReal &);
    PositiveUpperReal add(const PositiveUpperReal &, const PositiveUpperReal &);
    PositiveUpperReal mul(const PositiveUpperReal &, const PositiveUpperReal &);
    PositiveUpperReal div(const PositiveUpperReal &, const PositiveLowerReal &);
    PositiveLowerReal div(const PositiveLowerReal &, const PositiveUpperReal &);
};

PositiveLowerReal rec(const PositiveUpperReal &pur);
PositiveUpperReal add(const PositiveUpperReal &pur1, const PositiveUpperReal &pur2);
PositiveUpperReal mul(const PositiveUpperReal &pur1, const PositiveUpperReal &pur2);
PositiveUpperReal div(const PositiveUpperReal &pur1, const PositiveLowerReal &plr2);

class ValidatedUpperReal
    : public Handle<const ValidatedRealInterface> {
  public:
    ValidatedUpperReal(const DyadicUpperBound &);
    ValidatedUpperReal(const ValidatedReal &);
    DyadicUpperBound get() const;
    FloatDPUpperBound get(DoublePrecision) const;
    FloatMPUpperBound get(MultiplePrecision) const;
    friend OutputStream &operator<<(OutputStream &, const ValidatedUpperReal &);
};

} // namespace Ariadne

#endif
