/***************************************************************************
 *            numeric/naive_real.hpp
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

/** @file numeric/naive_real.hpp
 *  @brief Real numbers, including lower, upper and naive types, and validated versions.
 */

#ifndef ARIADNE_NAIVE_REAL_HPP
#define ARIADNE_NAIVE_REAL_HPP

#include "foundation/logical.decl.hpp"
#include "numeric/arithmetic.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/number.decl.hpp"
#include "utility/handle.hpp"

namespace Ariadne {

class Accuracy;

class Real;
class LowerReal;
class UpperReal;
class NaiveReal;

class PositiveNaiveReal;

class ApproximateReal;

class RealInterface;
using NaiveRealInterface = RealInterface;

//! \ingroup NumericModule
//! \brief %Real number type defined as limits of convergent sequences of (dyadic) rationals, without bounds on the convergence rate.
//! \details It is possible to compute arbitrarily accurate (dyadic) rational approximations, but the convergence rate is not known,
//! so at no point in the the computation can anything concrete be deduced about the value of the number.
//!
//! In principle useless for rigorous computation, but quickly-computed approximations to real numbers may be useful for preconditioning rigorous algorithms.
//! \sa Real, LowerReal, UpperReal
class NaiveReal
    : public Handle<const NaiveRealInterface>,
      public DeclareRealOperations<NaiveReal, PositiveNaiveReal>,
      public DeclareAnalyticFieldOperations<NaiveReal>,
      public DeclareLatticeOperations<NaiveReal, PositiveNaiveReal>,
      public DeclareComparisonOperations<NaiveReal, ApproximateKleenean>,
      public DefineFieldOperators<NaiveReal> {
  public:
    typedef NaiveRealInterface Interface;
    typedef ApproximateTag Paradigm;
    typedef NaiveReal NumericType;

  public:
    //!@{
    //! \name Constructors
    NaiveReal();                                        //!< Default constructor yields the integer \c 0 as a real number.
    explicit NaiveReal(SharedPointer<const Interface>); //!< Construct from any class implementing the real number interface.
    // explicit NaiveReal(ConvergentSequence<DyadicBounds> const&); //!< Construct from a sequence of dyadic bounds converging to a singleton intersection.
    //!@}

    //!@{
    //! \name Conversion operators
    NaiveReal(Dbl d);             //!< Construction from a builtin double-precision approximation.
    NaiveReal(const Integer &n);  //!< Construct from an integer.
    NaiveReal(const Dyadic &d);   //!< Construct from a dyadic number \a d=p/2<sup>q</sup> for integers \a p, \a q.
    NaiveReal(const Decimal &d);  //!< Construct from a decimal number, given by its decimal expansion.
    NaiveReal(const Rational &q); //!< Construct from a rational number.

    NaiveReal(const Real &r);       //!< Construct from an (effective) real number.
    NaiveReal(const LowerReal &lr); //!< Construct from a lower real number.
    NaiveReal(const UpperReal &ur); //!< Construct from an upper real number.

    //    operator Number<Tag>() const; //!< Convert to an effective number for computations.
    //!@}

    //!@{
    //! \name Computation of approximations to the number.
    //! Compute an approximation with no guarantees on the error. The error converges to \a 0 as \a eff approaches infinity.
    //! The time taken should be roughly polynomial in \a eff.
    ApproximateReal compute(Effort eff) const;
    //!@}

    //!@{
    //! \name Standard arithmetic operators
    friend NaiveReal operator+(const NaiveReal &r);                       //!< Unary plus.
    friend NaiveReal operator-(const NaiveReal &r);                       //!< Unary minus.
    friend NaiveReal operator+(const NaiveReal &r1, const NaiveReal &r2); //!< Plus.
    friend NaiveReal operator-(const NaiveReal &r1, const NaiveReal &r2); //!< Minus.
    friend NaiveReal operator*(const NaiveReal &r1, const NaiveReal &r2); //!< Times.
    friend NaiveReal operator/(const NaiveReal &r1, const NaiveReal &r2); //!< Divides.
    friend NaiveReal &operator+=(NaiveReal &r1, const NaiveReal &r2);     //!< Inplace plus.
    friend NaiveReal &operator-=(NaiveReal &r1, const NaiveReal &r2);     //!< Inplace minus.
    friend NaiveReal &operator*=(NaiveReal &r1, const NaiveReal &r2);     //!< Inplace times.
    friend NaiveReal &operator/=(NaiveReal &r1, const NaiveReal &r2);     //!< Inplace divides.
    //!@}

    //!@{
    //! \name Named arithmetical functions
    friend NaiveReal pos(const NaiveReal &r);                                            //!< Identity \a +r.
    friend NaiveReal neg(const NaiveReal &r);                                            //!< Negative \a -r.
    friend NaiveReal hlf(const NaiveReal &r);                                            //!< Half \a r÷2.
    friend NaiveReal sqr(const NaiveReal &r);                                            //!< Square \a r<sup>2</sup>.
    friend NaiveReal rec(const NaiveReal &r);                                            //!< Reciprocal \a 1/r.
    friend NaiveReal add(const NaiveReal &r1, const NaiveReal &r2);                      //!< \brief Sum \a r1+r2.
    friend NaiveReal sub(const NaiveReal &r1, const NaiveReal &r2);                      //!< \brief Difference \a r1-r2.
    friend NaiveReal mul(const NaiveReal &r1, const NaiveReal &r2);                      //!< \brief Product \a r1×r2.
    friend NaiveReal div(const NaiveReal &r1, const NaiveReal &r2);                      //!< \brief Quotient \a r1÷r2.
    friend NaiveReal fma(const NaiveReal &r1, const NaiveReal &r2, const NaiveReal &r3); //!< \brief Fused multiply-and-add \a r1×r2+r3.
    friend NaiveReal pow(const NaiveReal &r, Int n);                                     //!< \brief Power \a r<sup>n</sup>.
    //!@}

    //!@{
    //! \name Algebraic and transcendental functions
    friend NaiveReal sqrt(const NaiveReal &r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend NaiveReal exp(const NaiveReal &r);  //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend NaiveReal log(const NaiveReal &r);  //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend NaiveReal sin(const NaiveReal &r);  //!< The sine of \a r.
    friend NaiveReal cos(const NaiveReal &r);  //!< The cosine of \a r.
    friend NaiveReal tan(const NaiveReal &r);  //!< The tangent of \a r, sin(\a r)/cos(\a r).
    friend NaiveReal atan(const NaiveReal &r); //!< The arc-tangent of \a r.
    //!@}

    //!@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(const NaiveReal &);                //!< Absolute value \a |r|.
    friend NaiveReal min(const NaiveReal &r1, const NaiveReal &r2); //!< The mimimum of \a r1 and \a r2.
    friend NaiveReal max(const NaiveReal &r1, const NaiveReal &r2); //!< The maximum of \a r1 and \a r2.
    //!@}

    //!@{
    //! Operations based on the metric structure.
    friend PositiveNaiveReal dist(const NaiveReal &r1, const NaiveReal &r2);
    //< The distance |\a r <sub>1</sub>-\a r <sub>2</sub>| between \a r<sub>1</sub> and \a r<sub>2</sub>.
    //!@}

    //!@{
    //! \name Comparison operations and operators.
    friend ApproximateKleenean leq(const NaiveReal &r1, const NaiveReal &r2);
    //!< Returns \c likely if \a r1<r2, \c unlikely if \a r1\>r2, and \c indetermiate if \a r1==r2.
    //!@}

    //!@{
    //! \name Limit operations.
    friend NaiveReal limit(const ConvergentSequence<Dyadic> &qas);
    //!< Create a real number from an converging sequence of dyadic approximants.
    //!@}

    //!@{
    //! \name Input/output operations
    friend OutputStream &operator<<(OutputStream &os, const NaiveReal &r); //< Write to an output stream.
    //!@}
};

//! \ingroup UserNumericTypeSubModule
//! \brief Computable lower real numbers defined by conversion to concrete floats.
class PositiveNaiveReal : public NaiveReal {
  public:
    PositiveNaiveReal() : NaiveReal() {}
    template<ConvertibleTo<NaiveReal> Y>
    PositiveNaiveReal(const Positive<Y> &p) : NaiveReal(p) {}
    explicit PositiveNaiveReal(const NaiveReal &r) : NaiveReal(r) {}
    PositiveBounds<Dyadic> compute_get(Effort) const;

  public:
    PositiveNaiveReal max(const PositiveNaiveReal &, const PositiveNaiveReal &);
    PositiveNaiveReal max(const NaiveReal &, const PositiveNaiveReal &);
    PositiveNaiveReal max(const PositiveNaiveReal &, const NaiveReal &);

    PositiveNaiveReal min(const PositiveNaiveReal &, const PositiveNaiveReal &);

    PositiveNaiveReal rec(const PositiveNaiveReal &);
    PositiveNaiveReal add(const PositiveNaiveReal &, const PositiveNaiveReal &);
    PositiveNaiveReal mul(const PositiveNaiveReal &, const PositiveNaiveReal &);
    PositiveNaiveReal div(const PositiveNaiveReal &, const PositiveNaiveReal &);
    PositiveNaiveReal sqrt(const PositiveNaiveReal &);
    PositiveNaiveReal atan(const PositiveNaiveReal &);

  public:
    friend PositiveNaiveReal min(const PositiveNaiveReal &pr1, const PositiveNaiveReal &pr2);
    friend PositiveNaiveReal add(const PositiveNaiveReal &pr1, const PositiveNaiveReal &pr2);
    friend PositiveNaiveReal mul(const PositiveNaiveReal &pr1, const PositiveNaiveReal &pr2);
    friend PositiveNaiveReal div(const PositiveNaiveReal &pr1, const PositiveNaiveReal &pr2);
    friend PositiveNaiveReal rec(const PositiveNaiveReal &pr);
};

PositiveNaiveReal cast_positive(const NaiveReal &x);

} // namespace Ariadne

#endif
