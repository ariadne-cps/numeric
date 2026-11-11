/***************************************************************************
 *            numeric/float_factory.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file numeric/float_factory.hpp
 *  @brief Factories for creating floating-point numbers.
 */

#ifndef ARIADNE_FLOAT_FACTORY_HPP
#define ARIADNE_FLOAT_FACTORY_HPP

#include "float.decl.hpp"
#include "number.decl.hpp"

namespace Ariadne {

/** @brief Factory constructing compatible float. */
template<class PR>
class FloatFactory {
  protected:
    PR _pr;
    typedef RawFloatType<PR> F;

  public:
    typedef PR PrecisionType;
    typedef PR CharacteristicsType;
    /** @brief Construct a \c FloatFactory from pr. */
    FloatFactory(PR const& pr) : _pr(pr) {}
    /** @brief Return the working precision. */
    PR precision() const { return this->_pr; }
    /** @brief Return the characteristics required to construct compatible values. */
    PR characteristics() const { return this->_pr; }

  public:
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatApproximation<PR> create(const ApproximateNumber& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatLowerBound<PR> create(const ValidatedLowerNumber& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatUpperBound<PR> create(const ValidatedUpperNumber& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBounds<PR> create(const ValidatedNumber& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBounds<PR> create(const EffectiveNumber& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBounds<PR> create(const ExactNumber& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBounds<PR> create(const Real& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBounds<PR> create(const Rational& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBounds<PR> create(const Dyadic& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBounds<PR> create(const Integer& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    Float<PR> create(const ExactDouble& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    Float<PR> create(const Dyadic& y, ExactTag);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    Float<PR> create(const Integer& y, ExactTag);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    template<BuiltinSignedIntegral N>
    Float<PR> create(const N& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    template<BuiltinUnsignedIntegral M>
    PositiveFloat<PR> create(const M& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    template<BuiltinFloatingPoint D>
    FloatApproximation<PR> create(const D& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    PositiveFloatApproximation<PR> create(const PositiveApproximateNumber& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    PositiveFloatLowerBound<PR> create(const PositiveValidatedLowerNumber& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    PositiveFloatUpperBound<PR> create(const PositiveValidatedUpperNumber& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    PositiveFloatBounds<PR> create(const PositiveValidatedNumber& y);
};

/** @brief Factory constructing compatible float ball. */
template<class PR, class PRE>
class FloatBallFactory : public FloatFactory<PR> {
    PRE _pre;

  public:
    /** @brief Construct a \c FloatBallFactory from pr, pre. */
    FloatBallFactory(PR const& pr, PRE const& pre) : FloatFactory<PR>(pr), _pre(pre) {}
    using FloatFactory<PR>::create;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBall<PR, PRE> create(const Number<ValidatedTag>& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBall<PR, PRE> create(const Number<EffectiveTag>& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBall<PR, PRE> create(const Number<ExactTag>& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBall<PR, PRE> create(const Real& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBall<PR, PRE> create(const Rational& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBall<PR, PRE> create(const Dyadic& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    FloatBall<PR, PRE> create(const Integer& y);
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    PositiveFloatBall<PR, PRE> create(const PositiveNumber<ValidatedTag>& y);
};

template<class Y, class PR>
using ConcreteType = decltype(declval<FloatFactory<PR>>().create(declval<Y>()));
/** @brief Construct float from the supplied arguments. */
template<class Y, class PR>
inline decltype(auto) make_float(const Y& y, PR pr) { return float_factory(pr).create(y); }

/** @brief Return a floating-point factory with the supplied precision. */
template<class PR>
inline FloatFactory<PR> float_factory(PR pr) { return FloatFactory<PR>(pr); }
template<class F>
inline FloatFactory<PrecisionType<F>> factory(const Approximation<F>& flt);
template<class F>
inline FloatFactory<PrecisionType<F>> factory(const LowerBound<F>& flt);
template<class F>
inline FloatFactory<PrecisionType<F>> factory(const UpperBound<F>& flt);
template<class F>
inline FloatFactory<PrecisionType<F>> factory(const Bounds<F>& flt);
template<class F, class FE>
inline FloatBallFactory<PrecisionType<F>, PrecisionType<FE>> factory(const Ball<F, FE>& flt);
template<ARawFloat F>
inline FloatFactory<PrecisionType<F>> factory(const F& flt);

/** @brief Return a factory for compatible objects. */
template<ARawFloat F>
inline FloatFactory<PrecisionType<F>> factory(const F& flt) {
    return FloatFactory<PrecisionType<F>>(flt.precision());
}
template<class PR>
inline Float<PR> FloatFactory<PR>::create(const Dyadic& y, ExactTag) { return Float<PR>(y, _pr); }
template<class PR>
inline Float<PR> FloatFactory<PR>::create(const Integer& y, ExactTag) { return Float<PR>(y, _pr); }
template<class PR>
inline Float<PR> FloatFactory<PR>::create(const ExactDouble& y) { return Float<PR>(y, _pr); }

template<class PR>
template<BuiltinSignedIntegral N>
inline Float<PR> FloatFactory<PR>::create(const N& y) { return Float<PR>(y, _pr); }
template<class PR>
/** @brief Construct a compatible object from the supplied value and characteristics. */
template<BuiltinUnsignedIntegral M>
inline Positive<Float<PR>> FloatFactory<PR>::create(const M& y) { return Positive<Float<PR>>(y, _pr); }

/*
template<class F> inline FloatFactory<PR> factory(Approximation<F> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F> inline FloatFactory<PR> factory(LowerBound<F> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F> inline FloatFactory<PR> factory(UpperBound<F> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F> inline FloatFactory<PR> factory(Bounds<F> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F, class FE> inline FloatFactory<PR> factory(FloatBall<PR,PRE> const& flt) { return FloatFactory<PR>(flt.precision()); }
template<class F> inline FloatFactory<PR> factory(F const& flt) { return FloatFactory<PR>(flt.precision()); }

template<class F> inline Approximation<F> FloatFactory<PR>::(Number<ApproximateTag> const& y) { return Approximation<F>(y,_pr); }
template<class F> inline LowerBound<F> FloatFactory<PR>::(Number<ValidatedLowerTag> const& y) { return LowerBound<F>(y,_pr); }
template<class F> inline UpperBound<F> FloatFactory<PR>::(Number<ValidatedUpperTag> const& y) { return UpperBound<F>(y,_pr); }

template<class F> inline Bounds<F> FloatFactory<PR>::(Number<ValidatedTag> const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Number<EffectiveTag> const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Number<ExactTag> const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Real const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Rational const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Dyadic const& y) { return Bounds<F>(y,_pr); }
template<class F> inline Bounds<F> FloatFactory<PR>::(Integer const& y) { return Bounds<F>(y,_pr); }

template<class F> inline F FloatFactory<PR>::(Dyadic const& y, ExactTag) { return F(y,_pr); }
template<class F> inline F FloatFactory<PR>::(Integer const& y, ExactTag) { return F(y,_pr); }

template<class F> inline template<IsBuiltinSignedIntegral N> F FloatFactory<PR>::(N const& y) { return F(y,_pr); }
template<class F> inline template<IsBuiltinUnsignedIntegral M> Positive<F> FloatFactory<PR>::(M const& y) { return Positive<F>(y,_pr); }
template<class F> inline template<IsBuiltinFloatingPoint D> Approximation<F> FloatFactory<PR>::(D const& y) { return Approximation<F>(F(y,_pr)); }
*/

} // namespace Ariadne

#endif
