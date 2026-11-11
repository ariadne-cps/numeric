/***************************************************************************
 *            numeric/float_factory.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file numeric/float_factory.hpp
 *  @brief Factories for creating floating-point numbers.
 */

#ifndef ARIADNE_FLOAT_FACTORY_HPP
#define ARIADNE_FLOAT_FACTORY_HPP

#include "float.decl.hpp"
#include "number.decl.hpp"

namespace Ariadne {

template<class PR>
class FloatFactory {
  protected:
    PR _pr;
    typedef RawFloatType<PR> F;

  public:
    typedef PR PrecisionType;
    typedef PR CharacteristicsType;
    FloatFactory(PR const& pr) : _pr(pr) {}
    PR precision() const { return this->_pr; }
    PR characteristics() const { return this->_pr; }

  public:
    FloatApproximation<PR> create(const ApproximateNumber& y);
    FloatLowerBound<PR> create(const ValidatedLowerNumber& y);
    FloatUpperBound<PR> create(const ValidatedUpperNumber& y);
    FloatBounds<PR> create(const ValidatedNumber& y);
    FloatBounds<PR> create(const EffectiveNumber& y);
    FloatBounds<PR> create(const ExactNumber& y);
    FloatBounds<PR> create(const Real& y);
    FloatBounds<PR> create(const Rational& y);
    FloatBounds<PR> create(const Dyadic& y);
    FloatBounds<PR> create(const Integer& y);
    Float<PR> create(const ExactDouble& y);
    Float<PR> create(const Dyadic& y, ExactTag);
    Float<PR> create(const Integer& y, ExactTag);
    template<BuiltinSignedIntegral N>
    Float<PR> create(const N& y);
    template<BuiltinUnsignedIntegral M>
    PositiveFloat<PR> create(const M& y);
    template<BuiltinFloatingPoint D>
    FloatApproximation<PR> create(const D& y);
    PositiveFloatApproximation<PR> create(const PositiveApproximateNumber& y);
    PositiveFloatLowerBound<PR> create(const PositiveValidatedLowerNumber& y);
    PositiveFloatUpperBound<PR> create(const PositiveValidatedUpperNumber& y);
    PositiveFloatBounds<PR> create(const PositiveValidatedNumber& y);
};

template<class PR, class PRE>
class FloatBallFactory : public FloatFactory<PR> {
    PRE _pre;

  public:
    FloatBallFactory(PR const& pr, PRE const& pre) : FloatFactory<PR>(pr), _pre(pre) {}
    using FloatFactory<PR>::create;
    FloatBall<PR, PRE> create(const Number<ValidatedTag>& y);
    FloatBall<PR, PRE> create(const Number<EffectiveTag>& y);
    FloatBall<PR, PRE> create(const Number<ExactTag>& y);
    FloatBall<PR, PRE> create(const Real& y);
    FloatBall<PR, PRE> create(const Rational& y);
    FloatBall<PR, PRE> create(const Dyadic& y);
    FloatBall<PR, PRE> create(const Integer& y);
    PositiveFloatBall<PR, PRE> create(const PositiveNumber<ValidatedTag>& y);
};

template<class Y, class PR>
using ConcreteType = decltype(declval<FloatFactory<PR>>().create(declval<Y>()));
template<class Y, class PR>
inline decltype(auto) make_float(const Y& y, PR pr) { return float_factory(pr).create(y); }

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
