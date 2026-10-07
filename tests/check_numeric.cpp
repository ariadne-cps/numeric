/***************************************************************************
 *            check_numeric.cpp
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

#include <type_traits>
#include <utility>

#include "foundation/logical.hpp"

#include "numeric/concepts.hpp"
#include "numeric/operators.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include "numeric/real.hpp"
#include "numeric/number.hpp"
#include "numeric/floats.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"

using namespace Ariadne;

namespace {

template<class... TS> struct Types { };

struct NoResult { };

template<class OP, class A1, class A2, class=void>
struct SafeTypeTrait { using Type=NoResult; };

template<class OP, class A1, class A2>
struct SafeTypeTrait<OP,A1,A2,std::void_t<decltype(OP()(std::declval<A1>(),std::declval<A2>()))>> {
    using Type=decltype(OP()(std::declval<A1>(),std::declval<A2>()));
};

template<class OP, class A1, class A2>
using SafeType=typename SafeTypeTrait<OP,A1,A2>::Type;

template<class Expected, class OP, class A1, class A2s> struct CheckRow;

template<class OP, class A1>
struct CheckRow<Types<>,OP,A1,Types<>> {
    static constexpr bool value=true;
};

template<class E, class... ES, class OP, class A1, class A2, class... A2S>
struct CheckRow<Types<E,ES...>,OP,A1,Types<A2,A2S...>> {
    using R=SafeType<OP,A1,A2>;
    static_assert(Same<R,NoResult> or Same<R,E>);
    static constexpr bool value=CheckRow<Types<ES...>,OP,A1,Types<A2S...>>::value;
};

template<class ExpectedRows, class OP, class A1s, class A2s> struct CheckMatrix;

template<class OP, class A2s>
struct CheckMatrix<Types<>,OP,Types<>,A2s> {
    static constexpr bool value=true;
};

template<class Row, class... Rows, class OP, class A1, class... A1S, class A2s>
struct CheckMatrix<Types<Row,Rows...>,OP,Types<A1,A1S...>,A2s> {
    static_assert(CheckRow<Row,OP,A1,A2s>::value);
    static constexpr bool value=CheckMatrix<Types<Rows...>,OP,Types<A1S...>,A2s>::value;
};

using B=bool;
using Nat=unsigned int;
using Int=int;
using Dbl=double;
using Z=Integer;
using Q=Rational;
using R=Real;

using ExN=Number<ExactTag>;
using EfN=Number<EffectiveTag>;
using VaN=Number<ValidatedTag>;
using UpN=ValidatedUpperNumber;
using LoN=ValidatedLowerNumber;
using ApN=Number<ApproximateTag>;

using ExF=FloatDP;
using MeF=FloatDPBall;
using BoF=FloatDPBounds;
using UpF=FloatDPUpperBound;
using LoF=FloatDPLowerBound;
using ApF=FloatDPApproximation;

using EfF=decltype(std::declval<ExN>()+std::declval<ExF>());
using VaF=decltype(std::declval<MeF>()+std::declval<BoF>());
using ApD=decltype(std::declval<double>()+std::declval<ApN>());

using NumericTypes=Types<
    Nat,Int,Dbl,Z,Q,R,
    ExN,EfN,VaN,UpN,LoN,ApN,
    ExF,MeF,BoF,UpF,LoF,ApF>;

using ExpectedWeakerTable=
    Types<
        Types<Nat,Nat,Dbl, Z,Q,R, ExN,EfN,VaN,UpN,LoN,ApN, EfF,MeF,BoF,UpF,LoF,ApF>,
        Types<Nat,Int,Dbl, Z,Q,R, ExN,EfN,VaN,UpN,LoN,ApN, EfF,MeF,BoF,UpF,LoF,ApF>,
        Types<Dbl,Dbl,Dbl, ApD,ApD,ApD, ApD,ApD,ApD,ApD,ApD,ApD, ApF,ApF,ApF,ApF,ApF,ApF>,
        Types<Z,Z,ApD, Z,Q,R, ExN,EfN,VaN,UpN,LoN,ApN, EfF,MeF,BoF,UpF,LoF,ApF>,
        Types<Q,Q,ApD, Q,Q,R, ExN,EfN,VaN,UpN,LoN,ApN, EfF,MeF,BoF,UpF,LoF,ApF>,
        Types<R,R,ApD, R,R,R, EfN,EfN,VaN,UpN,LoN,ApN, EfF,MeF,BoF,UpF,LoF,ApF>,
        Types<ExN,ExN,ApN, ExN,ExN,EfN, ExN,EfN,VaN,UpN,LoN,ApN, EfF,MeF,BoF,UpF,LoF,ApF>,
        Types<EfN,EfN,ApN, EfN,EfN,EfN, EfN,EfN,VaN,UpN,LoN,ApN, EfF,MeF,BoF,UpF,LoF,ApF>,
        Types<VaN,VaN,ApN, VaN,VaN,VaN, VaN,VaN,VaN,UpN,LoN,ApN, EfF,MeF,BoF,UpF,LoF,ApF>,
        Types<UpN,UpN,ApN, UpN,UpN,UpN, UpN,UpN,UpN,UpN,ApN,ApN, UpF,UpF,UpF,UpF,ApF,ApF>,
        Types<LoN,LoN,ApN, LoN,LoN,LoN, LoN,LoN,LoN,ApN,LoN,ApN, LoF,LoF,LoF,ApF,LoF,ApF>,
        Types<ApN,ApN,ApN, ApN,ApN,ApN, ApN,ApN,ApN,ApN,ApN,ApN, ApF,ApF,ApF,ApF,ApF,ApF>,
        Types<EfF,EfF,ApF, EfF,EfF,EfF, EfF,EfF,EfF,UpF,LoF,ApF, EfF,MeF,BoF,UpF,LoF,ApF>,
        Types<MeF,MeF,ApF, MeF,MeF,MeF, MeF,MeF,MeF,UpF,LoF,ApF, MeF,MeF,VaF,UpF,LoF,ApF>,
        Types<BoF,BoF,ApF, BoF,BoF,BoF, BoF,BoF,BoF,UpF,LoF,ApF, BoF,VaF,BoF,UpF,LoF,ApF>,
        Types<UpF,UpF,ApF, UpF,UpF,UpF, UpF,UpF,UpF,UpF,ApF,ApF, UpF,UpF,UpF,UpF,ApF,ApF>,
        Types<LoF,LoF,ApF, LoF,LoF,LoF, LoF,LoF,LoF,ApF,LoF,ApF, LoF,LoF,LoF,ApF,LoF,ApF>,
        Types<ApF,ApF,ApF, ApF,ApF,ApF, ApF,ApF,ApF,ApF,ApF,ApF, ApF,ApF,ApF,ApF,ApF,ApF>
    >;

// For supported mixed operations, verify the result type predicted by the
// weakening table. Unsupported combinations are not treated as API contracts.
static_assert(CheckMatrix<ExpectedWeakerTable,Plus,NumericTypes,NumericTypes>::value);
static_assert(CheckMatrix<ExpectedWeakerTable,Times,NumericTypes,NumericTypes>::value);

// Conversion contracts that are part of the current public API.
static_assert(Convertible<Nat,Z>);
static_assert(Convertible<Int,Z>);
static_assert(Convertible<Z,Q>);
static_assert(Convertible<Q,R>);
static_assert(not Constructible<Q,Dbl>);
static_assert(not Constructible<R,Dbl>);
static_assert(Convertible<ExF,ExN>);
static_assert(Convertible<BoF,VaN>);
static_assert(Convertible<UpF,UpN>);
static_assert(Convertible<LoF,LoN>);
static_assert(Convertible<ApF,ApN>);
static_assert(Constructible<UpN,Dyadic>);
static_assert(Constructible<LoN,Dyadic>);
static_assert(Constructible<ApN,Dbl>);

// Directed subtraction is not the same weakening rule as addition. Check the
// important cross-family cases explicitly.
static_assert(Same<DifferenceType<ExN,UpN>,LoN>);
static_assert(Same<DifferenceType<LoN,UpN>,LoN>);
static_assert(Same<DifferenceType<LoN,VaN>,LoN>);
static_assert(Same<DifferenceType<Integer,UpF>,LoF>);
static_assert(Same<DifferenceType<Dyadic,UpF>,LoF>);
static_assert(Same<DifferenceType<UpF,Integer>,UpF>);
static_assert(Same<DifferenceType<UpF,Dyadic>,UpF>);

// Representative mixed-operation result types that historically regressed.
static_assert(Same<decltype(Add()(std::declval<UpF>(),std::declval<UpF>())),UpF>);
static_assert(Same<decltype(Sub()(std::declval<UpF>(),std::declval<LoF>())),UpF>);
static_assert(Same<decltype(BinaryElementaryOperator(Add())(std::declval<UpF>(),std::declval<UpF>())),ApF>);
static_assert(Same<decltype(BinaryElementaryOperator(Sub())(std::declval<UpF>(),std::declval<LoF>())),ApF>);

} // namespace
