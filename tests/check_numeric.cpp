/***************************************************************************
 *            check_numeric.cpp
 *
 *  Compile-time checks for cross-family numeric API relationships.
 ****************************************************************************/

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

template<class Expected, class OP, class A1, class A2s> struct CheckRow;

template<class OP, class A1>
struct CheckRow<Types<>,OP,A1,Types<>> {
    static constexpr bool value=true;
};

template<class E, class... ES, class OP, class A1, class A2, class... A2S>
struct CheckRow<Types<E,ES...>,OP,A1,Types<A2,A2S...>> {
    static_assert(Same<SafeType<OP,A1,A2>,E>);
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
using UpN=Number<ValidatedUpperTag>;
using LoN=Number<ValidatedLowerTag>;
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

static_assert(CheckMatrix<ExpectedWeakerTable,Plus,NumericTypes,NumericTypes>::value);
static_assert(CheckMatrix<ExpectedWeakerTable,Times,NumericTypes,NumericTypes>::value);

// The old checker encoded these as runtime checks. Keep the same conversion
// contract, but make it a compile-time property of the public API.
static_assert(Convertible<Nat,FloatDPError>);
static_assert((not Constructible<FloatDPError,Int>));
static_assert((Constructible<FloatDPError,Dbl> and not Convertible<Dbl,FloatDPError>));
static_assert(Convertible<Nat,Nat>);
static_assert(Convertible<Int,Nat>);
static_assert(Convertible<Dbl,Nat>);
static_assert((not Constructible<Nat,Z>));
static_assert((not Constructible<Int,Q>));
static_assert((not Constructible<Int,R>));
static_assert(Convertible<Nat,Int>);
static_assert(Convertible<Int,Int>);
static_assert(Convertible<Dbl,Int>);
static_assert((not Constructible<Int,Z>));
static_assert((not Constructible<Int,Q>));
static_assert((not Constructible<Int,R>));
static_assert(Convertible<Nat,Dbl>);
static_assert(Convertible<Int,Dbl>);
static_assert(Convertible<Dbl,Dbl>);
static_assert((not Constructible<Dbl,Z>));
static_assert((not Constructible<Dbl,Q>));
static_assert((not Constructible<Dbl,R>));
static_assert(Convertible<Nat,Z>);
static_assert(Convertible<Int,Z>);
static_assert((not Constructible<Z,Dbl>));
static_assert(Convertible<Z,Z>);
static_assert((not Constructible<Z,Q>));
static_assert((not Constructible<Z,R>));
static_assert(Convertible<Nat,Q>);
static_assert(Convertible<Int,Q>);
static_assert((Constructible<Q,Dbl> and not Convertible<Dbl,Q>));
static_assert(Convertible<Z,Q>);
static_assert(Convertible<Q,Q>);
static_assert((not Constructible<Q,R>));
static_assert(Convertible<Nat,R>);
static_assert(Convertible<Int,R>);
static_assert((Constructible<R,Dbl> and not Convertible<Dbl,R>));
static_assert(Convertible<Z,R>);
static_assert(Convertible<Q,R>);
static_assert(Convertible<R,R>);
static_assert(Convertible<ExN,ExN>);
static_assert((not Constructible<ExN,EfN>));
static_assert((not Constructible<ExN,VaN>));
static_assert((not Constructible<ExN,UpN>));
static_assert((not Constructible<ExN,LoN>));
static_assert((not Constructible<ExN,ApN>));
static_assert(Convertible<ExN,EfN>);
static_assert(Convertible<EfN,EfN>);
static_assert((not Constructible<EfN,VaN>));
static_assert((not Constructible<EfN,UpN>));
static_assert((not Constructible<EfN,LoN>));
static_assert((not Constructible<EfN,ApN>));
static_assert(Convertible<ExN,VaN>);
static_assert(Convertible<EfN,VaN>);
static_assert(Convertible<VaN,VaN>);
static_assert((not Constructible<VaN,UpN>));
static_assert((not Constructible<VaN,LoN>));
static_assert((not Constructible<VaN,ApN>));
static_assert(Convertible<ExN,UpN>);
static_assert(Convertible<EfN,UpN>);
static_assert(Convertible<VaN,UpN>);
static_assert(Convertible<UpN,UpN>);
static_assert((not Constructible<UpN,LoN>));
static_assert((not Constructible<UpN,ApN>));
static_assert(Convertible<ExN,LoN>);
static_assert(Convertible<EfN,LoN>);
static_assert(Convertible<VaN,LoN>);
static_assert((not Constructible<LoN,UpN>));
static_assert(Convertible<LoN,LoN>);
static_assert((not Constructible<LoN,ApN>));
static_assert(Convertible<ExN,ApN>);
static_assert(Convertible<EfN,ApN>);
static_assert(Convertible<VaN,ApN>);
static_assert(Convertible<UpN,ApN>);
static_assert(Convertible<LoN,ApN>);
static_assert(Convertible<ApN,ApN>);
static_assert(Convertible<Nat,ExN>);
static_assert(Convertible<Int,ExN>);
static_assert((not Constructible<ExN,Dbl>));
static_assert(Convertible<Z,ExN>);
static_assert(Convertible<Q,ExN>);
static_assert((not Constructible<ExN,R>));
static_assert(Convertible<Nat,EfN>);
static_assert(Convertible<Int,EfN>);
static_assert((not Constructible<EfN,Dbl>));
static_assert(Convertible<Z,EfN>);
static_assert(Convertible<Q,EfN>);
static_assert(Convertible<R,EfN>);
static_assert(Convertible<Nat,VaN>);
static_assert(Convertible<Int,VaN>);
static_assert((not Constructible<VaN,Dbl>));
static_assert(Convertible<Z,VaN>);
static_assert(Convertible<Q,VaN>);
static_assert(Convertible<R,VaN>);
static_assert(Convertible<Nat,UpN>);
static_assert(Convertible<Int,UpN>);
static_assert((not Constructible<UpN,Dbl>));
static_assert(Convertible<Z,UpN>);
static_assert(Convertible<Q,UpN>);
static_assert(Convertible<R,UpN>);
static_assert(Convertible<Nat,LoN>);
static_assert(Convertible<Int,LoN>);
static_assert((not Constructible<LoN,Dbl>));
static_assert(Convertible<Z,LoN>);
static_assert(Convertible<Q,LoN>);
static_assert(Convertible<R,LoN>);
static_assert(Convertible<Nat,ApN>);
static_assert(Convertible<Int,ApN>);
static_assert(Convertible<Dbl,ApN>);
static_assert(Convertible<Z,ApN>);
static_assert(Convertible<Q,ApN>);
static_assert(Convertible<R,ApN>);
static_assert((not Constructible<Nat,ExN>));
static_assert((not Constructible<Nat,EfN>));
static_assert((not Constructible<Nat,VaN>));
static_assert((not Constructible<Nat,UpN>));
static_assert((not Constructible<Nat,LoN>));
static_assert((not Constructible<Nat,ApN>));
static_assert((not Constructible<Int,ExN>));
static_assert((not Constructible<Int,EfN>));
static_assert((not Constructible<Int,VaN>));
static_assert((not Constructible<Int,UpN>));
static_assert((not Constructible<Int,LoN>));
static_assert((not Constructible<Int,ApN>));
static_assert((not Constructible<Dbl,ExN>));
static_assert((not Constructible<Dbl,EfN>));
static_assert((not Constructible<Dbl,VaN>));
static_assert((not Constructible<Dbl,UpN>));
static_assert((not Constructible<Dbl,LoN>));
static_assert((not Constructible<Dbl,ApN>));
static_assert((not Constructible<Z,ExN>));
static_assert((not Constructible<Z,EfN>));
static_assert((not Constructible<Z,VaN>));
static_assert((not Constructible<Z,UpN>));
static_assert((not Constructible<Z,LoN>));
static_assert((not Constructible<Z,ApN>));
static_assert((not Constructible<Q,ExN>));
static_assert((not Constructible<Q,EfN>));
static_assert((not Constructible<Q,VaN>));
static_assert((not Constructible<Q,UpN>));
static_assert((not Constructible<Q,LoN>));
static_assert((not Constructible<Q,ApN>));
static_assert((not Constructible<R,ExN>));
static_assert((not Constructible<R,EfN>));
static_assert((not Constructible<R,VaN>));
static_assert((not Constructible<R,UpN>));
static_assert((not Constructible<R,LoN>));
static_assert((not Constructible<R,ApN>));
static_assert(Convertible<ExF,ExF>);
static_assert((not Constructible<ExF,MeF>));
static_assert((not Constructible<ExF,BoF>));
static_assert((not Constructible<ExF,UpF>));
static_assert((not Constructible<ExF,LoF>));
static_assert((not Constructible<ExF,ApF>));
static_assert(Convertible<ExF,MeF>);
static_assert(Convertible<MeF,MeF>);
static_assert(Convertible<BoF,MeF>);
static_assert((not Constructible<MeF,UpF>));
static_assert((not Constructible<MeF,LoF>));
static_assert((not Constructible<MeF,ApF>));
static_assert(Convertible<ExF,BoF>);
static_assert(Convertible<MeF,BoF>);
static_assert(Convertible<BoF,BoF>);
static_assert((not Constructible<BoF,UpF>));
static_assert((not Constructible<BoF,LoF>));
static_assert((not Constructible<BoF,ApF>));
static_assert(Convertible<ExF,UpF>);
static_assert(Convertible<MeF,UpF>);
static_assert(Convertible<BoF,UpF>);
static_assert(Convertible<UpF,UpF>);
static_assert((not Constructible<UpF,LoF>));
static_assert((not Constructible<UpF,ApF>));
static_assert(Convertible<ExF,LoF>);
static_assert(Convertible<MeF,LoF>);
static_assert(Convertible<BoF,LoF>);
static_assert((not Constructible<LoF,UpF>));
static_assert(Convertible<LoF,LoF>);
static_assert((not Constructible<LoF,ApF>));
static_assert(Convertible<ExF,ApF>);
static_assert(Convertible<MeF,ApF>);
static_assert(Convertible<BoF,ApF>);
static_assert(Convertible<UpF,ApF>);
static_assert(Convertible<LoF,ApF>);
static_assert(Convertible<ApF,ApF>);
static_assert((not Constructible<ExF,ExN>));
static_assert((not Constructible<ExF,EfN>));
static_assert((not Constructible<ExF,VaN>));
static_assert((not Constructible<ExF,UpN>));
static_assert((not Constructible<ExF,LoN>));
static_assert((not Constructible<ExF,ApN>));
static_assert((Constructible<MeF,ExN> and not Convertible<ExN,MeF>));
static_assert((Constructible<MeF,EfN> and not Convertible<EfN,MeF>));
static_assert((Constructible<MeF,VaN> and not Convertible<VaN,MeF>));
static_assert((not Constructible<MeF,UpN>));
static_assert((not Constructible<MeF,LoN>));
static_assert((not Constructible<MeF,ApN>));
static_assert((Constructible<BoF,ExN> and not Convertible<ExN,BoF>));
static_assert((Constructible<BoF,EfN> and not Convertible<EfN,BoF>));
static_assert((Constructible<BoF,VaN> and not Convertible<VaN,BoF>));
static_assert((not Constructible<BoF,UpN>));
static_assert((not Constructible<BoF,LoN>));
static_assert((not Constructible<BoF,ApN>));
static_assert((Constructible<UpF,ExN> and not Convertible<ExN,UpF>));
static_assert((Constructible<UpF,EfN> and not Convertible<EfN,UpF>));
static_assert((Constructible<UpF,VaN> and not Convertible<VaN,UpF>));
static_assert((Constructible<UpF,UpN> and not Convertible<UpN,UpF>));
static_assert((not Constructible<UpF,LoN>));
static_assert((not Constructible<UpF,ApN>));
static_assert((Constructible<LoF,ExN> and not Convertible<ExN,LoF>));
static_assert((Constructible<LoF,EfN> and not Convertible<EfN,LoF>));
static_assert((Constructible<LoF,VaN> and not Convertible<VaN,LoF>));
static_assert((not Constructible<LoF,UpN>));
static_assert((Constructible<LoF,LoN> and not Convertible<LoN,LoF>));
static_assert((not Constructible<LoF,ApN>));
static_assert((Constructible<ApF,ExN> and not Convertible<ExN,ApF>));
static_assert((Constructible<ApF,EfN> and not Convertible<EfN,ApF>));
static_assert((Constructible<ApF,VaN> and not Convertible<VaN,ApF>));
static_assert((Constructible<ApF,UpN> and not Convertible<UpN,ApF>));
static_assert((Constructible<ApF,LoN> and not Convertible<LoN,ApF>));
static_assert((Constructible<ApF,ApN> and not Convertible<ApN,ApF>));
static_assert(Convertible<ExF,ExN>);
static_assert((not Constructible<ExN,MeF>));
static_assert((not Constructible<ExN,BoF>));
static_assert((not Constructible<ExN,UpF>));
static_assert((not Constructible<ExN,LoF>));
static_assert((not Constructible<ExN,ApF>));
static_assert(Convertible<ExF,EfN>);
static_assert((not Constructible<EfN,MeF>));
static_assert((not Constructible<EfN,BoF>));
static_assert((not Constructible<EfN,UpF>));
static_assert((not Constructible<EfN,LoF>));
static_assert((not Constructible<EfN,ApF>));
static_assert(Convertible<ExF,VaN>);
static_assert(Convertible<MeF,VaN>);
static_assert(Convertible<BoF,VaN>);
static_assert((not Constructible<VaN,UpF>));
static_assert((not Constructible<VaN,LoF>));
static_assert((not Constructible<VaN,ApF>));
static_assert(Convertible<ExF,UpN>);
static_assert(Convertible<MeF,UpN>);
static_assert(Convertible<BoF,UpN>);
static_assert(Convertible<UpF,UpN>);
static_assert((not Constructible<UpN,LoF>));
static_assert((not Constructible<UpN,ApF>));
static_assert(Convertible<ExF,LoN>);
static_assert(Convertible<MeF,LoN>);
static_assert(Convertible<BoF,LoN>);
static_assert((not Constructible<LoN,UpF>));
static_assert(Convertible<LoF,LoN>);
static_assert((not Constructible<LoN,ApF>));
static_assert(Convertible<ExF,ApN>);
static_assert(Convertible<MeF,ApN>);
static_assert(Convertible<BoF,ApN>);
static_assert(Convertible<UpF,ApN>);
static_assert(Convertible<LoF,ApN>);
static_assert(Convertible<ApF,ApN>);
static_assert((Constructible<ExF,Dbl> and not Convertible<Dbl,ExF>));
static_assert(Convertible<Int,ExF>);
static_assert((Constructible<ExF,Z> and not Convertible<Z,ExF>));
static_assert((not Constructible<ExF,Q>));
static_assert((not Constructible<ExF,R>));
static_assert((Constructible<MeF,Dbl> and not Convertible<Dbl,MeF>));
static_assert(Convertible<Int,MeF>);
static_assert((Constructible<MeF,Z> and not Convertible<Z,MeF>));
static_assert((Constructible<MeF,Q> and not Convertible<Q,MeF>));
static_assert((Constructible<MeF,R> and not Convertible<R,MeF>));
static_assert((Constructible<BoF,Dbl> and not Convertible<Dbl,BoF>));
static_assert(Convertible<Int,BoF>);
static_assert((Constructible<BoF,Z> and not Convertible<Z,BoF>));
static_assert((Constructible<BoF,Q> and not Convertible<Q,BoF>));
static_assert((Constructible<BoF,R> and not Convertible<R,BoF>));
static_assert((Constructible<UpF,Dbl> and not Convertible<Dbl,UpF>));
static_assert(Convertible<Int,UpF>);
static_assert((Constructible<UpF,Z> and not Convertible<Z,UpF>));
static_assert((Constructible<UpF,Q> and not Convertible<Q,UpF>));
static_assert((Constructible<UpF,R> and not Convertible<R,UpF>));
static_assert((Constructible<LoF,Dbl> and not Convertible<Dbl,LoF>));
static_assert(Convertible<Int,LoF>);
static_assert((Constructible<LoF,Z> and not Convertible<Z,LoF>));
static_assert((Constructible<LoF,Q> and not Convertible<Q,LoF>));
static_assert((Constructible<LoF,R> and not Convertible<R,LoF>));
static_assert(Convertible<Dbl,ApF>);
static_assert(Convertible<Int,ApF>);
static_assert((Constructible<ApF,Z> and not Convertible<Z,ApF>));
static_assert((Constructible<ApF,Q> and not Convertible<Q,ApF>));
static_assert((Constructible<ApF,R> and not Convertible<R,ApF>));
static_assert((not Constructible<Dbl,ExF>));
static_assert((not Constructible<Dbl,MeF>));
static_assert((not Constructible<Dbl,BoF>));
static_assert((not Constructible<Dbl,UpF>));
static_assert((not Constructible<Dbl,LoF>));
static_assert((not Constructible<Dbl,ApF>));
static_assert((not Constructible<Int,ExF>));
static_assert((not Constructible<Int,MeF>));
static_assert((not Constructible<Int,BoF>));
static_assert((not Constructible<Int,UpF>));
static_assert((not Constructible<Int,LoF>));
static_assert((not Constructible<Int,ApF>));
static_assert((not Constructible<Z,ExF>));
static_assert((not Constructible<Z,MeF>));
static_assert((not Constructible<Z,BoF>));
static_assert((not Constructible<Z,UpF>));
static_assert((not Constructible<Z,LoF>));
static_assert((not Constructible<Z,ApF>));
static_assert((Constructible<Q,ExF> and not Convertible<ExF,Q>));
static_assert((not Constructible<Q,MeF>));
static_assert((not Constructible<Q,BoF>));
static_assert((not Constructible<Q,UpF>));
static_assert((not Constructible<Q,LoF>));
static_assert((not Constructible<Q,ApF>));
static_assert((Constructible<R,ExF> and not Convertible<ExF,R>));
static_assert((not Constructible<R,MeF>));
static_assert((not Constructible<R,BoF>));
static_assert((not Constructible<R,UpF>));
static_assert((not Constructible<R,LoF>));
static_assert((not Constructible<R,ApF>));

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
