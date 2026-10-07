/***************************************************************************
 *            test_float_approximation.cpp
 *
 *  Copyright  2006-20  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "numeric/builtin.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/float_approximation.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_error.hpp"

#include "utility/test.hpp"
#include "test_floats.hpp"

using namespace Ariadne;
using namespace std;

namespace {

template<class PR>
constexpr bool check_float_approximation_concept()
{
    using PRE=DoublePrecision;
    using A=FloatApproximation<PR>;
    return requires(PR pr, PRE pre, Nat m, Int n, double ad, ExactDouble ed, Integer z, Dyadic w, Rational q, A ax, A rx) {
        A(0u,pr); A(0,pr); A(0.0,pr); A(Integer(0),pr); A(Dyadic(0),pr); A(Rational(0),pr);
        rx=A(RawFloat<PR>(pr)); rx=A(Float<PR>(pr)); A(FloatBall<PR,PRE>(pr,pre)); A(FloatBounds<PR>(pr));
        rx=A(FloatUpperBound<PR>(pr)); rx=A(FloatLowerBound<PR>(pr)); rx=A(A(pr));
        rx=m; rx=n; rx=ad; rx=ed; rx=z; rx=w; rx=q;
        rx=nul(ax); rx=pos(ax); rx=neg(ax); rx=hlf(ax); rx=sqr(ax); rx=rec(ax);
        rx=operator+(ax); rx=operator-(ax); rx=operator+(ax,ax); rx=operator-(ax,ax); rx=operator*(ax,ax); rx=operator/(ax,ax);
        rx=add(ax,ax); rx=sub(ax,ax); rx=mul(ax,ax); rx=div(ax,ax); rx=pow(ax,m); rx=pow(ax,n);
        rx=max(ax,ax); rx=min(ax,ax); rx=abs(ax);
        rx=sqrt(ax); rx=exp(ax); rx=log(ax); rx=sin(ax); rx=cos(ax); rx=tan(ax); rx=atan(ax);
    };
}

static_assert(check_float_approximation_concept<DoublePrecision>());
static_assert(check_float_approximation_concept<MultiplePrecision>());

} // namespace


template<class PR>
class TestFloatApproximation
    : public TestFloats<PR>
{
    typedef FloatType<ApproximateTag,PR> FloatApproximationType;

    using TestFloats<PR>::m,TestFloats<PR>::n,TestFloats<PR>::ad,TestFloats<PR>::ed;
    using TestFloats<PR>::z,TestFloats<PR>::w,TestFloats<PR>::d,TestFloats<PR>::q;
  private:
    PR precision;
  public:
    TestFloatApproximation(PR prec) : precision(prec) { };
    Void test();
  private:
    Void test_conversions();
    Void test_arithmetic();
    Void test_comparison();
};


template<class PR> Void
TestFloatApproximation<PR>::test()
{
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_comparison());
}

template<class PR> Void
TestFloatApproximation<PR>::test_conversions()
{
    PR pr=precision;
    ARIADNE_TEST_EQUALS(cast_integer(FloatApproximation<PR>(Dyadic(-3,1u),pr)),Integer(-2));
    ARIADNE_TEST_EQUALS(cast_integer(FloatApproximation<PR>(Dyadic(2),pr)),Integer(2));
    ARIADNE_TEST_EQUALS(cast_integer(FloatApproximation<PR>(Dyadic(3,1u),pr)),Integer(2));
    ARIADNE_TEST_EQUALS(cast_integer(FloatApproximation<PR>(Dyadic(7,2u),pr)),Integer(2));

    FloatLowerBound<PR> lower(Dyadic(5,2u),pr);
    FloatUpperBound<PR> upper(Dyadic(5,2u),pr);
    FloatError<PR> error(1u,pr);
    ARIADNE_TEST_EQUALS(FloatApproximation<PR>(lower).raw(),lower.raw());
    ARIADNE_TEST_EQUALS(FloatApproximation<PR>(upper).raw(),upper.raw());
    ARIADNE_TEST_EQUALS(FloatApproximation<PR>(error).raw(),error.raw());
}

template<class PR> Void
TestFloatApproximation<PR>::test_arithmetic()
{
    RawFloat<PR> third(Rational(1,3),near,precision);
    RawFloat<PR> fifth(Rational(1,5),near,precision);
    RawFloat<PR> seventh(Rational(1,7),near,precision);

    ARIADNE_TEST_SAME(max(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(third));
    ARIADNE_TEST_SAME(min(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(fifth));
    ARIADNE_TEST_SAME(abs(FloatApproximation<PR>(neg(third))),FloatApproximation<PR>(third));

    ARIADNE_TEST_SAME(mag(FloatApproximation<PR>(neg(third))),PositiveFloatApproximation<PR>(third));
    ARIADNE_TEST_SAME(mig(FloatApproximation<PR>(neg(third))),PositiveFloatApproximation<PR>(third));

    ARIADNE_TEST_SAME(operator+(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(add(near,third,fifth)));
    ARIADNE_TEST_SAME(operator-(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(sub(near,third,fifth)));
    ARIADNE_TEST_SAME(operator*(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(mul(near,third,fifth)));
    ARIADNE_TEST_SAME(operator/(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(div(near,third,fifth)));

    ARIADNE_TEST_SAME(nul(FloatApproximation<PR>(third)),FloatApproximation<PR>(nul(third)));
    ARIADNE_TEST_SAME(pos(FloatApproximation<PR>(third)),FloatApproximation<PR>(pos(third)));
    ARIADNE_TEST_SAME(neg(FloatApproximation<PR>(third)),FloatApproximation<PR>(neg(third)));
    ARIADNE_TEST_SAME(hlf(FloatApproximation<PR>(third)),FloatApproximation<PR>(hlf(third)));
    ARIADNE_TEST_SAME(add(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(add(near,third,fifth)));
    ARIADNE_TEST_SAME(sub(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(sub(near,third,fifth)));
    ARIADNE_TEST_SAME(mul(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(mul(near,third,fifth)));
    ARIADNE_TEST_SAME(div(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(div(near,third,fifth)));
    ARIADNE_TEST_SAME(fma(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth),FloatApproximation<PR>(seventh)),
                          FloatApproximation<PR>(fma(near,third,fifth,seventh)));
    ARIADNE_TEST_SAME(pow(FloatApproximation<PR>(fifth),3),FloatApproximation<PR>(pow(near,fifth,3)));
}

template<class PR> Void
TestFloatApproximation<PR>::test_comparison()
{
    PR pr=precision;
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)==FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)!=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)< FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)> FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)<=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-2,pr)>=FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));

    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)==FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)!=FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)< FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)> FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)<=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(FloatApproximation<PR>(-1,pr)>=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
}


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestFloatApproximation<DoublePrecision>(dp).test();
    TestFloatApproximation<MultiplePrecision>(MultiplePrecision(128_bits)).test();

    return ARIADNE_TEST_FAILURES;
}

