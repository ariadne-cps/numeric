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
#include "numeric/floatmp.hpp"
#include "numeric/float_approximation.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_error.hpp"

#include "ariadne-testing.hpp"
#include "test_floats.hpp"

using namespace Ariadne;
using namespace std;

namespace {

template<class PR>
constexpr bool check_float_approximation_concept()
{
    using A=FloatApproximation<PR>;
    using F=RawFloat<PR>;
    return requires(PR pr, Nat m, Int n, double d, Integer z, Dyadic w, Rational q, A x, A y, F f) {
        A(pr); A(0u,pr); A(0,pr); A(d,pr); A(z,pr); A(w,pr); A(q,pr); A(f);
        x=m; x=n; x=d; x=z; x=w; x=q; x=f;
        +x; -x; x+y; x-y; x*y; x/y; x+=y; x-=y; x*=y; x/=y;
        nul(x); pos(x); neg(x); hlf(x); sqr(x); rec(x);
        add(x,y); sub(x,y); mul(x,y); div(x,y); pow(x,m); pow(x,n);
        max(x,y); min(x,y); abs(x);
        sqrt(x); exp(x); log(x); sin(x); cos(x); tan(x); atan(x);
        x.precision(); x.characteristics(); x.raw(); x.generic();
    };
}

static_assert(check_float_approximation_concept<DoublePrecision>());
static_assert(check_float_approximation_concept<MultiplePrecision>());

} // namespace


template<class PR>
class TestFloatApproximation
    : public TestFloats<PR>, public Ariadne::Testing::ContextFixture
{
    typedef FloatType<ApproximateTag,PR> FloatApproximationType;

    using TestFloats<PR>::m,TestFloats<PR>::n,TestFloats<PR>::ad,TestFloats<PR>::ed;
    using TestFloats<PR>::z,TestFloats<PR>::w,TestFloats<PR>::d,TestFloats<PR>::q;
  private:
    PR precision;
  public:
    TestFloatApproximation(PR prec, Ariadne::Testing::Context& ctx)
        : Ariadne::Testing::ContextFixture(ctx), precision(prec) { }
    Void test_header_api();
    Void test_conversions();
    Void test_arithmetic();
    Void test_comparison();
};


template<class PR> Void
TestFloatApproximation<PR>::test_header_api()
{
    using F=RawFloat<PR>;
    using A=FloatApproximation<PR>;

    PR pr=precision;

    A zero(pr);
    ARIADNE_TESTING_EXPECT_EQUAL(zero.raw(),F(0,pr));

    A from_approximate_double(ApproximateDouble(1.5),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(from_approximate_double.get_d(),1.5);

    A from_exact_double(ExactDouble(1.0),pr);
    A from_twoexp(TwoExp(1),pr);
    A from_integer(Integer(1),pr);
    A from_decimal(Decimal(String("1.5")),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(from_exact_double.raw(),F(1,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(from_twoexp.raw(),F(2,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(from_integer.raw(),F(1,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(from_decimal.get_d(),1.5);

    A copied(from_integer,pr);
    ARIADNE_TESTING_EXPECT_EQUAL(copied.raw(),from_integer.raw());

    if constexpr (Same<PR,DoublePrecision>) {
        // Preserve the cross-backend conversion using the native type.
        FloatMPApproximation other(1,MultiplePrecision(128_bits));
        A cross(other,pr);
        ARIADNE_TESTING_EXPECT_EQUAL(cross.raw(),F(1,pr));
    } else {
        FloatDPApproximation other(1,dp);
        A cross(other,pr);
        ARIADNE_TESTING_EXPECT_EQUAL(cross.raw(),F(1,pr));
    }

    FloatLowerBound<PR> lower(1,pr);
    FloatUpperBound<PR> upper(1,pr);
    FloatBounds<PR> bounds(1,1,pr);
    F raw_one(1,pr);

    A assigned(pr);
    assigned=lower;
    ARIADNE_TESTING_EXPECT_EQUAL(assigned.raw(),F(1,pr));
    assigned=upper;
    ARIADNE_TESTING_EXPECT_EQUAL(assigned.raw(),F(1,pr));
    assigned=bounds;
    ARIADNE_TESTING_EXPECT_EQUAL(assigned.raw(),F(1,pr));
    assigned=raw_one;
    ARIADNE_TESTING_EXPECT_EQUAL(assigned.raw(),F(1,pr));

    ApproximateNumber generic=assigned.generic();
    assigned=generic;
    ARIADNE_TESTING_EXPECT_EQUAL(assigned.raw(),F(1,pr));

    A created=assigned.create(generic);
    ARIADNE_TESTING_EXPECT_EQUAL(created.raw(),F(1,pr));

    ARIADNE_TESTING_EXPECT_EQUAL(assigned.characteristics(),pr);
    F raw_copy=static_cast<F>(assigned);
    ARIADNE_TESTING_EXPECT_EQUAL(raw_copy,F(1,pr));
    ApproximateDouble approximate_double=static_cast<ApproximateDouble>(assigned);
    ARIADNE_TESTING_EXPECT_EQUAL(approximate_double.get_d(),1.0);
    ARIADNE_TESTING_EXPECT_EQUAL(assigned.get_d(),1.0);

    A positive(1,pr);
    A negative(-1,pr);
    ARIADNE_TESTING_PROVE(tanh(positive).raw()>F(0,pr));
    ARIADNE_TESTING_PROVE(tanh(negative).raw()<F(0,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(asin(zero).raw(),F(0,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(acos(A(1,pr)).raw(),F(0,pr));

    ARIADNE_TESTING_EXPECT_SAME(assigned.pm(zero),assigned);

    A generic_two(2,pr);
    Integer generic_one(1);
    ARIADNE_TESTING_EXPECT_EQUAL(sub(generic_two,generic_one).raw(),F(1,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(mul(generic_two,generic_one).raw(),F(2,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(div(generic_two,generic_one).raw(),F(2,pr));

    ARIADNE_TESTING_EXPECT_EQUAL(add(generic_one,generic_two).raw(),F(3,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(sub(generic_one,generic_two).raw(),F(-1,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(mul(generic_one,generic_two).raw(),F(2,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(div(generic_one,generic_two).raw(),F(Dyadic(1,1u),pr));

    ARIADNE_TESTING_EXPECT_EQUAL(max(generic_two,generic_one).raw(),F(2,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(min(generic_two,generic_one).raw(),F(1,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(max(generic_one,generic_two).raw(),F(2,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(min(generic_one,generic_two).raw(),F(1,pr));

    ARIADNE_TESTING_EXPECT_SAME(generic_two!=generic_one,ApproximateKleenean(true));
    ARIADNE_TESTING_EXPECT_SAME(generic_one==generic_two,ApproximateKleenean(false));
    ARIADNE_TESTING_EXPECT_SAME(generic_one!=generic_two,ApproximateKleenean(true));
    ARIADNE_TESTING_EXPECT_SAME(generic_one<=generic_two,ApproximateKleenean(true));
    ARIADNE_TESTING_EXPECT_SAME(generic_one>=generic_two,ApproximateKleenean(false));
    ARIADNE_TESTING_EXPECT_SAME(generic_one>generic_two,ApproximateKleenean(false));
}

template<class PR> Void
TestFloatApproximation<PR>::test_conversions()
{
    PR pr=precision;
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatApproximation<PR>(Dyadic(-3,1u),pr)),Integer(-2));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatApproximation<PR>(Dyadic(2),pr)),Integer(2));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatApproximation<PR>(Dyadic(3,1u),pr)),Integer(2));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatApproximation<PR>(Dyadic(7,2u),pr)),Integer(2));

    FloatLowerBound<PR> lower(Dyadic(5,2u),pr);
    FloatUpperBound<PR> upper(Dyadic(5,2u),pr);
    FloatError<PR> error(1u,pr);
    ARIADNE_TESTING_EXPECT_EQUAL(FloatApproximation<PR>(lower).raw(),lower.raw());
    ARIADNE_TESTING_EXPECT_EQUAL(FloatApproximation<PR>(upper).raw(),upper.raw());
    ARIADNE_TESTING_EXPECT_EQUAL(FloatApproximation<PR>(error).raw(),error.raw());
}

template<class PR> Void
TestFloatApproximation<PR>::test_arithmetic()
{
    RawFloat<PR> third(Rational(1,3),near,precision);
    RawFloat<PR> fifth(Rational(1,5),near,precision);
    RawFloat<PR> seventh(Rational(1,7),near,precision);

    ARIADNE_TESTING_EXPECT_SAME(max(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(third));
    ARIADNE_TESTING_EXPECT_SAME(min(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(fifth));
    ARIADNE_TESTING_EXPECT_SAME(abs(FloatApproximation<PR>(neg(third))),FloatApproximation<PR>(third));

    ARIADNE_TESTING_EXPECT_SAME(mag(FloatApproximation<PR>(neg(third))),PositiveFloatApproximation<PR>(third));
    ARIADNE_TESTING_EXPECT_SAME(mig(FloatApproximation<PR>(neg(third))),PositiveFloatApproximation<PR>(third));

    ARIADNE_TESTING_EXPECT_SAME(operator+(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(add(near,third,fifth)));
    ARIADNE_TESTING_EXPECT_SAME(operator-(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(sub(near,third,fifth)));
    ARIADNE_TESTING_EXPECT_SAME(operator*(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(mul(near,third,fifth)));
    ARIADNE_TESTING_EXPECT_SAME(operator/(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(div(near,third,fifth)));

    ARIADNE_TESTING_EXPECT_SAME(nul(FloatApproximation<PR>(third)),FloatApproximation<PR>(nul(third)));
    ARIADNE_TESTING_EXPECT_SAME(pos(FloatApproximation<PR>(third)),FloatApproximation<PR>(pos(third)));
    ARIADNE_TESTING_EXPECT_SAME(neg(FloatApproximation<PR>(third)),FloatApproximation<PR>(neg(third)));
    ARIADNE_TESTING_EXPECT_SAME(hlf(FloatApproximation<PR>(third)),FloatApproximation<PR>(hlf(third)));
    ARIADNE_TESTING_EXPECT_SAME(add(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(add(near,third,fifth)));
    ARIADNE_TESTING_EXPECT_SAME(sub(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(sub(near,third,fifth)));
    ARIADNE_TESTING_EXPECT_SAME(mul(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(mul(near,third,fifth)));
    ARIADNE_TESTING_EXPECT_SAME(div(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth)),FloatApproximation<PR>(div(near,third,fifth)));
    ARIADNE_TESTING_EXPECT_SAME(fma(FloatApproximation<PR>(third),FloatApproximation<PR>(fifth),FloatApproximation<PR>(seventh)),
                          FloatApproximation<PR>(fma(near,third,fifth,seventh)));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatApproximation<PR>(fifth),3),FloatApproximation<PR>(pow(near,fifth,3)));
}

template<class PR> Void
TestFloatApproximation<PR>::test_comparison()
{
    PR pr=precision;
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-2,pr)==FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-2,pr)!=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-2,pr)< FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-2,pr)> FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-2,pr)<=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-2,pr)>=FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));

    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-1,pr)==FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-1,pr)!=FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-1,pr)< FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-1,pr)> FloatApproximation<PR>(-1,pr),ApproximateKleenean(false));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-1,pr)<=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
    ARIADNE_TESTING_EXPECT_SAME(FloatApproximation<PR>(-1,pr)>=FloatApproximation<PR>(-1,pr),ApproximateKleenean(true));
}


int main(int argc, char** argv) {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    Ariadne::Testing::Suite suite("float_approximation");
    auto register_precision = [&suite]<class PR>(std::string variant, PR precision) {
        using Fixture = TestFloatApproximation<PR>;
        auto add = [&](char const* phase, auto method) {
            suite.add(variant + " / " + phase,
                [precision, method](Ariadne::Testing::Context& ctx) {
                    Fixture fixture(precision,ctx);
                    (fixture.*method)();
                }, Ariadne::Testing::CaseMetadata{{"float-approximation"},variant});
        };
        add("header_api", &Fixture::test_header_api);
        add("conversions", &Fixture::test_conversions);
        add("arithmetic", &Fixture::test_arithmetic);
        add("comparison", &Fixture::test_comparison);
    };
    register_precision("FloatDP",dp);
    register_precision("FloatMP-128",MultiplePrecision(128_bits));
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
