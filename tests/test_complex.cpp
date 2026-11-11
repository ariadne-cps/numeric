/***************************************************************************
 *            test_complex.cpp
 *
 *  Copyright  2019-20  Pieter Collins
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

#include "foundation/logical.hpp"
#include "numeric/builtin.hpp"
#include "numeric/complex.hpp"
#include "numeric/decimal.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/float_approximation.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_error.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include "numeric/real.hpp"
#include "utility/module.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

template<class X>
constexpr bool check_concept() {
    return requires(X x, Complex<X> z, Complex<X> w) {
        Complex<X>(x);
        Complex<X>(x, x);
        +z;
        -z;
        z + w;
        z - w;
        z * w;
        z / w;
        add(z, w);
        sub(z, w);
        mul(z, w);
        div(z, w);
        sqr(z);
        rec(z);
        pow(z, 2);
        pow(z, 2u);
        sqrt(z);
        exp(z);
        log(z);
        abs(z);
        mag(z);
        mig(z);
        arg(z);
        conj(z);
        dist(z, w);
        z == w;
        z != w;
    };
}

static_assert(check_concept<FloatDPBounds>());

} // namespace

template<class X>
class TestComplex : public Ariadne::Testing::ContextFixture {
    X _one;
    typedef Complex<X> ComplexType;

  public:
    TestComplex(Ariadne::Testing::Context& context, const X& one)
        : ContextFixture(context), _one(one) {}

  public:
    void test_constructors();
    void test_conversions();
    void test_arithmetic();
    void test_polar();
    void test_transcendental();
    void test_comparison();
    void test_fourier();
};

template<class X>
void TestComplex<X>::test_conversions() {
    //    auto one_dp=one.get(dp);
    //    ARIADNE_TESTING_INSPECT("one_dp", one_dp);
}

template<class X>
void TestComplex<X>::test_constructors() {
    ARIADNE_TESTING_EXPECT_NO_THROW(Positive<Real>(Real(1)));
    /*
        ARIADNE_TEST_CONSTRUCT(ComplexType,zz,);
        ARIADNE_TESTING_EXPECT_EQUAL(zz.get(dp),0);
        ARIADNE_TEST_CONSTRUCT(ComplexType,zr,(1));
        ARIADNE_TESTING_EXPECT_EQUAL(zr.get(dp),1);
        ARIADNE_TEST_CONSTRUCT(ComplexType,zc,(4,3));
        ARIADNE_TESTING_EXPECT_EQUAL(zc.get(dp).real_part(),4);
    */
}

template<class X>
void TestComplex<X>::test_arithmetic() {
    FloatDPApproximation::set_output_places(18);
    X x(2.5_dyadic * _one);
    X y(4.0_dyadic * _one);

    using Constants::i;

    ARIADNE_TESTING_EXPECT_EQUAL(x + i * y, ComplexType(x, y));

    ARIADNE_TESTING_EXPECT_EQUAL(x, 2.5_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(y, 4.0_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(+x, 2.5_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(-x, -2.5_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(x + y, 6.5_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(x - y, -1.5_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(x * y, 10.0_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(x / y, 0.625_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(add(x, y), 6.5_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(sub(x, y), -1.5_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(mul(x, y), 10.0_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(div(x, y), 0.625_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(pow(x, 3u), 15.625_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(pow(x, 3), 15.625_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(pos(x), +2.5_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(neg(x), -2.5_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(x), 1.25_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(sqr(x), 6.25_dy);
    ARIADNE_TESTING_EXPECT_EQUAL(rec(y), 0.25_dy);
}

template<class X>
void TestComplex<X>::test_transcendental() {
    using Constants::i;
    Dyadic eps{FloatDP::eps(dp)};
    Dyadic tol = 1024 * eps;
    Complex<X> x(3.0_dyadic * _one, 4.0_dyadic * _one);

    ARIADNE_TESTING_EXPECT_WITHIN(sqrt(sqr(x)), x, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(sqr(sqrt(x)), x, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(log(exp(x)) + 2 * pi * i, x, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(exp(log(x)), x, tol);

    auto minimum_abs = mig(x);
    ARIADNE_TESTING_EXPECT_EQUAL(minimum_abs.raw(), FloatDP(5, dp));

    ARIADNE_TESTING_PROVE(possibly(sqrt(Complex<X>(4 * _one)) == 2.0_dy));
    ARIADNE_TESTING_PROVE(possibly(log(Complex<X>(1 * _one)) == 0.0_dy));
    ARIADNE_TESTING_PROVE(possibly(exp(Complex<X>(0 * _one)) == 1.0_dy));
    ARIADNE_TESTING_PROVE(possibly(log(exp(Complex<X>(1 * _one))) == 1.0_dy));
}

template<class X>
void TestComplex<X>::test_polar() {
    Dyadic eps{FloatDP::eps(dp)};
    Dyadic tol = 16 * eps;
    X zero = nul(_one);
    Complex<X> i = Complex<X>(zero, _one);
    ARIADNE_TESTING_EXPECT_EXCEPTION(atan2(zero, zero));
    ARIADNE_TESTING_EXPECT_WITHIN(arg(exp(-3 * i)), -3, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(arg(exp(-2 * i)), -2, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(arg(exp(-1 * i)), -1, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(arg(exp(0 * i)), 0, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(arg(exp(1 * i)), 1, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(arg(exp(2 * i)), 2, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(arg(exp(3 * i)), 3, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(arg(exp(-pi / 2 * i)), -pi / 2, tol);
    ARIADNE_TESTING_EXPECT_WITHIN(arg(exp(+pi / 2 * i)), +pi / 2, tol);
}

template<class X>
void TestComplex<X>::test_comparison() {
    Effort effort(0);
    Complex<X> one = _one;
    Complex<X> e = exp(one);
    //    Complex<X> elogpi=e*log(pi);
    //    ARIADNE_TEST_SAME(e,e);
    //    ARIADNE_TESTING_EXPECT_EQUAL(log(e),one);

    ARIADNE_TESTING_INSPECT("pi.get(dp)", pi.get(dp));
    ARIADNE_TESTING_INSPECT("e*log(pi).get(dp)", e * log(pi).get(dp));

    //    ARIADNE_TESTING_PROVE(not check(pi==elogpi,effort));
    //    ARIADNE_TESTING_PROVE(check(pi!=elogpi,effort));
}

template<class X>
void TestComplex<X>::test_fourier() {
    ARIADNE_TESTING_EXPECT_NO_THROW(([&] {
        auto pr = dp;

        SizeType N = 32;
        Array<X> cs(N, [&](SizeType i) {if (i==0) return X(1,pr); return X(1,pr)/(i*i); });

        using Constants::i;

        X x = X(1, pr) / 2;

        Complex<X> r(0, 0, pr);
        for (SizeType k = 0; k != 32; ++k) {
            r = r + cs[k] * exp(i * k * x);
        }
    }()));
}

int main(int argc, char** argv) {
    FloatDPApproximation::set_output_places(18);
    Ariadne::Testing::Suite suite("complex");
    auto add_complex_case = [&](std::string name, void (TestComplex<FloatDPBounds>::*method)()) {
        suite.add(std::move(name), [method](Ariadne::Testing::Context& context) {
            TestComplex<FloatDPBounds> fixture(context, FloatDP(1, dp));
            (fixture.*method)();
        });
    };
    add_complex_case("constructors", &TestComplex<FloatDPBounds>::test_constructors);
    add_complex_case("conversions", &TestComplex<FloatDPBounds>::test_conversions);
    add_complex_case("arithmetic", &TestComplex<FloatDPBounds>::test_arithmetic);
    add_complex_case("transcendental", &TestComplex<FloatDPBounds>::test_transcendental);
    add_complex_case("polar", &TestComplex<FloatDPBounds>::test_polar);
    add_complex_case("comparison", &TestComplex<FloatDPBounds>::test_comparison);
    add_complex_case("fourier", &TestComplex<FloatDPBounds>::test_fourier);
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
