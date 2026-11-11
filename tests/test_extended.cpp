/***************************************************************************
 *            test_extended.cpp
 *
 *  Copyright  2026
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric. If not, see <https://www.gnu.org/licenses/>.
 */

#include "numeric/extended.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/rational.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

class TestExtended : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    void test_integer_comparison();
    void test_rational_extended();
    void test_dyadic_extended();
};

void TestExtended::test_integer_comparison()
{
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(1,1),Comparison::EQUAL);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(1,2),Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(2,1),Comparison::GREATER);
}

void TestExtended::test_rational_extended()
{
    const Rational pinf=Rational::inf(Sign::POSITIVE);
    const Rational ninf=Rational::inf(Sign::NEGATIVE);
    const Rational nan=Rational::nan();

    ARIADNE_TESTING_PROVE(is_nan(nan+Rational(2)));
    ARIADNE_TESTING_PROVE(is_nan(Rational(2)+nan));
    ARIADNE_TESTING_PROVE(is_nan(nan+nan));
    ARIADNE_TESTING_EXPECT_EQUAL(Rational(-2)+pinf,pinf);

    ARIADNE_TESTING_EXPECT_EQUAL(pinf-Rational(2),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL(pinf-ninf,pinf);
    ARIADNE_TESTING_PROVE(is_nan(nan-Rational(2)));

    ARIADNE_TESTING_PROVE(is_inf(sqr(pinf)));
    ARIADNE_TESTING_PROVE(is_nan(sqr(nan)));
    ARIADNE_TESTING_PROVE(is_nan(rec(nan)));
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(pinf),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(ninf),ninf);
    ARIADNE_TESTING_PROVE(is_nan(hlf(nan)));

    ARIADNE_TESTING_PROVE(is_nan(max(Rational(0),nan)));
    ARIADNE_TESTING_EXPECT_EQUAL(max(pinf,Rational(-2)),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL(max(Rational(-2),ninf),Rational(-2));

    ARIADNE_TESTING_PROVE(is_nan(min(nan,Rational(0))));
    ARIADNE_TESTING_PROVE(is_nan(min(Rational(0),nan)));
    ARIADNE_TESTING_EXPECT_EQUAL(min(ninf,Rational(2)),ninf);
    ARIADNE_TESTING_EXPECT_EQUAL(min(Rational(2),ninf),ninf);
    ARIADNE_TESTING_EXPECT_EQUAL(min(pinf,Rational(2)),Rational(2));
    ARIADNE_TESTING_EXPECT_EQUAL(min(Rational(2),pinf),Rational(2));

    ARIADNE_TESTING_PROVE(not eq(Rational(0),pinf));
    ARIADNE_TESTING_PROVE(not eq(pinf,Rational(0)));
    ARIADNE_TESTING_PROVE(not eq(nan,nan));
    ARIADNE_TESTING_PROVE(not eq(pinf,nan));

    ARIADNE_TESTING_EXPECT_EQUAL(cmp(nan,Rational(0)),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(Rational(0),nan),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(pinf,nan),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(pinf,ninf),Comparison::GREATER);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(ninf,pinf),Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(nan,Int(0)),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(pinf,Int(0)),Comparison::GREATER);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(ninf,Int(0)),Comparison::LESS);
}

void TestExtended::test_dyadic_extended()
{
    const Dyadic pinf=Dyadic::inf(Sign::POSITIVE);
    const Dyadic ninf=Dyadic::inf(Sign::NEGATIVE);
    const Dyadic nan=Dyadic::nan();

    ARIADNE_TESTING_PROVE(is_nan(nan+Dyadic(2)));
    ARIADNE_TESTING_PROVE(is_nan(Dyadic(2)+nan));
    ARIADNE_TESTING_PROVE(is_nan(nan+nan));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(-2)+pinf,pinf);

    ARIADNE_TESTING_EXPECT_EQUAL(pinf-Dyadic(2),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL(pinf-ninf,pinf);
    ARIADNE_TESTING_PROVE(is_nan(nan-Dyadic(2)));

    ARIADNE_TESTING_EXPECT_EQUAL(pow(pinf,Nat(2u)),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL(pow(ninf,Nat(2u)),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL(pow(ninf,Nat(3u)),ninf);
    ARIADNE_TESTING_PROVE(is_nan(pow(nan,Nat(2u))));

    ARIADNE_TESTING_PROVE(is_nan(min(nan,Dyadic(0))));
    ARIADNE_TESTING_PROVE(is_nan(min(Dyadic(0),nan)));
    ARIADNE_TESTING_EXPECT_EQUAL(min(ninf,Dyadic(2)),ninf);
    ARIADNE_TESTING_EXPECT_EQUAL(min(Dyadic(2),ninf),ninf);
    ARIADNE_TESTING_EXPECT_EQUAL(min(pinf,Dyadic(2)),Dyadic(2));
    ARIADNE_TESTING_EXPECT_EQUAL(min(Dyadic(2),pinf),Dyadic(2));

    ARIADNE_TESTING_EXPECT_EQUAL(cmp(nan,Dyadic(0)),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(Dyadic(0),nan),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(pinf,nan),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(pinf,ninf),Comparison::GREATER);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(ninf,pinf),Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(nan,Int(0)),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(pinf,Int(0)),Comparison::GREATER);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(ninf,Int(0)),Comparison::LESS);
}

} // namespace

int main(int argc, char** argv)
{
    Ariadne::Testing::Suite suite("extended");
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestExtended, test_integer_comparison);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestExtended, test_rational_extended);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestExtended, test_dyadic_extended);
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
