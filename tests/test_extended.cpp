/***************************************************************************
 *            test_extended.cpp
 *
 *  Tests for common extended-number semantics.
 ****************************************************************************/

#include "numeric/extended.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/rational.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

void test_integer_comparison(Ariadne::Testing::Context& ctx)
{
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(1,1),Comparison::EQUAL);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(1,2),Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(2,1),Comparison::GREATER);
}

void test_rational_extended(Ariadne::Testing::Context& ctx)
{
    const Rational pinf=Rational::inf(Sign::POSITIVE);
    const Rational ninf=Rational::inf(Sign::NEGATIVE);
    const Rational nan=Rational::nan();

    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(nan+Rational(2)));
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(Rational(2)+nan));
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(nan+nan));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Rational(-2)+pinf,pinf);

    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,pinf-Rational(2),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,pinf-ninf,pinf);
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(nan-Rational(2)));

    ARIADNE_TESTING_PROVE_CTX(ctx,is_inf(sqr(pinf)));
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(sqr(nan)));
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(rec(nan)));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,hlf(pinf),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,hlf(ninf),ninf);
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(hlf(nan)));

    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(max(Rational(0),nan)));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,max(pinf,Rational(-2)),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,max(Rational(-2),ninf),Rational(-2));

    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(min(nan,Rational(0))));
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(min(Rational(0),nan)));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,min(ninf,Rational(2)),ninf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,min(Rational(2),ninf),ninf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,min(pinf,Rational(2)),Rational(2));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,min(Rational(2),pinf),Rational(2));

    ARIADNE_TESTING_PROVE_CTX(ctx,not eq(Rational(0),pinf));
    ARIADNE_TESTING_PROVE_CTX(ctx,not eq(pinf,Rational(0)));
    ARIADNE_TESTING_PROVE_CTX(ctx,not eq(nan,nan));
    ARIADNE_TESTING_PROVE_CTX(ctx,not eq(pinf,nan));

    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(nan,Rational(0)),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(Rational(0),nan),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(pinf,nan),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(pinf,ninf),Comparison::GREATER);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(ninf,pinf),Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(nan,Int(0)),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(pinf,Int(0)),Comparison::GREATER);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(ninf,Int(0)),Comparison::LESS);
}

void test_dyadic_extended(Ariadne::Testing::Context& ctx)
{
    const Dyadic pinf=Dyadic::inf(Sign::POSITIVE);
    const Dyadic ninf=Dyadic::inf(Sign::NEGATIVE);
    const Dyadic nan=Dyadic::nan();

    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(nan+Dyadic(2)));
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(Dyadic(2)+nan));
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(nan+nan));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Dyadic(-2)+pinf,pinf);

    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,pinf-Dyadic(2),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,pinf-ninf,pinf);
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(nan-Dyadic(2)));

    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,pow(pinf,Nat(2u)),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,pow(ninf,Nat(2u)),pinf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,pow(ninf,Nat(3u)),ninf);
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(pow(nan,Nat(2u))));

    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(min(nan,Dyadic(0))));
    ARIADNE_TESTING_PROVE_CTX(ctx,is_nan(min(Dyadic(0),nan)));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,min(ninf,Dyadic(2)),ninf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,min(Dyadic(2),ninf),ninf);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,min(pinf,Dyadic(2)),Dyadic(2));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,min(Dyadic(2),pinf),Dyadic(2));

    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(nan,Dyadic(0)),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(Dyadic(0),nan),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(pinf,nan),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(pinf,ninf),Comparison::GREATER);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(ninf,pinf),Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(nan,Int(0)),Comparison::INCOMPARABLE);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(pinf,Int(0)),Comparison::GREATER);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,cmp(ninf,Int(0)),Comparison::LESS);
}

} // namespace

int main(int argc, char** argv)
{
    Ariadne::Testing::Suite suite("extended");
    suite.add("integer comparison", [](Ariadne::Testing::Context& ctx) {
        test_integer_comparison(ctx);
    });
    suite.add("rational extended", [](Ariadne::Testing::Context& ctx) {
        test_rational_extended(ctx);
    });
    suite.add("dyadic extended", [](Ariadne::Testing::Context& ctx) {
        test_dyadic_extended(ctx);
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
