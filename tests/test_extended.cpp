/***************************************************************************
 *            test_extended.cpp
 *
 *  Tests for common extended-number semantics.
 ****************************************************************************/

#include "numeric/extended.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/rational.hpp"

#include "utility/test.hpp"

using namespace Ariadne;

namespace {

void test_integer_comparison()
{
    ARIADNE_TEST_EQUALS(cmp(1,1),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(1,2),Comparison::LESS);
    ARIADNE_TEST_EQUALS(cmp(2,1),Comparison::GREATER);
}

void test_rational_extended()
{
    const Rational pinf=Rational::inf(Sign::POSITIVE);
    const Rational ninf=Rational::inf(Sign::NEGATIVE);
    const Rational nan=Rational::nan();

    ARIADNE_TEST_ASSERT(is_nan(nan+Rational(2)));
    ARIADNE_TEST_ASSERT(is_nan(Rational(2)+nan));
    ARIADNE_TEST_ASSERT(is_nan(nan+nan));
    ARIADNE_TEST_EQUALS(Rational(-2)+pinf,pinf);

    ARIADNE_TEST_EQUALS(pinf-Rational(2),pinf);
    ARIADNE_TEST_EQUALS(pinf-ninf,pinf);
    ARIADNE_TEST_ASSERT(is_nan(nan-Rational(2)));

    ARIADNE_TEST_ASSERT(is_inf(sqr(pinf)));
    ARIADNE_TEST_ASSERT(is_nan(sqr(nan)));
    ARIADNE_TEST_ASSERT(is_nan(rec(nan)));

    ARIADNE_TEST_ASSERT(is_nan(max(Rational(0),nan)));
    ARIADNE_TEST_EQUALS(max(pinf,Rational(-2)),pinf);
    ARIADNE_TEST_EQUALS(max(Rational(-2),ninf),Rational(-2));

    ARIADNE_TEST_ASSERT(is_nan(min(nan,Rational(0))));
    ARIADNE_TEST_ASSERT(is_nan(min(Rational(0),nan)));
    ARIADNE_TEST_EQUALS(min(ninf,Rational(2)),ninf);
    ARIADNE_TEST_EQUALS(min(Rational(2),ninf),ninf);
    ARIADNE_TEST_EQUALS(min(pinf,Rational(2)),Rational(2));
    ARIADNE_TEST_EQUALS(min(Rational(2),pinf),Rational(2));

    ARIADNE_TEST_ASSERT(not eq(Rational(0),pinf));
    ARIADNE_TEST_ASSERT(not eq(pinf,Rational(0)));
    ARIADNE_TEST_ASSERT(not eq(nan,nan));
    ARIADNE_TEST_ASSERT(not eq(pinf,nan));

    ARIADNE_TEST_EQUALS(cmp(nan,Rational(0)),Comparison::INCOMPARABLE);
    ARIADNE_TEST_EQUALS(cmp(Rational(0),nan),Comparison::INCOMPARABLE);
    ARIADNE_TEST_EQUALS(cmp(pinf,nan),Comparison::INCOMPARABLE);
    ARIADNE_TEST_EQUALS(cmp(pinf,ninf),Comparison::GREATER);
    ARIADNE_TEST_EQUALS(cmp(ninf,pinf),Comparison::LESS);
}

void test_dyadic_extended()
{
    const Dyadic pinf=Dyadic::inf(Sign::POSITIVE);
    const Dyadic ninf=Dyadic::inf(Sign::NEGATIVE);
    const Dyadic nan=Dyadic::nan();

    ARIADNE_TEST_ASSERT(is_nan(nan+Dyadic(2)));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic(2)+nan));
    ARIADNE_TEST_ASSERT(is_nan(nan+nan));
    ARIADNE_TEST_EQUALS(Dyadic(-2)+pinf,pinf);

    ARIADNE_TEST_EQUALS(pinf-Dyadic(2),pinf);
    ARIADNE_TEST_EQUALS(pinf-ninf,pinf);
    ARIADNE_TEST_ASSERT(is_nan(nan-Dyadic(2)));

    ARIADNE_TEST_EQUALS(pow(pinf,Nat(2u)),pinf);
    ARIADNE_TEST_EQUALS(pow(ninf,Nat(2u)),pinf);
    ARIADNE_TEST_EQUALS(pow(ninf,Nat(3u)),ninf);
    ARIADNE_TEST_ASSERT(is_nan(pow(nan,Nat(2u))));

    ARIADNE_TEST_ASSERT(is_nan(min(nan,Dyadic(0))));
    ARIADNE_TEST_ASSERT(is_nan(min(Dyadic(0),nan)));
    ARIADNE_TEST_EQUALS(min(ninf,Dyadic(2)),ninf);
    ARIADNE_TEST_EQUALS(min(Dyadic(2),ninf),ninf);
    ARIADNE_TEST_EQUALS(min(pinf,Dyadic(2)),Dyadic(2));
    ARIADNE_TEST_EQUALS(min(Dyadic(2),pinf),Dyadic(2));

    ARIADNE_TEST_EQUALS(cmp(nan,Dyadic(0)),Comparison::INCOMPARABLE);
    ARIADNE_TEST_EQUALS(cmp(Dyadic(0),nan),Comparison::INCOMPARABLE);
    ARIADNE_TEST_EQUALS(cmp(pinf,nan),Comparison::INCOMPARABLE);
    ARIADNE_TEST_EQUALS(cmp(pinf,ninf),Comparison::GREATER);
    ARIADNE_TEST_EQUALS(cmp(ninf,pinf),Comparison::LESS);
}

} // namespace

int main()
{
    ARIADNE_TEST_CALL(test_integer_comparison());
    ARIADNE_TEST_CALL(test_rational_extended());
    ARIADNE_TEST_CALL(test_dyadic_extended());
    return ARIADNE_TEST_FAILURES;
}
