/***************************************************************************
 *            test_decimal.cpp
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

#include <sstream>

#include "numeric/decimal.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/integer.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"

#include "utility/test.hpp"

using namespace Ariadne;

namespace {

constexpr bool check_concept()
{
    return requires(Decimal x, Decimal y, Nat m) {
        Decimal(); Decimal(0,0u); Decimal(0.0); Decimal("0");
        x=+x; x=-x; x=x+y; x=x-y; x=x*y; x/y;
        x+=y; nul(x); sqr(x); hlf(x); abs(x); max(x,y); min(x,y);
        x.literal();
        Rational(x);
    };
}

static_assert(check_concept());

void test_decimal()
{
    Decimal d0;
    ARIADNE_TEST_EQUALS(d0,Decimal(0,0u));
    ARIADNE_TEST_CONSTRUCT(Decimal,d1,(23,1u));
    ARIADNE_TEST_CONSTRUCT(Decimal,d2,(-42,2u));
    ARIADNE_TEST_EQUALS(7_decimal,Decimal(7,0u));
    ARIADNE_TEST_EQUALS(8_dec,Decimal(8,0u));
    ARIADNE_TEST_EQUALS(Decimal(0.0),Decimal(0,0u));
    ARIADNE_TEST_EQUALS(Decimal(-3.14),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal(10000000000.0),Decimal(10000000000_z,0u));
    ARIADNE_TEST_EQUALS(Decimal("-3.14"),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal("+3.14"),Decimal(314,2u));
    ARIADNE_TEST_FAIL(Decimal("1.2.3"));
    ARIADNE_TEST_FAIL(Decimal("1a"));
    ARIADNE_TEST_FAIL(Decimal("/1"));
    ARIADNE_TEST_EQUALS(Decimal("3.141592653589793238462643383279"),Decimal(3141592653589793_z,15u)+Decimal(238462643383279_z,30u));
    ARIADNE_TEST_EQUALS(Decimal("3.141592653589793238462643383279"),Decimal("3141592653589793238462643383279"_z,30u));
    ARIADNE_TEST_FAIL(Decimal d(0.33333333333));
    ARIADNE_TEST_EQUALS(Decimal(Dyadic(7,3u)),Decimal(875,3u));
    ARIADNE_TEST_EQUALS(Decimal(Dyadic(140,1u)),Decimal(70,0u));
    ARIADNE_TEST_FAIL(static_cast<void>(Decimal(Dyadic::inf())));
    ARIADNE_TEST_EQUALS(Decimal("-3.14"),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal("-3.1400"),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal("-0.0031400"),Decimal(-314,5u));
    ARIADNE_TEST_EQUALS(Decimal("3140"),Decimal(3140,0u));
    ARIADNE_TEST_EQUALS(Decimal("3.14159265")*Decimal("3.14159265"),Decimal("9.8696043785340225"));
    ARIADNE_TEST_PRINT(Decimal("3.14159265")*Decimal("3.14159265"));
    ARIADNE_TEST_EQUALS(Rational(d1),Rational(23,10));
    ARIADNE_TEST_EQUALS(Rational(d2),Rational(-21,50));
    ARIADNE_TEST_EQUAL(Rational(+d2),+Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(-d2),-Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(d1+d2),Rational(d1)+Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(d1-d2),Rational(d1)-Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(d1*d2),Rational(d1)*Rational(d2));
    ARIADNE_TEST_EQUAL(d1/d2,Rational(d1)/Rational(d2));
    ARIADNE_TEST_EQUAL(Decimal(1,2u)/Decimal(1,1u),Rational(1,10));

    Decimal accumulated(1);
    accumulated+=Decimal("2.5");
    ARIADNE_TEST_EQUALS(accumulated,Decimal("3.5"));
    ARIADNE_TEST_EQUALS(nul(d1),Decimal(0));
    ARIADNE_TEST_EQUALS(sqr(Decimal("1.5")),Decimal("2.25"));
    ARIADNE_TEST_EQUALS(hlf(Decimal(24,1u)),Decimal(12,1u));
    ARIADNE_TEST_EQUALS(hlf(Decimal(3)),Decimal(15,1u));
    ARIADNE_TEST_EQUALS(abs(Decimal("-2.5")),Decimal("2.5"));
    ARIADNE_TEST_EQUALS(max(Decimal("-2.5"),Decimal("1.5")),Decimal("1.5"));
    ARIADNE_TEST_EQUALS(min(Decimal("-2.5"),Decimal("1.5")),Decimal("-2.5"));

    Decimal comparison_low("1.25");
    Decimal comparison_high("2.5");
    ARIADNE_TEST_ASSERT(comparison_low<=comparison_high);
    ARIADNE_TEST_ASSERT(comparison_high>=comparison_low);
    ARIADNE_TEST_ASSERT(comparison_high>comparison_low);

    PositiveDecimal positive_default;
    PositiveDecimal positive_value(Decimal("1.25"));
    PositiveDecimal positive_cast=cast_positive(Decimal("2.5"));
    ARIADNE_TEST_EQUALS(positive_default,Decimal(0));
    ARIADNE_TEST_EQUALS(positive_value,Decimal("1.25"));
    ARIADNE_TEST_EQUALS(positive_cast,Decimal("2.5"));
    ARIADNE_TEST_FAIL(PositiveDecimal(Decimal("-1")));

    DecimalBounds point_bounds(Decimal("1.25"));
    ARIADNE_TEST_EQUALS(point_bounds.lower(),Decimal("1.25"));
    ARIADNE_TEST_EQUALS(point_bounds.upper(),Decimal("1.25"));
    std::ostringstream bounds_stream;
    bounds_stream << point_bounds;
    ARIADNE_TEST_ASSERT(not bounds_stream.str().empty());

    ARIADNE_TEST_EQUALS(Decimal("3.014").literal(),String("3.014"));
    ARIADNE_TEST_EQUALS(Decimal(3).literal(),String("3."));
    ARIADNE_TEST_EQUALS(class_name<DecimalBounds>(),String("DecimalBounds"));

    ExactNumber exact_decimal=static_cast<ExactNumber>(Decimal("1.25"));
    ARIADNE_TEST_EQUALS(exact_decimal.class_name(),String("Rational"));

}

} // namespace

int main()
{
    ARIADNE_TEST_CALL(test_decimal());
    return ARIADNE_TEST_FAILURES;
}
