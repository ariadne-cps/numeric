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
#include <fstream>
#include <stdexcept>
#include <set>

#include "numeric/decimal.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/integer.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"

#include "ariadne-testing.hpp"

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

// Coefficients and scales identify each input independently of its rendering.
void test_decimal_oracle(Ariadne::Testing::Context& ctx)
{
    auto oracle_scope=ctx.scope("oracle", "decimal_exact.tsv");
    std::ifstream data(ARIADNE_TEST_DATA_DIR "/decimal_exact.tsv");
    if (!data) throw std::runtime_error("Missing decimal_exact.tsv");
    std::string line;
    std::set<unsigned> ids;
    while (std::getline(data,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        unsigned id,aq,bq,sq,pq;
        std::string ap,bp,sp,pp,extra;
        int unequal;
        if (!(fields>>id>>ap>>aq>>bp>>bq>>sp>>sq>>pp>>pq>>unequal)
            || (fields>>extra) || !ids.insert(id).second)
            throw std::runtime_error("Malformed decimal_exact.tsv row");
        const Decimal a(Integer(String(ap)),aq),b(Integer(String(bp)),bq);
        auto row_scope=ctx.scope("row", id);
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,a!=b,Boolean(unequal!=0));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,a+b,Decimal(Integer(String(sp)),sq));
        ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,a*b,Decimal(Integer(String(pp)),pq));
    }
    if (ids.size()!=64u || *ids.begin()!=0u || *ids.rbegin()!=63u)
        throw std::runtime_error("Incomplete decimal_exact.tsv");
}

void test_decimal(Ariadne::Testing::Context& ctx)
{
    Decimal d0;
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,d0,Decimal(0,0u));
    Decimal d1(23,1u);
    ctx.inspect("d1",d1);
    Decimal d2(-42,2u);
    ctx.inspect("d2",d2);
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,7_decimal,Decimal(7,0u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,8_dec,Decimal(8,0u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal(0.0),Decimal(0,0u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal(-3.14),Decimal(-314,2u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal(10000000000.0),Decimal(10000000000_z,0u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("-3.14"),Decimal(-314,2u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("+3.14"),Decimal(314,2u));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(ctx,Decimal("1.2.3"));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(ctx,Decimal("1a"));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(ctx,Decimal("/1"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("3.141592653589793238462643383279"),Decimal(3141592653589793_z,15u)+Decimal(238462643383279_z,30u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("3.141592653589793238462643383279"),Decimal("3141592653589793238462643383279"_z,30u));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(ctx,Decimal d(0.33333333333));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal(Dyadic(7,3u)),Decimal(875,3u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal(Dyadic(140,1u)),Decimal(70,0u));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(ctx,static_cast<void>(Decimal(Dyadic::inf())));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(ctx,static_cast<void>(Decimal(-Dyadic::inf())));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(ctx,static_cast<void>(Decimal(Dyadic::nan())));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("-3.14"),Decimal(-314,2u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("-3.1400"),Decimal(-314,2u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("-0.0031400"),Decimal(-314,5u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("3140"),Decimal(3140,0u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("3.14159265")*Decimal("3.14159265"),Decimal("9.8696043785340225"));
    ctx.inspect("Decimal pi approximation squared", Decimal("3.14159265")*Decimal("3.14159265"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Rational(d1),Rational(23,10));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Rational(d2),Rational(-21,50));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Rational(+d2),+Rational(d2));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Rational(-d2),-Rational(d2));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Rational(d1+d2),Rational(d1)+Rational(d2));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Rational(d1-d2),Rational(d1)-Rational(d2));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Rational(d1*d2),Rational(d1)*Rational(d2));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,d1/d2,Rational(d1)/Rational(d2));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal(1,2u)/Decimal(1,1u),Rational(1,10));

    Decimal accumulated(1);
    accumulated+=Decimal("2.5");
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,accumulated,Decimal("3.5"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,nul(d1),Decimal(0));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,sqr(Decimal("1.5")),Decimal("2.25"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,hlf(Decimal(24,1u)),Decimal(12,1u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,hlf(Decimal(3)),Decimal(15,1u));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,abs(Decimal("-2.5")),Decimal("2.5"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,max(Decimal("-2.5"),Decimal("1.5")),Decimal("1.5"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,min(Decimal("-2.5"),Decimal("1.5")),Decimal("-2.5"));

    Decimal comparison_low("1.25");
    Decimal comparison_high("2.5");
    ARIADNE_TESTING_PROVE_CTX(ctx,comparison_low<=comparison_high);
    ARIADNE_TESTING_PROVE_CTX(ctx,comparison_high>=comparison_low);
    ARIADNE_TESTING_PROVE_CTX(ctx,comparison_high>comparison_low);

    PositiveDecimal positive_default;
    PositiveDecimal positive_value(Decimal("1.25"));
    PositiveDecimal positive_cast=cast_positive(Decimal("2.5"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,positive_default,Decimal(0));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,positive_value,Decimal("1.25"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,positive_cast,Decimal("2.5"));
    ARIADNE_TESTING_EXPECT_EXCEPTION_CTX(ctx,PositiveDecimal(Decimal("-1")));

    DecimalBounds point_bounds(Decimal("1.25"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,point_bounds.lower(),Decimal("1.25"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,point_bounds.upper(),Decimal("1.25"));
    std::ostringstream bounds_stream;
    bounds_stream << point_bounds;
    ARIADNE_TESTING_PROVE_CTX(ctx,not bounds_stream.str().empty());

    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal("3.014").literal(),String("3.014"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,Decimal(3).literal(),String("3."));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,class_name<DecimalBounds>(),String("DecimalBounds"));

    ExactNumber exact_decimal=static_cast<ExactNumber>(Decimal("1.25"));
    ARIADNE_TESTING_EXPECT_EQUAL_CTX(ctx,exact_decimal.class_name(),String("Rational"));

}

} // namespace

int main()
{
    Ariadne::Testing::Suite suite("decimal");
    suite.add("decimal semantics", [](Ariadne::Testing::Context& ctx) {
        test_decimal(ctx);
    });
    suite.add("decimal exact oracle", [](Ariadne::Testing::Context& ctx) {
        test_decimal_oracle(ctx);
    });
    return Ariadne::Testing::run_suite(suite);
}
