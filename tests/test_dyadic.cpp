/***************************************************************************
 *            test_dyadic.cpp
 *
 *  Copyright  2013-20  Pieter Collins
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


#include "utility/string.hpp"
#include "utility/module.hpp"

#include "numeric/twoexp.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/builtin.hpp"
#include "numeric/integer.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "foundation/logical.hpp"

#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/float_approximation.hpp"
#include <bit>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <set>

#include "ariadne-testing.hpp"

using namespace std;
using namespace Ariadne;

namespace {

constexpr bool check_concept()
{
    return requires(unsigned int m, unsigned long int lm, int n, long int ln, Integer z, Dyadic w, Dyadic w2, Boolean b) {
        w=Dyadic(); w=Dyadic(m); w=Dyadic(lm); w=Dyadic(n); w=Dyadic(ln); w=Dyadic(z); w=Dyadic(w);
        w2=Dyadic();
        w=m; w=lm; w=n; w=ln; w=z; w=w2;
        w=+w; w=-w;
        w=w+w; w=w-w; w=w*w;
        w=1; w+=Dyadic(2);
        w=w+n; w=w-n; w=w*n;
        w=n+w; w=n-w; w=n*w;
        w=w+z; w=w-z; w=w*z;
        w=z+w; w=z-w; w=z*w;
        w=max(w,w); w=min(w,w); w=abs(w);
        w=pos(w); w=neg(w); w=sqr(w); w=hlf(w);
        w=1.5_q2; w=1.5_dy; w=1.5_dyadic;
        b=(w==w); b=(w!=w); b=(w<=w); b=(w>=w); b=(w<w); b=(w>w);
        b=(w==n); b=(w!=n); b=(w<=n); b=(w>=n); b=(w<n); b=(w>n);
        b=(n==w); b=(n!=w); b=(n<=w); b=(n>=w); b=(n<w); b=(n>w);
        b=(w==z); b=(w!=z); b=(w<=z); b=(w>=z); b=(w<z); b=(w>z);
        b=(z==w); b=(z!=w); b=(z<=w); b=(z>=w); b=(z<w); b=(z>w);
    };
}

static_assert(check_concept());

} // namespace



class TestDyadic : public Ariadne::Testing::ContextFixture
{
  public:
    using ContextFixture::ContextFixture;
    void test_literal();
    void test_conversions();
    void test_arithmetic();
    void test_bounds();
    void test_ball();
    void test_comparisons();
    void test_infinity();
    void test_native_wrappers();
    void test_dyadic_add_sub_oracles();
    void test_dyadic_conversion_oracle();
    void test_dyadic_public_contracts();
    void test_dyadic_mixed_bounds_oracles();
    void test_dyadic_bounds_power_oracle();
};

const Writer<Dyadic> fraction_write=FractionWriter();

void TestDyadic::test_literal() {
    std::ostringstream diagnostic_stream;
    Dyadic q(3.25_q2);
    ARIADNE_TESTING_INSPECT("constructed q",q);
    ARIADNE_TESTING_EXPECT_EQUAL(q,Dyadic(13,2u));
    ARIADNE_TESTING_EXPECT_EQUAL(3.25_q2,Dyadic(13,2u));
    ARIADNE_TESTING_EXPECT_EQUAL(-11.375_q2,Dyadic(-91,3u));
    ARIADNE_TESTING_EXPECT_EQUAL(0.375_q2,Dyadic(3,3u));

    ARIADNE_TESTING_EXPECT_EQUAL(_2^3,Dyadic(8));
    ARIADNE_TESTING_EXPECT_EQUAL(_2^-3,Dyadic(1,3u));
    ARIADNE_TESTING_EXPECT_EQUAL(+(_2^-4),Dyadic(1,4u));
    ARIADNE_TESTING_EXPECT_EQUAL(-(_2^-4),Dyadic(-1,4u));
    ARIADNE_TESTING_EXPECT_EQUAL(5*(_2^-3),Dyadic(5,3u));
    ARIADNE_TESTING_EXPECT_EQUAL(5/(_2^3),Dyadic(5,3u));
    ARIADNE_TESTING_EXPECT_EQUAL(5/(_2^-3),Dyadic(40));
// The following should not compile, since ^ binds less tightly than -,*,/
//    ARIADNE_TESTING_EXPECT_EQUAL(-_2^-4,Dyadic(-1,4u));
//    ARIADNE_TESTING_EXPECT_EQUAL(5*_2^-4,Dyadic(5,4u));
//    ARIADNE_TESTING_EXPECT_EQUAL(5/_2^4,Dyadic(5,4u));

    DecimalWriter decimal_write;

    ARIADNE_TESTING_INSPECT("q", q);
    ARIADNE_TESTING_EXPECT_NO_THROW(Dyadic::set_default_writer(fraction_write));
    ARIADNE_TESTING_INSPECT("q", q);
    ARIADNE_TESTING_EXPECT_NO_THROW(Dyadic::set_default_writer(fraction_write));
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<q);
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<decimal_write(q));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<fraction_write(q));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_NO_THROW(Dyadic::set_default_writer(decimal_write));
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<q);
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();

    RepresentationWriter<Dyadic> representation_write;
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<representation_write(q));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();

    ARIADNE_TESTING_EXPECT_EQUAL(to_string(fraction_write(q)),"13/2^2");
    ARIADNE_TESTING_EXPECT_EQUAL(to_string(decimal_write(q)),"3.25");
    ARIADNE_TESTING_EXPECT_EQUAL(to_string(representation_write(q)),"Dyadic(13,2u)");
    ARIADNE_TESTING_EXPECT_EQUAL(q.literal(),String("3.25"));
    ARIADNE_TESTING_EXPECT_NO_THROW(Dyadic::default_writer());
    ARIADNE_TESTING_EXPECT_EQUAL(1.5_bin,Dyadic(3,1u));

    ScientificWriter scientific_write;
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<scientific_write(Dyadic(-8)));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<scientific_write(Dyadic(1,2u)));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<scientific_write(Dyadic::inf()));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<decimal_write(Dyadic::inf()));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<decimal_write(-Dyadic::inf()));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<decimal_write(Dyadic::nan()));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_EQUAL(to_string(fraction_write(Dyadic(3))),"3");
    ARIADNE_TESTING_EXPECT_EQUAL(to_string(fraction_write(Dyadic(3,1u))),"3/2");
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<fraction_write(Dyadic::inf()));
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();

    StringStream two_exp_stream;
    two_exp_stream << TwoExp(3);
    ARIADNE_TESTING_EXPECT_EQUAL(two_exp_stream.str(),String("2^3"));

    ARIADNE_TESTING_EXPECT_EXCEPTION(Dyadic(String("not-a-dyadic")));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Dyadic(String("0.1")));
}

void TestDyadic::test_conversions() {
    Dyadic assigned;
    Dyadic assign_source(7,2u);
    assigned=assign_source;
    ARIADNE_TESTING_EXPECT_EQUAL(assigned,assign_source);

    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(Integer(-3)),Dyadic(-3,0u));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(Dyadic(-13)),Dyadic(-13));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(Dyadic(-13,3u)),Dyadic(-13,3u));
    Natural huge_exponent{Integer{String("9223372036854775808")}};
    ARIADNE_TESTING_EXPECT_EXCEPTION((Dyadic{Integer(1),huge_exponent}));
    Natural valid_exponent{63u};
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(Integer(1),valid_exponent),Dyadic(1,63u));
    ARIADNE_TESTING_EXPECT_EXCEPTION((PositiveDyadic{Dyadic(-1)}));

    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(-11,2u)),Integer(-3));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(-10,2u)),Integer(-3));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(-9,2u)),Integer(-2));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(-2,2u)),Integer(-1));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(-1,2u)),Integer(0));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(0,2u)),Integer(0));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(1,2u)),Integer(0));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(2,2u)),Integer(1));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(9,2u)),Integer(2));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(10,2u)),Integer(3));
    ARIADNE_TESTING_EXPECT_EQUAL(round(Dyadic(11,2u)),Integer(3));
    ARIADNE_TESTING_EXPECT_EQUAL(floor(Dyadic(-5,1u)),Integer(-3));
    ARIADNE_TESTING_EXPECT_EQUAL(ceil(Dyadic(-5,1u)),Integer(-2));
    ARIADNE_TESTING_EXPECT_EXCEPTION(round(Dyadic::inf()));
    ARIADNE_TESTING_EXPECT_EXCEPTION(floor(Dyadic::inf()));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ceil(Dyadic::inf()));

    ARIADNE_TESTING_EXPECT_EQUAL(static_cast<Int>(3)/two,Dyadic(3,1u));
    ARIADNE_TESTING_EXPECT_EQUAL(Integer(3)/two,Dyadic(3,1u));
    ARIADNE_TESTING_EXPECT_EQUAL(ExactDouble(3)/two,Dyadic(3,1u));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(3)/two,Dyadic(3,1u));

    ARIADNE_TESTING_EXPECT_EQUAL(make_dyadic(7u),Dyadic(7));
    ARIADNE_TESTING_EXPECT_EQUAL(make_dyadic(1ULL<<40),Dyadic(TwoExp(40)));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<Dyadic>(),String("Dyadic"));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<DyadicBounds>(),String("DyadicBounds"));
}

void TestDyadic::test_arithmetic() {
    Dyadic accumulated=1;
    accumulated+=Dyadic(2);
    ARIADNE_TESTING_EXPECT_EQUAL(accumulated,Dyadic(3));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(3,2u)+Dyadic(-5,3u),Dyadic(1,3u));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(3,2u)-Dyadic(-5,3u),Dyadic(11,3u));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(3,2u)*Dyadic(-5,3u),Dyadic(-15,5u));
    ARIADNE_TESTING_EXPECT_EQUAL(pow(Dyadic(2),static_cast<Int>(3)),Rational(8));
    ARIADNE_TESTING_EXPECT_EQUAL(pow(Dyadic(2),static_cast<Int>(-3)),Rational(1,8));
}

void TestDyadic::test_bounds() {
    std::ostringstream diagnostic_stream;
    DyadicBounds positive(1,3);
    DyadicBounds negative(-3,-1);
    DyadicBounds crossing(-2,3);
    DyadicBounds singleton(2,2);

    ARIADNE_TESTING_EXPECT_EQUAL((+positive).lower_raw(),Dyadic(1));
    ARIADNE_TESTING_EXPECT_EQUAL((-positive).lower_raw(),Dyadic(-3));
    ARIADNE_TESTING_EXPECT_EQUAL((positive+DyadicBounds(2,4)).upper_raw(),Dyadic(7));
    ARIADNE_TESTING_EXPECT_EQUAL((positive-DyadicBounds(2,4)).lower_raw(),Dyadic(-3));
    ARIADNE_TESTING_EXPECT_EQUAL(add(positive,DyadicBounds(2,4)).lower_raw(),Dyadic(3));
    ARIADNE_TESTING_EXPECT_EQUAL(sub(positive,DyadicBounds(2,4)).upper_raw(),Dyadic(1));
    ARIADNE_TESTING_EXPECT_EQUAL(nul(positive).lower_raw(),Dyadic(0));
    ARIADNE_TESTING_EXPECT_EQUAL(pos(positive).upper_raw(),Dyadic(3));
    ARIADNE_TESTING_EXPECT_EQUAL(neg(positive).upper_raw(),Dyadic(-1));
    ARIADNE_TESTING_EXPECT_EQUAL(sqr(positive).lower_raw(),Dyadic(1));
    ARIADNE_TESTING_EXPECT_EQUAL(sqr(negative).upper_raw(),Dyadic(9));
    ARIADNE_TESTING_EXPECT_EQUAL(sqr(crossing).lower_raw(),Dyadic(0));
    ARIADNE_TESTING_EXPECT_EQUAL(sqr(crossing).upper_raw(),Dyadic(9));
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(positive).upper_raw(),Dyadic(3,1u));
    ARIADNE_TESTING_EXPECT_EQUAL(abs(crossing).lower_raw(),Dyadic(0));
    ARIADNE_TESTING_EXPECT_EQUAL(max(positive,DyadicBounds(2,4)).lower_raw(),Dyadic(2));
    ARIADNE_TESTING_EXPECT_EQUAL(min(positive,DyadicBounds(2,4)).upper_raw(),Dyadic(3));

    ARIADNE_TESTING_PROVE(definitely(sgn(positive)));
    ARIADNE_TESTING_PROVE(definitely(!sgn(negative)));
    ARIADNE_TESTING_PROVE(is_indeterminate(sgn(crossing)));

    DyadicBounds left(1,2);
    DyadicBounds right(3,4);
    DyadicBounds overlap(2,3);
    ARIADNE_TESTING_PROVE(definitely(singleton==singleton));
    ARIADNE_TESTING_PROVE(definitely(!(left==right)));
    ARIADNE_TESTING_PROVE(is_indeterminate(left==overlap));
    ARIADNE_TESTING_PROVE(definitely(left!=right));
    ARIADNE_TESTING_PROVE(definitely(right!=left));
    ARIADNE_TESTING_PROVE(definitely(!(singleton!=singleton)));
    ARIADNE_TESTING_PROVE(is_indeterminate(left!=overlap));

    ARIADNE_TESTING_PROVE(consistent(left,overlap));
    ARIADNE_TESTING_PROVE(not consistent(left,right));
    ARIADNE_TESTING_PROVE(inconsistent(left,right));
    ARIADNE_TESTING_PROVE(inconsistent(right,left));
    ARIADNE_TESTING_PROVE(not inconsistent(left,overlap));
    ARIADNE_TESTING_PROVE(refines(DyadicBounds(2,2),DyadicBounds(1,3)));
    ARIADNE_TESTING_PROVE(not refines(DyadicBounds(0,2),DyadicBounds(1,3)));
    ARIADNE_TESTING_PROVE(not refines(DyadicBounds(2,4),DyadicBounds(1,3)));
    ARIADNE_TESTING_EXPECT_EQUAL(refinement(DyadicBounds(1,3),DyadicBounds(2,4)).lower_raw(),Dyadic(2));
    ARIADNE_TESTING_EXPECT_EQUAL(coarsening(DyadicBounds(1,3),DyadicBounds(2,4)).upper_raw(),Dyadic(4));

    ARIADNE_TESTING_EXPECT_EQUAL(pow(crossing,Nat(2u)).upper_raw(),Dyadic(9));
    ARIADNE_TESTING_EXPECT_EQUAL(pow(crossing,Nat(3u)).lower_raw(),Dyadic(-8));

    DyadicLowerBound lower(Dyadic(1));
    DyadicLowerBound lower_from_bounds(positive);
    DyadicUpperBound upper(Dyadic(3));
    DyadicUpperBound upper_from_bounds(positive);
    DyadicApproximation approximation(Dyadic(2));
    ARIADNE_TESTING_EXPECT_EQUAL(lower.raw(),Dyadic(1));
    ARIADNE_TESTING_EXPECT_EQUAL(lower_from_bounds.raw(),Dyadic(1));
    ARIADNE_TESTING_EXPECT_EQUAL(upper.raw(),Dyadic(3));
    ARIADNE_TESTING_EXPECT_EQUAL(upper_from_bounds.raw(),Dyadic(3));
    ARIADNE_TESTING_EXPECT_EQUAL(approximation.raw(),Dyadic(2));
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<lower<<" "<<upper);
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
    ARIADNE_TESTING_EXPECT_NO_THROW((PositiveDyadicBounds{positive}));
    ARIADNE_TESTING_EXPECT_NO_THROW((PositiveDyadicLowerBound{lower}));
    ARIADNE_TESTING_EXPECT_NO_THROW((PositiveDyadicUpperBound{upper}));
    ARIADNE_TESTING_EXPECT_NO_THROW((PositiveDyadicApproximation{approximation}));
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<positive);
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
}

void TestDyadic::test_ball() {
    std::ostringstream diagnostic_stream;
    DyadicBall exact(Dyadic(2));
    DyadicBall uncertain(Dyadic(3),Dyadic(1));
    DyadicBall from_bounds(DyadicBounds(1,5));
    ARIADNE_TESTING_EXPECT_EQUAL(exact.value(),Dyadic(2));
    ARIADNE_TESTING_EXPECT_EQUAL(exact.error(),Dyadic(0));
    ARIADNE_TESTING_EXPECT_EQUAL(uncertain.value_raw(),Dyadic(3));
    ARIADNE_TESTING_EXPECT_EQUAL(uncertain.error_raw(),Dyadic(1));
    ARIADNE_TESTING_EXPECT_EQUAL(from_bounds.value(),Dyadic(3));
    ARIADNE_TESTING_EXPECT_EQUAL(from_bounds.error(),Dyadic(2));

    DyadicBounds converted=static_cast<DyadicBounds>(uncertain);
    ARIADNE_TESTING_EXPECT_EQUAL(converted.lower_raw(),Dyadic(2));
    ARIADNE_TESTING_EXPECT_EQUAL(converted.upper_raw(),Dyadic(4));

    ARIADNE_TESTING_EXPECT_EQUAL((exact+uncertain).value_raw(),Dyadic(5));
    ARIADNE_TESTING_EXPECT_EQUAL((uncertain-exact).value_raw(),Dyadic(1));
    ARIADNE_TESTING_EXPECT_EQUAL((exact*uncertain).value_raw(),Dyadic(6));
    ARIADNE_TESTING_EXPECT_EQUAL(abs(DyadicBall(Dyadic(-3),Dyadic(1))).value_raw(),Dyadic(3));
    DyadicBall crossing_abs=abs(DyadicBall(Dyadic(0),Dyadic(2)));
    ARIADNE_TESTING_EXPECT_EQUAL(crossing_abs.value_raw(),Dyadic(1));
    ARIADNE_TESTING_EXPECT_EQUAL(crossing_abs.error_raw(),Dyadic(1));

    ARIADNE_TESTING_PROVE(definitely(DyadicBall(Dyadic(1),Dyadic(0)) < DyadicBall(Dyadic(3),Dyadic(0))));
    ARIADNE_TESTING_PROVE(definitely(!(DyadicBall(Dyadic(3),Dyadic(0)) < DyadicBall(Dyadic(1),Dyadic(0)))));
    ARIADNE_TESTING_PROVE(is_indeterminate(DyadicBall(Dyadic(2),Dyadic(1)) < DyadicBall(Dyadic(3),Dyadic(1))));
    ARIADNE_TESTING_PROVE(inconsistent(DyadicBall(Dyadic(1),Dyadic(0)),DyadicBall(Dyadic(3),Dyadic(0))));
    ARIADNE_TESTING_PROVE(not inconsistent(DyadicBall(Dyadic(2),Dyadic(1)),DyadicBall(Dyadic(3),Dyadic(1))));
    ARIADNE_TESTING_PROVE(refines(DyadicBall(Dyadic(2),Dyadic(1)),DyadicBall(Dyadic(2),Dyadic(2))));
    ARIADNE_TESTING_PROVE(not refines(DyadicBall(Dyadic(2),Dyadic(2)),DyadicBall(Dyadic(2),Dyadic(1))));
    ARIADNE_TESTING_EXPECT_NO_THROW(refinement(DyadicBall(Dyadic(2),Dyadic(1)),DyadicBall(Dyadic(2),Dyadic(2))));
    ARIADNE_TESTING_EXPECT_NO_THROW(coarsening(DyadicBall(Dyadic(2),Dyadic(1)),DyadicBall(Dyadic(3),Dyadic(1))));
    ARIADNE_TESTING_EXPECT_NO_THROW(diagnostic_stream<<uncertain);
    ARIADNE_TESTING_INSPECT("stream output", diagnostic_stream.str());
    diagnostic_stream.str("");
    diagnostic_stream.clear();
}

void TestDyadic::test_comparisons() {
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<,Dyadic(-4,3u),Dyadic(-1,2u));
}

// D48: native exact bounds and ball arithmetic, without GMP storage.
void TestDyadic::test_native_wrappers() {
    const DyadicBounds positive(Dyadic(1),Dyadic(3));
    const DyadicBounds negative(Dyadic(-3),Dyadic(-1));
    const DyadicBounds crossing(Dyadic(-2),Dyadic(3));
    ARIADNE_TESTING_PROVE((positive+DyadicBounds(Dyadic(2),Dyadic(4))).upper_raw()
          ==Dyadic(7));
    ARIADNE_TESTING_PROVE((positive-negative).lower_raw()==Dyadic(2));
    ARIADNE_TESTING_PROVE((positive*negative).lower_raw()==Dyadic(-9));
    ARIADNE_TESTING_PROVE((positive*negative).upper_raw()==Dyadic(-1));
    ARIADNE_TESTING_PROVE(sqr(crossing).lower_raw()==Dyadic(0));
    ARIADNE_TESTING_PROVE(sqr(crossing).upper_raw()==Dyadic(9));
    ARIADNE_TESTING_PROVE(pow(crossing,Nat(3u)).lower_raw()==Dyadic(-8));
    ARIADNE_TESTING_PROVE(pow(crossing,Nat(2u)).lower_raw()==Dyadic(0));
    ARIADNE_TESTING_PROVE(abs(crossing).upper_raw()==Dyadic(3));
    ARIADNE_TESTING_PROVE(refines(DyadicBounds(Dyadic(2),Dyadic(2)),positive));
    const DyadicBall ball(Dyadic(3),Dyadic(1));
    const DyadicBounds enclosure=static_cast<DyadicBounds>(ball);
    ARIADNE_TESTING_PROVE(enclosure.lower_raw()==Dyadic(2) &&
          enclosure.upper_raw()==Dyadic(4));
    ARIADNE_TESTING_PROVE((ball*DyadicBall(Dyadic(2),Dyadic(0))).value_raw()
          ==Dyadic(6));
    const DyadicLowerBound lower(positive);
    const DyadicUpperBound upper(positive);
    const DyadicApproximation approximation(Dyadic(2));
    ARIADNE_TESTING_PROVE(lower.raw()==Dyadic(1) && upper.raw()==Dyadic(3)
          && approximation.raw()==Dyadic(2));
}


void TestDyadic::test_infinity() {

    // Verify the new canonical zero representation without GMP internals.
    const Dyadic z(0);
    ARIADNE_TESTING_PROVE(z.coefficient()==Integer(0));
    ARIADNE_TESTING_PROVE(z.denominator_exponent()==0u);
    ARIADNE_TESTING_PROVE(Dyadic(0,17u).denominator_exponent()==0u);

    ARIADNE_TESTING_INSPECT("Dyadic::nan()", Dyadic::nan());
    ARIADNE_TESTING_INSPECT("Dyadic::inf()", Dyadic::inf());
    ARIADNE_TESTING_INSPECT("Dyadic::inf(Sign(0))", Dyadic::inf(Sign(0)));
    ARIADNE_TESTING_INSPECT("Dyadic::inf(Sign(+1))", Dyadic::inf(Sign(+1)));
    ARIADNE_TESTING_INSPECT("Dyadic::inf(Sign(-1))", Dyadic::inf(Sign(-1)));ARIADNE_TESTING_PROVE(is_nan(Dyadic::nan()));
    ARIADNE_TESTING_PROVE(is_inf(Dyadic::inf()));
    ARIADNE_TESTING_PROVE(is_inf(Dyadic::inf(Sign(+1))));
    ARIADNE_TESTING_PROVE(is_inf(Dyadic::inf(Sign(-1))));
    ARIADNE_TESTING_PROVE(not is_inf(Dyadic::nan()));
    ARIADNE_TESTING_PROVE(not is_inf(Dyadic(0)));
    ARIADNE_TESTING_PROVE(is_finite(Dyadic(0)));
    ARIADNE_TESTING_PROVE(is_zero(Dyadic(0)));
    ARIADNE_TESTING_PROVE(not is_zero(Dyadic(1)));
    ARIADNE_TESTING_PROVE(not is_zero(Dyadic::inf()));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf(Sign(+1)),Dyadic::inf());
    ARIADNE_TESTING_PROVE(Dyadic::inf(Sign(+1))>Dyadic(0));
    ARIADNE_TESTING_PROVE(Dyadic::inf(Sign(-1))<Dyadic(0));

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator==,Dyadic::inf(),Dyadic::inf());
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<,Dyadic(-1,2u),Dyadic::inf());
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<,Dyadic(0),Dyadic::inf());
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<,Dyadic(3,1u),Dyadic::inf());
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<,-Dyadic::inf(),Dyadic::inf());
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<,-Dyadic::inf(),Dyadic(1,2u));
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<,-Dyadic::inf(),Dyadic(-3,1u));

    ExactDouble double_max(std::numeric_limits<double>::max());
    ExactDouble double_inf(std::numeric_limits<double>::infinity());
    ExactDouble double_nan(double_inf.get_d()*0.0);

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<,Dyadic(double_max),Dyadic(double_inf));

    ARIADNE_TESTING_PROVE(Dyadic::inf(Sign::POSITIVE).get_d()==std::numeric_limits<double>::infinity());
    ARIADNE_TESTING_PROVE(Dyadic::inf(Sign::NEGATIVE).get_d()==-std::numeric_limits<double>::infinity());
    ARIADNE_TESTING_PROVE(isnan(Dyadic::inf(Sign::ZERO).get_d()));

    ARIADNE_TESTING_PROVE(Dyadic(+double_inf)==Dyadic::inf(Sign(+1)));
    ARIADNE_TESTING_PROVE(Dyadic(-double_inf)==Dyadic::inf(Sign(-1)));
    ARIADNE_TESTING_PROVE(is_nan(Dyadic(double_nan)));

    ARIADNE_TESTING_PROVE(is_nan(-Dyadic::nan()));
    ARIADNE_TESTING_EXPECT_EQUAL(-Dyadic::inf(Sign::POSITIVE),Dyadic::inf(Sign::NEGATIVE));
    ARIADNE_TESTING_EXPECT_EQUAL(-Dyadic::inf(Sign::NEGATIVE),Dyadic::inf(Sign::POSITIVE));
    ARIADNE_TESTING_EXPECT_EQUAL(sgn(Dyadic::inf()), Sign::POSITIVE);
    ARIADNE_TESTING_EXPECT_EQUAL(sgn(-Dyadic::inf()), Sign::NEGATIVE);

    ARIADNE_TESTING_PROVE(is_nan(Dyadic::inf()+(-Dyadic::inf())));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf()+Dyadic::inf(),Dyadic::inf());
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf()+Dyadic(-2),Dyadic::inf());

    ARIADNE_TESTING_PROVE(is_nan(Dyadic::inf()-Dyadic::inf()));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(2)-Dyadic::inf(),-Dyadic::inf());

    ARIADNE_TESTING_PROVE(is_nan(Dyadic::inf(Sign(+1))*Dyadic(0)));
    ARIADNE_TESTING_PROVE(is_nan(Dyadic::inf(Sign(-1))*Dyadic(0)));
    ARIADNE_TESTING_PROVE(is_nan(Dyadic(0)*Dyadic::inf(Sign(+1))));
    ARIADNE_TESTING_PROVE(is_nan(Dyadic(0)*Dyadic::inf(Sign(-1))));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf(Sign(+1))*Dyadic::inf(Sign(+1)),Dyadic::inf(Sign(+1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf(Sign(+1))*Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(-1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf(Sign(-1))*Dyadic::inf(Sign(+1)),Dyadic::inf(Sign(-1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf(Sign(-1))*Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(+1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(+2)*Dyadic::inf(Sign(+1)),Dyadic::inf(Sign(+1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(+2)*Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(-1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(-2)*Dyadic::inf(Sign(+1)),Dyadic::inf(Sign(-1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(-2)*Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(+1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf(Sign(+1))*Dyadic(+2),Dyadic::inf(Sign(+1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf(Sign(+1))*Dyadic(-2),Dyadic::inf(Sign(-1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf(Sign(-1))*Dyadic(+2),Dyadic::inf(Sign(-1)));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf(Sign(-1))*Dyadic(-2),Dyadic::inf(Sign(+1)));

    ARIADNE_TESTING_PROVE(is_nan(hlf(Dyadic::nan())));
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(Dyadic::inf()),Dyadic::inf());
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(-Dyadic::inf()),-Dyadic::inf());

    ARIADNE_TESTING_PROVE(is_nan(abs(Dyadic::nan())));
    ARIADNE_TESTING_EXPECT_EQUAL(abs(Dyadic::inf()),Dyadic::inf());
    ARIADNE_TESTING_EXPECT_EQUAL(abs(-Dyadic::inf()),Dyadic::inf());

    ARIADNE_TESTING_PROVE(is_nan(max(Dyadic::nan(),Dyadic(0))));
    ARIADNE_TESTING_PROVE(is_nan(max(Dyadic::nan(),Dyadic::inf())));
    ARIADNE_TESTING_EXPECT_EQUAL(max(Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(-1))),Dyadic::inf(Sign(-1)));
    ARIADNE_TESTING_EXPECT_EQUAL(max(Dyadic::inf(Sign(-1)),Dyadic(-2)),Dyadic(-2));
    ARIADNE_TESTING_EXPECT_EQUAL(max(Dyadic(-2),Dyadic::inf(Sign(+1))),Dyadic::inf(Sign(+1)));

    ARIADNE_TESTING_PROVE(DyadicBounds(1,3)<=3);
    ARIADNE_TESTING_PROVE(DyadicBounds(1,3)>=1);
    ARIADNE_TESTING_PROVE(DyadicBounds(1,3)< 4);
    ARIADNE_TESTING_PROVE(DyadicBounds(1,3)> 0);
    ARIADNE_TESTING_PROVE(is_indeterminate(DyadicBounds(1,3)<=2));
    ARIADNE_TESTING_PROVE(is_indeterminate(DyadicBounds(1,3)>=2));
    ARIADNE_TESTING_PROVE(is_indeterminate(DyadicBounds(1,3)< 3));
    ARIADNE_TESTING_PROVE(is_indeterminate(DyadicBounds(1,3)> 1));
    ARIADNE_TESTING_PROVE(!(DyadicBounds(1,3)<=0));
    ARIADNE_TESTING_PROVE(!(DyadicBounds(1,3)>=4));
    ARIADNE_TESTING_PROVE(not(DyadicBounds(1,3)< 1));
    ARIADNE_TESTING_PROVE(not(DyadicBounds(1,3)> 3));
}


void TestDyadic::test_dyadic_add_sub_oracles() {
    for (unsigned limbs:{1u,5u,16u,17u,31u,32u,64u,65u}) {
        const std::string name="dyadic_add_sub_"+std::to_string(limbs)+".tsv";
        std::ifstream data(std::string(ARIADNE_TEST_DATA_DIR)+"/"+name);
        if (!data) throw std::runtime_error("Missing "+name);
        std::string line;
        std::set<unsigned> ids;
        while (std::getline(data,line)) {
            if (line.empty() || line[0]=='#') continue;
            std::istringstream fields(line);
            unsigned id;
            std::size_t ae,be,se,de;
            std::string an,bn,sn,dn,extra;
            if (!(fields>>id>>an>>ae>>bn>>be>>sn>>se>>dn>>de) || (fields>>extra)
                || !ids.insert(id).second)
                throw std::runtime_error("Malformed "+name);
            const Dyadic a(Integer(String(an)),ae),b(Integer(String(bn)),be);
            const Dyadic sum_expected(Integer(String(sn)),se),difference_expected(Integer(String(dn)),de);
            [[maybe_unused]] auto row_scope = this->ctx().scope(name + " row", id);
            ARIADNE_TESTING_EXPECT_EQUAL(a+b,sum_expected);
            ARIADNE_TESTING_EXPECT_EQUAL(b+a,sum_expected);
            ARIADNE_TESTING_EXPECT_EQUAL(a-b,difference_expected);
            ARIADNE_TESTING_EXPECT_EQUAL(b-a,-difference_expected);
            // Canonicalization is observable through the public accessors.
            const auto sum=a+b,difference=a-b;
            ARIADNE_TESTING_EXPECT_EQUAL(sum.coefficient().literal(),String(sn));
            ARIADNE_TESTING_EXPECT_EQUAL(sum.denominator_exponent(),se);
            ARIADNE_TESTING_EXPECT_EQUAL(difference.coefficient().literal(),String(dn));
            ARIADNE_TESTING_EXPECT_EQUAL(difference.denominator_exponent(),de);
        }
        if (ids.size()!=256u || *ids.begin()!=0u || *ids.rbegin()!=255u)
            throw std::runtime_error("Incomplete "+name);
    }
}

void TestDyadic::test_dyadic_conversion_oracle() {
    std::ifstream data(ARIADNE_TEST_DATA_DIR "/dyadic_conversion.tsv");
    if (!data) throw std::runtime_error("Missing dyadic_conversion.tsv");
    std::string line;
    std::set<unsigned> ids;
    const MultiplePrecision precision(2048_bits);
    while (std::getline(data,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        unsigned id,exponent;
        std::uint64_t expected;
        int exact;
        std::string coefficient,denominator,floor_expected,ceil_expected,extra;
        if (!(fields>>id>>coefficient>>exponent>>expected>>exact>>denominator>>floor_expected>>ceil_expected) || (fields>>extra)
            || !ids.insert(id).second)
            throw std::runtime_error("Malformed dyadic_conversion.tsv");
        const Dyadic source(Integer(String(coefficient)),static_cast<std::size_t>(exponent));
        [[maybe_unused]] auto row_scope = this->ctx().scope("dyadic conversion row", id);
        ARIADNE_TESTING_EXPECT_EQUAL(std::bit_cast<std::uint64_t>(source.get_d()),expected);
        if (exponent==0u) ARIADNE_TESTING_EXPECT_EQUAL(source.literal(),String(coefficient));
        ARIADNE_TESTING_EXPECT_EQUAL(floor(source),Integer(String(floor_expected)));
        ARIADNE_TESTING_EXPECT_EQUAL(ceil(source),Integer(String(ceil_expected)));
        ARIADNE_TESTING_EXPECT_EQUAL(hlf(source),Dyadic(Integer(String(coefficient)),static_cast<std::size_t>(exponent)+1u));
        ARIADNE_TESTING_EXPECT_EQUAL(fma(source,Dyadic(1),Dyadic(0)),source);
        const Rational exact_value{Integer(String(coefficient)),Integer(String(denominator))};
        const Rational reciprocal{Integer(String(denominator)),Integer(String(coefficient))};
        ARIADNE_TESTING_EXPECT_EQUAL(div(source,Integer(1)),exact_value);
        ARIADNE_TESTING_EXPECT_EQUAL(div(Integer(1),source),reciprocal);
        const FloatMP mp(source,FloatMP::ROUND_TO_NEAREST,precision);
        ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(mp),source);
        const Approximation<Dyadic> approximate(source);
        ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(approximate.get(precision).raw()),source);
        if (exact!=0) {
            const double value=std::bit_cast<double>(expected);
            const FloatDP fp(ExactDouble(value),dp);
            ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(fp),source);
            ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(ExactDouble(value)),source);
            ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(value),source);
            ARIADNE_TESTING_EXPECT_EQUAL(std::bit_cast<std::uint64_t>(approximate.get(dp).raw().get_d()),expected);
        }
    }
    if (ids.size()!=168u || *ids.begin()!=0u || *ids.rbegin()!=167u)
        throw std::runtime_error("Incomplete dyadic_conversion.tsv");
}

void TestDyadic::test_dyadic_public_contracts() {
    const std::size_t largest=std::numeric_limits<std::size_t>::max();
    const Dyadic tiny(Integer(1),largest);
    ARIADNE_TESTING_EXPECT_EXCEPTION(tiny.exponent());
    ARIADNE_TESTING_EXPECT_EXCEPTION(tiny.literal());
    ARIADNE_TESTING_EXPECT_EXCEPTION(tiny*Dyadic(1,1u));
    ARIADNE_TESTING_EXPECT_EXCEPTION(hlf(tiny));
    const Integer outside_size=Integer::from_magnitude(shift_left(BigUnsigned(1u),
        std::numeric_limits<std::size_t>::digits),false);
    const Natural too_large(outside_size);
    ARIADNE_TESTING_EXPECT_EXCEPTION((Dyadic(Integer(1),too_large)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Dyadic(1)/3);

    const PositiveDyadic default_positive;
    const Dyadic positive(3,1u),negative(-3,1u);
    const PositiveDyadic copied_positive(positive);
    ARIADNE_TESTING_EXPECT_EQUAL(default_positive,Dyadic(0));
    ARIADNE_TESTING_EXPECT_EQUAL(copied_positive,positive);
    ARIADNE_TESTING_EXPECT_EQUAL(cast_positive(positive),positive);
    ARIADNE_TESTING_EXPECT_EXCEPTION(cast_positive(negative));
    ARIADNE_TESTING_PROVE(positive.kind()==Dyadic::Kind::Finite);
    ARIADNE_TESTING_PROVE(Dyadic::inf().kind()==Dyadic::Kind::PositiveInfinity);
    ARIADNE_TESTING_PROVE((-Dyadic::inf()).kind()==Dyadic::Kind::NegativeInfinity);
    ARIADNE_TESTING_PROVE(Dyadic::nan().kind()==Dyadic::Kind::NaN);
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::nan().literal(),String("NaN"));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic::inf().literal(),String("inf"));
    ARIADNE_TESTING_EXPECT_EQUAL((-Dyadic::inf()).literal(),String("-inf"));
    ARIADNE_TESTING_PROVE(is_nan(Dyadic(Rational::nan())));
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(Rational::inf()),Dyadic::inf());
    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(-Rational::inf()),-Dyadic::inf());
    ARIADNE_TESTING_PROVE(is_nan(Dyadic::inf(Sign::ZERO)));
    ARIADNE_TESTING_PROVE(is_nan(-Dyadic::nan()));
    ARIADNE_TESTING_PROVE(is_nan(Dyadic(0)*Dyadic::inf()));
    ARIADNE_TESTING_PROVE(is_nan(Dyadic::inf()*Dyadic(0)));
    ARIADNE_TESTING_PROVE(is_nan(Dyadic::nan()*positive));
    ARIADNE_TESTING_PROVE(is_nan(positive*Dyadic::nan()));
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(Int(1),positive),Comparison::LESS);
    ARIADNE_TESTING_PROVE(positive!=negative);
    ARIADNE_TESTING_PROVE(!eq(positive,Dyadic::nan()));
    ARIADNE_TESTING_PROVE(!eq(positive,Dyadic(3)));
    ARIADNE_TESTING_PROVE(!lt(positive,Dyadic::nan()));
    ARIADNE_TESTING_PROVE(!lt(Dyadic::nan(),positive));
    ARIADNE_TESTING_PROVE(is_nan(max(positive,Dyadic::nan())));

    ARIADNE_TESTING_EXPECT_EQUAL(Dyadic(String("1.5")),positive);
    ARIADNE_TESTING_EXPECT_EQUAL(dy_(1.5L),positive);
    std::ostringstream stream;
    stream << Approximation<Dyadic>(positive);
    ARIADNE_TESTING_EXPECT_EQUAL(stream.str(),std::string("1.5"));
}

void TestDyadic::test_dyadic_mixed_bounds_oracles() {
    for (bool logical:{false,true}) {
        const std::string name=logical?"dyadic_mixed_comparisons.tsv":"dyadic_mixed_bounds.tsv";
        std::ifstream data(std::string(ARIADNE_TEST_DATA_DIR)+"/"+name);
        if (!data) throw std::runtime_error("Missing "+name);
        std::string line;
        std::set<unsigned> ids;
        while (std::getline(data,line)) {
            if (line.empty() || line[0]=='#') continue;
            std::istringstream fields(line);
            unsigned id,reverse;
            Int dl,du,dd,ql,qu,qd;
            std::string op,extra;
            if (!(fields>>id>>op>>reverse>>dl>>du>>dd>>ql>>qu>>qd)
                || !ids.insert(id).second || dd!=2 || qd!=3)
                throw std::runtime_error("Malformed "+name);
            const DyadicBounds d(Dyadic(dl,1u),Dyadic(du,1u));
            const RationalBounds q(Rational(ql,qd),Rational(qu,qd));
            [[maybe_unused]] auto row_scope = this->ctx().scope(name + " row", id);
            if (!logical) {
                Int ln,ld,un,ud;
                if (!(fields>>ln>>ld>>un>>ud) || (fields>>extra))
                    throw std::runtime_error("Malformed interval result");
                const auto evaluate=[&]<class A,class B>(A const& a,B const& b)->RationalBounds {
                    if (op=="add") return add(a,b);
                    if (op=="sub") return sub(a,b);
                    if (op=="mul") return mul(a,b);
                    if (op=="div") return div(a,b);
                    if (op=="max") return max(a,b);
                    if (op=="min") return min(a,b);
                    throw std::runtime_error("Unknown mixed interval operation");
                };
                const auto result=reverse?evaluate(q,d):evaluate(d,q);
                ARIADNE_TESTING_EXPECT_EQUAL(result.lower(),Rational(ln,ld));
                ARIADNE_TESTING_EXPECT_EQUAL(result.upper(),Rational(un,ud));
                if (op=="div" && ql%qd==0 && qu%qd==0) {
                    const DyadicBounds other(Dyadic(ql/qd),Dyadic(qu/qd));
                    const auto direct=reverse?div(other,d):div(d,other);
                    ARIADNE_TESTING_EXPECT_EQUAL(direct.lower(),Rational(ln,ld));
                    ARIADNE_TESTING_EXPECT_EQUAL(direct.upper(),Rational(un,ud));
                }
            } else {
                int expected;
                if (!(fields>>expected) || (fields>>extra))
                    throw std::runtime_error("Malformed logical result");
                const auto evaluate=[&]<class A,class B>(A const& a,B const& b)->ValidatedKleenean {
                    if (op=="eq") return a==b;
                    if (op=="ne") return a!=b;
                    if (op=="lt") return a<b;
                    if (op=="le") return a<=b;
                    if (op=="gt") return a>b;
                    if (op=="ge") return a>=b;
                    throw std::runtime_error("Unknown mixed interval comparison");
                };
                const auto result=reverse?evaluate(q,d):evaluate(d,q);
                const int truth=definitely(result)?1:definitely(!result)?-1:0;
                ARIADNE_TESTING_EXPECT_EQUAL(truth,expected);
            }
        }
        const unsigned expected_rows=logical?300u:285u;
        if (ids.size()!=expected_rows || *ids.begin()!=0u || *ids.rbegin()!=expected_rows-1u)
            throw std::runtime_error("Incomplete "+name);
    }
}

void TestDyadic::test_dyadic_bounds_power_oracle() {
    std::ifstream data(ARIADNE_TEST_DATA_DIR "/dyadic_bounds_power.tsv");
    if (!data) throw std::runtime_error("Missing dyadic_bounds_power.tsv");
    std::string line;
    std::set<unsigned> ids;
    while (std::getline(data,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        unsigned id,degree;
        Int l,u,ln,ld,un,ud;
        std::string extra;
        if (!(fields>>id>>l>>u>>degree>>ln>>ld>>un>>ud) || (fields>>extra)
            || !ids.insert(id).second)
            throw std::runtime_error("Malformed dyadic_bounds_power.tsv");
        const DyadicBounds input(Dyadic(l,1u),Dyadic(u,1u));
        for (auto result:{pow(input,static_cast<Nat>(degree)),pow(input,static_cast<Int>(degree))}) {
            ARIADNE_TESTING_EXPECT_EQUAL(Rational(result.lower()),Rational(ln,ld));
            ARIADNE_TESTING_EXPECT_EQUAL(Rational(result.upper()),Rational(un,ud));
        }
    }
    if (ids.size()!=25u || *ids.begin()!=0u || *ids.rbegin()!=24u)
        throw std::runtime_error("Incomplete dyadic_bounds_power.tsv");
    ARIADNE_TESTING_EXPECT_EXCEPTION(pow(DyadicBounds(Dyadic(1)),Int(-1)));
}

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("dyadic");
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_literal);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_conversions);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_arithmetic);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_bounds);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_ball);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_comparisons);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_infinity);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_native_wrappers);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_dyadic_add_sub_oracles);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_dyadic_conversion_oracle);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_dyadic_public_contracts);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_dyadic_mixed_bounds_oracles);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, test_dyadic_bounds_power_oracle);
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
