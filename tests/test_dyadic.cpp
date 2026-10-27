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

#include "macro/testing.hpp"

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



class TestDyadic
{
  public:
    void test();
  private:
    void test_literal();
    void test_conversions();
    void test_arithmetic();
    void test_bounds();
    void test_ball();
    void test_comparisons();
    void test_infinity();
    void test_native_wrappers();
};

void TestDyadic::test()
{
    ARIADNE_TEST_CALL(test_literal());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_bounds());
    ARIADNE_TEST_CALL(test_ball());
    ARIADNE_TEST_CALL(test_comparisons());
    ARIADNE_TEST_CALL(test_infinity());
    ARIADNE_TEST_CALL(test_native_wrappers());
}

const Writer<Dyadic> fraction_write=FractionWriter();

void TestDyadic::test_literal() {
    ARIADNE_TEST_CONSTRUCT(Dyadic,q,(3.25_q2));
    ARIADNE_TEST_EQUALS(q,Dyadic(13,2u));
    ARIADNE_TEST_EQUALS(3.25_q2,Dyadic(13,2u));
    ARIADNE_TEST_EQUALS(-11.375_q2,Dyadic(-91,3u));
    ARIADNE_TEST_EQUALS(0.375_q2,Dyadic(3,3u));

    ARIADNE_TEST_EQUALS(_2^3,Dyadic(8));
    ARIADNE_TEST_EQUALS(_2^-3,Dyadic(1,3u));
    ARIADNE_TEST_EQUALS(+(_2^-4),Dyadic(1,4u));
    ARIADNE_TEST_EQUALS(-(_2^-4),Dyadic(-1,4u));
    ARIADNE_TEST_EQUALS(5*(_2^-3),Dyadic(5,3u));
    ARIADNE_TEST_EQUALS(5/(_2^3),Dyadic(5,3u));
    ARIADNE_TEST_EQUALS(5/(_2^-3),Dyadic(40));
// The following should not compile, since ^ binds less tightly than -,*,/
//    ARIADNE_TEST_EQUALS(-_2^-4,Dyadic(-1,4u));
//    ARIADNE_TEST_EQUALS(5*_2^-4,Dyadic(5,4u));
//    ARIADNE_TEST_EQUALS(5/_2^4,Dyadic(5,4u));

    DecimalWriter decimal_write;

    ARIADNE_TEST_PRINT(q);
    ARIADNE_TEST_EXECUTE(Dyadic::set_default_writer(fraction_write));
    ARIADNE_TEST_PRINT(q);
    ARIADNE_TEST_EXECUTE(Dyadic::set_default_writer(fraction_write));
    ARIADNE_TEST_EXECUTE(std::cout<<q);
    ARIADNE_TEST_EXECUTE(std::cout<<decimal_write(q));
    ARIADNE_TEST_EXECUTE(std::cout<<fraction_write(q));
    ARIADNE_TEST_EXECUTE(Dyadic::set_default_writer(decimal_write));
    ARIADNE_TEST_EXECUTE(std::cout<<q);

    RepresentationWriter<Dyadic> representation_write;
    ARIADNE_TEST_EXECUTE(std::cout<<representation_write(q));

    ARIADNE_TEST_EQUALS(to_string(fraction_write(q)),"13/2^2");
    ARIADNE_TEST_EQUALS(to_string(decimal_write(q)),"3.25");
    ARIADNE_TEST_EQUALS(to_string(representation_write(q)),"Dyadic(13,2u)");
    ARIADNE_TEST_EQUALS(q.literal(),String("3.25"));
    ARIADNE_TEST_EXECUTE(Dyadic::default_writer());
    ARIADNE_TEST_EQUALS(1.5_bin,Dyadic(3,1u));

    ScientificWriter scientific_write;
    ARIADNE_TEST_EXECUTE(std::cout<<scientific_write(Dyadic(-8)));
    ARIADNE_TEST_EXECUTE(std::cout<<scientific_write(Dyadic(1,2u)));
    ARIADNE_TEST_EXECUTE(std::cout<<scientific_write(Dyadic::inf()));
    ARIADNE_TEST_EXECUTE(std::cout<<decimal_write(Dyadic::inf()));
    ARIADNE_TEST_EXECUTE(std::cout<<decimal_write(-Dyadic::inf()));
    ARIADNE_TEST_EXECUTE(std::cout<<decimal_write(Dyadic::nan()));
    ARIADNE_TEST_EQUALS(to_string(fraction_write(Dyadic(3))),"3");
    ARIADNE_TEST_EQUALS(to_string(fraction_write(Dyadic(3,1u))),"3/2");
    ARIADNE_TEST_EXECUTE(std::cout<<fraction_write(Dyadic::inf()));

    StringStream two_exp_stream;
    two_exp_stream << TwoExp(3);
    ARIADNE_TEST_EQUALS(two_exp_stream.str(),String("2^3"));

    ARIADNE_TEST_FAIL(Dyadic(String("not-a-dyadic")));
    ARIADNE_TEST_FAIL(Dyadic(String("0.1")));
}

void TestDyadic::test_conversions() {
    Dyadic assigned;
    Dyadic assign_source(7,2u);
    assigned=assign_source;
    ARIADNE_TEST_EQUAL(assigned,assign_source);

    ARIADNE_TEST_EQUAL(Dyadic(Integer(-3)),Dyadic(-3,0u));
    ARIADNE_TEST_EQUAL(Dyadic(Dyadic(-13)),Dyadic(-13));
    ARIADNE_TEST_EQUAL(Dyadic(Dyadic(-13,3u)),Dyadic(-13,3u));
    Natural huge_exponent{Integer{String("9223372036854775808")}};
    ARIADNE_TEST_FAIL((Dyadic{Integer(1),huge_exponent}));
    Natural valid_exponent{63u};
    ARIADNE_TEST_EQUALS(Dyadic(Integer(1),valid_exponent),Dyadic(1,63u));
    ARIADNE_TEST_FAIL((PositiveDyadic{Dyadic(-1)}));

    ARIADNE_TEST_EQUAL(round(Dyadic(-11,2u)),Integer(-3));
    ARIADNE_TEST_EQUAL(round(Dyadic(-10,2u)),Integer(-3));
    ARIADNE_TEST_EQUAL(round(Dyadic(-9,2u)),Integer(-2));
    ARIADNE_TEST_EQUAL(round(Dyadic(-2,2u)),Integer(-1));
    ARIADNE_TEST_EQUAL(round(Dyadic(-1,2u)),Integer(0));
    ARIADNE_TEST_EQUAL(round(Dyadic(0,2u)),Integer(0));
    ARIADNE_TEST_EQUAL(round(Dyadic(1,2u)),Integer(0));
    ARIADNE_TEST_EQUAL(round(Dyadic(2,2u)),Integer(1));
    ARIADNE_TEST_EQUAL(round(Dyadic(9,2u)),Integer(2));
    ARIADNE_TEST_EQUAL(round(Dyadic(10,2u)),Integer(3));
    ARIADNE_TEST_EQUAL(round(Dyadic(11,2u)),Integer(3));
    ARIADNE_TEST_EQUAL(floor(Dyadic(-5,1u)),Integer(-3));
    ARIADNE_TEST_EQUAL(ceil(Dyadic(-5,1u)),Integer(-2));
    ARIADNE_TEST_FAIL(round(Dyadic::inf()));
    ARIADNE_TEST_FAIL(floor(Dyadic::inf()));
    ARIADNE_TEST_FAIL(ceil(Dyadic::inf()));

    ARIADNE_TEST_EQUAL(static_cast<Int>(3)/two,Dyadic(3,1u));
    ARIADNE_TEST_EQUAL(Integer(3)/two,Dyadic(3,1u));
    ARIADNE_TEST_EQUAL(ExactDouble(3)/two,Dyadic(3,1u));
    ARIADNE_TEST_EQUAL(Dyadic(3)/two,Dyadic(3,1u));

    ARIADNE_TEST_EQUAL(make_dyadic(7u),Dyadic(7));
    ARIADNE_TEST_EQUAL(make_dyadic(1ULL<<40),Dyadic(TwoExp(40)));
    ARIADNE_TEST_EQUALS(class_name<Dyadic>(),String("Dyadic"));
    ARIADNE_TEST_EQUALS(class_name<DyadicBounds>(),String("DyadicBounds"));
}

void TestDyadic::test_arithmetic() {
    Dyadic accumulated=1;
    accumulated+=Dyadic(2);
    ARIADNE_TEST_EQUALS(accumulated,Dyadic(3));
    ARIADNE_TEST_EQUAL(Dyadic(3,2u)+Dyadic(-5,3u),Dyadic(1,3u));
    ARIADNE_TEST_EQUAL(Dyadic(3,2u)-Dyadic(-5,3u),Dyadic(11,3u));
    ARIADNE_TEST_EQUAL(Dyadic(3,2u)*Dyadic(-5,3u),Dyadic(-15,5u));
    ARIADNE_TEST_EQUAL(pow(Dyadic(2),static_cast<Int>(3)),Rational(8));
    ARIADNE_TEST_EQUAL(pow(Dyadic(2),static_cast<Int>(-3)),Rational(1,8));
}

void TestDyadic::test_bounds() {
    DyadicBounds positive(1,3);
    DyadicBounds negative(-3,-1);
    DyadicBounds crossing(-2,3);
    DyadicBounds singleton(2,2);

    ARIADNE_TEST_EQUALS((+positive).lower_raw(),Dyadic(1));
    ARIADNE_TEST_EQUALS((-positive).lower_raw(),Dyadic(-3));
    ARIADNE_TEST_EQUALS((positive+DyadicBounds(2,4)).upper_raw(),Dyadic(7));
    ARIADNE_TEST_EQUALS((positive-DyadicBounds(2,4)).lower_raw(),Dyadic(-3));
    ARIADNE_TEST_EQUALS(add(positive,DyadicBounds(2,4)).lower_raw(),Dyadic(3));
    ARIADNE_TEST_EQUALS(sub(positive,DyadicBounds(2,4)).upper_raw(),Dyadic(1));
    ARIADNE_TEST_EQUALS(nul(positive).lower_raw(),Dyadic(0));
    ARIADNE_TEST_EQUALS(pos(positive).upper_raw(),Dyadic(3));
    ARIADNE_TEST_EQUALS(neg(positive).upper_raw(),Dyadic(-1));
    ARIADNE_TEST_EQUALS(sqr(positive).lower_raw(),Dyadic(1));
    ARIADNE_TEST_EQUALS(sqr(negative).upper_raw(),Dyadic(9));
    ARIADNE_TEST_EQUALS(sqr(crossing).lower_raw(),Dyadic(0));
    ARIADNE_TEST_EQUALS(sqr(crossing).upper_raw(),Dyadic(9));
    ARIADNE_TEST_EQUALS(hlf(positive).upper_raw(),Dyadic(3,1u));
    ARIADNE_TEST_EQUALS(abs(crossing).lower_raw(),Dyadic(0));
    ARIADNE_TEST_EQUALS(max(positive,DyadicBounds(2,4)).lower_raw(),Dyadic(2));
    ARIADNE_TEST_EQUALS(min(positive,DyadicBounds(2,4)).upper_raw(),Dyadic(3));

    ARIADNE_TEST_ASSERT(definitely(sgn(positive)));
    ARIADNE_TEST_ASSERT(definitely(!sgn(negative)));
    ARIADNE_TEST_ASSERT(is_indeterminate(sgn(crossing)));

    DyadicBounds left(1,2);
    DyadicBounds right(3,4);
    DyadicBounds overlap(2,3);
    ARIADNE_TEST_ASSERT(definitely(singleton==singleton));
    ARIADNE_TEST_ASSERT(definitely(!(left==right)));
    ARIADNE_TEST_ASSERT(is_indeterminate(left==overlap));
    ARIADNE_TEST_ASSERT(definitely(left!=right));
    ARIADNE_TEST_ASSERT(definitely(right!=left));
    ARIADNE_TEST_ASSERT(definitely(!(singleton!=singleton)));
    ARIADNE_TEST_ASSERT(is_indeterminate(left!=overlap));

    ARIADNE_TEST_ASSERT(consistent(left,overlap));
    ARIADNE_TEST_ASSERT(not consistent(left,right));
    ARIADNE_TEST_ASSERT(inconsistent(left,right));
    ARIADNE_TEST_ASSERT(inconsistent(right,left));
    ARIADNE_TEST_ASSERT(not inconsistent(left,overlap));
    ARIADNE_TEST_ASSERT(refines(DyadicBounds(2,2),DyadicBounds(1,3)));
    ARIADNE_TEST_ASSERT(not refines(DyadicBounds(0,2),DyadicBounds(1,3)));
    ARIADNE_TEST_ASSERT(not refines(DyadicBounds(2,4),DyadicBounds(1,3)));
    ARIADNE_TEST_EQUALS(refinement(DyadicBounds(1,3),DyadicBounds(2,4)).lower_raw(),Dyadic(2));
    ARIADNE_TEST_EQUALS(coarsening(DyadicBounds(1,3),DyadicBounds(2,4)).upper_raw(),Dyadic(4));

    ARIADNE_TEST_EQUALS(pow(crossing,Nat(2u)).upper_raw(),Dyadic(9));
    ARIADNE_TEST_EQUALS(pow(crossing,Nat(3u)).lower_raw(),Dyadic(-8));

    DyadicLowerBound lower(Dyadic(1));
    DyadicLowerBound lower_from_bounds(positive);
    DyadicUpperBound upper(Dyadic(3));
    DyadicUpperBound upper_from_bounds(positive);
    DyadicApproximation approximation(Dyadic(2));
    ARIADNE_TEST_EQUALS(lower.raw(),Dyadic(1));
    ARIADNE_TEST_EQUALS(lower_from_bounds.raw(),Dyadic(1));
    ARIADNE_TEST_EQUALS(upper.raw(),Dyadic(3));
    ARIADNE_TEST_EQUALS(upper_from_bounds.raw(),Dyadic(3));
    ARIADNE_TEST_EQUALS(approximation.raw(),Dyadic(2));
    ARIADNE_TEST_EXECUTE(std::cout<<lower<<" "<<upper);
    ARIADNE_TEST_EXECUTE((PositiveDyadicBounds{positive}));
    ARIADNE_TEST_EXECUTE((PositiveDyadicLowerBound{lower}));
    ARIADNE_TEST_EXECUTE((PositiveDyadicUpperBound{upper}));
    ARIADNE_TEST_EXECUTE((PositiveDyadicApproximation{approximation}));
    ARIADNE_TEST_EXECUTE(std::cout<<positive);
}

void TestDyadic::test_ball() {
    DyadicBall exact(Dyadic(2));
    DyadicBall uncertain(Dyadic(3),Dyadic(1));
    DyadicBall from_bounds(DyadicBounds(1,5));
    ARIADNE_TEST_EQUALS(exact.value(),Dyadic(2));
    ARIADNE_TEST_EQUALS(exact.error(),Dyadic(0));
    ARIADNE_TEST_EQUALS(uncertain.value_raw(),Dyadic(3));
    ARIADNE_TEST_EQUALS(uncertain.error_raw(),Dyadic(1));
    ARIADNE_TEST_EQUALS(from_bounds.value(),Dyadic(3));
    ARIADNE_TEST_EQUALS(from_bounds.error(),Dyadic(2));

    DyadicBounds converted=static_cast<DyadicBounds>(uncertain);
    ARIADNE_TEST_EQUALS(converted.lower_raw(),Dyadic(2));
    ARIADNE_TEST_EQUALS(converted.upper_raw(),Dyadic(4));

    ARIADNE_TEST_EQUALS((exact+uncertain).value_raw(),Dyadic(5));
    ARIADNE_TEST_EQUALS((uncertain-exact).value_raw(),Dyadic(1));
    ARIADNE_TEST_EQUALS((exact*uncertain).value_raw(),Dyadic(6));
    ARIADNE_TEST_EQUALS(abs(DyadicBall(Dyadic(-3),Dyadic(1))).value_raw(),Dyadic(3));
    DyadicBall crossing_abs=abs(DyadicBall(Dyadic(0),Dyadic(2)));
    ARIADNE_TEST_EQUALS(crossing_abs.value_raw(),Dyadic(1));
    ARIADNE_TEST_EQUALS(crossing_abs.error_raw(),Dyadic(1));

    ARIADNE_TEST_ASSERT(definitely(DyadicBall(Dyadic(1),Dyadic(0)) < DyadicBall(Dyadic(3),Dyadic(0))));
    ARIADNE_TEST_ASSERT(definitely(!(DyadicBall(Dyadic(3),Dyadic(0)) < DyadicBall(Dyadic(1),Dyadic(0)))));
    ARIADNE_TEST_ASSERT(is_indeterminate(DyadicBall(Dyadic(2),Dyadic(1)) < DyadicBall(Dyadic(3),Dyadic(1))));
    ARIADNE_TEST_ASSERT(inconsistent(DyadicBall(Dyadic(1),Dyadic(0)),DyadicBall(Dyadic(3),Dyadic(0))));
    ARIADNE_TEST_ASSERT(not inconsistent(DyadicBall(Dyadic(2),Dyadic(1)),DyadicBall(Dyadic(3),Dyadic(1))));
    ARIADNE_TEST_ASSERT(refines(DyadicBall(Dyadic(2),Dyadic(1)),DyadicBall(Dyadic(2),Dyadic(2))));
    ARIADNE_TEST_ASSERT(not refines(DyadicBall(Dyadic(2),Dyadic(2)),DyadicBall(Dyadic(2),Dyadic(1))));
    ARIADNE_TEST_EXECUTE(refinement(DyadicBall(Dyadic(2),Dyadic(1)),DyadicBall(Dyadic(2),Dyadic(2))));
    ARIADNE_TEST_EXECUTE(coarsening(DyadicBall(Dyadic(2),Dyadic(1)),DyadicBall(Dyadic(3),Dyadic(1))));
    ARIADNE_TEST_EXECUTE(std::cout<<uncertain);
}

void TestDyadic::test_comparisons() {
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Dyadic(-4,3u),Dyadic(-1,2u));
}

// D48: native exact bounds and ball arithmetic, without GMP storage.
void TestDyadic::test_native_wrappers() {
    const DyadicBounds positive(Dyadic(1),Dyadic(3));
    const DyadicBounds negative(Dyadic(-3),Dyadic(-1));
    const DyadicBounds crossing(Dyadic(-2),Dyadic(3));
    ARIADNE_TEST_ASSERT((positive+DyadicBounds(Dyadic(2),Dyadic(4))).upper_raw()
          ==Dyadic(7));
    ARIADNE_TEST_ASSERT((positive-negative).lower_raw()==Dyadic(2));
    ARIADNE_TEST_ASSERT((positive*negative).lower_raw()==Dyadic(-9));
    ARIADNE_TEST_ASSERT((positive*negative).upper_raw()==Dyadic(-1));
    ARIADNE_TEST_ASSERT(sqr(crossing).lower_raw()==Dyadic(0));
    ARIADNE_TEST_ASSERT(sqr(crossing).upper_raw()==Dyadic(9));
    ARIADNE_TEST_ASSERT(pow(crossing,Nat(3u)).lower_raw()==Dyadic(-8));
    ARIADNE_TEST_ASSERT(pow(crossing,Nat(2u)).lower_raw()==Dyadic(0));
    ARIADNE_TEST_ASSERT(abs(crossing).upper_raw()==Dyadic(3));
    ARIADNE_TEST_ASSERT(refines(DyadicBounds(Dyadic(2),Dyadic(2)),positive));
    const DyadicBall ball(Dyadic(3),Dyadic(1));
    const DyadicBounds enclosure=static_cast<DyadicBounds>(ball);
    ARIADNE_TEST_ASSERT(enclosure.lower_raw()==Dyadic(2) &&
          enclosure.upper_raw()==Dyadic(4));
    ARIADNE_TEST_ASSERT((ball*DyadicBall(Dyadic(2),Dyadic(0))).value_raw()
          ==Dyadic(6));
    const DyadicLowerBound lower(positive);
    const DyadicUpperBound upper(positive);
    const DyadicApproximation approximation(Dyadic(2));
    ARIADNE_TEST_ASSERT(lower.raw()==Dyadic(1) && upper.raw()==Dyadic(3)
          && approximation.raw()==Dyadic(2));
}


void TestDyadic::test_infinity() {

    // Verify the new canonical zero representation without GMP internals.
    const Dyadic z(0);
    ARIADNE_TEST_ASSERT(z.coefficient()==Integer(0));
    ARIADNE_TEST_ASSERT(z.denominator_exponent()==0u);
    ARIADNE_TEST_ASSERT(Dyadic(0,17u).denominator_exponent()==0u);

    ARIADNE_TEST_PRINT(Dyadic::nan());
    ARIADNE_TEST_PRINT(Dyadic::inf());
    ARIADNE_TEST_PRINT(Dyadic::inf(Sign(0)));
    ARIADNE_TEST_PRINT(Dyadic::inf(Sign(+1)));
    ARIADNE_TEST_PRINT(Dyadic::inf(Sign(-1)))

    ARIADNE_TEST_ASSERT(is_nan(Dyadic::nan()));
    ARIADNE_TEST_ASSERT(is_inf(Dyadic::inf()));
    ARIADNE_TEST_ASSERT(is_inf(Dyadic::inf(Sign(+1))));
    ARIADNE_TEST_ASSERT(is_inf(Dyadic::inf(Sign(-1))));
    ARIADNE_TEST_ASSERT(not is_inf(Dyadic::nan()));
    ARIADNE_TEST_ASSERT(not is_inf(Dyadic(0)));
    ARIADNE_TEST_ASSERT(is_finite(Dyadic(0)));
    ARIADNE_TEST_ASSERT(is_zero(Dyadic(0)));
    ARIADNE_TEST_ASSERT(not is_zero(Dyadic(1)));
    ARIADNE_TEST_ASSERT(not is_zero(Dyadic::inf()));
    ARIADNE_TEST_EQUALS(Dyadic::inf(Sign(+1)),Dyadic::inf());
    ARIADNE_TEST_ASSERT(Dyadic::inf(Sign(+1))>Dyadic(0));
    ARIADNE_TEST_ASSERT(Dyadic::inf(Sign(-1))<Dyadic(0));

    ARIADNE_TEST_BINARY_PREDICATE(operator==,Dyadic::inf(),Dyadic::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Dyadic(-1,2u),Dyadic::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Dyadic(0),Dyadic::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Dyadic(3,1u),Dyadic::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-Dyadic::inf(),Dyadic::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-Dyadic::inf(),Dyadic(1,2u));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-Dyadic::inf(),Dyadic(-3,1u));

    ExactDouble double_max(std::numeric_limits<double>::max());
    ExactDouble double_inf(std::numeric_limits<double>::infinity());
    ExactDouble double_nan(double_inf.get_d()*0.0);

    ARIADNE_TEST_BINARY_PREDICATE(operator<,Dyadic(double_max),Dyadic(double_inf));

    ARIADNE_TEST_ASSERT(Dyadic::inf(Sign::POSITIVE).get_d()==std::numeric_limits<double>::infinity());
    ARIADNE_TEST_ASSERT(Dyadic::inf(Sign::NEGATIVE).get_d()==-std::numeric_limits<double>::infinity());
    ARIADNE_TEST_ASSERT(isnan(Dyadic::inf(Sign::ZERO).get_d()));

    ARIADNE_TEST_ASSERT(Dyadic(+double_inf)==Dyadic::inf(Sign(+1)));
    ARIADNE_TEST_ASSERT(Dyadic(-double_inf)==Dyadic::inf(Sign(-1)));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic(double_nan)));

    ARIADNE_TEST_ASSERT(is_nan(-Dyadic::nan()));
    ARIADNE_TEST_EQUAL(-Dyadic::inf(Sign::POSITIVE),Dyadic::inf(Sign::NEGATIVE));
    ARIADNE_TEST_EQUAL(-Dyadic::inf(Sign::NEGATIVE),Dyadic::inf(Sign::POSITIVE));
    ARIADNE_TEST_EQUALS(sgn(Dyadic::inf()), Sign::POSITIVE);
    ARIADNE_TEST_EQUALS(sgn(-Dyadic::inf()), Sign::NEGATIVE);

    ARIADNE_TEST_ASSERT(is_nan(Dyadic::inf()+(-Dyadic::inf())));
    ARIADNE_TEST_EQUALS(Dyadic::inf()+Dyadic::inf(),Dyadic::inf());
    ARIADNE_TEST_EQUALS(Dyadic::inf()+Dyadic(-2),Dyadic::inf());

    ARIADNE_TEST_ASSERT(is_nan(Dyadic::inf()-Dyadic::inf()));
    ARIADNE_TEST_EQUALS(Dyadic(2)-Dyadic::inf(),-Dyadic::inf())

    ARIADNE_TEST_ASSERT(is_nan(Dyadic::inf(Sign(+1))*Dyadic(0)));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic::inf(Sign(-1))*Dyadic(0)));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic(0)*Dyadic::inf(Sign(+1))));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic(0)*Dyadic::inf(Sign(-1))));
    ARIADNE_TEST_EQUALS(Dyadic::inf(Sign(+1))*Dyadic::inf(Sign(+1)),Dyadic::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Dyadic::inf(Sign(+1))*Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Dyadic::inf(Sign(-1))*Dyadic::inf(Sign(+1)),Dyadic::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Dyadic::inf(Sign(-1))*Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Dyadic(+2)*Dyadic::inf(Sign(+1)),Dyadic::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Dyadic(+2)*Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Dyadic(-2)*Dyadic::inf(Sign(+1)),Dyadic::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Dyadic(-2)*Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Dyadic::inf(Sign(+1))*Dyadic(+2),Dyadic::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Dyadic::inf(Sign(+1))*Dyadic(-2),Dyadic::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Dyadic::inf(Sign(-1))*Dyadic(+2),Dyadic::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Dyadic::inf(Sign(-1))*Dyadic(-2),Dyadic::inf(Sign(+1)));

    ARIADNE_TEST_ASSERT(is_nan(hlf(Dyadic::nan())));
    ARIADNE_TEST_EQUALS(hlf(Dyadic::inf()),Dyadic::inf());
    ARIADNE_TEST_EQUALS(hlf(-Dyadic::inf()),-Dyadic::inf());

    ARIADNE_TEST_ASSERT(is_nan(abs(Dyadic::nan())));
    ARIADNE_TEST_EQUALS(abs(Dyadic::inf()),Dyadic::inf());
    ARIADNE_TEST_EQUALS(abs(-Dyadic::inf()),Dyadic::inf());

    ARIADNE_TEST_ASSERT(is_nan(max(Dyadic::nan(),Dyadic(0))));
    ARIADNE_TEST_ASSERT(is_nan(max(Dyadic::nan(),Dyadic::inf())));
    ARIADNE_TEST_EQUALS(max(Dyadic::inf(Sign(-1)),Dyadic::inf(Sign(-1))),Dyadic::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(max(Dyadic::inf(Sign(-1)),Dyadic(-2)),Dyadic(-2));
    ARIADNE_TEST_EQUALS(max(Dyadic(-2),Dyadic::inf(Sign(+1))),Dyadic::inf(Sign(+1)));

    ARIADNE_TEST_ASSERT(DyadicBounds(1,3)<=3);
    ARIADNE_TEST_ASSERT(DyadicBounds(1,3)>=1);
    ARIADNE_TEST_ASSERT(DyadicBounds(1,3)< 4);
    ARIADNE_TEST_ASSERT(DyadicBounds(1,3)> 0);
    ARIADNE_TEST_ASSERT(is_indeterminate(DyadicBounds(1,3)<=2));
    ARIADNE_TEST_ASSERT(is_indeterminate(DyadicBounds(1,3)>=2));
    ARIADNE_TEST_ASSERT(is_indeterminate(DyadicBounds(1,3)< 3));
    ARIADNE_TEST_ASSERT(is_indeterminate(DyadicBounds(1,3)> 1));
    ARIADNE_TEST_ASSERT(!(DyadicBounds(1,3)<=0));
    ARIADNE_TEST_ASSERT(!(DyadicBounds(1,3)>=4));
    ARIADNE_TEST_ASSERT(not(DyadicBounds(1,3)< 1));
    ARIADNE_TEST_ASSERT(not(DyadicBounds(1,3)> 3));
}


void test_dyadic_add_sub_oracles() {
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
            std::cout << name << " row " << id << std::endl;
            ARIADNE_TEST_EQUALS(a+b,sum_expected);
            ARIADNE_TEST_EQUALS(b+a,sum_expected);
            ARIADNE_TEST_EQUALS(a-b,difference_expected);
            ARIADNE_TEST_EQUALS(b-a,-difference_expected);
            // Canonicalization is observable through the public accessors.
            const auto sum=a+b,difference=a-b;
            ARIADNE_TEST_EQUALS(sum.coefficient().literal(),String(sn));
            ARIADNE_TEST_EQUALS(sum.denominator_exponent(),se);
            ARIADNE_TEST_EQUALS(difference.coefficient().literal(),String(dn));
            ARIADNE_TEST_EQUALS(difference.denominator_exponent(),de);
        }
        if (ids.size()!=256u || *ids.begin()!=0u || *ids.rbegin()!=255u)
            throw std::runtime_error("Incomplete "+name);
    }
}

void test_dyadic_conversion_oracle() {
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
        std::cout << "Dyadic conversion row " << id << std::endl;
        ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(source.get_d()),expected);
        if (exponent==0u) ARIADNE_TEST_EQUALS(source.literal(),String(coefficient));
        ARIADNE_TEST_EQUALS(floor(source),Integer(String(floor_expected)));
        ARIADNE_TEST_EQUALS(ceil(source),Integer(String(ceil_expected)));
        ARIADNE_TEST_EQUALS(hlf(source),Dyadic(Integer(String(coefficient)),static_cast<std::size_t>(exponent)+1u));
        ARIADNE_TEST_EQUALS(fma(source,Dyadic(1),Dyadic(0)),source);
        const Rational exact_value{Integer(String(coefficient)),Integer(String(denominator))};
        const Rational reciprocal{Integer(String(denominator)),Integer(String(coefficient))};
        ARIADNE_TEST_EQUALS(div(source,Integer(1)),exact_value);
        ARIADNE_TEST_EQUALS(div(Integer(1),source),reciprocal);
        const FloatMP mp(source,FloatMP::ROUND_TO_NEAREST,precision);
        ARIADNE_TEST_EQUALS(Dyadic(mp),source);
        const Approximation<Dyadic> approximate(source);
        ARIADNE_TEST_EQUALS(Dyadic(approximate.get(precision).raw()),source);
        if (exact!=0) {
            const double value=std::bit_cast<double>(expected);
            const FloatDP fp(ExactDouble(value),dp);
            ARIADNE_TEST_EQUALS(Dyadic(fp),source);
            ARIADNE_TEST_EQUALS(Dyadic(ExactDouble(value)),source);
            ARIADNE_TEST_EQUALS(Dyadic(value),source);
            ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(approximate.get(dp).raw().get_d()),expected);
        }
    }
    if (ids.size()!=168u || *ids.begin()!=0u || *ids.rbegin()!=167u)
        throw std::runtime_error("Incomplete dyadic_conversion.tsv");
}

void test_dyadic_public_contracts() {
    const std::size_t largest=std::numeric_limits<std::size_t>::max();
    const Dyadic tiny(Integer(1),largest);
    ARIADNE_TEST_FAIL(tiny.exponent());
    ARIADNE_TEST_FAIL(tiny.literal());
    ARIADNE_TEST_FAIL(tiny*Dyadic(1,1u));
    ARIADNE_TEST_FAIL(hlf(tiny));
    const Integer outside_size=Integer::from_magnitude(shift_left(BigUnsigned(1u),
        std::numeric_limits<std::size_t>::digits),false);
    const Natural too_large(outside_size);
    ARIADNE_TEST_FAIL((Dyadic(Integer(1),too_large)));
    ARIADNE_TEST_FAIL(Dyadic(1)/3);

    const PositiveDyadic default_positive;
    const Dyadic positive(3,1u),negative(-3,1u);
    const PositiveDyadic copied_positive(positive);
    ARIADNE_TEST_EQUALS(default_positive,Dyadic(0));
    ARIADNE_TEST_EQUALS(copied_positive,positive);
    ARIADNE_TEST_EQUALS(cast_positive(positive),positive);
    ARIADNE_TEST_FAIL(cast_positive(negative));
    ARIADNE_TEST_ASSERT(positive.kind()==Dyadic::Kind::Finite);
    ARIADNE_TEST_ASSERT(Dyadic::inf().kind()==Dyadic::Kind::PositiveInfinity);
    ARIADNE_TEST_ASSERT((-Dyadic::inf()).kind()==Dyadic::Kind::NegativeInfinity);
    ARIADNE_TEST_ASSERT(Dyadic::nan().kind()==Dyadic::Kind::NaN);
    ARIADNE_TEST_EQUALS(Dyadic::nan().literal(),String("NaN"));
    ARIADNE_TEST_EQUALS(Dyadic::inf().literal(),String("inf"));
    ARIADNE_TEST_EQUALS((-Dyadic::inf()).literal(),String("-inf"));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic(Rational::nan())));
    ARIADNE_TEST_EQUALS(Dyadic(Rational::inf()),Dyadic::inf());
    ARIADNE_TEST_EQUALS(Dyadic(-Rational::inf()),-Dyadic::inf());
    ARIADNE_TEST_ASSERT(is_nan(Dyadic::inf(Sign::ZERO)));
    ARIADNE_TEST_ASSERT(is_nan(-Dyadic::nan()));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic(0)*Dyadic::inf()));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic::inf()*Dyadic(0)));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic::nan()*positive));
    ARIADNE_TEST_ASSERT(is_nan(positive*Dyadic::nan()));
    ARIADNE_TEST_EQUALS(cmp(Int(1),positive),Comparison::LESS);
    ARIADNE_TEST_ASSERT(positive!=negative);
    ARIADNE_TEST_ASSERT(!eq(positive,Dyadic::nan()));
    ARIADNE_TEST_ASSERT(!eq(positive,Dyadic(3)));
    ARIADNE_TEST_ASSERT(!lt(positive,Dyadic::nan()));
    ARIADNE_TEST_ASSERT(!lt(Dyadic::nan(),positive));
    ARIADNE_TEST_ASSERT(is_nan(max(positive,Dyadic::nan())));

    ARIADNE_TEST_EQUALS(Dyadic(String("1.5")),positive);
    ARIADNE_TEST_EQUALS(dy_(1.5L),positive);
    std::ostringstream stream;
    stream << Approximation<Dyadic>(positive);
    ARIADNE_TEST_EQUALS(stream.str(),std::string("1.5"));
}

void test_dyadic_mixed_bounds_oracles() {
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
            std::cout << name << " row " << id << std::endl;
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
                ARIADNE_TEST_EQUALS(result.lower(),Rational(ln,ld));
                ARIADNE_TEST_EQUALS(result.upper(),Rational(un,ud));
                if (op=="div" && ql%qd==0 && qu%qd==0) {
                    const DyadicBounds other(Dyadic(ql/qd),Dyadic(qu/qd));
                    const auto direct=reverse?div(other,d):div(d,other);
                    ARIADNE_TEST_EQUALS(direct.lower(),Rational(ln,ld));
                    ARIADNE_TEST_EQUALS(direct.upper(),Rational(un,ud));
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
                ARIADNE_TEST_EQUALS(truth,expected);
            }
        }
        const unsigned expected_rows=logical?300u:285u;
        if (ids.size()!=expected_rows || *ids.begin()!=0u || *ids.rbegin()!=expected_rows-1u)
            throw std::runtime_error("Incomplete "+name);
    }
}

void test_dyadic_bounds_power_oracle() {
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
            ARIADNE_TEST_EQUALS(Rational(result.lower()),Rational(ln,ld));
            ARIADNE_TEST_EQUALS(Rational(result.upper()),Rational(un,ud));
        }
    }
    if (ids.size()!=25u || *ids.begin()!=0u || *ids.rbegin()!=24u)
        throw std::runtime_error("Incomplete dyadic_bounds_power.tsv");
    ARIADNE_TEST_FAIL(pow(DyadicBounds(Dyadic(1)),Int(-1)));
}

int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    ARIADNE_TEST_CLASS(Dyadic,TestDyadic());
    ARIADNE_TEST_CALL(test_dyadic_add_sub_oracles());
    ARIADNE_TEST_CALL(test_dyadic_conversion_oracle());
    ARIADNE_TEST_CALL(test_dyadic_public_contracts());
    ARIADNE_TEST_CALL(test_dyadic_mixed_bounds_oracles());
    ARIADNE_TEST_CALL(test_dyadic_bounds_power_oracle());

    return ARIADNE_TEST_FAILURES;
}
