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

#include <iomanip>

#include "utility/test.hpp"

using namespace std;
using namespace Ariadne;


class TestDyadic
{
  public:
    void test();
  private:
    void test_concept();
    void test_literal();
    void test_conversions();
    void test_arithmetic();
    void test_bounds();
    void test_ball();
    void test_comparisons();
    void test_infinity();
};

void TestDyadic::test()
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_literal());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_bounds());
    ARIADNE_TEST_CALL(test_ball());
    ARIADNE_TEST_CALL(test_comparisons());
    ARIADNE_TEST_CALL(test_infinity());
}

void TestDyadic::test_concept() {
    unsigned int m=1; unsigned long int lm=1; int n=-2; long int ln=-2; Integer z=-5;
    Dyadic w, w2; Boolean b;

    w=Dyadic(); w=Dyadic(m); w=Dyadic(lm); w=Dyadic(n); w=Dyadic(ln); w=Dyadic(z); w=Dyadic(z);
    w2=Dyadic();
    w=m; w=lm; w=n; w=ln; w=z; w=w2;

    w=+w; w=-w;
    w=w+w; w=w-w; w=w*w;

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
    mpf_t raw_mpf;
    mpf_init_set_ui(raw_mpf,7u);
    ARIADNE_TEST_EQUAL(Dyadic(raw_mpf),Dyadic(7));
    mpf_clear(raw_mpf);

    Dyadic assigned;
    Dyadic assign_source(7,2u);
    assigned=assign_source;
    ARIADNE_TEST_EQUAL(assigned,assign_source);

    ARIADNE_TEST_EQUAL(Dyadic(Integer(-3)),Dyadic(-3,0u));
    ARIADNE_TEST_EQUAL(Dyadic(Dyadic(-13)),Dyadic(-13));
    ARIADNE_TEST_EQUAL(Dyadic(Dyadic(-13,3u)),Dyadic(-13,3u));

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
    ARIADNE_TEST_EQUALS(abs(DyadicBall(Dyadic(0),Dyadic(2))).value_raw(),Dyadic(2));

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

void TestDyadic::test_infinity() {

    mpf_t mpf; static_assert(std::is_signed<decltype(mpf[0]._mp_exp)>::value,"");

    mp_exp_t mp_exp_t_min = std::numeric_limits<mp_exp_t>::min();
    ARIADNE_TEST_PRINT(mp_exp_t_min);
    Dyadic z(0);
    assert(z._mpf[0]._mp_size==0 && z._mpf[0]._mp_exp==0);

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


int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    ARIADNE_TEST_CLASS(Dyadic,TestDyadic());

    return ARIADNE_TEST_FAILURES;
}
