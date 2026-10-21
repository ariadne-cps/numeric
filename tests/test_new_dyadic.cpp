/***************************************************************************
 * test_new_dyadic.cpp
 *
 * Parallel port of tests/test_dyadic.cpp. Each numbered section follows
 * the legacy Dyadic test's literal/conversion/arithmetic/comparison/
 * infinity coverage. Dyadic remains the reference and is not modified.
 ****************************************************************************/
#include "numeric/dyadic.hpp"
#include "numeric/new_dyadic.hpp"
#include "numeric/decimal.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/rational.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace Ariadne;

namespace {
void check(bool condition,char const* description) {
    if (!condition) { throw std::runtime_error(description); }
}
template<class Action> void must_fail(Action&& action,char const* description) {
    bool thrown=false;
    try { action(); } catch (std::exception const&) { thrown=true; }
    check(thrown,description);
}
template<class A,class B>
void same_value(A const& actual,B const& expected,char const* description) {
    check(actual.to_rational()==Rational(expected),description);
}

// TestDyadic::test_literal: exact dyadic decimal literals, powers of two,
// negative numbers and rejection of non-dyadic decimals.
void test_literal() {
    check(NewDyadic(String("3.25"))==NewDyadic(13,2u),"3.25");
    check(NewDyadic(String("-11.375"))==NewDyadic(-91,3u),"-11.375");
    check(NewDyadic(String("0.375"))==NewDyadic(3,3u),"0.375");
    check(NewDyadic(TwoExp(3))==NewDyadic(8),"2^3");
    check(NewDyadic(TwoExp(-3))==NewDyadic(1,3u),"2^-3");
    check(NewDyadic(TwoExp(-4))==NewDyadic(1,4u),"2^-4");
    check(NewDyadic(String("3.25")).literal()==String("3.25"),
          "exact decimal formatting");
    must_fail([] { (void)NewDyadic(String("not-a-dyadic")); },
              "reject invalid literal");
    must_fail([] { (void)NewDyadic(String("0.1")); },
              "reject non-dyadic literal");
}

// TestDyadic::test_conversions: copy, integral/decimal powers,
// round/floor/ceil, non-finite errors. Preserve each original integer
// rounding expectation, including negative and exact-half boundaries.
void test_conversions() {
    NewDyadic assigned;
    NewDyadic source(7,2u);
    assigned=source;
    check(assigned==source,"copy assignment");
    check(NewDyadic(Integer(-3))==NewDyadic(-3,0u),"Integer");
    check(NewDyadic(NewDyadic(-13))==NewDyadic(-13),"copy constructor");
    check(NewDyadic(NewDyadic(-13,3u))==NewDyadic(-13,3u),
          "copy fractional");
    check(NewDyadic(Integer(12),3u)==NewDyadic(3,1u),
          "canonical normalization");
    check(NewDyadic(Integer(-12),3u)==NewDyadic(-3,1u),
          "negative normalization");
    check(NewDyadic(0,100u)==NewDyadic(0),"canonical zero");

    struct RoundCase { int numerator; unsigned exponent; int expected; };
    constexpr RoundCase rounds[]={
        {-11,2u,-3},{-10,2u,-3},{-9,2u,-2},
        {-2,2u,-1},{-1,2u,0},{0,2u,0},{1,2u,0},
        {2,2u,1},{9,2u,2},{10,2u,3},{11,2u,3}
    };
    for (auto const& x:rounds) {
        check(round(NewDyadic(x.numerator,x.exponent))==
              Integer(x.expected),"legacy round case");
        check(round(NewDyadic(x.numerator,x.exponent))==
              round(Dyadic(x.numerator,x.exponent)),
              "legacy round oracle");
    }
    check(floor(NewDyadic(-5,1u))==Integer(-3),"negative floor");
    check(ceil(NewDyadic(-5,1u))==Integer(-2),"negative ceil");
    must_fail([] { (void)round(NewDyadic::inf()); },"round infinity");
    must_fail([] { (void)floor(NewDyadic::inf()); },"floor infinity");
    must_fail([] { (void)ceil(NewDyadic::inf()); },"ceil infinity");
    check(NewDyadic(7).mantissa()==Integer(7),"mantissa");
    check(NewDyadic(7,2u).exponent()==2,"exponent");

    for (unsigned width:{53u,128u,256u,1024u,2048u,4096u}) {
        const BigUnsigned high=shift_left(BigUnsigned(1u),width-1u);
        const Integer a=Integer::from_magnitude(
            add(high,BigUnsigned(13u)),true);
        const Integer c=Integer::from_magnitude(
            add(high,BigUnsigned(31u)),false);
        const unsigned e=width/5u+1u;
        const NewDyadic x(a,e),y(c,e+1u);
        const Dyadic old_x(a,e),old_y(c,e+1u);
        same_value(x,old_x,"wide constructor x");
        same_value(y,old_y,"wide constructor y");
        check(floor(x)==floor(old_x),"wide floor");
        check(ceil(x)==ceil(old_x),"wide ceil");
        check(round(x)==round(old_x),"wide round");
        check(x.get_d()==x.to_rational().get_d(),
              "native Rational binary64 conversion");
    }
}

// TestDyadic::test_arithmetic: use the exact values, signed powers,
// named arithmetic and in-place operators from the original test.
void test_arithmetic() {
    NewDyadic accumulated=1;
    accumulated+=NewDyadic(2);
    check(accumulated==NewDyadic(3),"compound plus");
    check(NewDyadic(3,2u)+NewDyadic(-5,3u)==NewDyadic(1,3u),
          "addition");
    check(NewDyadic(3,2u)-NewDyadic(-5,3u)==NewDyadic(11,3u),
          "subtraction");
    check(NewDyadic(3,2u)*NewDyadic(-5,3u)==NewDyadic(-15,5u),
          "multiplication");
    check(pow(NewDyadic(2),static_cast<Int>(3))==Rational(8),
          "positive signed power");
    check(pow(NewDyadic(2),static_cast<Int>(-3))==Rational(1,8),
          "negative signed power");
    check(sqr(NewDyadic(-3,2u))==NewDyadic(9,4u),"positive square");
    check(hlf(NewDyadic(3))==NewDyadic(3,1u),"half");
    check(abs(NewDyadic(-3,1u))==NewDyadic(3,1u),"positive abs");
    check(fma(NewDyadic(3),NewDyadic(2),NewDyadic(1))==
          NewDyadic(7),"fused multiply-add");
    NewDyadic m(-3,2u);
    m*=NewDyadic(2);
    check(m==NewDyadic(-3,1u),"compound multiply");
    m-=NewDyadic(-3,1u);
    check(m==NewDyadic(0),"compound subtract");
    check(nul(NewDyadic(2))==NewDyadic(0),"nul");
    check(pos(NewDyadic(2))==NewDyadic(2),"pos");
    check(neg(NewDyadic(2))==NewDyadic(-2),"neg");

    for (unsigned width:{128u,256u,512u,1024u,2048u,4096u}) {
        const auto high=shift_left(BigUnsigned(1u),width-1u);
        const Integer a=Integer::from_magnitude(add(high,BigUnsigned(13u)),false);
        const Integer b=Integer::from_magnitude(add(high,BigUnsigned(31u)),true);
        const unsigned e=width/5u+1u;
        const NewDyadic x(a,e),y(b,e+1u);
        const Dyadic old_x(a,e),old_y(b,e+1u);
        same_value(x+y,old_x+old_y,"wide add");
        same_value(x-y,old_x-old_y,"wide sub");
        same_value(x*y,old_x*old_y,"wide mul");
        same_value(pow(x,Nat(3u)),pow(old_x,Nat(3u)),"wide cube");
    }
}

// TestDyadic::test_comparisons: same signed comparison and ordering
// scenarios, plus a compile-time signature check for Boolean/Comparison.
void test_comparisons() {
    static_assert(std::is_same_v<
        decltype(std::declval<NewDyadic>()<std::declval<NewDyadic>()),
        Boolean>);
    static_assert(std::is_same_v<
        decltype(cmp(std::declval<NewDyadic>(),std::declval<NewDyadic>())),
        Comparison>);
    check(NewDyadic(-4,3u)<NewDyadic(-1,2u), "negative comparison");
    check(cmp(NewDyadic(-4,3u),NewDyadic(-1,2u))==Comparison::LESS,
          "negative cmp");
    check(eq(NewDyadic(2,2u),NewDyadic(1,1u)),"canonical equality");
    check(lt(NewDyadic(-2),NewDyadic(1)),"lt");
    check(min(NewDyadic(-2),NewDyadic(1))==NewDyadic(-2),"min");
    check(max(NewDyadic(-2),NewDyadic(1))==NewDyadic(1),"max");
}

// TestDyadic::test_infinity: public special-value contracts, sign,
// mixed arithmetic, half, abs, min/max and double conversion.
void test_infinity() {
    const NewDyadic positive=NewDyadic::inf();
    const NewDyadic negative=NewDyadic::inf(Sign::NEGATIVE);
    const NewDyadic missing=NewDyadic::nan();
    check(is_nan(missing),"NaN");
    check(is_inf(positive) && is_inf(negative),"infinities");
    check(!is_inf(missing),"NaN not infinity");
    check(is_finite(NewDyadic(0)) && is_zero(NewDyadic(0)),"zero");
    check(!is_zero(positive),"infinity not zero");
    check(positive==NewDyadic::inf(Sign::POSITIVE),"positive infinity");
    check(positive>NewDyadic(0),"positive infinity ordering");
    check(negative<NewDyadic(0),"negative infinity ordering");
    check(negative<positive,"infinity ordering");
    check(sgn(positive)==Sign::POSITIVE,"positive sign");
    check(sgn(negative)==Sign::NEGATIVE,"negative sign");
    check(is_nan(positive+negative),"opposite infinities");
    check(positive+positive==positive,"same infinity addition");
    check(is_nan(positive-positive),"infinity subtraction");
    check(is_nan(positive*NewDyadic(0)),"zero times infinity");
    check(positive*negative==negative,"signed infinity multiply");
    check(negative*negative==positive,"negative negative multiply");
    check(is_nan(hlf(missing)),"half NaN");
    check(hlf(positive)==positive && hlf(negative)==negative,
          "half infinity");
    check(is_nan(abs(missing)),"absolute NaN");
    check(abs(negative)==positive,"absolute infinity");
    check(is_nan(max(missing,NewDyadic(0))),"NaN max");
    check(max(negative,NewDyadic(-2))==NewDyadic(-2),"max negative");
    check(min(negative,NewDyadic(-2))==negative,"min negative");
    check(std::isinf(positive.get_d()) && positive.get_d()>0,"double infinity");
    check(std::isinf(negative.get_d()) && negative.get_d()<0,"double -infinity");
    check(std::isnan(missing.get_d()),"double NaN");
}
}

int main() {
    try {
        test_literal();
        test_conversions();
        test_arithmetic();
        test_comparisons();
        test_infinity();
        std::cout<<"NewDyadic test port (literal/conversions/arithmetic/"
                    "comparisons/infinity) PASS\n";
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"NewDyadic test port FAILED: "<<error.what()<<'\n';
        return 1;
    }
}
