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
#include "numeric/floatap.hpp"
#include <bit>
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
// R34: bit-for-bit legacy binary64 and exact cross-exponent cmp.
// R35: ported Dyadic arithmetic oracle over virtual-shift fast path
// and both fallback boundaries. Test magnitude reversal, both signs,
// exact cancellation, carry, borrow and +=/-= alias safety.
// R37: a normalized finite coefficient is odd whenever the
// denominator exponent is nonzero. Thus every positive exponent
// means a non-integral value, and floor/ceil must round away from
// truncation in the appropriate direction. Exercise both signs,
// full limb boundaries, and values that become integral on input.
void test_canonical_fraction_rounding() {
    for (unsigned bits:{1u,2u,53u,63u,64u,65u,127u,128u,
                        129u,255u,256u,257u,1024u,2048u,4096u,4097u}) {
        const BigUnsigned high=shift_left(BigUnsigned(1u),bits-1u);
        for (unsigned tail:{0u,1u,2u,3u,7u,63u}) {
            const Integer magnitude=Integer::from_magnitude(
                add(high,BigUnsigned(tail)),false);
            for (unsigned negative=0u;negative<2u;++negative) {
                const Integer numerator=negative!=0u ? -magnitude : magnitude;
                for (unsigned exponent:{0u,1u,2u,7u,63u,64u,
                                        bits/2u,bits,bits+1u}) {
                    const NewDyadic native(numerator,exponent);
                    const Dyadic legacy(numerator,exponent);
                    check(floor(native)==floor(legacy),
                          "R37 canonical fractional floor");
                    check(ceil(native)==ceil(legacy),
                          "R37 canonical fractional ceil");
                    check(floor(native)<=ceil(native),
                          "R37 floor/ceil order");
                    const bool integer=native.denominator_exponent()==0u;
                    check(static_cast<bool>(floor(native)==ceil(native))==integer,
                          "R37 exponent/whole-value invariant");
                }
            }
        }
    }
    for (unsigned exponent:{0u,1u,2u,63u,64u,65u,256u,4096u}) {
        const NewDyadic n(Integer(0),exponent);
        check(n.denominator_exponent()==0u,"R37 canonical zero");
        check(floor(n)==Integer(0) && ceil(n)==Integer(0),
              "R37 zero floor/ceil");
        if (exponent!=0u) {
            const BigUnsigned exact=shift_left(BigUnsigned(1u),exponent);
            const Integer exact_integer=Integer::from_magnitude(exact,false);
            const NewDyadic q(exact_integer,exponent);
            check(q.denominator_exponent()==0u,"R37 integral factor reduction");
            check(floor(q)==Integer(1) && ceil(q)==Integer(1),
                  "R37 exact whole power of two");
        }
    }
}

// R38: test direct owned add/sub on both sides of the ARM64
// 32-limb and x86-64 17-limb crossover. The original Dyadic is
// the independent oracle, not an input-specific benchmark expectation.
void test_owned_shifted_arithmetic() {
    for (unsigned bits:{960u,961u,1023u,1024u,1025u,
                        1087u,1088u,1089u,
                        1983u,1984u,1985u,
                        2047u,2048u,2049u,
                        4095u,4096u,4097u,8192u}) {
        const BigUnsigned high=shift_left(BigUnsigned(1u),bits-1u);
        for (unsigned low:{0u,1u,7u,63u}) {
            const Integer ax=Integer::from_magnitude(
                add(high,BigUnsigned(11u+low)),false);
            const Integer by=Integer::from_magnitude(
                add(high,BigUnsigned(13u+low)),false);
            for (unsigned signs=0u;signs<4u;++signs) {
                const Integer aa=(signs&1u)!=0u ? -ax : ax;
                const Integer bb=(signs&2u)!=0u ? -by : by;
                for (unsigned gap:{0u,1u,2u,7u,31u,63u,64u}) {
                    const unsigned ea=bits/5u+1u,eb=ea+gap;
                    const NewDyadic a(aa,ea),b(bb,eb);
                    const Dyadic ref_a(aa,ea),ref_b(bb,eb);
                    check((a+b).to_rational()==Rational(ref_a+ref_b),
                          "R38 owned shifted add oracle");
                    check((a-b).to_rational()==Rational(ref_a-ref_b),
                          "R38 owned shifted sub oracle");
                    check((b-a).to_rational()==Rational(ref_b-ref_a),
                          "R38 owned shifted reversed sub oracle");
                    NewDyadic inplace(a);
                    inplace+=b;
                    check(inplace.to_rational()==Rational(ref_a+ref_b),
                          "R38 owned shifted += alias");
                    inplace=a;
                    inplace-=b;
                    check(inplace.to_rational()==Rational(ref_a-ref_b),
                          "R38 owned shifted -= alias");
                }
            }
        }
    }
}

void test_aligned_add_sub() {
    for (unsigned bits:{53u,64u,65u,127u,128u,255u,256u,257u,
                        320u,511u,512u,1024u,2048u,4096u,4097u}) {
        const BigUnsigned leading=shift_left(BigUnsigned(1u),bits-1u);
        for (unsigned tail:{0u,1u,63u}) {
            const Integer magnitude=Integer::from_magnitude(
                add(leading,BigUnsigned(tail)),false);
            for (unsigned signs=0u;signs<4u;++signs) {
                const Integer a=(signs&1u)!=0u ? -magnitude : magnitude;
                const Integer b=(signs&2u)!=0u ? -magnitude : magnitude;
                for (unsigned gap:{0u,1u,7u,31u,63u,64u,65u}) {
                    const unsigned ea=bits/5u+1u,eb=ea+gap;
                    const NewDyadic x(a,ea),y(b,eb);
                    const Dyadic old_x(a,ea),old_y(b,eb);
                    const Rational expected_sum=Rational(old_x+old_y);
                    const Rational expected_difference=Rational(old_x-old_y);
                    check((x+y).to_rational()==expected_sum,
                          "R35 sum versus legacy Dyadic");
                    check((x-y).to_rational()==expected_difference,
                          "R35 difference versus legacy Dyadic");
                    check((y-x).to_rational()==Rational(old_y-old_x),
                          "R35 reversed subtraction");
                    NewDyadic inplace(x);
                    inplace+=y;
                    check(inplace.to_rational()==expected_sum,
                          "R35 in-place addition");
                    inplace=x;
                    inplace-=y;
                    check(inplace.to_rational()==expected_difference,
                          "R35 in-place subtraction");
                }
            }
        }
    }
    // R36: exercise the 16/17-limb x86-64 crossover and
    // 64/65-limb ARM64 crossover with the same exact Dyadic oracle.
    // Both branches must be semantically interchangeable.
    for (unsigned bits:{959u,960u,961u,1023u,1024u,1025u,
                        1087u,1088u,1089u,
                        4031u,4032u,4033u,
                        4095u,4096u,4097u,4159u,4160u,4161u}) {
        const BigUnsigned leading=shift_left(BigUnsigned(1u),bits-1u);
        const BigUnsigned right=add(leading,BigUnsigned(7u));
        for (unsigned signs=0u;signs<4u;++signs) {
            const Integer lhs=Integer::from_magnitude(
                add(leading,BigUnsigned(3u)),(signs&1u)!=0u);
            const Integer rhs=Integer::from_magnitude(
                right,(signs&2u)!=0u);
            for (unsigned gap:{0u,1u,63u,64u}) {
                const unsigned exp=bits/5u+1u;
                const NewDyadic x(lhs,exp),y(rhs,exp+gap);
                const Dyadic ref_x(lhs,exp),ref_y(rhs,exp+gap);
                check((x+y).to_rational()==Rational(ref_x+ref_y),
                      "R36 arch dispatch boundary addition");
                check((x-y).to_rational()==Rational(ref_x-ref_y),
                      "R36 arch dispatch boundary subtraction");
                check((y-x).to_rational()==Rational(ref_y-ref_x),
                      "R36 arch dispatch boundary reverse subtraction");
            }
        }
    }
    for (unsigned gap:{0u,1u,63u,64u,65u}) {
        const NewDyadic zero(0),one(Integer(1),gap);
        check(zero+one==one,"R35 zero plus nonzero");
        check(zero-one==-one,"R35 zero minus nonzero");
        check(one-zero==one,"R35 nonzero minus zero");
    }
}

void test_binary64_and_scaled_comparison() {
    const auto check_double=[](Integer const& coefficient,unsigned exponent) {
        const Dyadic old_value(coefficient,exponent);
        const NewDyadic new_value(coefficient,exponent);
        if (std::bit_cast<std::uint64_t>(old_value.get_d())!=
            std::bit_cast<std::uint64_t>(new_value.get_d())) {
            throw std::runtime_error("legacy binary64 mismatch at exponent "+
                                     std::to_string(exponent));
        }
    };
    for (unsigned width:{1u,2u,3u,52u,53u,54u,63u,64u,65u,
                         127u,128u,129u,255u,256u,257u,
                         1023u,1024u,1025u,4095u,4096u,4097u}) {
        const BigUnsigned leading=shift_left(BigUnsigned(1u),width-1u);
        for (unsigned tail:{0u,1u,3u,7u,127u}) {
            const BigUnsigned magnitude=add(leading,BigUnsigned(tail));
            const Integer pos=Integer::from_magnitude(magnitude,false);
            const Integer neg=Integer::from_magnitude(magnitude,true);
            for (unsigned exponent:{0u,1u,2u,width/2u,width,width+1u,
                                    width+1021u,width+1022u,
                                    width+1072u,width+1073u,width+1074u}) {
                check_double(pos,exponent);
                check_double(neg,exponent);
            }
        }
    }
    for (unsigned width:{1u,53u,64u,65u,127u,128u,129u,
                         255u,256u,257u,1024u,2048u,4096u}) {
        const BigUnsigned leading=shift_left(BigUnsigned(1u),width-1u);
        for (unsigned i=0u;i<24u;++i) {
            const Integer a=Integer::from_magnitude(
                add(leading,BigUnsigned(1u+2u*i)),i%4u>=2u);
            const Integer b=Integer::from_magnitude(
                add(leading,BigUnsigned(7u+5u*i)),
                i%4u==1u || i%4u==2u);
            for (unsigned ea:{0u,1u,width/2u,width,width+1u}) {
                for (unsigned eb:{0u,1u,width/2u,width,width+1u}) {
                    const NewDyadic x(a,ea),y(b,eb);
                    const Dyadic old_x(a,ea),old_y(b,eb);
                    if (cmp(x,y)!=cmp(old_x,old_y)) {
                        throw std::runtime_error("legacy comparison mismatch");
                    }
                }
            }
        }
    }
}

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

// D42: exact absolute value of large canonical signed coefficients.
// Preserve the input, denominator exponent and special-value behavior.
void test_wide_abs() {
    for (unsigned width:{1u,53u,128u,256u,1024u,2048u,4096u,8192u}) {
        const BigUnsigned high=shift_left(BigUnsigned(1u),width-1u);
        for (unsigned tail:{1u,3u,15u}) {
            const BigUnsigned magnitude=add(high,BigUnsigned(tail));
            for (unsigned exponent:{0u,1u,7u,63u,width/5u+1u}) {
                for (unsigned sign=0u;sign<2u;++sign) {
                    const Integer coefficient=Integer::from_magnitude(
                        magnitude,sign!=0u);
                    const NewDyadic value(coefficient,exponent);
                    const NewDyadic original(value);
                    const Dyadic legacy(coefficient,exponent);
                    const Positive<NewDyadic> actual=abs(value);
                    same_value(actual,abs(legacy),"D42 legacy abs");
                    check(actual.denominator_exponent()==
                          value.denominator_exponent(),"D42 abs exponent");
                    check(sgn(actual)!=Sign::NEGATIVE,"D42 abs sign");
                    check(value==original,"D42 abs input preserved");
                }
            }
        }
    }
    check(abs(NewDyadic(0))==NewDyadic(0),"D42 zero abs");
    check(abs(NewDyadic::inf(Sign::NEGATIVE))==
          NewDyadic::inf(Sign::POSITIVE),"D42 negative infinity abs");
    check(abs(NewDyadic::inf(Sign::POSITIVE))==
          NewDyadic::inf(Sign::POSITIVE),"D42 positive infinity abs");
    check(is_nan(abs(NewDyadic::nan())),"D42 NaN abs");
}

// D47: native FloatAP conversion, preserving exactness.
void test_floatap_native_bridge() {
    const ArbitraryPrecision precision(256_bits);
    for (unsigned width:{1u,53u,127u,128u,191u}) {
        const BigUnsigned magnitude=add(shift_left(BigUnsigned(1u),width-1u),BigUnsigned(13u));
        for (unsigned exponent:{0u,1u,7u,63u,127u}) {
            for (bool negative:{false,true}) {
                const Integer coefficient=Integer::from_magnitude(magnitude,negative);
                const NewDyadic source(coefficient,exponent);
                const FloatAP value(source,precision);
                const NewDyadic restored=static_cast<NewDyadic>(value);
                check(restored.to_rational()==source.to_rational(),"D47 exact roundtrip");
                check(restored.to_rational()==static_cast<Rational>(value),"D47 Rational oracle");
                FloatAP assigned(precision);
                assigned=source;
                check(static_cast<NewDyadic>(assigned).to_rational()==source.to_rational(),
                      "D47 exact assignment");
            }
        }
    }
    const BigUnsigned wide=add(shift_left(BigUnsigned(1u),256u),BigUnsigned(1u));
    const NewDyadic inexact(Integer::from_magnitude(wide,false));
    must_fail([&] { (void)FloatAP(inexact,precision); },"D47 reject inexact");
    const FloatAP rounded(inexact,FloatAP::ROUND_TO_NEAREST,precision);
    check(is_finite(static_cast<NewDyadic>(rounded)),"D47 explicit rounding");
    check(is_nan(static_cast<NewDyadic>(FloatAP::nan(precision))),"D47 NaN");
    check(is_inf(static_cast<NewDyadic>(FloatAP::inf(Sign::NEGATIVE,precision))),"D47 infinity");
}

// D48: native exact bounds and ball arithmetic, without GMP storage.
void test_native_wrappers() {
    const NewDyadicBounds positive(NewDyadic(1),NewDyadic(3));
    const NewDyadicBounds negative(NewDyadic(-3),NewDyadic(-1));
    const NewDyadicBounds crossing(NewDyadic(-2),NewDyadic(3));
    check((positive+NewDyadicBounds(NewDyadic(2),NewDyadic(4))).upper_raw()
          ==NewDyadic(7),"D48 bounds add");
    check((positive-negative).lower_raw()==NewDyadic(2),"D48 bounds sub");
    check((positive*negative).lower_raw()==NewDyadic(-9),"D48 bounds mul lower");
    check((positive*negative).upper_raw()==NewDyadic(-1),"D48 bounds mul upper");
    check(sqr(crossing).lower_raw()==NewDyadic(0),"D48 square crossing lower");
    check(sqr(crossing).upper_raw()==NewDyadic(9),"D48 square crossing upper");
    check(pow(crossing,Nat(3u)).lower_raw()==NewDyadic(-8),"D48 odd power");
    check(pow(crossing,Nat(2u)).lower_raw()==NewDyadic(0),"D48 even power");
    check(abs(crossing).upper_raw()==NewDyadic(3),"D48 bounds abs");
    check(refines(NewDyadicBounds(NewDyadic(2),NewDyadic(2)),positive),
          "D48 bounds refinement");
    const NewDyadicBall ball(NewDyadic(3),NewDyadic(1));
    const NewDyadicBounds enclosure=static_cast<NewDyadicBounds>(ball);
    check(enclosure.lower_raw()==NewDyadic(2) &&
          enclosure.upper_raw()==NewDyadic(4),"D48 ball enclosure");
    check((ball*NewDyadicBall(NewDyadic(2),NewDyadic(0))).value_raw()
          ==NewDyadic(6),"D48 ball multiplication");
    const NewDyadicLowerBound lower(positive);
    const NewDyadicUpperBound upper(positive);
    const NewDyadicApproximation approximation(NewDyadic(2));
    check(lower.raw()==NewDyadic(1) && upper.raw()==NewDyadic(3)
          && approximation.raw()==NewDyadic(2),"D48 wrapper values");
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
        test_binary64_and_scaled_comparison();
        test_aligned_add_sub();
        test_owned_shifted_arithmetic();
        test_canonical_fraction_rounding();
        test_wide_abs();
        test_floatap_native_bridge();
        test_native_wrappers();
        test_infinity();
        std::cout<<"NewDyadic test port (literal/conversions/arithmetic/"
                    "comparisons/infinity) PASS\n";
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"NewDyadic test port FAILED: "<<error.what()<<'\n';
        return 1;
    }
}
