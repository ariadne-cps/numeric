/***************************************************************************
 *            test_float.cpp
 *
 *  Copyright  2006-20  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

#include "numeric/builtin.hpp"
#include "numeric/double.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/number.hpp"
#include "numeric/upper_number.hpp"
#include "numeric/lower_number.hpp"
#include "numeric/float.decl.hpp"


#include "numeric/float_ball.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_error.hpp"
#include "numeric/float_literals.hpp"
#include "numeric/positive.hpp"
#include "numeric/rounded_float.hpp"
#include "numeric/casts.hpp"

#include "numeric/floatmp.hpp"
#include "utility/uniform_array.hpp"
#include "macro/testing.hpp"
#include "test_floats.hpp"

using namespace Ariadne;
using namespace std;

#ifndef ARIADNE_TEST_DATA_DIR
#error "Provide the absolute test data directory from CMake"
#endif

namespace {

// These exact directed-rounding fixtures are generated offline using signed
// integer division only. The test must never use another float implementation
// as an oracle for FloatMP.
Void test_floatmp_directed_bound_oracles() {
    constexpr char filename[]=ARIADNE_TEST_DATA_DIR "/floatmp_directed_bounds.tsv";
    std::ifstream data(filename);
    if (!data) throw std::runtime_error(std::string("Missing FloatMP oracle: ")+filename);
    std::set<std::string> unique_keys;
    std::string line;
    Nat rows=0u;
    while (std::getline(data,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        std::size_t bits=0;
        std::string n,d,lo_n,lo_d,hi_n,hi_d,extra;
        if (!(fields>>bits>>n>>d>>lo_n>>lo_d>>hi_n>>hi_d) || (fields>>extra))
            throw std::runtime_error("Malformed FloatMP directed-bound oracle row");
        if (!unique_keys.emplace(std::to_string(bits)+":"+n+"/"+d).second)
            throw std::runtime_error("Duplicate FloatMP directed-bound oracle input");
        const MultiplePrecision prec(bits);
        const Rational input{Integer(String(n)),Integer(String(d))};
        const Rational lower_expected{Integer(String(lo_n)),Integer(String(lo_d))};
        const Rational upper_expected{Integer(String(hi_n)),Integer(String(hi_d))};
        ARIADNE_TEST_ASSERT(lower_expected<=input && input<=upper_expected);
        const LowerBound<FloatMP> direct_lower(input,prec);
        const UpperBound<FloatMP> direct_upper(input,prec);
        const ValidatedLowerNumber exact_generic_lower(input);
        const ValidatedUpperNumber exact_generic_upper(input);
        const LowerBound<FloatMP> generic_lower(exact_generic_lower,prec);
        const UpperBound<FloatMP> generic_upper(exact_generic_upper,prec);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(direct_lower.raw()),lower_expected);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(direct_upper.raw()),upper_expected);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(generic_lower.raw()),lower_expected);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(generic_upper.raw()),upper_expected);
        ARIADNE_TEST_EQUALS(direct_lower.precision(),prec);
        ARIADNE_TEST_EQUALS(direct_upper.precision(),prec);
        ++rows;
    }
    if (rows!=90u) throw std::runtime_error("Unexpected FloatMP directed-bound oracle row count");
}


template<class PR>
constexpr bool check_float_concept()
{
    using PRE=DoublePrecision;
    using RawFloatType=RawFloat<PR>;
    using FloatType=Float<PR>;
    using FloatBoundsType=FloatBounds<PR>;
    return requires(PR pr, PRE pre, Boolean b, Nat m, Int n, Integer z, Dyadic w, ExactDouble d, TwoExp t,
                    RawFloatType f, FloatType vx, FloatType rx, FloatBoundsType rbx, FloatBall<PR,PRE> rmx) {
        FloatType(m,pr); FloatType(n,pr); FloatType(z,pr); FloatType(w,pr); FloatType(d,pr); FloatType(t,pr); FloatType(pr); FloatType(f);
        rx=m; rx=n; rx=z; rx=w; rx=d; rx=t;
        rx=operator+(vx); rx=operator-(vx); rx=t*vx; rx=vx*t; rx=vx/t;
        rbx=operator+(vx,vx); rbx=operator-(vx,vx); rbx=operator*(vx,vx); rbx=operator/(vx,vx);
        rx=nul(vx); rx=pos(vx); rx=neg(vx); rx=hlf(vx); rx=mul(vx,t); rx=div(vx,t);
        rmx=add(vx,vx,pre); rmx=sub(vx,vx,pre); rmx=mul(vx,vx,pre); rmx=div(vx,vx,pre);
        rbx=add(vx,vx); rbx=sub(vx,vx); rbx=mul(vx,vx); rbx=div(vx,vx);
        rbx=sqr(vx); rbx=rec(vx); rbx=pow(vx,m); rbx=pow(vx,n);
        rx=max(vx,vx); rx=min(vx,vx); rx=abs(vx);
        b=(vx==vx); b=(vx!=vx); b=(vx<=vx); b=(vx>=vx); b=(vx<vx); b=(vx>vx);
    };
}

static_assert(check_float_concept<DoublePrecision>());
static_assert(check_float_concept<MultiplePrecision>());


} // namespace


Void test_float_literals()
{
    FloatDP half(ExactDouble(0.5),dp);

    auto error=0.5_error;
    auto exact=0.5_exact;
    auto near_value=0.5_near;
    auto upper=0.5_upper;
    auto lower=0.5_lower;
    auto approximation=0.5_approx;

    ARIADNE_TEST_EQUALS(error.raw(),half);
    ARIADNE_TEST_EQUALS(exact,half);
    ARIADNE_TEST_EQUALS(near_value.value_raw(),half);
    ARIADNE_TEST_EQUALS(near_value.error_raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(upper.raw(),half);
    ARIADNE_TEST_EQUALS(lower.raw(),half);
    ARIADNE_TEST_EQUALS(approximation.raw(),half);

    DecimalPlaces places(3u);
    DecimalPrecision precision(5u);
    ARIADNE_TEST_EQUALS(static_cast<unsigned int>(places),3u);
    ARIADNE_TEST_EQUALS(static_cast<unsigned int>(precision),5u);

    if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
        ARIADNE_TEST_FAIL(operator""_error(0.1L));
        ARIADNE_TEST_FAIL(operator""_exact(0.1L));

        long double near_input=0.1L;
        auto near_inexact=operator""_near(near_input);
        long double near_value_raw=static_cast<long double>(near_inexact.value_raw().get_d());
        long double near_error_raw=static_cast<long double>(near_inexact.error_raw().get_d());
        ARIADNE_TEST_ASSERT(std::abs(near_value_raw-near_input)<=near_error_raw);

        long double upper_input=0.3L;
        auto upper_inexact=operator""_upper(upper_input);
        ARIADNE_TEST_ASSERT(static_cast<long double>(upper_inexact.raw().get_d())>=upper_input);

        long double lower_input=0.1L;
        auto lower_inexact=operator""_lower(lower_input);
        ARIADNE_TEST_ASSERT(static_cast<long double>(lower_inexact.raw().get_d())<=lower_input);
    }
}


template<class PR>
class TestFloat
{
    using PRE=DoublePrecision;
    typedef RawFloat<PR> RawFloatType;
    typedef FloatBounds<PR> FloatBoundsType;
    typedef FloatBall<PR,PRE> FloatBallType;
    typedef Float<PR> FloatType;
  private:
    PR precision;
  public:
    TestFloat(PR prec) : precision(prec) { }
    Void test();
  private:
    Void test_conversions();
    Void test_operations();
    Void test_predicates();
};

template<class PR> Void
TestFloat<PR>::test()
{
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_operations());
    ARIADNE_TEST_CALL(test_predicates());
}

template<class PR> Void
TestFloat<PR>::test_conversions()
{
    ARIADNE_TEST_EQUALS(cast_integer(Float<PR>(Dyadic(2),precision)),Integer(2));
    ARIADNE_TEST_FAIL(cast_integer(Float<PR>(Dyadic(7,2u),precision)));

    Positive<Float<PR>> positive_value(2u,precision);
    ARIADNE_TEST_EQUALS(cast_unsigned(positive_value),Float<PR>(2u,precision));

    FloatError<PR> error(1u,precision);
    ARIADNE_TEST_EQUALS(cast_exact(error),error.raw());

    FloatApproximation<PR> approximation(2u,precision);
    ARIADNE_TEST_EQUALS(cast_exact(approximation),approximation.raw());

    if constexpr (Same<PR,DoublePrecision>) {
        RawFloatDP raw(2u,dp);
        FloatDPApproximation approximate(2u,dp);
        ARIADNE_TEST_EQUALS(cast_raw(raw),raw);
        ARIADNE_TEST_EQUALS(cast_raw(approximate),approximate.raw());
        ARIADNE_TEST_EQUALS(cast_approximate(raw).raw(),raw);
        ARIADNE_TEST_EQUALS(cast_approximate(approximate).raw(),approximate.raw());

        ARIADNE_TEST_EXECUTE((RoundedFloatDP{approximate}));
        ARIADNE_TEST_EQUALS(FloatDP(String("1.25"),dp),Dyadic(5,2u));
        ARIADNE_TEST_EQUALS(FloatDP(TwoExp(10),dp),Dyadic(TwoExp(10)));
        ARIADNE_TEST_FAIL((FloatDP{TwoExp(1024),dp}));

        Dyadic too_precise{Integer(String("9007199254740993"))};
        ARIADNE_TEST_FAIL((FloatDP{too_precise,dp}));

        FloatDP exact_assignment(0,dp);
        ARIADNE_TEST_EQUALS((exact_assignment=Dyadic(2)),Dyadic(2));
        ARIADNE_TEST_FAIL(exact_assignment=too_precise);
        ARIADNE_TEST_EXECUTE(exact_assignment=Dyadic::nan());
        ARIADNE_TEST_ASSERT(is_nan(exact_assignment));
        ARIADNE_TEST_EXECUTE(exact_assignment=Dyadic::inf(Sign::POSITIVE));
        ARIADNE_TEST_ASSERT(is_inf(exact_assignment));
        ARIADNE_TEST_EXECUTE(exact_assignment=Dyadic::inf(Sign::NEGATIVE));
        ARIADNE_TEST_ASSERT(is_inf(exact_assignment));

        ARIADNE_TEST_EQUALS(integer_cast<Nat>(FloatDP(3u,dp)),Nat(3u));
        ARIADNE_TEST_EQUALS(integer_cast<Int>(FloatDP(-3,dp)),Int(-3));

        FloatDP resized(1,dp);
        resized.set_precision(dp);
        ARIADNE_TEST_EQUALS(resized.precision(),dp);

        Float32 f32(FloatDP(1.25_x,dp),FloatDP::ROUND_TO_NEAREST);
        ARIADNE_TEST_EQUALS(static_cast<FloatDP>(f32),FloatDP(1.25_x,dp));

        ARIADNE_TEST_EQUALS(class_name<double>(),String("double"));
        ARIADNE_TEST_EQUALS(class_name<ApproximateDouble>(),String("ApproximateDouble"));
        ARIADNE_TEST_EQUALS(class_name<ExactDouble>(),String("ExactDouble"));
        ARIADNE_TEST_EQUALS(class_name<Rounded<FloatDP>>(),String("Rounded<FloatDP>"));
    } else if constexpr (Same<PR,MultiplePrecision>) {
        PR pr=precision;

        ARIADNE_TEST_ASSERT(PR(64_bits)<=PR(128_bits));

        FloatType from_string(String("1.25"),pr);
        ARIADNE_TEST_EQUALS(from_string,Dyadic(5,2u));

        const std::uint64_t max_unsigned=std::numeric_limits<std::uint64_t>::max();
        const FloatType from_max_unsigned(max_unsigned,pr);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(from_max_unsigned),Rational(Integer(max_unsigned)));

        FloatType parsed(7,pr);
        std::istringstream empty_input("");
        ARIADNE_TEST_EXECUTE(empty_input >> parsed);
        ARIADNE_TEST_ASSERT(empty_input.fail());
        ARIADNE_TEST_EQUALS(parsed,Dyadic(7));

        std::istringstream invalid_input("2e+");
        ARIADNE_TEST_EXECUTE(invalid_input >> parsed);
        ARIADNE_TEST_ASSERT(invalid_input.fail());
        ARIADNE_TEST_EQUALS(parsed,Dyadic(7));

        Dyadic too_wide{Integer(String("18446744073709551617"))};
        ARIADNE_TEST_FAIL(FloatType(too_wide,PR(64_bits)));
        ARIADNE_TEST_ASSERT(is_nan(FloatType(Dyadic::nan(),pr)));

        FloatType exact_assignment(0,PR(64_bits));
        ARIADNE_TEST_EQUALS((exact_assignment=Dyadic(2)),Dyadic(2));
        ARIADNE_TEST_FAIL(exact_assignment=too_wide);

        // Regression: UniformArray exposes CheckedAssignable<FloatMP>.
        // Its precision-preserving T assignment must hide FloatMP's base
        // copy assignment rather than compete with a by-value overload.
        UniformArray<FloatMP> checked_values(1u,pr);
        checked_values[0u]=FloatMP(3,pr);
        ARIADNE_TEST_EQUALS(checked_values[0u],FloatMP(3,pr));
        ARIADNE_TEST_EQUALS(checked_values[0u].characteristics(),pr);

        FloatType zero(pr);
        ARIADNE_TEST_ASSERT(is_zero(zero));
        ARIADNE_TEST_ASSERT(is_inf(FloatType::inf(pr)));
        ARIADNE_TEST_ASSERT(is_nan(FloatType::nan(pr)));

        FloatType resized(Dyadic(3,2u),pr);
        resized.set_precision(PR(96_bits));
        ARIADNE_TEST_EQUALS(resized.precision(),PR(96_bits));

        FloatType scaled(3,pr);
        scaled*=TwoExp(2);
        ARIADNE_TEST_EQUALS(scaled,Dyadic(12));
        scaled/=TwoExp(1);
        ARIADNE_TEST_EQUALS(scaled,Dyadic(6));
        ARIADNE_TEST_ASSERT(same(scaled,FloatType(6,pr)));
        ARIADNE_TEST_EQUALS(integer_cast<Nat>(FloatType(7u,pr)),Nat(7u));

        if constexpr (Same<PR,MultiplePrecision>) {
            // The FloatDP conversion must preserve IEEE exceptional values,
            // including the sign of zero, without converting them to Rational.
            ARIADNE_TEST_ASSERT(std::isnan(FloatDP(FloatType::nan(pr),
                FloatDP::ROUND_TO_NEAREST,dp).get_d()));

            const FloatDP positive_infinity(FloatType::inf(Sign::POSITIVE,pr),
                FloatDP::ROUND_TO_NEAREST,dp);
            const FloatDP negative_infinity(FloatType::inf(Sign::NEGATIVE,pr),
                FloatDP::ROUND_TO_NEAREST,dp);
            ARIADNE_TEST_ASSERT(std::isinf(positive_infinity.get_d()));
            ARIADNE_TEST_ASSERT(not std::signbit(positive_infinity.get_d()));
            ARIADNE_TEST_ASSERT(std::isinf(negative_infinity.get_d()));
            ARIADNE_TEST_ASSERT(std::signbit(negative_infinity.get_d()));

            const FloatType negative_zero=FloatType::zero(
                static_cast<unsigned long>(pr.bits()),true);
            const FloatType positive_zero=FloatType::zero(
                static_cast<unsigned long>(pr.bits()));
            const FloatDP converted_negative_zero(negative_zero,FloatDP::ROUND_TO_NEAREST,dp);
            const FloatDP converted_positive_zero(positive_zero,FloatDP::ROUND_TO_NEAREST,dp);
            ARIADNE_TEST_ASSERT(converted_negative_zero.get_d()==0.0);
            ARIADNE_TEST_ASSERT(std::signbit(converted_negative_zero.get_d()));
            ARIADNE_TEST_ASSERT(converted_positive_zero.get_d()==0.0);
            ARIADNE_TEST_ASSERT(not std::signbit(converted_positive_zero.get_d()));

            const FloatType third(Rational(1,3),FloatType::ROUND_TO_NEAREST,pr);
            ARIADNE_TEST_EQUALS(FloatDP(third,FloatDP::ROUND_DOWNWARD,dp),
                FloatDP(static_cast<Rational>(third),FloatDP::ROUND_DOWNWARD,dp));
            ARIADNE_TEST_EQUALS(FloatDP(third,FloatDP::ROUND_UPWARD,dp),
                FloatDP(static_cast<Rational>(third),FloatDP::ROUND_UPWARD,dp));
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Preserve the public zero-interval representation,
            // rather than the old generic formatter's "{0.:0.}" spelling.
            std::ostringstream ap_zero_stream;
            ap_zero_stream << Bounds<FloatMP>(0,pr);
            ARIADNE_TEST_EQUALS(ap_zero_stream.str(),String("0.0[:]"));

            // Verify sign-sensitive interval arithmetic in the FloatMP
            // specialization, including all endpoint sign combinations.
            using B=Bounds<FloatMP>;
            const B positive(1,2,pr);
            const B negative(-4,-2,pr);
            const B crossing(-2,3,pr);
            const B singleton(2,pr);
            const B another_singleton(2,pr);
            const B other_singleton(3,pr);

            const B positive_square=sqr(positive);
            const B negative_square=sqr(negative);
            const B crossing_square=sqr(crossing);
            ARIADNE_TEST_EQUALS(positive_square.lower_raw(),FloatType(1,pr));
            ARIADNE_TEST_EQUALS(positive_square.upper_raw(),FloatType(4,pr));
            ARIADNE_TEST_EQUALS(negative_square.lower_raw(),FloatType(4,pr));
            ARIADNE_TEST_EQUALS(negative_square.upper_raw(),FloatType(16,pr));
            ARIADNE_TEST_EQUALS(crossing_square.lower_raw(),FloatType(0,pr));
            ARIADNE_TEST_EQUALS(crossing_square.upper_raw(),FloatType(9,pr));

            ARIADNE_TEST_ASSERT(definitely(eq(singleton,another_singleton)));
            ARIADNE_TEST_ASSERT(not possibly(eq(singleton,other_singleton)));
            ARIADNE_TEST_ASSERT(not definitely(eq(positive,crossing)));
            ARIADNE_TEST_ASSERT(possibly(eq(positive,crossing)));
            // Touching intervals are consistent but not identical singletons.
            // This exercises the false outcome of the second conjunct in
            // Bounds<F>::eq for the FloatMP instantiation.
            const B touching_left(1,2,pr);
            const B touching_right(0,1,pr);
            ARIADNE_TEST_ASSERT(possibly(eq(touching_left,touching_right)));
            ARIADNE_TEST_ASSERT(not definitely(eq(touching_left,touching_right)));
            ARIADNE_TEST_ASSERT(not possibly(eq(negative,positive)));
            ARIADNE_TEST_ASSERT(not same(positive,negative));
            ARIADNE_TEST_ASSERT(not same(positive,B(1,3,pr)));

            const B products[]={
                mul(positive,positive),mul(positive,negative),mul(positive,crossing),
                mul(negative,positive),mul(negative,negative),mul(negative,crossing),
                mul(crossing,positive),mul(crossing,negative),mul(crossing,crossing)
            };
            const Int lower[]={1,-8,-4,-8,4,-12,-4,-12,-6};
            const Int upper[]={4,-2,6,-2,16,8,6,8,9};
            for (Nat i=0u; i<9u; ++i) {
                ARIADNE_TEST_EQUALS(products[i].lower_raw(),FloatType(lower[i],pr));
                ARIADNE_TEST_EQUALS(products[i].upper_raw(),FloatType(upper[i],pr));
            }

            const B division_negative=div(positive,negative);
            ARIADNE_TEST_ASSERT(division_negative.lower_raw()<=FloatType(Dyadic(-1,1u),pr));
            ARIADNE_TEST_ASSERT(division_negative.upper_raw()>=FloatType(Dyadic(-1,2u),pr));
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Verify the concrete FloatMP wrapper specializations, including
            // rigorous enclosure when converting from generic real numbers.
            ARIADNE_TEST_EQUALS(class_name<Bounds<FloatMP>>(),String("FloatMPBounds"));
            ARIADNE_TEST_EQUALS(class_name<LowerBound<FloatMP>>(),String("FloatMPLowerBound"));
            ARIADNE_TEST_EQUALS(class_name<UpperBound<FloatMP>>(),String("FloatMPUpperBound"));

            const Rational exact_third(1,3);
            const Decimal exact_decimal(Integer(125),3u);
            const FloatMP decimal_nearest(exact_decimal,FloatMP::ROUND_TO_NEAREST,pr);
            ARIADNE_TEST_EQUALS(static_cast<Rational>(decimal_nearest),Rational(1,8));
            const Bounds<FloatMP> real_enclosure(Real(exact_third),pr);
            ARIADNE_TEST_ASSERT(cmp(real_enclosure.lower_raw(),exact_third)!=Comparison::GREATER);
            ARIADNE_TEST_ASSERT(cmp(real_enclosure.upper_raw(),exact_third)!=Comparison::LESS);

            // The direct and polymorphic FloatMP conversion paths are tested
            // against a frozen, independently computed Rational TSV oracle.
            // No external floating-point backend is used as an oracle.
        }

        // Exercise native precision controls and decimal formatting.
        // The independent exact arithmetic raw-pointer path has no counterpart in FloatMP.
        const MultiplePrecision from_bits(53_bits);
        ARIADNE_TEST_EQUALS(from_bits.bits(),53u);

        std::ostringstream precision_repr;
        repr(precision_repr,MultiplePrecision(96_bits));
        ARIADNE_TEST_EQUALS(precision_repr.str(),String("MultiplePrecision(96)"));

        FloatType source(Dyadic(3,2u),pr);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(source),Rational(3,4));

        const MultiplePrecision old_default=FloatType::get_default_precision();
        FloatType::set_default_precision(MultiplePrecision(96_bits));
        ARIADNE_TEST_EQUALS(FloatType::get_default_precision(),MultiplePrecision(96_bits));
        FloatType::set_default_precision(old_default);

        resized=Integer(2);
        ARIADNE_TEST_EQUALS(resized,Dyadic(2));
    }
}

template<class PR> Void
TestFloat<PR>::test_operations()
{
    PR pr=precision;
    DoublePrecision pre;
    FloatType vr(pr);

    Nat m(5);
    Int n(-3);
    TwoExp t(-12);
    Dyadic w(-3,1u);
    FloatType vx(w,pr);
    Dyadic w1(-3,1u), w2(5,2u);
    FloatType vx1(w1,pr), vx2(w2,pr);

    ARIADNE_TEST_EQUALS(vx,w);

    FloatType mixed_float_one(1,pr);
    Real mixed_real_two(2);

    ARIADNE_TEST_ASSERT(not possibly(mixed_real_two==mixed_float_one));
    ARIADNE_TEST_ASSERT(definitely(mixed_real_two!=mixed_float_one));
    ARIADNE_TEST_ASSERT(not possibly(mixed_real_two< mixed_float_one));
    ARIADNE_TEST_ASSERT(definitely(mixed_real_two> mixed_float_one));
    ARIADNE_TEST_ASSERT(not possibly(mixed_real_two<=mixed_float_one));
    ARIADNE_TEST_ASSERT(definitely(mixed_real_two>=mixed_float_one));

    ARIADNE_TEST_ASSERT(not possibly(mixed_float_one==mixed_real_two));
    ARIADNE_TEST_ASSERT(definitely(mixed_float_one!=mixed_real_two));
    ARIADNE_TEST_ASSERT(definitely(mixed_float_one< mixed_real_two));
    ARIADNE_TEST_ASSERT(not possibly(mixed_float_one> mixed_real_two));
    ARIADNE_TEST_ASSERT(definitely(mixed_float_one<=mixed_real_two));
    ARIADNE_TEST_ASSERT(not possibly(mixed_float_one>=mixed_real_two));

    FloatBoundsType mixed_three(3,pr);
    FloatBoundsType mixed_minus_one(-1,pr);
    FloatBoundsType mixed_two(2,pr);
    FloatBoundsType mixed_one(1,pr);
    FloatBoundsType mixed_half(Dyadic(1,1u),pr);

    ARIADNE_TEST_SAME(mixed_float_one+mixed_real_two,mixed_three);
    ARIADNE_TEST_SAME(mixed_float_one-mixed_real_two,mixed_minus_one);
    ARIADNE_TEST_SAME(mixed_float_one*mixed_real_two,mixed_two);
    ARIADNE_TEST_SAME(mixed_float_one/mixed_real_two,mixed_half);
    ARIADNE_TEST_SAME(mixed_real_two+mixed_float_one,mixed_three);
    ARIADNE_TEST_SAME(mixed_real_two-mixed_float_one,mixed_one);
    ARIADNE_TEST_SAME(mixed_real_two*mixed_float_one,mixed_two);
    ARIADNE_TEST_SAME(mixed_real_two/mixed_float_one,mixed_two);

    ARIADNE_TEST_SAME(add(mixed_float_one,mixed_real_two),mixed_three);
    ARIADNE_TEST_SAME(sub(mixed_float_one,mixed_real_two),mixed_minus_one);
    ARIADNE_TEST_SAME(mul(mixed_float_one,mixed_real_two),mixed_two);
    ARIADNE_TEST_SAME(div(mixed_float_one,mixed_real_two),mixed_half);
    ARIADNE_TEST_SAME(add(mixed_real_two,mixed_float_one),mixed_three);
    ARIADNE_TEST_SAME(sub(mixed_real_two,mixed_float_one),mixed_one);
    ARIADNE_TEST_SAME(mul(mixed_real_two,mixed_float_one),mixed_two);
    ARIADNE_TEST_SAME(div(mixed_real_two,mixed_float_one),mixed_two);

    ARIADNE_TEST_SAME(max(mixed_float_one,mixed_real_two),mixed_two);
    ARIADNE_TEST_SAME(min(mixed_float_one,mixed_real_two),mixed_one);
    ARIADNE_TEST_SAME(max(mixed_real_two,mixed_float_one),mixed_two);
    ARIADNE_TEST_SAME(min(mixed_real_two,mixed_float_one),mixed_one);

    ARIADNE_TEST_EQUALS(FloatType(RawFloatType(1.25_x,pr)),1.25_dy);

    ARIADNE_TEST_EQUALS(FloatType(3u,pr),3.0_dy);
    ARIADNE_TEST_EQUALS(FloatType(-5,pr),-5.0_dy);
    ARIADNE_TEST_EQUALS(FloatType(ExactDouble(1.25),pr),1.25_dy);
    ARIADNE_TEST_EQUALS(FloatType(TwoExp(-3),pr),0.125_dy);
    ARIADNE_TEST_EQUALS(FloatType(Integer(-23),pr),-23.0_dy);
    ARIADNE_TEST_EQUALS(FloatType(Dyadic(-23,3u),pr),-2.875_dy);


    ARIADNE_TEST_EQUALS((vr=3u),3.0_dy);
    ARIADNE_TEST_EQUALS((vr=-5),-5.0_dy);
    ARIADNE_TEST_EQUALS((vr=Integer(-23)),-23.0_dy);
    ARIADNE_TEST_EQUALS((vr=TwoExp(-3)),0.125_dy);
    ARIADNE_TEST_EQUALS((vr=Dyadic(-23,3u)),-2.875_dy);

    ARIADNE_TEST_EQUALS(FloatType(w,pr).operator Dyadic(),w);
    ARIADNE_TEST_EQUALS(FloatType(w,pr).operator Rational(),w);

    ARIADNE_TEST_EQUALS(nul(vx),nul(w));
    ARIADNE_TEST_EQUALS(pos(vx),pos(w));
    ARIADNE_TEST_EQUALS(neg(vx),neg(w));
    ARIADNE_TEST_EQUALS(hlf(vx),hlf(w));

    ARIADNE_TEST_EQUALS(mul(vx,t),mul(w,t));
    ARIADNE_TEST_EQUALS(div(vx,t),w/t);

    ARIADNE_TEST_BINARY_PREDICATE(models,sqr(vx),sqr(w));
    ARIADNE_TEST_BINARY_PREDICATE(models,rec(vx),rec(w));

//    ARIADNE_TEST_BINARY_PREDICATE(models,sqr(vx,pre),sqr(w));
//    ARIADNE_TEST_BINARY_PREDICATE(models,rec(vx,pre),rec(w));

    ARIADNE_TEST_EQUALS(max(vx1,vx2),max(w1,w2));
    ARIADNE_TEST_EQUALS(min(vx1,vx2),min(w1,w2));
    ARIADNE_TEST_EQUALS(abs(vx),abs(w));
    ARIADNE_TEST_EQUALS(mig(vx).raw(),abs(w));
    ARIADNE_TEST_EQUALS(mag(vx).raw(),abs(w));
    ARIADNE_TEST_SAME(mag(vx),PositiveFloatUpperBound<PR>(abs(w),pr));
    ARIADNE_TEST_SAME(mig(vx),PositiveFloatLowerBound<PR>(abs(w),pr));

    ARIADNE_TEST_EQUALS(t*vx,t*w);
    ARIADNE_TEST_EQUALS(vx*t,w*t);
    ARIADNE_TEST_EQUALS(vx/t,w/t);

    ARIADNE_TEST_BINARY_PREDICATE(models,vx1+vx2,w1+w2);
    ARIADNE_TEST_BINARY_PREDICATE(models,vx1-vx2,w1-w2);
    ARIADNE_TEST_BINARY_PREDICATE(models,vx1*vx2,w1*w2);
    ARIADNE_TEST_BINARY_PREDICATE(models,vx1/vx2,w1/w2);

    ARIADNE_TEST_BINARY_PREDICATE(models,add(vx1,vx2),add(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,sub(vx1,vx2),sub(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,mul(vx1,vx2),mul(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,div(vx1,vx2),div(w1,w2));

    ARIADNE_TEST_BINARY_PREDICATE(models,add(vx1,vx2,pre),add(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,sub(vx1,vx2,pre),sub(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,mul(vx1,vx2,pre),mul(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,div(vx1,vx2,pre),div(w1,w2));


    ARIADNE_TEST_SAME(add(vx1,vx2,pre),FloatBallType(add(w1,w2),pr,pre));

    ARIADNE_TEST_BINARY_PREDICATE(models,pow(vx,m),pow(w,m));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(vx,n),pow(Rational(w),n));

//    friend Bounds<F> med(F const& x1, F const& x2);
//    friend Bounds<F> rad(F const& x1, F const& x2);

    ARIADNE_TEST_BINARY_PREDICATE(models,sqr(sqrt(abs(vx))),abs(w));
    ARIADNE_TEST_BINARY_PREDICATE(models,log(exp(vx)),w);
    ARIADNE_TEST_BINARY_PREDICATE(models,exp(log(abs(vx))),abs(w));
    ARIADNE_TEST_BINARY_PREDICATE(models,atan(tan(vx)),w);
    ARIADNE_TEST_BINARY_PREDICATE(models,tan(atan(vx)),w);

    ARIADNE_TEST_SAME(sin(vx),sin(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(cos(vx),cos(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(tan(vx),tan(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(tanh(vx),tanh(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(atan(vx),atan(FloatBounds<PR>(vx)));

    FloatBoundsType tanh_wide(-100,100,pr);
    FloatBoundsType tanh_positive(1,2,pr);
    FloatBoundsType tanh_negative(-2,-1,pr);
    auto wide_image=tanh(tanh_wide);
    auto positive_image=tanh(tanh_positive);
    auto negative_image=tanh(tanh_negative);
    ARIADNE_TEST_ASSERT(wide_image.lower_raw()>=FloatType(-1,pr));
    ARIADNE_TEST_ASSERT(wide_image.upper_raw()<=FloatType(1,pr));
    ARIADNE_TEST_ASSERT(wide_image.lower_raw()<FloatType(0,pr));
    ARIADNE_TEST_ASSERT(wide_image.upper_raw()>FloatType(0,pr));
    ARIADNE_TEST_ASSERT(positive_image.lower_raw()>FloatType(0,pr));
    ARIADNE_TEST_ASSERT(positive_image.upper_raw()<FloatType(1,pr));
    ARIADNE_TEST_ASSERT(negative_image.lower_raw()>FloatType(-1,pr));
    ARIADNE_TEST_ASSERT(negative_image.upper_raw()<FloatType(0,pr));

    FloatBallType tanh_ball(FloatBoundsType(-1,1,pr),pre);
    auto tanh_ball_image=tanh(tanh_ball);
    ARIADNE_TEST_ASSERT(tanh_ball_image.lower_raw()<=FloatType(0,pr));
    ARIADNE_TEST_ASSERT(tanh_ball_image.upper_raw()>=FloatType(0,pr));

    FloatLowerBound<PR> tanh_lower(-2,pr);
    FloatUpperBound<PR> tanh_upper(2,pr);
    ARIADNE_TEST_ASSERT(tanh(tanh_lower).raw()>FloatType(-1,pr));
    ARIADNE_TEST_ASSERT(tanh(tanh_upper).raw()<FloatType(1,pr));

//    ARIADNE_TEST_SAME(shft(vx,n),shft(w,n));

    if constexpr (Same<PR,DoublePrecision>) {
        FloatDP one(1,dp);
        FloatDP half(Dyadic(1,1u),dp);
        CurrentRoundingMode current;
        FloatDP::set_rounding_to_nearest();

        ARIADNE_TEST_EQUALS(pow_rnd(FloatDP(2,dp),Int(-2)),FloatDP(Dyadic(1,2u),dp));
        ARIADNE_TEST_EXECUTE(sqrt_rnd(one));
        ARIADNE_TEST_EXECUTE(exp_rnd(one));
        ARIADNE_TEST_EXECUTE(log_rnd(one));
        ARIADNE_TEST_EXECUTE(sin_rnd(one));
        ARIADNE_TEST_EXECUTE(cos_rnd(one));
        ARIADNE_TEST_EXECUTE(tan_rnd(one));
        ARIADNE_TEST_EXECUTE(atan_rnd(one));

        ARIADNE_TEST_EXECUTE(sqr(current,one));
        ARIADNE_TEST_EXECUTE(rec(current,one));
        ARIADNE_TEST_EXECUTE(add(current,one,half));
        ARIADNE_TEST_EXECUTE(sub(current,one,half));
        ARIADNE_TEST_EXECUTE(mul(current,one,half));
        ARIADNE_TEST_EXECUTE(div(current,one,half));
        ARIADNE_TEST_EXECUTE(fma(current,one,half,one));
        ARIADNE_TEST_EXECUTE(pow(current,one,Nat(2u)));
        ARIADNE_TEST_EXECUTE(pow(current,one,Int(-2)));
        ARIADNE_TEST_EXECUTE(sqrt(current,one));
        ARIADNE_TEST_EXECUTE(exp(current,one));
        ARIADNE_TEST_EXECUTE(log(current,one));
        ARIADNE_TEST_EXECUTE(sin(current,one));
        ARIADNE_TEST_EXECUTE(cos(current,one));
        ARIADNE_TEST_EXECUTE(tan(current,one));
        ARIADNE_TEST_EXECUTE(asin(current,half));
        ARIADNE_TEST_EXECUTE(acos(current,half));
        ARIADNE_TEST_EXECUTE(atan(current,one));
        ARIADNE_TEST_EXECUTE(FloatDP::pi(current,dp));

        ARIADNE_TEST_EQUALS(abs(FloatDP::ROUND_TO_NEAREST,FloatDP(-2,dp)),FloatDP(2,dp));
        ARIADNE_TEST_EQUALS(mag(FloatDP::ROUND_TO_NEAREST,FloatDP(-2,dp)),FloatDP(2,dp));
        ARIADNE_TEST_ASSERT(same(FloatDP(2,dp),FloatDP(2,dp)));
        ARIADNE_TEST_ASSERT(not same(FloatDP(2,dp),FloatDP(3,dp)));

        FloatDP scaled(2,dp);
        scaled*=TwoExp(2);
        ARIADNE_TEST_EQUALS(scaled,FloatDP(8,dp));
        scaled/=TwoExp(3);
        ARIADNE_TEST_EQUALS(scaled,FloatDP(1,dp));

        ARIADNE_TEST_EXECUTE(pos_opp(one));
        ARIADNE_TEST_EXECUTE(neg_opp(one));
        ARIADNE_TEST_EXECUTE(sqr_opp(one));
        ARIADNE_TEST_EXECUTE(rec_opp(one));
        ARIADNE_TEST_EXECUTE(add_opp(one,half));
        ARIADNE_TEST_EXECUTE(sub_opp(one,half));
        ARIADNE_TEST_EXECUTE(mul_opp(one,half));
        ARIADNE_TEST_EXECUTE(div_opp(one,half));

        Dyadic downward_target=Dyadic(1)+Dyadic(3,54u);
        FloatDP downward_value(downward_target,FloatDP::ROUND_DOWNWARD,dp);
        ARIADNE_TEST_ASSERT(Dyadic(downward_value)<=downward_target);

        Dyadic negative_downward_target=-Dyadic(1)-Dyadic(3,54u);
        FloatDP negative_downward_value(negative_downward_target,FloatDP::ROUND_DOWNWARD,dp);
        ARIADNE_TEST_ASSERT(Dyadic(negative_downward_value)<=negative_downward_target);

        ARIADNE_TEST_EXECUTE(FloatDP::pi(FloatDP::ROUND_UPWARD,dp));
        ARIADNE_TEST_EXECUTE(FloatDP::pi(FloatDP::ROUND_DOWNWARD,dp));
        ARIADNE_TEST_EXECUTE(FloatDP::pi(FloatDP::ROUND_TO_NEAREST,dp));
        ARIADNE_TEST_FAIL(FloatDP::pi(static_cast<FloatDP::RoundingModeType>(0xffffu),dp));

        ARIADNE_TEST_ASSERT(not one.literal().empty());
        ARIADNE_TEST_ASSERT(not one.literal(FloatDP::ROUND_UPWARD).empty());
        std::ostringstream dp_places_stream;
        write(dp_places_stream,one,DecimalPlaces(3u),FloatDP::ROUND_DOWNWARD);
        ARIADNE_TEST_ASSERT(not dp_places_stream.str().empty());
        std::ostringstream dp_repr_stream;
        repr(dp_repr_stream,one,FloatDP::ROUND_UPWARD);
        ARIADNE_TEST_ASSERT(not dp_repr_stream.str().empty());
        FloatDP::set_output_places(7u);
        ARIADNE_TEST_EQUALS(FloatDP::output_places,Nat(7u));
        ARIADNE_TEST_FAIL(one.literal(FloatDP::ROUND_TOWARD_ZERO));

        DoublePrecision p1;
        DoublePrecision p2;
        ARIADNE_TEST_EXECUTE(max(p1,p2));
        ARIADNE_TEST_ASSERT(p1<=p2);
        std::ostringstream precision_repr;
        repr(precision_repr,p1);
        ARIADNE_TEST_EQUALS(precision_repr.str(),std::string("DoublePrecision()"));

        FloatDP noinit(NoInit{});
        noinit.raw()=FloatDP(2,dp);
        ARIADNE_TEST_EQUALS(noinit,FloatDP(2,dp));
        ARIADNE_TEST_EXECUTE((FloatDP{one,dp}));
        ARIADNE_TEST_ASSERT(is_inf(FloatDP::inf(dp)));
        ARIADNE_TEST_ASSERT(is_finite(one));
        ARIADNE_TEST_ASSERT(is_zero(FloatDP(0,dp)));

        std::istringstream dp_input("+2.5:");
        FloatDP parsed(dp);
        dp_input >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatDP(Dyadic(5,1u),dp));
        ARIADNE_TEST_EQUALS(dp_input.peek(),static_cast<int>(':'));

        std::istringstream dp_negative("-.5");
        dp_negative >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatDP(Dyadic(-1,1u),dp));

        std::istringstream dp_exp("2E+3");
        dp_exp >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatDP(2000,dp));

        std::istringstream dp_exp_delimited("2e3:");
        dp_exp_delimited >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatDP(2000,dp));
        ARIADNE_TEST_EQUALS(dp_exp_delimited.peek(),static_cast<int>(':'));

        std::istringstream dp_exp_negative("2e-3");
        ARIADNE_TEST_EXECUTE(dp_exp_negative >> parsed);

        std::istringstream dp_overflow("1e9999");
        ARIADNE_TEST_EXECUTE(dp_overflow >> parsed);
        ARIADNE_TEST_ASSERT(dp_overflow.fail());

        std::istringstream dp_no_digits(".");
        ARIADNE_TEST_EXECUTE(dp_no_digits >> parsed);
        ARIADNE_TEST_ASSERT(dp_no_digits.fail());

        std::istringstream dp_bad_exp("2e+");
        ARIADNE_TEST_EXECUTE(dp_bad_exp >> parsed);
        ARIADNE_TEST_ASSERT(dp_bad_exp.fail());

        std::istringstream dp_prefailed("1");
        dp_prefailed.setstate(std::ios::failbit);
        ARIADNE_TEST_EXECUTE(dp_prefailed >> parsed);
    }

    if constexpr (Same<PR,MultiplePrecision>) {
        FloatType a(2,pr);
        FloatType b(1,pr);
        FloatDP dpa(1,dp);
        ExactDouble ed(0.5);
        CurrentRoundingMode current;

        FloatType::set_rounding_toward_zero();
        ARIADNE_TEST_EQUALS(fma(current,a,b,b),FloatType(3,pr));
        ARIADNE_TEST_EQUALS(pow(current,a,Nat(2u)),FloatType(4,pr));
        ARIADNE_TEST_EXECUTE(tan(current,b));
        ARIADNE_TEST_EXECUTE(asin(current,b));
        ARIADNE_TEST_EXECUTE(acos(current,b));
        ARIADNE_TEST_EXECUTE(FloatType::pi(current,pr));

        ARIADNE_TEST_EQUALS(add(current,a,dpa),FloatType(3,pr));
        ARIADNE_TEST_EQUALS(sub(current,a,dpa),FloatType(1,pr));
        ARIADNE_TEST_EQUALS(mul(current,a,dpa),FloatType(2,pr));
        ARIADNE_TEST_EQUALS(div(current,a,dpa),FloatType(2,pr));
        ARIADNE_TEST_EQUALS(add(current,dpa,a),FloatType(3,pr));
        ARIADNE_TEST_EQUALS(sub(current,dpa,a),FloatType(-1,pr));
        ARIADNE_TEST_EQUALS(mul(current,dpa,a),FloatType(2,pr));
        ARIADNE_TEST_EQUALS(div(current,dpa,a),FloatType(Dyadic(1,1u),pr));

        ARIADNE_TEST_EQUALS(add(current,a,ed),FloatType(Dyadic(5,1u),pr));
        ARIADNE_TEST_EQUALS(sub(current,a,ed),FloatType(Dyadic(3,1u),pr));
        ARIADNE_TEST_EQUALS(mul(current,a,ed),FloatType(1,pr));
        ARIADNE_TEST_EQUALS(div(current,a,ed),FloatType(4,pr));
        ARIADNE_TEST_EQUALS(add(current,ed,a),FloatType(Dyadic(5,1u),pr));
        ARIADNE_TEST_EQUALS(sub(current,ed,a),FloatType(Dyadic(-3,1u),pr));
        ARIADNE_TEST_EQUALS(mul(current,ed,a),FloatType(1,pr));
        ARIADNE_TEST_EQUALS(div(current,ed,a),FloatType(Dyadic(1,2u),pr));

        const auto nearest=FloatType::ROUND_TO_NEAREST;
        ARIADNE_TEST_EQUALS(abs(nearest,FloatType(-2,pr)),FloatType(2,pr));
        ARIADNE_TEST_EQUALS(mag(nearest,FloatType(-2,pr)),FloatType(2,pr));

        FloatType inplace(2,pr);
        ARIADNE_TEST_EQUALS(iadd(nearest,inplace,FloatType(1,pr)),FloatType(3,pr));
        ARIADNE_TEST_EQUALS(isub(nearest,inplace,FloatType(1,pr)),FloatType(2,pr));
        ARIADNE_TEST_EQUALS(imul(nearest,inplace,FloatType(3,pr)),FloatType(6,pr));
        ARIADNE_TEST_EQUALS(idiv(nearest,inplace,FloatType(2,pr)),FloatType(3,pr));

        ARIADNE_TEST_ASSERT(not a.literal().empty());
        ARIADNE_TEST_ASSERT(not a.literal(nearest).empty());
        std::ostringstream places_stream;
        write(places_stream,a,DecimalPlaces(3u),nearest);
        ARIADNE_TEST_ASSERT(not places_stream.str().empty());
        std::ostringstream repr_stream;
        repr(repr_stream,a,nearest);
        ARIADNE_TEST_ASSERT(not repr_stream.str().empty());

        ARIADNE_TEST_EQUALS(print(a,DecimalPlaces(0u),nearest),String("2."));
        ARIADNE_TEST_EQUALS(print(FloatType::nan(pr),DecimalPrecision(3u),nearest),String("nan"));

        std::istringstream float_input("\t\n+2.5:");
        FloatType parsed(pr);
        float_input >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatType(Dyadic(5,1u),pr));
        ARIADNE_TEST_EQUALS(float_input.peek(),static_cast<int>(':'));

        FloatType::set_output_places(7u);
        ARIADNE_TEST_EQUALS(FloatType::output_places,Nat(7u));
        ARIADNE_TEST_EQUALS(class_name<Rounded<FloatType>>(),String("Rounded<")+class_name<FloatType>()+String(">"));

        if constexpr (Same<PR,MultiplePrecision>) {
            ARIADNE_TEST_EXECUTE(FloatDP(a,FloatDP::ROUND_TO_NEAREST,dp));
            ARIADNE_TEST_EXECUTE(FloatDP(a,FloatDP::ROUND_TOWARD_ZERO,dp));
            const auto bad_mode=static_cast<FloatDP::RoundingModeType>(0xffffu);
            ARIADNE_TEST_FAIL(FloatDP(a,bad_mode,dp));
            ARIADNE_TEST_FAIL(FloatDP(FloatType(0,pr),bad_mode,dp));
            ARIADNE_TEST_FAIL(FloatDP(FloatType::nan(pr),bad_mode,dp));
            ARIADNE_TEST_FAIL(FloatDP(FloatType::inf(pr),bad_mode,dp));
        }
        FloatType::set_rounding_to_nearest();
    }
}

template<class PR> Void
TestFloat<PR>::test_predicates()
{
    PR pr=precision;
    ARIADNE_TEST_BINARY_PREDICATE(eq,Float<PR>(-1,pr),Float<PR>(-1,pr));
    ARIADNE_TEST_BINARY_PREDICATE(not eq,Float<PR>(-1,pr),Float<PR>(-2,pr));
    ARIADNE_TEST_BINARY_PREDICATE(lt,Float<PR>(-2,pr),Float<PR>(-1,pr));
    ARIADNE_TEST_BINARY_PREDICATE(not lt,Float<PR>(-2,pr),Float<PR>(-2,pr));
    ARIADNE_TEST_BINARY_PREDICATE(not lt,Float<PR>(-1,pr),Float<PR>(-2,pr));

    Dyadic w(3,1u);
    Integer zl(1), zu(2);
    Dyadic wl(5,2u), wu(7,2u);
    Rational ql(4,3), qu(5,3);

    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),w),cmp(w,w));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),zl),cmp(w,zl));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),zu),cmp(w,zu));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),wl),cmp(w,wl));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),wu),cmp(w,wu));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),ql),cmp(w,ql));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),qu),cmp(w,qu));

    if constexpr (Same<PR,DoublePrecision>) {
        Rational zero(0);
        ARIADNE_TEST_EQUALS(cmp(FloatDP::inf(Sign::POSITIVE,dp),zero),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(FloatDP::inf(Sign::NEGATIVE,dp),zero),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(zero,FloatDP::inf(Sign::POSITIVE,dp)),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(zero,FloatDP::inf(Sign::NEGATIVE,dp)),Comparison::GREATER);
    }

    if constexpr (Same<PR,MultiplePrecision>) {
        Float<PR> x(Dyadic(3,1u),pr);
        Float<PR> one(1,pr);
        Float<PR> two_native(2,pr);

        ARIADNE_TEST_EQUALS(cmp(two_native,one),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(two_native,Nat(1u)),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(one,Nat(2u)),Comparison::LESS);

        ARIADNE_TEST_EQUALS(cmp(Int(1),two_native),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(Int(2),two_native),Comparison::EQUAL);
        ARIADNE_TEST_EQUALS(cmp(Int(3),two_native),Comparison::GREATER);

        ARIADNE_TEST_EQUALS(cmp(ExactDouble(1.0),two_native),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(ExactDouble(2.0),two_native),Comparison::EQUAL);
        ARIADNE_TEST_EQUALS(cmp(ExactDouble(3.0),two_native),Comparison::GREATER);

        ARIADNE_TEST_EQUALS(cmp(Nat(1u),x),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(Nat(2u),x),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(Nat(3u),Float<PR>(3,pr)),Comparison::EQUAL);

        FloatDP dp_low(Dyadic(1),dp);
        FloatDP dp_equal(Dyadic(3,1u),dp);
        FloatDP dp_high(Dyadic(2),dp);
        ARIADNE_TEST_EQUALS(cmp(x,dp_low),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(x,dp_equal),Comparison::EQUAL);
        ARIADNE_TEST_EQUALS(cmp(x,dp_high),Comparison::LESS);
    }
}


Void test_builtin_runtime()
{
    ApproximateDouble default_approximation;
    ApproximateDouble a(2.0);
    ApproximateDouble b(-3.0);
    ApproximateDouble literal=1.25_a;
    ARIADNE_TEST_EQUALS(default_approximation.get_d(),0.0);
    ARIADNE_TEST_EQUALS(static_cast<double>(a),2.0);
    ARIADNE_TEST_EQUALS(literal.get_d(),1.25);

    std::ostringstream approximate_stream;
    approximate_stream << a;
    ARIADNE_TEST_ASSERT(not approximate_stream.str().empty());

    ARIADNE_TEST_EQUALS(nul(a).get_d(),0.0);
    ARIADNE_TEST_EQUALS(pos(a).get_d(),2.0);
    ARIADNE_TEST_EQUALS(neg(a).get_d(),-2.0);
    ARIADNE_TEST_EQUALS(add(a,b).get_d(),-1.0);
    ARIADNE_TEST_EQUALS(sub(a,b).get_d(),5.0);
    ARIADNE_TEST_EQUALS(mul(a,b).get_d(),-6.0);
    ARIADNE_TEST_EQUALS(div(a,b).get_d(),-2.0/3.0);

    ARIADNE_TEST_EQUALS(abs(a).get_d(),2.0);
    ARIADNE_TEST_EQUALS(abs(b).get_d(),3.0);
    ARIADNE_TEST_EQUALS(max(a,b).get_d(),2.0);
    ARIADNE_TEST_EQUALS(max(b,a).get_d(),2.0);
    ARIADNE_TEST_EQUALS(min(a,b).get_d(),-3.0);
    ARIADNE_TEST_EQUALS(min(b,a).get_d(),-3.0);
    ARIADNE_TEST_EQUALS(mag(b).get_d(),3.0);
    ARIADNE_TEST_EQUALS(mig(b).get_d(),3.0);

    ARIADNE_TEST_EQUALS((+a).get_d(),2.0);
    ARIADNE_TEST_EQUALS((-a).get_d(),-2.0);
    ARIADNE_TEST_EQUALS((a+b).get_d(),-1.0);
    ARIADNE_TEST_EQUALS((a-b).get_d(),5.0);
    ARIADNE_TEST_EQUALS((a*b).get_d(),-6.0);
    ARIADNE_TEST_EQUALS((a/b).get_d(),-2.0/3.0);
    ARIADNE_TEST_EQUALS((a*TwoExp(2)).get_d(),8.0);
    ARIADNE_TEST_EQUALS((a/TwoExp(1)).get_d(),1.0);

    ApproximateDouble inplace(8.0);
    inplace+=ApproximateDouble(2.0);
    ARIADNE_TEST_EQUALS(inplace.get_d(),10.0);
    inplace-=ApproximateDouble(3.0);
    ARIADNE_TEST_EQUALS(inplace.get_d(),7.0);
    inplace*=ApproximateDouble(2.0);
    ARIADNE_TEST_EQUALS(inplace.get_d(),14.0);
    inplace/=ApproximateDouble(7.0);
    ARIADNE_TEST_EQUALS(inplace.get_d(),2.0);

    ARIADNE_TEST_ASSERT(same(a,ApproximateDouble(2.0)));
    ARIADNE_TEST_ASSERT(not same(a,b));
    ARIADNE_TEST_SAME(a==ApproximateDouble(2.0),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(a!=b,ApproximateKleenean(true));
    ARIADNE_TEST_SAME(a>b,ApproximateKleenean(true));
    ARIADNE_TEST_SAME(b<a,ApproximateKleenean(true));
    ARIADNE_TEST_SAME(a>=ApproximateDouble(2.0),ApproximateKleenean(true));
    ARIADNE_TEST_SAME(b<=a,ApproximateKleenean(true));

    ExactDouble exact_default;
    ExactDouble exact_two(2);
    ExactDouble exact_three(3u);
    ExactDouble exact_signed_ll(4LL);
    ExactDouble exact_unsigned_ll(5ULL);
    ARIADNE_TEST_EQUALS(exact_default.get_d(),0.0);
    ARIADNE_TEST_EQUALS(exact_two.get_d(),2.0);
    ARIADNE_TEST_EQUALS(exact_three.get_d(),3.0);
    ARIADNE_TEST_EQUALS(exact_signed_ll.get_d(),4.0);
    ARIADNE_TEST_EQUALS(exact_unsigned_ll.get_d(),5.0);

    ARIADNE_TEST_FAIL((ExactDouble{9007199254740993LL}));
    ARIADNE_TEST_FAIL((ExactDouble{std::numeric_limits<long long>::max()}));
    ARIADNE_TEST_FAIL((ExactDouble{9007199254740993ULL}));
    ARIADNE_TEST_FAIL((ExactDouble{std::numeric_limits<unsigned long long>::max()}));

    ExactDouble exact_nan(std::numeric_limits<double>::quiet_NaN());
    ARIADNE_TEST_ASSERT(std::isnan(exact_nan.get_d()));
    if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
        ARIADNE_TEST_FAIL((ExactDouble{0.1L}));
        ARIADNE_TEST_FAIL(exact(0.1L));
    }

    ARIADNE_TEST_ASSERT(std::isinf(ExactDouble::infinity().get_d()));
    ARIADNE_TEST_ASSERT(std::isinf(ExactDouble::inf().get_d()));
    ARIADNE_TEST_ASSERT(std::isnan(ExactDouble::nan().get_d()));
    ARIADNE_TEST_EQUALS(nul(exact_two).get_d(),0.0);
    ARIADNE_TEST_EQUALS(abs(ExactDouble(-2)).get_d(),2.0);
    ARIADNE_TEST_ASSERT(is_finite(exact_two));
    ARIADNE_TEST_ASSERT(not is_finite(ExactDouble::inf()));

    ARIADNE_TEST_EQUALS((ExactDouble(2)*TwoExp(3)).get_d(),16.0);
    ARIADNE_TEST_EQUALS((ExactDouble(2)/TwoExp(1)).get_d(),1.0);
    ARIADNE_TEST_EQUALS(cmp(ExactDouble(2),ExactDouble(2)),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(ExactDouble(1),ExactDouble(2)),Comparison::LESS);
    ARIADNE_TEST_EQUALS(cmp(ExactDouble(3),ExactDouble(2)),Comparison::GREATER);

    ARIADNE_TEST_ASSERT(ExactDouble(2)==ExactDouble(2));
    ARIADNE_TEST_ASSERT(ExactDouble(2)!=ExactDouble(3));
    ARIADNE_TEST_ASSERT(ExactDouble(3)>=ExactDouble(2));
    ARIADNE_TEST_ASSERT(ExactDouble(2)<=ExactDouble(3));
    ARIADNE_TEST_ASSERT(ExactDouble(3)>ExactDouble(2));
    ARIADNE_TEST_ASSERT(ExactDouble(2)<ExactDouble(3));

    ARIADNE_TEST_EQUALS(exact(0.625L).get_d(),0.625);
    ARIADNE_TEST_EQUALS(cast_exact(ApproximateDouble(1.5)).get_d(),1.5);
    ARIADNE_TEST_EQUALS(cast_positive(ExactDouble(2)).get_d(),2.0);
    ARIADNE_TEST_EQUALS(Positive<ExactDouble>(ExactDouble(3)).get_d(),3.0);
    ARIADNE_TEST_EQUALS(ApproximateDouble(ExactDouble(1.25)).get_d(),1.25);
}

Void test_double_runtime()
{
    set_default_builtin_rounding();
    ARIADNE_TEST_EQUALS(get_builtin_rounding_mode(),ROUND_UPWARD);
    set_builtin_rounding_to_nearest();

    ARIADNE_TEST_ASSERT(texp(1.0)>2.7 && texp(1.0)<2.8);
    ARIADNE_TEST_EQUALS(add_opp(1.0,2.0),3.0);
    ARIADNE_TEST_EQUALS(neg_rec_rnd(2.0),-0.5);
    ARIADNE_TEST_EQUALS(neg_rec_opp(2.0),-0.5);

    ARIADNE_TEST_EQUALS(pow_rnd(-2.0,Int(-2)),0.25);
    ARIADNE_TEST_EQUALS(pow_rnd(-2.0,Nat(2u)),4.0);
    ARIADNE_TEST_FAIL(pow_rnd(0.0,Int(-1)));

    const double infinity_value=std::numeric_limits<double>::infinity();
    const double nan=std::numeric_limits<double>::quiet_NaN();
    ARIADNE_TEST_ASSERT(std::isnan(sqrt_rnd(nan)));
    ARIADNE_TEST_ASSERT(std::isinf(sqrt_rnd(infinity_value)));
    ARIADNE_TEST_EQUALS(sqrt_rnd(0.0),0.0);
    ARIADNE_TEST_FAIL(sqrt_rnd(-1.0));
    ARIADNE_TEST_ASSERT(std::isnan(exp_rnd(nan)));
    ARIADNE_TEST_ASSERT(std::isinf(exp_rnd(infinity_value)));
    ARIADNE_TEST_EQUALS(exp_rnd(-infinity_value),0.0);
    ARIADNE_TEST_ASSERT(std::isinf(log_rnd(0.0)));
    ARIADNE_TEST_ASSERT(std::signbit(log_rnd(0.0)));
    ARIADNE_TEST_FAIL(log_rnd(-1.0));
    ARIADNE_TEST_FAIL(atan_rnd_series(0.5));

    set_builtin_rounding_toward_zero();
    ARIADNE_TEST_EXECUTE(pi_rnd());
    ARIADNE_TEST_EXECUTE(pi_opp());
    set_builtin_rounding_to_nearest();

    const double sin_inputs[] = {-4.0,-2.8,-2.0,-1.2,-0.5,0.1,1.0,1.7,2.5,4.0};
    for(double x : sin_inputs) {
        ARIADNE_TEST_ASSERT(std::abs(sin_rnd(x)-std::sin(x))<1e-12);
    }
    ARIADNE_TEST_ASSERT(std::isnan(sin_rnd(infinity_value)));

    ARIADNE_TEST_ASSERT(std::abs(cos_rnd(0.1)-std::cos(0.1))<1e-12);
    ARIADNE_TEST_ASSERT(std::isnan(cos_rnd(infinity_value)));

    ARIADNE_TEST_ASSERT(std::abs(tan_rnd(-4.0)-std::tan(-4.0))<1e-11);
    ARIADNE_TEST_ASSERT(std::isnan(tan_rnd(infinity_value)));

    set_builtin_rounding_toward_zero();
    ARIADNE_TEST_EXECUTE(atan_rnd(-0.5));
    set_builtin_rounding_to_nearest();
}

// D46: exact FloatMP -> Dyadic conversion without mpf_t intermediary.
Void test_floatmp_exact_dyadic_export() {
    const auto precision=MultiplePrecision(256_bits);
    for (unsigned width:{1u,53u,127u,128u,191u}) {
        const BigUnsigned magnitude=add(
            shift_left(BigUnsigned(1u),width-1u),BigUnsigned(13u));
        for (unsigned denominator_exponent:{0u,1u,7u,63u,127u}) {
            for (bool negative:{false,true}) {
                const Integer coefficient=Integer::from_magnitude(
                    magnitude,negative);
                const Dyadic source(coefficient,denominator_exponent);
                const FloatMP exact(source,precision);
                const Dyadic restored=static_cast<Dyadic>(exact);
                ARIADNE_TEST_EQUALS(Rational(restored),Rational(source));
                ARIADNE_TEST_EQUALS(Rational(restored),
                                    static_cast<Rational>(exact));
            }
        }
    }
    ARIADNE_TEST_ASSERT(is_nan(static_cast<Dyadic>(
        FloatMP::nan(precision))));
    ARIADNE_TEST_ASSERT(is_inf(static_cast<Dyadic>(
        FloatMP::inf(Sign::POSITIVE,precision))));
    ARIADNE_TEST_EQUALS(static_cast<Dyadic>(FloatMP(0,precision)),Dyadic(0));
}

// Native-storage regression uses exact binary significand/exponent facts,
// not another floating implementation or GMP as the test oracle.
// Frozen mathematical expectations for direct native ratio rounding.
// In particular, the huge common factor exercises a non-reduced
// numerator/denominator without constructing a canonical Rational.
Void test_floatmp_native_ratio_rounding() {
    const MultiplePrecision p(7_bits);
    const auto nearest=FloatMP::ROUND_TO_NEAREST;
    const auto round_up=FloatMP::ROUND_UPWARD;
    const auto round_down=FloatMP::ROUND_DOWNWARD;
    const FloatMP thirds=FloatMP::from_ratio(Integer(1),Integer(3),nearest,p);
    const FloatMP thirds_up=FloatMP::from_ratio(Integer(1),Integer(3),round_up,p);
    const FloatMP negative_thirds=FloatMP::from_ratio(Integer(-1),Integer(3),round_down,p);
    const FloatMP negative_denominator=FloatMP::from_ratio(
        Integer(1),Integer(-3),round_up,p);
    ARIADNE_TEST_ASSERT(compare(thirds.significand(),BigUnsigned(85u))==0);
    ARIADNE_TEST_ASSERT(compare(thirds_up.significand(),BigUnsigned(86u))==0);
    ARIADNE_TEST_ASSERT(compare(negative_thirds.significand(),BigUnsigned(86u))==0);
    ARIADNE_TEST_ASSERT(compare(negative_denominator.significand(),BigUnsigned(85u))==0);
    ARIADNE_TEST_ASSERT(thirds.exponent()==-8);
    ARIADNE_TEST_ASSERT(thirds_up.exponent()==-8);
    ARIADNE_TEST_ASSERT(negative_thirds.exponent()==-8);
    ARIADNE_TEST_ASSERT(!thirds.negative() && !thirds_up.negative());
    ARIADNE_TEST_ASSERT(negative_thirds.negative() && negative_denominator.negative());

    const BigUnsigned common=shift_left(BigUnsigned(1u),256u);
    const Integer numerator=Integer::from_magnitude(common,false);
    const Integer denominator=Integer::from_magnitude(
        mul(common,BigUnsigned(3u)),false);
    const FloatMP nonreduced=FloatMP::from_ratio(numerator,denominator,nearest,p);
    ARIADNE_TEST_ASSERT(compare(nonreduced.significand(),BigUnsigned(85u))==0);
    ARIADNE_TEST_ASSERT(nonreduced.exponent()==-8);
    const FloatMP exact=FloatMP::from_ratio(Integer(5),Integer(8),nearest,p);
    ARIADNE_TEST_ASSERT(compare(exact.significand(),BigUnsigned(80u))==0);
    ARIADNE_TEST_ASSERT(exact.exponent()==-7);

    bool threw=false;
    try { (void)FloatMP::from_ratio(Integer(1),Integer(0),nearest,p); }
    catch(std::domain_error const&) { threw=true; }
    ARIADNE_TEST_ASSERT(threw);
}

// Decimal-format oracle: exact binary rationals with frozen expected
// strings. No GMP or other floating implementation generates the answers.
Void test_native_decimal_format() {
    const MultiplePrecision precision(53_bits);
    const FloatMP one_eighth(Rational(Integer(1),Integer(8)),
                             FloatMP::ROUND_TO_NEAREST,precision);
    const FloatMP negative_eighth(Rational(Integer(-1),Integer(8)),
                                  FloatMP::ROUND_TO_NEAREST,precision);
    const FloatMP three_eighths(Rational(Integer(3),Integer(8)),
                                FloatMP::ROUND_TO_NEAREST,precision);
    ARIADNE_TEST_EQUALS(print(one_eighth,DecimalPlaces(2u),
                              FloatMP::ROUND_TO_NEAREST),String("0.12"));
    ARIADNE_TEST_EQUALS(print(one_eighth,DecimalPlaces(2u),
                              FloatMP::ROUND_UPWARD),String("0.13"));
    ARIADNE_TEST_EQUALS(print(one_eighth,DecimalPlaces(3u),
                              FloatMP::ROUND_TO_NEAREST),String("0.125"));
    ARIADNE_TEST_EQUALS(print(negative_eighth,DecimalPlaces(2u),
                              FloatMP::ROUND_TO_NEAREST),String("-0.12"));
    ARIADNE_TEST_EQUALS(print(negative_eighth,DecimalPlaces(2u),
                              FloatMP::ROUND_DOWNWARD),String("-0.13"));
    ARIADNE_TEST_EQUALS(print(three_eighths,DecimalPlaces(2u),
                              FloatMP::ROUND_TO_NEAREST),String("0.38"));
    ARIADNE_TEST_EQUALS(print(three_eighths,DecimalPlaces(2u),
                              FloatMP::ROUND_TOWARD_ZERO),String("0.37"));
    ARIADNE_TEST_EQUALS(print(FloatMP(2,precision),DecimalPlaces(0u),
                              FloatMP::ROUND_TO_NEAREST),String("2."));
    ARIADNE_TEST_EQUALS(print(FloatMP::zero(53u,true),DecimalPlaces(2u),
                              FloatMP::ROUND_TO_NEAREST),String("-0.00"));
    ARIADNE_TEST_EQUALS(print(FloatMP::infinity(53u,true),DecimalPlaces(2u),
                              FloatMP::ROUND_TO_NEAREST),String("-inf"));
    ARIADNE_TEST_EQUALS(print(FloatMP::nan(53u),DecimalPlaces(2u),
                              FloatMP::ROUND_TO_NEAREST),String("nan"));

    // The separate binary64 formatter also preserves exact halfway
    // rounding and the signed-zero payload.
    StringStream dp_nearest_stream,dp_upward_stream,dp_tie_stream,dp_zero_stream;
    write(dp_nearest_stream,FloatDP(ExactDouble(0.125),dp),DecimalPlaces(2u),FloatDP::ROUND_TO_NEAREST);
    write(dp_upward_stream,FloatDP(ExactDouble(0.125),dp),DecimalPlaces(2u),FloatDP::ROUND_UPWARD);
    write(dp_tie_stream,FloatDP(ExactDouble(0.375),dp),DecimalPlaces(2u),FloatDP::ROUND_TO_NEAREST);
    write(dp_zero_stream,FloatDP(ExactDouble(-0.0),dp),DecimalPlaces(2u),FloatDP::ROUND_TO_NEAREST);
    ARIADNE_TEST_EQUALS(dp_nearest_stream.str(),String("0.12"));
    ARIADNE_TEST_EQUALS(dp_upward_stream.str(),String("0.13"));
    ARIADNE_TEST_EQUALS(dp_tie_stream.str(),String("0.38"));
    ARIADNE_TEST_EQUALS(dp_zero_stream.str(),String("-0.00"));
}

Void test_floatmp_native_significand_storage() {
    for (std::size_t bits:{7u,53u,128u,256u}) {
        const MultiplePrecision precision(bits);
        const BigUnsigned leading=shift_left(BigUnsigned(1u),bits-1u);
        const BigUnsigned all_ones=sub(shift_left(BigUnsigned(1u),bits),
                                      BigUnsigned(1u));
        const FloatMP one(std::int64_t{1},precision);
        const FloatMP negative_one(std::int64_t{-1},precision);
        ARIADNE_TEST_ASSERT(compare(one.significand(),leading)==0);
        ARIADNE_TEST_ASSERT(compare(negative_one.significand(),leading)==0);
        ARIADNE_TEST_ASSERT(!one.negative() && negative_one.negative());
        ARIADNE_TEST_ASSERT(one.exponent()==-static_cast<std::int64_t>(bits-1u));
        ARIADNE_TEST_ASSERT(negative_one.exponent()==one.exponent());
        ARIADNE_TEST_EQUALS(one.get_d(),1.0);
        ARIADNE_TEST_EQUALS(negative_one.get_d(),-1.0);

        // Immediately below +1 and immediately above -1, the p-bit
        // significand wraps to 2^p-1 with the exponent decremented.
        const FloatMP before_one=next(RoundDownward{},one);
        const FloatMP after_negative_one=next(RoundUpward{},negative_one);
        ARIADNE_TEST_ASSERT(compare(before_one.significand(),all_ones)==0);
        ARIADNE_TEST_ASSERT(compare(after_negative_one.significand(),all_ones)==0);
        ARIADNE_TEST_ASSERT(before_one.exponent()==one.exponent()-1);
        ARIADNE_TEST_ASSERT(after_negative_one.exponent()==one.exponent()-1);

        const FloatMP after_one=next(RoundUpward{},one);
        const FloatMP before_negative_one=next(RoundDownward{},negative_one);
        const BigUnsigned incremented=add(leading,BigUnsigned(1u));
        ARIADNE_TEST_ASSERT(compare(after_one.significand(),incremented)==0);
        ARIADNE_TEST_ASSERT(compare(before_negative_one.significand(),incremented)==0);
        ARIADNE_TEST_ASSERT(after_one.exponent()==one.exponent());
        ARIADNE_TEST_ASSERT(before_negative_one.exponent()==one.exponent());

        // Storage copy/assignment and canonical signed zero are explicit.
        FloatMP copied=one;
        ARIADNE_TEST_ASSERT(compare(copied.significand(),leading)==0);
        copied=negative_one;
        ARIADNE_TEST_ASSERT(copied.negative());
        ARIADNE_TEST_ASSERT(compare(copied.significand(),leading)==0);
        const FloatMP negative_zero=FloatMP::zero(
            static_cast<unsigned long>(bits),true);
        ARIADNE_TEST_ASSERT(negative_zero.is_zero() && negative_zero.negative());
        ARIADNE_TEST_ASSERT(negative_zero.significand().is_zero());

        const FloatMP smallest=FloatMP::min(precision);
        const FloatMP largest=FloatMP::max(precision);
        ARIADNE_TEST_ASSERT(compare(smallest.significand(),leading)==0);
        ARIADNE_TEST_ASSERT(compare(largest.significand(),all_ones)==0);
        ARIADNE_TEST_ASSERT(smallest.exponent()==
            std::numeric_limits<FloatMP::Exponent>::min());
        ARIADNE_TEST_ASSERT(largest.exponent()==
            std::numeric_limits<FloatMP::Exponent>::max());
    }
}

Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);
    
    test_float_literals();
    test_floatmp_directed_bound_oracles();
    test_builtin_runtime();
    test_double_runtime();
    TestFloat<DoublePrecision>(dp).test();
        TestFloat<MultiplePrecision>(MultiplePrecision(128_bits)).test();
    test_floatmp_exact_dyadic_export();
    test_floatmp_native_significand_storage();
    test_floatmp_native_ratio_rounding();
    test_native_decimal_format();


    return ARIADNE_TEST_FAILURES;
}

