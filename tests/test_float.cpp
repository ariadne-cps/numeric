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
#include <array>
#include "numeric/detail/checked_arithmetic.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <bit>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>
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

namespace Ariadne { int abslog10floor(double); }

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
    ARIADNE_TEST_EQUALS(abslog10floor(0.0),std::numeric_limits<int>::min());
    ARIADNE_TEST_EQUALS(abslog10floor(-0.0),std::numeric_limits<int>::min());
    ARIADNE_TEST_EQUALS(abslog10floor(infinity_value),0);
    ARIADNE_TEST_EQUALS(abslog10floor(nan),0);
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

    set_builtin_rounding_upward();
    const double positive_atan_upper=atan_rnd(1.0);
    set_builtin_rounding_downward();
    ARIADNE_TEST_EQUALS(atan_rnd(-1.0),-positive_atan_upper);
    ARIADNE_TEST_EQUALS(get_builtin_rounding_mode(),ROUND_DOWNWARD);
    set_builtin_rounding_to_nearest();
    const double positive_atan_nearest=atan_rnd(1.0);
    ARIADNE_TEST_EQUALS(atan_rnd(-1.0),-positive_atan_nearest);
    ARIADNE_TEST_EQUALS(get_builtin_rounding_mode(),ROUND_TO_NEAREST);
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
                const Dyadic restored=exact.operator Dyadic();
                const Dyadic constructed(exact);
                ARIADNE_TEST_EQUALS(Rational(constructed),Rational(source));
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
    ARIADNE_TEST_ASSERT(is_nan(FloatMP::nan(precision).operator Dyadic()));
    for(bool negative:{false,true}) {
        const auto infinite=FloatMP::infinity(256u,negative).operator Dyadic();
        ARIADNE_TEST_ASSERT(is_inf(infinite));
        ARIADNE_TEST_EQUALS(std::signbit(infinite.get_d()),negative);
        ARIADNE_TEST_EQUALS(FloatMP::zero(256u,negative).operator Dyadic(),Dyadic(0));
    }
    const FloatMP positive_scale(TwoExp(512),precision);
    ARIADNE_TEST_EQUALS(positive_scale.operator Dyadic(),Dyadic(TwoExp(512)));

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
    for (auto value : {FloatDP::nan(dp),FloatDP::inf(Sign::POSITIVE,dp),
                       FloatDP::inf(Sign::NEGATIVE,dp)}) {
        StringStream stream;
        write(stream,value,DecimalPlaces(2u),FloatDP::ROUND_TO_NEAREST);
        const String expected=is_nan(value) ? String("nan") :
                              value.get_d()<0.0 ? String("-inf") : String("inf");
        ARIADNE_TEST_EQUALS(stream.str(),expected);
    }
    const auto old_mode=FloatDP::get_rounding_mode();
    FloatDP::set_rounding_to_nearest();
    const FloatDP zero(ExactDouble(-0.0),dp);
    ARIADNE_TEST_ASSERT(std::signbit(tanh(CurrentRoundingMode{},zero).get_d()));
    FloatDP::set_rounding_mode(old_mode);
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

// A failed explicitly rounded operation must preserve its caller's mode.
void test_floatdp_rounding_exception_safety() {
    RoundingSession<FloatDP> restore;
    for (auto mode:{FloatDP::ROUND_DOWNWARD,FloatDP::ROUND_TO_NEAREST}) {
        restore.set(mode);
        ARIADNE_TEST_FAIL(log(FloatDP::ROUND_UPWARD,FloatDP(-0.5_x,dp)));
        ARIADNE_TEST_EQUALS(FloatDP::get_rounding_mode(),mode);
        restore.set(FloatDP::ROUND_UPWARD);
        restore.set(mode);
        ARIADNE_TEST_FAIL(pow(FloatDP::ROUND_UPWARD,FloatDP(0,dp),Int(-1)));
        ARIADNE_TEST_EQUALS(FloatDP::get_rounding_mode(),mode);
        restore.set(FloatDP::ROUND_UPWARD);
    }
}

// Invalid raw mode integers must be rejected before touching the hardware.
void test_floatdp_rounding_mode_contract() {
    const auto initial=FloatDP::get_rounding_mode();
    const auto invalid=std::numeric_limits<FloatDP::RoundingModeType>::max();
    ARIADNE_TEST_FAIL(FloatDP::set_rounding_mode(invalid));
    ARIADNE_TEST_EQUALS(FloatDP::get_rounding_mode(),initial);
    ARIADNE_TEST_FAIL(sin(invalid,FloatDP(0.5_x,dp)));
    ARIADNE_TEST_EQUALS(FloatDP::get_rounding_mode(),initial);
    ARIADNE_TEST_FAIL(cos(invalid,FloatDP(2,dp)));
    ARIADNE_TEST_EQUALS(FloatDP::get_rounding_mode(),initial);
}

void test_floatdp_sine_outside_certificate() {
    const auto initial=FloatDP::get_rounding_mode();
    const FloatDP input(0.75_x,dp);
    for (auto mode:{FloatDP::ROUND_TO_NEAREST,FloatDP::ROUND_DOWNWARD,
                    FloatDP::ROUND_UPWARD,FloatDP::ROUND_TOWARD_ZERO}) {
        const auto sine_value=sin(mode,input);
        // For 0<x<1, x-x^3/6 < sin(x) < x. Both bounds are exact dyadics.
        ARIADNE_TEST_ASSERT(sine_value>FloatDP(0.6796875_x,dp));
        ARIADNE_TEST_ASSERT(sine_value<input);
        ARIADNE_TEST_EQUALS(FloatDP::get_rounding_mode(),initial);
    }
}

void test_floatdp_unary_wrappers() {
    std::ifstream data(ARIADNE_TEST_DATA_DIR "/floatdp_wrappers.tsv");
    if (!data) throw std::runtime_error("Missing floatdp_wrappers.tsv");
    std::string line;
    std::set<std::uint64_t> inputs;
    while (std::getline(data,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        std::uint64_t input,zero,positive,negative,half;
        std::string extra;
        if (!(fields>>std::hex>>input>>zero>>positive>>negative>>half) || (fields>>extra)
            || !inputs.insert(input).second)
            throw std::runtime_error("Malformed floatdp_wrappers.tsv");
        const ExactDouble source(std::bit_cast<double>(input));
        for (auto mode:{ROUND_TO_NEAREST,ROUND_DOWNWARD,ROUND_UPWARD,ROUND_TOWARD_ZERO}) {
            const FloatDP x(source,mode,dp);
            ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(nul(mode,x).get_d()),zero);
            ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(pos(mode,x).get_d()),positive);
            ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(neg(mode,x).get_d()),negative);
            ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(hlf(mode,x).get_d()),half);
        }
    }
    if (inputs.size()!=8u) throw std::runtime_error("Incomplete floatdp_wrappers.tsv");
}

void test_floatmp_minmax_oracle() {
    std::ifstream data(ARIADNE_TEST_DATA_DIR "/floatmp_minmax.tsv");
    if (!data) throw std::runtime_error("Missing floatmp_minmax.tsv");
    std::string line;
    std::set<unsigned> ids;
    const auto saved=FloatMP::get_rounding_mode();
    while (std::getline(data,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        unsigned id,ap,bp,p;
        int negative;
        char mode,kind;
        std::string op,a,b,n,d,extra;
        if (!(fields>>id>>op>>mode>>ap>>a>>bp>>b>>p>>kind>>n>>d>>negative)
            || (fields>>extra) || !ids.insert(id).second)
            throw std::runtime_error("Malformed floatmp_minmax.tsv");
        const auto rounding=mode=='N'?FloatMP::ROUND_TO_NEAREST:mode=='D'?FloatMP::ROUND_DOWNWARD:
                            mode=='U'?FloatMP::ROUND_UPWARD:FloatMP::ROUND_TOWARD_ZERO;
        const auto input=[](std::string const& source,unsigned bits) {
            const MultiplePrecision precision(static_cast<std::size_t>(bits));
            if (source=="nan") return FloatMP::nan(precision);
            if (source=="inf" || source=="-inf") return FloatMP::inf(source[0]=='-'?Sign::NEGATIVE:Sign::POSITIVE,precision);
            if (source=="+0" || source=="-0") return FloatMP::zero(static_cast<unsigned long>(bits),source[0]=='-');
            return FloatMP(Rational(String(source)),FloatMP::ROUND_TO_NEAREST,precision);
        };
        const FloatMP left=input(a,ap),right=input(b,bp);
        const FloatMP actual=op=="min"?min(rounding,left,right):max(rounding,left,right);
        std::cout << "FloatMP min/max oracle row " << id << std::endl;
        ARIADNE_TEST_EQUALS(actual.precision().bits(),static_cast<std::size_t>(p));
        ARIADNE_TEST_ASSERT(FloatMP::get_rounding_mode()==saved);
        if (kind=='N') {
            ARIADNE_TEST_ASSERT(actual.kind()==FloatMP::Kind::NaN);
            ARIADNE_TEST_ASSERT(is_nan(actual));
        } else if (kind=='I') {
            ARIADNE_TEST_ASSERT(actual.kind()==FloatMP::Kind::Infinity);
            ARIADNE_TEST_ASSERT(is_inf(actual));
            ARIADNE_TEST_EQUALS(actual.negative(),negative!=0);
        } else if (kind=='F') {
            ARIADNE_TEST_ASSERT(actual.kind()==FloatMP::Kind::Finite);
            ARIADNE_TEST_ASSERT(is_finite(actual));
            ARIADNE_TEST_EQUALS(static_cast<Rational>(actual),(Rational(Integer(String(n)),Integer(String(d)))));
            ARIADNE_TEST_EQUALS(actual.negative(),negative!=0);
        } else throw std::runtime_error("Unknown FloatMP min/max actual kind");
    }
    if (ids.size()!=648u || *ids.begin()!=0u || *ids.rbegin()!=647u)
        throw std::runtime_error("Incomplete floatmp_minmax.tsv");
}

void test_floatmp_construction_oracle() {
    std::ifstream data(ARIADNE_TEST_DATA_DIR "/floatmp_construction.tsv");
    if (!data) throw std::runtime_error("Missing floatmp_construction.tsv");
    std::string line;
    std::set<unsigned> ids;
    while (std::getline(data,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        unsigned id,bits;
        int n,d;
        std::string en,ed,extra;
        if (!(fields>>id>>n>>d>>bits>>en>>ed) || (fields>>extra) || !ids.insert(id).second)
            throw std::runtime_error("Malformed floatmp_construction.tsv");
        const MultiplePrecision precision(static_cast<std::size_t>(bits));
        const Rational expected{Integer(String(en)),Integer(String(ed))};
        std::cout << "FloatMP construction oracle row " << id << std::endl;
        const FloatMP from_rational(Rational(n,d),precision);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(from_rational),expected);
        ARIADNE_TEST_EQUALS(from_rational.precision(),precision);
        const FloatMP empty(precision);
        ARIADNE_TEST_ASSERT(empty.is_zero() && !empty.negative());
        if (d==1) {
            const FloatMP from_integer(Integer(n),precision);
            ARIADNE_TEST_EQUALS(static_cast<Rational>(from_integer),expected);
            for (const auto mode:{FloatMP::ROUND_TO_NEAREST,FloatMP::ROUND_DOWNWARD,
                                  FloatMP::ROUND_UPWARD,FloatMP::ROUND_TOWARD_ZERO}) {
                ARIADNE_TEST_EQUALS(static_cast<Rational>(FloatMP(n,mode,precision)),expected);
                if (n>=0) {
                    const auto natural=static_cast<unsigned>(n);
                    ARIADNE_TEST_EQUALS(static_cast<Rational>(FloatMP(natural,mode,precision)),expected);
                }
            }
        }
        FloatMP self_move=from_rational;
        auto* alias=&self_move;
        self_move=std::move(*alias);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(self_move),expected);
        ARIADNE_TEST_EQUALS(self_move.precision(),precision);
    }
    if (ids.size()!=42u || *ids.begin()!=0u || *ids.rbegin()!=41u)
        throw std::runtime_error("Incomplete floatmp_construction.tsv");
}

void test_floatmp_precision_contracts() {
    const MultiplePrecision low(2u),high(7u);
    ARIADNE_TEST_ASSERT(low<high && !(high<low));
    ARIADNE_TEST_ASSERT(high>low && !(low>high));
    ARIADNE_TEST_ASSERT(high>=low && low>=low && !(low>=high));
    ARIADNE_TEST_EQUALS(max(low,high),high);
    ARIADNE_TEST_EQUALS(max(high,low),high);
    std::ostringstream rendered;
    rendered << high;
    ARIADNE_TEST_EQUALS(rendered.str(),std::string("mp(7)"));
    ARIADNE_TEST_FAIL(FloatMP(0ul));
    ARIADNE_TEST_FAIL(FloatMP(1ul));
    ARIADNE_TEST_FAIL((FloatMP(MultiplePrecision(1u))));
    if constexpr (std::numeric_limits<unsigned long>::digits>std::numeric_limits<FloatMP::Exponent>::digits) {
        ARIADNE_TEST_FAIL(static_cast<void>(FloatMP(std::numeric_limits<unsigned long>::max())));
    }
    ARIADNE_TEST_EQUALS(floatmp_exponent_from_bit_length(0u),std::int64_t(0));
    if constexpr (std::numeric_limits<std::size_t>::digits>std::numeric_limits<std::int64_t>::digits) {
        const auto boundary=static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max());
        ARIADNE_TEST_EQUALS(floatmp_exponent_from_bit_length(boundary),std::numeric_limits<std::int64_t>::max());
        ARIADNE_TEST_FAIL(floatmp_exponent_from_bit_length(boundary+1u));
    }
}

void test_floatmp_exponent_limit_oracle() {
    std::ifstream data(ARIADNE_TEST_DATA_DIR "/floatmp_exponent_limits.tsv");
    if (!data) throw std::runtime_error("Missing floatmp_exponent_limits.tsv");
    std::string line;
    std::set<unsigned> ids;
    while (std::getline(data,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        unsigned id,bits;
        char mode;
        FloatMP::Exponent expected_exponent;
        std::string op,a,b,status,magnitude,extra;
        if (!(fields>>id>>op>>mode>>bits>>a>>b>>status>>magnitude>>expected_exponent)
            || (fields>>extra) || !ids.insert(id).second)
            throw std::runtime_error("Malformed floatmp_exponent_limits.tsv");
        const MultiplePrecision precision(static_cast<std::size_t>(bits));
        const auto input=[&](std::string const& token) {
            if (token=="min") return FloatMP::min(precision);
            if (token=="max") return FloatMP::max(precision);
            throw std::runtime_error("Unknown native exponent source");
        };
        const auto rounding=mode=='N'?FloatMP::ROUND_TO_NEAREST:mode=='D'?FloatMP::ROUND_DOWNWARD:
                            mode=='U'?FloatMP::ROUND_UPWARD:FloatMP::ROUND_TOWARD_ZERO;
        const FloatMP left=input(a);
        bool overflow=false;
        std::cout << "FloatMP exponent limit row " << id << std::endl;
        try {
            const auto evaluate=[&]() {
                if (op=="shift") return shft(left,static_cast<Int>(std::stoi(b)));
                const FloatMP right=input(b);
                if (op=="div") return div(rounding,left,right);
                if (op=="mul") return mul(rounding,left,right);
                throw std::runtime_error("Unknown native exponent operation");
            };
            const FloatMP actual=evaluate();
            ARIADNE_TEST_EQUALS(status,std::string("finite"));
            ARIADNE_TEST_ASSERT(actual.kind()==FloatMP::Kind::Finite);
            ARIADNE_TEST_ASSERT(!actual.negative());
            ARIADNE_TEST_EQUALS(actual.precision(),precision);
            ARIADNE_TEST_EQUALS(actual.exponent(),expected_exponent);
            ARIADNE_TEST_EQUALS(Integer::from_magnitude(actual.significand(),false),Integer(String(magnitude)));
        } catch (std::overflow_error const&) { overflow=true; }
        ARIADNE_TEST_EQUALS(overflow,status=="overflow");
    }
    if (ids.size()!=76u || *ids.begin()!=0u || *ids.rbegin()!=75u)
        throw std::runtime_error("Incomplete floatmp_exponent_limits.tsv");
}

Void test_floatmp_special_arithmetic() {
    const auto old=FloatMP::get_rounding_mode();
    FloatMP::set_rounding_downward();
    const FloatMP::Rounding current(CurrentRoundingMode{});
    ARIADNE_TEST_ASSERT(current==FloatMP::ROUND_DOWNWARD);
    FloatMP::set_rounding_mode(old);
    const MultiplePrecision precision(128_bits);
    const FloatMP zero=FloatMP::zero(128u),one(1,precision),two_value(2,precision);
    const FloatMP pinf=FloatMP::infinity(128u),ninf=FloatMP::infinity(128u,true);
    const FloatMP nan=FloatMP::nan(128u);
    ARIADNE_TEST_ASSERT(is_nan(neg(nan)));
    ARIADNE_TEST_ASSERT(is_nan(add(near,nan,one)));
    ARIADNE_TEST_ASSERT(is_nan(add(near,one,nan)));
    ARIADNE_TEST_ASSERT(is_nan(add(near,pinf,ninf)));
    ARIADNE_TEST_EQUALS(add(near,pinf,pinf),pinf);
    ARIADNE_TEST_EQUALS(add(near,one,ninf),ninf);
    ARIADNE_TEST_EQUALS(add(near,ninf,one),ninf);
    ARIADNE_TEST_ASSERT(is_nan(mul(near,nan,one)));
    ARIADNE_TEST_ASSERT(is_nan(mul(near,one,nan)));
    ARIADNE_TEST_ASSERT(is_nan(mul(near,zero,pinf)));
    ARIADNE_TEST_ASSERT(is_nan(mul(near,ninf,zero)));
    ARIADNE_TEST_EQUALS(mul(near,one,ninf),ninf);
    ARIADNE_TEST_EQUALS(mul(near,ninf,-one),pinf);
    ARIADNE_TEST_ASSERT(is_nan(div(near,nan,one)));
    ARIADNE_TEST_ASSERT(is_nan(div(near,one,nan)));
    ARIADNE_TEST_ASSERT(is_nan(div(near,zero,zero)));
    ARIADNE_TEST_ASSERT(is_nan(div(near,pinf,ninf)));
    ARIADNE_TEST_EQUALS(div(near,ninf,one),ninf);
    ARIADNE_TEST_EQUALS(div(near,-one,zero),ninf);
    ARIADNE_TEST_ASSERT(div(near,zero,-one).negative());
    ARIADNE_TEST_ASSERT(div(near,one,ninf).negative());
    ARIADNE_TEST_ASSERT(is_nan(fma(near,nan,one,one)));
    ARIADNE_TEST_ASSERT(is_nan(fma(near,one,nan,one)));
    ARIADNE_TEST_ASSERT(is_nan(fma(near,one,one,nan)));
    ARIADNE_TEST_ASSERT(is_nan(fma(near,zero,pinf,one)));
    ARIADNE_TEST_ASSERT(is_nan(fma(near,pinf,zero,one)));
    ARIADNE_TEST_ASSERT(is_nan(fma(near,pinf,one,ninf)));
    ARIADNE_TEST_EQUALS(fma(near,pinf,one,pinf),pinf);
    ARIADNE_TEST_EQUALS(fma(near,one,ninf,one),ninf);
    ARIADNE_TEST_EQUALS(fma(near,one,one,ninf),ninf);
    ARIADNE_TEST_EQUALS(fma(near,zero,one,two_value),two_value);
    ARIADNE_TEST_EQUALS(fma(near,one,zero,two_value),two_value);
    for (auto value : {nan,pinf,ninf,FloatMP::zero(128u,true)}) {
        const FloatMP narrowed(value,near,MultiplePrecision(8_bits));
        ARIADNE_TEST_EQUALS(narrowed.precision(),MultiplePrecision(8_bits));
        ARIADNE_TEST_ASSERT(narrowed.is_nan()==value.is_nan());
        ARIADNE_TEST_ASSERT(narrowed.is_inf()==value.is_inf());
        ARIADNE_TEST_ASSERT(narrowed.negative()==value.negative());
        value.set_precision(precision); // unchanged precision is an explicit no-op.
        value.set_precision(MultiplePrecision(8_bits));
        ARIADNE_TEST_EQUALS(value.precision(),MultiplePrecision(8_bits));
        ARIADNE_TEST_ASSERT(value.is_nan()==narrowed.is_nan());
        ARIADNE_TEST_ASSERT(value.is_inf()==narrowed.is_inf());
        ARIADNE_TEST_ASSERT(value.negative()==narrowed.negative());
    }
    ARIADNE_TEST_ASSERT(is_nan(static_cast<Dyadic>(nan)));
    ARIADNE_TEST_ASSERT(is_inf(static_cast<Dyadic>(pinf)));
    ARIADNE_TEST_ASSERT(is_inf(static_cast<Dyadic>(ninf)));
    ARIADNE_TEST_EQUALS(static_cast<Dyadic>(zero),Dyadic(0));
    ARIADNE_TEST_EQUALS(static_cast<Dyadic>(FloatMP(256,MultiplePrecision(8_bits))),Dyadic(256));
    ARIADNE_TEST_FAIL(cast_integer(nan));
    ARIADNE_TEST_FAIL(cast_integer(pinf));
    StringStream implicit_repr,explicit_repr;
    repr(implicit_repr,two_value);
    repr(explicit_repr,two_value,FloatMP::get_rounding_mode());
    ARIADNE_TEST_EQUALS(implicit_repr.str(),explicit_repr.str());
    // Call the raw binary64 adapter, which is distinct from tagged tanh.
    ARIADNE_TEST_EQUALS(tanh_rnd(FloatDP(0,dp)),FloatDP(0,dp));
}

// Exact small integers isolate precision selection from numerical rounding.
Void test_floatmp_mixed_precision_arithmetic() {
    for (auto sizes : {std::pair{8u,16u},std::pair{16u,8u},std::pair{8u,8u},std::pair{8u,129u},std::pair{129u,8u}}) {
        const MultiplePrecision pa(std::size_t{sizes.first});
        const MultiplePrecision pb(std::size_t{sizes.second});
        const MultiplePrecision expected(std::size_t{std::min(sizes.first,sizes.second)});
        const FloatMP a(2,pa), b(4,pb);
        const auto check=[&](FloatMP const& result, int value) {
            ARIADNE_TEST_EQUALS(result.precision(),expected);
            ARIADNE_TEST_EQUALS(static_cast<Rational>(result),Rational(value));
        };
        check(add(near,a,b),6);
        check(mul(near,a,b),8);
        ARIADNE_TEST_EQUALS(div(near,a,b).precision(),expected);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(div(near,a,b)),Rational(1,2));
        check(med(near,a,b),3);
        check(rad(near,a,b),1);
        for (std::size_t c_bits : {4u,32u}) {
            const FloatMP c(1,MultiplePrecision(c_bits));
            const auto result=fma(near,a,b,c);
            ARIADNE_TEST_EQUALS(result.precision(),MultiplePrecision(std::min(expected.bits(),c_bits)));
            ARIADNE_TEST_EQUALS(static_cast<Rational>(result),Rational(9));
        }
    }
    const FloatMP value_two(2,MultiplePrecision(8_bits));
    ARIADNE_TEST_ASSERT(value_two>=Rational(1));
    ARIADNE_TEST_ASSERT(value_two>=Rational(2));
    ARIADNE_TEST_ASSERT(!(value_two>=Rational(3)));
    ARIADNE_TEST_ASSERT(!(FloatMP::nan(8u)>=Rational(2)));
}

// Exponent limits are scalar metadata: these neighbors need no huge allocation.
void test_floatmp_neighbor_limits() {
    for (unsigned bits:{2u,65u}) {
        const MultiplePrecision precision(std::size_t{bits});
        const auto smallest=FloatMP::min(precision);
        const auto largest=FloatMP::max(precision);
        const auto positive_infinity=FloatMP::inf(precision);
        const auto negative_infinity=FloatMP::inf(Sign::NEGATIVE,precision);
        const auto quiet_nan=FloatMP::nan(precision);
        ARIADNE_TEST_ASSERT(is_nan(next(up,quiet_nan)));
        ARIADNE_TEST_ASSERT(is_nan(next(down,quiet_nan)));
        ARIADNE_TEST_EQUALS(next(up,positive_infinity),positive_infinity);
        ARIADNE_TEST_EQUALS(next(down,negative_infinity),negative_infinity);
        ARIADNE_TEST_EQUALS(next(down,positive_infinity),largest);
        ARIADNE_TEST_EQUALS(next(up,negative_infinity),-largest);
        ARIADNE_TEST_EQUALS(next(up,largest),positive_infinity);
        ARIADNE_TEST_EQUALS(next(down,-largest),negative_infinity);
        for (bool negative:{false,true}) {
            const auto zero=FloatMP::zero(bits,negative);
            ARIADNE_TEST_EQUALS(next(up,zero),smallest);
            ARIADNE_TEST_EQUALS(next(down,zero),-smallest);
        }
        const auto negative_zero=next(up,-smallest);
        const auto positive_zero=next(down,smallest);
        ARIADNE_TEST_ASSERT(negative_zero.is_zero() && negative_zero.negative());
        ARIADNE_TEST_ASSERT(positive_zero.is_zero() && !positive_zero.negative());
        const auto first_above_min=next(up,smallest);
        ARIADNE_TEST_EQUALS(first_above_min.exponent(),smallest.exponent());
        ARIADNE_TEST_ASSERT(compare(first_above_min.significand(),
                                  add(smallest.significand(),BigUnsigned(1u)))==0);
        ARIADNE_TEST_EQUALS(next(down,first_above_min),smallest);
        ARIADNE_TEST_EQUALS(next(up,-first_above_min),-smallest);
        const auto last_below_max=next(down,largest);
        ARIADNE_TEST_EQUALS(next(up,last_below_max),largest);
        ARIADNE_TEST_EQUALS(next(down,-last_below_max),-largest);

        // At a power of two, the predecessor has twice as fine a spacing.
        const auto predecessor=next(down,FloatMP(1,precision));
        ARIADNE_TEST_EQUALS(predecessor.exponent(),-static_cast<FloatMP::Exponent>(bits));
        ARIADNE_TEST_ASSERT(compare(predecessor.significand(),
            sub(shift_left(BigUnsigned(1u),bits),BigUnsigned(1u)))==0);
        ARIADNE_TEST_EQUALS(next(up,predecessor),FloatMP(1,precision));
        ARIADNE_TEST_EQUALS(next(up,FloatMP(-1,precision)),-predecessor);
        ARIADNE_TEST_EQUALS(next(down,-predecessor),FloatMP(-1,precision));
    }
}

void test_floatmp_ratio_ties_and_carry() {
    const MultiplePrecision precision(2_bits);
    // n/3 exercises scalar normalization with zero and negative shifts.
    struct Row { int numerator; int nearest; int lower; int upper; };
    for (auto row:{Row{7,2,2,3},Row{15,4,4,6},Row{21,8,6,8},Row{100,32,32,48}}) {
        for (bool negative:{false,true}) {
            const Integer numerator(negative ? -row.numerator : row.numerator);
            const auto value=[&](auto mode) {
                return static_cast<Rational>(FloatMP::from_ratio(numerator,Integer(3),mode,precision));
            };
            ARIADNE_TEST_EQUALS(value(near),Rational(negative ? -row.nearest : row.nearest));
            ARIADNE_TEST_EQUALS(value(down),Rational(negative ? -row.upper : row.lower));
            ARIADNE_TEST_EQUALS(value(up),Rational(negative ? -row.lower : row.upper));
            ARIADNE_TEST_EQUALS(value(toward_zero),Rational(negative ? -row.lower : row.lower));
        }
    }
    // Common factors preserve exact halfway values while forcing multi-limb division.
    const Integer two65=pow(Integer(2),Nat(65u));
    for (int extra:{1,3}) {
        const Integer numerator=3*(two65+Integer(extra));
        const auto rounded_value=FloatMP::from_ratio(numerator,Integer(3),near,MultiplePrecision(65_bits));
        ARIADNE_TEST_EQUALS(static_cast<Rational>(rounded_value),Rational(two65+Integer(extra==1 ? 0 : 4)));
    }
    // 11/3 * 2^MAX rounds from the largest 2-bit finite significand into infinity.
    auto source=FloatMP::max(MultiplePrecision(4_bits));
    for (unsigned step=0u;step<4u;++step) source=next(down,source);
    const FloatMP denominator(3,precision);
    ARIADNE_TEST_EQUALS(div(down,source,denominator),FloatMP::max(precision));
    ARIADNE_TEST_ASSERT(is_inf(div(near,source,denominator)));

    // The same rounding carry at 65-bit precision uses the general ratio kernel.
    // Scale (2^130-2^65-1) to exponent MAX through supported arithmetic.
    const Integer wide_denominator=two65-Integer(1);
    FloatMP wide_source(wide_denominator*two65-Integer(1),MultiplePrecision(130_bits));
    const auto scale=rec(near,FloatMP::min(MultiplePrecision(2_bits)));
    imul(near,wide_source,scale);
    const FloatMP wide_divisor(wide_denominator,MultiplePrecision(65_bits));
    ARIADNE_TEST_ASSERT(div(down,wide_source,wide_divisor)==FloatMP::max(MultiplePrecision(65_bits)));
    ARIADNE_TEST_ASSERT(is_inf(div(near,wide_source,wide_divisor)));
    ARIADNE_TEST_ASSERT(is_inf(div(down,-wide_source,wide_divisor)));
}

void test_floatmp_negative_power_specials() {
    const MultiplePrecision precision(65_bits);
    ARIADNE_TEST_ASSERT(is_nan(pow(near,FloatMP::nan(precision),Int(-1))));
    for (bool negative:{false,true}) {
        for (Int exponent:{-1,-2}) {
            const bool expected_negative=negative && exponent==-1;
            const auto zero_power=pow(near,FloatMP::zero(65u,negative),exponent);
            const auto infinite_power=pow(near,FloatMP::infinity(65u,negative),exponent);
            ARIADNE_TEST_ASSERT(is_inf(zero_power));
            ARIADNE_TEST_EQUALS(zero_power.negative(),expected_negative);
            ARIADNE_TEST_ASSERT(infinite_power.is_zero());
            ARIADNE_TEST_EQUALS(infinite_power.negative(),expected_negative);
        }
    }
    const auto tiny=FloatMP::min(precision);
    ARIADNE_TEST_EQUALS(floor(-tiny),FloatMP(-1,precision));
    ARIADNE_TEST_EQUALS(ceil(tiny),FloatMP(1,precision));
    ARIADNE_TEST_ASSERT(round(tiny).is_zero());
    ARIADNE_TEST_ASSERT(round(-tiny).is_zero() && round(-tiny).negative());
}

void test_floatmp_decimal_extreme_scales() {
    const MultiplePrecision precision(2_bits);
    const FloatMP small(TwoExp(-2000),precision);
    // Exact integer certificates: 85*2^2000 < 10^604 < 90*2^2000.
    // Thus 2^-2000 lies between 8.5 and 9 units at decimal place 603.
    const std::string prefix="0."+std::string(602u,'0');
    ARIADNE_TEST_EQUALS(print(small,DecimalPrecision(1u),near),prefix+"9");
    ARIADNE_TEST_EQUALS(print(small,DecimalPrecision(1u),down),prefix+"8");
    ARIADNE_TEST_EQUALS(print(small,DecimalPrecision(1u),up),prefix+"9");
    ARIADNE_TEST_EQUALS(print(-small,DecimalPrecision(1u),up),"-"+prefix+"8");
    ARIADNE_TEST_EQUALS(print(-small,DecimalPrecision(1u),down),"-"+prefix+"9");

    // Reuse the independently generated exact decimal oracle for 2^2047.
    std::ifstream data(ARIADNE_TEST_DATA_DIR "/big_unsigned_b132_decimal_structured.tsv");
    if(!data) throw std::runtime_error("Missing structured decimal oracle");
    std::string line,large_expected;
    while(std::getline(data,line)) {
        if(line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        unsigned ordinal;
        std::string source_decimal,source_hex,canonical;
        if(!(fields>>ordinal>>source_decimal>>source_hex>>canonical))
            throw std::runtime_error("Malformed structured decimal oracle");
        if(source_hex=="8"+std::string(511u,'0')) large_expected=canonical+".";
    }
    if(large_expected.empty()) throw std::runtime_error("Missing power-of-two decimal oracle");
    const FloatMP large(TwoExp(2047),precision);
    ARIADNE_TEST_EQUALS(print(large,DecimalPrecision(1u),near),large_expected);
    ARIADNE_TEST_EQUALS(print(-large,DecimalPrecision(1u),near),"-"+large_expected);
    RoundingSession<FloatMP> ambient;
    for(auto ambient_mode:{FloatMP::ROUND_TO_NEAREST,FloatMP::ROUND_DOWNWARD,
                           FloatMP::ROUND_UPWARD,FloatMP::ROUND_TOWARD_ZERO}) {
        ambient.set(ambient_mode);
        ARIADNE_TEST_EQUALS(print(small,DecimalPrecision(1u),near),prefix+"9");
        ARIADNE_TEST_EQUALS(print(-small,DecimalPrecision(1u),near),"-"+prefix+"9");
        ARIADNE_TEST_ASSERT(FloatMP::get_rounding_mode()==ambient_mode);
    }
    // Reject an unrepresentable decimal-place count before attempting storage.
    bool places_overflow=false;
    try { (void)print(FloatMP::min(precision),DecimalPrecision(1u),near); }
    catch(std::overflow_error const&) { places_overflow=true; }
    ARIADNE_TEST_ASSERT(places_overflow);
}

// Scalar boundaries exercise metadata failure paths without allocating giant values.
void test_checked_precision_arithmetic() {
    using Detail::checked_unsigned_add;
    using Detail::checked_unsigned_cast;
    constexpr char message[]="test precision overflow";
    const auto largest=std::numeric_limits<unsigned long>::max();
    ARIADNE_TEST_EQUALS(checked_unsigned_add(0ul,0ul,message),0ul);
    ARIADNE_TEST_EQUALS(checked_unsigned_add(largest,0ul,message),largest);
    ARIADNE_TEST_EQUALS(checked_unsigned_add(0ul,largest,message),largest);
    ARIADNE_TEST_EQUALS(checked_unsigned_add(largest-5ul,5ul,message),largest);
    ARIADNE_TEST_FAIL(checked_unsigned_add(largest-5ul,6ul,message));
    // A narrowing destination models the portability check independently of this ABI.
    ARIADNE_TEST_EQUALS(checked_unsigned_cast<std::uint8_t>(255ul,message),std::uint8_t(255u));
    ARIADNE_TEST_FAIL(checked_unsigned_cast<std::uint8_t>(256ul,message));
    ARIADNE_TEST_EQUALS(checked_unsigned_cast<unsigned long>(std::uint8_t(255u),message),255ul);
    ARIADNE_TEST_ASSERT(Detail::unsigned_fits<unsigned long>(std::uint8_t(255u)));
    const auto check_product_bounds=[&]<class CountType>() {
        const auto limit=std::numeric_limits<CountType>::max();
        ARIADNE_TEST_EQUALS(Detail::compare_unsigned_sum(limit,CountType{1u},limit),1);
        ARIADNE_TEST_EQUALS(Detail::compare_unsigned_sum(limit,CountType{0u},limit),0);
        ARIADNE_TEST_EQUALS(Detail::compare_unsigned_sum(CountType{0u},limit,limit),0);
        ARIADNE_TEST_EQUALS(Detail::compare_unsigned_sum(CountType{0u},CountType{0u},limit),-1);
        ARIADNE_TEST_EQUALS(Detail::compare_unsigned_sum(limit,CountType{0u},CountType{0u}),1);

        ARIADNE_TEST_EQUALS(Detail::checked_unsigned_multiply(limit,CountType{0u},message),CountType{0u});
        ARIADNE_TEST_EQUALS(Detail::checked_unsigned_multiply(CountType{0u},limit,message),CountType{0u});
        ARIADNE_TEST_EQUALS(Detail::checked_unsigned_multiply(limit,CountType{1u},message),limit);
        ARIADNE_TEST_EQUALS(Detail::checked_unsigned_multiply(limit/CountType{2u},CountType{2u},message),limit-CountType{1u});
        ARIADNE_TEST_FAIL(Detail::checked_unsigned_multiply(limit/CountType{2u}+CountType{1u},CountType{2u},message));
    };
    check_product_bounds.template operator()<unsigned long>();
    check_product_bounds.template operator()<std::size_t>();
}

// Public IEEE boundary behavior, using exact constants and elementary pi bounds.
void test_floatmp_transcendental_specials() {
    const MultiplePrecision precision(128_bits);
    using Unary=FloatMP(*)(FloatMP const&);
    const std::array<Unary,5u> bounded_domain={
        +[](FloatMP const& x){return sin(near,x);},
        +[](FloatMP const& x){return cos(near,x);},
        +[](FloatMP const& x){return tan(near,x);},
        +[](FloatMP const& x){return asin(near,x);},
        +[](FloatMP const& x){return acos(near,x);}};
    const auto quiet_nan=FloatMP::nan(precision);
    for(auto operation:bounded_domain) {
        ARIADNE_TEST_ASSERT(is_nan(operation(quiet_nan)));
        ARIADNE_TEST_ASSERT(is_nan(operation(FloatMP::inf(precision))));
        ARIADNE_TEST_ASSERT(is_nan(operation(FloatMP::inf(Sign::NEGATIVE,precision))));
    }
    ARIADNE_TEST_ASSERT(is_nan(atan(near,quiet_nan)));
    const std::array<Unary,4u> odd_functions={
        +[](FloatMP const& x){return sin(near,x);},
        +[](FloatMP const& x){return tan(near,x);},
        +[](FloatMP const& x){return asin(near,x);},
        +[](FloatMP const& x){return atan(near,x);}};
    for(bool negative:{false,true}) {
        const auto zero=FloatMP::zero(128u,negative);
        for(auto operation:odd_functions) {
            const auto zero_image=operation(zero);
            ARIADNE_TEST_ASSERT(zero_image.is_zero());
            ARIADNE_TEST_EQUALS(zero_image.negative(),negative);
            ARIADNE_TEST_EQUALS(zero_image.precision(),precision);
        }
        ARIADNE_TEST_EQUALS(cos(near,zero),FloatMP(1,precision));
        const auto right_angle=acos(near,zero);
        ARIADNE_TEST_ASSERT(right_angle>FloatMP(1,precision));
        ARIADNE_TEST_ASSERT(right_angle<FloatMP(2,precision));
        const auto infinite_argument=FloatMP::infinity(128u,negative);
        const auto infinite_atan=atan(near,infinite_argument);
        ARIADNE_TEST_ASSERT(abs(infinite_atan)>FloatMP(1,precision));
        ARIADNE_TEST_ASSERT(abs(infinite_atan)<FloatMP(2,precision));
        ARIADNE_TEST_EQUALS(infinite_atan.negative(),negative);
    }
    for(int sign:{-1,1}) {
        const auto boundary_asin=asin(near,FloatMP(sign,precision));
        ARIADNE_TEST_ASSERT(abs(boundary_asin)>FloatMP(1,precision));
        ARIADNE_TEST_ASSERT(abs(boundary_asin)<FloatMP(2,precision));
        ARIADNE_TEST_EQUALS(boundary_asin.negative(),sign<0);
        ARIADNE_TEST_ASSERT(is_nan(asin(near,FloatMP(sign*2,precision))));
        ARIADNE_TEST_ASSERT(is_nan(acos(near,FloatMP(sign*2,precision))));
    }
    ARIADNE_TEST_ASSERT(acos(near,FloatMP(1,precision)).is_zero());
    const auto straight_angle=acos(near,FloatMP(-1,precision));
    ARIADNE_TEST_ASSERT(straight_angle>FloatMP(3,precision));
    ARIADNE_TEST_ASSERT(straight_angle<FloatMP(4,precision));
    // An exact positive binary scale also exercises rational conversion
    // with a positive exponent, rather than a fractional significand.
    const auto large_log=log(near,FloatMP(1024,MultiplePrecision(8_bits)));
    ARIADNE_TEST_ASSERT(large_log>FloatMP(6,MultiplePrecision(8_bits)));
    ARIADNE_TEST_ASSERT(large_log<FloatMP(7,MultiplePrecision(8_bits)));
    const FloatMP::Rounding typed_zero{Ariadne::Rounding(toward_zero)};
    ARIADNE_TEST_ASSERT(typed_zero==FloatMP::ROUND_TOWARD_ZERO);
}

// Exact power-of-two inputs force zero coefficients at the first working
// precision. Taylor inequalities determine the final neighboring p-bit
// answer independently of any floating-point implementation or oracle.
void test_floatmp_tiny_transcendentals() {
    const MultiplePrecision precision(128_bits);
    const FloatMP tiny(TwoExp(-256),precision);
    const BigUnsigned leading(std::array<BigUnsigned::Limb,2u>{0u,0x8000000000000000ULL});
    const BigUnsigned below(std::array<BigUnsigned::Limb,2u>{0xffffffffffffffffULL,0xffffffffffffffffULL});
    const BigUnsigned above(std::array<BigUnsigned::Limb,2u>{1u,0x8000000000000000ULL});
    const auto check=[&](FloatMP const& value,BigUnsigned const& significand,
                         FloatMP::Exponent exponent,bool negative) {
        ARIADNE_TEST_ASSERT(value.is_finite() && !value.is_zero());
        ARIADNE_TEST_ASSERT(compare(value.significand(),significand)==0);
        ARIADNE_TEST_EQUALS(value.exponent(),exponent);
        ARIADNE_TEST_EQUALS(value.negative(),negative);
        ARIADNE_TEST_EQUALS(value.precision(),precision);
    };
    using Unary=FloatMP(*)(FloatMP::RoundingModeType,FloatMP const&);
    const std::array<Unary,3u> below_magnitude={
        +[](FloatMP::RoundingModeType mode,FloatMP const& x){return sin(mode,x);},
        +[](FloatMP::RoundingModeType mode,FloatMP const& x){return atan(mode,x);},
        +[](FloatMP::RoundingModeType mode,FloatMP const& x){return tanh(mode,x);}};
    const std::array<Unary,2u> above_magnitude={
        +[](FloatMP::RoundingModeType mode,FloatMP const& x){return tan(mode,x);},
        +[](FloatMP::RoundingModeType mode,FloatMP const& x){return asin(mode,x);}};
    for(bool negative:{false,true}) {
        const FloatMP argument=negative ? -tiny : tiny;
        for(auto mode:{FloatMP::ROUND_TO_NEAREST,FloatMP::ROUND_DOWNWARD,
                       FloatMP::ROUND_UPWARD,FloatMP::ROUND_TOWARD_ZERO}) {
            const bool toward_zero_result=mode==FloatMP::ROUND_TOWARD_ZERO ||
                (mode==FloatMP::ROUND_DOWNWARD && !negative) ||
                (mode==FloatMP::ROUND_UPWARD && negative);
            // sin, atan and tanh have magnitude just below |x|;
            // their error is O(x^3), far below half an ulp at x=2^-256.
            for(auto operation:below_magnitude)
                check(operation(mode,argument),toward_zero_result ? below : leading,
                      toward_zero_result ? -384 : -383,negative);
            const bool away_from_zero=(mode==FloatMP::ROUND_UPWARD && !negative) ||
                (mode==FloatMP::ROUND_DOWNWARD && negative);
            // tan and asin have magnitude just above |x| by O(x^3).
            for(auto operation:above_magnitude)
                check(operation(mode,argument),away_from_zero ? above : leading,-383,negative);
            // cos(x)<1 by O(x^2), and exp(x) is just above/below 1.
            const bool cos_below=mode==FloatMP::ROUND_DOWNWARD || mode==FloatMP::ROUND_TOWARD_ZERO;
            check(cos(mode,argument),cos_below ? below : leading,cos_below ? -128 : -127,false);
            const bool exp_below=negative &&
                (mode==FloatMP::ROUND_DOWNWARD || mode==FloatMP::ROUND_TOWARD_ZERO);
            const bool exp_above=!negative && mode==FloatMP::ROUND_UPWARD;
            check(exp(mode,argument),exp_below ? below : exp_above ? above : leading,
                  exp_below ? -128 : -127,false);
        }
    }
}

void test_floatmp_ieee_conversion_boundaries() {
    const MultiplePrecision precision(128_bits);
    const auto largest=FloatMP::max(precision);
    const auto smallest=FloatMP::min(precision);
    const FloatMP outside_double(TwoExp(2048),precision);
    const FloatMP near_two(Dyadic(Integer(0x3fffffffffffffULL),53u),precision);
    std::array<BigUnsigned::Limb,16u> boundary_words{};
    boundary_words.back()=0xfffffffffffffc00ULL; // 2^1024-2^970, the overflow halfway point.
    const FloatMP overflow_halfway(Integer::from_magnitude(BigUnsigned(boundary_words),false),precision);
    struct Row {
        FloatMP::Rounding mode;
        std::uint64_t overflow_positive,overflow_negative;
        std::uint64_t underflow_positive,underflow_negative;
        std::uint64_t carry_positive,carry_negative;
    };
    const std::array<Row,4u> rows={
        Row{FloatMP::ROUND_TO_NEAREST,0x7ff0000000000000ULL,0xfff0000000000000ULL,
            0u,0x8000000000000000ULL,0x4000000000000000ULL,0xc000000000000000ULL},
        Row{FloatMP::ROUND_DOWNWARD,0x7fefffffffffffffULL,0xfff0000000000000ULL,
            0u,0x8000000000000001ULL,0x3fffffffffffffffULL,0xc000000000000000ULL},
        Row{FloatMP::ROUND_UPWARD,0x7ff0000000000000ULL,0xffefffffffffffffULL,
            1u,0x8000000000000000ULL,0x4000000000000000ULL,0xbfffffffffffffffULL},
        Row{FloatMP::ROUND_TOWARD_ZERO,0x7fefffffffffffffULL,0xffefffffffffffffULL,
            0u,0x8000000000000000ULL,0x3fffffffffffffffULL,0xbfffffffffffffffULL}};
    RoundingSession<FloatMP> ambient;
    for(auto const& row:rows) {
        ambient.set(row.mode);
        for(auto const& value:{largest,outside_double,overflow_halfway}) {
            ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(value.get_d()),row.overflow_positive);
            ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>((-value).get_d()),row.overflow_negative);
        }
        ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(smallest.get_d()),row.underflow_positive);
        ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>((-smallest).get_d()),row.underflow_negative);
        ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(near_two.get_d()),row.carry_positive);
        ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>((-near_two).get_d()),row.carry_negative);
        const FloatMP short_one(1,MultiplePrecision(2_bits));
        ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(short_one.get_d()),0x3ff0000000000000ULL);
    }
}

void test_floatmp_scalar_ratio_limits() {
    const auto largest_word=std::numeric_limits<std::uint64_t>::max();
    const MultiplePrecision precision(64_bits);
    const Rational source{Integer(largest_word-1u),Integer(largest_word)};
    // 2^64*(d-1)/d=2^64-1-1/d, so the floor is 2^64-2;
    // nearest and upward retain 2^64-1. A final increment never wraps.
    for(auto mode:{FloatMP::ROUND_TO_NEAREST,FloatMP::ROUND_DOWNWARD,
                   FloatMP::ROUND_UPWARD,FloatMP::ROUND_TOWARD_ZERO}) {
        const auto image=FloatMP(source,mode,precision);
        const auto expected=(mode==FloatMP::ROUND_DOWNWARD || mode==FloatMP::ROUND_TOWARD_ZERO)
            ? largest_word-1u : largest_word;
        ARIADNE_TEST_EQUALS(image.significand().limb(0u),expected);
        ARIADNE_TEST_EQUALS(image.exponent(),FloatMP::Exponent(-64));
    }
    const MultiplePrecision short_precision(2_bits);
    const Rational short_source(11,3);
    ARIADNE_TEST_ASSERT(FloatMP(short_source,near,short_precision)==FloatMP(4,short_precision));
    ARIADNE_TEST_ASSERT(FloatMP(short_source,up,short_precision)==FloatMP(4,short_precision));
    ARIADNE_TEST_ASSERT(FloatMP(short_source,down,short_precision)==FloatMP(3,short_precision));
}

void test_floatmp_computational_limits() {
    const MultiplePrecision precision(128_bits);
    const auto expect_overflow=[](auto operation) {
        bool raised=false;
        try { (void)operation(); }
        catch(std::overflow_error const&) { raised=true; }
        ARIADNE_TEST_ASSERT(raised);
    };
    // These compact native inputs hit resource guards BEFORE a Taylor
    // series or an enormous exact rational power can be constructed.
    expect_overflow([&]{return exp(near,FloatMP(TwoExp(2050),precision));});
    expect_overflow([&]{return sin(near,FloatMP(TwoExp(4097),precision));});
    expect_overflow([&]{return pow(near,FloatMP::min(precision),Nat(2));});
}

void test_floatmp_exact_series_thresholds() {
    // With p=128 exp_series receives work=p+40=168. This input
    // makes its first fixed-point term EXACTLY the stopping threshold.
    // 0<exp(y)-1<2y and 0<1-exp(-y)<y put both inside half an ulp.
    const MultiplePrecision exp_precision(128_bits);
    const FloatMP exp_input(Dyadic(1,168u),exp_precision);
    const FloatMP exp_one(1,exp_precision);
    for(auto mode:{FloatMP::ROUND_TO_NEAREST,FloatMP::ROUND_DOWNWARD,
                   FloatMP::ROUND_UPWARD,FloatMP::ROUND_TOWARD_ZERO}) {
        const auto positive_expected=mode==FloatMP::ROUND_UPWARD
            ? next(RoundUpward{},exp_one) : exp_one;
        const auto negative_expected=(mode==FloatMP::ROUND_DOWNWARD ||
                                      mode==FloatMP::ROUND_TOWARD_ZERO)
            ? next(RoundDownward{},exp_one) : exp_one;
        ARIADNE_TEST_ASSERT(exp(mode,exp_input)==positive_expected);
        ARIADNE_TEST_ASSERT(exp(mode,neg(exp_input))==negative_expected);
    }
    // p=129 gives the cosine kernel work=161. For y=2^-80,
    // its first omitted term Q*y*y/2 equals 2^guard exactly.
    // 0<1-cos(y)<y*y/2 is far below half the spacing beneath 1.
    const MultiplePrecision cos_precision(129_bits);
    const FloatMP cos_input(Dyadic(1,80u),cos_precision);
    const FloatMP cos_one(1,cos_precision);
    for(auto mode:{FloatMP::ROUND_TO_NEAREST,FloatMP::ROUND_DOWNWARD,
                   FloatMP::ROUND_UPWARD,FloatMP::ROUND_TOWARD_ZERO}) {
        const auto expected=(mode==FloatMP::ROUND_DOWNWARD ||
                             mode==FloatMP::ROUND_TOWARD_ZERO)
            ? next(RoundDownward{},cos_one) : cos_one;
        ARIADNE_TEST_ASSERT(cos(mode,cos_input)==expected);
        ARIADNE_TEST_ASSERT(cos(mode,neg(cos_input))==expected);
    }
}

void test_floatmp_public_edge_paths() {
    const MultiplePrecision precision(128_bits);
    const FloatMP one(1,precision);
    // A short subnormal reaches the general binary64 converter while
    // retaining every source bit: there are no discarded bits to round.
    const FloatMP exact_subnormal(Dyadic(1,1023u),MultiplePrecision(2_bits));
    for(auto mode:{FloatMP::ROUND_TO_NEAREST,FloatMP::ROUND_DOWNWARD,
                   FloatMP::ROUND_UPWARD,FloatMP::ROUND_TOWARD_ZERO}) {
        RoundingSession<FloatMP> rounding_scope;
        rounding_scope.set(mode);
        ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(exact_subnormal.get_d()),
                            std::uint64_t{0x0008000000000000});
    }
    ARIADNE_TEST_ASSERT(ceil(one)==one);
    ARIADNE_TEST_ASSERT(round(FloatMP(Dyadic(9,3u),precision))==one);
    ARIADNE_TEST_ASSERT(abs(near,one)==one);
    ARIADNE_TEST_ASSERT(abs(near,neg(one))==one);

    // Reading a malformed exponent preserves the old value. Uppercase E
    // and nondigit mantissa prefixes are valid, and punctuation is retained.
    FloatMP parsed(7,precision);
    std::istringstream malformed("2e3.4:");
    malformed>>parsed;
    ARIADNE_TEST_ASSERT(malformed.fail());
    ARIADNE_TEST_ASSERT(parsed==FloatMP(7,precision));
    std::istringstream valid("-.5E+2:");
    valid>>parsed;
    ARIADNE_TEST_ASSERT(!valid.fail());
    ARIADNE_TEST_ASSERT(parsed==FloatMP(-50,precision));
    ARIADNE_TEST_EQUALS(valid.peek(),int(':'));

    const MultiplePrecision two_bits(2_bits),three_bits(3_bits);
    const FloatMP broad(8,three_bits),narrow(1,two_bits);
    ARIADNE_TEST_EQUALS(cmp(broad,narrow),Comparison::GREATER);
    ARIADNE_TEST_EQUALS(cmp(narrow,broad),Comparison::LESS);
    ARIADNE_TEST_EQUALS(cmp(FloatMP(1,three_bits),narrow),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(FloatMP(Dyadic(7,2u),three_bits),
                            FloatMP(Dyadic(3,1u),two_bits)),Comparison::GREATER);

    for(bool negative:{false,true}) {
        const auto signed_zero=FloatMP::zero(128u,negative);
        const auto shifted_zero=shft(signed_zero,Int(3));
        ARIADNE_TEST_ASSERT(shifted_zero.is_zero());
        ARIADNE_TEST_EQUALS(shifted_zero.negative(),negative);
        const auto infinite_source=FloatMP::infinity(128u,negative);
        const auto shifted_infinity=shft(infinite_source,Int(-3));
        ARIADNE_TEST_ASSERT(shifted_infinity.is_inf());
        ARIADNE_TEST_EQUALS(shifted_infinity.negative(),negative);
        const auto squared_infinity=pow(near,infinite_source,Nat(2));
        ARIADNE_TEST_ASSERT(squared_infinity.is_inf() && !squared_infinity.negative());
        const auto cubed_infinity=pow(near,infinite_source,Nat(3));
        ARIADNE_TEST_ASSERT(cubed_infinity.is_inf());
        ARIADNE_TEST_EQUALS(cubed_infinity.negative(),negative);
    }
    ARIADNE_TEST_ASSERT(shft(FloatMP::nan(precision),Int(3)).is_nan());

    // 0<y-tanh(y)<y^3, far below half the spacing at this y. Unlike
    // an exact power of two, its predecessor keeps the same binade.
    const FloatMP tiny(Dyadic(3,256u),precision);
    const auto tiny_predecessor=next(RoundDownward{},tiny);
    ARIADNE_TEST_ASSERT(tanh(near,tiny)==tiny);
    ARIADNE_TEST_ASSERT(tanh(down,tiny)==tiny_predecessor);
    ARIADNE_TEST_ASSERT(tanh(up,tiny)==tiny);
    ARIADNE_TEST_ASSERT(tanh(FloatMP::ROUND_TOWARD_ZERO,tiny)==tiny_predecessor);
    ARIADNE_TEST_ASSERT(tanh(near,neg(tiny))==neg(tiny));
    ARIADNE_TEST_ASSERT(tanh(down,neg(tiny))==neg(tiny));
    ARIADNE_TEST_ASSERT(tanh(up,neg(tiny))==neg(tiny_predecessor));
    ARIADNE_TEST_ASSERT(tanh(FloatMP::ROUND_TOWARD_ZERO,neg(tiny))==neg(tiny_predecessor));

    // sum(k=0..3,1/k!)<e<sum(k=0..4,1/k!)+1/100:
    // 8/3<e<1631/600<30/11. Raising the positive bounds to the
    // fourth power proves 48<exp(4)<56. At p=2 the enclosing
    // representable values are 48 and 64, with midpoint 56.
    const FloatMP four(4,two_bits);
    ARIADNE_TEST_ASSERT(exp(near,four)==FloatMP(48,two_bits));
    ARIADNE_TEST_ASSERT(exp(down,four)==FloatMP(48,two_bits));
    ARIADNE_TEST_ASSERT(exp(up,four)==FloatMP(64,two_bits));
    ARIADNE_TEST_ASSERT(exp(FloatMP::ROUND_TOWARD_ZERO,four)==FloatMP(48,two_bits));
}

// Native exponent limits are represented as metadata, never enormous
// Rational inputs. The fixture's exact results come from a standalone oracle.
void test_floatmp_log_native_exponents() {
    std::ifstream data(ARIADNE_TEST_DATA_DIR "/floatmp_log_native.tsv");
    if (!data) throw std::runtime_error("Missing floatmp_log_native.tsv");
    std::set<std::tuple<unsigned,char,std::string>> keys;
    std::string line;
    while(std::getline(data,line)) {
        if(line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        unsigned bits;
        char mode;
        std::string source,numerator,denominator,extra;
        if(!(fields>>bits>>mode>>source>>numerator>>denominator) || (fields>>extra) ||
            (mode!='N' && mode!='D' && mode!='U' && mode!='Z') ||
            (source!="min" && source!="max") ||
            !keys.emplace(bits,mode,source).second) {
            throw std::runtime_error("Malformed floatmp_log_native.tsv");
        }
        const MultiplePrecision precision(static_cast<std::size_t>(bits));
        const auto input=source=="min" ? FloatMP::min(precision) : FloatMP::max(precision);
        const auto rounding=mode=='N'?FloatMP::ROUND_TO_NEAREST:
            mode=='D'?FloatMP::ROUND_DOWNWARD:
            mode=='U'?FloatMP::ROUND_UPWARD:FloatMP::ROUND_TOWARD_ZERO;
        const Rational expected{Integer(String(numerator)),Integer(String(denominator))};
        std::cout<<"log native "<<source<<" precision "<<bits<<" mode "<<mode<<std::endl;
        const auto image=log(rounding,input);
        ARIADNE_TEST_ASSERT(image.is_finite());
        ARIADNE_TEST_EQUALS(image.precision(),precision);
        ARIADNE_TEST_EQUALS(static_cast<Rational>(image),expected);
    }
    if(keys.size()!=32u) throw std::runtime_error("Incomplete native logarithm oracle");
}

void test_floatmp_elementary_specials() {
    const MultiplePrecision precision(128_bits);
    const auto quiet_nan=FloatMP::nan(precision);
    ARIADNE_TEST_ASSERT(is_nan(sqrt(near,quiet_nan)));
    ARIADNE_TEST_ASSERT(is_nan(sqrt(near,FloatMP(-1,precision))));
    ARIADNE_TEST_ASSERT(is_nan(log(near,quiet_nan)));
    ARIADNE_TEST_ASSERT(is_nan(log(near,FloatMP(-1,precision))));
    ARIADNE_TEST_ASSERT(is_nan(exp(near,quiet_nan)));
    using Integral=FloatMP(*)(FloatMP const&);
    const std::array<Integral,3u> integral_operations={
        +[](FloatMP const& x){return floor(x);},
        +[](FloatMP const& x){return ceil(x);},
        +[](FloatMP const& x){return round(x);}};
    for(auto operation:integral_operations) {
        ARIADNE_TEST_ASSERT(is_nan(operation(quiet_nan)));
        const FloatMP scaled_integer(4,MultiplePrecision(2_bits));
        ARIADNE_TEST_ASSERT(operation(scaled_integer)==scaled_integer);
        for(bool negative:{false,true}) {
            const auto infinity_image=operation(FloatMP::infinity(128u,negative));
            ARIADNE_TEST_ASSERT(infinity_image.is_inf());
            ARIADNE_TEST_EQUALS(infinity_image.negative(),negative);
            const auto zero_image=operation(FloatMP::zero(128u,negative));
            ARIADNE_TEST_ASSERT(zero_image.is_zero());
            ARIADNE_TEST_EQUALS(zero_image.negative(),negative);
        }
    }
    for(bool negative:{false,true}) {
        const auto infinite_source=FloatMP::infinity(128u,negative);
        const auto root_image=sqrt(near,infinite_source);
        const auto log_image=log(near,infinite_source);
        const auto exp_image=exp(near,infinite_source);
        if(negative) {
            ARIADNE_TEST_ASSERT(root_image.is_nan());
            ARIADNE_TEST_ASSERT(log_image.is_nan());
            ARIADNE_TEST_ASSERT(exp_image.is_zero() && !exp_image.negative());
        } else {
            ARIADNE_TEST_ASSERT(root_image.is_inf() && !root_image.negative());
            ARIADNE_TEST_ASSERT(log_image.is_inf() && !log_image.negative());
            ARIADNE_TEST_ASSERT(exp_image.is_inf() && !exp_image.negative());
        }
        const auto zero=FloatMP::zero(128u,negative);
        const auto root_zero=sqrt(near,zero);
        ARIADNE_TEST_ASSERT(root_zero.is_zero());
        ARIADNE_TEST_EQUALS(root_zero.negative(),negative);
        const auto log_zero=log(near,zero);
        ARIADNE_TEST_ASSERT(log_zero.is_inf() && log_zero.negative());
        ARIADNE_TEST_ASSERT(exp(near,zero)==FloatMP(1,precision));
    }
}

Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);
    
    test_checked_precision_arithmetic();
    test_floatmp_decimal_extreme_scales();
    test_floatmp_special_arithmetic();
    test_floatmp_neighbor_limits();
    test_floatmp_transcendental_specials();
    test_floatmp_tiny_transcendentals();
    test_floatmp_ieee_conversion_boundaries();
    test_floatmp_elementary_specials();
    test_floatmp_log_native_exponents();
    test_floatmp_public_edge_paths();
    test_floatmp_exact_series_thresholds();
    test_floatmp_scalar_ratio_limits();
    test_floatmp_computational_limits();
    test_floatmp_ratio_ties_and_carry();
    test_floatmp_negative_power_specials();
    test_floatmp_mixed_precision_arithmetic();
    test_floatmp_exponent_limit_oracle();
    test_floatmp_precision_contracts();
    test_floatmp_construction_oracle();
    test_floatmp_minmax_oracle();
    test_floatdp_unary_wrappers();
    test_floatdp_rounding_exception_safety();
    test_floatdp_rounding_mode_contract();
    test_floatdp_sine_outside_certificate();
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

