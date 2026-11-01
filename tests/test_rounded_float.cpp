/***************************************************************************
 *            test_rounded_float.cpp
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
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <map>
#include <fenv.h>

#include "numeric/rounding.hpp"
#include "numeric/floats.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/rounded_float.hpp"
#include "numeric/numeric.hpp"
#include "numeric/dyadic.hpp"


#include "ariadne-testing.hpp"

// The original test's EQUAL/COMPARE checks used decide() on
// ApproximateKleenean, not definitely(). Preserve that explicit
// approximate decision without altering the framework's generic evidence.
namespace {
template<class T>
auto rounded_test_verdict(T const& predicate) {
    if constexpr (std::is_same_v<std::remove_cvref_t<T>, Ariadne::ApproximateKleenean>) {
        using Ariadne::decide;
        return static_cast<bool>(decide(predicate));
    } else {
        return Ariadne::Testing::Detail::evidence(predicate);
    }
}
} // namespace

// Snapshot both operand values while their expression-owned temporaries are
// still alive; evaluate actual then expected once, before comparison/rendering.
#define ARIADNE_ROUNDED_DECIDED_RELATION(actual, relation, expected) \
    this->ctx().check_equal_snapshot([&] { \
        return ([&](auto const& rounded_actual) { \
            return ([&](auto const& rounded_expected) { \
                auto verdict = (rounded_actual relation rounded_expected); \
                return Ariadne::Testing::Detail::EqualityObservation{ \
                    Ariadne::Testing::Detail::evidence(rounded_test_verdict(verdict)), \
                    Ariadne::Testing::Detail::render(rounded_actual), \
                    Ariadne::Testing::Detail::render(rounded_expected)}; \
            }((expected))); \
        }((actual))); \
    }, #actual " " #relation " " #expected)

#define ARIADNE_ROUNDED_DECIDED_EQUAL(actual, expected) \
    ARIADNE_ROUNDED_DECIDED_RELATION(actual, ==, expected)
#define ARIADNE_ROUNDED_DECIDED_NOT_EQUAL(actual, expected) \
    ARIADNE_ROUNDED_DECIDED_RELATION(actual, !=, expected)
#define ARIADNE_ROUNDED_DECIDED_LESS(actual, expected) \
    ARIADNE_ROUNDED_DECIDED_RELATION(actual, <, expected)
#define ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(actual, expected) \
    ARIADNE_ROUNDED_DECIDED_RELATION(actual, <=, expected)
#define ARIADNE_ROUNDED_DECIDED_GREATER(actual, expected) \
    ARIADNE_ROUNDED_DECIDED_RELATION(actual, >, expected)
#define ARIADNE_ROUNDED_DECIDED_GREATER_EQUAL(actual, expected) \
    ARIADNE_ROUNDED_DECIDED_RELATION(actual, >=, expected)

// For old ASSERTs an approximate comparison produced only a warning, never
// a decisive assertion. Here we explicitly verify its approximate decision:
// a false decision fails rather than being silently accepted.
#define ARIADNE_ROUNDED_CHECK_APPROXIMATE(expression) \
    this->ctx().check([&] { \
        decltype(auto) predicate = (expression); \
        return rounded_test_verdict(predicate); \
    }, #expression " [approximate decision if applicable]")



using namespace std;
using namespace Ariadne;

namespace {
// These high-precision bounds are independent certified oracles; the
// FloatDP directed results must strictly enclose them, as before.
struct CertifiedRounding { Rational down; Rational up; };
using CertifiedKey=std::pair<std::string,Rational>;
std::map<CertifiedKey,CertifiedRounding> const& certified_rounding_data() {
#ifndef ARIADNE_TEST_DATA_DIR
#error "Configure ARIADNE_TEST_DATA_DIR"
#endif
    static const auto data=[] {
        std::ifstream input(ARIADNE_TEST_DATA_DIR "/rounded_float_directed.tsv");
        if(!input) throw std::runtime_error("Missing certified rounding TSV");
        std::map<CertifiedKey,CertifiedRounding> rows;
        std::string line;
        while(std::getline(input,line)) {
            if(line.empty() || line[0]=='#') continue;
            std::istringstream row(line);
            std::string kind,n,d,ln,ld,un,ud,extra;
            if(!(row>>kind>>n>>d>>ln>>ld>>un>>ud) || (row>>extra))
                throw std::runtime_error("Malformed certified rounding TSV row");
            Rational x{Integer(String(n)),Integer(String(d))};
            Rational lo{Integer(String(ln)),Integer(String(ld))};
            Rational hi{Integer(String(un)),Integer(String(ud))};
            if(lo>hi || !rows.emplace(CertifiedKey{kind,x},CertifiedRounding{lo,hi}).second)
                throw std::runtime_error("Duplicate or invalid certified rounding row");
        }
        if(rows.size()!=24u) throw std::runtime_error("Wrong certified rounding TSV row count");
        return rows;
    }();
    return data;
}
CertifiedRounding const& certified_rounding(std::string const& kind,Rational const& x) {
    auto const& rows=certified_rounding_data();
    auto it=rows.find(CertifiedKey{kind,x});
    if(it==rows.end()) throw std::runtime_error("Missing certified rounding input");
    return it->second;
}

template<class FLT>
constexpr bool check_rounded_float_concept()
{
    using R=Rounded<FLT>;
    using PR=typename FLT::PrecisionType;
    return requires(PR pr, Int n, Nat m, ExactDouble d, R x, R y) {
        R(pr); R(n,pr); R(m,pr); R(d,pr); R(x,pr);
        x=n; x=m; x=d; x=y;
        max(x,y); min(x,y); abs(x); nul(x); pos(x); neg(x); hlf(x);
        add(x,y); sub(x,y); mul(x,y); div(x,y); fma(x,y,x); pow(x,n); pow(x,m);
        sqr(x); rec(x); sqrt(x); exp(x); log(x); sin(x); cos(x); tan(x); asin(x); acos(x); atan(x);
        +x; -x; x+y; x-y; x*y; x/y; x-=y; x*=y; x/=y;
        R::set_rounding_to_nearest(); R::set_rounding_downward(); R::set_rounding_upward(); R::set_rounding_toward_zero();
        R::get_rounding_mode(); x.precision(); x.raw();
    };
}

static_assert(check_rounded_float_concept<FloatDP>());
static_assert(check_rounded_float_concept<FloatMP>());

class TestNumericCasts : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    void test_numeric_casts();
};

Void TestNumericCasts::test_numeric_casts()
{
    // A heuristic decision is distinct from a proof. Check both outcomes:
    // accepting an approximate 'true' must not also accept 'false'.
    const ApproximateKleenean likely_true(true), likely_false(false);
    ARIADNE_TESTING_EXPECT_EQUAL( rounded_test_verdict(likely_true), true);
    ARIADNE_TESTING_EXPECT_EQUAL( rounded_test_verdict(likely_false), false);
    ARIADNE_TESTING_EXPECT_EQUAL(
        Testing::Detail::evidence(likely_true), Testing::Evidence::Indeterminate);
    ARIADNE_TESTING_EXPECT_EQUAL(
        Testing::Detail::evidence(likely_false), Testing::Evidence::Indeterminate);
    FloatDP dp_two(2,dp);
    const MultiplePrecision mp(128_bits);
    FloatMP ap_two(2,mp);
    Dyadic dy_two(2);
    Real real_two(2);
    FloatDPBounds bounds_two(2,dp);
    FloatDPApproximation approximation_two(2,dp);
    double double_two=2.0;

    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<Int>(dp_two),Int(2));
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<Int>(ap_two),Int(2));
    // The public numeric_cast API must match the legacy backend on
    // positive, negative and non-integral finite inputs.
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<Int>(FloatMP(-2,mp)),Int(-2));
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<Int>(FloatMP(ExactDouble(2.5),mp)),Int(2));
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<Int>(FloatMP(ExactDouble(-2.5),mp)),Int(-2));

    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<double>(dp_two),2.0);
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<double>(dy_two),2.0);
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<double>(real_two),2.0);
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<double>(bounds_two),2.0);
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<double>(approximation_two),2.0);

    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<float>(double_two),2.0f);
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<float>(dp_two),2.0f);
    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<float>(real_two),2.0f);

    ARIADNE_TESTING_EXPECT_EQUAL( numeric_cast<FloatDP>(dp_two),dp_two);

    Real real_from_dp=numeric_cast<Real>(dp_two);
    Real real_from_bounds=numeric_cast<Real>(bounds_two);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE( models,real_from_dp.get(dp),dp_two);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE( models,real_from_bounds.get(dp),dp_two);

    FloatDPBall ball_from_real=numeric_cast<FloatDPBall>(real_two);
    FloatDPBounds bounds_from_real=numeric_cast<FloatDPBounds>(real_two);
    FloatDPApproximation approximation_from_real=numeric_cast<FloatDPApproximation>(real_two);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE( models,ball_from_real,dp_two);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE( models,bounds_from_real,dp_two);
    ARIADNE_TESTING_EXPECT_EQUAL( approximation_from_real.raw(),dp_two);
}

} // namespace


template<class FLT>
class TestRounded : public Ariadne::Testing::ContextFixture
{
    using FloatType = FLT;
    using RoundedFloatType = Rounded<FLT>;
    using PR=typename FLT::PrecisionType;
    PR precision;

  public:
    TestRounded(Ariadne::Testing::Context& context, PR prec);
    Void test_header_api();
    Void test_class();
    Void test_conversion_from_to();
    Void test_comparison();
    Void test_arithmetic();
    Void test_function();
    Void test_sine();
    Void test_cosine();
    Void test_arctan();
};




template<class FLT>
TestRounded<FLT>::TestRounded(Ariadne::Testing::Context& context, PR pr)
    : ContextFixture(context), precision(pr)
{
}



// Test that the type implements all operations of
// the RoundedFloatType concept without testing correctness
template<class FLT> Void
TestRounded<FLT>::test_header_api()
{
    RoundedFloatType::set_rounding_to_nearest();

    RoundedFloatType two_value(2,precision);
    RoundedFloatType four(4,precision);
    RoundedFloatType quotient=four;
    quotient/=two_value;
    ARIADNE_ROUNDED_DECIDED_EQUAL(quotient,two_value);

    if constexpr (Same<FloatType,FloatDP>) {
        auto saved_rounding=RoundedFloatType::get_rounding_mode();
        ARIADNE_TESTING_EXPECT_NO_THROW(RoundedFloatType::set_rounding_toward_zero());
        ARIADNE_TESTING_EXPECT_NO_THROW(RoundedFloatType::set_rounding_mode(saved_rounding));

        FloatDP raw_two(2,dp);
        RoundedFloatType from_raw_precision(raw_two,dp);
        ARIADNE_ROUNDED_DECIDED_EQUAL(from_raw_precision,ExactDouble(2));

        RoundedFloatType assigned(0,dp);
        ARIADNE_TESTING_EXPECT_NO_THROW(assigned=raw_two);
        ARIADNE_ROUNDED_DECIDED_EQUAL(assigned,ExactDouble(2));

        FloatDP raw_copy=static_cast<FloatDP>(assigned);
        ARIADNE_ROUNDED_DECIDED_EQUAL(raw_copy,raw_two);
        ARIADNE_ROUNDED_DECIDED_EQUAL(assigned.precision(),dp);
        ARIADNE_ROUNDED_DECIDED_EQUAL(assigned.characteristics(),dp);

        RoundedFloatType nan_value(FloatDP::nan(dp));
        ARIADNE_ROUNDED_CHECK_APPROXIMATE(is_nan(nan_value));

        RoundedFloatType positive(0.5_x,dp);
        RoundedFloatType negative(-0.5_x,dp);
        RoundedFloatType one(1,dp);

        ARIADNE_TESTING_EXPECT_NO_THROW(round(positive));
        ARIADNE_TESTING_EXPECT_NO_THROW(nul(positive));
        ARIADNE_TESTING_EXPECT_NO_THROW(pos(positive));
        ARIADNE_TESTING_EXPECT_NO_THROW(rec(positive));
        ARIADNE_TESTING_EXPECT_NO_THROW(fma(positive,one,one));
        ARIADNE_TESTING_EXPECT_NO_THROW(tan(positive));
        ARIADNE_TESTING_EXPECT_NO_THROW(tanh(positive));
        ARIADNE_TESTING_EXPECT_NO_THROW(tanh(negative));
        ARIADNE_TESTING_EXPECT_NO_THROW(asin(positive));
        ARIADNE_TESTING_EXPECT_NO_THROW(acos(positive));
        ARIADNE_TESTING_EXPECT_NO_THROW(mag(negative));
        ARIADNE_TESTING_EXPECT_NO_THROW(mig(negative));

        ARIADNE_TESTING_EXPECT_NO_THROW(positive+one);
        ARIADNE_TESTING_EXPECT_NO_THROW(positive*one);

        RoundedFloatType compound(2,dp);
        ARIADNE_TESTING_EXPECT_NO_THROW(compound+=one);
        ARIADNE_TESTING_EXPECT_NO_THROW(compound-=one);
        ARIADNE_TESTING_EXPECT_NO_THROW(compound*=one);

        ExactDouble exact_one(1);
        ARIADNE_TESTING_EXPECT_NO_THROW(positive+exact_one);
        ARIADNE_TESTING_EXPECT_NO_THROW(positive-exact_one);
        ARIADNE_TESTING_EXPECT_NO_THROW(exact_one+positive);
        ARIADNE_TESTING_EXPECT_NO_THROW(exact_one-positive);

        ARIADNE_TESTING_EXPECT_NO_THROW(compound+=exact_one);
        ARIADNE_TESTING_EXPECT_NO_THROW(compound-=exact_one);
        ARIADNE_TESTING_EXPECT_NO_THROW(compound*=exact_one);
        ARIADNE_TESTING_EXPECT_NO_THROW(compound/=exact_one);

        ARIADNE_TESTING_EXPECT_NO_THROW(max(positive,exact_one));
        ARIADNE_TESTING_EXPECT_NO_THROW(min(positive,exact_one));
        ARIADNE_TESTING_EXPECT_NO_THROW(max(exact_one,positive));
        ARIADNE_TESTING_EXPECT_NO_THROW(min(exact_one,positive));

        ARIADNE_ROUNDED_CHECK_APPROXIMATE(same(positive,RoundedFloatType(0.5_x,dp)));
    }
}


template<class FLT> Void
TestRounded<FLT>::test_class()
{
    ARIADNE_TESTING_NOTE(ARIADNE_PRETTY_FUNCTION);
    // Construct from an Int
    RoundedFloatType f1(2,precision);
    ARIADNE_TESTING_INSPECT("constructed f1", f1);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f1,2);
    // Construct from a double
    RoundedFloatType f2(1.25_x,precision);
    ARIADNE_TESTING_INSPECT("constructed f2", f2);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f2,1.25_x);
    // Copy constructor
    RoundedFloatType f3(f2);
    ARIADNE_TESTING_INSPECT("constructed f3", f3);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f3,f2);


    // Assign from an Int
    ARIADNE_TESTING_EXPECT_NO_THROW(f1=3);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f1,3);
    // Assign from a double
    ARIADNE_TESTING_EXPECT_NO_THROW(f2=2.25_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f2,2.25_x);
    // Copy assignment
    ARIADNE_TESTING_EXPECT_NO_THROW(f3=f2);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f3,f2);
    // Self-assignment
    ARIADNE_TESTING_EXPECT_NO_THROW(f3=f3);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f3,f2);

    if(not (precision==FloatType::get_default_precision())) {
        PR low_precision=min(precision,FloatType::get_default_precision());
        PR high_precision=max(precision,FloatType::get_default_precision());
        RoundedFloatType f4(2.25_x,high_precision);
    ARIADNE_TESTING_INSPECT("constructed f4", f4);
        RoundedFloatType f5(3.75_x,low_precision);
    ARIADNE_TESTING_INSPECT("constructed f5", f5);
        ARIADNE_TESTING_EXPECT_NO_THROW(f5=f4);
        ARIADNE_ROUNDED_DECIDED_EQUAL(f5.precision(),high_precision);
        ARIADNE_ROUNDED_DECIDED_EQUAL(f5,f4);
        RoundedFloatType f6(low_precision);
    ARIADNE_TESTING_INSPECT("constructed f6", f6);
        ARIADNE_TESTING_EXPECT_NO_THROW(f6=std::move(f4));
        ARIADNE_ROUNDED_DECIDED_EQUAL(f6.precision(),high_precision);
        ARIADNE_ROUNDED_DECIDED_EQUAL(f6,f5);
    }


}



template<class FLT> Void
TestRounded<FLT>::test_conversion_from_to()
{
    // Convert from dyadic
    RoundedFloatType f1(Dyadic(5,2u),precision);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f1,1.25_x);

    // Convert from self
    RoundedFloatType f2(f1,precision);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f1,f2);
    RoundedFloatType f3(f1,precision);
    ARIADNE_ROUNDED_DECIDED_EQUAL(f1,f3);

    // Convert from integers
    Int n;
    n=std::numeric_limits<Int>::min();
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(RoundedFloatType(n,precision)==n);
    n=std::numeric_limits<Int>::max();
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(RoundedFloatType(n,precision)==n);
    Nat m = std::numeric_limits<Nat>::max();
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(RoundedFloatType(m,precision)==m);

    // Convert from a rational
    Int num=1; Int den=3;
    Rational q(1,3);
    RoundedFloatType::set_rounding_downward();
    ARIADNE_ROUNDED_DECIDED_LESS(Rational(RoundedFloatType(q,precision).raw()), q);
    RoundedFloatType::set_rounding_upward();
    ARIADNE_ROUNDED_DECIDED_GREATER(Rational(RoundedFloatType(Rational(num,den),precision).raw()), q);

    // Convert from a negative rational
    num=-2; den=5;
    RoundedFloatType::set_rounding_downward();
    ARIADNE_ROUNDED_DECIDED_LESS(Rational(RoundedFloatType(q,precision).raw()), q);
    RoundedFloatType::set_rounding_upward();
    ARIADNE_ROUNDED_DECIDED_GREATER(Rational(RoundedFloatType(q,precision).raw()), q);
}

template<class FLT> Void
TestRounded<FLT>::test_comparison()
{
    ARIADNE_TESTING_NOTE(ARIADNE_PRETTY_FUNCTION);

    RoundedFloatType f0(3.00_x,precision);
    ARIADNE_TESTING_INSPECT("constructed f0", f0);
    RoundedFloatType f1(1.25_x,precision);
    ARIADNE_TESTING_INSPECT("constructed f1", f1);
    RoundedFloatType f2(-1.25_x,precision);
    ARIADNE_TESTING_INSPECT("constructed f2", f2);
    RoundedFloatType f3(-2.25_x,precision);
    ARIADNE_TESTING_INSPECT("constructed f3", f3);
    RoundedFloatType f4(1.25_x,precision);
    ARIADNE_TESTING_INSPECT("constructed f4", f4);

    // Test comparison of two equal numbers
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1==f4); ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1!=f4));
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1<=f4); ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1> f4));
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1>=f4); ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1< f4));


    // Test comparison of two different numbers
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1==f2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1!=f2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1<=f2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1> f2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1>=f2); ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1< f2));

    // Test comparison of two negative numbers
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f2==f3)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f2!=f3);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f2<=f3)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f2> f3);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f2>=f3); ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f2< f3));

    // Test comparison with in int
    Int i2=1;
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1==i2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1!=i2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1<=i2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1> i2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1>=i2); ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1< i2));

    Int i1=1;
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(i1==f2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(i1!=f2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(i1<=f2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(i1> f2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(i1>=f2); ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(i1< f2));

    // Test comparison with an exact double
    ExactDouble x2=1.0_x;
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1==x2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1!=x2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1<=x2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1> x2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1>=x2); ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(f1< x2));

    ExactDouble x1=1.0_x;
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(x1==f2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(x1!=f2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(x1<=f2)); ARIADNE_ROUNDED_CHECK_APPROXIMATE(x1> f2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(x1>=f2); ARIADNE_ROUNDED_CHECK_APPROXIMATE(!(x1< f2));

}

template<class FLT> Void
TestRounded<FLT>::test_arithmetic()
{
    ARIADNE_TESTING_NOTE(ARIADNE_PRETTY_FUNCTION);

    static bool full_precision_warning = false;
    RoundedFloatType::set_rounding_to_nearest();
    const RoundedFloatType pi_near=RoundedFloatType::pi(precision);
    const RoundedFloatType four   =RoundedFloatType(4,precision);
    RoundedFloatType::set_rounding_upward();
    if( cast_exact(pi_near * 4.0_x) != cast_exact(mul(pi_near,four)) ) {
        if(not full_precision_warning) {
            full_precision_warning=true;
            this->ctx().warning("Mixed RoundedFloatType/BuiltinTag operations may not use full precision.");
        }
    }

    static const PR pr = precision;

    //static const double eps = std::numeric_limits<double>::epsilon();
    static const RoundedFloatType eps = FloatType::eps(pr);
    static const RoundedFloatType one = FloatType(1,pr);

    // Set next_up some variables
    RoundedFloatType f1(1.25_x,precision); RoundedFloatType f2(2.25_x,precision); RoundedFloatType f3(-3.25_x,precision);
    RoundedFloatType f4(precision); RoundedFloatType f5(precision);

    // Minimum (this should always remain exact)
    f4=min(f1,f2);
    ARIADNE_TESTING_NOTE("min(" << f1 << "," << f2 << ") = " << f4);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f4==f1);
    f4=min(f1,f3);
    ARIADNE_TESTING_NOTE("min(" << f1 << "," << f3 << ") = " << f4);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f4==f3);

    // Maximum (this should always remain exact)
    f4=max(f1,f2);
    ARIADNE_TESTING_NOTE("max(" << f1 << "," << f2 << ") = " << f4);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f4==f2);
    f4=max(f1,f3);
    ARIADNE_TESTING_NOTE("max(" << f1 << "," << f3 << ") = " << f4);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f4==f1);

    // Absolute value (this should always remain exact)
    f4=abs(f1);
    ARIADNE_TESTING_NOTE("abs(" << f1 << ") = " << f4);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f4==1.25_x);
    f5=abs(f3);
    ARIADNE_TESTING_NOTE("abs(" << f3 << ") = " << f5);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f5==3.25_x);

    // Median (this should remain exact here)
    f3=med(f1,f2);
    ARIADNE_TESTING_NOTE(f1 << " <= med(" << f1 << "," << f2 << ")=" << f3 << " <= " << f2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f1<=f3); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f3<=f2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f3==1.75_x);

    // Negation (this should always remain exact)
    f3=neg(f1);
    ARIADNE_TESTING_NOTE("neg(" << f1 << ") = " << f3);
    f3=-f1;
    ARIADNE_TESTING_NOTE("- " << f1 << " = " << f3);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f3==-1.25_x);

    // Addition (this should remain exact here)
    f3=add(f1,f2);
    f4=add(f1,f2);
    ARIADNE_TESTING_NOTE(f3 << " <= " << f1 << " + " << f2 << " <= " << f4);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f3<=3.5_x); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f4>=3.5_x);
    // Addition
    ARIADNE_TESTING_EXPECT_NO_THROW(RoundedFloatType::set_rounding_upward());
    ARIADNE_ROUNDED_DECIDED_GREATER(add(one,eps/2), 1.0_x);
    ARIADNE_TESTING_EXPECT_NO_THROW(RoundedFloatType::set_rounding_downward());
    ARIADNE_ROUNDED_DECIDED_GREATER_EQUAL(add(one,eps/2), 1.0_x);
    if(cast_exact(add(one,eps/2)) != 1.0_x) {
        this->ctx().warning("Results of floating-point operations stored to higher-pr in registers than memory.");
        ARIADNE_ROUNDED_DECIDED_GREATER(add(one,eps/(1<<11)), 1.0_x);
    }
    ARIADNE_ROUNDED_DECIDED_EQUAL(add(one,eps/(1<<12)), RoundedFloatType(1.0_x,pr));
    RoundedFloatType one_add_down_half_epsilon(add(one,eps/2));
    ARIADNE_TESTING_INSPECT("constructed one_add_down_half_epsilon", one_add_down_half_epsilon);
    ARIADNE_ROUNDED_DECIDED_EQUAL(one_add_down_half_epsilon, RoundedFloatType(1.0_x,pr));
    ARIADNE_ROUNDED_DECIDED_EQUAL(one,1.0_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(one, 1.0_x);


    // Subtraction (this should remain exact here)
    RoundedFloatType::set_rounding_downward();
    f3=sub(f1,f2);
    RoundedFloatType::set_rounding_upward();
    f4=sub(f1,f2);
    ARIADNE_TESTING_NOTE(f3 << " <= " << f1 << " - " << f2 << " <= " << f4);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f3<=-1); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f4>=-1);

    // Multiplication (this should remain exact here)
    RoundedFloatType::set_rounding_downward();
    f3=mul(f1,f2);
    RoundedFloatType::set_rounding_upward();
    f4=mul(f1,f2);
    ARIADNE_TESTING_NOTE(f3 << " <= " << f1 << " * " << f2 << " <= " << f4);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f3<=2.8125_x); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f4>=2.8125_x);

    // Division (not exact; should catch errors here)
    RoundedFloatType::set_rounding_downward();
    f3=div(f1,f2);
    RoundedFloatType::set_rounding_upward();
    f4=div(f1,f2);
    RoundedFloatType::set_rounding_to_nearest();
    f5=div(f1,f2);
    ARIADNE_TESTING_NOTE(f3 << " <= " << f1 << " / " << f2 << " <= " << f4);

    RoundedFloatType five(5,precision);
    RoundedFloatType nine(9,precision);

    if (Same<FloatType,FloatDP>) {
        RoundedFloatType expected_five_ninths_up(0.55555555555555558023_pr,pr);
        RoundedFloatType::set_rounding_downward();
        ARIADNE_ROUNDED_DECIDED_LESS(five/nine, expected_five_ninths_up);
        ARIADNE_ROUNDED_DECIDED_LESS(div(five,nine), expected_five_ninths_up);
        RoundedFloatType five_divby_nine_down=five;
        five_divby_nine_down/=nine;
        ARIADNE_ROUNDED_DECIDED_LESS(five_divby_nine_down, expected_five_ninths_up);
    }

    RoundedFloatType::set_rounding_downward();
    RoundedFloatType five_ninths_down = div(five,nine);
    RoundedFloatType five_down = mul(five_ninths_down,nine);
    RoundedFloatType::set_rounding_upward();
    RoundedFloatType five_ninths_up = div(five,nine);
    RoundedFloatType five_up = mul(five_ninths_up,nine);
    ARIADNE_ROUNDED_DECIDED_LESS(five_ninths_down, five_ninths_up);
    ARIADNE_ROUNDED_DECIDED_LESS(five_down, five);
    ARIADNE_ROUNDED_DECIDED_GREATER(five_up, five);


    // Power (not exact; should catch errors here)
    RoundedFloatType::set_rounding_downward();
    f3=pow(f1,3);
    RoundedFloatType::set_rounding_upward();
    f4=pow(f1,3);
    ARIADNE_TESTING_NOTE(f3 << " <= pow(" << f1 << ",3) <= " << f4);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(f3<=1.953125_x); ARIADNE_ROUNDED_CHECK_APPROXIMATE(f4>=1.953125_x);

    RoundedFloatType::set_rounding_downward();
    f3=pow(f1,-2);
    RoundedFloatType::set_rounding_upward();
    f4=pow(f1,-2);
    ARIADNE_TESTING_NOTE(f3 << " <= pow(" << f1 << ",-2) <= " << f4);
    //ARIADNE_ROUNDED_CHECK_APPROXIMATE(Rational(f3)<Rational(16,25));
    //ARIADNE_ROUNDED_CHECK_APPROXIMATE(Rational(f4)>Rational(16,25));

    // Floor and ceiling
    f2=RoundedFloatType(-3.25_x,precision); f3=RoundedFloatType(-2,precision);

    ARIADNE_ROUNDED_CHECK_APPROXIMATE(floor(f1)==1); ARIADNE_ROUNDED_CHECK_APPROXIMATE(ceil(f1)==2);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(floor(f2)==-4); ARIADNE_ROUNDED_CHECK_APPROXIMATE(ceil(f2)==-3);
    ARIADNE_ROUNDED_CHECK_APPROXIMATE(floor(f3)==-2); ARIADNE_ROUNDED_CHECK_APPROXIMATE(ceil(f3)==-2);

    // Check interval conversions
    RoundedFloatType z(0,pr); RoundedFloatType o(1,pr); RoundedFloatType t(3,pr);

    RoundedFloatType::set_rounding_downward();
    RoundedFloatType odtd=div(o,t);
    RoundedFloatType::set_rounding_upward();
    RoundedFloatType odtu=div(o,t);
    ARIADNE_TESTING_NOTE(odtd << " <= 1/3 <= " << odtu);
    ARIADNE_ROUNDED_DECIDED_LESS(odtd, odtu);
    // Regression test to catch errors when RoundedFloatType result is not assigned to a variable

    ARIADNE_ROUNDED_DECIDED_EQUAL(med(RoundedFloatType(2,pr),RoundedFloatType(3,pr)), 2.5_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(rad(RoundedFloatType(2,pr),RoundedFloatType(3,pr)), 0.5_x);
}


template<class FLT> Void
TestRounded<FLT>::test_function()
{
    ARIADNE_TESTING_NOTE(ARIADNE_PRETTY_FUNCTION);

    // Stream precision is configured in main; NOTE diagnostics are independent.

    // Set up some variables
    RoundedFloatType x(precision); RoundedFloatType ra(precision); RoundedFloatType rl(precision); RoundedFloatType ru(precision);

    RoundedFloatType zero=RoundedFloatType(0.0_x,precision);
    RoundedFloatType one=RoundedFloatType(1.0_x,precision);
    RoundedFloatType two_value=RoundedFloatType(2.0_x,precision);
    RoundedFloatType half=RoundedFloatType(0.5_x,precision);

    RoundedFloatType::set_rounding_mode(upward);
    ARIADNE_ROUNDED_DECIDED_EQUAL(exp(zero), one);
    ARIADNE_ROUNDED_DECIDED_EQUAL(log(one), zero);
    ARIADNE_ROUNDED_DECIDED_GREATER(log(exp(one)), one);
    ARIADNE_ROUNDED_DECIDED_GREATER(log(exp(-one)), -one);
    ARIADNE_ROUNDED_DECIDED_GREATER(log(exp(two_value)), two_value);
    ARIADNE_ROUNDED_DECIDED_GREATER(log(exp(-two_value)), -two_value);
    ARIADNE_ROUNDED_DECIDED_GREATER(mul(exp(two_value),exp(-two_value)), one);
    ARIADNE_ROUNDED_DECIDED_GREATER(exp(log(two_value)), two_value);
    ARIADNE_ROUNDED_DECIDED_GREATER(add(log(two_value),log(half)), zero);
    RoundedFloatType::set_rounding_mode(downward);
    ARIADNE_ROUNDED_DECIDED_EQUAL(exp(zero), one);
    ARIADNE_ROUNDED_DECIDED_EQUAL(log(one), zero);
    ARIADNE_ROUNDED_DECIDED_LESS(log(exp(one)), one);
    ARIADNE_ROUNDED_DECIDED_LESS(log(exp(-one)), -one);
    ARIADNE_ROUNDED_DECIDED_LESS(log(exp(two_value)), two_value);
    ARIADNE_ROUNDED_DECIDED_LESS(log(exp(-two_value)), -two_value);
    ARIADNE_ROUNDED_DECIDED_LESS(mul(exp(two_value),exp(-two_value)), one);
    ARIADNE_ROUNDED_DECIDED_LESS(exp(log(two_value)), two_value);
    ARIADNE_ROUNDED_DECIDED_LESS(add(log(two_value),log(half)), zero);

    if constexpr (Same<FLT,FloatDP>) {
        std::vector<ExactDouble> wexp={-2.0_x,-1.5_x,-1.0_x,-0.75_x,-0.5_x,-0.375_x,0.375_x,0.5_x,0.75_x,1.0_x,1.5_x,2.0_x};
        for (auto w : wexp) {
            ARIADNE_TESTING_INSPECT("w",w);
            const CertifiedRounding expected=certified_rounding("exp",Rational(w));
            RoundedFloatType::set_rounding_mode(upward);
            ARIADNE_ROUNDED_DECIDED_GREATER(Rational(cast_exact(exp(RoundedFloatType(w,dp)))), expected.up);
            RoundedFloatType::set_rounding_mode(downward);
            ARIADNE_ROUNDED_DECIDED_LESS(Rational(cast_exact(exp(RoundedFloatType(w,dp)))), expected.down);
        }

        std::vector<ExactDouble> wlog={3.0_x,2.0_x,0.75_x,0.625_x,0.5_x,0.375_x};
        for (auto w : wlog) {
            ARIADNE_TESTING_INSPECT("w",w);
            const CertifiedRounding expected=certified_rounding("log",Rational(w));
            RoundedFloatType::set_rounding_mode(upward);
            ARIADNE_ROUNDED_DECIDED_GREATER(Rational(cast_exact(log(RoundedFloatType(w,dp)))), expected.up);
            RoundedFloatType::set_rounding_mode(downward);
            ARIADNE_ROUNDED_DECIDED_LESS(Rational(cast_exact(log(RoundedFloatType(w,dp)))), expected.down);
        }
    }

}

template<class FLT> Void
TestRounded<FLT>::test_sine()
{
    static const RoundedFloatType one(1.0_x,precision);
    RoundedFloatType::set_rounding_mode(upward);
    RoundedFloatType sin_rnd_up_one=sin(one);
    RoundedFloatType::set_rounding_mode(downward);
    RoundedFloatType sin_rnd_down_one=sin(one);
    RoundedFloatType::set_rounding_mode(to_nearest);
    RoundedFloatType sin_rnd_approx_one=sin(one);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(sin_rnd_down_one, sin_rnd_approx_one);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(sin_rnd_approx_one, sin_rnd_up_one);
    ARIADNE_ROUNDED_DECIDED_LESS(sin_rnd_down_one, sin_rnd_up_one);

}

template<class FLT> Void
TestRounded<FLT>::test_cosine()
{
    //3.14159265358979323846264338327950288419716939937510
    const RoundedFloatType pi_down=RoundedFloatType(FloatType::pi(down,precision));
    const RoundedFloatType pi_up  =RoundedFloatType(FloatType::pi(up,precision));

    static const RoundedFloatType three=RoundedFloatType(3,precision);
    static const RoundedFloatType zero =RoundedFloatType(0,precision);

    RoundedFloatType::set_rounding_mode(downward);
    static const RoundedFloatType third_pi_down=div(pi_down,three);
    RoundedFloatType::set_rounding_mode(upward);
    static const RoundedFloatType third_pi_up  =div(pi_up  ,three);

    RoundedFloatType::set_rounding_mode(upward);
    ARIADNE_ROUNDED_DECIDED_GREATER(cos(hlf(pi_down)), 0.0_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(cos(zero),1.0_x);
    ARIADNE_ROUNDED_DECIDED_GREATER(cos(third_pi_down), 0.5_x);
    ARIADNE_ROUNDED_DECIDED_GREATER(sqr(cos(pi_down/4)), 0.5_x);
    ARIADNE_ROUNDED_DECIDED_GREATER(cos(pi_down/2), 0.0_x);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(cos(pi_up/2), 0.0_x);
    ARIADNE_ROUNDED_DECIDED_GREATER(cos(pi_down), -1.0_x);
    ARIADNE_ROUNDED_DECIDED_GREATER(cos(pi_up), -1.0_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(cos(2*pi_down), 1.0_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(cos(2*pi_up), 1.0_x);
    ARIADNE_ROUNDED_DECIDED_GREATER(cos(3*pi_down), -1.0_x);
    ARIADNE_ROUNDED_DECIDED_GREATER(cos(3*pi_up), -1.0_x);

    RoundedFloatType::set_rounding_mode(downward);
    ARIADNE_ROUNDED_DECIDED_EQUAL(cos(zero),1.0_x);
    ARIADNE_ROUNDED_DECIDED_LESS(cos(third_pi_up), 0.5_x);
    ARIADNE_ROUNDED_DECIDED_LESS(sqr(cos(pi_up/4)), 0.5_x);
    ARIADNE_ROUNDED_DECIDED_GREATER_EQUAL(cos(pi_down/2), 0.0_x);
    ARIADNE_ROUNDED_DECIDED_LESS(cos(pi_up/2), 0.0_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(cos(pi_down), -1.0_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(cos(pi_up), -1.0_x);
    ARIADNE_ROUNDED_DECIDED_LESS(cos(2*pi_down), 1.0_x);
    ARIADNE_ROUNDED_DECIDED_LESS(cos(2*pi_up), 1.0_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(cos(3*pi_down), -1.0_x);
    ARIADNE_ROUNDED_DECIDED_EQUAL(cos(3*pi_up), -1.0_x);
}



template<class FLT> Void
TestRounded<FLT>::test_arctan()
{
    // _pr must round to nearest even regardless of the ambient FP mode.
    if constexpr(std::is_same_v<FloatType,FloatDP>) {
        auto saved_rounding=FloatDP::get_rounding_mode();
        for(auto mode : {FloatDP::ROUND_TO_NEAREST,FloatDP::ROUND_DOWNWARD,
                         FloatDP::ROUND_UPWARD,FloatDP::ROUND_TOWARD_ZERO}) {
            FloatDP::set_rounding_mode(mode);
            ARIADNE_ROUNDED_DECIDED_EQUAL(3.1415926535897931_pr,0x1.921fb54442d18p+1_x);
            ARIADNE_ROUNDED_DECIDED_EQUAL(FloatDP::get_rounding_mode(),mode);
        }
        FloatDP::set_rounding_mode(saved_rounding);
    }

    //pi~=3.14159265358979323846264338327950288419716939937510
    const RoundedFloatType pi_down(FloatType::pi(down,precision));
    const RoundedFloatType pi_near(FloatType::pi(near,precision));
    const RoundedFloatType pi_up  (FloatType::pi(up,precision));

    if constexpr(std::is_same<FloatType,FloatDP>::value) {
        ARIADNE_ROUNDED_DECIDED_EQUAL(pi_down.raw(),3.1415926535897931_pr);
        ARIADNE_ROUNDED_DECIDED_EQUAL(pi_near.raw(),3.1415926535897931_pr);
        ARIADNE_ROUNDED_DECIDED_EQUAL(pi_up.raw(),  3.1415926535897936_pr);
    }

    RoundedFloatType three=RoundedFloatType(3,precision);
    RoundedFloatType::set_rounding_mode(downward);
    const RoundedFloatType sqrt_three_down=sqrt(three);
    RoundedFloatType::set_rounding_mode(upward);
    const RoundedFloatType sqrt_three_up  =sqrt(three);

    const RoundedFloatType zero(0,precision);
    const RoundedFloatType one(1,precision);
    const RoundedFloatType four(4,precision);
    const RoundedFloatType six(6,precision);
    const RoundedFloatType eps(FloatType::eps(precision));

    RoundedFloatType::set_rounding_mode(upward);
    ARIADNE_ROUNDED_DECIDED_EQUAL(atan(zero),zero);
    ARIADNE_ROUNDED_DECIDED_GREATER_EQUAL(mul(atan(+one),four), pi_up);
    ARIADNE_ROUNDED_DECIDED_GREATER_EQUAL(mul(atan(-one),four), -pi_down);
    ARIADNE_ROUNDED_DECIDED_GREATER_EQUAL(mul(atan(+sqrt_three_up),three), pi_up);
    ARIADNE_ROUNDED_DECIDED_GREATER_EQUAL(mul(atan(-sqrt_three_down),three), -pi_down);
    ARIADNE_ROUNDED_DECIDED_GREATER_EQUAL(mul(atan(+1/sqrt_three_down),six), pi_up);
    ARIADNE_ROUNDED_DECIDED_GREATER_EQUAL(mul(atan(-1/sqrt_three_up),six), -pi_down);

    RoundedFloatType::set_rounding_mode(downward);
    ARIADNE_ROUNDED_DECIDED_EQUAL(atan(zero),zero);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(mul(atan(+one),four), pi_down);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(mul(atan(-one),four), -pi_up);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(mul(atan(+sqrt_three_down),three), pi_down);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(mul(atan(-sqrt_three_up),three), -pi_up);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(mul(atan(+1/sqrt_three_up),six), pi_down);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(mul(atan(-1/sqrt_three_down),six), -pi_up);

    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(sub(mul(atan(one),four),pi_up), 4*eps);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(sub(mul(atan(sqrt_three_up),three),pi_up), 4*eps);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(sub(mul(atan(1/sqrt_three_down),six),pi_up), 4*eps);

    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(atan(one)*4-pi_up, 4*eps);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(atan(sqrt_three_up)*3-pi_up, 4*eps);
    ARIADNE_ROUNDED_DECIDED_LESS_EQUAL(atan(1/sqrt_three_down)*6-pi_up, 4*eps);

    // Regression test
    ARIADNE_ROUNDED_DECIDED_LESS(pi_down/4, pi_up/4);


    if constexpr (Same<FLT,FloatDP>) {
        Dyadic e(1,48u);

        // Same six boundary inputs, now tested against certified directed results.
        std::vector<Dyadic> tests({1+e,1,1-e,-1+e,-1,-1-e});
        for (Dyadic w : tests) {
            ARIADNE_TESTING_INSPECT("w",w);
            const CertifiedRounding expected=certified_rounding("atan",Rational(w));
            RoundedFloatType::set_rounding_mode(upward);
            ARIADNE_ROUNDED_DECIDED_GREATER(Rational(cast_exact(atan(RoundedFloatType(w,dp)))), expected.up);
            RoundedFloatType::set_rounding_mode(downward);
            ARIADNE_ROUNDED_DECIDED_LESS(Rational(cast_exact(atan(RoundedFloatType(w,dp)))), expected.down);
        }
    }
}


template<class FLT>
void add_rounded_cases(Ariadne::Testing::Suite& suite,
                       std::string const& variant,
                       typename FLT::PrecisionType precision) {
    using Fixture=TestRounded<FLT>;
    auto add=[&](std::string const& name, auto method) {
        suite.add(variant+" / "+name, [precision,method](Ariadne::Testing::Context& context) {
            Fixture fixture(context,precision);
            (fixture.*method)();
        });
    };
    add("header api", &Fixture::test_header_api);
    add("class", &Fixture::test_class);
    add("conversion from to", &Fixture::test_conversion_from_to);
    add("comparison", &Fixture::test_comparison);
    add("arithmetic", &Fixture::test_arithmetic);
    add("function", &Fixture::test_function);
    add("sine", &Fixture::test_sine);
    add("cosine", &Fixture::test_cosine);
    add("arctan", &Fixture::test_arctan);
}

int main(int argc, char** argv) {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);
    Ariadne::Testing::Suite suite("rounded floats");
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestNumericCasts, test_numeric_casts);
    add_rounded_cases<FloatDP>(suite,"FloatDP",dp);
    add_rounded_cases<FloatMP>(suite,"FloatMP",MultiplePrecision(128_bits));
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
