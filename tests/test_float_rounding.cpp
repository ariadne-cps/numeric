/***************************************************************************
 *            test_float_rounding.cpp
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

#include <bit>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <fenv.h>

#include "numeric/rounding.hpp"
#include "numeric/floats.hpp"
#include "numeric/numeric.hpp"
#include "numeric/dyadic.hpp"

#include "numeric/floatmp.hpp"

#include "macro/testing.hpp"

using namespace std;
using namespace Ariadne;

namespace {

template<class PR>
constexpr bool check_float_rounding_concept()
{
    using F=RawFloat<PR>;
    return requires(PR pr, Int n, Nat m, ExactDouble d, F x, F y) {
        F(pr); F(n,pr); F(m,pr); F(d,pr);
        x=n; x=m; x=d; x=y;
        max(x,y); min(x,y); abs(x); nul(x); pos(x); neg(x); hlf(x);
        add(near,x,y); add(down,x,y); add(up,x,y);
        sub(near,x,y); sub(down,x,y); sub(up,x,y);
        mul(near,x,y); mul(down,x,y); mul(up,x,y);
        div(near,x,y); div(down,x,y); div(up,x,y);
        sqr(approx,x); rec(approx,x); sqrt(approx,x); exp(approx,x); log(approx,x);
        sin(approx,x); cos(approx,x); tan(approx,x); tanh(approx,x); asin(approx,x); acos(approx,x); atan(approx,x);
        pow(approx,x,n); pow(approx,x,m);
        +x; -x; x+y; x-y; x*y; x/y;
        F::set_rounding_to_nearest(); F::set_rounding_downward(); F::set_rounding_upward(); F::set_rounding_toward_zero();
        F::get_rounding_mode(); x.precision();
    };
}

static_assert(check_float_rounding_concept<DoublePrecision>());
static_assert(check_float_rounding_concept<MultiplePrecision>());

#ifndef ARIADNE_TEST_DATA_DIR
#error "Provide the absolute test data directory from CMake"
#endif

inline constexpr char ARIADNE_FLOATMP_ROUNDING_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_rounding.tsv";
inline constexpr char ARIADNE_FLOATMP_ADD_SUB_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_add_sub.tsv";
inline constexpr char ARIADNE_FLOATMP_EXP_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_exp.tsv";
inline constexpr char ARIADNE_FLOATMP_POW_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_pow.tsv";
inline constexpr char ARIADNE_FLOATMP_SQRT_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_sqrt.tsv";
inline constexpr char ARIADNE_FLOATMP_LOG_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_log.tsv";
inline constexpr char ARIADNE_FLOATMP_CENTERED_EXP_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_centered_exp.tsv";
inline constexpr char ARIADNE_FLOATMP_ASIN_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_asin.tsv";
inline constexpr char ARIADNE_FLOATMP_ACOS_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_acos.tsv";
inline constexpr char ARIADNE_FLOATMP_ATAN_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_atan.tsv";
inline constexpr char ARIADNE_FLOATMP_SIN_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_sin.tsv";
inline constexpr char ARIADNE_FLOATMP_COS_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_cos.tsv";
inline constexpr char ARIADNE_FLOATMP_TAN_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_tan.tsv";
inline constexpr char ARIADNE_FLOATMP_GET_D_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/floatmp_get_d.tsv";
// Exact references, keyed by precision, mode, and the rational input.
// Ordered lookup avoids quadratic duplicate checking and linear scans as
// the frozen oracle dataset grows. The comparator uses exact Rational order.
struct FloatMPGoldenKey {
    std::size_t bits;
    char mode;
    Rational input;
};
struct FloatMPGoldenKeyLess {
    bool operator()(FloatMPGoldenKey const& a,
                    FloatMPGoldenKey const& b) const {
        if (a.bits != b.bits) return a.bits < b.bits;
        if (a.mode != b.mode) return a.mode < b.mode;
        return static_cast<bool>(a.input < b.input);
    }
};
using FloatMPGoldenTable = std::map<FloatMPGoldenKey,Rational,FloatMPGoldenKeyLess>;
FloatMPGoldenTable load_floatmp_golden_rows(char const* path,
                                            std::size_t expected_rows) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error(
        std::string("Missing FloatMP oracle: ")+path);
    FloatMPGoldenTable rows;
    std::string line;
    unsigned long line_number=0ul;
    while (std::getline(file,line)) {
        ++line_number;
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        std::size_t bits=0u;
        std::string tag,in_num,in_den,out_num,out_den,extra;
        if (!(fields>>bits>>tag>>in_num>>in_den>>out_num>>out_den)
            || (fields>>extra)) {
            throw std::runtime_error("Malformed FloatMP oracle line "
                                     +std::to_string(line_number));
        }
        if (bits<2u || tag.size()!=1u ||
            (tag[0]!='N' && tag[0]!='D' && tag[0]!='U' && tag[0]!='Z')) {
            throw std::runtime_error("Invalid FloatMP oracle mode/precision at line "
                                     +std::to_string(line_number));
        }
        const Rational input{Integer(String(in_num)),Integer(String(in_den))};
        const Rational expected{Integer(String(out_num)),Integer(String(out_den))};
        if (!rows.emplace(FloatMPGoldenKey{bits,tag[0],input},
                          expected).second) {
            throw std::runtime_error("Duplicate FloatMP oracle line "
                                     +std::to_string(line_number));
        }
    }
    if (rows.size()!=expected_rows) throw std::runtime_error(
        std::string("Unexpected FloatMP oracle row count in ")+path+
        ": "+std::to_string(rows.size())+
        " (expected "+std::to_string(expected_rows)+")");
    return rows;
}
Rational const& find_floatmp_golden(FloatMPGoldenTable const& rows,
                                    std::size_t bits,
                                    FloatMP::RoundingModeType mode,
                                    Rational const& input) {
    char tag='?';
    if (mode==FloatMP::ROUND_TO_NEAREST) tag='N';
    else if (mode==FloatMP::ROUND_DOWNWARD) tag='D';
    else if (mode==FloatMP::ROUND_UPWARD) tag='U';
    else if (mode==FloatMP::ROUND_TOWARD_ZERO) tag='Z';
    else throw std::invalid_argument("Unknown FloatMP rounding mode");
    auto const found=rows.find(FloatMPGoldenKey{bits,tag,input});
    if (found!=rows.end()) return found->second;
    throw std::runtime_error("Missing FloatMP oracle: bits="
                             +std::to_string(bits)+" mode="+tag);
}
Rational const& floatmp_golden(std::size_t bits,
                               FloatMP::RoundingModeType mode,
                               Rational const& input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_ROUNDING_ORACLE_PATH,4794u);
    return find_floatmp_golden(rows,bits,mode,input);
}
Rational const& floatmp_add_sub_golden(std::size_t bits,
                                       FloatMP::RoundingModeType mode,
                                       Rational const& input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_ADD_SUB_ORACLE_PATH,1780u);
    return find_floatmp_golden(rows,bits,mode,input);
}
Rational const& floatmp_exp_golden(std::size_t bits,
                                   FloatMP::RoundingModeType mode,
                                   Rational const& input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_EXP_ORACLE_PATH,320u);
    return find_floatmp_golden(rows,bits,mode,input);
}
Rational const& floatmp_sqrt_golden(std::size_t bits,
                                    FloatMP::RoundingModeType mode,
                                    Rational const& input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_SQRT_ORACLE_PATH,300u);
    return find_floatmp_golden(rows,bits,mode,input);
}
Rational const& floatmp_log_golden(std::size_t bits,
                                   FloatMP::RoundingModeType mode,
                                   Rational const& input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_LOG_ORACLE_PATH,452u);
    return find_floatmp_golden(rows,bits,mode,input);
}
 
Rational const& floatmp_centered_exp_golden(std::size_t bits,
                                            FloatMP::RoundingModeType mode,
                                            Rational const& input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_CENTERED_EXP_ORACLE_PATH,568u);
    return find_floatmp_golden(rows,bits,mode,input);
}

// Inverse-trigonometric references are keyed by the ORIGINAL Rational.
// Nearest/even conversion of that source was computed offline, independent
// of FloatMP, and cannot silently follow a broken FloatMP constructor.
Rational const& floatmp_asin_golden(std::size_t bits,
                                    FloatMP::RoundingModeType mode,
                                    Rational const& original_input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_ASIN_ORACLE_PATH,728u);
    return find_floatmp_golden(rows,bits,mode,original_input);
}
Rational const& floatmp_acos_golden(std::size_t bits,
                                    FloatMP::RoundingModeType mode,
                                    Rational const& original_input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_ACOS_ORACLE_PATH,728u);
    return find_floatmp_golden(rows,bits,mode,original_input);
}

Rational const& floatmp_atan_golden(std::size_t bits,
                                    FloatMP::RoundingModeType mode,
                                    Rational const& original_input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_ATAN_ORACLE_PATH,1340u);
    return find_floatmp_golden(rows,bits,mode,original_input);
}

// Trigonometric fixtures use the original signed Rational. Offline
// input conversion to a p-bit binary value was independently rounded
// nearest/even, so no FloatMP constructor error can shift the key.
Rational const& floatmp_sin_golden(std::size_t bits,
                                   FloatMP::RoundingModeType mode,
                                   Rational const& original_input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_SIN_ORACLE_PATH,1080u);
    return find_floatmp_golden(rows,bits,mode,original_input);
}
Rational const& floatmp_cos_golden(std::size_t bits,
                                   FloatMP::RoundingModeType mode,
                                   Rational const& original_input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_COS_ORACLE_PATH,1080u);
    return find_floatmp_golden(rows,bits,mode,original_input);
}

Rational const& floatmp_tan_golden(std::size_t bits,
                                   FloatMP::RoundingModeType mode,
                                   Rational const& original_input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_TAN_ORACLE_PATH,1208u);
    return find_floatmp_golden(rows,bits,mode,original_input);
}

// Expected values are unsigned binary64 bit-patterns stored as exact
// integral Rationals (denominator 1) in the existing golden file format.
// Original source keys are intentionally independent of FloatMP::get_d
// and of any FloatMP constructor's result.
Rational const& floatmp_get_d_bits_golden(std::size_t bits,
                                          FloatMP::RoundingModeType mode,
                                          Rational const& original_input) {
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_FLOATMP_GET_D_ORACLE_PATH,272u);
    return find_floatmp_golden(rows,bits,mode,original_input);
}

// Integer power lookup uses the original exact Rational and exponent
// as its key. The offline data performs the nearest-even input rounding
// independently, so the test does not derive its oracle from FloatMP.
struct FloatMPPowerGoldenKey {
    std::size_t bits;
    char mode;
    Rational input;
    std::size_t degree;
};
struct FloatMPPowerGoldenKeyLess {
    bool operator()(FloatMPPowerGoldenKey const& a,
                    FloatMPPowerGoldenKey const& b) const {
        if (a.bits!=b.bits) return a.bits<b.bits;
        if (a.mode!=b.mode) return a.mode<b.mode;
        if (a.input!=b.input) return static_cast<bool>(a.input<b.input);
        return a.degree<b.degree;
    }
};
using FloatMPPowerGoldenTable =
    std::map<FloatMPPowerGoldenKey,Rational,FloatMPPowerGoldenKeyLess>;

FloatMPPowerGoldenTable load_floatmp_power_golden_rows(
        char const* path=ARIADNE_FLOATMP_POW_ORACLE_PATH,
        std::size_t expected_rows=2400u) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error(
        std::string("Missing FloatMP integer power oracle: ")+path);
    FloatMPPowerGoldenTable rows;
    std::string line;
    std::size_t line_number=0u;
    while (std::getline(file,line)) {
        ++line_number;
        if (line.empty() || line[0]=='#') continue;
        std::istringstream fields(line);
        std::size_t bits=0u;
        std::size_t degree=0u;
        std::string mode,in_num,in_den,out_num,out_den,extra;
        if (!(fields>>bits>>mode>>in_num>>in_den>>degree>>out_num>>out_den)
            || (fields>>extra)) {
            throw std::runtime_error("Malformed FloatMP power fixture line "
                                     +std::to_string(line_number));
        }
        if (bits<2u || mode.size()!=1u ||
            (mode[0]!='N' && mode[0]!='D' &&
             mode[0]!='U' && mode[0]!='Z')) {
            throw std::runtime_error("Invalid FloatMP power fixture metadata at line "
                                     +std::to_string(line_number));
        }
        const Rational input{Integer(String(in_num)),Integer(String(in_den))};
        const Rational expected{Integer(String(out_num)),Integer(String(out_den))};
        if (!rows.emplace(FloatMPPowerGoldenKey{bits,mode[0],input,degree},
                          expected).second) {
            throw std::runtime_error("Duplicate FloatMP power fixture line "
                                     +std::to_string(line_number));
        }
    }
    if (rows.size()!=expected_rows) throw std::runtime_error(
        "Unexpected FloatMP integer power oracle row count: "
        +std::to_string(rows.size()));
    return rows;
}
Rational const& floatmp_power_golden(std::size_t bits,
                                     FloatMP::RoundingModeType mode,
                                     Rational const& input,
                                     Nat degree) {
    static const auto rows=load_floatmp_power_golden_rows();
    char tag='?';
    if (mode==FloatMP::ROUND_TO_NEAREST) tag='N';
    else if (mode==FloatMP::ROUND_DOWNWARD) tag='D';
    else if (mode==FloatMP::ROUND_UPWARD) tag='U';
    else if (mode==FloatMP::ROUND_TOWARD_ZERO) tag='Z';
    else throw std::invalid_argument("Unknown FloatMP integer power rounding mode");
    auto const found=rows.find(FloatMPPowerGoldenKey{
        bits,tag,input,static_cast<std::size_t>(degree)});
    if (found!=rows.end()) return found->second;
    throw std::runtime_error("Missing FloatMP power fixture: bits="
                             +std::to_string(bits)+" mode="+tag+
                             " degree="+std::to_string(degree));
}



// Diagnostic-only wall-clock checkpoints. They do not alter test inputs,
// assertions, execution order, or the arithmetic implementation.
class FloatMPPhaseTiming {
    const char* _name;
    std::size_t _enclosing_bits;
    std::chrono::steady_clock::time_point _begin;
  public:
    FloatMPPhaseTiming(const char* name,std::size_t enclosing_bits)
      : _name(name),_enclosing_bits(enclosing_bits),
        _begin(std::chrono::steady_clock::now()) {
        std::cerr << "[floatmp-timing] BEGIN " << _name
                  << " suite_bits=" << _enclosing_bits << '\n';
    }
    ~FloatMPPhaseTiming() {
        const auto elapsed=std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-_begin).count();
        std::cerr << "[floatmp-timing] END " << _name
                  << " suite_bits=" << _enclosing_bits
                  << " elapsed_ms=" << elapsed << '\n';
    }
    FloatMPPhaseTiming(FloatMPPhaseTiming const&)=delete;
    FloatMPPhaseTiming& operator=(FloatMPPhaseTiming const&)=delete;
};


} // namespace


template<class PR>
class TestFloatRounding
{
    using Float = RawFloat<PR>;
    PR precision;

  public:
    TestFloatRounding(PR prec);
    Void test();
  private:
    Void test_class();
    Void test_limits();
    Void test_conversion_from_to();
    Void test_conversion_between();
    Void test_stream();
    Void test_comparison();
    Void test_rounding_mode();
    Void test_rounding_session();
    Void test_double_rounding();
    Void test_arithmetic();
    Void test_cosine();
    Void test_arctan();
    Void test_precision_regressions();
    Void test_short_power_and_tanh();
    Void test_function();
};


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    ARIADNE_TEST_CLASS(FloatDP,TestFloatRounding<DoublePrecision>(dp));
    ARIADNE_TEST_CLASS(FloatMP,TestFloatRounding<MultiplePrecision>(MultiplePrecision(64_bits)));
    ARIADNE_TEST_CLASS(FloatMP,TestFloatRounding<MultiplePrecision>(MultiplePrecision(192_bits)));


    return ARIADNE_TEST_FAILURES;
}


template<class PR>
TestFloatRounding<PR>::TestFloatRounding(PR pr)
    : precision(pr)
{
}

template<class PR> Void
TestFloatRounding<PR>::test()
{
    ARIADNE_TEST_CALL(test_rounding_mode());
    ARIADNE_TEST_CALL(test_rounding_session());
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_limits());
    ARIADNE_TEST_CALL(test_conversion_from_to());
    ARIADNE_TEST_CALL(test_conversion_between());
    ARIADNE_TEST_CALL(test_stream());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_double_rounding());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_function());
    ARIADNE_TEST_CALL(test_cosine());
    ARIADNE_TEST_CALL(test_arctan());
    ARIADNE_TEST_CALL(test_precision_regressions());
    ARIADNE_TEST_CALL(test_short_power_and_tanh());
}


// Test that the type implements all operations of
// the Float concept without testing correctness
template<class PR> Void
TestFloatRounding<PR>::test_class()
{
    cout << ARIADNE_PRETTY_FUNCTION << endl;
    // Construct from an Int
    ARIADNE_TEST_CONSTRUCT(Float,f1,(2,precision));
    ARIADNE_TEST_EQUALS(f1,2);
    // Construct from an ExactDouble
    ARIADNE_TEST_CONSTRUCT(Float,f2,(1.25_x,precision));
    ARIADNE_TEST_EQUALS(f2,1.25_x);
    // Copy constructor
    ARIADNE_TEST_CONSTRUCT(Float,f3,(f2));
    ARIADNE_TEST_EQUAL(f3,f2);


    // Assign from an Int
    ARIADNE_TEST_EXECUTE(f1=3);
    ARIADNE_TEST_EQUALS(f1,3);
    // Assign from a double
    ARIADNE_TEST_EXECUTE(f2=2.25_x);
    ARIADNE_TEST_EQUALS(f2,2.25_x);
    // Copy assignment
    ARIADNE_TEST_EXECUTE(f3=f2);
    ARIADNE_TEST_EQUAL(f3,f2);
    // Self-assignment
    ARIADNE_TEST_EXECUTE(f3=f3);
    ARIADNE_TEST_EQUAL(f3,f2);

    if(not (precision==Float::get_default_precision())) {
        PR low_precision=min(precision,Float::get_default_precision());
        PR high_precision=max(precision,Float::get_default_precision());
        ARIADNE_TEST_CONSTRUCT(Float,f4,(2.25_x,high_precision));
        ARIADNE_TEST_CONSTRUCT(Float,f5,(3.75_x,low_precision));
        ARIADNE_TEST_EXECUTE(f5=f4);
        ARIADNE_TEST_EQUAL(f5.precision(),high_precision);
        ARIADNE_TEST_EQUAL(f5,f4);
        ARIADNE_TEST_CONSTRUCT(Float,f6,(low_precision));
        ARIADNE_TEST_EXECUTE(f6=std::move(f4));
        ARIADNE_TEST_EQUAL(f6.precision(),high_precision);
        ARIADNE_TEST_EQUAL(f6,f5);
    }


}


template<class PR> Void
TestFloatRounding<PR>::test_limits()
{
    ExactDouble pinf_d = cast_exact(std::numeric_limits<double>::infinity());
    ExactDouble ninf_d = cast_exact(-std::numeric_limits<double>::infinity());
    ExactDouble nan_d = cast_exact(std::numeric_limits<double>::quiet_NaN());

    Float zero=Float(0,precision);
    Float one=Float(1,precision);
    Float two_value=Float(2,precision);

    Float min=Float::min(precision);
    Float eps=Float::eps(precision);
    Float max=Float::max(precision);

    Float inf_=Float::inf(precision);
    Float pinf=Float::inf(Sign::POSITIVE,precision);
    Float ninf=Float::inf(Sign::NEGATIVE,precision);
    Float zinf=Float::inf(Sign::ZERO,precision);
    Float nan=Float::nan(precision);

    ARIADNE_TEST_PRINT(eps);
    ARIADNE_TEST_COMPARE(add(down,one,eps),>,one);
    ARIADNE_TEST_COMPARE(add(down,one,div(up,eps,two_value)),==,one);

    ARIADNE_TEST_ASSERT(min>zero);
    ARIADNE_TEST_ASSERT(max<inf_);

    ARIADNE_TEST_ASSERT(add(down,inf_,min)==inf_);
    ARIADNE_TEST_ASSERT(add(up,max,min)==inf_);

    if constexpr (Same<PR,MultiplePrecision>) {
        // The MP exponent gap spans the entire signed 64-bit range.
        // Directed rounding must respond to the tiny operand's sign.
        ARIADNE_TEST_ASSERT(add(near,max,min)==max);
        ARIADNE_TEST_ASSERT(add(down,max,min)==max);
        ARIADNE_TEST_ASSERT(add(up,max,-min)==max);
        ARIADNE_TEST_ASSERT(add(down,max,-min)==next(down,max));
        ARIADNE_TEST_ASSERT(add(down,-max,-min)==-inf_);
        ARIADNE_TEST_ASSERT(add(up,-max,-min)==-max);
        ARIADNE_TEST_ASSERT(add(up,-max,min)==next(up,-max));
        // In-place addition retains the destination precision even when
        // the dominant source operand has fewer significant bits.
        const MultiplePrecision wider(precision.bits()*2u);
        FloatMP wide_tiny=FloatMP::min(wider);
        const FloatMP wide_high(max,wider);
        // Rescaling FloatMP at INT64_MAX must not construct a gigantic Rational.
        ARIADNE_TEST_ASSERT(wide_high.precision()==wider);
        ARIADNE_TEST_ASSERT(wide_high==max);
        ARIADNE_TEST_ASSERT(FloatMP(wide_high,precision)==max);
        ARIADNE_TEST_ASSERT(FloatMP(-max,wider)==-wide_high);
        ARIADNE_TEST_ASSERT(iadd(up,wide_tiny,max)==next(up,wide_high));
    }

    ARIADNE_TEST_ASSERT(is_nan(zinf));
    ARIADNE_TEST_EQUAL(pinf,inf_);
    ARIADNE_TEST_EQUAL(ninf,-inf_);

    ARIADNE_TEST_ASSERT(is_nan(nan));
    ARIADNE_TEST_ASSERT(not(nan==nan));
    ARIADNE_TEST_ASSERT(not(nan==one));
    ARIADNE_TEST_ASSERT(not(one==nan));
    ARIADNE_TEST_ASSERT((nan!=nan));
    ARIADNE_TEST_ASSERT((nan!=one));
    ARIADNE_TEST_ASSERT((one!=nan));
    ARIADNE_TEST_ASSERT(not(nan<=nan));
    ARIADNE_TEST_ASSERT(not(nan<=one));
    ARIADNE_TEST_ASSERT(not(one<=nan));
    ARIADNE_TEST_ASSERT(not(nan>=nan));
    ARIADNE_TEST_ASSERT(not(nan>=one));
    ARIADNE_TEST_ASSERT(not(one>=nan));
    ARIADNE_TEST_ASSERT(not(nan< nan));
    ARIADNE_TEST_ASSERT(not(nan< one));
    ARIADNE_TEST_ASSERT(not(one< nan));
    ARIADNE_TEST_ASSERT(not(nan> nan));
    ARIADNE_TEST_ASSERT(not(nan> one));
    ARIADNE_TEST_ASSERT(not(one> nan));

    ARIADNE_TEST_EQUALS(Float(pinf_d,precision).get_d(),pinf_d.get_d());
    ARIADNE_TEST_EQUALS(Float(ninf_d,precision).get_d(),ninf_d.get_d());
    ARIADNE_TEST_ASSERT(is_nan(Float(nan_d,precision)));

    ARIADNE_TEST_EQUALS(Float(Dyadic::inf(Sign::POSITIVE),precision).get_d(),pinf_d.get_d());
    ARIADNE_TEST_EQUALS(Float(Dyadic::inf(Sign::NEGATIVE),precision).get_d(),ninf_d.get_d());
    ARIADNE_TEST_ASSERT(is_nan(Float(Dyadic::inf(Sign::ZERO),precision)));

    ARIADNE_TEST_EQUALS(Float(Dyadic::inf(Sign::POSITIVE),upward,precision).get_d(),pinf_d.get_d());
    ARIADNE_TEST_EQUALS(Float(Dyadic::inf(Sign::NEGATIVE),upward,precision).get_d(),ninf_d.get_d());
    ARIADNE_TEST_ASSERT(is_nan(Float(Dyadic::inf(Sign::ZERO),upward,precision)));

    ARIADNE_TEST_EQUALS(Float(Rational::inf(Sign::POSITIVE),downward,precision).get_d(),pinf_d.get_d());
    ARIADNE_TEST_EQUALS(Float(Rational::inf(Sign::NEGATIVE),downward,precision).get_d(),ninf_d.get_d());
    ARIADNE_TEST_ASSERT(is_nan(Float(Rational::inf(Sign::ZERO),downward,precision)));
}


template<class PR> Void
TestFloatRounding<PR>::test_conversion_from_to()
{
    // Convert from dyadic
    Float f1(Dyadic(5,2u),precision);
    ARIADNE_TEST_EQUAL(f1,1.25_x);

    // Convert from self
    Float f2(f1,downward,precision);
    ARIADNE_TEST_EQUAL(f1,f2);
    Float f3(f1,upward,precision);
    ARIADNE_TEST_EQUAL(f1,f3);

    // Convert from integers
    Int n;
    n=std::numeric_limits<Int>::min();
    ARIADNE_TEST_ASSERT(Float(n,precision)==n);
    n=std::numeric_limits<Int>::max();
    ARIADNE_TEST_ASSERT(Float(n,precision)==n);
    Nat nn = std::numeric_limits<Nat>::max();
    ARIADNE_TEST_ASSERT(Float(nn,precision)==nn);

    // Convert to a rational
    ARIADNE_TEST_EQUAL(Rational(Float(1.0_x,precision)),Rational(1));
    ARIADNE_TEST_EQUAL(Rational(Float(2.25_x,precision)),Rational(9,4));

    // Convert from a rational
    Int num=1; Int den=3;
    Rational q(1,3);
    ARIADNE_TEST_COMPARE(Rational(Float(q,downward,precision)),<,q);
    ARIADNE_TEST_COMPARE(Rational(Float(Rational(num,den),upward,precision)),>,q);
    ARIADNE_TEST_COMPARE(Float(Rational(num,den),downward,precision),<,Float(Rational(num,den),upward,precision));

    // Convert from a negative rational
    num=-2; den=5;
    ARIADNE_TEST_COMPARE(Rational(Float(q,downward,precision)),<,q);
    ARIADNE_TEST_COMPARE(Rational(Float(q,upward,precision)),>,q);
    ARIADNE_TEST_COMPARE(Float(q,downward,precision),<,Float(q,upward,precision));

    // Convert into a dyadic
    Dyadic d1(f1);
    ARIADNE_TEST_EQUAL(d1,f1);
    ARIADNE_TEST_EQUALS(Dyadic(Float::inf(Sign::POSITIVE,precision)),Dyadic(std::numeric_limits<double>::infinity()));
    ARIADNE_TEST_EQUALS(Dyadic(Float::inf(Sign::NEGATIVE,precision)),Dyadic(-std::numeric_limits<double>::infinity()));
    ARIADNE_TEST_ASSERT(is_nan(Dyadic(Float::inf(Sign::ZERO,precision))));
}

template<> Void
TestFloatRounding<DoublePrecision>::test_conversion_between()
{
    // Keep cross-precision conversion coverage with the native backend.
    // A 100-bit directed 1/3 lies strictly between the corresponding
    // 53-bit rounded values and the exact rational 1/3.
    const MultiplePrecision source_precision(100_bits);
    const Rational exact(1,3);
    const FloatMP xap1(exact,downward,source_precision);
    const FloatDP xdp1(xap1,downward,DoublePrecision());
    ARIADNE_TEST_COMPARE(Rational(xdp1),<,Rational(xap1));
    const FloatMP xap2(exact,upward,source_precision);
    const FloatDP xdp2(xap2,upward,DoublePrecision());
    ARIADNE_TEST_COMPARE(Rational(xdp2),>,Rational(xap2));
    ARIADNE_TEST_COMPARE(xdp2,>,xdp1);
}

template<class PR> Void
TestFloatRounding<PR>::test_conversion_between()
{
    // The DoublePrecision-specific FloatMP-to-FloatDP conversion test
    // is specialized above; this case checks FloatDP-to-FloatMP.
    const FloatDP source(Rational(1,3),upward,DoublePrecision());
    const Float converted(source,PR(100_bits));
    ARIADNE_TEST_EQUALS(converted,source);

    const FloatDP exact_source(1.25_x,dp);
    const Float exact_converted(exact_source,PR(100_bits));
    ARIADNE_TEST_EQUALS(exact_converted,exact_source);
}

template<class PR> Void
TestFloatRounding<PR>::test_precision_regressions()
{
    // Exercise the native non-dyadic ratio path through the public FloatMP
    // API; the expected 1/3 mantissas follow exactly from integer division
    // of 2^(p+1) by three, independently of FloatMP.
    if constexpr (Same<PR,MultiplePrecision>) {
        FloatMPPhaseTiming timing_native_ratio("native_ratio",precision.bits());
        for (Bits bits : {2_bits,7_bits,31_bits,32_bits,33_bits,53_bits,
                          64_bits,65_bits,127_bits,128_bits,129_bits,
                          255_bits,256_bits,257_bits,1023_bits,1024_bits,
                          1100_bits,1536_bits}) {
            const MultiplePrecision p(bits);
            const Integer scale=pow(Integer(2),Nat(p.bits()+1u));
            const Integer floor_mantissa=quot(scale,Integer(3));
            const Integer rem_mantissa=scale-3*floor_mantissa;
            const bool nearest_increment=rem_mantissa*2>3;
            for (bool negative_input : {false,true}) {
                const Rational input=negative_input ? Rational(-1,3)
                                                     : Rational(1,3);
                const auto check=[&](FloatMP::RoundingModeType ap_mode,
                                     bool increment) {
                    const Integer magnitude=floor_mantissa+
                                             Integer(increment ? 1 : 0);
                    const Rational expected(negative_input ? -magnitude : magnitude,
                                             scale);
                    const FloatMP native(input,ap_mode,p);
                    ARIADNE_TEST_EQUALS(static_cast<Rational>(native),expected);
                    ARIADNE_TEST_EQUALS(static_cast<Rational>(native),
                                        floatmp_golden(p.bits(),ap_mode,input));
                };
                check(FloatMP::ROUND_TO_NEAREST,
                      nearest_increment);
                check(FloatMP::ROUND_DOWNWARD,
                      negative_input);
                check(FloatMP::ROUND_UPWARD,
                      !negative_input);
                check(FloatMP::ROUND_TOWARD_ZERO,
                      false);
            }

            // This non-dyadic pair forces multiprecision quotient/remainder.
            const Rational large_ratio(
                pow(Integer(2),Nat(1800u))+Integer(7),
                pow(Integer(2),Nat(1500u))+Integer(11));
            const auto verify=[&](auto ap_mode) {
                const FloatMP native(large_ratio,ap_mode,p);
                ARIADNE_TEST_EQUALS(static_cast<Rational>(native),
                                    floatmp_golden(p.bits(),ap_mode,large_ratio));
            };
            verify(FloatMP::ROUND_TO_NEAREST);
            verify(FloatMP::ROUND_DOWNWARD);
            verify(FloatMP::ROUND_UPWARD);
            verify(FloatMP::ROUND_TOWARD_ZERO);
        }
    }
    if constexpr (Same<PR,MultiplePrecision>) {
        // Exact halfway cases at p=2 (ties to even) and the signed
        // directed-rounding counterparts. All literals are independent
        // of both floating-point backends.
        const MultiplePrecision p(2_bits);
        for (Rational const& input : {Rational(5,4),Rational(7,4),
                                     Rational(-5,4),Rational(-7,4)}) {
            const Rational nearest_expected=
                input==Rational(5,4) ? Rational(1) :
                input==Rational(7,4) ? Rational(2) :
                input==Rational(-5,4) ? Rational(-1) : Rational(-2);
            const auto check=[&](FloatMP::RoundingModeType ap_mode,
                                 Rational const& exact) {
                const FloatMP candidate(input,ap_mode,p);
                ARIADNE_TEST_EQUALS(static_cast<Rational>(candidate),exact);
                ARIADNE_TEST_EQUALS(static_cast<Rational>(candidate),
                                    floatmp_golden(p.bits(),ap_mode,input));
            };
            check(FloatMP::ROUND_TO_NEAREST,
                  nearest_expected);
            const Rational lower=input==Rational(5,4) ? Rational(1) :
                input==Rational(7,4) ? Rational(3,2) :
                input==Rational(-5,4) ? Rational(-3,2) : Rational(-2);
            const Rational upper=input==Rational(5,4) ? Rational(3,2) :
                input==Rational(7,4) ? Rational(2) :
                input==Rational(-5,4) ? Rational(-1) : Rational(-3,2);
            check(FloatMP::ROUND_DOWNWARD,lower);
            check(FloatMP::ROUND_UPWARD,upper);
            check(FloatMP::ROUND_TOWARD_ZERO,
                  input<0 ? upper : lower);
        }
    }
    if constexpr (Same<PR,MultiplePrecision>) {
        // Broad, exact non-dyadic quotients across the limb boundaries.
        // Stored exact integer-rounding references replace backend cross-checks;
        // independently expected 1/3 and halfway cases remain unchanged.
        for (Bits bits : {53_bits,128_bits,256_bits,1024_bits,1536_bits}) {
            const MultiplePrecision p(bits);
            const Integer power=pow(Integer(2),Nat(p.bits()));
            for (Int offset : {Int(1),Int(7),Int(19),Int(257)}) {
                const Rational positive_ratio(
                    power+Integer(offset),
                    2*power-Integer(offset+2));
                for (bool negative_input : {false,true}) {
                    const Rational exact_input=negative_input
                        ? -positive_ratio : positive_ratio;
                    const auto verify=[&](auto ap_mode) {
                            const FloatMP candidate(exact_input,ap_mode,p);
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(candidate),
                                floatmp_golden(p.bits(),ap_mode,exact_input));
                        };
                    verify(FloatMP::ROUND_TO_NEAREST);
                    verify(FloatMP::ROUND_DOWNWARD);
                    verify(FloatMP::ROUND_UPWARD);
                    verify(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }
    }
    if constexpr (Same<PR,MultiplePrecision>) {
        // Sparse multi-limb denominators remain sparse after Knuth
        // normalization; verify exact rounding for each sign and mode.
        for (Bits bits : {128_bits,256_bits,1024_bits,1536_bits}) {
            const MultiplePrecision p(bits);
            const Integer power=pow(Integer(2),Nat(p.bits()));
            for (Int offset : {Int(1),Int(17),Int(257)}) {
                const Rational positive_ratio(
                    power+Integer(offset),2*power+Integer(offset+2));
                for (bool negative_input : {false,true}) {
                    const Rational input=negative_input
                        ? -positive_ratio : positive_ratio;
                    const auto verify=[&](auto ap_mode) {
                            const FloatMP candidate(input,ap_mode,p);
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(candidate),
                                floatmp_golden(p.bits(),ap_mode,input));
                        };
                    verify(FloatMP::ROUND_TO_NEAREST);
                    verify(FloatMP::ROUND_DOWNWARD);
                    verify(FloatMP::ROUND_UPWARD);
                    verify(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }
    }
    if constexpr (Same<PR,MultiplePrecision>) {
        // Fused multiply/subtract carry can include a borrow from the
        // previous word. Stress runs of zero and all-one low limbs,
        // quotient-correction boundaries, signs and all rounding modes.
        // Stored exact integer-rounding references replace backend cross-checks;
        // exact 1/3 and halfway tests above remain in place.
        for (Bits bits : {65_bits,128_bits,129_bits,255_bits,256_bits,
                          257_bits,511_bits,512_bits,1023_bits,1024_bits,
                          1025_bits,1536_bits}) {
            const MultiplePrecision p(bits);
            const Integer power=pow(Integer(2),Nat(p.bits()));
            for (Int offset : {Int(1),Int(3),Int(257)}) {
                const Integer x=power-Integer(offset);
                const Integer y=power+Integer(offset+2);
                for (bool reverse : {false,true}) {
                    const Rational positive_ratio=reverse
                        ? Rational(y,x) : Rational(x,y);
                    for (bool negative_input : {false,true}) {
                        const Rational input=negative_input
                            ? -positive_ratio : positive_ratio;
                        const auto verify=[&](auto ap_mode) {
                            const FloatMP candidate(input,ap_mode,p);
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(candidate),
                                floatmp_golden(p.bits(),ap_mode,input));
                        };
                        verify(FloatMP::ROUND_TO_NEAREST);
                        verify(FloatMP::ROUND_DOWNWARD);
                        verify(FloatMP::ROUND_UPWARD);
                        verify(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }
    }
    // Preserve exact-rounding regression coverage for FloatMP with the
    // public arithmetic API over all required precisions and modes.
    if constexpr (Same<PR,MultiplePrecision>) {
        using F=RawFloat<PR>;
        for (Bits bits : {2_bits,7_bits,24_bits,53_bits,113_bits,256_bits,1024_bits}) {
            const PR pr(bits);
            for (auto mode : {F::ROUND_TO_NEAREST,F::ROUND_DOWNWARD,
                              F::ROUND_UPWARD,F::ROUND_TOWARD_ZERO}) {
                const F x(Rational(17,13),mode,pr);
                const F y(Rational(-11,19),mode,pr);
                const Rational a=static_cast<Rational>(x);
                const Rational b=static_cast<Rational>(y);
                const auto check=[&](F const& actual,Rational const& exact) {
                    ARIADNE_TEST_EQUALS(actual,F(exact,mode,pr));
                };
                check(add(mode,x,y),a+b);
                check(sub(mode,x,y),a-b);
                check(mul(mode,x,y),a*b);
                check(div(mode,x,y),a/b);
                check(fma(mode,x,y,x),a*b+a);
                check(pow(mode,x,Nat(7u)),a*a*a*a*a*a*a);

                ARIADNE_TEST_EQUALS(sqrt(mode,F(4,pr)),F(2,pr));
                if (bits>=4_bits) {
                    ARIADNE_TEST_EQUALS(sqrt(mode,F(9,pr)),F(3,pr));
                } else {
                    // At two bits, 9 rounds to 8 except in upward mode,
                    // where it rounds to 12, before sqrt is evaluated.
                    const F rounded_nine(Rational(9),mode,pr);
                    const bool round_upward=mode==F::ROUND_UPWARD;
                    ARIADNE_TEST_EQUALS(static_cast<Rational>(rounded_nine),
                                        Rational(round_upward?12:8));
                    const Int expected=round_upward?4:
                        (mode==F::ROUND_TO_NEAREST?3:2);
                    ARIADNE_TEST_EQUALS(sqrt(mode,rounded_nine),F(expected,pr));
                }

                const F two_value(2,pr);
                const Rational root_lower=static_cast<Rational>(
                    sqrt(F::ROUND_DOWNWARD,two_value));
                const Rational root_upper=static_cast<Rational>(
                    sqrt(F::ROUND_UPWARD,two_value));
                ARIADNE_TEST_ASSERT(root_lower*root_lower<=Rational(2));
                ARIADNE_TEST_ASSERT(root_upper*root_upper>=Rational(2));
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Fixed exact dyadic literals, independently generated from
            // rigorous rational enclosures, independently of FloatMP
            // conversion or FloatMP arithmetic.
            // For each row the nearest-even mantissa M has exponent e;
            // down/up use M plus their separately certified integer delta.
            struct CertifiedLiteral {
                char const* operation;
                unsigned bits;
                Int exponent;
                char const* nearest_mantissa;
                Int down_offset;
                Int up_offset;
            };
            static const CertifiedLiteral expected[] = {
                {"sqrt2",128u,-127,
                 "240615969168004511545033772477625056927",0,1},
                {"sqrt2",256u,-255,
                 "81877371507464127617551201542979628307507432471243237061821853600756754782485",0,1},
                {"sqrt2",1024u,-1023,
                 "127116100615364628366052028422228385109549339363277773668712774121205113882836385245579298214279956427891394679477132134293915806934836260631367995992135823438869876926245592272026536661973063905698010856540683205419763801114347235341470756348731874877833367845541044116835881896947807071342075615375081392715",-1,0},
                {"exp_half",128u,-127,
                 "280515388193368458015406427511040113881",-1,0},
                {"exp_half",256u,-255,
                 "95454440252185294983886966005124823429889078541787437355320193447306720442726",0,1},
                {"exp_half",1024u,-1023,
                 "148194745481954716082317888363455823456070480866279839338113084411204029676709818203132581631311711177728819167258549487238671063611984378104029593810429844208825070153844097581609333789019814401213030206345714196895102345447040710261659875317126576627766610707108418090765348768199197232529394147236509021090",-1,0},
                {"log2",128u,-128,
                 "235865763225513294137944142764154484399",0,1},
                {"log2",256u,-256,
                 "80260960185991308862233904206310070533990667611589946606122867505419956976172",-1,0},
                {"log2",1024u,-1024,
                 "124606592794178379834942157015327973260725777549818548919889062160856542278503538017240508556319343545235776417566655832902220557522478644110217075261513836839396089380897954160502582582011535456352891435689181660887803024528449058262732172678320723819450786384318433514902743298912695256684345870847729194438",0,1},
                {"pi",128u,-126,
                 "267257146016241686964920093290467695825",0,1},
                {"pi",256u,-254,
                 "90942894222941581070058735694432465663348344332098107489693037779484723616546",0,1},
                {"pi",1024u,-1022,
                 "141190488647306426951497947939679424974931379369199690146440126562466409291619758234410430492087736172907390888290563508264155720409054555597508847013932624303633293809435642998784554044791800749147419391451534200559935013628202594044220234864921212360009141120460241474543779193546148756924903705307226952602",0,1},
                {"sin_half",53u,-54,"8636562708039152",-1,0},
                {"sin_half",128u,-129,"326280114077167906689097124699132648243",0,1},
                {"sin_half",256u,-257,
                 "111027369497412508784314791148902648293994217576031053546226341313477387335312",-1,0},
                {"sin_half",1024u,-1025,
                 "172371999885288789128716840567144047953759349476369256334737458968483142997910650916441548737119101143602980846427532337793230109910643907421443819149435799871037387373869579663640966923121927526805195778868067048300347801194629542755640128151062422101982536195161049560508328652059747480969032691608841231624",-1,0},
                {"cos_half",53u,-53,"7904560997432656",-1,0},
                {"cos_half",128u,-128,"298625871328596996571466326832836172165",0,1},
                {"cos_half",256u,-256,
                 "101617118319522600545601981648807607350213579319835970884288805016705398675944",0,1},
                {"cos_half",1024u,-1024,
                 "157762414678520649592707644765600265788274071390691627946472084316739856684105644928097150639079359415661197561565092641738842940635149605067256924545601776964818603475456673753501326513258696538222222055589487086313973209608156611433960402938690734994439009100183662423673858021971513755127657214343043463501",0,1}
            };
            for (auto const& row : expected) {
                const MultiplePrecision p{static_cast<std::size_t>(row.bits)};
                const Integer mantissa(String(row.nearest_mantissa));
                const Integer denominator=pow(Integer(2),Nat(-row.exponent));
                const auto exact=[&](Int delta) {
                    return Rational(mantissa+Integer(delta),denominator);
                };
                const FloatMP operand_half(Rational(1,2),p);
                const FloatMP operand_two(2,p);
                const auto evaluate=[&](auto mode) -> FloatMP {
                    if (std::string(row.operation)=="sqrt2") return sqrt(mode,operand_two);
                    if (std::string(row.operation)=="exp_half") return exp(mode,operand_half);
                    if (std::string(row.operation)=="log2") return log(mode,operand_two);
                    if (std::string(row.operation)=="sin_half") return sin(mode,operand_half);
                    if (std::string(row.operation)=="cos_half") return cos(mode,operand_half);
                    return FloatMP::pi(mode,p);
                };
                const FloatMP nearest=evaluate(FloatMP::ROUND_TO_NEAREST);
                const FloatMP lower=evaluate(FloatMP::ROUND_DOWNWARD);
                const FloatMP upper=evaluate(FloatMP::ROUND_UPWARD);
                const FloatMP zero=evaluate(FloatMP::ROUND_TOWARD_ZERO);
                ARIADNE_TEST_EQUALS(static_cast<Rational>(nearest),exact(0));
                ARIADNE_TEST_EQUALS(static_cast<Rational>(lower),exact(row.down_offset));
                ARIADNE_TEST_EQUALS(static_cast<Rational>(upper),exact(row.up_offset));
                ARIADNE_TEST_EQUALS(static_cast<Rational>(zero),exact(row.down_offset));
                ARIADNE_TEST_EQUALS(nearest.precision().bits(),row.bits);
                ARIADNE_TEST_EQUALS(lower.precision().bits(),row.bits);
                ARIADNE_TEST_EQUALS(upper.precision().bits(),row.bits);
                ARIADNE_TEST_EQUALS(zero.precision().bits(),row.bits);
            }
        }


        if constexpr (Same<PR,MultiplePrecision>) {
            // Independent exact dyadic oracle values around +/-1/2:
            // both rational Taylor enclosure endpoints were rounded
            // offline to the SAME target literal for every mode.
            // This is an independent golden oracle, not a FloatMP-derived
            // expectation. Check both signs and both sides of 1/2.
            struct NearHalfTrigLiteral {
                char const* function;
                unsigned bits;
                unsigned offset_bits;
                Int offset_sign;
                bool negative_input;
                Int exponent;
                char const* nearest_mantissa;
                Int down_offset;
                Int up_offset;
                Int zero_offset;
            };
            static const NearHalfTrigLiteral expected[] = {
                {"sin",53u,50u,-1,false,-54,"8636562708039138",-1,0,-1},
                {"cos",53u,50u,-1,false,-53,"7904560997432659",0,1,0},
                {"sin",53u,50u,-1,true,-54,"-8636562708039138",0,1,1},
                {"cos",53u,50u,-1,true,-53,"7904560997432659",0,1,0},
                {"sin",53u,50u,1,false,-54,"8636562708039166",-1,0,-1},
                {"cos",53u,50u,1,false,-53,"7904560997432652",-1,0,-1},
                {"sin",53u,50u,1,true,-54,"-8636562708039166",0,1,1},
                {"cos",53u,50u,1,true,-53,"7904560997432652",-1,0,-1},
                {"sin",128u,125u,-1,false,-129,"326280114077167906689097124699132648229",0,1,0},
                {"cos",128u,125u,-1,false,-128,"298625871328596996571466326832836172169",0,1,0},
                {"sin",128u,125u,-1,true,-129,"-326280114077167906689097124699132648229",-1,0,0},
                {"cos",128u,125u,-1,true,-128,"298625871328596996571466326832836172169",0,1,0},
                {"sin",128u,125u,1,false,-129,"326280114077167906689097124699132648257",0,1,0},
                {"cos",128u,125u,1,false,-128,"298625871328596996571466326832836172162",-1,0,-1},
                {"sin",128u,125u,1,true,-129,"-326280114077167906689097124699132648257",-1,0,0},
                {"cos",128u,125u,1,true,-128,"298625871328596996571466326832836172162",-1,0,-1},
                {"sin",256u,253u,-1,false,-257,"111027369497412508784314791148902648293994217576031053546226341313477387335298",-1,0,-1},
                {"cos",256u,253u,-1,false,-256,"101617118319522600545601981648807607350213579319835970884288805016705398675948",-1,0,-1},
                {"sin",256u,253u,-1,true,-257,"-111027369497412508784314791148902648293994217576031053546226341313477387335298",0,1,1},
                {"cos",256u,253u,-1,true,-256,"101617118319522600545601981648807607350213579319835970884288805016705398675948",-1,0,-1},
                {"sin",256u,253u,1,false,-257,"111027369497412508784314791148902648293994217576031053546226341313477387335326",-1,0,-1},
                {"cos",256u,253u,1,false,-256,"101617118319522600545601981648807607350213579319835970884288805016705398675940",0,1,0},
                {"sin",256u,253u,1,true,-257,"-111027369497412508784314791148902648293994217576031053546226341313477387335326",0,1,1},
                {"cos",256u,253u,1,true,-256,"101617118319522600545601981648807607350213579319835970884288805016705398675940",0,1,0},
                {"sin",1024u,1021u,-1,false,-1025,"172371999885288789128716840567144047953759349476369256334737458968483142997910650916441548737119101143602980846427532337793230109910643907421443819149435799871037387373869579663640966923121927526805195778868067048300347801194629542755640128151062422101982536195161049560508328652059747480969032691608841231610",-1,0,-1},
                {"cos",1024u,1021u,-1,false,-1024,"157762414678520649592707644765600265788274071390691627946472084316739856684105644928097150639079359415661197561565092641738842940635149605067256924545601776964818603475456673753501326513258696538222222055589487086313973209608156611433960402938690734994439009100183662423673858021971513755127657214343043463505",0,1,0},
                {"sin",1024u,1021u,-1,true,-1025,"-172371999885288789128716840567144047953759349476369256334737458968483142997910650916441548737119101143602980846427532337793230109910643907421443819149435799871037387373869579663640966923121927526805195778868067048300347801194629542755640128151062422101982536195161049560508328652059747480969032691608841231610",0,1,1},
                {"cos",1024u,1021u,-1,true,-1024,"157762414678520649592707644765600265788274071390691627946472084316739856684105644928097150639079359415661197561565092641738842940635149605067256924545601776964818603475456673753501326513258696538222222055589487086313973209608156611433960402938690734994439009100183662423673858021971513755127657214343043463505",0,1,0},
                {"sin",1024u,1021u,1,false,-1025,"172371999885288789128716840567144047953759349476369256334737458968483142997910650916441548737119101143602980846427532337793230109910643907421443819149435799871037387373869579663640966923121927526805195778868067048300347801194629542755640128151062422101982536195161049560508328652059747480969032691608841231638",-1,0,-1},
                {"cos",1024u,1021u,1,false,-1024,"157762414678520649592707644765600265788274071390691627946472084316739856684105644928097150639079359415661197561565092641738842940635149605067256924545601776964818603475456673753501326513258696538222222055589487086313973209608156611433960402938690734994439009100183662423673858021971513755127657214343043463497",0,1,0},
                {"sin",1024u,1021u,1,true,-1025,"-172371999885288789128716840567144047953759349476369256334737458968483142997910650916441548737119101143602980846427532337793230109910643907421443819149435799871037387373869579663640966923121927526805195778868067048300347801194629542755640128151062422101982536195161049560508328652059747480969032691608841231638",0,1,1},
                {"cos",1024u,1021u,1,true,-1024,"157762414678520649592707644765600265788274071390691627946472084316739856684105644928097150639079359415661197561565092641738842940635149605067256924545601776964818603475456673753501326513258696538222222055589487086313973209608156611433960402938690734994439009100183662423673858021971513755127657214343043463497",0,1,0}
            };
            for (auto const& row : expected) {
                const MultiplePrecision p{static_cast<std::size_t>(row.bits)};
                const Rational delta(Integer(1),
                    pow(Integer(2),Nat(row.offset_bits)));
                const Rational positive = Rational(1,2) +
                    (row.offset_sign < 0 ? -delta : delta);
                const Rational exact_input =
                    row.negative_input ? -positive : positive;
                const FloatMP argument(exact_input,FloatMP::ROUND_TO_NEAREST,p);
                // All selected dyadic arguments are p-bit representable:
                // the checked function must not see a different input.
                ARIADNE_TEST_EQUALS(static_cast<Rational>(argument),exact_input);
                const Integer mantissa(String(row.nearest_mantissa));
                const Integer denominator =
                    pow(Integer(2),Nat(-row.exponent));
                const auto exact = [&](Int offset) {
                    return Rational(mantissa+Integer(offset),denominator);
                };
                const auto evaluate = [&](auto mode) {
                    return std::string(row.function) == "sin"
                        ? sin(mode,argument) : cos(mode,argument);
                };
                const FloatMP nearest=evaluate(FloatMP::ROUND_TO_NEAREST);
                const FloatMP lower=evaluate(FloatMP::ROUND_DOWNWARD);
                const FloatMP upper=evaluate(FloatMP::ROUND_UPWARD);
                const FloatMP zero=evaluate(FloatMP::ROUND_TOWARD_ZERO);
                ARIADNE_TEST_EQUALS(static_cast<Rational>(nearest),exact(0));
                ARIADNE_TEST_EQUALS(static_cast<Rational>(lower),exact(row.down_offset));
                ARIADNE_TEST_EQUALS(static_cast<Rational>(upper),exact(row.up_offset));
                ARIADNE_TEST_EQUALS(static_cast<Rational>(zero),exact(row.zero_offset));
                ARIADNE_TEST_EQUALS(nearest.precision().bits(),row.bits);
                ARIADNE_TEST_EQUALS(lower.precision().bits(),row.bits);
                ARIADNE_TEST_EQUALS(upper.precision().bits(),row.bits);
                ARIADNE_TEST_EQUALS(zero.precision().bits(),row.bits);
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Independent mathematically certified exact dyadic literals
            // of atan(+-1/2 +- 2^-(p-3)). No FloatMP is used
            // to derive the expected mantissas or directed offsets.
            struct AtanExpected {
                unsigned bits, offset_bits;
                Int offset_sign;
                bool negative_input;
                Int exponent;
                char const* nearest_mantissa;
                Int down_offset, up_offset, zero_offset;
            };
            static const AtanExpected expected[] = {
                {53u,50u,-1,false,-54,"8352332796508995",-1,0,-1},
                {53u,50u,-1,true,-54,"-8352332796508995",0,1,1},
                {53u,50u,1,false,-54,"8352332796509020",0,1,0},
                {53u,50u,1,true,-54,"-8352332796509020",-1,0,0},
                {128u,125u,-1,false,-129,"315542211616056235508655905213297364215",-1,0,-1},
                {128u,125u,-1,true,-129,"-315542211616056235508655905213297364215",0,1,1},
                {128u,125u,1,false,-129,"315542211616056235508655905213297364240",0,1,0},
                {128u,125u,1,true,-129,"-315542211616056235508655905213297364240",-1,0,0},
                {256u,253u,-1,false,-257,"107373450632179258931464468992502143670939545571605869416383036103061999339261",0,1,0},
                {256u,253u,-1,true,-257,"-107373450632179258931464468992502143670939545571605869416383036103061999339261",-1,0,0},
                {256u,253u,1,false,-257,"107373450632179258931464468992502143670939545571605869416383036103061999339287",-1,0,-1},
                {256u,253u,1,true,-257,"-107373450632179258931464468992502143670939545571605869416383036103061999339287",0,1,1},
                {1024u,1021u,-1,false,-1025,"166699224739215292881125319598038346462953795978162326636292450871431032119258242184401154484764005946607347345377978929188361233397565402058863543332761952867976377120522096751547155783154721876987164034590575540231000736940058769795865598096281237086736548576233083171219137660006199486594544307781066487450",0,1,0},
                {1024u,1021u,-1,true,-1025,"-166699224739215292881125319598038346462953795978162326636292450871431032119258242184401154484764005946607347345377978929188361233397565402058863543332761952867976377120522096751547155783154721876987164034590575540231000736940058769795865598096281237086736548576233083171219137660006199486594544307781066487450",-1,0,0},
                {1024u,1021u,1,false,-1025,"166699224739215292881125319598038346462953795978162326636292450871431032119258242184401154484764005946607347345377978929188361233397565402058863543332761952867976377120522096751547155783154721876987164034590575540231000736940058769795865598096281237086736548576233083171219137660006199486594544307781066487476",-1,0,-1},
                {1024u,1021u,1,true,-1025,"-166699224739215292881125319598038346462953795978162326636292450871431032119258242184401154484764005946607347345377978929188361233397565402058863543332761952867976377120522096751547155783154721876987164034590575540231000736940058769795865598096281237086736548576233083171219137660006199486594544307781066487476",0,1,1}
            };
            for (auto const& row : expected) {
                const MultiplePrecision p{static_cast<std::size_t>(row.bits)};
                const Rational delta(Integer(1),
                                     pow(Integer(2),Nat(row.offset_bits)));
                const Rational y=Rational(1,2) +
                    (row.offset_sign<0 ? -delta : delta);
                const Rational input=row.negative_input ? -y : y;
                const FloatMP argument(input,FloatMP::ROUND_TO_NEAREST,p);
                ARIADNE_TEST_EQUALS(static_cast<Rational>(argument),input);
                const Integer mantissa(String(row.nearest_mantissa));
                const Integer denominator=pow(Integer(2),Nat(-row.exponent));
                const auto target=[&](Int offset) {
                    return Rational(mantissa+Integer(offset),denominator);
                };
                const FloatMP nearest_result=atan(FloatMP::ROUND_TO_NEAREST,argument);
                const FloatMP downward_result=atan(FloatMP::ROUND_DOWNWARD,argument);
                const FloatMP upward_result=atan(FloatMP::ROUND_UPWARD,argument);
                const FloatMP toward_zero_result=atan(FloatMP::ROUND_TOWARD_ZERO,argument);
                ARIADNE_TEST_EQUALS(static_cast<Rational>(nearest_result),target(0));
                ARIADNE_TEST_EQUALS(static_cast<Rational>(downward_result),target(row.down_offset));
                ARIADNE_TEST_EQUALS(static_cast<Rational>(upward_result),target(row.up_offset));
                ARIADNE_TEST_EQUALS(static_cast<Rational>(toward_zero_result),target(row.zero_offset));
                ARIADNE_TEST_EQUALS(nearest_result.precision().bits(),row.bits);
                ARIADNE_TEST_EQUALS(downward_result.precision().bits(),row.bits);
                ARIADNE_TEST_EQUALS(upward_result.precision().bits(),row.bits);
                ARIADNE_TEST_EQUALS(toward_zero_result.precision().bits(),row.bits);
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Integer square-root remainder and halfway rounding:
            // freeze the correctly rounded results computed independently
            // from exact integer square root and midpoint comparisons.
            // Keep the same rounded inputs and all four rounding modes.
            for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                const Rational step(Integer(1),
                                    pow(Integer(2),Nat(mp.bits())));
                const Rational exact_values[]={
                    Rational(1),Rational(2),Rational(3),Rational(4),
                    Rational(5),Rational(9),Rational(15,16),
                    Rational(17,16),Rational(1,4),Rational(9,16),
                    Rational(1)-step,Rational(1)+step,
                    Rational(4)-step,Rational(4)+step
                };
                const auto check=[&](FloatMP const& input) {
                    const Rational exact_input=static_cast<Rational>(input);
                    const auto check_mode=[&](auto am) {
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(sqrt(am,input)),
                            floatmp_sqrt_golden(mp.bits(),am,exact_input));
                    };
                    check_mode(FloatMP::ROUND_TO_NEAREST);
                    check_mode(FloatMP::ROUND_DOWNWARD);
                    check_mode(FloatMP::ROUND_UPWARD);
                    check_mode(FloatMP::ROUND_TOWARD_ZERO);
                };
                for (Rational const& value : exact_values) {
                    check(FloatMP(value,FloatMP::ROUND_TO_NEAREST,mp));
                }
                const FloatMP scale_source(Rational(5,4),mp);
                for (Int shift : {Int(-2000),Int(-127),Int(127),Int(2000)}) {
                    check(shft(scale_source,shift));
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Equal-precision normalized values with distinct exponents:
            // compare against exact Rational ordering, including negative
            // values, reversal, same-exponent cases and mixed precision.
            // The min/max comparisons exercise the int64 exponent limits
            // without materializing impossibly large Rational values.
            for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                const FloatMP low(Rational(3,4),mp);
                const FloatMP high(Rational(9,8),mp);
                const FloatMP same_exp_a(Rational(5,4),mp);
                const FloatMP same_exp_b(Rational(7,4),mp);
                const FloatMP operands[]={
                    low,high,neg(low),neg(high),
                    same_exp_a,same_exp_b,neg(same_exp_a),neg(same_exp_b)
                };
                const auto check=[&](FloatMP const& lhs,FloatMP const& rhs) {
                    ARIADNE_TEST_ASSERT(
                        cmp(lhs,rhs)==cmp(static_cast<Rational>(lhs),
                                          static_cast<Rational>(rhs)));
                };
                for (FloatMP const& lhs : operands) {
                    for (FloatMP const& rhs : operands) {
                        check(lhs,rhs);
                    }
                }

                const FloatMP other_precision(Rational(5,4),
                                              MultiplePrecision(24_bits));
                check(same_exp_a,other_precision);
                check(other_precision,high);
                check(high,other_precision);

                const FloatMP smallest=FloatMP::min(mp);
                const FloatMP largest=FloatMP::max(mp);
                ARIADNE_TEST_ASSERT(largest>smallest);
                ARIADNE_TEST_ASSERT(neg(largest)<neg(smallest));
                ARIADNE_TEST_ASSERT(largest>neg(largest));
                ARIADNE_TEST_ASSERT(smallest>neg(smallest));
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Portable small-limb quotient/remainder, with exact golden
            // values on identical rounded binary inputs.
            // 7/53/64 exercise small and capacity-boundary cases;
            // 128/256 exercise wider native quotient/remainder.
            for (Bits bits : {7_bits,53_bits,64_bits,128_bits,256_bits}) {
                const MultiplePrecision mp(bits);
                const Integer base = pow(Integer(2),Nat(mp.bits()));
                for (unsigned i = 0; i < 32u; ++i) {
                    const Integer numerator=base+Integer(13u+2u*i);
                    const Integer denominator=2*base+Integer(17u+3u*i);
                    const Rational exact(numerator,denominator);
                    const Rational divisor(
                        base+Integer(31u+5u*i),
                        3*base+Integer(19u+7u*i));
                    for (Int sign : {Int(-1),Int(1)}) {
                        const Rational signed_exact=sign<0 ? -exact : exact;
                        const FloatMP input(signed_exact,
                            FloatMP::ROUND_TO_NEAREST,mp);
                        const FloatMP other(divisor,
                            FloatMP::ROUND_TO_NEAREST,mp);
                        // Divide the values actually rounded on input,
                        // not the original exact rationals.
                        const Rational exact_quotient=
                            static_cast<Rational>(input)/
                            static_cast<Rational>(other);
                        const auto check=[&](auto am) {
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(FloatMP(signed_exact,am,mp)),
                                floatmp_golden(mp.bits(),am,signed_exact));
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(div(am,input,other)),
                                floatmp_golden(mp.bits(),am,exact_quotient));
                        };
                        check(FloatMP::ROUND_TO_NEAREST);
                        check(FloatMP::ROUND_DOWNWARD);
                        check(FloatMP::ROUND_UPWARD);
                        check(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // The general non-unit-denominator path of from_ratio()
            // must round both signs correctly for division and Rational
            // construction, including exact and inexact quotients.
            for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                for (Int numerator : {Int(-17),Int(-3),Int(1),Int(17)}) {
                    for (Int denominator : {Int(-16),Int(-7),Int(-3),
                                            Int(3),Int(7),Int(16)}) {
                        const FloatMP a(numerator,mp);
                        const FloatMP b(denominator,mp);
                        const Rational ratio(numerator,denominator);
                        const auto check=[&](auto ap_mode) {
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(div(ap_mode,a,b)),
                                floatmp_golden(mp.bits(),ap_mode,ratio));
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(FloatMP(ratio,ap_mode,mp)),
                                floatmp_golden(mp.bits(),ap_mode,ratio));
                        };
                        check(FloatMP::ROUND_TO_NEAREST);
                        check(FloatMP::ROUND_DOWNWARD);
                        check(FloatMP::ROUND_UPWARD);
                        check(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // A dyadic Rational uses the direct binary-significand path.
            // Exercise both signs and exact/tie cases against frozen exact values.
            for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                const Integer power=pow(Integer(2),Nat(mp.bits()));
                const Integer numerators[]={
                    power-Integer(1), power+Integer(1),
                    Integer(1)-power, -power-Integer(1)
                };
                for (unsigned shift : {1u,7u,31u,130u}) {
                    const Integer denominator=pow(Integer(2),Nat(shift));
                    for (Integer const& numerator : numerators) {
                        const Rational exact(numerator,denominator);
                        const auto check=[&](auto am) {
                            const FloatMP result(exact,am,mp);
                            ARIADNE_TEST_EQUALS(static_cast<Rational>(result),
                                                floatmp_golden(mp.bits(),am,exact));
                        };
                        check(FloatMP::ROUND_TO_NEAREST);
                        check(FloatMP::ROUND_DOWNWARD);
                        check(FloatMP::ROUND_UPWARD);
                        check(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Differential checks for direct binary-integer rounding.
            // Construct values just above and below precision boundaries,
            // including both signs and significand-carry cases.
            for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                const Integer power_of_two=pow(Integer(2),Nat(mp.bits()));
                const Integer values[]={
                    power_of_two+Integer(1),power_of_two-Integer(1),
                    -power_of_two-Integer(1),-power_of_two+Integer(1),
                    Integer(3),Integer(-3)
                };
                for (Integer const& exact : values) {
                    const auto check=[&](auto ap_mode) {
                        const FloatMP a(exact,ap_mode,mp);
                        const Rational rounded_a=static_cast<Rational>(a);
                        ARIADNE_TEST_EQUALS(rounded_a,
                            floatmp_golden(mp.bits(),ap_mode,Rational(exact)));

                        const FloatMP b(-3,mp);
                        const FloatMP c(1,mp);
                        const Rational exact_b(-3);
                        const Rational exact_c(1);
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(add(ap_mode,a,b)),
                            floatmp_golden(mp.bits(),ap_mode,rounded_a+exact_b));
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(sub(ap_mode,a,b)),
                            floatmp_golden(mp.bits(),ap_mode,rounded_a-exact_b));
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(mul(ap_mode,a,b)),
                            floatmp_golden(mp.bits(),ap_mode,rounded_a*exact_b));
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(fma(ap_mode,a,b,c)),
                            floatmp_golden(mp.bits(),ap_mode,
                                           rounded_a*exact_b+exact_c));
                    };
                    check(FloatMP::ROUND_TO_NEAREST);
                    check(FloatMP::ROUND_DOWNWARD);
                    check(FloatMP::ROUND_UPWARD);
                    check(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Cancellation, ties and large exponent gaps: use independently
            // frozen exact rational values of the *rounded binary operands*.
            for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                const FloatMP one(1,mp);
                const FloatMP fraction(Rational(7,11),FloatMP::ROUND_TO_NEAREST,mp);
                const FloatMP tiny=shft(one,-static_cast<Int>(mp.bits()));
                const FloatMP values[]={
                    one,neg(one),fraction,neg(fraction),
                    tiny,neg(tiny),shft(one,Int(-200)),
                    shft(one,Int(200)),shft(one,Int(-1400)),
                    shft(one,Int(1400))
                };
                for (FloatMP const& x : values) {
                    for (FloatMP const& y : values) {
                        const Rational exact_x=static_cast<Rational>(x);
                        const Rational exact_y=static_cast<Rational>(y);
                        const Rational exact_sum=exact_x+exact_y;
                        const Rational exact_difference=exact_x-exact_y;
                        const auto verify=[&](auto am) {
                            ARIADNE_TEST_EQUALS(static_cast<Rational>(add(am,x,y)),
                                floatmp_add_sub_golden(mp.bits(),am,exact_sum));
                            ARIADNE_TEST_EQUALS(static_cast<Rational>(sub(am,x,y)),
                                floatmp_add_sub_golden(mp.bits(),am,exact_difference));
                        };
                        verify(FloatMP::ROUND_TO_NEAREST);
                        verify(FloatMP::ROUND_DOWNWARD);
                        verify(FloatMP::ROUND_UPWARD);
                        verify(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
                const FloatMP z(0,mp);
                ARIADNE_TEST_ASSERT(sub(FloatMP::ROUND_DOWNWARD,one,one).negative());
                ARIADNE_TEST_ASSERT(not sub(FloatMP::ROUND_UPWARD,one,one).negative());
                ARIADNE_TEST_ASSERT(add(FloatMP::ROUND_DOWNWARD,z,neg(z)).negative());
                ARIADNE_TEST_ASSERT(is_nan(sub(FloatMP::ROUND_TO_NEAREST,
                    FloatMP::inf(Sign::POSITIVE,mp),
                    FloatMP::inf(Sign::POSITIVE,mp))));
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Exercise precision-dependent exp reduction and reconstruction
            // against fixed exact binary outputs from rigorously certified
            // rational Taylor enclosures (independently of FloatMP).
            for (Bits bits : {53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                for (Dyadic const& argument :
                     {Dyadic(1,2u),Dyadic(-1,2u),Dyadic(3,2u),
                      Dyadic(-3,2u),Dyadic(1,1u),Dyadic(-1,1u),
                      Dyadic(1,0u),Dyadic(-1,0u),Dyadic(3,0u),
                      Dyadic(-3,0u),Dyadic(7,0u),Dyadic(-7,0u),
                      Dyadic(1,16u),Dyadic(-1,16u),
                      Dyadic(1,64u),Dyadic(-1,64u),
                      Dyadic(1023,10u),Dyadic(-1023,10u),
                      Dyadic(256,0u),Dyadic(-256,0u)}) {
                    const FloatMP a(argument,mp);
                    const Rational exact_input(argument);
                    const auto compare_exponential=[&](auto ap_mode) {
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(exp(ap_mode,a)),
                            floatmp_exp_golden(mp.bits(),ap_mode,exact_input));
                    };
                    compare_exponential(FloatMP::ROUND_TO_NEAREST);
                    compare_exponential(FloatMP::ROUND_DOWNWARD);
                    compare_exponential(FloatMP::ROUND_UPWARD);
                    compare_exponential(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // One-round integer power, with an independent exact oracle:
            // nearest-round the original input, form its exact integer power,
            // then round once in the requested direction. Fixture keys use
            // signed_input so no FloatMP-derived value becomes the oracle.
            FloatMPPhaseTiming timing_integer_pow("integer_pow",precision.bits());
            for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                const Rational magnitudes[] = {
                    Rational(0),Rational(1),Rational(1,2),
                    Rational(3,4),Rational(17,16),Rational(5,3),
                    Rational(3),Rational(131,128)
                };
                for (Rational const& magnitude : magnitudes) {
                    for (Int sign_factor : {Int(-1),Int(1)}) {
                        const Rational signed_input =
                            sign_factor < 0 ? -magnitude : magnitude;
                        const FloatMP ap_input(
                            signed_input,FloatMP::ROUND_TO_NEAREST,mp);
                        for (Nat power : {Nat(0u),Nat(1u),Nat(2u),
                                          Nat(3u),Nat(7u),Nat(8u),
                                          Nat(13u),Nat(31u)}) {
                            const auto check=[&](auto ap_mode) {
                                ARIADNE_TEST_EQUALS(
                                    static_cast<Rational>(
                                        pow(ap_mode,ap_input,power)),
                                    floatmp_power_golden(
                                        mp.bits(),ap_mode,signed_input,power));
                            };
                            check(FloatMP::ROUND_TO_NEAREST);
                            check(FloatMP::ROUND_DOWNWARD);
                            check(FloatMP::ROUND_UPWARD);
                            check(FloatMP::ROUND_TOWARD_ZERO);
                        }
                    }
                }
            }

            // Signed zero, infinities and NaN follow the preexisting
            // IEEE-style special-case rules; 0th power stays one.
            const MultiplePrecision mp(53_bits);
            const FloatMP minus_zero=FloatMP::zero(53ul,true);
            const FloatMP minus_infinity=
                FloatMP::inf(Sign::NEGATIVE,mp);
            const FloatMP plus_infinity=
                FloatMP::inf(Sign::POSITIVE,mp);
            ARIADNE_TEST_ASSERT(
                pow(FloatMP::ROUND_TO_NEAREST,minus_zero,Nat(7u)).negative());
            ARIADNE_TEST_ASSERT(
                not pow(FloatMP::ROUND_TO_NEAREST,minus_zero,Nat(8u)).negative());
            ARIADNE_TEST_ASSERT(
                pow(FloatMP::ROUND_TO_NEAREST,minus_infinity,Nat(7u)).negative());
            ARIADNE_TEST_ASSERT(
                not pow(FloatMP::ROUND_TO_NEAREST,minus_infinity,Nat(8u)).negative());
            ARIADNE_TEST_EQUALS(
                pow(FloatMP::ROUND_TO_NEAREST,plus_infinity,Nat(0u)),
                FloatMP(1,mp));
            ARIADNE_TEST_EQUALS(
                pow(FloatMP::ROUND_TO_NEAREST,FloatMP::nan(mp),Nat(0u)),
                FloatMP(1,mp));
            ARIADNE_TEST_ASSERT(is_nan(
                pow(FloatMP::ROUND_TO_NEAREST,FloatMP::nan(mp),Nat(7u))));
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Certified centered logarithm around +1/2: exactly rounded
            // golden values obtained from rigorous rational atanh series.
            // Preserve all input visits, cutoff probes and rounding modes.
            FloatMPPhaseTiming timing_centered_log("centered_log",precision.bits());
            for (Bits bits : {53_bits,7_bits,128_bits,256_bits,
                              1024_bits,53_bits}) {
                const MultiplePrecision mp(bits);
                const Nat ulp_power(static_cast<Nat>(mp.bits()));
                const Rational ulp(Integer(1),pow(Integer(2),ulp_power));
                const Nat cutoff_power(static_cast<Nat>(
                    (static_cast<unsigned long>(mp.bits())+33ul)/2ul));
                const Rational cutoff(Integer(1),pow(Integer(2),cutoff_power));
                const Rational values[]={
                    Rational(1,2), Rational(1,2)-ulp,
                    Rational(1,2)+ulp,
                    Rational(1,2)-Rational(4)*ulp,
                    Rational(1,2)+Rational(4)*ulp,
                    Rational(1,2)-cutoff,
                    Rational(1,2)+cutoff,
                    Rational(1,2)-Rational(2)*cutoff,
                    Rational(1,2)+Rational(2)*cutoff,
                    Rational(1,4),Rational(3,4),
                    Rational(1,3),Rational(2,3),
                    Rational(1),Rational(2)
                };
                for (Rational const& exact : values) {
                    const FloatMP ap_input(exact,
                        FloatMP::ROUND_TO_NEAREST,mp);
                    const Rational rounded_input=static_cast<Rational>(ap_input);
                    const auto check=[&](auto ap_mode) {
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(log(ap_mode,ap_input)),
                            floatmp_log_golden(mp.bits(),ap_mode,rounded_input));
                    };
                    check(FloatMP::ROUND_TO_NEAREST);
                    check(FloatMP::ROUND_DOWNWARD);
                    check(FloatMP::ROUND_UPWARD);
                    check(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Certified exponential oracle around BOTH +/-1/2.
            // Keep the exact window boundary, one step outside it,
            // all four directions, very low precision, high-precision
            // cache upgrades and revisiting a lower precision.
            FloatMPPhaseTiming timing_centered_exp("centered_exp",precision.bits());
            for (Bits bits : {53_bits,7_bits,128_bits,256_bits,
                              1024_bits,53_bits}) {
                const MultiplePrecision mp(bits);
                const Nat exponent_ulp(static_cast<Nat>(mp.bits()));
                const Rational ulp(Integer(1),
                                   pow(Integer(2),exponent_ulp));
                const Nat cutoff_exponent(static_cast<Nat>(
                    (static_cast<unsigned long>(mp.bits())+33ul)/2ul));
                const Rational cutoff(Integer(1),
                    pow(Integer(2),cutoff_exponent));
                const Rational magnitudes[] = {
                    Rational(1,2), Rational(1,2)-ulp,
                    Rational(1,2)+ulp,
                    Rational(1,2)-Rational(4)*ulp,
                    Rational(1,2)+Rational(4)*ulp,
                    Rational(1,2)-cutoff,
                    Rational(1,2)+cutoff,
                    Rational(1,2)-Rational(2)*cutoff,
                    Rational(1,2)+Rational(2)*cutoff,
                    Rational(1,4), Rational(3,4),
                    Rational(1,3), Rational(2,3),
                    Rational(1), Rational(2)
                };
                for (Rational const& magnitude : magnitudes) {
                    for (Int sign_factor : {Int(-1),Int(1)}) {
                        const Rational exact = sign_factor < 0
                            ? -magnitude : magnitude;
                        const FloatMP ap_input(exact,
                            FloatMP::ROUND_TO_NEAREST,mp);
                        const Rational rounded_input=
                            static_cast<Rational>(ap_input);
                        const auto check=[&](auto ap_mode) {
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(exp(ap_mode,ap_input)),
                                floatmp_centered_exp_golden(
                                    mp.bits(),ap_mode,rounded_input));
                        };
                        check(FloatMP::ROUND_TO_NEAREST);
                        check(FloatMP::ROUND_DOWNWARD);
                        check(FloatMP::ROUND_UPWARD);
                        check(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Certified arcsine/arccosine local Taylor enclosures at
            // both signs of 1/2. Preserve cache growth/reuse and the
            // exact inclusion-threshold and exterior-fallback probes.
            // The golden reference is keyed by the ORIGINAL source Rational.
            FloatMPPhaseTiming timing_centered_asin_acos("centered_asin_acos",precision.bits());
            for (Bits bits : {53_bits,7_bits,128_bits,256_bits,
                              1024_bits,53_bits}) {
                const MultiplePrecision mp(bits);
                const Nat offset_power(static_cast<Nat>(mp.bits()));
                const Rational ulp(Integer(1),
                                   pow(Integer(2),offset_power));
                // Initial work is bits+32: exercise the exact centered
                // inclusion threshold as well as the fallback exterior.
                const Nat threshold_power(
                    (static_cast<unsigned long>(mp.bits())+33ul)/2ul);
                const Rational threshold_offset(
                    Integer(1),pow(Integer(2),threshold_power));
                const Rational magnitudes[] = {
                    Rational(1,2), Rational(1,2)-ulp, Rational(1,2)+ulp,
                    Rational(1,2)-Rational(4)*ulp,
                    Rational(1,2)+Rational(4)*ulp,
                    Rational(1,2)-threshold_offset,
                    Rational(1,2)+threshold_offset,
                    Rational(1,2)-Rational(2)*threshold_offset,
                    Rational(1,2)+Rational(2)*threshold_offset,
                    Rational(1,4), Rational(3,4),
                    Rational(1,3), Rational(2,3)
                };
                for (Rational const& magnitude : magnitudes) {
                    for (Int input_sign : {Int(-1),Int(1)}) {
                        const Rational signed_value =
                            input_sign<0 ? -magnitude : magnitude;
                        const FloatMP ap_input(signed_value,
                            FloatMP::ROUND_TO_NEAREST,mp);
                        const auto check=[&](auto ap_mode) {
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(asin(ap_mode,ap_input)),
                                floatmp_asin_golden(mp.bits(),ap_mode,signed_value));
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(acos(ap_mode,ap_input)),
                                floatmp_acos_golden(mp.bits(),ap_mode,signed_value));
                        };
                        check(FloatMP::ROUND_TO_NEAREST);
                        check(FloatMP::ROUND_DOWNWARD);
                        check(FloatMP::ROUND_UPWARD);
                        check(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Certified atan(1/2) linearization on both sides of the
            // exact work-dependent crossover, including signs, directed
            // rounding, and cache precision upgrades/reuse. Reference
            // keyed by the ORIGINAL input, independently rounded offline.
            FloatMPPhaseTiming timing_centered_atan("centered_atan",precision.bits());
            for (Bits bits : {53_bits,7_bits,128_bits,256_bits,
                              1024_bits,53_bits}) {
                const MultiplePrecision mp(bits);
                const Nat exponent_ulp(static_cast<Nat>(mp.bits()));
                const Rational ulp(Integer(1),
                                   pow(Integer(2),exponent_ulp));
                const Nat cutoff_power(static_cast<Nat>(
                    (static_cast<unsigned long>(mp.bits())+33ul)/2ul));
                const Rational cutoff(Integer(1),
                    pow(Integer(2),cutoff_power));
                const Rational magnitudes[] = {
                    Rational(1,2), Rational(1,2)-ulp,
                    Rational(1,2)+ulp,
                    Rational(1,2)-Rational(4)*ulp,
                    Rational(1,2)+Rational(4)*ulp,
                    Rational(1,2)-cutoff, Rational(1,2)+cutoff,
                    Rational(1,2)-Rational(2)*cutoff,
                    Rational(1,2)+Rational(2)*cutoff,
                    Rational(1,4), Rational(3,4), Rational(5,7),
                    Rational(1,3), Rational(2,3), Rational(1),
                    Rational(2)
                };
                for (Rational const& magnitude : magnitudes) {
                    for (Int sign_factor : {Int(-1),Int(1)}) {
                        const Rational exact =
                            sign_factor < 0 ? -magnitude : magnitude;
                        const FloatMP ap_input(exact,
                            FloatMP::ROUND_TO_NEAREST,mp);
                        const auto check=[&](auto ap_mode) {
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(atan(ap_mode,ap_input)),
                                floatmp_atan_golden(mp.bits(),ap_mode,exact));
                        };
                        check(FloatMP::ROUND_TO_NEAREST);
                        check(FloatMP::ROUND_DOWNWARD);
                        check(FloatMP::ROUND_UPWARD);
                        check(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Single-recurrence atan certification at |x|=1/2 and |x|=2,
            // small arguments, both signs, and all directed modes.
            // Use an independently certified exact binary oracle, keyed
            // by the original Dyadic source before input rounding.
            FloatMPPhaseTiming timing_single_recurrence_atan("single_recurrence_atan",precision.bits());
            for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                for (Dyadic const& argument :
                     {Dyadic(1,64u),Dyadic(-1,64u),
                      Dyadic(255,9u),Dyadic(-255,9u),
                      Dyadic(1,1u),Dyadic(-1,1u),
                      Dyadic(257,9u),Dyadic(-257,9u),
                      Dyadic(1,0u),Dyadic(-1,0u),
                      Dyadic(1023,9u),Dyadic(-1023,9u),
                      Dyadic(2,0u),Dyadic(-2,0u),
                      Dyadic(1025,9u),Dyadic(-1025,9u)}) {
                    const FloatMP a(Rational(argument),FloatMP::ROUND_TO_NEAREST,mp);
                    const auto check_atan=[&](auto ap_mode) {
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(atan(ap_mode,a)),
                            floatmp_atan_golden(
                                mp.bits(),ap_mode,Rational(argument)));
                    };
                    check_atan(FloatMP::ROUND_TO_NEAREST);
                    check_atan(FloatMP::ROUND_DOWNWARD);
                    check_atan(FloatMP::ROUND_UPWARD);
                    check_atan(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Certified atan(1/2) cache and centered argument reduction:
            // probe either side of 1/4, 1/2 and 5/7, both signs,
            // different precision/cache ordering, and all rounding modes.
            // The independent golden key is the ORIGINAL source Rational.
            const Rational epsilon(Dyadic(1,40u));
            const Rational tiny(Dyadic(1,100u));
            const Rational values[]={
                Rational(0),Rational(1,8),
                Rational(1,4)-epsilon,Rational(1,4),Rational(1,4)+epsilon,
                Rational(1,3)-epsilon,Rational(1,3),Rational(1,3)+epsilon,
                Rational(1,2)-epsilon,
                Rational(1,2)-tiny,Rational(1,2),
                Rational(1,2)+tiny,Rational(1,2)+epsilon,
                Rational(5,7)-epsilon,Rational(5,7),Rational(5,7)+epsilon,
                Rational(3,4),Rational(1),Rational(2),Rational(3)
            };
            for (Bits bits : {53_bits,7_bits,128_bits,1024_bits,256_bits,53_bits}) {
                const MultiplePrecision mp(bits);
                for (Rational const& value : values) {
                    for (Int sign : {Int(-1),Int(1)}) {
                        const Rational exact=sign<0 ? -value : value;
                        const FloatMP a(exact,FloatMP::ROUND_TO_NEAREST,mp);
                        const auto verify=[&](auto am) {
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(atan(am,a)),
                                floatmp_atan_golden(mp.bits(),am,exact));
                        };
                        verify(FloatMP::ROUND_TO_NEAREST);
                        verify(FloatMP::ROUND_DOWNWARD);
                        verify(FloatMP::ROUND_UPWARD);
                        verify(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Certified asin/acos one-evaluation and shared enclosure:
            // retain endpoint choices, near +/-1 and zero probes, all
            // rounding modes, and independently rounded source-input keys.
            for (Bits bits : {53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                const Rational epsilon(Dyadic(1,40u));
                const Rational inputs[] = {
                    Rational(0), Rational(1,64), Rational(-1,64),
                    Rational(1,2), Rational(-1,2),
                    Rational(1,2)-epsilon, Rational(1,2)+epsilon,
                    epsilon-Rational(1,2), -Rational(1,2)-epsilon,
                    Rational(3,4), Rational(-3,4),
                    Rational(7,8), Rational(-7,8),
                    Rational(1) - epsilon, epsilon - Rational(1),
                    Rational(1), Rational(-1)
                };
                for (Rational const& input : inputs) {
                    const FloatMP a(input,FloatMP::ROUND_TO_NEAREST,mp);
                    const auto compare_inverse=[&](auto ap_mode) {
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(asin(ap_mode,a)),
                            floatmp_asin_golden(mp.bits(),ap_mode,input));
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(acos(ap_mode,a)),
                            floatmp_acos_golden(mp.bits(),ap_mode,input));
                    };
                    compare_inverse(FloatMP::ROUND_TO_NEAREST);
                    compare_inverse(FloatMP::ROUND_DOWNWARD);
                    compare_inverse(FloatMP::ROUND_UPWARD);
                    compare_inverse(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Certified sine/cosine around both +/-1/2 with cache/guard
            // threshold crossings at 96/128/192/224/256/1024 bits.
            // Revisit 128 bits for cache reuse, comparing the same test
            // calls to independently certified original-input fixtures.
            FloatMPPhaseTiming timing_centered_sin_cos("centered_sin_cos",precision.bits());
            for (Bits bits : {96_bits,128_bits,192_bits,224_bits,
                              256_bits,1024_bits,128_bits}) {
                const MultiplePrecision mp(bits);
                const Nat offset_power(static_cast<unsigned long>(mp.bits())/2ul+16ul);
                const Rational small_offset(Integer(1),pow(Integer(2),offset_power));
                for (Rational const& magnitude : {
                        Rational(1,2),Rational(1,2)-small_offset,
                        Rational(1,2)+small_offset,
                        Rational(1,2)-Rational(2)*small_offset,
                        Rational(1,2)+Rational(2)*small_offset,
                        Rational(1,4),Rational(3,4)}) {
                    for (Int sign : {Int(-1),Int(1)}) {
                        const Rational exact =
                            sign < 0 ? -magnitude : magnitude;
                        const FloatMP a(exact,
                                        FloatMP::ROUND_TO_NEAREST,mp);
                        const auto verify=[&](auto am) {
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(sin(am,a)),
                                floatmp_sin_golden(mp.bits(),am,exact));
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(cos(am,a)),
                                floatmp_cos_golden(mp.bits(),am,exact));
                        };
                        verify(FloatMP::ROUND_TO_NEAREST);
                        verify(FloatMP::ROUND_DOWNWARD);
                        verify(FloatMP::ROUND_UPWARD);
                        verify(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Local certified sin/cos/tan near +/-1/2, including
            // exact centers, small offsets and fallback arguments.
            // Preserve all rounding modes and the work>=256 dispatch
            // transition; use exact source-keyed independent oracles.
            FloatMPPhaseTiming timing_local_sin_cos("local_sin_cos",precision.bits());
            for (Bits bits : {53_bits,7_bits,128_bits,192_bits,224_bits,
                              256_bits,384_bits,1024_bits,53_bits}) {
                const MultiplePrecision mp(bits);
                const Rational unit(Integer(1),
                    pow(Integer(2),Nat(mp.bits())));
                const Rational half(1,2);
                const Rational nearby[]={
                    half, half-unit, half+unit,
                    half-Rational(8)*unit, half+Rational(8)*unit,
                    Rational(3,8),Rational(5,8),Rational(1,4),
                    Rational(3,4)
                };
                for (Rational const& positive : nearby) {
                    for (Int sign : {Int(-1),Int(1)}) {
                        const Rational value=sign<0 ? -positive : positive;
                        const FloatMP a(value,FloatMP::ROUND_TO_NEAREST,mp);
                        const auto verify=[&](auto am) {
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(sin(am,a)),
                                floatmp_sin_golden(mp.bits(),am,value));
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(cos(am,a)),
                                floatmp_cos_golden(mp.bits(),am,value));
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(tan(am,a)),
                                floatmp_tan_golden(mp.bits(),am,value));
                        };
                        verify(FloatMP::ROUND_TO_NEAREST);
                        verify(FloatMP::ROUND_DOWNWARD);
                        verify(FloatMP::ROUND_UPWARD);
                        verify(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Tangent-centered enclosure at +/-1/2. Preserve precision
            // ordering, cache reuse and upgrades, and both sides of the
            // exact |delta| <= 2^-ceil(work/2) threshold. Compare against
            // certified independently rounded source-keyed golden values.
            FloatMPPhaseTiming timing_centered_tan("centered_tan",precision.bits());
            for (Bits bits : {53_bits,7_bits,256_bits,128_bits,1024_bits,53_bits}) {
                const MultiplePrecision mp(bits);
                const Nat offset_denominator_power(static_cast<unsigned long>(mp.bits()) / 2ul + 16ul);
                const Rational step(Integer(1),pow(Integer(2),offset_denominator_power));
                const Rational half(1,2);
                const Rational values[]={
                    half, half-step, half+step,
                    half-Rational(2)*step, half+Rational(2)*step,
                    half-Rational(8)*step, half+Rational(8)*step,
                    Rational(1,4), Rational(3,4)
                };
                for (Rational const& magnitude : values) {
                    for (Int sign : {Int(-1),Int(1)}) {
                        const Rational exact = sign < 0 ? -magnitude : magnitude;
                        const FloatMP a(exact,FloatMP::ROUND_TO_NEAREST,mp);
                        const auto verify=[&](auto am) {
                            ARIADNE_TEST_EQUALS(
                                static_cast<Rational>(tan(am,a)),
                                floatmp_tan_golden(mp.bits(),am,exact));
                        };
                        verify(FloatMP::ROUND_TO_NEAREST);
                        verify(FloatMP::ROUND_DOWNWARD);
                        verify(FloatMP::ROUND_UPWARD);
                        verify(FloatMP::ROUND_TOWARD_ZERO);
                    }
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // The single-sqrt cosine enclosure is used by tangent for
            // |x|<=1. Preserve both signs, low precision and reconstruction
            // across |x|=1, with certified exact source-keyed references.
            FloatMPPhaseTiming timing_tangent_sqrt_cos("tangent_sqrt_cos",precision.bits());
            for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                for (Dyadic const& argument :
                     {Dyadic(1,64u),Dyadic(-1,64u),
                      Dyadic(1,1u),Dyadic(-1,1u),
                      Dyadic(3,2u),Dyadic(-3,2u),
                      Dyadic(1023,10u),Dyadic(-1023,10u),
                      Dyadic(1,0u),Dyadic(-1,0u),
                      Dyadic(1025,10u),Dyadic(-1025,10u),
                      Dyadic(3,1u),Dyadic(-3,1u)}) {
                    // Some source dyadics need more than seven significant
                    // bits; retain nearest/even input rounding instead of
                    // using the exact Dyadic constructor at that precision.
                    const Rational exact(argument);
                    const FloatMP a(exact,FloatMP::ROUND_TO_NEAREST,mp);
                    const auto check_tangent=[&](auto ap_mode) {
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(tan(ap_mode,a)),
                            floatmp_tan_golden(mp.bits(),ap_mode,exact));
                    };
                    check_tangent(FloatMP::ROUND_TO_NEAREST);
                    check_tangent(FloatMP::ROUND_DOWNWARD);
                    check_tangent(FloatMP::ROUND_UPWARD);
                    check_tangent(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Single-recurrence sine/cosine/tangent rounding (including
            // cosine reconstructed from sine). Keep the exact dyadic
            // inputs, both sides of |x|=1 and the very small arguments,
            // using independent correctly rounded rational references.
            for (Bits bits : {53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                for (Dyadic const& argument :
                     {Dyadic(1,2u),Dyadic(-1,2u),Dyadic(3,2u),
                      Dyadic(-3,2u),Dyadic(5,1u),Dyadic(-5,1u),
                      Dyadic(7,1u),Dyadic(-7,1u),Dyadic(2,1u),Dyadic(-2,1u),
                      Dyadic(7,3u),Dyadic(-7,3u),
                      Dyadic(1023,10u),Dyadic(-1023,10u),
                      Dyadic(31,5u),Dyadic(-31,5u),
                      Dyadic(17,4u),Dyadic(-17,4u),
                      Dyadic(1,64u),Dyadic(-1,64u),
                      Dyadic(1,128u),Dyadic(-1,128u)}) {
                    const FloatMP a(argument,mp);
                    const Rational source(argument);
                    const auto compare_trigonometric=[&](auto ap_mode) {
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(sin(ap_mode,a)),
                            floatmp_sin_golden(mp.bits(),ap_mode,source));
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(cos(ap_mode,a)),
                            floatmp_cos_golden(mp.bits(),ap_mode,source));
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(tan(ap_mode,a)),
                            floatmp_tan_golden(mp.bits(),ap_mode,source));
                    };
                    compare_trigonometric(FloatMP::ROUND_TO_NEAREST);
                    compare_trigonometric(FloatMP::ROUND_DOWNWARD);
                    compare_trigonometric(FloatMP::ROUND_UPWARD);
                    compare_trigonometric(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Check log range-reduction boundaries, cancellation near 1,
            // cached log(2) enclosures and all directed rounding modes
            // against independently certified original-input golden values.
            FloatMPPhaseTiming timing_log_boundaries("log_boundaries",precision.bits());
            {
                const MultiplePrecision cancellation_precision(128_bits);
                const Rational epsilon(Dyadic(1,98u));
                for(const auto& input:{Rational(1)-epsilon,Rational(1)+epsilon}) {
                    const FloatMP value(input,near,cancellation_precision);
                    for(auto mode:{FloatMP::ROUND_TO_NEAREST,FloatMP::ROUND_DOWNWARD,
                                   FloatMP::ROUND_UPWARD,FloatMP::ROUND_TOWARD_ZERO}) {
                        ARIADNE_TEST_EQUALS(static_cast<Rational>(log(mode,value)),
                                            floatmp_log_golden(128u,mode,input));
                    }
                }
            }

            for (Bits bits : {53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                const Rational epsilon(Dyadic(1,40u));
                const Rational inputs[] = {
                    Rational(1,2), Rational(3,4), Rational(3,2),
                    Rational(2), Rational(1,8), Rational(8),
                    Rational(1) - epsilon, Rational(1) + epsilon,
                    Rational(3,4) - epsilon, Rational(3,4) + epsilon,
                    Rational(3,2) - epsilon, Rational(3,2) + epsilon,
                    Rational(1,2) + epsilon
                };
                for (Rational const& input : inputs) {
                    const FloatMP ap_input(input,FloatMP::ROUND_TO_NEAREST,mp);
                    const auto check=[&](auto ap_mode) {
                        ARIADNE_TEST_EQUALS(
                            static_cast<Rational>(log(ap_mode,ap_input)),
                            floatmp_log_golden(mp.bits(),ap_mode,input));
                    };
                    check(FloatMP::ROUND_TO_NEAREST);
                    check(FloatMP::ROUND_DOWNWARD);
                    check(FloatMP::ROUND_UPWARD);
                    check(FloatMP::ROUND_TOWARD_ZERO);
                }
            }
        }

        // Transcendental directed-rounding checks at multiple precisions.
        // Exact mathematical identities supplement these enclosures.
        for (Bits bits : {7_bits,24_bits,53_bits,113_bits}) {
            const PR pr(bits);
            const F quarter(Dyadic(1,2u),pr);
            const F negative_quarter(Dyadic(-1,2u),pr);
            const F one(1,pr);
            const F two_value(2,pr);
            const auto enclosure=[&](auto operation,F const& x) {
                const F lower=operation(F::ROUND_DOWNWARD,x);
                const F upper=operation(F::ROUND_UPWARD,x);
                const F nearest=operation(F::ROUND_TO_NEAREST,x);
                const F toward_zero_result=operation(F::ROUND_TOWARD_ZERO,x);
                ARIADNE_TEST_ASSERT(lower<=nearest && nearest<=upper);
                ARIADNE_TEST_ASSERT(lower<=toward_zero_result && toward_zero_result<=upper);
            };
            enclosure([](auto mode,auto const& x){return sqrt(mode,x);},two_value);
            enclosure([](auto mode,auto const& x){return exp(mode,x);},quarter);
            enclosure([](auto mode,auto const& x){return exp(mode,x);},negative_quarter);
            enclosure([](auto mode,auto const& x){return log(mode,x);},two_value);
            enclosure([](auto mode,auto const& x){return sin(mode,x);},quarter);
            enclosure([](auto mode,auto const& x){return cos(mode,x);},negative_quarter);
            enclosure([](auto mode,auto const& x){return tan(mode,x);},quarter);
            enclosure([](auto mode,auto const& x){return atan(mode,x);},negative_quarter);
            enclosure([](auto mode,auto const& x){return asin(mode,x);},quarter);
            enclosure([](auto mode,auto const& x){return acos(mode,x);},quarter);
            ARIADNE_TEST_ASSERT(F::pi(F::ROUND_DOWNWARD,pr)<=
                                F::pi(F::ROUND_TO_NEAREST,pr));
            ARIADNE_TEST_ASSERT(F::pi(F::ROUND_TO_NEAREST,pr)<=
                                F::pi(F::ROUND_UPWARD,pr));
            ARIADNE_TEST_EQUALS(exp(F::ROUND_TO_NEAREST,F(0,pr)),one);
            ARIADNE_TEST_EQUALS(log(F::ROUND_TO_NEAREST,one),F(0,pr));
            ARIADNE_TEST_EQUALS(sin(F::ROUND_TO_NEAREST,F(0,pr)),F(0,pr));
            ARIADNE_TEST_EQUALS(cos(F::ROUND_TO_NEAREST,F(0,pr)),one);
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // Exercise cold and warm pi-cache paths, switching precisions
            // and rounding modes. pi = 4*atan(1): multiplication by the
            // exact power of two preserves p-bit correct rounding. The atan
            // source-keyed fixture was certified independently of FloatMP
            // (and its scaled pi values also checked with a Machin bound).
            for (Bits bits : {53_bits,256_bits,7_bits,1024_bits,128_bits,53_bits}) {
                const MultiplePrecision mp(bits);
                const auto check_pi=[&](auto ap_mode) {
                    const FloatMP first=FloatMP::pi(ap_mode,mp);
                    const FloatMP repeated=FloatMP::pi(ap_mode,mp);
                    const Rational reference=Rational(4) *
                        floatmp_atan_golden(mp.bits(),ap_mode,Rational(1));
                    ARIADNE_TEST_EQUALS(static_cast<Rational>(first),
                                        reference);
                    ARIADNE_TEST_SAME(first,repeated);
                };
                check_pi(FloatMP::ROUND_TO_NEAREST);
                check_pi(FloatMP::ROUND_UPWARD);
                check_pi(FloatMP::ROUND_DOWNWARD);
                check_pi(FloatMP::ROUND_TOWARD_ZERO);
                check_pi(FloatMP::ROUND_TO_NEAREST);
            }
        }

        if constexpr (Same<PR,MultiplePrecision>) {
            // FloatMP::get_d must implement correctly rounded IEEE 754
            // binary64 conversion from each stored p-bit value, including
            // midpoint ties, signed subnormals, overflows and signed zero.
            // Frozen bit-pattern references are calculated independently
            // from the ORIGINAL exact Rational and p-bit input rounding.
            const auto previous_ap=FloatMP::get_rounding_mode();
            for (Bits bits : {53_bits,128_bits,256_bits,1024_bits}) {
                const MultiplePrecision mp(bits);
                const FloatMP one(1,mp);
                const FloatMP negative_one(-1,mp);
                const FloatMP tie_even(Integer(String("9007199254740993")),mp);
                const FloatMP tie_odd(Integer(String("9007199254740995")),mp);
                const FloatMP near_overflow(
                    Integer(String("18014398509481983")),mp);
                const FloatMP values[]={
                    one,negative_one,
                    FloatMP(Rational(1,3),FloatMP::ROUND_TO_NEAREST,mp),
                    FloatMP(Rational(-17,13),FloatMP::ROUND_TO_NEAREST,mp),
                    tie_even,tie_odd,
                    shft(one,Int(-1200)),shft(negative_one,Int(-1200)),
                    shft(one,Int(-1076)),shft(one,Int(-1075)),
                    shft(negative_one,Int(-1075)),
                    shft(one,Int(-1074)),shft(one,Int(-1022)),
                    shft(one,Int(1023)),shft(one,Int(1024)),
                    shft(negative_one,Int(1024)),
                    shft(near_overflow,Int(970))
                };
                const Rational original_inputs[]={
                    Rational(1),Rational(-1),Rational(1,3),Rational(-17,13),
                    Rational(Integer(String("9007199254740993"))),
                    Rational(Integer(String("9007199254740995"))),
                    Rational(Dyadic(1,1200u)),Rational(Dyadic(-1,1200u)),
                    Rational(Dyadic(1,1076u)),Rational(Dyadic(1,1075u)),
                    Rational(Dyadic(-1,1075u)),
                    Rational(Dyadic(1,1074u)),Rational(Dyadic(1,1022u)),
                    Rational(Dyadic(TwoExp(static_cast<Int>(1023)))),
                    Rational(Dyadic(TwoExp(static_cast<Int>(1024)))),
                    -Rational(Dyadic(TwoExp(static_cast<Int>(1024)))),
                    Rational(Integer(String("18014398509481983"))) *
                        Rational(Dyadic(TwoExp(static_cast<Int>(970))))
                };
                static_assert(std::size(values)==std::size(original_inputs));
                const auto check=[&](auto ap_mode) {
                    FloatMP::set_rounding_mode(ap_mode);
                    for (std::size_t i=0u;i<std::size(values);++i) {
                        const Rational actual_bits{Integer(String(std::to_string(
                            std::bit_cast<std::uint64_t>(values[i].get_d()))))};
                        ARIADNE_TEST_EQUALS(actual_bits,
                            floatmp_get_d_bits_golden(
                                mp.bits(),ap_mode,original_inputs[i]));
                    }
                };
                check(FloatMP::ROUND_TO_NEAREST);
                check(FloatMP::ROUND_DOWNWARD);
                check(FloatMP::ROUND_UPWARD);
                check(FloatMP::ROUND_TOWARD_ZERO);

                const FloatMP negative_zero=FloatMP::zero(
                    static_cast<unsigned long>(mp.bits()),true);
                const FloatMP positive_zero=FloatMP::zero(
                    static_cast<unsigned long>(mp.bits()));
                ARIADNE_TEST_ASSERT(std::signbit(negative_zero.get_d()));
                ARIADNE_TEST_ASSERT(not std::signbit(positive_zero.get_d()));
                ARIADNE_TEST_ASSERT(std::isnan(FloatMP::nan(mp).get_d()));
                ARIADNE_TEST_ASSERT(std::isinf(FloatMP::inf(mp).get_d()));
            }
            FloatMP::set_rounding_mode(previous_ap);
        }

        // A deterministic comparison oracle with varied binary scales.
        unsigned long state=0xC0A4E5u;
        for (Bits bits : {7_bits,53_bits,128_bits,256_bits,1024_bits}) {
            const PR pr(bits);
            for (unsigned i=0;i<300u;++i) {
                state=state*1664525u+1013904223u;
                const long na=static_cast<long>(state%200001u)-100000L;
                state=state*1664525u+1013904223u;
                const long nb=static_cast<long>(state%200001u)-100000L;
                const unsigned sa=(i*17u)%200u;
                const unsigned sb=(i*29u)%200u;
                const Rational a=Rational(na)/Rational(Dyadic(TwoExp(static_cast<Int>(sa))));
                const Rational b=Rational(nb)/Rational(Dyadic(TwoExp(static_cast<Int>(sb))));
                const F fa(a,F::ROUND_TO_NEAREST,pr);
                const F fb(b,F::ROUND_TO_NEAREST,pr);
                const Rational exact_a=static_cast<Rational>(fa);
                const Rational exact_b=static_cast<Rational>(fb);
                ARIADNE_TEST_EQUALS(cmp(fa,fb),cmp(exact_a,exact_b));
            }
        }

        const PR pr(53_bits);
        const F positive_zero(0,pr);
        const F negative_zero=neg(positive_zero);
        ARIADNE_TEST_EQUALS(cmp(positive_zero,negative_zero),Comparison::EQUAL);
        ARIADNE_TEST_ASSERT(is_nan(add(near,F::inf(Sign::POSITIVE,pr),
                                       F::inf(Sign::NEGATIVE,pr))));
        ARIADNE_TEST_ASSERT(is_nan(mul(near,F(0,pr),
                                       F::inf(Sign::POSITIVE,pr))));
        ARIADNE_TEST_ASSERT(is_nan(div(near,F(0,pr),F(0,pr))));
    }
}


template<class PR> Void
TestFloatRounding<PR>::test_stream()
{
    cout << ARIADNE_PRETTY_FUNCTION << endl;

    stringstream ss("1.25 -2.25 42 2.375e1 2.35e1");
    Float f1(precision),f2(precision),f3(precision),f4(precision),f5(precision);
    ss >> f1;
    cout << f1 << endl;
    ARIADNE_TEST_EQUALS(f1,1.25_x);
    ss >> f2;
    cout << f2 << endl;
    ARIADNE_TEST_EQUALS(f2,-2.25_x);
    ss >> f3;
    cout << f3 << endl;
    ARIADNE_TEST_EQUALS(f3,42.0_x);
    try {
        ss >> f4;
        cout << f4 << endl;
        if(f4!=23.75_x) {
            ARIADNE_TEST_WARN("Cannot create Float<"<<class_name<PR>()<<"> from string literal in exponential form 2.375e1");
        }
    }
    catch(const std::exception& e) {
        cerr << e.what() << endl;
    }

    try {
        ss >> f5;
        cout << f5 << endl;
        if(f4!=23.5_x) {
            ARIADNE_TEST_WARN("Cannot create Float<"<<class_name<PR>()<<"> from string literal in exponential form 2.35e1");
        }
    }
    catch(const std::exception& e) {
        cerr << e.what() << endl;
    }

    stringstream delimited("1.25,2.5");
    Float fd(precision);
    char delimiter='\0';
    delimited >> fd;
    ARIADNE_TEST_EQUALS(fd,1.25_x);
    delimited >> delimiter;
    ARIADNE_TEST_EQUAL(delimiter,',');

    stringstream terminal("2");
    Float terminal_value(precision);
    terminal >> terminal_value;
    ARIADNE_TEST_EQUALS(terminal_value,2);
    ARIADNE_TEST_ASSERT(terminal.eof());
    ARIADNE_TEST_ASSERT(not terminal.fail());

    if constexpr (std::is_same_v<PR,DoublePrecision>) {
        typename Float::RoundingModeType old_rnd=Float::get_rounding_mode();
        auto test_rounding_independent_parse = [&]() {
            stringstream input_stream("1.25");
            Float fr(precision);
            input_stream >> fr;
            ARIADNE_TEST_EQUALS(fr,1.25_x);
        };

        Float::set_rounding_upward();
        test_rounding_independent_parse();
        Float::set_rounding_downward();
        test_rounding_independent_parse();
        Float::set_rounding_toward_zero();
        test_rounding_independent_parse();
        Float::set_rounding_to_nearest();
        test_rounding_independent_parse();

        Float::set_rounding_mode(old_rnd);
    }

}


template<class PR> Void
TestFloatRounding<PR>::test_comparison()
{
    cout << ARIADNE_PRETTY_FUNCTION << endl;

    ARIADNE_TEST_CONSTRUCT(Float,f0,(3.00_x,precision));
    ARIADNE_TEST_CONSTRUCT(Float,f1,(1.25_x,precision));
    ARIADNE_TEST_CONSTRUCT(Float,f2,(-1.25_x,precision));
    ARIADNE_TEST_CONSTRUCT(Float,f3,(-2.25_x,precision));
    ARIADNE_TEST_CONSTRUCT(Float,f4,(1.25_x,precision));

    // Test comparison of two equal numbers
    ARIADNE_TEST_ASSERT(f1==f4); ARIADNE_TEST_ASSERT(!(f1!=f4));
    ARIADNE_TEST_ASSERT(f1<=f4); ARIADNE_TEST_ASSERT(!(f1> f4));
    ARIADNE_TEST_ASSERT(f1>=f4); ARIADNE_TEST_ASSERT(!(f1< f4));


    // Test comparison of two different numbers
    ARIADNE_TEST_ASSERT(!(f1==f2)); ARIADNE_TEST_ASSERT(f1!=f2);
    ARIADNE_TEST_ASSERT(!(f1<=f2)); ARIADNE_TEST_ASSERT(f1> f2);
    ARIADNE_TEST_ASSERT(f1>=f2); ARIADNE_TEST_ASSERT(!(f1< f2));

    // Test comparison of two negative numbers
    ARIADNE_TEST_ASSERT(!(f2==f3)); ARIADNE_TEST_ASSERT(f2!=f3);
    ARIADNE_TEST_ASSERT(!(f2<=f3)); ARIADNE_TEST_ASSERT(f2> f3);
    ARIADNE_TEST_ASSERT(f2>=f3); ARIADNE_TEST_ASSERT(!(f2< f3));

    // Test comparison with in int
    Int i2=1;
    ARIADNE_TEST_ASSERT(!(f1==i2)); ARIADNE_TEST_ASSERT(f1!=i2);
    ARIADNE_TEST_ASSERT(!(f1<=i2)); ARIADNE_TEST_ASSERT(f1> i2);
    ARIADNE_TEST_ASSERT(f1>=i2); ARIADNE_TEST_ASSERT(!(f1< i2));

    Int i1=1;
    ARIADNE_TEST_ASSERT(!(i1==f2)); ARIADNE_TEST_ASSERT(i1!=f2);
    ARIADNE_TEST_ASSERT(!(i1<=f2)); ARIADNE_TEST_ASSERT(i1> f2);
    ARIADNE_TEST_ASSERT(i1>=f2); ARIADNE_TEST_ASSERT(!(i1< f2));

    // Test comparison with a double
    ExactDouble x2=1.0_x;
    ARIADNE_TEST_ASSERT(!(f1==x2)); ARIADNE_TEST_ASSERT(f1!=x2);
    ARIADNE_TEST_ASSERT(!(f1<=x2)); ARIADNE_TEST_ASSERT(f1> x2);
    ARIADNE_TEST_ASSERT(f1>=x2); ARIADNE_TEST_ASSERT(!(f1< x2));

    ExactDouble x1=1.0_x;
    ARIADNE_TEST_ASSERT(!(x1==f2)); ARIADNE_TEST_ASSERT(x1!=f2);
    ARIADNE_TEST_ASSERT(!(x1<=f2)); ARIADNE_TEST_ASSERT(x1> f2);
    ARIADNE_TEST_ASSERT(x1>=f2); ARIADNE_TEST_ASSERT(!(x1< f2));

    // Test comparison with an integer
    ARIADNE_TEST_CONSTRUCT(Integer,z0,(-3));
    ARIADNE_TEST_ASSERT(!(f1==z0)); ARIADNE_TEST_ASSERT(f1!=z0);
    ARIADNE_TEST_ASSERT(!(f1<=z0)); ARIADNE_TEST_ASSERT(f1> z0);
    ARIADNE_TEST_ASSERT(f1>=z0); ARIADNE_TEST_ASSERT(!(f1< z0));
    ARIADNE_TEST_ASSERT(!(z0==f1)); ARIADNE_TEST_ASSERT(z0<f1);

    // Test comparison with a dyadic
    ARIADNE_TEST_CONSTRUCT(Dyadic,w2,(-5,2u));
    ARIADNE_TEST_ASSERT(!(f1==w2)); ARIADNE_TEST_ASSERT(f1!=w2);
    ARIADNE_TEST_ASSERT(!(w2==f1)); ARIADNE_TEST_ASSERT(w2<f1);
    ARIADNE_TEST_ASSERT(!(f1<=w2)); ARIADNE_TEST_ASSERT(f1> w2);
    ARIADNE_TEST_ASSERT(f1>=w2); ARIADNE_TEST_ASSERT(!(f1< w2));
    ARIADNE_TEST_ASSERT(f2==w2); ARIADNE_TEST_ASSERT(!(f2!=w2));

    // Test comparison with a rational
    ARIADNE_TEST_CONSTRUCT(Rational,q2,(-5,4));
    ARIADNE_TEST_ASSERT(!(f1==q2)); ARIADNE_TEST_ASSERT(f1!=q2);
    ARIADNE_TEST_ASSERT(!(f1<=q2)); ARIADNE_TEST_ASSERT(f1> q2);
    ARIADNE_TEST_ASSERT(f1>=q2); ARIADNE_TEST_ASSERT(!(f1< q2));
    ARIADNE_TEST_ASSERT(f2==q2); ARIADNE_TEST_ASSERT(!(f2!=q2));

    if constexpr (Same<PR,MultiplePrecision>) {
        // Exercise the public Rational-left comparison API in all three
        // ordered outcomes, not merely the FloatMP-left overloads.
        const Rational equal(-5,4);
        const Rational lower(-9,4);
        const Rational higher(5,4);

        ARIADNE_TEST_ASSERT(equal==f2);
        ARIADNE_TEST_ASSERT(!(equal!=f2));
        ARIADNE_TEST_ASSERT(!(equal< f2));
        ARIADNE_TEST_ASSERT(!(equal> f2));
        ARIADNE_TEST_ASSERT(equal<=f2);
        ARIADNE_TEST_ASSERT(equal>=f2);

        ARIADNE_TEST_ASSERT(!(lower==f2));
        ARIADNE_TEST_ASSERT(lower!=f2);
        ARIADNE_TEST_ASSERT(lower< f2);
        ARIADNE_TEST_ASSERT(!(lower> f2));
        ARIADNE_TEST_ASSERT(lower<=f2);
        ARIADNE_TEST_ASSERT(!(lower>=f2));

        ARIADNE_TEST_ASSERT(!(higher==f2));
        ARIADNE_TEST_ASSERT(higher!=f2);
        ARIADNE_TEST_ASSERT(!(higher< f2));
        ARIADNE_TEST_ASSERT(higher> f2);
        ARIADNE_TEST_ASSERT(!(higher<=f2));
        ARIADNE_TEST_ASSERT(higher>=f2);

        // Both directions must reject unordered operands before trying
        // to materialize a nonfinite Rational from the FloatMP value.
        const FloatMP ap_nan=FloatMP::nan(precision);
        const Rational rational_nan=Rational::nan();
        ARIADNE_TEST_ASSERT(cmp(equal,ap_nan)==Comparison::INCOMPARABLE);
        ARIADNE_TEST_ASSERT(cmp(rational_nan,f2)==Comparison::INCOMPARABLE);
        ARIADNE_TEST_ASSERT(cmp(ap_nan,equal)==Comparison::INCOMPARABLE);
        ARIADNE_TEST_ASSERT(cmp(f2,rational_nan)==Comparison::INCOMPARABLE);

        // Compare FloatMP and FloatDP directly in both operand orders.
        // -1.25 and the two neighbors used below are exactly representable
        // by both backends, so expected order does not depend on rounding.
        const FloatDP dp_equal(-1.25_x,dp);
        const FloatDP dp_lower(-2.25_x,dp);
        const FloatDP dp_higher(1.25_x,dp);

        // FloatMP on the left: equal, less, greater.
        ARIADNE_TEST_ASSERT(f2==dp_equal);
        ARIADNE_TEST_ASSERT(!(f2!=dp_equal));
        ARIADNE_TEST_ASSERT(!(f2<dp_equal));
        ARIADNE_TEST_ASSERT(!(f2>dp_equal));
        ARIADNE_TEST_ASSERT(f2<=dp_equal);
        ARIADNE_TEST_ASSERT(f2>=dp_equal);

        ARIADNE_TEST_ASSERT(!(f2==dp_lower));
        ARIADNE_TEST_ASSERT(f2!=dp_lower);
        ARIADNE_TEST_ASSERT(!(f2<dp_lower));
        ARIADNE_TEST_ASSERT(f2>dp_lower);
        ARIADNE_TEST_ASSERT(!(f2<=dp_lower));
        ARIADNE_TEST_ASSERT(f2>=dp_lower);

        ARIADNE_TEST_ASSERT(!(f2==dp_higher));
        ARIADNE_TEST_ASSERT(f2!=dp_higher);
        ARIADNE_TEST_ASSERT(f2<dp_higher);
        ARIADNE_TEST_ASSERT(!(f2>dp_higher));
        ARIADNE_TEST_ASSERT(f2<=dp_higher);
        ARIADNE_TEST_ASSERT(!(f2>=dp_higher));

        // FloatDP on the left: exercise the distinct mixed overloads.
        ARIADNE_TEST_ASSERT(dp_equal==f2);
        ARIADNE_TEST_ASSERT(!(dp_equal!=f2));
        ARIADNE_TEST_ASSERT(!(dp_equal<f2));
        ARIADNE_TEST_ASSERT(!(dp_equal>f2));
        ARIADNE_TEST_ASSERT(dp_equal<=f2);
        ARIADNE_TEST_ASSERT(dp_equal>=f2);

        ARIADNE_TEST_ASSERT(!(dp_lower==f2));
        ARIADNE_TEST_ASSERT(dp_lower!=f2);
        ARIADNE_TEST_ASSERT(dp_lower<f2);
        ARIADNE_TEST_ASSERT(!(dp_lower>f2));
        ARIADNE_TEST_ASSERT(dp_lower<=f2);
        ARIADNE_TEST_ASSERT(!(dp_lower>=f2));

        ARIADNE_TEST_ASSERT(!(dp_higher==f2));
        ARIADNE_TEST_ASSERT(dp_higher!=f2);
        ARIADNE_TEST_ASSERT(!(dp_higher<f2));
        ARIADNE_TEST_ASSERT(dp_higher>f2);
        ARIADNE_TEST_ASSERT(!(dp_higher<=f2));
        ARIADNE_TEST_ASSERT(dp_higher>=f2);
    }

}

template<class PR> Void
TestFloatRounding<PR>::test_rounding_mode()
{
    Float three(3,precision);
    Float third_down(precision), third_near(precision), third_up(precision), third_rnd(precision);

    Float::set_rounding_mode(downward);
    third_up=rec(up,three);
    third_rnd=rec(rounded,three);
    if (third_rnd==third_up) {
        ARIADNE_TEST_NOTIFY("Performing upward-rounded operation changes rounding mode for future operations.");
    } else {
        ARIADNE_TEST_NOTIFY("Performing upward-rounded operation preserves rounding mode for future operations.");
    }

    Float::set_rounding_mode(upward);
    third_down=rec(down,three);
    third_rnd=rec(rounded,three);
    if (third_rnd==third_down) {
        ARIADNE_TEST_NOTIFY("Performing downward-rounded operation changes rounding mode for future operations.");
    } else {
        ARIADNE_TEST_NOTIFY("Performing downward-rounded operation preserves rounding mode for future operations.");
    }

    third_near=rec(near,three);
    if(third_near==third_down) {
        Float::set_rounding_mode(upward);
    } else {
        Float::set_rounding_mode(downward);
    }
    third_near=rec(near,three);
    third_rnd=rec(rounded,three);
    if (third_rnd==third_near) {
        ARIADNE_TEST_NOTIFY("Performing to-nearest-rounded operation changes rounding mode for future operations.");
    } else {
        ARIADNE_TEST_NOTIFY("Performing to-nearest-rounded operation preserves rounding mode for future operations.");
    }
}


template<class PR> Void
TestFloatRounding<PR>::test_rounding_session()
{
    cout << ARIADNE_PRETTY_FUNCTION << endl;

    using RoundingModeType=typename Float::RoundingModeType;

    const RoundingModeType original=Float::get_rounding_mode();
    const RoundingModeType down_mode=static_cast<RoundingModeType>(down);
    const RoundingModeType near_mode=static_cast<RoundingModeType>(near);
    const RoundingModeType up_mode=static_cast<RoundingModeType>(up);

    const Float one(1,precision);
    const Float three(3,precision);
    const Float third_near=div(near,one,three);
    const Float third_up=div(up,one,three);
    const Float third_down=div(down,one,three);

    Float::set_rounding_mode(down_mode);
    {
        RoundingSession<Float> session;
        ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==down_mode);

        CurrentRoundingMode current=session.set(near);
        ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==near_mode);
        ARIADNE_TEST_EQUAL(div(current,one,three),third_near);

        current=session.set(near);
        ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==near_mode);
        ARIADNE_TEST_EQUAL(div(current,one,three),third_near);

        current=session.set(up);
        ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==up_mode);
        ARIADNE_TEST_EQUAL(div(current,one,three),third_up);
    }
    ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==down_mode);

    {
        RoundingSession<Float> session;
        CurrentRoundingMode current=session.set(down);
        ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==down_mode);
        ARIADNE_TEST_EQUAL(div(current,one,three),third_down);
    }
    ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==down_mode);

    {
        RoundingSession<Float> outer;
        outer.set(up);
        ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==up_mode);
        {
            RoundingSession<Float> inner;
            inner.set(near);
            ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==near_mode);
        }
        ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==up_mode);
    }
    ARIADNE_TEST_ASSERT(Float::get_rounding_mode()==down_mode);

    Float::set_rounding_mode(original);
}


template<> Void
TestFloatRounding<DoublePrecision>::test_double_rounding()
{
    volatile double one   = 1;
    volatile double two_value = 2;
    volatile double three = 3;
    volatile double five  = 5;
    const double onethirddown    = 0.33333333333333331483;
    const double onethirdup      = 0.33333333333333337034;
    const double onethirdchop    = 0.33333333333333331483;
    const double onethirdnearest = 0.33333333333333331483;
    const double twofifthsdown   = 0.39999999999999996669;
    const double twofifthsup     = 0.40000000000000002220;
    const double twofifthschop   = 0.39999999999999996669;
    const double twofifthsnearest= 0.40000000000000002220;

    Ariadne::set_builtin_rounding_mode(downward);
    double onethirdrounddown=one/three;
    ARIADNE_TEST_EQUAL(onethirdrounddown, onethirddown);
    Ariadne::set_builtin_rounding_mode(upward);
    double onethirdroundup=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundup, onethirdup);
    Ariadne::set_builtin_rounding_mode(toward_zero);
    double onethirdroundchop=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundchop, onethirdchop);
    Ariadne::set_builtin_rounding_mode(to_nearest);
    double onethirdroundnearest=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundnearest, onethirdnearest);

    Ariadne::set_builtin_rounding_downward();
    double twofifthsrounddown=two_value/five;
    ARIADNE_TEST_EQUAL(twofifthsrounddown, twofifthsdown);
    Ariadne::set_builtin_rounding_upward();
    double twofifthsroundup=two_value/five;
    ARIADNE_TEST_EQUAL(twofifthsroundup, twofifthsup);
    Ariadne::set_builtin_rounding_toward_zero();
    double twofifthsroundchop=two_value/five;
    ARIADNE_TEST_EQUAL(twofifthsroundchop, twofifthschop);
    Ariadne::set_builtin_rounding_to_nearest();
    double twofifthsroundnearest=two_value/five;
    ARIADNE_TEST_EQUAL(twofifthsroundnearest, twofifthsnearest);
}

template<class PR> Void
TestFloatRounding<PR>::test_double_rounding()
{
}

template<class PR> Void
TestFloatRounding<PR>::test_arithmetic()
{
    cout << ARIADNE_PRETTY_FUNCTION << endl;

/*
    static bool full_precision_warning = false;
    const Float pi_near=Float::pi(near,precision);
    const Float four   =Float(4,precision);
    if( mul(up,pi_near,4.0_x) != mul(up,pi_near,four) ) {
        if(not full_precision_warning) {
            full_precision_warning=true;
            ARIADNE_TEST_WARN("Mixed Float/Builtin operations may not use full precision.");
        }
    }
*/

    const Float eps = Float::eps(precision);
    const Float one = Float(1.0_x,precision);
    const Float two_value = Float(2.0_x,precision);

    // Set next_up some variables
    Float f1(1.25_x,precision); Float f2(2.25_x,precision); Float f3(-3.25_x,precision);
    Float f4(precision); Float f5(precision);

    // Minimum (this should always remain exact)
    f4=min(f1,f2);
    cout << "min(" << f1 << "," << f2 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f1);
    f4=min(f1,f3);
    cout << "min(" << f1 << "," << f3 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f3);

    // Maximum (this should always remain exact)
    f4=max(f1,f2);
    cout << "max(" << f1 << "," << f2 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f2);
    f4=max(f1,f3);
    cout << "max(" << f1 << "," << f3 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f1);

    // Absolute value (this should always remain exact)
    f4=abs(f1);
    cout << "abs(" << f1 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==1.25_x);
    f5=abs(f3);
    cout << "abs(" << f3 << ") = " << f5 << endl;
    ARIADNE_TEST_ASSERT(f5==3.25_x);

    // Median (this should remain exact here)
    f3=med(near,f1,f2);
    cout << f1 << " <= med(" << f1 << "," << f2 << ")=" << f3 << " <= " << f2 << endl;
    ARIADNE_TEST_ASSERT(f1<=f3); ARIADNE_TEST_ASSERT(f3<=f2);
    ARIADNE_TEST_ASSERT(f3==1.75_x);

    // Negation (this should always remain exact)
    f3=neg(f1);
    cout << "neg(" << f1 << ") = " << f3 << endl;
    f3=-f1;
    cout << "- " << f1 << " = " << f3 << endl;
    ARIADNE_TEST_ASSERT(f3==-1.25_x);

    // Rounding
    f3=next(down,f1);
    f4=next(up,f1);
    cout << f3 << " < " << f1 << " < " << f4 << endl;
    ARIADNE_TEST_COMPARE(f3,<,f1); ARIADNE_TEST_COMPARE(f4,>,f1);

    // Addition (this should remain exact here)
    f3=add(down,f1,f2);
    f4=add(up,f1,f2);
    cout << f3 << " <= " << f1 << " + " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=3.5_x); ARIADNE_TEST_ASSERT(f4>=3.5_x);
    // Addition
    ARIADNE_TEST_COMPARE(add(up,one,hlf(eps)),>,1.0_x);
    ARIADNE_TEST_COMPARE(add(down,one,eps),>,1.0_x);
    if(add(down,one,hlf(eps)) != one) {
        ARIADNE_TEST_WARN("Results of floating-point operations stored to higher-precision in registers than memory.");
        ARIADNE_TEST_COMPARE(add(down,one,shft(eps,-11)),>,one);
    }
    ARIADNE_TEST_COMPARE(add(down,one,shft(eps,-12)),==,one);
    Float one_add_down_half_epsilon=add(down,one,hlf(eps));
    ARIADNE_TEST_COMPARE(one_add_down_half_epsilon,==,one);
    ARIADNE_TEST_EQUAL(one,1.0_x);
    ARIADNE_TEST_COMPARE(one,==,1.0_x);


    // Subtraction (this should remain exact here)
    f3=sub(down,f1,f2);
    f4=sub(up,f1,f2);
    cout << f3 << " <= " << f1 << " - " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=-1); ARIADNE_TEST_ASSERT(f4>=-1);

    // Multiplication (this should remain exact here)
    f3=mul(down,f1,f2);
    f4=mul(up,f1,f2);
    cout << f3 << " <= " << f1 << " * " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=2.8125_x); ARIADNE_TEST_ASSERT(f4>=2.8125_x);

    // Division (not exact; should catch errors here)
    f3=div(down,f1,f2);
    f4=div(up,f1,f2);
    f5=div(approx,f1,f2);
    cout << f3 << " <= " << f1 << " / " << f2 << " <= " << f4 << endl;
    Float five(5,precision);
    Float nine(9,precision);

    ExactDouble expected_five_ninths_up=0.55555555555555558023_x;
    ARIADNE_TEST_COMPARE(div(down,five,nine),<,expected_five_ninths_up);

    ARIADNE_TEST_COMPARE(div(down,Float(5.0_x,precision),Float(9.0_x,precision)),<,div(up,Float(5.0_x,precision),Float(9.0_x,precision)));
    ARIADNE_TEST_COMPARE(div(down,Float(5,precision),Float(9,precision)),<,div(up,Float(5,precision),Float(9,precision)));
    ARIADNE_TEST_COMPARE(div(down,five,nine),<,div(up,five,nine));
    Float five_ninths_down = div(down,five,nine);
    Float five_ninths_near = div(near,five,nine);
    Float five_ninths_up = div(up,five,nine);
    ARIADNE_TEST_COMPARE(five_ninths_down, < , five_ninths_up);
    ARIADNE_TEST_COMPARE(mul(down,five_ninths_down,nine), < , five);
    ARIADNE_TEST_COMPARE(mul(up,five_ninths_up,nine), > , five);

    Float::set_rounding_downward();
    ARIADNE_TEST_COMPARE(div(rounded,five,nine),==,five_ninths_down);
    Float::set_rounding_to_nearest();
    ARIADNE_TEST_COMPARE(div(rounded,five,nine),==,five_ninths_near);
    Float::set_rounding_upward();
    ARIADNE_TEST_COMPARE(div(rounded,five,nine),==,five_ninths_up);

    Float pi_down=Float::pi(down,precision);
    ARIADNE_TEST_EQUAL(div(down,pi_down,two_value),hlf(pi_down));

    // Power (not exact; should catch errors here)
    f3=pow(down,f1,3);
    f4=pow(up,f1,3);
    cout << f3 << " <= pow(" << f1 << ",3) <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=1.953125_x); ARIADNE_TEST_ASSERT(f4>=1.953125_x);

    f3=pow(down,f1,-2);
    f4=pow(up,f1,-2);
    cout << f3 << " <= pow(" << f1 << ",-2) <= " << f4 << endl;
    //ARIADNE_TEST_ASSERT(Rational(f3)<Rational(16,25));
    //ARIADNE_TEST_ASSERT(Rational(f4)>Rational(16,25));

    // Floor and ceiling
    f2=Float(-3.25_x,precision); f3=Float(-2,precision);

    ARIADNE_TEST_ASSERT(floor(f1)==1); ARIADNE_TEST_ASSERT(ceil(f1)==2);
    ARIADNE_TEST_ASSERT(floor(f2)==-4); ARIADNE_TEST_ASSERT(ceil(f2)==-3);
    ARIADNE_TEST_ASSERT(floor(f3)==-2); ARIADNE_TEST_ASSERT(ceil(f3)==-2);

    if constexpr (Same<PR,MultiplePrecision>) {
        // Independent exact rational/integer expectations for native
        // integral rounding, including ties away from zero.
        const FloatMP half(Dyadic(1,1u),precision);
        const FloatMP minus_half(Dyadic(-1,1u),precision);
        const FloatMP quarter(Dyadic(1,2u),precision);
        const FloatMP minus_quarter(Dyadic(-1,2u),precision);
        ARIADNE_TEST_EQUALS(round(half),FloatMP(1,precision));
        ARIADNE_TEST_EQUALS(round(minus_half),FloatMP(-1,precision));
        ARIADNE_TEST_EQUALS(round(quarter),FloatMP(0,precision));
        ARIADNE_TEST_ASSERT(is_zero(round(minus_quarter)));
        ARIADNE_TEST_ASSERT(round(minus_quarter).negative());
        ARIADNE_TEST_EQUALS(floor(quarter),FloatMP(0,precision));
        ARIADNE_TEST_EQUALS(ceil(quarter),FloatMP(1,precision));
        ARIADNE_TEST_EQUALS(floor(minus_quarter),FloatMP(-1,precision));
        ARIADNE_TEST_EQUALS(ceil(minus_quarter),FloatMP(0,precision));
        ARIADNE_TEST_ASSERT(not ceil(minus_quarter).negative());
        // Import of the most negative signed 64-bit integer must not
        // overflow by attempting to negate it in the signed domain.
        const std::int64_t min_i64=std::numeric_limits<std::int64_t>::min();
        const FloatMP min_i64_float(min_i64,MultiplePrecision(128_bits));
        ARIADNE_TEST_EQUALS(Rational(min_i64_float),Rational(Integer(min_i64)));
        // No huge shift allocation at the minimum signed exponent.
        const FloatMP tiny=FloatMP::min(precision);
        const FloatMP minus_tiny=-tiny;
        ARIADNE_TEST_EQUALS(floor(tiny),FloatMP(0,precision));
        ARIADNE_TEST_EQUALS(ceil(tiny),FloatMP(1,precision));
        ARIADNE_TEST_EQUALS(floor(minus_tiny),FloatMP(-1,precision));
        ARIADNE_TEST_EQUALS(ceil(minus_tiny),FloatMP(0,precision));
        ARIADNE_TEST_ASSERT(is_zero(round(minus_tiny)));
        ARIADNE_TEST_ASSERT(round(minus_tiny).negative());
    }

    // Conversion to integer types
    Int i3,i4;
    i3=integer_cast<Int>(floor(f1));
    i4=integer_cast<Int>(ceil(f1));
    cout << i3 << " < " << f1 << " < " << i4 << endl;
    ARIADNE_TEST_ASSERT(i3==1); ARIADNE_TEST_ASSERT(i4==2);
    i3=integer_cast<Int>(floor(f2));
    i4=integer_cast<Int>(ceil(f2));
    cout << i3 << " < " << f2 << " < " << i4 << endl;
    ARIADNE_TEST_ASSERT(i3==-4); ARIADNE_TEST_ASSERT(i4==-3);

    // Check interval conversions
    Float z(0,precision); Float o(1,precision); Float t(3,precision);

    Float odtd=div(down,o,t);
    Float odtu=div(up,o,t);
    cout << odtd << " <= 1/3 <= " << odtu << endl;
    ARIADNE_TEST_COMPARE(odtd,<,odtu);
    // Regression test to catch errors when Float result is not assigned to a variable
    cout << div(down,o,t) << " <= 1/3 <= " << div(up,o,t) << endl;
    ARIADNE_TEST_COMPARE(div(down,o,t),<,div(up,o,t));

    ARIADNE_TEST_COMPARE(med(near,Float(2,precision),Float(3,precision)),==,2.5_x);
    ARIADNE_TEST_COMPARE(rad(up,Float(2,precision),Float(3,precision)),>=,0.5_x);
    if constexpr (Same<Float,FloatDP>) {
        ARIADNE_TEST_COMPARE(rad(up,Float(2,precision),Float(3,precision)),<=,0.5000000000000002_pr);
    }

    // The following line should not compile
    // f5=f1+f2;

}


template<class PR> Void
TestFloatRounding<PR>::test_function()
{
    cout << ARIADNE_PRETTY_FUNCTION << endl;

    cout << setprecision(20);

    Float x(1,precision);
    Float ra=exp(near,x);
    Float rl=next(down,ra);
    Float ru=next(up,ra);
    ARIADNE_TEST_PRINT(rl);
    ARIADNE_TEST_PRINT(ru);
    ARIADNE_TEST_ASSERT(rl<ru);
    ARIADNE_TEST_ASSERT(2.71_dec<rl);
    ARIADNE_TEST_ASSERT(ru<2.72_dec);


    //The following don't work as rounded operators not exported.
    //test_inverse_pair("sin",&sin_down,&sin_up,&asin_down,&asin_up);
    //The following don't work as acos is decreasing
    //test_inverse_pair("cos",&cos_down,&cos_up,&acos_down,&acos_up);
    //test_inverse_pair("cos",&cos_up,&cos_down,&acos_down,&acos_up);
    //test_inverse_pair("tan",&tan_down,&tan_up,&atan_down,&atan_up);
    //test_inverse_pair("sinh",&sinh_down,&sinh_up,&asinh_down,&asinh_up);
    //test_inverse_pair("cosh",&cosh_down,&cosh_up,&acosh_down,&acosh_up);
    //test_inverse_pair("tanh",&tanh_down,&tanh_up,&atanh_down,&atanh_up);
}

template<class X> X twc(X x) { return shft(x,1); }
template<class X> X qtr(X x) { return shft(x,-2); }

template<class PR> Void
TestFloatRounding<PR>::test_cosine()
{
    //3.14159265358979323846264338327950288419716939937510
    const Float pi_down=Float::pi(down,precision);
    const Float pi_up  =Float::pi(up,precision);
    const Float three  =Float(3,precision);

    const Float third_pi_down=div(down,pi_down,three);
    const Float third_pi_up  =div(up  ,pi_up  ,three);

    const Float zero(0.0_x,precision);
    const Float one(1.0_x,precision);

    const Float three_pi_down=mul(down,pi_down,three);
    const Float three_pi_up  =mul(up  ,pi_up  ,three);

    ARIADNE_TEST_COMPARE(cos(up,hlf(pi_down)),>,0.0_x);
    ARIADNE_TEST_EQUAL(cos(up,zero),1.0_x);
    ARIADNE_TEST_COMPARE(cos(up,third_pi_down),>,0.5_x);
    ARIADNE_TEST_COMPARE(sqr(up,cos(up,qtr(pi_down))),>,0.5_x);
    ARIADNE_TEST_COMPARE(cos(up,hlf(pi_down)),>,0.0_x);
    ARIADNE_TEST_COMPARE(cos(up,hlf(pi_up)),<=,0.0_x);
    ARIADNE_TEST_COMPARE(cos(up,pi_down),>,-1.0_x);
    ARIADNE_TEST_COMPARE(cos(up,pi_up),>,-1.0_x);
    ARIADNE_TEST_COMPARE(cos(up,twc(pi_down)),==,1.0_x);
    ARIADNE_TEST_COMPARE(cos(up,twc(pi_up)),==,1.0_x);
    ARIADNE_TEST_COMPARE(cos(up,three_pi_down),>,-1.0_x);
    ARIADNE_TEST_COMPARE(cos(up,three_pi_up),>,-1.0_x);

    ARIADNE_TEST_EQUAL(cos(down,zero),1.0_x);
    ARIADNE_TEST_COMPARE(cos(down,third_pi_up),<,0.5_x);
    ARIADNE_TEST_COMPARE(sqr(down,cos(down,qtr(pi_up))),<,0.5_x);
    ARIADNE_TEST_COMPARE(cos(down,hlf(pi_down)),>=,0.0_x);
    ARIADNE_TEST_COMPARE(cos(down,hlf(pi_up)),<,0.0_x);
    ARIADNE_TEST_COMPARE(cos(down,pi_down),==,-1.0_x);
    ARIADNE_TEST_COMPARE(cos(down,pi_up),==,-1.0_x);
    ARIADNE_TEST_COMPARE(cos(down,twc(pi_down)),<,1.0_x);
    ARIADNE_TEST_COMPARE(cos(down,twc(pi_up)),<,1.0_x);
    ARIADNE_TEST_COMPARE(cos(down,three_pi_down),==,-1.0_x);
    ARIADNE_TEST_COMPARE(cos(down,three_pi_up),==,-1.0_x);

    Float sin_rnd_up_one=sin(up,one);
    Float sin_rnd_down_one=sin(down,one);
    Float sin_rnd_approx_one=sin(near,one);
    ARIADNE_TEST_COMPARE(sin_rnd_down_one,<=,sin_rnd_approx_one);
    ARIADNE_TEST_COMPARE(sin_rnd_approx_one,<=,sin_rnd_up_one);
    ARIADNE_TEST_COMPARE(sin_rnd_down_one,< ,sin_rnd_up_one);
}

template<> Void
TestFloatRounding<DoublePrecision>::test_arctan()
{
    DoublePrecision pr;

    static const FloatDP pi_down(3.1415926535897931_pr,pr);
    static const FloatDP pi_near(3.1415926535897931_pr,pr);
    static const FloatDP pi_up  (3.1415926535897936_pr,pr);

    static const FloatDP third_pi_down(1.0471975511965976_pr,pr);
    static const FloatDP third_pi_up  (1.0471975511965978_pr,pr);

    static const FloatDP quarter_pi_down(0.7853981633974483_pr,pr);
    //static const FloatDP quarter_pi_near(0.7853981633974483_pr,pr);
    static const FloatDP quarter_pi_up  (0.7853981633974484_pr,pr);

    static const FloatDP sixth_pi_down(0.52359877559829882_pr,pr);
    static const FloatDP sixth_pi_up  (0.52359877559829893_pr,pr);

    static const FloatDP atan_quarter_down(0.244978663126864143_pr,pr);
    //static const FloatDP atan_quarter_near(0.244978663126864154_pr,pr);
    static const FloatDP atan_quarter_up  (0.244978663126864171_pr,pr);

    static const FloatDP sqrt_three_down(1.732050807568877193_pr,pr);
    //static const FloatDP sqrt_three_near(1.732050807568877294_pr,pr);
    static const FloatDP sqrt_three_up  (1.732050807568877415_pr,pr);

    static const FloatDP sqrt_rec_three_down(0.57735026918962573_pr,pr);
    //static const FloatDP sqrt_three_near(0.57735026918962573_pr,pr);
    static const FloatDP sqrt_rec_three_up  (0.57735026918962584_pr,pr);

    static const FloatDP zero(0.0_x,pr);
    static const FloatDP one (1.0_x,pr);
    static const FloatDP quarter (0.25_x,pr);

    static const FloatDP eps=FloatDP::eps(pr);

    ARIADNE_TEST_EQUAL(atan(up,zero),zero);
    ARIADNE_TEST_COMPARE(atan(up,one),>=,quarter_pi_up);
    ARIADNE_TEST_COMPARE(atan(up,-one),>=,-quarter_pi_down);
    ARIADNE_TEST_COMPARE(atan(up,sqrt_three_up),>=,third_pi_up);
    ARIADNE_TEST_COMPARE(atan(up,-sqrt_three_down),>=,-third_pi_down);
    ARIADNE_TEST_COMPARE(atan(up,sqrt_rec_three_up),>=,sixth_pi_up);
    ARIADNE_TEST_COMPARE(atan(up,-sqrt_rec_three_down),>=,-sixth_pi_down);
    ARIADNE_TEST_COMPARE(atan(up,quarter),>=,atan_quarter_up);
    ARIADNE_TEST_COMPARE(atan(up,-quarter),>=,-atan_quarter_down);

    ARIADNE_TEST_EQUAL(atan(down,zero),zero);
    ARIADNE_TEST_COMPARE(atan(down,+one),<=,quarter_pi_down);
    ARIADNE_TEST_COMPARE(atan(down,-one),<=,-quarter_pi_up);
    ARIADNE_TEST_COMPARE(atan(down,+sqrt_three_down),<=,third_pi_down);
    ARIADNE_TEST_COMPARE(atan(down,-sqrt_three_up),<=,-third_pi_up);
    ARIADNE_TEST_COMPARE(atan(down,sqrt_rec_three_down),<=,sixth_pi_down);
    ARIADNE_TEST_COMPARE(atan(down,-sqrt_rec_three_up),<=,-sixth_pi_up);
    ARIADNE_TEST_COMPARE(atan(down,+quarter),<=,atan_quarter_down);
    ARIADNE_TEST_COMPARE(atan(down,-quarter),<=,-atan_quarter_up);

    ARIADNE_TEST_COMPARE(sub(down,atan(down,one),quarter_pi_up),<=,eps);
    ARIADNE_TEST_COMPARE(sub(down,atan(down,sqrt_three_up),third_pi_up),<=,eps);
    ARIADNE_TEST_COMPARE(sub(down,atan(down,sqrt_rec_three_down),sixth_pi_up),<=,eps);
    ARIADNE_TEST_COMPARE(sub(down,atan(down,quarter),atan_quarter_up),<=,hlf(eps));

}


template<class PR> Void
TestFloatRounding<PR>::test_arctan()
{
    //pi~=3.14159265358979323846264338327950288419716939937510
    const Float pi_down=Float::pi(down,precision);
    const Float pi_near=Float::pi(near,precision);
    const Float pi_up  =Float::pi(up,precision);

    Float three=Float(3,precision);
    const Float sqrt_three_down=sqrt(down,three);
    const Float sqrt_three_up  =sqrt(up,three);

    const Float sqrt_rec_three_down=sqrt(down,rec(down,three));
    const Float sqrt_rec_three_up  =sqrt(up,rec(up,three));

    const Float zero(0,precision);
    const Float one(1,precision);
    const Float four(4,precision);
    const Float six(6,precision);
    const Float eps=Float::eps(precision);

    ARIADNE_TEST_EQUAL(atan(up,zero),zero);
    ARIADNE_TEST_COMPARE(mul(up,atan(up,+one),four),>=,pi_up);
    ARIADNE_TEST_COMPARE(mul(up,atan(up,-one),four),>=,-pi_down);
    ARIADNE_TEST_COMPARE(mul(up,atan(up,+sqrt_three_up),three),>=,pi_up);
    ARIADNE_TEST_COMPARE(mul(up,atan(up,-sqrt_three_down),three),>=,-pi_down);
    ARIADNE_TEST_COMPARE(mul(up,atan(up,+sqrt_rec_three_up),six),>=,pi_up);
    ARIADNE_TEST_COMPARE(mul(up,atan(up,-sqrt_rec_three_down),six),>=,-pi_down);

    ARIADNE_TEST_EQUAL(atan(down,zero),zero);
    ARIADNE_TEST_COMPARE(mul(down,atan(down,+one),four),<=,pi_down);
    ARIADNE_TEST_COMPARE(mul(down,atan(down,-one),four),<=,-pi_up);
    ARIADNE_TEST_COMPARE(mul(down,atan(down,+sqrt_three_down),three),<=,pi_down);
    ARIADNE_TEST_COMPARE(mul(down,atan(down,-sqrt_three_up),three),<=,-pi_up);
    ARIADNE_TEST_COMPARE(mul(down,atan(down,+sqrt_rec_three_down),six),<=,pi_down);
    ARIADNE_TEST_COMPARE(mul(down,atan(down,-sqrt_rec_three_up),six),<=,-pi_up);

    ARIADNE_TEST_COMPARE(sub(down,mul(down,atan(down,one),four),pi_up),<=,shft(eps,2));
    ARIADNE_TEST_COMPARE(sub(down,mul(down,atan(down,sqrt_three_up),three),pi_up),<=,shft(eps,2));
    ARIADNE_TEST_COMPARE(sub(down,mul(down,atan(down,sqrt_rec_three_down),six),pi_up),<=,shft(eps,2));

    // Regression test
    ARIADNE_TEST_COMPARE(div(down,pi_down,4),<,div(up,pi_up,4));

}


template<class PR> Void TestFloatRounding<PR>::test_short_power_and_tanh() {
    if constexpr (Same<PR,MultiplePrecision>) {
        static const auto rows=load_floatmp_power_golden_rows(
            ARIADNE_TEST_DATA_DIR "/floatmp_pow_short.tsv",3680u);
        for (auto const& [key,expected]:rows) {
            const FloatMP value(key.input,near,MultiplePrecision(key.bits));
            const auto check=[&](auto mode) {
                ARIADNE_TEST_EQUALS(static_cast<Rational>(pow(mode,value,Nat(key.degree))),expected);
            };
            switch(key.mode) {
                case 'N': check(near); break;
                case 'D': check(down); break;
                case 'U': check(up); break;
                case 'Z': check(toward_zero); break;
                default: throw std::runtime_error("Invalid short-power oracle mode");
            }
        }
    }
    static const auto rows=load_floatmp_golden_rows(
        ARIADNE_TEST_DATA_DIR "/floatmp_tanh.tsv",1180u);
    for (auto const& [key,expected]:rows) {
        if constexpr (Same<PR,DoublePrecision>) { if(key.bits!=53u) continue; }
        const auto pr=[&]() {
            if constexpr (Same<PR,MultiplePrecision>) return PR(key.bits);
            else return PR{};
        }();
        const Float value(key.input,near,pr);
        const auto check=[&](auto mode) {
            ARIADNE_TEST_EQUALS(static_cast<Rational>(tanh(mode,value)),expected);
            const auto old=Float::get_rounding_mode();
            Float::set_rounding_mode(mode);
            ARIADNE_TEST_EQUALS(static_cast<Rational>(tanh(CurrentRoundingMode{},value)),expected);
            Float::set_rounding_mode(old);
        };
        switch(key.mode) {
            case 'N': check(near); break;
            case 'D': check(down); break;
            case 'U': check(up); break;
            case 'Z': check(toward_zero); break;
            default: throw std::runtime_error("Invalid tanh oracle mode");
        }
    }
    if constexpr (Same<PR,DoublePrecision>) {
        // Reuse independent frozen transcendental literals for the shared
        // binary64 adapter. Expected values remain offline exact certificates.
        // Recover nearest-even binary64 input exactly: Rational::get_d()
        // truncates toward zero, whereas the frozen oracles were generated
        // after nearest-even source conversion.
        const auto nearest_binary64=[](Rational const& source) {
            const double truncated=source.get_d();
            const Rational base{ExactDouble(truncated)};
            if (cmp(source,base)==Comparison::EQUAL) return truncated;
            const bool toward_negative=cmp(source,base)==Comparison::LESS;
            const double adjacent=std::nextafter(truncated,
                toward_negative ? -std::numeric_limits<double>::infinity()
                                :  std::numeric_limits<double>::infinity());
            const Rational midpoint=(base+Rational(ExactDouble(adjacent)))/Rational(2);
            const Comparison position=cmp(source,midpoint);
            if (position==Comparison::EQUAL) {
                return (std::bit_cast<std::uint64_t>(truncated)&1u)!=0u
                    ? adjacent : truncated;
            }
            const bool closer_to_adjacent=toward_negative
                ? position==Comparison::LESS : position==Comparison::GREATER;
            return closer_to_adjacent ? adjacent : truncated;
        };
        const auto verify=[&](char const* path,std::size_t count,
                              auto evaluate,bool check_nearest) {
            const auto golden=load_floatmp_golden_rows(path,count);
            const Rational window(Integer(1),pow(Integer(2),Nat(43u)));
            for (auto const& [key,expected]:golden) {
                if (key.bits!=53u || (key.mode=='N' && !check_nearest)) continue;
                const FloatDP input(ExactDouble(nearest_binary64(key.input)),dp);
                if (abs(abs(static_cast<Rational>(input))-Rational(1,2))>window) continue;
                const auto check=[&](auto mode) {
                    const auto old=FloatDP::get_rounding_mode();
                    ARIADNE_TEST_EQUALS(static_cast<Rational>(evaluate(mode,input)),expected);
                    ARIADNE_TEST_ASSERT(FloatDP::get_rounding_mode()==old);
                };
                switch(key.mode) {
                    case 'N': check(near); break;
                    case 'D': check(down); break;
                    case 'U': check(up); break;
                    case 'Z': check(toward_zero); break;
                    default: throw std::runtime_error("Invalid shared transcendental oracle mode");
                }
            }
        };
        verify(ARIADNE_FLOATMP_CENTERED_EXP_ORACLE_PATH,568u,
            [](auto mode,FloatDP value) { return exp(mode,value); },false);
        verify(ARIADNE_FLOATMP_LOG_ORACLE_PATH,452u,
            [](auto mode,FloatDP value) { return log(mode,value); },false);
        verify(ARIADNE_FLOATMP_SIN_ORACLE_PATH,1080u,
            [](auto mode,FloatDP value) { return sin(mode,value); },true);
        verify(ARIADNE_FLOATMP_COS_ORACLE_PATH,1080u,
            [](auto mode,FloatDP value) { return cos(mode,value); },true);
        verify(ARIADNE_FLOATMP_TAN_ORACLE_PATH,1208u,
            [](auto mode,FloatDP value) { return tan(mode,value); },false);
        verify(ARIADNE_FLOATMP_ATAN_ORACLE_PATH,1340u,
            [](auto mode,FloatDP value) { return atan(mode,value); },false);
        verify(ARIADNE_FLOATMP_ASIN_ORACLE_PATH,728u,
            [](auto mode,FloatDP value) { return asin(mode,value); },true);
        verify(ARIADNE_FLOATMP_ACOS_ORACLE_PATH,728u,
            [](auto mode,FloatDP value) { return acos(mode,value); },true);
    }
    const Float negative_zero(ExactDouble(-0.0),precision);
    ARIADNE_TEST_ASSERT(std::signbit(negative_zero.get_d()));
    ARIADNE_TEST_ASSERT(std::signbit(tanh(near,negative_zero).get_d()));
    for (auto mode : {Float::ROUND_TO_NEAREST,Float::ROUND_DOWNWARD,
                      Float::ROUND_UPWARD,Float::ROUND_TOWARD_ZERO}) {
        const Float converted(ExactDouble(-0.0),mode,precision);
        ARIADNE_TEST_ASSERT(std::signbit(converted.get_d()));
        ARIADNE_TEST_ASSERT(std::signbit(tanh(mode,converted).get_d()));
    }
    if constexpr (Same<PR,MultiplePrecision>) {
        const FloatDP source(ExactDouble(-0.0),dp);
        ARIADNE_TEST_ASSERT(Float(source,precision).negative());
        ARIADNE_TEST_ASSERT(Float(source,up,precision).negative());
    }
    ARIADNE_TEST_ASSERT(is_nan(tanh(near,Float::nan(precision))));
    ARIADNE_TEST_EQUALS(tanh(near,Float::inf(Sign::POSITIVE,precision)),Float(1,precision));
    ARIADNE_TEST_EQUALS(tanh(near,Float::inf(Sign::NEGATIVE,precision)),Float(-1,precision));
    if constexpr (Same<PR,DoublePrecision>) {
        const Float least(ExactDouble(std::numeric_limits<double>::denorm_min()),precision);
        ARIADNE_TEST_EQUALS(tanh(near,least),least);
        ARIADNE_TEST_EQUALS(tanh(up,least),least);
        ARIADNE_TEST_EQUALS(tanh(down,least),Float(0,precision));
        ARIADNE_TEST_EQUALS(tanh(toward_zero,least),Float(0,precision));
        ARIADNE_TEST_EQUALS(tanh(near,-least),-least);
        ARIADNE_TEST_EQUALS(tanh(down,-least),-least);
        ARIADNE_TEST_ASSERT(std::signbit(tanh(up,-least).get_d()));
        ARIADNE_TEST_EQUALS(tanh(up,-least),negative_zero);
    }
}
