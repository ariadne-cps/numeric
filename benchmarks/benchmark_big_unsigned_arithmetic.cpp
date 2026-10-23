/***************************************************************************
 *            benchmark_big_unsigned_arithmetic.cpp
 *
 *  Copyright  2026  Ariadne contributors
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

/*! \file benchmark_big_unsigned_arithmetic.cpp
 *  \brief Integrated GMP-differential benchmark for unsigned primitives.
 */

#include "numeric/big_unsigned.hpp"
#include "benchmark_command_line.hpp"
#include "logging/command_line.hpp"
#include "logging/logging.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <gmp.h>

namespace {

using namespace Ariadne;
using Clock=std::chrono::steady_clock;
using Limb=BigUnsigned::Limb;

struct GmpInteger {
    mpz_t value;
    GmpInteger() { mpz_init(value); }
    ~GmpInteger() { mpz_clear(value); }
    GmpInteger(GmpInteger const&)=delete;
    GmpInteger& operator=(GmpInteger const&)=delete;
};

template<class T> void keep_alive(T const& value)
{
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&value) : "memory");
#else
    static void const* volatile sink=nullptr;
    sink=&value;
#endif
}

template<class F> double median_ns(F&& work,unsigned iterations,unsigned samples)
{
    std::vector<double> durations;
    durations.reserve(samples);
    for (unsigned sample=0u;sample<=samples;++sample) {
        const auto begin=Clock::now();
        for (unsigned i=0u;i<iterations;++i) {
            work();
        }
        const auto end=Clock::now();
        if (sample!=0u) {
            durations.push_back(
                std::chrono::duration<double,std::nano>(end-begin).count()
                /static_cast<double>(iterations));
        }
    }
    std::sort(durations.begin(),durations.end());
    return durations[durations.size()/2u];
}

void import_gmp(mpz_ptr dst,BigUnsigned const& value)
{
    auto view=value.limbs();
    mpz_import(dst,view.size(),-1,sizeof(Limb),0,0,view.data());
}

void verify(BigUnsigned const& native,mpz_srcptr expected)
{
    GmpInteger actual;
    import_gmp(actual.value,native);
    if (mpz_cmp(actual.value,expected)!=0) {
        throw std::runtime_error("BigUnsigned differs from GMP");
    }
}

struct Selection {
    unsigned iterations=20000u;
    unsigned samples=31u;
    std::vector<unsigned> sizes{53u,128u,256u,1024u,2048u};
    bool summary_only=false;
};

List<ArgumentParser> parsers(Selection& config)
{
    List<ArgumentParser> result;
    result.append(make_valued_argument_parser(
        "i","iterations","Positive iterations per sample",
        [&config](String const& value) -> VoidFunction {
            const unsigned count=BenchmarkCommandLine::parse_unsigned<unsigned>(value);
            if (count==0u) throw std::invalid_argument("iterations must be positive");
            return [&config,count] { config.iterations=count; };
        }));
    result.append(make_valued_argument_parser(
        "n","samples","Odd number of samples, at least three",
        [&config](String const& value) -> VoidFunction {
            const unsigned count=BenchmarkCommandLine::parse_unsigned<unsigned>(value);
            if (count<3u || count%2u==0u) {
                throw std::invalid_argument("invalid sample count");
            }
            return [&config,count] { config.samples=count; };
        }));
    result.append(make_valued_argument_parser(
        "b","bits","Distinct comma-separated operand bit sizes",
        [&config](String const& value) -> VoidFunction {
            const std::string text(value);
            if (text.empty() || text.back()==',') {
                throw std::invalid_argument("invalid bit-size list");
            }
            std::vector<unsigned> parsed;
            std::istringstream input(text);
            std::string token;
            while (std::getline(input,token,',')) {
                const unsigned bits=BenchmarkCommandLine::parse_unsigned<unsigned>(
                    String(token));
                if (bits==0u || std::find(parsed.begin(),parsed.end(),bits)!=parsed.end()) {
                    throw std::invalid_argument("invalid or duplicate operand bit size");
                }
                parsed.push_back(bits);
            }
            if (parsed.empty()) throw std::invalid_argument("empty size list");
            return [&config,parsed] { config.sizes=parsed; };
        }));
    result.append(make_unvalued_argument_parser(
        "","summary-only","Print aggregate statistics only",
        [&config] { config.summary_only=true; }));
    return result;
}

struct Result {
    unsigned bits;
    std::string operation;
    double owning_native=0.0;
    double owning_gmp=0.0;
};

struct Aggregate {
    unsigned count=0u;
    double logs=0.0;
    double native_sum=0.0;
    double gmp_sum=0.0;
    double worst=0.0;
    std::string worst_operation;

    void add(Result const& r) {
        const double native=r.owning_native;
        const double gmp=r.owning_gmp;
        if (!(native>0.0 && gmp>0.0)) return;
        const double ratio=native/gmp;
        ++count;
        logs+=std::log(ratio);
        native_sum+=native;
        gmp_sum+=gmp;
        if (ratio>worst) {
            worst=ratio;
            worst_operation=r.operation+" @"+std::to_string(r.bits);
        }
    }

    void print(char const* group) const {
        if (count==0u) return;
        LOGGING_PRINTLN(group << "  ops=" << count
                        << " geomean=" << std::fixed << std::setprecision(3)
                        << std::exp(logs/count) << "x"
                        << " total_time=" << native_sum/gmp_sum << "x"
                        << " worst=" << worst << "x (" << worst_operation << ")")
    }
};

std::vector<Result> measure(unsigned bits,unsigned iterations,unsigned samples)
{
    // Equal-size unsigned inputs; subtraction is always nonnegative.
    std::mt19937_64 engine(0xB1BADDULL+bits);
    const std::size_t words=(static_cast<std::size_t>(bits)+63u)/64u;
    std::vector<Limb> aa(words),bb(words);
    for (auto& limb:aa) limb=engine();
    for (auto& limb:bb) limb=engine();
    const unsigned top=bits%64u;
    if (top!=0u) {
        const Limb mask=(Limb{1}<<top)-1u;
        aa.back()&=mask;
        bb.back()&=mask;
    }
    aa.back()|=Limb{1}<<((bits-1u)%64u);
    bb.back()|=Limb{1}<<((bits-1u)%64u);
    const BigUnsigned a(aa),b(bb);
    const bool swapped=compare(a,b)<0;
    BigUnsigned const& greater=swapped ? b : a;
    BigUnsigned const& smaller=swapped ? a : b;
    GmpInteger ga,gb,expected;
    import_gmp(ga.value,greater);
    import_gmp(gb.value,smaller);
    std::vector<Result> result;
    result.reserve(9u);

    mpz_add(expected.value,ga.value,gb.value);
    verify(add(greater,smaller),expected.value);
    Result add_row{bits,"add"};
    add_row.owning_native=median_ns([&] {
        auto value=add(greater,smaller);
        keep_alive(value);
    },iterations,samples);
    add_row.owning_gmp=median_ns([&] {
        GmpInteger value;
        mpz_add(value.value,ga.value,gb.value);
        keep_alive(value.value);
    },iterations,samples);
    result.push_back(add_row);

    mpz_sub(expected.value,ga.value,gb.value);
    verify(sub(greater,smaller),expected.value);
    Result sub_row{bits,"sub"};
    sub_row.owning_native=median_ns([&] {
        auto value=sub(greater,smaller);
        keep_alive(value);
    },iterations,samples);
    sub_row.owning_gmp=median_ns([&] {
        GmpInteger value;
        mpz_sub(value.value,ga.value,gb.value);
        keep_alive(value.value);
    },iterations,samples);
    result.push_back(sub_row);

    // Compare is already allocation-free on both implementations.
    const int native_comparison=compare(a,b);
    const int gmp_comparison=swapped ? -mpz_cmp(ga.value,gb.value)
                                     : mpz_cmp(ga.value,gb.value);
    if ((native_comparison<0)!=(gmp_comparison<0) ||
        (native_comparison>0)!=(gmp_comparison>0)) {
        throw std::runtime_error("comparison differs from GMP");
    }
    Result cmp_row{bits,"compare"};
    cmp_row.owning_native=median_ns([&] {
        const int value=compare(a,b);
        keep_alive(value);
    },iterations,samples);
    cmp_row.owning_gmp=median_ns([&] {
        const int value=swapped ? mpz_cmp(gb.value,ga.value)
                                : mpz_cmp(ga.value,gb.value);
        keep_alive(value);
    },iterations,samples);
    result.push_back(cmp_row);

    // Both shift implementations own their result.
    static_assert(std::numeric_limits<mp_bitcnt_t>::max()>=65u);
    for (unsigned shift_bits:{1u,65u}) {
        mpz_mul_2exp(expected.value,ga.value,
                      static_cast<mp_bitcnt_t>(shift_bits));
        verify(shift_left(greater,shift_bits),expected.value);
        Result left{bits,"shl_"+std::to_string(shift_bits)};
        left.owning_native=median_ns([&] {
            auto value=shift_left(greater,shift_bits);
            keep_alive(value);
        },iterations,samples);
        left.owning_gmp=median_ns([&] {
            GmpInteger value;
            mpz_mul_2exp(value.value,ga.value,
                         static_cast<mp_bitcnt_t>(shift_bits));
            keep_alive(value.value);
        },iterations,samples);
        result.push_back(left);

        mpz_fdiv_q_2exp(expected.value,ga.value,
                        static_cast<mp_bitcnt_t>(shift_bits));
        verify(shift_right(greater,shift_bits),expected.value);
        Result right{bits,"shr_"+std::to_string(shift_bits)};
        right.owning_native=median_ns([&] {
            auto value=shift_right(greater,shift_bits);
            keep_alive(value);
        },iterations,samples);
        right.owning_gmp=median_ns([&] {
            GmpInteger value;
            mpz_fdiv_q_2exp(value.value,ga.value,
                            static_cast<mp_bitcnt_t>(shift_bits));
            keep_alive(value.value);
        },iterations,samples);
        result.push_back(right);
    }

    // Only owning APIs exist for isqrt/sqrtrem. Both compared paths
    // allocate and destroy their result(s), with common immutable input
    // and identical iterations/samples. Verification stays OUTSIDE timing.
    mpz_sqrt(expected.value,ga.value);
    verify(isqrt(greater),expected.value);
    const unsigned sqrt_iterations=std::max(
        1u,iterations/std::max(1u,bits/64u));
    Result sqrt_row{bits,"isqrt"};
    sqrt_row.owning_native=median_ns([&] {
        auto value=isqrt(greater);
        keep_alive(value);
    },sqrt_iterations,samples);
    sqrt_row.owning_gmp=median_ns([&] {
        GmpInteger value;
        mpz_sqrt(value.value,ga.value);
        keep_alive(value.value);
    },sqrt_iterations,samples);
    verify(isqrt(greater),expected.value);
    result.push_back(sqrt_row);

    // Match two owned native results against two owned GMP mpz_t values.
    // The result pair is checked before/after the measured loops.
    GmpInteger expected_remainder;
    mpz_sqrtrem(expected.value,expected_remainder.value,ga.value);
    const auto check_sqrtrem=[&] {
        const auto [root,remainder]=sqrtrem(greater);
        verify(root,expected.value);
        verify(remainder,expected_remainder.value);
    };
    check_sqrtrem();
    Result sqrtrem_row{bits,"sqrtrem"};
    sqrtrem_row.owning_native=median_ns([&] {
        auto value=sqrtrem(greater);
        keep_alive(value);
    },sqrt_iterations,samples);
    sqrtrem_row.owning_gmp=median_ns([&] {
        GmpInteger root,remainder;
        mpz_sqrtrem(root.value,remainder.value,ga.value);
        keep_alive(root.value);
        keep_alive(remainder.value);
    },sqrt_iterations,samples);
    check_sqrtrem();
    result.push_back(sqrtrem_row);

    // Check all operations once again after measurement.
    mpz_add(expected.value,ga.value,gb.value);
    verify(add(greater,smaller),expected.value);
    mpz_sub(expected.value,ga.value,gb.value);
    verify(sub(greater,smaller),expected.value);
    for (unsigned shift_bits:{1u,65u}) {
        mpz_mul_2exp(expected.value,ga.value,
                     static_cast<mp_bitcnt_t>(shift_bits));
        verify(shift_left(greater,shift_bits),expected.value);
        mpz_fdiv_q_2exp(expected.value,ga.value,
                        static_cast<mp_bitcnt_t>(shift_bits));
        verify(shift_right(greater,shift_bits),expected.value);
    }
    return result;
}

void print_results(std::vector<Result> const& rows,bool summary_only)
{
    Aggregate owning;
    if (!summary_only) {
        LOGGING_PRINTLN("BigUnsigned arithmetic vs GMP (ns; owning results)")
        LOGGING_PRINTLN(std::string(66u,'='))
        LOGGING_PRINTLN(std::left << std::setw(10) << "Operation"
                        << std::setw(7) << "Bits"
                        << std::setw(16) << "Native (ns)"
                        << std::setw(16) << "GMP (ns)"
                        << std::setw(14) << "Native/GMP")
        LOGGING_PRINTLN(std::string(66u,'-'))
    }
    for (auto const& r:rows) {
        owning.add(r);
        if (!summary_only) {
            const auto show=[](double ns) {
                if (ns==0.0) return std::string("-");
                std::ostringstream os;
                os << std::fixed << std::setprecision(3) << ns;
                return os.str();
            };
            LOGGING_PRINTLN(std::left << std::setw(10) << r.operation
                            << std::setw(7) << r.bits
                            << std::setw(16) << show(r.owning_native)
                            << std::setw(16) << show(r.owning_gmp)
                            << std::setw(14) << show(r.owning_native/r.owning_gmp))
        }
    }
    LOGGING_PRINTLN(std::string(66u,'='))
    LOGGING_PRINTLN("BigUnsigned / GMP slowdown (x; >1 means slower)")
    owning.print("owning");
}

} // namespace

int main(int argc,const char* argv[])
{
    BenchmarkCommandLine::initialize_logging();
    try {
        Selection config;
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution(
            "numeric-big-unsigned-arithmetic",parsers(config)));
        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }
        LOGGING_PRINTLN("BigUnsigned arithmetic: " << config.sizes.size()
                        << " sizes | iterations=" << config.iterations
                        << " | samples=" << config.samples)
        std::vector<Result> data;
        for (unsigned bits:config.sizes) {
            LOGGING_PRINTLN("Measuring " << bits << " bits ...")
            auto rows=measure(bits,config.iterations,config.samples);
            data.insert(data.end(),rows.begin(),rows.end());
            LOGGING_PRINTLN("Completed " << bits << " bits")
        }
        print_results(data,config.summary_only);
        return 0;
    } catch (std::exception const& error) {
        LOGGING_PRINTLN("benchmark_big_unsigned_arithmetic: " << error.what())
        return 1;
    }
}
