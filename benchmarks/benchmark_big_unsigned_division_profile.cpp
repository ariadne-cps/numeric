/***************************************************************************
 *            benchmark_big_unsigned_division_profile.cpp
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

/*! \file benchmark_big_unsigned_division_profile.cpp
 *  \brief B56 diagnostic: compare exact division by divisor shape.
 *
 *  This benchmark is deliberately separate from the immutable official
 *  84-case Integer/GMP benchmark. All operand pools are deterministic;
 *  independent GMP verification occurs entirely outside timed loops.
 *  Native and GMP owning lifetimes are matched, as are the reusable
 *  joint-result buffers in the preallocated comparison. Preallocated
 *  joint-result timings are NOT additive to the owning timings.
 */

#include "benchmark_exact_arithmetic_report.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>

namespace {

using namespace Ariadne;
using namespace AriadneBenchmark;
using Limb=BigUnsigned::Limb;
using Clock=std::chrono::steady_clock;

constexpr std::size_t pool_size=32u;

struct Gmp {
    mpz_t value;
    Gmp() { mpz_init(value); }
    ~Gmp() { mpz_clear(value); }
    Gmp(Gmp const&)=delete;
    Gmp& operator=(Gmp const&)=delete;
    // B56.1: make fixture pools movable without copying GMP internals.
    // Initialize our own object, then swap GMP-owned storage. A moved
    // source remains a valid, independently destructible mpz_t.
    Gmp(Gmp&& other) noexcept {
        mpz_init(value);
        mpz_swap(value,other.value);
    }
    Gmp& operator=(Gmp&& other) noexcept {
        if (this!=&other) {
            mpz_swap(value,other.value);
        }
        return *this;
    }
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

template<class Function>
double median_ns(Function&& function,unsigned iterations,unsigned samples)
{
    std::vector<double> measurements;
    measurements.reserve(samples);
    for (unsigned sample=0u;sample<=samples;++sample) {
        const auto begin=Clock::now();
        for (unsigned i=0u;i<iterations;++i) {
            function(i);
        }
        const auto end=Clock::now();
        if (sample!=0u) {
            measurements.push_back(
                std::chrono::duration<double,std::nano>(end-begin).count()/
                static_cast<double>(iterations));
        }
    }
    std::sort(measurements.begin(),measurements.end());
    return measurements[samples/2u];
}

void import_gmp(mpz_ptr out,BigUnsigned const& number)
{
    const auto view=number.limbs();
    mpz_import(out,view.size(),-1,sizeof(Limb),0,0,view.data());
}

void verify(std::vector<Limb> const& words,mpz_srcptr expected,
            std::string const& context)
{
    Gmp got;
    mpz_import(got.value,words.size(),-1,sizeof(Limb),0,0,words.data());
    if (mpz_cmp(got.value,expected)!=0) {
        throw std::runtime_error("preallocated division differs from GMP: "+context);
    }
}

void verify(BigUnsigned const& native,mpz_srcptr expected,
            std::string const& context)
{
    Gmp got;
    import_gmp(got.value,native);
    if (mpz_cmp(got.value,expected)!=0) {
        throw std::runtime_error("division differs from GMP: "+context);
    }
}

enum class Shape {
    limb_small,
    limb_wide,
    short_two,
    short_eight,
    quarter_dense,
    half_dense,
    half_sparse,
    long_dense,
    near_less,
    near_greater,
    shifted_half,
    power_of_two
};

constexpr std::array<Shape,12u> shapes{
    Shape::limb_small,Shape::limb_wide,Shape::short_two,
    Shape::short_eight,Shape::quarter_dense,Shape::half_dense,
    Shape::half_sparse,Shape::long_dense,Shape::near_less,
    Shape::near_greater,Shape::shifted_half,Shape::power_of_two
};

char const* name(Shape shape)
{
    switch (shape) {
    case Shape::limb_small: return "limb_small";
    case Shape::limb_wide: return "limb_wide";
    case Shape::short_two: return "short_2";
    case Shape::short_eight: return "short_8";
    case Shape::quarter_dense: return "quarter_dense";
    case Shape::half_dense: return "half_dense";
    case Shape::half_sparse: return "half_sparse";
    case Shape::long_dense: return "long_dense";
    case Shape::near_less: return "near_less_q0";
    case Shape::near_greater: return "near_greater_q1";
    case Shape::shifted_half: return "shifted_half";
    case Shape::power_of_two: return "power_of_two";
    }
    throw std::logic_error("unknown division shape");
}

struct Fixture {
    BigUnsigned dividend,divisor;
    Gmp a,b,prealloc_q,prealloc_r;
    std::vector<Limb> quotient_words,remainder_words,scratch_words;
};

std::array<Fixture,pool_size> make_pool(unsigned bits,Shape shape)
{
    // Whole-limb input sizes make the case groups directly comparable.
    if (bits<512u || bits%64u!=0u) {
        throw std::invalid_argument(
            "division-profile --bits requires multiples of 64 >=512");
    }
    const std::size_t n=bits/64u;
    std::array<Fixture,pool_size> result;
    std::mt19937_64 rng(0xB560D1A6ULL+
                        static_cast<std::uint64_t>(bits)*1009u+
                        static_cast<std::uint64_t>(shape)*65537u);
    for (std::size_t k=0u;k<result.size();++k) {
        std::size_t m=n/2u;
        switch (shape) {
        case Shape::limb_small:
        case Shape::limb_wide: m=1u; break;
        case Shape::short_two: m=2u; break;
        case Shape::short_eight: m=8u; break;
        case Shape::quarter_dense: m=n/4u; break;
        case Shape::long_dense: m=n-1u; break;
        case Shape::near_less:
        case Shape::near_greater: m=n; break;
        default: break;
        }
        std::vector<Limb> a(n),b(m);
        for (auto& digit:a) digit=rng();
        for (auto& digit:b) digit=rng();
        a.back()|=Limb{1}<<63u;
        b.back()|=Limb{1}<<63u;
        switch (shape) {
        case Shape::limb_small:
            b[0u]=static_cast<Limb>(3u+2u*(k%31u));
            break;
        case Shape::limb_wide:
            b[0u]|=Limb{1}<<63u;
            break;
        case Shape::half_sparse:
            std::fill(b.begin(),b.end(),Limb{0});
            b.front()=static_cast<Limb>(3u+k*2u);
            b[m/2u]=Limb{11};
            b.back()=Limb{1}<<63u;
            break;
        case Shape::near_less:
        case Shape::near_greater:
            b=a;
            if (shape==Shape::near_less) {
                a[0u]&=~Limb{1}; // even => adding one cannot overflow
                b[0u]=a[0u]+1u; // a<b; exact quotient 0
            } else {
                a[0u]|=Limb{1}; // odd => subtraction cannot underflow
                b[0u]=a[0u]-1u; // a>b; exact quotient 1
            }
            break;
        case Shape::shifted_half:
            std::fill_n(b.begin(),m/4u,Limb{0});
            break;
        case Shape::power_of_two:
            std::fill(b.begin(),b.end(),Limb{0});
            b.back()=Limb{1}<<63u;
            break;
        default: break;
        }
        auto& fixture=result[k];
        fixture.dividend=BigUnsigned(a);
        fixture.divisor=BigUnsigned(b);
        import_gmp(fixture.a.value,fixture.dividend);
        import_gmp(fixture.b.value,fixture.divisor);
        if (mpz_sgn(fixture.b.value)==0) {
            throw std::logic_error("division benchmark zero divisor");
        }
        fixture.quotient_words.resize(fixture.dividend.size());
        fixture.remainder_words.resize(
            std::min(fixture.dividend.size(),fixture.divisor.size()));
        fixture.scratch_words.resize(BigUnsigned::divrem_scratch_limbs(
            fixture.dividend.size(),fixture.divisor.size()));
        // Avoid any GMP result growth in the reusable-output benchmark.
        mpz_realloc2(fixture.prealloc_q.value,2u*bits+64u);
        mpz_realloc2(fixture.prealloc_r.value,2u*bits+64u);
    }
    return result;
}

void verify_pool(std::array<Fixture,pool_size>& pool,
                 unsigned bits,Shape shape)
{
    for (std::size_t k=0u;k<pool.size();++k) {
        auto& test=pool[k];
        Gmp q,r;
        mpz_tdiv_qr(q.value,r.value,test.a.value,test.b.value);
        const std::string context=std::to_string(bits)+" bits "+
            name(shape)+" case "+std::to_string(k);
        verify(quot(test.dividend,test.divisor),q.value,context+" quot");
        verify(rem(test.dividend,test.divisor),r.value,context+" rem");
        const auto both=divrem(test.dividend,test.divisor);
        verify(both.first,q.value,context+" divrem quot");
        verify(both.second,r.value,context+" divrem rem");
        BigUnsigned::divrem_into(
            test.dividend.limbs(),test.divisor.limbs(),
            test.quotient_words,test.remainder_words,test.scratch_words);
        verify(test.quotient_words,q.value,context+" prealloc quot");
        verify(test.remainder_words,r.value,context+" prealloc rem");
        if (shape==Shape::near_less && mpz_sgn(q.value)!=0) {
            throw std::runtime_error("near_less did not produce quotient 0");
        }
        if (shape==Shape::near_greater && mpz_cmp_ui(q.value,1ul)!=0) {
            throw std::runtime_error("near_greater did not produce quotient 1");
        }
    }
}

struct Result {
    unsigned bits=0u;
    std::string shape,operation;
    double native_ns=0.0,gmp_ns=0.0;
};

std::vector<Result> measure(unsigned bits,Shape shape,
                            unsigned iterations,unsigned samples)
{
    auto pool=make_pool(bits,shape);
    verify_pool(pool,bits,shape);
    auto select=[&](unsigned i)->Fixture& {
        return pool[static_cast<std::size_t>(i)%pool.size()];
    };
    std::vector<Result> result;
    result.reserve(4u);
    const auto record=[&](char const* operation,auto native,auto gmp) {
        Result entry;
        entry.bits=bits;
        entry.shape=name(shape);
        entry.operation=operation;
        entry.native_ns=median_ns(native,iterations,samples);
        entry.gmp_ns=median_ns(gmp,iterations,samples);
        if (entry.native_ns<=0.0 || entry.gmp_ns<=0.0) {
            throw std::runtime_error("nonpositive division timing");
        }
        result.push_back(std::move(entry));
    };

    record("quot",
        [&](unsigned i) {
            auto answer=quot(select(i).dividend,select(i).divisor);
            keep_alive(answer);
        },
        [&](unsigned i) {
            auto& item=select(i);
            Gmp answer;
            mpz_tdiv_q(answer.value,item.a.value,item.b.value);
            keep_alive(answer.value);
        });
    record("rem",
        [&](unsigned i) {
            auto answer=rem(select(i).dividend,select(i).divisor);
            keep_alive(answer);
        },
        [&](unsigned i) {
            auto& item=select(i);
            Gmp answer;
            mpz_tdiv_r(answer.value,item.a.value,item.b.value);
            keep_alive(answer.value);
        });
    record("divrem",
        [&](unsigned i) {
            auto answer=divrem(select(i).dividend,select(i).divisor);
            keep_alive(answer);
        },
        [&](unsigned i) {
            auto& item=select(i);
            Gmp q,r;
            mpz_tdiv_qr(q.value,r.value,item.a.value,item.b.value);
            keep_alive(q.value);
            keep_alive(r.value);
        });
    record("prealloc_divrem",
        [&](unsigned i) {
            auto& item=select(i);
            BigUnsigned::divrem_into(
                item.dividend.limbs(),item.divisor.limbs(),
                item.quotient_words,item.remainder_words,item.scratch_words);
            keep_alive(item.quotient_words);
            keep_alive(item.remainder_words);
        },
        [&](unsigned i) {
            auto& item=select(i);
            mpz_tdiv_qr(item.prealloc_q.value,item.prealloc_r.value,
                        item.a.value,item.b.value);
            keep_alive(item.prealloc_q.value);
            keep_alive(item.prealloc_r.value);
        });
    verify_pool(pool,bits,shape); // Guard correctness after all timed calls.
    return result;
}

void print_results(std::vector<Result> const& rows,bool summary_only)
{
    if (!summary_only) {
        LOGGING_PRINTLN("BigUnsigned division profile (ns; matched modes)")
        LOGGING_PRINTLN(std::string(92u,'='))
        LOGGING_PRINTLN(std::left << std::setw(18) << "Shape"
                        << std::setw(9) << "Bits"
                        << std::setw(18) << "Operation"
                        << std::right << std::setw(15) << "Integer"
                        << std::setw(15) << "GMP"
                        << std::setw(13) << "Ratio")
        LOGGING_PRINTLN(std::string(92u,'-'))
        for (auto const& row:rows) {
            LOGGING_PRINTLN(std::left << std::setw(18) << row.shape
                            << std::setw(9) << row.bits
                            << std::setw(18) << row.operation
                            << std::right << std::fixed << std::setprecision(3)
                            << std::setw(15) << row.native_ns
                            << std::setw(15) << row.gmp_ns
                            << std::setw(12) << row.native_ns/row.gmp_ns << "x")
        }
        LOGGING_PRINTLN(std::string(92u,'='))
    }
    LOGGING_PRINTLN("Worst native/GMP slowdown by operation (diagnostic categories)")
    for (char const* operation:{"quot","rem","divrem","prealloc_divrem"}) {
        Result const* worst=nullptr;
        for (auto const& row:rows) {
            if (row.operation!=operation) continue;
            if (worst==nullptr ||
                    row.native_ns/row.gmp_ns >
                    worst->native_ns/worst->gmp_ns) {
                worst=&row;
            }
        }
        if (worst!=nullptr) {
            LOGGING_PRINTLN(std::left << std::setw(18) << operation
                            << std::setw(18) << worst->shape
                            << worst->bits << " bits  "
                            << std::fixed << std::setprecision(3)
                            << worst->native_ns/worst->gmp_ns << "x")
        }
    }
    LOGGING_PRINTLN("GMP verification PASS; all fixture generation, "
                    "oracles and comparisons are outside timed loops.")
    LOGGING_PRINTLN("prealloc_divrem uses reusable native spans and "
                    "reusable GMP result buffers; it is not additive "
                    "to owning-operation timings.")
}

} // namespace

int main(int argc,const char* argv[])
{
    BenchmarkCommandLine::initialize_logging();
    try {
        ExactArithmeticSelection selection;
        selection.sizes={1024u,2048u,4096u};
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution(
            "numeric-big-unsigned-division-profile",
            exact_arithmetic_parsers(selection,false)));
        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }
        LOGGING_PRINTLN("B56 division profile: " << selection.sizes.size()
                        << " sizes x " << shapes.size() << " divisor shapes"
                        << " x 32 exact GMP-verified input pairs"
                        << " | iterations=" << selection.iterations
                        << " | samples=" << selection.samples)
        std::vector<Result> rows;
        for (unsigned bits:selection.sizes) {
            if (bits<512u || bits%64u!=0u) {
                throw std::invalid_argument(
                    "division profile requires --bits multiples of 64 >=512");
            }
            LOGGING_PRINTLN("Measuring " << bits << " bits ...")
            for (auto shape:shapes) {
                auto measured=measure(
                    bits,shape,selection.iterations,selection.samples);
                rows.insert(rows.end(),measured.begin(),measured.end());
            }
            LOGGING_PRINTLN("Completed " << bits << " bits")
        }
        print_results(rows,selection.summary_only);
        return 0;
    } catch (std::exception const& error) {
        LOGGING_PRINTLN("benchmark_big_unsigned_division_profile: "
                        << error.what())
        return 1;
    }
}
