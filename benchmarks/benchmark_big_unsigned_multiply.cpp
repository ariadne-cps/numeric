/***************************************************************************
 *            benchmark_big_unsigned_multiply.cpp
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
#include <sstream>
#include <string>
#include <random>
#include <stdexcept>
#include <vector>
#include <gmp.h>

namespace {

using namespace Ariadne;
using Clock = std::chrono::steady_clock;

struct GmpInteger {
    mpz_t value;
    GmpInteger() { mpz_init(value); }
    ~GmpInteger() { mpz_clear(value); }
    GmpInteger(GmpInteger const&) = delete;
    GmpInteger& operator=(GmpInteger const&) = delete;
};

template<class Function>
double median_ns(Function&& function,
                 unsigned iterations,unsigned samples)
{
    std::vector<double> times;
    times.reserve(samples);
    function();
    for (unsigned sample=0u; sample<samples; ++sample) {
        const auto start = Clock::now();
        for (unsigned i=0u; i<iterations; ++i) {
            function();
        }
        const auto finish = Clock::now();
        const auto elapsed =
            std::chrono::duration<double,std::nano>(finish-start).count();
        times.push_back(elapsed/static_cast<double>(iterations));
    }
    std::sort(times.begin(),times.end());
    return times[times.size()/2u];
}

void import_gmp(mpz_ptr target,BigUnsigned const& x)
{
    const auto digits = x.limbs();
    mpz_import(target,digits.size(),-1,sizeof(BigUnsigned::Limb),
               0,0,digits.data());
}

void verify_gmp(BigUnsigned const& actual,mpz_srcptr expected)
{
    GmpInteger converted;
    import_gmp(converted.value,actual);
    if (mpz_cmp(converted.value,expected)!=0) {
        throw std::runtime_error("BigUnsigned multiply differs from GMP");
    }
}

struct MultiplicationDatum {
    unsigned bits;
    double prealloc_portable_ns;
    double prealloc_selected_ns;
    double gmp_ns;
    double new_selected_ns;
    double prealloc_comba_ns=0.0; // 0 means the experiment was disabled
    double prealloc_fused_ns=0.0; // 0 means the experiment was disabled
    double prealloc_schoolbook_ns=0.0; // 0 means the experiment was disabled
    double new_mpz_ns=0.0; // 0 means the experiment was disabled
    double prealloc_karatsuba_ns=0.0; // 0 means the experiment was disabled
};

struct MultiplicationSelection {
    unsigned iterations=20000u;
    unsigned samples=31u;
    std::vector<unsigned> sizes{53u,128u,256u,1024u,2048u};
    bool summary_only=false;
    bool details=false;
    bool include_comba=false;
    bool include_fused=false;
    bool include_schoolbook=false;
    bool include_new_mpz=false;
    bool include_karatsuba=false;
    double max_slowdown=std::numeric_limits<double>::infinity();
};

List<ArgumentParser> multiplication_parsers(MultiplicationSelection& selection)
{
    List<ArgumentParser> parsers;
    parsers.append(make_valued_argument_parser(
        "i","iterations","Number of operations in each measured sample (positive)",
        [&selection](String const& value) -> VoidFunction {
            const unsigned count=BenchmarkCommandLine::parse_unsigned<unsigned>(value);
            if (count==0u) { throw std::invalid_argument("iterations must be positive"); }
            return [&selection,count] { selection.iterations=count; };
        }));
    parsers.append(make_valued_argument_parser(
        "n","samples","Odd number of samples, at least three",
        [&selection](String const& value) -> VoidFunction {
            const unsigned count=BenchmarkCommandLine::parse_unsigned<unsigned>(value);
            if (count<3u || count%2u==0u) {
                throw std::invalid_argument("invalid sample count");
            }
            return [&selection,count] { selection.samples=count; };
        }));
    parsers.append(make_valued_argument_parser(
        "b","bits","Comma-separated operand sizes in bits (positive)",
        [&selection](String const& value) -> VoidFunction {
            const std::string text(value);
            if (text.empty() || text.back()==',') {
                throw std::invalid_argument("invalid size list");
            }
            std::vector<unsigned> parsed;
            std::istringstream input(text);
            std::string token;
            while (std::getline(input,token,',')) {
                const unsigned bits=
                    BenchmarkCommandLine::parse_unsigned<unsigned>(String(token));
                if (bits==0u ||
                    std::find(parsed.begin(),parsed.end(),bits)!=parsed.end()) {
                    throw std::invalid_argument("invalid or duplicate operand size");
                }
                parsed.push_back(bits);
            }
            if (parsed.empty()) { throw std::invalid_argument("empty size list"); }
            return [&selection,parsed] { selection.sizes=parsed; };
        }));
    parsers.append(make_unvalued_argument_parser(
        "","summary-only","Print only aggregate slowdown statistics",
        [&selection] { selection.summary_only=true; }));
    parsers.append(make_unvalued_argument_parser(
        "","details","Print per-size native/GMP slowdown ratios",
        [&selection] { selection.details=true; }));
    parsers.append(make_unvalued_argument_parser(
        "","include-comba","Include experimental x86-64 prealloc_comba kernel",
        [&selection] { selection.include_comba=true; }));
    parsers.append(make_unvalued_argument_parser(
        "","include-fused","Include explicitly forced x86-64 prealloc_fused kernel",
        [&selection] { selection.include_fused=true; }));
    parsers.append(make_unvalued_argument_parser(
        "","include-schoolbook","Include forced B3 prealloc_schoolbook kernel",
        [&selection] { selection.include_schoolbook=true; }));
    parsers.append(make_unvalued_argument_parser(
        "","include-new-mpz",
        "Include GMP multiplication with new result allocation and release",
        [&selection] { selection.include_new_mpz=true; }));
    parsers.append(make_unvalued_argument_parser(
        "","include-karatsuba",
        "Include experimental Karatsuba with reusable preallocated scratch",
        [&selection] { selection.include_karatsuba=true; }));
    parsers.append(make_valued_argument_parser(
        "","max-slowdown","Fail if worst prealloc_selected/GMP slowdown exceeds this ratio",
        [&selection](String const& value) -> VoidFunction {
            const std::string text(value);
            std::size_t used=0u;
            const double ratio=std::stod(text,&used);
            if (used!=text.size() || !std::isfinite(ratio) || ratio<=0.0) {
                throw std::invalid_argument("invalid positive slowdown ratio");
            }
            return [&selection,ratio] { selection.max_slowdown=ratio; };
        }));
    return parsers;
}

MultiplicationDatum measure(unsigned bits,unsigned iterations,unsigned samples,
                            bool include_comba,bool include_fused,
                            bool include_schoolbook,bool include_new_mpz,
                            bool include_karatsuba)
{
    std::mt19937_64 engine(0xB3B3B3ULL+bits);
    const std::size_t n=(bits+63u)/64u;
    std::vector<BigUnsigned::Limb> aa(n),bb(n);
    for (auto& v:aa) v=engine();
    for (auto& v:bb) v=engine();
    const auto top=static_cast<unsigned>(bits%64u);
    if (top!=0u) {
        aa.back() &= (std::uint64_t{1}<<top)-1u;
        bb.back() &= (std::uint64_t{1}<<top)-1u;
    }
    aa.back() |= std::uint64_t{1}<<((bits-1u)%64u);
    bb.back() |= std::uint64_t{1}<<((bits-1u)%64u);
    const BigUnsigned a(aa),b(bb);
    BigUnsigned::Scratch scratch;
    const auto result_size = a.size()+b.size();
    auto output=scratch.acquire(result_size);
    BigUnsigned::Scratch karatsuba_workspace;
    const auto karatsuba_work=include_karatsuba
        ? karatsuba_workspace.acquire(
            BigUnsigned::karatsuba_scratch_limbs(a.size(),b.size()))
        : BigUnsigned::MutableLimbView{};
    GmpInteger ga,gb,gr;
    import_gmp(ga.value,a);
    import_gmp(gb.value,b);
    mpz_mul(gr.value,ga.value,gb.value);
    BigUnsigned::mul_into_portable(a.limbs(),b.limbs(),output);
    verify_gmp(BigUnsigned(output),gr.value);
    BigUnsigned::mul_into(a.limbs(),b.limbs(),output);
    verify_gmp(BigUnsigned(output),gr.value);
    verify_gmp(mul(a,b),gr.value);
    if (include_comba) {
        BigUnsigned::mul_into_comba(a.limbs(),b.limbs(),output);
        verify_gmp(BigUnsigned(output),gr.value);
    }
    if (include_fused) {
        BigUnsigned::mul_into_fused(a.limbs(),b.limbs(),output);
        verify_gmp(BigUnsigned(output),gr.value);
    }
    if (include_schoolbook) {
        BigUnsigned::mul_into_schoolbook(a.limbs(),b.limbs(),output);
        verify_gmp(BigUnsigned(output),gr.value);
    }
    if (include_karatsuba) {
        BigUnsigned::mul_into_karatsuba(
            a.limbs(),b.limbs(),output,karatsuba_work);
        verify_gmp(BigUnsigned(output),gr.value);
    }
    if (include_new_mpz) {
        GmpInteger fresh;
        mpz_mul(fresh.value,ga.value,gb.value);
        if (mpz_cmp(fresh.value,gr.value)!=0) {
            throw std::runtime_error("GMP new-result multiplication mismatch");
        }
    }

    const auto portable=median_ns([&] {
        BigUnsigned::mul_into_portable(a.limbs(),b.limbs(),output);
    },iterations,samples);
    const auto native=median_ns([&] {
        BigUnsigned::mul_into(a.limbs(),b.limbs(),output);
    },iterations,samples);
    const auto raw=median_ns([&] {
        mpz_mul(gr.value,ga.value,gb.value);
    },iterations,samples);
    const auto owned=median_ns([&] {
        const BigUnsigned product=mul(a,b);
        if (product.is_zero()) std::abort();
    },iterations,samples);
    // Comba is measured only after the original four B3 measurements,
    // preserving their order and the meaning of their existing baselines.
    double comba=0.0;
    if (include_comba) {
        comba=median_ns([&] {
            BigUnsigned::mul_into_comba(a.limbs(),b.limbs(),output);
        },iterations,samples);
    }
    // The optional fused kernel comes after the frozen B3 measurements
    // (and after optional Comba); all default timings are unchanged.
    double fused=0.0;
    if (include_fused) {
        fused=median_ns([&] {
            BigUnsigned::mul_into_fused(a.limbs(),b.limbs(),output);
        },iterations,samples);
    }
    double schoolbook=0.0;
    if (include_schoolbook) {
        schoolbook=median_ns([&] {
            BigUnsigned::mul_into_schoolbook(a.limbs(),b.limbs(),output);
        },iterations,samples);
    }
    // New GMP destination lifecycle is intentionally opt-in. It includes
    // mpz_init, mpz_mul and mpz_clear in every timed iteration, unlike the
    // original "mpz" case that reuses its destination. The nonzero guard
    // mirrors "new_selected", and no correctness oracle runs in the loop.
    double new_mpz=0.0;
    if (include_new_mpz) {
        new_mpz=median_ns([&] {
            GmpInteger result;
            mpz_mul(result.value,ga.value,gb.value);
            if (mpz_sgn(result.value)==0) std::abort();
        },iterations,samples);
    }
    // Karatsuba is measured last and only when explicitly requested.
    // All buffers (including recursive scratch) were allocated above.
    double karatsuba=0.0;
    if (include_karatsuba) {
        karatsuba=median_ns([&] {
            BigUnsigned::mul_into_karatsuba(
                a.limbs(),b.limbs(),output,karatsuba_work);
        },iterations,samples);
    }
    // Verify after timing as well, with no exact checks inside timed code.
    verify_gmp(BigUnsigned(output),gr.value);
    verify_gmp(mul(a,b),gr.value);
    if (include_new_mpz) {
        GmpInteger fresh;
        mpz_mul(fresh.value,ga.value,gb.value);
        if (mpz_cmp(fresh.value,gr.value)!=0) {
            throw std::runtime_error("GMP new-result multiplication mismatch");
        }
    }
    return {bits,portable,native,raw,owned,comba,fused,schoolbook,
            new_mpz,karatsuba};
}

void print_measurements(std::vector<MultiplicationDatum> const& data,
                        unsigned iterations,unsigned samples,bool summary_only)
{
    if (summary_only) { return; }
    LOGGING_PRINTLN("BigUnsigned multiplication benchmark | iterations="
                    << iterations << " | samples=" << samples)
    LOGGING_PRINTLN(std::string(79u,'='))
    LOGGING_PRINTLN(std::left << std::setw(18) << "Backend"
                    << "  " << std::setw(24) << "Operation"
                    << "  " << std::setw(6) << "Bits"
                    << "  " << std::right << std::setw(15) << "Median (ns)")
    LOGGING_PRINTLN(std::string(79u,'-'))
    for (auto const& row : data) {
        const auto print_row=[&](char const* backend,double nanoseconds) {
            LOGGING_PRINTLN(std::left << std::setw(18) << backend
                            << "  " << std::setw(24) << "multiply"
                            << "  " << std::setw(6) << row.bits
                            << "  " << std::right << std::setw(15)
                            << std::fixed << std::setprecision(3)
                            << nanoseconds)
        };
        print_row("prealloc_portable",row.prealloc_portable_ns);
        print_row("prealloc_selected",row.prealloc_selected_ns);
        print_row("mpz",row.gmp_ns);
        print_row("new_selected",row.new_selected_ns);
        if (row.prealloc_comba_ns>0.0) {
            print_row("prealloc_comba",row.prealloc_comba_ns);
        }
        if (row.prealloc_fused_ns>0.0) {
            print_row("prealloc_fused",row.prealloc_fused_ns);
        }
        if (row.prealloc_schoolbook_ns>0.0) {
            print_row("prealloc_schoolbook",row.prealloc_schoolbook_ns);
        }
        if (row.new_mpz_ns>0.0) {
            print_row("new_mpz",row.new_mpz_ns);
        }
        if (row.prealloc_karatsuba_ns>0.0) {
            print_row("prealloc_karatsuba",row.prealloc_karatsuba_ns);
        }
        LOGGING_PRINTLN(std::string(79u,'-'))
    }
    LOGGING_PRINTLN(std::string(79u,'='))
}

struct MultiplicationAggregate {
    std::size_t count=0u;
    double log_sum=0.0;
    double reference_sum=0.0;
    double candidate_sum=0.0;
    double worst_ratio=0.0;
    unsigned worst_bits=0u;

    void add(unsigned bits,double reference,double candidate) {
        if (!std::isfinite(reference) || !std::isfinite(candidate)
            || reference<=0.0 || candidate<=0.0) {
            throw std::runtime_error("nonpositive or nonfinite multiply timing");
        }
        const double ratio=candidate/reference;
        ++count;
        log_sum+=std::log(ratio);
        reference_sum+=reference;
        candidate_sum+=candidate;
        if (ratio>worst_ratio) {
            worst_ratio=ratio;
            worst_bits=bits;
        }
    }

    void print(char const* name) const {
        if (count==0u || reference_sum<=0.0) {
            throw std::runtime_error("empty multiplication benchmark aggregate");
        }
        LOGGING_PRINTLN(std::left << std::setw(18) << name
                        << "  ops=" << std::setw(3) << count
                        << "  geomean=" << std::fixed << std::setprecision(3)
                        << std::exp(log_sum/static_cast<double>(count)) << "x"
                        << "  total_time=" << candidate_sum/reference_sum << "x"
                        << "  worst=" << worst_ratio << "x ("
                        << worst_bits << " bits)")
    }
};

void print_slowdowns(std::vector<MultiplicationDatum> const& data,
                     bool details,double max_slowdown)
{
    MultiplicationAggregate portable,selected,owned,comba,fused,schoolbook;
    MultiplicationAggregate new_mpz,owned_vs_new_mpz,karatsuba;
    if (details) {
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("Detailed comparisons (candidate / GMP)")
        LOGGING_PRINTLN(std::string(79u,'-'))
        LOGGING_PRINTLN(std::left << std::setw(18) << "Backend"
                        << std::setw(8) << "Bits"
                        << std::setw(16) << "GMP (ns)"
                        << std::setw(16) << "Native (ns)" << "Slowdown")
    }
    for (auto const& row : data) {
        portable.add(row.bits,row.gmp_ns,row.prealloc_portable_ns);
        selected.add(row.bits,row.gmp_ns,row.prealloc_selected_ns);
        owned.add(row.bits,row.gmp_ns,row.new_selected_ns);
        if (row.prealloc_comba_ns>0.0) {
            comba.add(row.bits,row.gmp_ns,row.prealloc_comba_ns);
        }
        if (row.prealloc_fused_ns>0.0) {
            fused.add(row.bits,row.gmp_ns,row.prealloc_fused_ns);
        }
        if (row.prealloc_schoolbook_ns>0.0) {
            schoolbook.add(row.bits,row.gmp_ns,row.prealloc_schoolbook_ns);
        }
        if (row.new_mpz_ns>0.0) {
            new_mpz.add(row.bits,row.gmp_ns,row.new_mpz_ns);
            owned_vs_new_mpz.add(
                row.bits,row.new_mpz_ns,row.new_selected_ns);
        }
        if (row.prealloc_karatsuba_ns>0.0) {
            karatsuba.add(
                row.bits,row.gmp_ns,row.prealloc_karatsuba_ns);
        }
        if (details) {
            const auto print_detail=[&](char const* backend,double candidate) {
                LOGGING_PRINTLN(std::left << std::setw(18) << backend
                                << std::setw(8) << row.bits
                                << std::setw(16) << std::fixed
                                << std::setprecision(3) << row.gmp_ns
                                << std::setw(16) << candidate
                                << candidate/row.gmp_ns << "x")
            };
            print_detail("prealloc_portable",row.prealloc_portable_ns);
            print_detail("prealloc_selected",row.prealloc_selected_ns);
            print_detail("new_selected",row.new_selected_ns);
            if (row.prealloc_comba_ns>0.0) {
                print_detail("prealloc_comba",row.prealloc_comba_ns);
            }
            if (row.prealloc_fused_ns>0.0) {
                print_detail("prealloc_fused",row.prealloc_fused_ns);
            }
            if (row.prealloc_schoolbook_ns>0.0) {
                print_detail("prealloc_schoolbook",row.prealloc_schoolbook_ns);
            }
            if (row.new_mpz_ns>0.0) {
                print_detail("new_mpz",row.new_mpz_ns);
            }
            if (row.prealloc_karatsuba_ns>0.0) {
                print_detail("prealloc_karatsuba",row.prealloc_karatsuba_ns);
            }
        }
    }
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN(std::string(79u,'='))
    LOGGING_PRINTLN("BigUnsigned multiplication / GMP slowdown (x; >1 means slower)")
    LOGGING_PRINTLN(std::string(79u,'-'))
    portable.print("prealloc_portable");
    selected.print("prealloc_selected");
    owned.print("new_selected");
    if (comba.count!=0u) {
        comba.print("prealloc_comba");
    }
    if (fused.count!=0u) {
        fused.print("prealloc_fused");
    }
    if (schoolbook.count!=0u) {
        schoolbook.print("prealloc_schoolbook");
    }
    if (new_mpz.count!=0u) {
        new_mpz.print("new_mpz");
    }
    if (karatsuba.count!=0u) {
        karatsuba.print("prealloc_karatsuba");
    }
    LOGGING_PRINTLN(std::string(79u,'='))
    if (owned_vs_new_mpz.count!=0u) {
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("new_selected / new_mpz slowdown (x; both create a result)")
        LOGGING_PRINTLN(std::string(79u,'-'))
        owned_vs_new_mpz.print("new_selected");
        LOGGING_PRINTLN(std::string(79u,'='))
    }
    LOGGING_PRINTLN("geomean weights each operand size equally; "
                    "total_time weights by measured cost.")
    if (std::isfinite(max_slowdown) &&
        selected.worst_ratio>max_slowdown) {
        throw std::runtime_error(
            "prealloc_selected slowdown exceeds --max-slowdown: "+
            std::to_string(selected.worst_ratio)+"x > "+
            std::to_string(max_slowdown)+"x");
    }
}

} // namespace

int main(int argc,const char* argv[])
{
    BenchmarkCommandLine::initialize_logging();
    try {
        MultiplicationSelection selection;
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution(
            "numeric-big-unsigned-multiplication",
            multiplication_parsers(selection)));

        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--max-slowdown","--scheduler","-s","--theme","-t",
             "--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }
        if ((selection.include_comba || selection.include_fused ||
                selection.include_schoolbook) &&
                !BigUnsigned::has_native_multiply()) {
            throw std::invalid_argument(
                "--include-comba, --include-fused and --include-schoolbook require x86-64");
        }
        LOGGING_PRINTLN("BigUnsigned multiplication - running "
                        << selection.sizes.size() << " operand-size groups"
                        << " | native available="
                        << BigUnsigned::has_native_multiply()
                        << " | selected kernel="
                        << (BigUnsigned::selected_multiply_uses_fused()
                            ? "fused" : "B3 schoolbook/portable"))
        LOGGING_PRINTLN("Iterations: " << selection.iterations
                        << "  |  Samples: " << selection.samples)

        std::vector<MultiplicationDatum> data;
        data.reserve(selection.sizes.size());
        for (unsigned bits : selection.sizes) {
            LOGGING_PRINTLN("Measuring " << bits << " bits ...")
            data.push_back(measure(bits,selection.iterations,selection.samples,
                                   selection.include_comba,
                                   selection.include_fused,
                                   selection.include_schoolbook,
                                   selection.include_new_mpz,
                                   selection.include_karatsuba));
            LOGGING_PRINTLN("Completed " << bits << " bits")
        }
        print_measurements(data,selection.iterations,selection.samples,
                           selection.summary_only);
        print_slowdowns(data,selection.details,selection.max_slowdown);
        return 0;
    } catch (std::exception const& error) {
        LOGGING_PRINTLN("benchmark_big_unsigned_multiply: " << error.what())
        return 1;
    }
}
