/***************************************************************************
 *            benchmark_integer_decimal_profile.cpp
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

/*! \file benchmark_integer_decimal_profile.cpp
 *  \brief Independent diagnostic of the native decimal conversion stages.
 */

// B44: do not modify production arithmetic or the official 84-case benchmark.
// Separate component replays are non-additive proxies, not CPU attribution.
#include "numeric/big_unsigned.hpp"
#include "benchmark_exact_arithmetic_report.hpp"

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace Ariadne;
using namespace AriadneBenchmark;

constexpr std::size_t cutoff_limbs=16u;

struct DecimalPower {
    BigUnsigned value;
    std::size_t width;
};

struct Split {
    BigUnsigned dividend;
    BigUnsigned divisor;
};

struct Leaf {
    BigUnsigned value;
    std::size_t minimum_width;
    std::string decimal;
};

struct ProfileCase {
    BigUnsigned value;
    std::vector<Split> splits;
    std::vector<Leaf> leaves;
    std::size_t digit_capacity=0u;
    std::string expected;
};

// Same B22 power-selection arithmetic as the production formatter.
std::vector<DecimalPower> make_powers(BigUnsigned const& value)
{
    std::vector<DecimalPower> powers;
    if (value.is_zero() || value.size()<=cutoff_limbs) return powers;
    powers.push_back({BigUnsigned(std::uint64_t{1000000000u}),9u});
    const std::size_t midpoint=value.bit_length()/2u;
    while (true) {
        auto const& last=powers.back();
        const std::size_t bits=last.value.bit_length();
        if (bits>(midpoint+1u)/2u) break;
        if (last.width>std::numeric_limits<std::size_t>::max()/2u) {
            throw std::length_error("decimal profile power width overflow");
        }
        BigUnsigned squared=square(last.value);
        if (squared.bit_length()>midpoint) break;
        const std::size_t width=last.width*2u;
        powers.push_back({std::move(squared),width});
    }
    return powers;
}

// Obtain exactly the production recursion's division operands and leaves
// outside every timed loop; GMP verifies the reconstructed output.
void collect(BigUnsigned const& value,std::vector<DecimalPower> const& powers,
             std::size_t minimum_width,ProfileCase& entry)
{
    if (value.size()>cutoff_limbs) {
        const std::size_t midpoint=value.bit_length()/2u;
        std::size_t count=powers.size();
        while (count!=0u &&
                powers[count-1u].value.bit_length()>midpoint) {
            --count;
        }
        if (count!=0u) {
            auto const& divisor=powers[count-1u];
            entry.splits.push_back({value,divisor.value});
            auto [high,low]=divrem(value,divisor.value);
            const std::size_t high_width=minimum_width>divisor.width
                ? minimum_width-divisor.width : 0u;
            collect(high,powers,high_width,entry);
            collect(low,powers,divisor.width,entry);
            return;
        }
    }
    // The production to_decimal directly calls to_decimal_blocks at
    // <=16 limbs, including the architecture-selected wide ARM64 path.
    entry.leaves.push_back({value,minimum_width,to_decimal(value)});
}

std::string replay_append(ProfileCase const& entry)
{
    std::string output;
    if (entry.digit_capacity<=output.max_size()) {
        output.reserve(entry.digit_capacity);
    }
    for (auto const& leaf:entry.leaves) {
        if (leaf.decimal.size()<leaf.minimum_width) {
            output.append(leaf.minimum_width-leaf.decimal.size(),'0');
        }
        output+=leaf.decimal;
    }
    return output;
}

std::vector<ProfileCase> prepare(unsigned bits)
{
    auto const specs=integer_specifications(bits); // Exact official pool.
    std::vector<ProfileCase> pool;
    pool.reserve(specs.size());
    for (auto const& spec:specs) {
        RawInteger oracle(spec.a);
        mpz_abs(oracle.value,oracle.value);
        ProfileCase entry;
        entry.expected=decimal(oracle.value); // Independent GMP oracle.
        entry.value=from_decimal(entry.expected);
        entry.digit_capacity=entry.value.bit_length()/3u+1u;
        auto const powers=make_powers(entry.value);
        collect(entry.value,powers,0u,entry);
        if (to_decimal(entry.value)!=entry.expected ||
                (entry.value.size()>cutoff_limbs &&
                 replay_append(entry)!=entry.expected)) {
            throw std::runtime_error(
                "decimal profile differs from GMP at "+
                std::to_string(bits)+" bits");
        }
        pool.push_back(std::move(entry));
    }
    return pool;
}

template<class F>
double measure(F&& operation,unsigned iterations,unsigned samples)
{
    return exact_arithmetic_median_ns(
        std::forward<F>(operation),iterations,samples);
}

void report(char const* name,double nanoseconds,double full)
{
    std::ostringstream formatted;
    formatted << std::left << std::setw(28) << name
              << std::right << std::setw(12) << std::fixed
              << std::setprecision(3) << nanoseconds << " ns"
              << "  (" << std::setprecision(1) << 100.0*nanoseconds/full
              << "% of full; non-additive)";
    LOGGING_PRINTLN(formatted.str())
}

void benchmark(unsigned bits,unsigned iterations,unsigned samples)
{
    auto const pool=prepare(bits);
    const auto select=[&](unsigned i) -> ProfileCase const& {
        return pool[i%pool.size()];
    };
    const double full=measure([&](unsigned i) {
        return to_decimal(select(i).value);
    },iterations,samples);

    LOGGING_PRINTLN(" ") 
    LOGGING_PRINTLN(bits << " bits | 32 GMP-verified official inputs"
                    << " | iterations=" << iterations
                    << " | samples=" << samples)
    report("production to_decimal",full,full);
    if (pool.front().value.size()<=cutoff_limbs) {
        LOGGING_PRINTLN("Power table / recursive divrem / append: not used "
                        "(<=16 limbs).")
        return;
    }
    const double powers=measure([&](unsigned i) {
        return make_powers(select(i).value);
    },iterations,samples);
    const double divisions=measure([&](unsigned i) {
        std::size_t checksum=0u;
        for (auto const& split:select(i).splits) {
            auto [quotient,remainder]=divrem(split.dividend,split.divisor);
            checksum+=quotient.size()+remainder.size();
        }
        return checksum;
    },iterations,samples);
    const double leaves=measure([&](unsigned i) {
        std::size_t checksum=0u;
        for (auto const& leaf:select(i).leaves) {
            auto value=to_decimal(leaf.value);
            checksum+=value.size();
        }
        return checksum;
    },iterations,samples);
    const double appending=measure([&](unsigned i) {
        return replay_append(select(i));
    },iterations,samples);
    report("power-table construction",powers,full);
    report("recursive divrem replay",divisions,full);
    report("production leaf conversion",leaves,full);
    report("append / zero-padding replay",appending,full);
    LOGGING_PRINTLN("Diagnostic replays are independently timed. "
                    "Do not add their percentages.")
}

} // namespace

int main(int argc,const char* argv[])
{
    try {
        BenchmarkCommandLine::initialize_logging();
        ExactArithmeticSelection selection;
        selection.sizes={1024u,2048u,4096u};
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution(
            "numeric-integer-decimal-profile",
            exact_arithmetic_parsers(selection,false)));
        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }
        for (unsigned bits:selection.sizes) {
            benchmark(bits,selection.iterations,selection.samples);
        }
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("Numerical verification PASS (GMP oracle "
                       "and full production conversion).")
        return 0;
    } catch (std::exception const& error) {
        std::cerr << "Integer decimal profile FAILED: "
                  << error.what() << '\n';
        return 1;
    }
}
