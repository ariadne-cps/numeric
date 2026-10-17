/***************************************************************************
 *            benchmark_compare_multiple_precision_backends.cpp
 *
 *  Copyright  2026  Luca Geretti
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

// Same monomorphic test code is instantiated for each implementation: no virtual dispatch.
#include "numeric/floatmp.hpp"
#include "numeric/floatap.hpp"
#include "numeric/rational.hpp"
#include "logging/logging.hpp"
#include "logging/command_line.hpp"
#include "benchmark_command_line.hpp"
#include <algorithm>
#include <bit>
#include <type_traits>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <limits>
#include <map>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Ariadne;

#include "benchmark_backend_workload.hpp"
using namespace AriadneBenchmark;

// All reports use the Ariadne logger, following the integration benchmarks.
void print_results(std::vector<BenchmarkDatum> const& data,
                   unsigned iterations, unsigned samples,
                   bool summary_only) {
    if (summary_only) return;

    constexpr int backend_width=10;
    constexpr int operation_width=24;
    constexpr int timing_width=15;
    constexpr int table_width=backend_width+operation_width+timing_width+12;
    LOGGING_PRINTLN("FloatMP / FloatAP benchmark | iterations=" << iterations
                    << " | samples=" << samples)
    LOGGING_PRINTLN(std::string(static_cast<std::size_t>(table_width), '='))
    LOGGING_PRINTLN(std::left << std::setw(backend_width) << "Backend"
                    << "  " << std::setw(operation_width) << "Operation"
                    << "  " << std::setw(6) << "Bits"
                    << "  " << std::right << std::setw(timing_width) << "Median (ns)")
    LOGGING_PRINTLN(std::string(static_cast<std::size_t>(table_width), '-'))
    unsigned previous_bits=0;
    for (auto const& datum : data) {
        if (previous_bits != 0u && previous_bits != datum.bits) {
            LOGGING_PRINTLN(std::string(static_cast<std::size_t>(table_width), '-'))
        }
        previous_bits=datum.bits;
        LOGGING_PRINTLN(std::left << std::setw(backend_width) << datum.backend
                        << "  " << std::setw(operation_width) << datum.operation
                        << "  " << std::setw(6) << datum.bits
                        << "  " << std::right << std::setw(timing_width)
                        << std::fixed << std::setprecision(3) << datum.nanoseconds)
    }
    LOGGING_PRINTLN(std::string(static_cast<std::size_t>(table_width), '='))
}

void print_slowdown_summary(std::vector<BenchmarkDatum> const& data,
bool details, double max_slowdown) {
    struct Aggregate {
        std::size_t count = 0;
        double logarithmic_sum = 0.0;
        double mpfr_sum = 0.0;
        double floatap_sum = 0.0;
        double worst_ratio = 0.0;
        std::string worst_operation;
        unsigned worst_bits = 0;

        void add(std::string const& name, unsigned bits, double reference, double candidate) {
            if (!(std::isfinite(reference) && std::isfinite(candidate))
            || reference <= 0.0 || candidate <= 0.0) {
                throw std::runtime_error("nonpositive/nonfinite benchmark timing");
            }
            const double ratio = candidate / reference;
            ++count;
            logarithmic_sum += std::log(ratio);
            mpfr_sum += reference;
            floatap_sum += candidate;
            if (ratio > worst_ratio) {
                worst_ratio = ratio;
                worst_operation = name;
                worst_bits = bits;
            }
        }
        void print(std::string const& name) const {
            if (count == 0 || mpfr_sum <= 0.0) {
                throw std::runtime_error("cannot summarize empty benchmark results");
            }
            LOGGING_PRINTLN(std::left << std::setw(12) << name
                << "  ops=" << std::setw(3) << count
                << "  geomean=" << std::fixed << std::setprecision(3)
                << std::exp(logarithmic_sum / static_cast<double>(count)) << "x"
                << "  total_time=" << floatap_sum / mpfr_sum << "x"
                << "  worst=" << worst_ratio << "x (" << worst_operation
                << ", " << worst_bits << " bits)")
        }
    };

    using Key = std::pair < unsigned, std::string >;
    std::map< Key, double > reference, candidate;
    for (auto const& datum : data) {
        auto& collection = datum.backend == "mpfr" ? reference : candidate;
        if (!collection.emplace(Key { datum.bits, datum.operation }, datum.nanoseconds).second) {
            throw std::runtime_error("duplicate benchmark measurement");
        }
    }
    if (reference.size() != candidate.size() || reference.empty()) {
        throw std::runtime_error("FloatAP and FloatMP benchmark rows are not paired");
    }
    std::map< unsigned, Aggregate > per_precision;
    Aggregate overall;
    if (details) {
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("Detailed comparisons (FloatAP / MPFR)")
        LOGGING_PRINTLN(std::string(79u,'-'))
        LOGGING_PRINTLN(std::left << std::setw(24) << "Operation"
                        << std::setw(8) << "Bits"
                        << std::setw(16) << "MPFR (ns)"
                        << std::setw(16) << "FloatAP (ns)" << "Slowdown")
    }
    for (auto const& [key, reference_ns] : reference) {
        auto it = candidate.find(key);
        if (it == candidate.end()) {
            throw std::runtime_error("FloatAP benchmark operation missing: " + key.second);
        }
        auto& precision_summary = per_precision[key.first];
        precision_summary.add(key.second, key.first, reference_ns, it->second);
        overall.add(key.second, key.first, reference_ns, it->second);
        if (details) {
            LOGGING_PRINTLN(std::left << std::setw(24) << key.second
                            << std::setw(8) << key.first
                            << std::setw(16) << std::fixed << std::setprecision(3)
                            << reference_ns << std::setw(16) << it->second
                            << std::setprecision(3) << it->second / reference_ns << "x")
        }
    }
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN(std::string(78u, '='))
    LOGGING_PRINTLN("FloatAP / FloatMP slowdown (x; >1 means slower)")
    LOGGING_PRINTLN(std::string(78u, '-'))
    for (auto const& [bits, summary] : per_precision) {
        summary.print(std::to_string(bits) + " bits");
    }
    overall.print("ALL");
    LOGGING_PRINTLN(std::string(78u, '='))
    LOGGING_PRINTLN("geomean weighs each (operation, precision) equally; "
                    "total_time weights by measured cost.")
    if (std::isfinite(max_slowdown) && overall.worst_ratio > max_slowdown) {
        throw std::runtime_error(
        "FloatAP slowdown exceeds --max-slowdown: "
        + std::to_string(overall.worst_ratio) + "x > "
        + std::to_string(max_slowdown) + "x");
    }
}


// This verification is deliberately performed only AFTER ALL stopwatch
// measurements. Comparing output values inside median_ns would contaminate
// the operation timings, particularly for transcendental functions.
// The independent expected-literal tests remain the correctness oracle:
// MPFR agreement alone is additional evidence, not a mathematical proof.
namespace {
struct CheckedValue {
    enum class Kind { exact_rational, comparison, binary64 };
    Kind kind = Kind::exact_rational;
    Rational exact;
    Comparison relation = Comparison::EQUAL;
    std::uint64_t binary64 = 0u;
};

template<class Value>
CheckedValue normalize_checked(Value const& value, unsigned bits,
                               char const* name, unsigned index) {
    CheckedValue result;
    if constexpr (std::is_same_v<Value,Comparison>) {
        result.kind = CheckedValue::Kind::comparison;
        result.relation = value;
    } else if constexpr (std::is_same_v<Value,double>) {
        result.kind = CheckedValue::Kind::binary64;
        result.binary64 = std::bit_cast<std::uint64_t>(value);
    } else {
        if constexpr (std::is_same_v<Value,FloatMP> ||
                      std::is_same_v<Value,FloatAP>) {
            if (static_cast<unsigned>(value.precision().bits()) != bits) {
                std::ostringstream error;
                error << "unexpected output precision: " << bits
                      << " bits, " << name << ", input " << index;
                throw std::runtime_error(error.str());
            }
        }
        // Conversion to Rational is *exact* (no double, no tolerance).
        result.exact = static_cast<Rational>(value);
    }
    return result;
}

void assert_same(CheckedValue const& mpfr, CheckedValue const& ap,
                 unsigned bits, std::string const& operation, unsigned index) {
    bool same = mpfr.kind == ap.kind;
    if (same) {
        switch (mpfr.kind) {
            case CheckedValue::Kind::exact_rational:
                same = mpfr.exact == ap.exact;
                break;
            case CheckedValue::Kind::comparison:
                same = mpfr.relation == ap.relation;
                break;
            case CheckedValue::Kind::binary64:
                same = mpfr.binary64 == ap.binary64;
                break;
        }
    }
    if (!same) {
        std::ostringstream error;
        error << "FloatMP/FloatAP value mismatch at " << bits
              << " bits, operation " << operation << ", input " << index;
        if (mpfr.kind == CheckedValue::Kind::exact_rational &&
            ap.kind == CheckedValue::Kind::exact_rational) {
            error << ": MPFR=" << mpfr.exact << "; FloatAP=" << ap.exact;
        } else if (mpfr.kind == CheckedValue::Kind::binary64 &&
                   ap.kind == CheckedValue::Kind::binary64) {
            error << ": MPFR double bits=0x" << std::hex << mpfr.binary64
                  << "; FloatAP double bits=0x" << ap.binary64;
        } else if (mpfr.kind == CheckedValue::Kind::comparison &&
                   ap.kind == CheckedValue::Kind::comparison) {
            error << ": comparison outcomes differ";
        } else {
            error << ": result kinds differ";
        }
        throw std::runtime_error(error.str());
    }
}
} // namespace

void verify_benchmark_values(std::vector<unsigned> const& precisions) {
    constexpr unsigned count = BackendPool<FloatMP>::pool_size;
    std::size_t verified_inputs = 0;
    std::size_t verified_operations = 0;
    for (unsigned bits : precisions) {
        // These use EXACTLY the same rational pool constructor as timing.
        const BackendPool<FloatMP> reference(bits);
        const BackendPool<FloatAP> candidate(bits);
        for (unsigned i = 0; i < count; ++i) {
            if (reference.exact_a[i] != candidate.exact_a[i] ||
                reference.exact_b[i] != candidate.exact_b[i]) {
                throw std::runtime_error("benchmark exact input pools diverge");
            }
            assert_same(normalize_checked(reference.a[i],bits,"input_a",i),
                        normalize_checked(candidate.a[i],bits,"input_a",i),
                        bits,"input_a",i);
            assert_same(normalize_checked(reference.b[i],bits,"input_b",i),
                        normalize_checked(candidate.b[i],bits,"input_b",i),
                        bits,"input_b",i);
            verified_inputs += 2u;
        }

        // Compare each of the 27 operation definitions at ALL 32 indices;
        // never test only the first benchmark input or merely compare times.
        std::map<std::string,std::vector<CheckedValue>> expected;
        for_each_backend_operation(reference,[&](char const* name, auto&& operation) {
            auto [it, inserted] = expected.try_emplace(name);
            if (!inserted) {
                throw std::runtime_error(std::string("duplicate operation: ")+name);
            }
            auto& outputs = it->second;
            outputs.reserve(count);
            for (unsigned i = 0; i < count; ++i) {
                outputs.push_back(normalize_checked(operation(i),bits,name,i));
            }
        });
        std::size_t operations_this_precision = 0;
        for_each_backend_operation(candidate,[&](char const* name, auto&& operation) {
            const auto found = expected.find(name);
            if (found == expected.end() || found->second.size() != count) {
                throw std::runtime_error(std::string("missing MPFR operation: ")+name);
            }
            for (unsigned i = 0; i < count; ++i) {
                assert_same(found->second[i],
                            normalize_checked(operation(i),bits,name,i),
                            bits,name,i);
                ++verified_operations;
            }
            ++operations_this_precision;
        });
        if (operations_this_precision != expected.size() ||
            operations_this_precision != 27u) {
            throw std::runtime_error("benchmark operation count changed");
        }
    }
    LOGGING_PRINTLN("Numerical verification PASS: "
                    << verified_inputs << " paired inputs and "
                    << verified_operations
                    << " paired operation results (exact values, "
                       "comparison outcomes and binary64 bits).")
}


struct ComparisonBenchmarkSelection {
    unsigned iterations=1000u;
    unsigned samples=11u;
    std::vector<unsigned> precisions{53u,128u,256u,1024u};
    bool summary_only=false;
    bool details=false;
    double max_slowdown=std::numeric_limits<double>::infinity();
};

List<ArgumentParser> make_comparison_argument_parsers(ComparisonBenchmarkSelection& selection) {
    List<ArgumentParser> parsers;

    parsers.append(make_valued_argument_parser(
        "i","iterations","Number of operations in each measured sample (positive)",
        [&selection](String const& value) -> VoidFunction {
            const unsigned count=BenchmarkCommandLine::parse_unsigned<unsigned>(value);
            if(count==0u) throw std::invalid_argument("iterations must be positive");
            return [&selection,count] { selection.iterations=count; };
        }));

    parsers.append(make_valued_argument_parser(
        "n","samples","Odd number of samples, at least three",
        [&selection](String const& value) -> VoidFunction {
            const unsigned count=BenchmarkCommandLine::parse_unsigned<unsigned>(value);
            if(count<3u || count%2u==0u) throw std::invalid_argument("invalid sample count");
            return [&selection,count] { selection.samples=count; };
        }));

    parsers.append(make_valued_argument_parser(
        "b","bits","Comma-separated precisions (53,128,256,1024)",
        [&selection](String const& value) -> VoidFunction {
            const std::string text(value);
            if(text.empty() || text.back()==',') throw std::invalid_argument("invalid precision list");
            std::vector<unsigned> parsed;
            std::istringstream input(text);
            std::string token;
            while(std::getline(input,token,',')) {
                const unsigned bits=BenchmarkCommandLine::parse_unsigned<unsigned>(String(token));
                if(bits!=53u && bits!=128u && bits!=256u && bits!=1024u) {
                    throw std::invalid_argument("unsupported precision");
                }
                if(std::find(parsed.begin(),parsed.end(),bits)!=parsed.end()) {
                    throw std::invalid_argument("duplicate precision");
                }
                parsed.push_back(bits);
            }
            if(parsed.empty()) throw std::invalid_argument("empty precision list");
            return [&selection,parsed] { selection.precisions=parsed; };
        }));

    parsers.append(make_unvalued_argument_parser(
        "","summary-only","Print only aggregate slowdowns",
        [&selection] { selection.summary_only=true; }));

    parsers.append(make_unvalued_argument_parser(
        "","details","Print each matched MPFR/FloatAP slowdown",
        [&selection] { selection.details=true; }));

    parsers.append(make_valued_argument_parser(
        "","max-slowdown","Fail if worst FloatAP/MPFR slowdown exceeds this ratio",
        [&selection](String const& value) -> VoidFunction {
            const std::string text(value);
            std::size_t used=0u;
            const double ratio=std::stod(text,&used);
            if(used!=text.size() || !std::isfinite(ratio) || ratio<=0.0) {
                throw std::invalid_argument("invalid positive slowdown ratio");
            }
            return [&selection,ratio] { selection.max_slowdown=ratio; };
        }));

    return parsers;
}

int main(int argc, const char* argv[]) {
    // Verbosity defaults to 0, which mutes all LOGGING_PRINTLN calls.
    // Follow the integration benchmarks and enable level-1 reports.
    BenchmarkCommandLine::initialize_logging();
    try {
        ComparisonBenchmarkSelection selection;
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution(
            "numeric-multiple-precision-benchmark",
            make_comparison_argument_parsers(selection)));

        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b","--max-slowdown",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }

        const unsigned iterations=selection.iterations;
        const unsigned samples=selection.samples;
        const auto& precisions=selection.precisions;
        const bool summary_only=selection.summary_only;
        const bool details=selection.details;
        const double max_slowdown=selection.max_slowdown;

        LOGGING_PRINTLN("FloatMP / FloatAP benchmark - running "
                        << precisions.size() << " precision groups")
        LOGGING_PRINTLN("Iterations: " << iterations << "  |  Samples: " << samples)
        LOGGING_PRINTLN(std::string(78u, '='))
        std::vector<BenchmarkDatum> measurements;
        for (unsigned bits : precisions) {
            LOGGING_PRINTLN("Measuring " << bits << " bits ...")
            run_backend<FloatMP>("mpfr",bits,iterations,samples,measurements);
            run_backend<FloatAP>("floatap",bits,iterations,samples,measurements);
            LOGGING_PRINTLN("Completed " << bits << " bits")
            LOGGING_PRINTLN(std::string(78u, '-'))
        }
        print_results(measurements,iterations,samples,summary_only);
        // No correctness work is performed inside median_ns. Check all
        // outputs before publishing the aggregate slowdown conclusions.
        verify_benchmark_values(precisions);
        print_slowdown_summary(measurements,details,max_slowdown);
        return 0;
    }
    catch (std::exception const& e) {
        LOGGING_PRINTLN("benchmark_compare_multiple_precision_backends: " << e.what())
        return 1;
    }
}

