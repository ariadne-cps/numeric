/***************************************************************************
 *            benchmark_compare_double_precision_backends.cpp
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
#include "numeric/floatdp.hpp"
#include "numeric/floatap.hpp"
#include "numeric/rational.hpp"
#include "logging/logging.hpp"
#include "logging/command_line.hpp"
#include "benchmark_command_line.hpp"
#include <algorithm>
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
    LOGGING_PRINTLN("FloatDP / FloatAP benchmark | iterations=" << iterations
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
        double floatdp_sum = 0.0;
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
            floatdp_sum += reference;
            floatap_sum += candidate;
            if (ratio > worst_ratio) {
                worst_ratio = ratio;
                worst_operation = name;
                worst_bits = bits;
            }
        }
        void print(std::string const& name) const {
            if (count == 0 || floatdp_sum <= 0.0) {
                throw std::runtime_error("cannot summarize empty benchmark results");
            }
            LOGGING_PRINTLN(std::left << std::setw(12) << name
                << "  ops=" << std::setw(3) << count
                << "  geomean=" << std::fixed << std::setprecision(3)
                << std::exp(logarithmic_sum / static_cast<double>(count)) << "x"
                << "  total_time=" << floatap_sum / floatdp_sum << "x"
                << "  worst=" << worst_ratio << "x (" << worst_operation
                << ", " << worst_bits << " bits)")
        }
    };

    using Key = std::pair < unsigned, std::string >;
    std::map< Key, double > reference, candidate;
    for (auto const& datum : data) {
        auto& collection = datum.backend == "floatdp" ? reference : candidate;
        if (!collection.emplace(Key { datum.bits, datum.operation }, datum.nanoseconds).second) {
            throw std::runtime_error("duplicate benchmark measurement");
        }
    }
    if (reference.size() != candidate.size() || reference.empty()) {
        throw std::runtime_error("FloatAP and FloatDP benchmark rows are not paired");
    }
    std::map< unsigned, Aggregate > per_precision;
    Aggregate overall;
    if (details) {
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("Detailed comparisons (FloatAP / FloatDP)")
        LOGGING_PRINTLN(std::string(79u,'-'))
        LOGGING_PRINTLN(std::left << std::setw(24) << "Operation"
                        << std::setw(8) << "Bits"
                        << std::setw(16) << "FloatDP (ns)"
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
    LOGGING_PRINTLN("FloatAP / FloatDP slowdown (x; >1 means slower)")
    LOGGING_PRINTLN(std::string(78u, '-'))
    for (auto const& [bits, summary] : per_precision) {
        summary.print(std::to_string(bits) + " bits");
    }
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


struct ComparisonBenchmarkSelection {
    unsigned iterations=1000u;
    unsigned samples=11u;
    std::vector<unsigned> precisions{53u};
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
        "b","bits","Precision (FloatDP supports 53 bits only)",
        [&selection](String const& value) -> VoidFunction {
            const std::string text(value);
            if(text.empty() || text.back()==',') throw std::invalid_argument("invalid precision list");
            std::vector<unsigned> parsed;
            std::istringstream input(text);
            std::string token;
            while(std::getline(input,token,',')) {
                const unsigned bits=BenchmarkCommandLine::parse_unsigned<unsigned>(String(token));
                if(bits!=53u) {
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
        "","details","Print each matched FloatDP/FloatAP slowdown",
        [&selection] { selection.details=true; }));

    parsers.append(make_valued_argument_parser(
        "","max-slowdown","Fail if worst FloatAP/FloatDP slowdown exceeds this ratio",
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
            "numeric-double-precision-benchmark",
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

        LOGGING_PRINTLN("FloatDP / FloatAP benchmark - running "
                        << precisions.size() << " precision groups")
        LOGGING_PRINTLN("Iterations: " << iterations << "  |  Samples: " << samples)
        LOGGING_PRINTLN(std::string(78u, '='))
        std::vector<BenchmarkDatum> measurements;
        for (unsigned bits : precisions) {
            LOGGING_PRINTLN("Measuring " << bits << " bits ...")
            run_backend<FloatDP>("floatdp",bits,iterations,samples,measurements);
            run_backend<FloatAP>("floatap",bits,iterations,samples,measurements);
            LOGGING_PRINTLN("Completed " << bits << " bits")
            LOGGING_PRINTLN(std::string(78u, '-'))
        }
        print_results(measurements,iterations,samples,summary_only);
        print_slowdown_summary(measurements,details,max_slowdown);
        return 0;
    }
    catch (std::exception const& e) {
        LOGGING_PRINTLN("benchmark_compare_double_precision_backends: " << e.what())
        return 1;
    }
}

