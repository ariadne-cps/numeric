/***************************************************************************
 *            benchmarks/benchmark_absolute_workload.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @file benchmarks/benchmark_absolute_workload.hpp
 *  @brief CLI, timing and report helpers for single-backend numeric benchmarks.
 */

#ifndef ARIADNE_NUMERIC_BENCHMARK_ABSOLUTE_WORKLOAD_HPP
#define ARIADNE_NUMERIC_BENCHMARK_ABSOLUTE_WORKLOAD_HPP

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ariadne-benchmarking.hpp"

namespace AriadneBenchmark {

namespace Bench = Ariadne::Benchmarking;
using Ariadne::Logger; // Required by LOGGING_PRINTLN's unqualified lookup.

struct AbsoluteSelection {
    unsigned iterations = 20000u;
    unsigned samples = 31u;
    std::vector<unsigned> bits{53u, 128u, 256u, 1024u, 2048u};
    bool summary_only = false;
};

// Preserve the established CLI contract: samples are odd and >=3, while
// operand widths are positive, unique and visited in the specified order.
inline Ariadne::List<Ariadne::ArgumentParser> absolute_parsers(AbsoluteSelection& s) {
    using namespace Ariadne;
    List<ArgumentParser> parsers;
    parsers.append(make_valued_argument_parser(
        "i", "iterations", "Positive operations per measured sample",
        [&s](const String& input) -> VoidFunction {
            const unsigned count = Bench::parse_unsigned<unsigned>(input);
            if (count == 0u)
                throw std::invalid_argument("iterations must be positive");
            return [&s, count] { s.iterations = count; };
        }));
    parsers.append(make_valued_argument_parser(
        "n", "samples", "Odd number of samples, at least three",
        [&s](const String& input) -> VoidFunction {
            const unsigned count = Bench::parse_unsigned<unsigned>(input);
            if (count < 3u || count % 2u == 0u)
                throw std::invalid_argument("samples must be odd and at least three");
            return [&s, count] { s.samples = count; };
        }));
    parsers.append(make_valued_argument_parser(
        "b", "bits", "Distinct operand bit widths, each at least two",
        [&s](const String& input) -> VoidFunction {
            const std::string text(input);
            if (text.empty() || text.back() == ',')
                throw std::invalid_argument("invalid bit-width list");
            std::istringstream stream(text);
            std::string token;
            std::vector<unsigned> values;
            while (std::getline(stream, token, ',')) {
                const unsigned value = Bench::parse_unsigned<unsigned>(
                    String(token));
                if (value < 2u || std::find(values.begin(), values.end(), value) != values.end())
                    throw std::invalid_argument("invalid or repeated bit width");
                values.push_back(value);
            }
            if (values.empty())
                throw std::invalid_argument("empty bit-width list");
            return [&s, values] { s.bits = values; };
        }));
    parsers.append(make_unvalued_argument_parser(
        "", "summary-only", "Print aggregate medians only",
        [&s] { s.summary_only = true; }));
    return parsers;
}

struct AbsoluteDatum {
    std::string operation;
    unsigned bits;
    Bench::Statistics statistics;
};

// The operation closure contains the same per-iteration call as the old
// benchmark. Only this closure is timed; setup, invariants and reports are not.
// With normalization=iterations, a sample is nanoseconds per operation.
template<class Operation>
void record_absolute(std::vector<AbsoluteDatum>& rows,
                     const Bench::Runner& runner,
                     const char* name, unsigned bits,
                     unsigned iterations, unsigned samples,
                     Operation&& operation) {
    Bench::Samples timings = runner.run(name, samples, [] {}, [&] {
            for(unsigned i=0u;i<iterations;++i) {
                auto value=operation(i);
                Bench::keep_alive(value);
            } }, [] {}, static_cast<double>(iterations));
    auto stats = Bench::summarize(timings);
    if (!std::isfinite(stats.median) || stats.median <= 0.0)
        throw std::runtime_error(std::string("invalid timing for ") + name);
    rows.push_back({name, bits, std::move(stats)});
}

inline void report_absolute(const char* type,
                            const AbsoluteSelection& s,
                            const Bench::SummaryFields& fields,
                            const std::vector<AbsoluteDatum>& rows) {
    Bench::Report report(std::string(type) + " operation times", fields);
    report.set_dimension_columns({
        {"Operation"},
        {"Bits", 0u, Bench::TableAlignment::Right}
    });
    auto layout = Bench::TableOptions::ruled();
    layout.group_columns = {1u};
    report.set_table_options(layout);

    // Validate that no bit-width group silently loses an operation, and
    // retain the established sum-of-medians aggregate (not a measured total).
    std::map<unsigned, std::pair<double, std::size_t>> totals;
    for (const auto& row: rows) {
        auto& group = totals[row.bits];
        group.first += row.statistics.median;
        ++group.second;
        if (!s.summary_only)
            report.add_row({row.operation, std::to_string(row.bits)}, row.statistics);
    }
    if (!s.summary_only) {
        LOGGING_PRINTLN(type << " | iterations=" << s.iterations
                             << " | samples=" << s.samples)
        report.print();
    }
    Bench::Table summary({
                             {"Group"},
                             {"Bits", 0u, Bench::TableAlignment::Right},
                             {"Operations", 0u, Bench::TableAlignment::Right},
                             {"Sum of medians (ns)", 0u, Bench::TableAlignment::Right}
    },
                         Bench::TableOptions::ruled());
    double grand_total = 0.0;
    for (const auto& [bits, total]: totals) {
        grand_total += total.first;
        summary.add_row({"TOTAL", std::to_string(bits),
                         std::to_string(total.second),
                         Bench::Table::number(total.first)});
    }
    summary.add_row({"TOTAL ALL", "-", "-",
                     Bench::Table::number(grand_total)});
    LOGGING_PRINTLN(summary.str())
    LOGGING_PRINTLN("Totals sum per-operation median durations (not wall-clock totals).")
}

inline bool acquire_absolute(AbsoluteSelection& selection,
                             Bench::CommandLine& cli,
                             const char* title, int argc, const char* argv[]) {
    using namespace Ariadne;
    cli.add(CommandLineContribution(title, absolute_parsers(selection)));
    // The shared command-line parser already handles named options and help.
    return cli.acquire(argc, argv);
}

} // namespace AriadneBenchmark

#endif // ARIADNE_NUMERIC_BENCHMARK_ABSOLUTE_WORKLOAD_HPP
