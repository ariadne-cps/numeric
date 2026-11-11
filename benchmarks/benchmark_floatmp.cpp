/***************************************************************************
 *            benchmark_floatmp.cpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric. If not, see <https://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "benchmark_backend_workload.hpp"
#include "numeric/floatmp.hpp"

#include "ariadne-benchmarking.hpp"

using namespace Ariadne;
using namespace AriadneBenchmark;

namespace {

struct Selection {
    unsigned iterations = 20000u;
    unsigned samples = 31u;
    std::vector<unsigned> precisions{53u, 128u, 256u, 1024u};
};

List<ArgumentParser> make_argument_parsers(Selection& selection) {
    List<ArgumentParser> parsers;
    parsers.append(make_valued_argument_parser(
        "i", "iterations", "Number of operations in each measured sample (positive)",
        [&selection](const String& value) -> VoidFunction {
            const unsigned count = Benchmarking::parse_unsigned<unsigned>(value);
            if (count == 0u)
                throw std::invalid_argument("iterations must be positive");
            return [&selection, count] { selection.iterations = count; };
        }));
    parsers.append(make_valued_argument_parser(
        "n", "samples", "Odd number of samples, at least three",
        [&selection](const String& value) -> VoidFunction {
            const unsigned count = Benchmarking::parse_unsigned<unsigned>(value);
            if (count < 3u || count % 2u == 0u)
                throw std::invalid_argument("invalid sample count");
            return [&selection, count] { selection.samples = count; };
        }));
    parsers.append(make_valued_argument_parser(
        "b", "bits", "Comma-separated bit precisions (each >=2)",
        [&selection](const String& value) -> VoidFunction {
            const std::string text(value);
            if (text.empty() || text.back() == ',')
                throw std::invalid_argument("invalid precision list");
            std::vector<unsigned> parsed;
            std::istringstream input(text);
            std::string token;
            while (std::getline(input, token, ',')) {
                const unsigned bits = Benchmarking::parse_unsigned<unsigned>(String(token));
                if (bits < 2u)
                    throw std::invalid_argument("precision must be at least 2 bits");
                if (std::find(parsed.begin(), parsed.end(), bits) != parsed.end())
                    throw std::invalid_argument("duplicate precision");
                parsed.push_back(bits);
            }
            if (parsed.empty())
                throw std::invalid_argument("empty precision list");
            return [&selection, parsed] { selection.precisions = parsed; };
        }));
    return parsers;
}

double print_results_and_totals(const std::vector<BenchmarkDatum>& data,
                                const std::vector<unsigned>& precisions,
                                unsigned iterations, unsigned samples,
                                const Ariadne::Benchmarking::SummaryFields& selected) {
    namespace Bench = Ariadne::Benchmarking;
    if (data.size() != precisions.size() * backend_operation_count)
        throw std::runtime_error("benchmark operation count changed");

    Bench::TableOptions style = Bench::TableOptions::ruled();
    style.group_columns = {1u};
    Bench::Report report("FloatMP benchmark", selected);
    report.set_dimension_columns({
        {"Operation"},
        {"Bits", 0u, Bench::TableAlignment::Right}
    });
    report.set_table_options(style);
    Bench::Table totals({
                            {"Group"},
                            {"Bits", 0u, Bench::TableAlignment::Right},
                            {"Sum of medians (ns)", 0u, Bench::TableAlignment::Right}
    },
                        Bench::TableOptions::ruled());
    double grand_total = 0.0;
    std::size_t offset = 0u;
    for (unsigned bits: precisions) {
        double precision_total = 0.0;
        std::size_t count = 0u;
        while (offset < data.size() && data[offset].bits == bits) {
            const auto& datum = data[offset++];
            if (datum.backend != "floatmp")
                throw std::runtime_error("unexpected benchmark backend");
            report.add_row({datum.operation, std::to_string(bits)}, datum.statistics);
            precision_total += datum.nanoseconds;
            ++count;
        }
        if (count != backend_operation_count)
            throw std::runtime_error("benchmark operation count changed for precision " + std::to_string(bits));
        grand_total += precision_total;
        totals.add_row({"TOTAL", std::to_string(bits),
                        Bench::Table::number(precision_total)});
    }
    if (offset != data.size())
        throw std::runtime_error("unexpected benchmark rows after precision groups");
    totals.add_row({"TOTAL ALL", "-", Bench::Table::number(grand_total)});
    LOGGING_PRINTLN("FloatMP benchmark | iterations=" << iterations
                                                      << " | samples=" << samples)
    report.print();
    LOGGING_PRINTLN(totals.str())
    LOGGING_PRINTLN("Totals are sums of per-operation median nanoseconds.")
    return grand_total;
}

} // namespace

int main(int argc, const char* argv[]) {
    try {
        Selection selection;
        Benchmarking::CommandLine command_line;
        command_line.add(CommandLineContribution("numeric-floatmp-benchmark",
                                                 make_argument_parsers(selection)));
        if (!command_line.acquire(argc, argv))
            return 0;

        Benchmarking::Runner runner(command_line, selection.precisions.size());
        std::vector<BenchmarkDatum> measurements;
        measurements.reserve(selection.precisions.size() * backend_operation_count);
        for (const unsigned bits: selection.precisions) {
            LOGGING_PRINTLN("Measuring " << bits << " bits ...")
            run_backend<FloatMP>("floatmp", bits, selection.iterations, selection.samples,
                                 runner, measurements);
            LOGGING_PRINTLN("Completed " << bits << " bits")
        }

        print_results_and_totals(measurements, selection.precisions,
                                 selection.iterations, selection.samples,
                                 command_line.summary_fields());
        return 0;
    } catch (const std::exception& e) {
        LOGGING_PRINTLN("benchmark_floatmp: " << e.what())
        return 1;
    }
}
