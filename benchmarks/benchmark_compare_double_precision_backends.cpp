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
#include "numeric/floatmp.hpp"
#include "numeric/rational.hpp"
#include "ariadne-benchmarking.hpp"
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
#include <utility>
#include <string>
#include <vector>

using namespace Ariadne;

#include "benchmark_backend_workload.hpp"
using namespace AriadneBenchmark;

// ComparisonResults owns all samples, validates pairing and computes
// geometric means, summed-time ratios and worst slowdowns. Report views
// select columns without changing the underlying measurements.
void print_comparison(
        Ariadne::Benchmarking::ComparisonResults const& results,
        bool summary_only, bool details, double max_slowdown,
        Ariadne::Benchmarking::SummaryFields const& selected) {
    namespace Bench=Ariadne::Benchmarking;
    Bench::ComparisonReports reports;
    if(!summary_only) {
        std::vector<std::string> columns={"implementation","operation"};
        for(auto field:selected)
            columns.push_back(std::string(Bench::summary_field_id(field)));
        reports.add({"measurements",Bench::ComparisonReportKind::Measurements,
                     std::move(columns)});
    }
    if(details) {
        reports.add({"pairs",Bench::ComparisonReportKind::Pairs,
                     {"operation","reference_ns","candidate_ns","slowdown"}});
    }
    auto aggregate=Bench::ComparisonReportDefinition(
        "aggregate slowdowns",Bench::ComparisonReportKind::Aggregates,
        {"candidate","operations","geomean","total_time",
         "worst","worst_operation"});
    aggregate.by_configuration=false;
    reports.add(std::move(aggregate));
    reports.print_all(results,"floatdp");
    for(auto const& row:results.aggregates("floatdp",false)) {
        if(std::isfinite(max_slowdown) && row.worst>max_slowdown) {
            throw std::runtime_error("FloatMP slowdown exceeds --max-slowdown: "
                +std::to_string(row.worst)+"x > "
                +std::to_string(max_slowdown)+"x");
        }
    }
}

struct ComparisonBenchmarkSelection {
    unsigned iterations=20000u;
    unsigned samples=31u;
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
            const unsigned count=Benchmarking::parse_unsigned<unsigned>(value);
            if(count==0u) throw std::invalid_argument("iterations must be positive");
            return [&selection,count] { selection.iterations=count; };
        }));

    parsers.append(make_valued_argument_parser(
        "n","samples","Odd number of samples, at least three",
        [&selection](String const& value) -> VoidFunction {
            const unsigned count=Benchmarking::parse_unsigned<unsigned>(value);
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
                const unsigned bits=Benchmarking::parse_unsigned<unsigned>(String(token));
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
        "","details","Print each matched FloatDP/FloatMP slowdown",
        [&selection] { selection.details=true; }));

    parsers.append(make_valued_argument_parser(
        "","max-slowdown","Fail if worst FloatMP/FloatDP slowdown exceeds this ratio",
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
    try {
        ComparisonBenchmarkSelection selection;
        Benchmarking::CommandLine command_line;
        command_line.add(CommandLineContribution(
            "numeric-double-precision-benchmark",
            make_comparison_argument_parsers(selection)));

        if (!command_line.acquire(argc,argv)) return 0;

        const unsigned iterations=selection.iterations;
        const unsigned samples=selection.samples;
        const auto& precisions=selection.precisions;
        const bool summary_only=selection.summary_only;
        const bool details=selection.details;
        const double max_slowdown=selection.max_slowdown;

        // FloatDP has a single supported precision (53 bits), so
        // configuration progress messages add no information here.
        LOGGING_PRINTLN("FloatDP / FloatMP comparison | iterations="
                        <<iterations<<" | samples="<<samples)
        // FloatDP is fixed to 53 bits. Register this axis even when it has
        // only one value, so comparison rows retain their configuration.
        Benchmarking::Parameter<unsigned> precision_axis(
            "bits","Floating-point precision",precisions);
        Benchmarking::Matrix matrix;
        matrix.add(precision_axis);
        Benchmarking::Runner runner(command_line,matrix.combinations());
        Benchmarking::ComparisonBattery<FloatDP,FloatMP> battery(
            {"floatdp","floatmp"});
        matrix.for_each([&](Benchmarking::Configuration const& configuration) {
            unsigned const bits=configuration.value(precision_axis);
            auto workload=[&]<class Backend>(auto&& record) {
                BackendPool<Backend> pool(bits);
                for_each_backend_operation(pool,[&](char const* name,auto&& operation) {
                    // The record callback invokes Runner synchronously while
                    // the pool and operation closure are still alive.
                    record(name,[] {},
                        [&] {
                            for(unsigned i=0u;i<iterations;++i) {
                                auto value=operation(i);
                                Benchmarking::keep_alive(value);
                            }
                        },
                        [] {},
                        static_cast<double>(iterations));
                });
            };
            battery.run(configuration,runner,samples,workload);
        });
        print_comparison(battery.results(),summary_only,details,max_slowdown,
                         command_line.summary_fields());
        return 0;
    }
    catch (std::exception const& e) {
        LOGGING_PRINTLN("benchmark_compare_double_precision_backends: " << e.what())
        return 1;
    }
}

