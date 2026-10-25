/***************************************************************************
 *            benchmark_floatap.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

#include "numeric/floatap.hpp"
#include "logging/logging.hpp"
#include "logging/command_line.hpp"
#include "benchmark_command_line.hpp"
#include "benchmark_backend_workload.hpp"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Ariadne;
using namespace AriadneBenchmark;

namespace {

struct Selection {
    unsigned iterations=20000u;
    unsigned samples=31u;
    std::vector<unsigned> precisions{53u,128u,256u,1024u};
};

List<ArgumentParser> make_argument_parsers(Selection& selection) {
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
        "b","bits","Comma-separated bit precisions (each >=2)",
        [&selection](String const& value) -> VoidFunction {
            const std::string text(value);
            if(text.empty() || text.back()==',') throw std::invalid_argument("invalid precision list");
            std::vector<unsigned> parsed;
            std::istringstream input(text);
            std::string token;
            while(std::getline(input,token,',')) {
                const unsigned bits=BenchmarkCommandLine::parse_unsigned<unsigned>(String(token));
                if(bits<2u) throw std::invalid_argument("precision must be at least 2 bits");
                if(std::find(parsed.begin(),parsed.end(),bits)!=parsed.end())
                    throw std::invalid_argument("duplicate precision");
                parsed.push_back(bits);
            }
            if(parsed.empty()) throw std::invalid_argument("empty precision list");
            return [&selection,parsed] { selection.precisions=parsed; };
        }));
    return parsers;
}

double print_results_and_totals(std::vector<BenchmarkDatum> const& data,
                                std::vector<unsigned> const& precisions,
                                unsigned iterations,
                                unsigned samples) {
    constexpr int operation_width=24;
    constexpr int timing_width=15;
    constexpr int table_width=operation_width+timing_width+12;

    if(data.size()!=precisions.size()*backend_operation_count)
        throw std::runtime_error("benchmark operation count changed");

    LOGGING_PRINTLN("FloatAP benchmark | iterations=" << iterations
                    << " | samples=" << samples)
    LOGGING_PRINTLN(std::string(static_cast<std::size_t>(table_width),'='))
    LOGGING_PRINTLN(std::left << std::setw(operation_width) << "Operation"
                    << "  " << std::setw(6) << "Bits"
                    << "  " << std::right << std::setw(timing_width) << "Median (ns)")
    LOGGING_PRINTLN(std::string(static_cast<std::size_t>(table_width),'-'))

    double all_precisions_total_ns=0.0;
    std::size_t offset=0u;
    for(const unsigned bits:precisions) {
        double precision_total_ns=0.0;
        std::size_t count=0u;
        while(offset<data.size() && data[offset].bits==bits) {
            const auto& datum=data[offset];
            if(datum.backend!="floatap")
                throw std::runtime_error("unexpected benchmark backend");
            LOGGING_PRINTLN(std::left << std::setw(operation_width) << datum.operation
                            << "  " << std::setw(6) << datum.bits
                            << "  " << std::right << std::setw(timing_width)
                            << std::fixed << std::setprecision(3) << datum.nanoseconds)
            precision_total_ns+=datum.nanoseconds;
            ++count;
            ++offset;
        }
        if(count!=backend_operation_count)
            throw std::runtime_error("benchmark operation count changed for precision "
                                     +std::to_string(bits));
        all_precisions_total_ns+=precision_total_ns;
        LOGGING_PRINTLN(std::string(static_cast<std::size_t>(table_width),'-'))
        LOGGING_PRINTLN(std::left << std::setw(operation_width) << "TOTAL"
                        << "  " << std::setw(6) << bits
                        << "  " << std::right << std::setw(timing_width)
                        << std::fixed << std::setprecision(3) << precision_total_ns)
        LOGGING_PRINTLN(std::string(static_cast<std::size_t>(table_width),'='))
    }
    if(offset!=data.size())
        throw std::runtime_error("unexpected benchmark rows after precision totals");

    LOGGING_PRINTLN(std::left << std::setw(operation_width) << "TOTAL ALL"
                    << "  " << std::setw(6) << "-"
                    << "  " << std::right << std::setw(timing_width)
                    << std::fixed << std::setprecision(3) << all_precisions_total_ns)
    LOGGING_PRINTLN("Totals are sums of the per-operation median nanoseconds.")
    return all_precisions_total_ns;
}

} // namespace

int main(int argc, const char* argv[]) {
    BenchmarkCommandLine::initialize_logging();
    try {
        Selection selection;
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution("numeric-floatap-benchmark",
                                                  make_argument_parsers(selection)));
        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if(!command_line.acquire(arguments))
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;

        std::vector<BenchmarkDatum> measurements;
        measurements.reserve(selection.precisions.size()*backend_operation_count);
        for(const unsigned bits:selection.precisions)
            run_backend<FloatAP>("floatap",bits,selection.iterations,selection.samples,
                                 measurements);

        print_results_and_totals(measurements,selection.precisions,
                                 selection.iterations,selection.samples);
        return 0;
    } catch(std::exception const& e) {
        LOGGING_PRINTLN("benchmark_floatap: " << e.what())
        return 1;
    }
}
