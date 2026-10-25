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
#include <chrono>
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
        "i","iterations","Number of repetitions of the full operation battery",
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

double total_time_ns(unsigned bits, unsigned iterations, unsigned samples) {
    BackendPool<FloatAP> pool(bits);
    std::vector<double> measurements;
    measurements.reserve(samples);

    for(unsigned sample=0u; sample<=samples; ++sample) {
        const auto start=std::chrono::steady_clock::now();
        std::size_t operation_groups=0u;
        for_each_backend_operation(pool,[&](char const*, auto&& operation) {
            ++operation_groups;
            for(unsigned i=0u; i<iterations; ++i) {
                auto result=operation(i);
                keep_alive(result);
            }
        });
        const auto stop=std::chrono::steady_clock::now();
        if(operation_groups!=backend_operation_count)
            throw std::runtime_error("benchmark operation count changed");
        if(sample!=0u) {
            measurements.push_back(
                std::chrono::duration<double,std::nano>(stop-start).count());
        }
    }
    std::sort(measurements.begin(),measurements.end());
    return measurements[samples/2u];
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

        double total_ns=0.0;
        for(const unsigned bits:selection.precisions)
            total_ns+=total_time_ns(bits,selection.iterations,selection.samples);

        LOGGING_PRINTLN("FloatAP total time (all operations): "
                        << total_ns/1.0e6 << " ms")
        return 0;
    } catch(std::exception const& e) {
        LOGGING_PRINTLN("benchmark_floatap: " << e.what())
        return 1;
    }
}
