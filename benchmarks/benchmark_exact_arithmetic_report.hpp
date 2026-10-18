/***************************************************************************
 *            benchmark_exact_arithmetic_report.hpp
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

/*! \file benchmark_exact_arithmetic_report.hpp
 *  \brief Reporting and correctness checks for exact arithmetic benchmarks.
 */

#ifndef ARIADNE_NUMERIC_BENCHMARK_EXACT_ARITHMETIC_REPORT_HPP
#define ARIADNE_NUMERIC_BENCHMARK_EXACT_ARITHMETIC_REPORT_HPP

#include "benchmark_exact_arithmetic_workload.hpp"
#include "benchmark_command_line.hpp"

#include "logging/logging.hpp"
#include "logging/command_line.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace AriadneBenchmark {

inline std::string exact_value(std::string const& value) {
    return value;
}

struct ExactArithmeticSelection {
    unsigned iterations=1000u;
    unsigned samples=11u;
    std::vector<unsigned> sizes{53u,128u,256u,1024u};
    bool summary_only=false;
    bool details=false;
    double max_slowdown=std::numeric_limits<double>::infinity();
};

inline List<ArgumentParser> exact_arithmetic_parsers(
        ExactArithmeticSelection& selection,bool comparison)
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
        "","summary-only","Print only aggregate statistics",
        [&selection] { selection.summary_only=true; }));
    parsers.append(make_unvalued_argument_parser(
        "","details","Print individual timings and matched ratios",
        [&selection] { selection.details=true; }));
    if (comparison) {
        parsers.append(make_valued_argument_parser(
            "","max-slowdown","Fail if worst native/GMP slowdown exceeds this ratio",
            [&selection](String const& value) -> VoidFunction {
                const std::string text(value);
                std::size_t used=0u;
                const double ratio=std::stod(text,&used);
                if (used!=text.size() || !std::isfinite(ratio) || ratio<=0.0) {
                    throw std::invalid_argument("invalid positive slowdown ratio");
                }
                return [&selection,ratio] { selection.max_slowdown=ratio; };
            }));
    }
    return parsers;
}

inline void print_exact_measurements(std::vector<ExactArithmeticDatum> const& data,
                                     char const* title,unsigned iterations,
                                     unsigned samples,bool summary_only)
{
    if (summary_only) { return; }
    LOGGING_PRINTLN(title << " benchmark | iterations=" << iterations
                    << " | samples=" << samples)
    LOGGING_PRINTLN(std::string(79u,'='))
    LOGGING_PRINTLN(std::left << std::setw(10) << "Backend"
                    << "  " << std::setw(24) << "Operation"
                    << "  " << std::setw(6) << "Bits"
                    << "  " << std::right << std::setw(15) << "Median (ns)")
    LOGGING_PRINTLN(std::string(79u,'-'))
    for (auto const& datum : data) {
        LOGGING_PRINTLN(std::left << std::setw(10) << datum.backend
                        << "  " << std::setw(24) << datum.operation
                        << "  " << std::setw(6) << datum.bits
                        << "  " << std::right << std::setw(15)
                        << std::fixed << std::setprecision(3)
                        << datum.nanoseconds)
    }
    LOGGING_PRINTLN(std::string(79u,'='))
}

struct ExactArithmeticAggregate {
    std::size_t count=0u;
    double log_sum=0.0;
    double reference_sum=0.0;
    double candidate_sum=0.0;
    double worst_ratio=0.0;
    std::string worst_operation;

    void add(std::string const& name,double reference,double candidate) {
        if (!(std::isfinite(reference) && std::isfinite(candidate)) ||
            reference<=0.0 || candidate<=0.0) {
            throw std::runtime_error("nonpositive or nonfinite benchmark timing");
        }
        const double ratio=candidate/reference;
        ++count;
        log_sum+=std::log(ratio);
        reference_sum+=reference;
        candidate_sum+=candidate;
        if (ratio>worst_ratio) {
            worst_ratio=ratio;
            worst_operation=name;
        }
    }

    void print(std::string const& label) const {
        if (count==0u || reference_sum<=0.0) {
            throw std::runtime_error("empty arithmetic benchmark aggregate");
        }
        LOGGING_PRINTLN(std::left << std::setw(12) << label
                        << "  ops=" << std::setw(3) << count
                        << "  geomean=" << std::fixed << std::setprecision(3)
                        << std::exp(log_sum/static_cast<double>(count)) << "x"
                        << "  total_time=" << candidate_sum/reference_sum << "x"
                        << "  worst=" << worst_ratio << "x ("
                        << worst_operation << ")")
    }
};

inline void print_exact_comparison(
        std::vector<ExactArithmeticDatum> const& data,
        char const* native_backend,char const* reference_backend,
        char const* title,bool details,double maximum_slowdown)
{
    using Key=std::pair<unsigned,std::string>;
    std::map<Key,double> reference,candidate;
    for (auto const& datum : data) {
        auto& target=datum.backend==reference_backend ? reference : candidate;
        if (datum.backend!=reference_backend && datum.backend!=native_backend) {
            throw std::runtime_error("unrecognized exact benchmark backend");
        }
        if (!target.emplace(Key{datum.bits,datum.operation},datum.nanoseconds).second) {
            throw std::runtime_error("duplicate exact benchmark timing");
        }
    }
    if (reference.empty() || reference.size()!=candidate.size()) {
        throw std::runtime_error("unmatched arithmetic benchmark timings");
    }
    std::map<unsigned,ExactArithmeticAggregate> grouped;
    ExactArithmeticAggregate overall;
    if (details) {
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("Detailed comparisons (native / GMP)")
        LOGGING_PRINTLN(std::string(79u,'-'))
        LOGGING_PRINTLN(std::left << std::setw(24) << "Operation"
                        << std::setw(8) << "Bits"
                        << std::setw(16) << "GMP (ns)"
                        << std::setw(16) << "Native (ns)" << "Slowdown")
    }
    for (auto const& [key,reference_ns] : reference) {
        const auto found=candidate.find(key);
        if (found==candidate.end()) {
            throw std::runtime_error("missing native operation: "+key.second);
        }
        const std::string name=key.second+", "+std::to_string(key.first)+" bits";
        grouped[key.first].add(name,reference_ns,found->second);
        overall.add(name,reference_ns,found->second);
        if (details) {
            LOGGING_PRINTLN(std::left << std::setw(24) << key.second
                            << std::setw(8) << key.first
                            << std::setw(16) << std::fixed << std::setprecision(3)
                            << reference_ns << std::setw(16) << found->second
                            << std::setprecision(3) << found->second/reference_ns << "x")
        }
    }
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN(std::string(79u,'='))
    LOGGING_PRINTLN(title << " / GMP slowdown (x; >1 means slower)")
    LOGGING_PRINTLN(std::string(79u,'-'))
    for (auto const& [bits,aggregate] : grouped) {
        aggregate.print(std::to_string(bits)+" bits");
    }
    overall.print("ALL");
    LOGGING_PRINTLN(std::string(79u,'='))
    LOGGING_PRINTLN("geomean weights each (operation, size) equally; "
                    "total_time weights by measured cost.")
    if (std::isfinite(maximum_slowdown) &&
        overall.worst_ratio>maximum_slowdown) {
        throw std::runtime_error(
            "native slowdown exceeds --max-slowdown: "+
            std::to_string(overall.worst_ratio)+"x > "+
            std::to_string(maximum_slowdown)+"x");
    }
}

inline void print_exact_absolute(
        std::vector<ExactArithmeticDatum> const& data,char const* title,
        bool details)
{
    std::map<unsigned,std::pair<double,std::size_t>> totals;
    for (auto const& datum : data) {
        if (!(std::isfinite(datum.nanoseconds)) || datum.nanoseconds<=0.0) {
            throw std::runtime_error("invalid absolute benchmark time");
        }
        auto& group=totals[datum.bits];
        group.first+=datum.nanoseconds;
        ++group.second;
        if (details) {
            LOGGING_PRINTLN(std::left << std::setw(24) << datum.operation
                            << std::setw(8) << datum.bits
                            << std::fixed << std::setprecision(3)
                            << datum.nanoseconds << " ns")
        }
    }
    LOGGING_PRINTLN(std::string(79u,'='))
    LOGGING_PRINTLN(title << " absolute timing summary")
    LOGGING_PRINTLN(std::string(79u,'-'))
    for (auto const& [bits,total] : totals) {
        LOGGING_PRINTLN(std::left << std::setw(12)
                        << std::to_string(bits)+" bits"
                        << "  ops=" << total.second
                        << "  sum_of_medians=" << std::fixed
                        << std::setprecision(3) << total.first << " ns")
    }
    LOGGING_PRINTLN(std::string(79u,'='))
}

template<ExactArithmeticKind Kind,class Operations>
void measure_exact_backend(char const* backend,unsigned bits,
                           unsigned iterations,unsigned samples,
                           std::vector<ExactArithmeticDatum>& data)
{
    const auto pool=make_exact_pool<Kind,Operations>(bits);
    for_each_exact_operation<Kind,Operations>(
        pool,[&](char const* name,auto&& operation) {
            const double time=exact_arithmetic_median_ns(
                operation,iterations,samples);
            data.push_back({backend,name,bits,time});
        });
}

template<ExactArithmeticKind Kind,class NativeOperations,class GmpOperations>
void verify_exact_results(std::vector<unsigned> const& sizes) {
    std::size_t checked_inputs=0u;
    std::size_t checked_operations=0u;
    for (unsigned bits : sizes) {
        const auto gmp=make_exact_pool<Kind,GmpOperations>(bits);
        const auto native=make_exact_pool<Kind,NativeOperations>(bits);
        constexpr unsigned pool_size=32u;
        for (unsigned i=0u; i<pool_size; ++i) {
            if (exact_value(gmp.a[i])!=exact_value(native.a[i]) ||
                exact_value(gmp.b[i])!=exact_value(native.b[i])) {
                throw std::runtime_error("input mismatch at "+
                                         std::to_string(bits)+" bits, index "+
                                         std::to_string(i));
            }
            checked_inputs+=2u;
        }
        std::map<std::string,std::vector<std::string>> reference;
        for_each_exact_operation<Kind,GmpOperations>(
            gmp,[&](char const* name,auto&& operation) {
                auto [it,inserted]=reference.try_emplace(name);
                if (!inserted) {
                    throw std::runtime_error(std::string("duplicate GMP operation: ")+name);
                }
                for (unsigned i=0u; i<pool_size; ++i) {
                    it->second.push_back(exact_value(operation(i)));
                }
            });
        std::size_t count=0u;
        for_each_exact_operation<Kind,NativeOperations>(
            native,[&](char const* name,auto&& operation) {
                const auto found=reference.find(name);
                if (found==reference.end() || found->second.size()!=pool_size) {
                    throw std::runtime_error(std::string("missing GMP operation: ")+name);
                }
                for (unsigned i=0u; i<pool_size; ++i) {
                    const std::string actual=exact_value(operation(i));
                    if (actual!=found->second[i]) {
                        std::ostringstream error;
                        error << "Numerical mismatch: " << bits << " bits, "
                              << name << ", index " << i
                              << "; GMP=" << found->second[i]
                              << "; native=" << actual;
                        throw std::runtime_error(error.str());
                    }
                    ++checked_operations;
                }
                ++count;
            });
        constexpr std::size_t expected_operations=
            Kind==ExactArithmeticKind::Integer ? 12u : 15u;
        if (reference.size()!=count || count!=expected_operations) {
            throw std::runtime_error("unexpected exact benchmark operation count");
        }
    }
    LOGGING_PRINTLN("Numerical verification PASS: "
                    << checked_inputs << " paired inputs and "
                    << checked_operations << " paired operation results.")
}

template<ExactArithmeticKind Kind,bool Comparison>
int run_exact_arithmetic_benchmark(int argc,const char* argv[]) {
    BenchmarkCommandLine::initialize_logging();
    constexpr bool integer=Kind==ExactArithmeticKind::Integer;
    const char* title=integer ? "Integer" : "Rational";
    const char* native_name=integer ? "integer" : "rational";
    const char* gmp_name=integer ? "mpz" : "mpq";

    try {
        ExactArithmeticSelection selection;
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution(
            integer ? (Comparison ? "numeric-integer-comparison" :
                                    "numeric-integer-absolute") :
                      (Comparison ? "numeric-rational-comparison" :
                                    "numeric-rational-absolute"),
            exact_arithmetic_parsers(selection,Comparison)));

        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--max-slowdown","--scheduler","-s","--theme","-t",
             "--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }
        LOGGING_PRINTLN(title << (Comparison ? " vs GMP" : " absolute")
                        << " benchmark | iterations=" << selection.iterations
                        << " | samples=" << selection.samples)
        std::vector<ExactArithmeticDatum> data;
        for (unsigned bits : selection.sizes) {
            LOGGING_PRINTLN("Measuring " << bits << " bits ...")
            if constexpr (integer) {
                if constexpr (Comparison) {
                    measure_exact_backend<Kind,GmpIntegerOperations>(
                        gmp_name,bits,selection.iterations,selection.samples,data);
                }
                measure_exact_backend<Kind,IntegerOperations>(
                    native_name,bits,selection.iterations,selection.samples,data);
            } else {
                if constexpr (Comparison) {
                    measure_exact_backend<Kind,GmpRationalOperations>(
                        gmp_name,bits,selection.iterations,selection.samples,data);
                }
                measure_exact_backend<Kind,RationalOperations>(
                    native_name,bits,selection.iterations,selection.samples,data);
            }
            LOGGING_PRINTLN("Completed " << bits << " bits")
        }
        print_exact_measurements(data,title,selection.iterations,
                                 selection.samples,selection.summary_only);
        // Match the FloatAP comparison: numerical verification is entirely
        // outside measurement loops and precedes slowdown reporting.
        if constexpr (Comparison) {
            if constexpr (integer) {
                verify_exact_results<Kind,IntegerOperations,GmpIntegerOperations>(
                    selection.sizes);
            } else {
                verify_exact_results<Kind,RationalOperations,GmpRationalOperations>(
                    selection.sizes);
            }
            print_exact_comparison(data,native_name,gmp_name,title,
                                   selection.details,selection.max_slowdown);
        } else {
            print_exact_absolute(data,title,selection.details);
        }
        return 0;
    }
    catch (std::exception const& e) {
        LOGGING_PRINTLN("exact arithmetic benchmark: " << e.what())
        return 1;
    }
}

} // namespace AriadneBenchmark

#endif // ARIADNE_NUMERIC_BENCHMARK_EXACT_ARITHMETIC_REPORT_HPP
