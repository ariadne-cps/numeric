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

template< class T > inline void keep_alive(T const& value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(& value) : "memory");
#else
    static void const* volatile sink = nullptr;
    sink =& value;
#endif
}

template< class Operation >
double median_ns(Operation && operation, unsigned iterations, unsigned samples) {
    std::vector<double> measurements;
    measurements.reserve(samples);
    for (unsigned sample = 0; sample <= samples;++sample) {
        auto start = std::chrono::steady_clock::now();
        for (unsigned i = 0; i < iterations;++i) {
            auto result = operation(i);
            keep_alive(result);
        }
        auto stop = std::chrono::steady_clock::now();
        if (sample != 0) {
            // unmeasured warm-up
            measurements.push_back(
            std::chrono::duration < double, std::nano >(stop - start).count() / iterations);
        }
    }
    std::sort(measurements.begin(), measurements.end());
    return measurements[samples / 2];
}

struct BenchmarkDatum {
    std::string backend;
    std::string operation;
    unsigned bits;
    double nanoseconds;
};

template< class Float >
void run_backend(char const* backend, unsigned bits, unsigned iterations,
unsigned samples, std::vector<BenchmarkDatum>& data) {
    using Precision = typename Float::PrecisionType;
    Precision precision = bits == 53 ? Precision(53_bits)
    : bits == 128 ? Precision(128_bits)
    : bits == 256 ? Precision(256_bits) : Precision(1024_bits);
    constexpr unsigned pool_size = 32;
    std::vector<Rational> exact_a, exact_b;
    std::vector<Float> a, b;
    exact_a.reserve(pool_size);
    exact_b.reserve(pool_size);
    a.reserve(pool_size);
    b.reserve(pool_size);

    // Nontrivial significands: avoid benchmarking only exact small integers.
    const Integer base = pow(Integer(2), Nat(bits));
    for (unsigned i = 0; i < pool_size;++i) {
        exact_a.emplace_back(base + Integer(13u + 2u* i), 2* base + Integer(17u + 3u* i));
        exact_b.emplace_back(base + Integer(31u + 5u* i), 3* base + Integer(19u + 7u* i));
        a.emplace_back(exact_a.back(), Float::ROUND_TO_NEAREST, precision);
        b.emplace_back(exact_b.back(), Float::ROUND_TO_NEAREST, precision);
    }

    auto report =[& ](char const* name, auto && operation) {
        const double ns = median_ns(operation, iterations, samples);
        data.push_back( { backend, name, bits, ns });

    };
    auto index =[](unsigned i) {
        return i % pool_size;
    };
    report("construct_rational",[& ](unsigned i) {
        return Float(exact_a[index(i)], Float::ROUND_TO_NEAREST, precision);
    }
    );
    report("copy",[& ](unsigned i) {
        return Float(a[index(i)]);
    }
    );
    report("compare",[& ](unsigned i) {
        return cmp(a[index(i)], b[index(i)]);
    }
    );
    report("add_nearest",[& ](unsigned i) {
        return add(Float::ROUND_TO_NEAREST, a[index(i)], b[index(i)]);
    }
    );
    report("add_down",[& ](unsigned i) {
        return add(Float::ROUND_DOWNWARD, a[index(i)], b[index(i)]);
    }
    );
    report("add_up",[& ](unsigned i) {
        return add(Float::ROUND_UPWARD, a[index(i)], b[index(i)]);
    }
    );
    report("sub_nearest",[& ](unsigned i) {
        return sub(Float::ROUND_TO_NEAREST, a[index(i)], b[index(i)]);
    }
    );
    report("mul_nearest",[& ](unsigned i) {
        return mul(Float::ROUND_TO_NEAREST, a[index(i)], b[index(i)]);
    }
    );
    report("mul_down",[& ](unsigned i) {
        return mul(Float::ROUND_DOWNWARD, a[index(i)], b[index(i)]);
    }
    );
    report("mul_up",[& ](unsigned i) {
        return mul(Float::ROUND_UPWARD, a[index(i)], b[index(i)]);
    }
    );
    report("div_nearest",[& ](unsigned i) {
        return div(Float::ROUND_TO_NEAREST, a[index(i)], b[index(i)]);
    }
    );
    report("div_down",[& ](unsigned i) {
        return div(Float::ROUND_DOWNWARD, a[index(i)], b[index(i)]);
    }
    );
    report("div_up",[& ](unsigned i) {
        return div(Float::ROUND_UPWARD, a[index(i)], b[index(i)]);
    }
    );
    report("fma_nearest",[& ](unsigned i) {
        return fma(Float::ROUND_TO_NEAREST, a[index(i)], b[index(i)], a[index(i)]);
    }
    );
    // Operations shared by both concrete backends are timed identically.
    report("sqrt_nearest",[& ](unsigned i) {
        return sqrt(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("exp_nearest",[& ](unsigned i) {
        return exp(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("log_nearest",[& ](unsigned i) {
        return log(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("sin_nearest",[& ](unsigned i) {
        return sin(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("cos_nearest",[& ](unsigned i) {
        return cos(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("tan_nearest",[& ](unsigned i) {
        return tan(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("atan_nearest",[& ](unsigned i) {
        return atan(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("asin_nearest", [&](unsigned i) {
        return asin(Float::ROUND_TO_NEAREST, a[index(i)]);
    });
    report("acos_nearest", [&](unsigned i) {
        return acos(Float::ROUND_TO_NEAREST, a[index(i)]);
    });
    report("pi_nearest", [&](unsigned) {
        return Float::pi(Float::ROUND_TO_NEAREST, precision);
    });
    report("pow7_nearest",[& ](unsigned i) {
        return pow(Float::ROUND_TO_NEAREST, a[index(i)], Nat(7u));
    }
    );
    report("to_double",[& ](unsigned i) {
        return a[index(i)].get_d();
    }
    );
    report("to_rational",[& ](unsigned i) {
        return static_cast<Rational>(a[index(i)]);
    }
    );
}

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

int main(int argc, char* * argv) {
    // Verbosity defaults to 0, which mutes all LOGGING_PRINTLN calls.
    // Follow the integration benchmarks and enable level-1 reports.
    Logger::instance().configuration().set_verbosity(1u);
    try {
        unsigned iterations = 300;
        unsigned samples = 7;
        bool summary_only = false;
        bool details = false;
        std::vector<unsigned> precisions {53u,128u,256u,1024u};
        double max_slowdown = std::numeric_limits<double>::infinity();
        unsigned positional = 0;

        for (int i = 1; i < argc;++i) {
            const std::string arg(argv[i]);
            if (arg == "--iterations" || arg == "--samples" || arg == "--bits") {
                if (++i >= argc) throw std::invalid_argument(arg + " requires a value");
                const std::string value(argv[i]);
                if (arg == "--bits") {
                    precisions.clear();
                    std::istringstream input(value);
                    std::string token;
                    while (std::getline(input,token,',')) {
                        std::size_t consumed=0;
                        const unsigned long bits=std::stoul(token,&consumed);
                        if (consumed!=token.size() || (bits!=53u && bits!=128u && bits!=256u && bits!=1024u)) {
                            throw std::invalid_argument("--bits supports 53,128,256,1024");
                        }
                        if (std::find(precisions.begin(),precisions.end(),static_cast<unsigned>(bits))!=precisions.end()) {
                            throw std::invalid_argument("duplicate --bits precision");
                        }
                        precisions.push_back(static_cast<unsigned>(bits));
                    }
                    if (precisions.empty()) throw std::invalid_argument("--bits requires precisions");
                } else {
                    std::size_t consumed=0;
                    const unsigned long parsed=std::stoul(value,&consumed);
                    if (consumed!=value.size() || parsed>std::numeric_limits<unsigned>::max()) {
                        throw std::invalid_argument("invalid count: " + value);
                    }
                    if (arg == "--iterations") iterations=static_cast<unsigned>(parsed);
                    else samples=static_cast<unsigned>(parsed);
                }
            }
            else if (arg == "--summary-only") {
                summary_only = true;
            }
            else if (arg == "--details") {
                details = true;
            }
            else if (arg == "--max-slowdown") {
                if (++i >= argc) {
                    throw std::invalid_argument("--max-slowdown requires a positive ratio");
                }
                const std::string value(argv[i]);
                std::size_t read = 0;
                max_slowdown = std::stod(value,& read);
                if (read != value.size() || !std::isfinite(max_slowdown)
                || max_slowdown <= 0.0) {
                    throw std::invalid_argument("--max-slowdown must be positive and finite");
                }
            }
            else if (arg == "--help") {
                LOGGING_PRINTLN("Usage: benchmark_compare_multiple_precision_backends "
                    "[iterations] [odd_samples] [--iterations N] [--samples ODD_N] "
                    "[--bits 53,128,256,1024] [--summary-only] [--details] "
                    "[--max-slowdown RATIO]")
                return 0;
            }
            else if (arg.starts_with("--")) {
                throw std::invalid_argument("unknown option: " + arg);
            } else {
                std::size_t read = 0;
                unsigned long parsed = std::stoul(arg,& read);
                if (read != arg.size() || parsed > std::numeric_limits<unsigned>::max()) {
                    throw std::invalid_argument("invalid iterations/sample count: " + arg);
                }
                if (positional == 0) {
                    iterations = static_cast<unsigned>(parsed);
                }
                else if (positional == 1) {
                    samples = static_cast<unsigned>(parsed);
                } else {
                    throw std::invalid_argument("too many positional arguments");
                }
                ++positional;
            }
        }
        if (iterations == 0 || samples < 3 || samples % 2 == 0) {
            throw std::invalid_argument("iterations must be positive; samples must be odd and >=3");
        }

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
        print_slowdown_summary(measurements,details,max_slowdown);
        return 0;
    }
    catch (std::exception const& e) {
        LOGGING_PRINTLN("benchmark_compare_multiple_precision_backends: " << e.what())
        return 1;
    }
}

