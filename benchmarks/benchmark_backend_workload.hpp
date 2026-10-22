/***************************************************************************
 *            benchmark_backend_workload.hpp
 *
 *  Shared operation pool for FloatMP/FloatAP and FloatDP/FloatAP.
 ****************************************************************************/

#ifndef ARIADNE_NUMERIC_BENCHMARK_BACKEND_WORKLOAD_HPP
#define ARIADNE_NUMERIC_BENCHMARK_BACKEND_WORKLOAD_HPP

#include "numeric/floatdp.hpp"
#include "numeric/rational.hpp"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace AriadneBenchmark {
using namespace Ariadne;

// Keep both comparison executables on the same 27-operation battery
// and identical, deterministic 32-element exact-rational input pool.
// Changing the workload here deliberately changes BOTH benchmarks.
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

// Both the timer and its post-measurement numerical checker use the
// *same* input construction and the *same* 27 operation callables.
// A new operation cannot silently enter only one of the two paths.
template<class Float> struct BackendPool {
    using Precision = typename Float::PrecisionType;
    static constexpr unsigned pool_size = 32u;

    static Precision precision_at(unsigned bits) {
        if constexpr (std::is_same_v<Precision,DoublePrecision>) {
            if (bits != 53u) {
                throw std::invalid_argument("FloatDP is fixed to 53 bits");
            }
            return Precision{};
        } else {
            // Construct the requested precision, not a predefined fallback.
            // The benchmark may exercise 2048 bits (or any supported width).
            return Precision(bits);
        }
    }

    Precision precision;
    std::vector<Rational> exact_a, exact_b;
    std::vector<Float> a, b;

    explicit BackendPool(unsigned bits) : precision(precision_at(bits)) {
        exact_a.reserve(pool_size);
        exact_b.reserve(pool_size);
        a.reserve(pool_size);
        b.reserve(pool_size);
        const Integer base = pow(Integer(2), Nat(bits));
        // emplace_back perfectly forwards arguments as references.
        // Passing FloatDP::ROUND_TO_NEAREST directly would ODR-use its
        // static const integral member (no separate linked definition).
        // Pass a locally materialized mode value instead.
        const typename Float::RoundingModeType nearest_mode =
            Float::ROUND_TO_NEAREST;
        for (unsigned i = 0; i < pool_size; ++i) {
            exact_a.emplace_back(base + Integer(13u + 2u*i),
                                 2*base + Integer(17u + 3u*i));
            exact_b.emplace_back(base + Integer(31u + 5u*i),
                                 3*base + Integer(19u + 7u*i));
            a.emplace_back(exact_a.back(),nearest_mode,precision);
            b.emplace_back(exact_b.back(),nearest_mode,precision);
        }
    }
};

template<class Float, class Visit>
void for_each_backend_operation(BackendPool<Float> const& pool, Visit&& report) {
    constexpr unsigned pool_size = BackendPool<Float>::pool_size;
    auto const& precision = pool.precision;
    auto const& exact_a = pool.exact_a;
    auto const& a = pool.a;
    auto const& b = pool.b;
    auto index = [](unsigned i) { return i % pool_size; };
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

template<class Float>
void run_backend(char const* backend, unsigned bits, unsigned iterations,
                 unsigned samples, std::vector<BenchmarkDatum>& data) {
    const BackendPool<Float> pool(bits);
    auto report = [&](char const* name, auto&& operation) {
        const double ns = median_ns(operation, iterations, samples);
        data.push_back({backend,name,bits,ns});
    };
    for_each_backend_operation(pool,report);
}


} // namespace AriadneBenchmark
#endif // ARIADNE_NUMERIC_BENCHMARK_BACKEND_WORKLOAD_HPP
