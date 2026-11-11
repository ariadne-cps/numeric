/***************************************************************************
 *            benchmarks/benchmark_backend_workload.hpp
 *
 *  Copyright  2026
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

/** @file benchmarks/benchmark_backend_workload.hpp
 *  @brief Shared benchmark workloads for single- and multiple-precision numeric backends.
 */

#ifndef ARIADNE_NUMERIC_BENCHMARK_BACKEND_WORKLOAD_HPP
#define ARIADNE_NUMERIC_BENCHMARK_BACKEND_WORKLOAD_HPP

#include "numeric/floatdp.hpp"
#include "numeric/rational.hpp"
#include <algorithm>
#include "ariadne-benchmarking.hpp"
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace AriadneBenchmark {
using namespace Ariadne;

inline constexpr unsigned backend_operation_count=62u;

// Keep both comparison executables on the same operation battery
// and identical, deterministic 32-element exact-rational input pool.
// Changing the workload here deliberately changes all backend benchmarks.
template<class T> inline void keep_alive(T const& value) {
    Ariadne::Benchmarking::keep_alive(value);
}

template<class Operation>
Ariadne::Benchmarking::Statistics measure_operation(
        Ariadne::Benchmarking::Runner const& runner,
        char const* name, Operation&& operation,
        unsigned iterations, unsigned samples) {
    // Preserve the original loop and its 32-entry deterministic operand pool.
    // The runner now owns warm-ups, timing boundaries, and normalized samples.
    auto values = runner.run(name, samples,
        [] {},
        [&] {
            for (unsigned i=0u; i<iterations; ++i) {
                auto result = operation(i);
                Ariadne::Benchmarking::keep_alive(result);
            }
        },
        [] {},
        static_cast<double>(iterations));
    return Ariadne::Benchmarking::summarize(values);
}

struct BenchmarkDatum {
    std::string backend;
    std::string operation;
    unsigned bits;
    double nanoseconds;
    Ariadne::Benchmarking::Statistics statistics;
};

// The same deterministic input construction and operation callables are
// used by the standalone FloatMP and paired FloatDP/FloatMP drivers.
// The benchmark clock is owned by Runner, not by this workload header.
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
            // Construct precisely the requested bit count: the comparison
            // CLI accepts every integer precision >=2, not a fixed preset list.
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
    report("sub_down",[& ](unsigned i) {
        return sub(Float::ROUND_DOWNWARD, a[index(i)], b[index(i)]);
    }
    );
    report("sub_up",[& ](unsigned i) {
        return sub(Float::ROUND_UPWARD, a[index(i)], b[index(i)]);
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
    report("fma_down",[& ](unsigned i) {
        return fma(Float::ROUND_DOWNWARD, a[index(i)], b[index(i)], a[index(i)]);
    }
    );
    report("fma_up",[& ](unsigned i) {
        return fma(Float::ROUND_UPWARD, a[index(i)], b[index(i)], a[index(i)]);
    }
    );
    // Operations shared by both concrete backends are timed identically.
    report("sqrt_nearest",[& ](unsigned i) {
        return sqrt(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("sqrt_down",[& ](unsigned i) {
        return sqrt(Float::ROUND_DOWNWARD, a[index(i)]);
    }
    );
    report("sqrt_up",[& ](unsigned i) {
        return sqrt(Float::ROUND_UPWARD, a[index(i)]);
    }
    );
    report("exp_nearest",[& ](unsigned i) {
        return exp(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("exp_down",[& ](unsigned i) {
        return exp(Float::ROUND_DOWNWARD, a[index(i)]);
    }
    );
    report("exp_up",[& ](unsigned i) {
        return exp(Float::ROUND_UPWARD, a[index(i)]);
    }
    );
    report("log_nearest",[& ](unsigned i) {
        return log(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("log_down",[&](unsigned i) {
        return log(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("log_up",[&](unsigned i) {
        return log(Float::ROUND_UPWARD,a[index(i)]);
    });
    report("sin_nearest",[& ](unsigned i) {
        return sin(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("sin_down",[&](unsigned i) {
        return sin(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("sin_up",[&](unsigned i) {
        return sin(Float::ROUND_UPWARD,a[index(i)]);
    });
    report("cos_nearest",[& ](unsigned i) {
        return cos(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("cos_down",[&](unsigned i) {
        return cos(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("cos_up",[&](unsigned i) {
        return cos(Float::ROUND_UPWARD,a[index(i)]);
    });
    report("tan_nearest",[& ](unsigned i) {
        return tan(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("tan_down",[&](unsigned i) {
        return tan(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("tan_up",[&](unsigned i) {
        return tan(Float::ROUND_UPWARD,a[index(i)]);
    });
    report("tanh_nearest",[&](unsigned i) {
        return tanh(Float::ROUND_TO_NEAREST,a[index(i)]);
    });
    report("tanh_down",[&](unsigned i) {
        return tanh(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("tanh_up",[&](unsigned i) {
        return tanh(Float::ROUND_UPWARD,a[index(i)]);
    });
    // Like tanh, each activation is one scalar operation per rounding mode.
    // Any sign-specific numerical algorithm belongs inside the implementation.
    report("sigmoid_nearest",[&](unsigned i) {
        return sigmoid(Float::ROUND_TO_NEAREST,a[index(i)]);
    });
    report("sigmoid_down",[&](unsigned i) {
        return sigmoid(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("sigmoid_up",[&](unsigned i) {
        return sigmoid(Float::ROUND_UPWARD,a[index(i)]);
    });
    report("softplus_nearest",[&](unsigned i) {
        return softplus(Float::ROUND_TO_NEAREST,a[index(i)]);
    });
    report("softplus_down",[&](unsigned i) {
        return softplus(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("softplus_up",[&](unsigned i) {
        return softplus(Float::ROUND_UPWARD,a[index(i)]);
    });
    report("atan_nearest",[& ](unsigned i) {
        return atan(Float::ROUND_TO_NEAREST, a[index(i)]);
    }
    );
    report("atan_down",[&](unsigned i) {
        return atan(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("atan_up",[&](unsigned i) {
        return atan(Float::ROUND_UPWARD,a[index(i)]);
    });
    report("asin_nearest", [&](unsigned i) {
        return asin(Float::ROUND_TO_NEAREST, a[index(i)]);
    });
    report("asin_down",[&](unsigned i) {
        return asin(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("asin_up",[&](unsigned i) {
        return asin(Float::ROUND_UPWARD,a[index(i)]);
    });
    report("acos_nearest", [&](unsigned i) {
        return acos(Float::ROUND_TO_NEAREST, a[index(i)]);
    });
    report("acos_down",[&](unsigned i) {
        return acos(Float::ROUND_DOWNWARD,a[index(i)]);
    });
    report("acos_up",[&](unsigned i) {
        return acos(Float::ROUND_UPWARD,a[index(i)]);
    });
    report("pi_nearest", [&](unsigned) {
        return Float::pi(Float::ROUND_TO_NEAREST, precision);
    });
    report("pi_down",[&](unsigned) {
        return Float::pi(Float::ROUND_DOWNWARD,precision);
    });
    report("pi_up",[&](unsigned) {
        return Float::pi(Float::ROUND_UPWARD,precision);
    });
    report("pow7_nearest",[& ](unsigned i) {
        return pow(Float::ROUND_TO_NEAREST, a[index(i)], Nat(7u));
    }
    );
    report("pow7_down",[&](unsigned i) {
        return pow(Float::ROUND_DOWNWARD,a[index(i)],Nat(7u));
    });
    report("pow7_up",[&](unsigned i) {
        return pow(Float::ROUND_UPWARD,a[index(i)],Nat(7u));
    });
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
                 unsigned samples, Ariadne::Benchmarking::Runner const& runner,
                 std::vector<BenchmarkDatum>& data) {
    const BackendPool<Float> pool(bits);
    auto report = [&](char const* name, auto&& operation) {
        auto stats = measure_operation(runner,name,operation,iterations,samples);
        data.push_back({backend,name,bits,stats.median,stats});
    };
    for_each_backend_operation(pool,report);
}


} // namespace AriadneBenchmark
#endif // ARIADNE_NUMERIC_BENCHMARK_BACKEND_WORKLOAD_HPP
