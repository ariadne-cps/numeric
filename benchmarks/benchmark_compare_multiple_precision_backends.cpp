// Part of Ariadne Numeric. GPL-3.0-or-later.
// Same monomorphic test code is instantiated for each implementation: no virtual dispatch.
#include "numeric/floatmp.hpp"
#include "numeric/rational.hpp"
#ifdef ARIADNE_NUMERIC_ENABLE_FLOATVP
#include "numeric/floatvp.hpp"
#endif
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <limits>
#include <map>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Ariadne;

template<class T> inline void keep_alive(T const& value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&value) : "memory");
#else
    static void const* volatile sink = nullptr;
    sink = &value;
#endif
}

template<class Operation>
double median_ns(Operation&& operation, unsigned iterations, unsigned samples) {
    std::vector<double> measurements;
    measurements.reserve(samples);
    for (unsigned sample=0; sample<=samples; ++sample) {
        auto start=std::chrono::steady_clock::now();
        for (unsigned i=0; i<iterations; ++i) {
            auto result=operation(i);
            keep_alive(result);
        }
        auto stop=std::chrono::steady_clock::now();
        if (sample!=0) { // unmeasured warm-up
            measurements.push_back(
                std::chrono::duration<double,std::nano>(stop-start).count()/iterations);
        }
    }
    std::sort(measurements.begin(),measurements.end());
    return measurements[samples/2];
}

struct BenchmarkDatum {
    std::string backend;
    std::string operation;
    unsigned bits;
    double nanoseconds;
};

template<class Float>
void run_backend(char const* backend, unsigned bits, unsigned iterations,
                 unsigned samples, bool print_csv, std::vector<BenchmarkDatum>& data) {
    using Precision=typename Float::PrecisionType;
    Precision precision = bits==53 ? Precision(53_bits)
                        : bits==128 ? Precision(128_bits)
                        : bits==256 ? Precision(256_bits) : Precision(1024_bits);
    constexpr unsigned pool_size=32;
    std::vector<Rational> exact_a,exact_b;
    std::vector<Float> a,b;
    exact_a.reserve(pool_size); exact_b.reserve(pool_size);
    a.reserve(pool_size); b.reserve(pool_size);

    // Nontrivial significands: avoid benchmarking only exact small integers.
    const Integer base=pow(Integer(2),Nat(bits));
    for (unsigned i=0; i<pool_size; ++i) {
        exact_a.emplace_back(base+Integer(13u+2u*i),2*base+Integer(17u+3u*i));
        exact_b.emplace_back(base+Integer(31u+5u*i),3*base+Integer(19u+7u*i));
        a.emplace_back(exact_a.back(),Float::ROUND_TO_NEAREST,precision);
        b.emplace_back(exact_b.back(),Float::ROUND_TO_NEAREST,precision);
    }

    auto report=[&](char const* name,auto&& operation) {
        const double ns=median_ns(operation,iterations,samples);
        data.push_back({backend,name,bits,ns});
        if (print_csv) {
            std::cout<<backend<<','<<name<<','<<bits<<','<<iterations<<','<<samples<<','
                     <<std::fixed<<std::setprecision(6)<<ns<<'\n';
        }
    };
    auto index=[](unsigned i){return i%pool_size;};
    report("construct_rational",[&](unsigned i) {
        return Float(exact_a[index(i)],Float::ROUND_TO_NEAREST,precision);
    });
    report("copy",[&](unsigned i){return Float(a[index(i)]);});
    report("compare",[&](unsigned i){return cmp(a[index(i)],b[index(i)]);});
    report("add_nearest",[&](unsigned i){return add(Float::ROUND_TO_NEAREST,a[index(i)],b[index(i)]);});
    report("add_down",[&](unsigned i){return add(Float::ROUND_DOWNWARD,a[index(i)],b[index(i)]);});
    report("add_up",[&](unsigned i){return add(Float::ROUND_UPWARD,a[index(i)],b[index(i)]);});
    report("sub_nearest",[&](unsigned i){return sub(Float::ROUND_TO_NEAREST,a[index(i)],b[index(i)]);});
    report("mul_nearest",[&](unsigned i){return mul(Float::ROUND_TO_NEAREST,a[index(i)],b[index(i)]);});
    report("mul_down",[&](unsigned i){return mul(Float::ROUND_DOWNWARD,a[index(i)],b[index(i)]);});
    report("mul_up",[&](unsigned i){return mul(Float::ROUND_UPWARD,a[index(i)],b[index(i)]);});
    report("div_nearest",[&](unsigned i){return div(Float::ROUND_TO_NEAREST,a[index(i)],b[index(i)]);});
    report("div_down",[&](unsigned i){return div(Float::ROUND_DOWNWARD,a[index(i)],b[index(i)]);});
    report("div_up",[&](unsigned i){return div(Float::ROUND_UPWARD,a[index(i)],b[index(i)]);});
    report("fma_nearest",[&](unsigned i){return fma(Float::ROUND_TO_NEAREST,a[index(i)],b[index(i)],a[index(i)]);});
    report("sqrt_nearest",[&](unsigned i){return sqrt(Float::ROUND_TO_NEAREST,a[index(i)]);});
    report("exp_nearest",[&](unsigned i){return exp(Float::ROUND_TO_NEAREST,a[index(i)]);});
    report("log_nearest",[&](unsigned i){return log(Float::ROUND_TO_NEAREST,a[index(i)]);});
    report("sin_nearest",[&](unsigned i){return sin(Float::ROUND_TO_NEAREST,a[index(i)]);});
    report("cos_nearest",[&](unsigned i){return cos(Float::ROUND_TO_NEAREST,a[index(i)]);});
    report("tan_nearest",[&](unsigned i){return tan(Float::ROUND_TO_NEAREST,a[index(i)]);});
    report("atan_nearest",[&](unsigned i){return atan(Float::ROUND_TO_NEAREST,a[index(i)]);});
    report("pow7_nearest",[&](unsigned i){return pow(Float::ROUND_TO_NEAREST,a[index(i)],Nat(7u));});
    report("to_double",[&](unsigned i){return a[index(i)].get_d();});
    report("to_rational",[&](unsigned i){return static_cast<Rational>(a[index(i)]);});
}

#ifdef ARIADNE_NUMERIC_ENABLE_FLOATVP
void print_slowdown_summary(std::vector<BenchmarkDatum> const& data,
                            bool details, double max_slowdown) {
    struct Aggregate {
        std::size_t count=0;
        double logarithmic_sum=0.0;
        double mpfr_sum=0.0;
        double floatvp_sum=0.0;
        double worst_ratio=0.0;
        std::string worst_operation;
        unsigned worst_bits=0;

        void add(std::string const& name, unsigned bits, double reference, double candidate) {
            if (!(std::isfinite(reference) && std::isfinite(candidate))
                || reference<=0.0 || candidate<=0.0) {
                throw std::runtime_error("nonpositive/nonfinite benchmark timing");
            }
            const double ratio=candidate/reference;
            ++count;
            logarithmic_sum+=std::log(ratio);
            mpfr_sum+=reference;
            floatvp_sum+=candidate;
            if (ratio>worst_ratio) {
                worst_ratio=ratio;
                worst_operation=name;
                worst_bits=bits;
            }
        }
        void print(std::string const& name) const {
            if (count==0 || mpfr_sum<=0.0) {
                throw std::runtime_error("cannot summarize empty benchmark results");
            }
            std::cerr<<std::left<<std::setw(12)<<name
                     <<" ops="<<std::setw(3)<<count
                     <<" geomean="<<std::fixed<<std::setprecision(3)
                     <<std::exp(logarithmic_sum/static_cast<double>(count))<<"x"
                     <<" total_time="<<floatvp_sum/mpfr_sum<<"x"
                     <<" worst="<<worst_ratio<<"x ("<<worst_operation
                     <<", "<<worst_bits<<" bits)\n";
        }
    };

    using Key=std::pair<unsigned,std::string>;
    std::map<Key,double> reference, candidate;
    for (auto const& datum : data) {
        auto& collection=datum.backend=="mpfr" ? reference : candidate;
        if (!collection.emplace(Key{datum.bits,datum.operation},datum.nanoseconds).second) {
            throw std::runtime_error("duplicate benchmark measurement");
        }
    }
    if (reference.size()!=candidate.size() || reference.empty()) {
        throw std::runtime_error("FloatVP and FloatMP benchmark rows are not paired");
    }
    std::map<unsigned,Aggregate> per_precision;
    Aggregate overall;
    if (details) {
        std::cerr<<"\nFloatVP / FloatMP detailed comparisons:\n"
                 <<"operation,bits,mpfr_ns,floatvp_ns,slowdown_floatvp_over_mpfr\n";
    }
    for (auto const& [key,reference_ns] : reference) {
        auto it=candidate.find(key);
        if (it==candidate.end()) {
            throw std::runtime_error("FloatVP benchmark operation missing: "+key.second);
        }
        auto& precision_summary=per_precision[key.first];
        precision_summary.add(key.second,key.first,reference_ns,it->second);
        overall.add(key.second,key.first,reference_ns,it->second);
        if (details) {
            std::cerr<<key.second<<','<<key.first<<','
                     <<std::fixed<<std::setprecision(6)<<reference_ns<<','
                     <<it->second<<','<<std::setprecision(3)
                     <<it->second/reference_ns<<'\n';
        }
    }
    std::cerr<<"\nFloatVP / FloatMP slowdown (x; >1 means slower):\n";
    for (auto const& [bits,summary] : per_precision) {
        summary.print(std::to_string(bits)+" bits");
    }
    overall.print("ALL");
    std::cerr<<"geomean weighs each (operation, precision) equally; "
                "total_time weights by measured cost.\n";
    if (std::isfinite(max_slowdown) && overall.worst_ratio>max_slowdown) {
        throw std::runtime_error(
            "FloatVP slowdown exceeds --max-slowdown: "
            +std::to_string(overall.worst_ratio)+"x > "
            +std::to_string(max_slowdown)+"x");
    }
}
#endif

int main(int argc,char** argv) {
    try {
        unsigned iterations=300;
        unsigned samples=7;
        bool summary_only=false;
        [[maybe_unused]] bool details=false;
        double max_slowdown=std::numeric_limits<double>::infinity();
        unsigned positional=0;

        for (int i=1; i<argc; ++i) {
            const std::string arg(argv[i]);
            if (arg=="--summary-only") {
                summary_only=true;
            } else if (arg=="--details") {
                details=true;
            } else if (arg=="--max-slowdown") {
                if (++i>=argc) {
                    throw std::invalid_argument("--max-slowdown requires a positive ratio");
                }
                const std::string value(argv[i]);
                std::size_t read=0;
                max_slowdown=std::stod(value,&read);
                if (read!=value.size() || !std::isfinite(max_slowdown)
                    || max_slowdown<=0.0) {
                    throw std::invalid_argument("--max-slowdown must be positive and finite");
                }
            } else if (arg=="--help") {
                std::cout<<"usage: benchmark_compare_multiple_precision_backends "
                            "[iterations] [odd_samples] "
                            "[--summary-only] [--details] [--max-slowdown RATIO]\n";
                return 0;
            } else if (arg.starts_with("--")) {
                throw std::invalid_argument("unknown option: "+arg);
            } else {
                std::size_t read=0;
                unsigned long parsed=std::stoul(arg,&read);
                if (read!=arg.size() || parsed>std::numeric_limits<unsigned>::max()) {
                    throw std::invalid_argument("invalid iterations/sample count: "+arg);
                }
                if (positional==0) {
                    iterations=static_cast<unsigned>(parsed);
                } else if (positional==1) {
                    samples=static_cast<unsigned>(parsed);
                } else {
                    throw std::invalid_argument("too many positional arguments");
                }
                ++positional;
            }
        }
        if (iterations==0 || samples<3 || samples%2==0) {
            throw std::invalid_argument("iterations must be positive; samples must be odd and >=3");
        }
#ifndef ARIADNE_NUMERIC_ENABLE_FLOATVP
        if (std::isfinite(max_slowdown)) {
            throw std::invalid_argument("--max-slowdown requires BUILD_FLOATVP_COMPARISON=ON");
        }
#endif

        std::vector<BenchmarkDatum> measurements;
        if (!summary_only) {
            std::cout<<"backend,operation,bits,iterations,samples,median_ns\n";
        }
        for (unsigned bits : {53u,128u,256u,1024u}) {
            run_backend<FloatMP>("mpfr",bits,iterations,samples,!summary_only,measurements);
#ifdef ARIADNE_NUMERIC_ENABLE_FLOATVP
            run_backend<FloatVP>("floatvp",bits,iterations,samples,!summary_only,measurements);
#endif
        }
#ifdef ARIADNE_NUMERIC_ENABLE_FLOATVP
        print_slowdown_summary(measurements,details,max_slowdown);
#else
        std::cerr<<"FloatVP slowdown unavailable: BUILD_FLOATVP_COMPARISON is OFF. "
                    "Only FloatMP reference timings were recorded.\n";
#endif
        return 0;
    } catch (std::exception const& e) {
        std::cerr<<"benchmark_compare_multiple_precision_backends: "<<e.what()<<'\n';
        return 1;
    }
}
