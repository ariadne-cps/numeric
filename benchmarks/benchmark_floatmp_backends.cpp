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

template<class Float>
void run_backend(char const* backend, unsigned bits, unsigned iterations, unsigned samples) {
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
        std::cout<<backend<<','<<name<<','<<bits<<','<<iterations<<','<<samples<<','
                 <<std::fixed<<std::setprecision(3)
                 <<median_ns(operation,iterations,samples)<<'\n';
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

int main(int argc,char** argv) {
    try {
        unsigned iterations=argc>=2 ? static_cast<unsigned>(std::stoul(argv[1])) : 300;
        unsigned samples=argc>=3 ? static_cast<unsigned>(std::stoul(argv[2])) : 7;
        if (argc>3 || iterations==0 || samples<3 || samples%2==0) {
            std::cerr<<"usage: benchmark_floatmp_backends [iterations>0] [odd_samples>=3]\n";
            return 2;
        }
        std::cout<<"backend,operation,bits,iterations,samples,median_ns\n";
        for (unsigned bits : {53u,128u,256u,1024u}) {
            run_backend<FloatMP>("mpfr",bits,iterations,samples);
#ifdef ARIADNE_NUMERIC_ENABLE_FLOATVP
            run_backend<FloatVP>("floatvp",bits,iterations,samples);
#endif
        }
    } catch (std::exception const& e) {
        std::cerr<<"benchmark_floatmp_backends: "<<e.what()<<'\n';
        return 1;
    }
}
