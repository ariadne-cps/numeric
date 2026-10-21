/***************************************************************************
 * benchmark_dyadic.cpp
 * Absolute timing of the native dyadic arithmetic backend.
 ****************************************************************************/
#include "numeric/new_dyadic.hpp"
#include "benchmark_command_line.hpp"
#include "logging/logging.hpp"
#include "logging/command_line.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace Ariadne;

template<class T> void consume(T const& value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&value) : "memory");
#else
    static void const* volatile sink=nullptr;
    sink=&value;
#endif
}

struct Selection {
    unsigned iterations=20000u;
    unsigned samples=31u;
    std::vector<unsigned> sizes{53u,128u,256u,1024u,2048u,4096u};
    bool summary_only=false;
};

List<ArgumentParser> benchmark_parsers(Selection& s) {
    List<ArgumentParser> parsers;
    parsers.append(make_valued_argument_parser(
        "i","iterations","Number of operations in each measured sample (positive)",
        [&s](String const& v)->VoidFunction {
            const unsigned n=BenchmarkCommandLine::parse_unsigned<unsigned>(v);
            if(n==0u) throw std::invalid_argument("iterations must be positive");
            return [&s,n] { s.iterations=n; };
        }));
    parsers.append(make_valued_argument_parser(
        "n","samples","Odd number of samples, at least three",
        [&s](String const& v)->VoidFunction {
            const unsigned n=BenchmarkCommandLine::parse_unsigned<unsigned>(v);
            if(n<3u || n%2u==0u) throw std::invalid_argument("invalid sample count");
            return [&s,n] { s.samples=n; };
        }));
    parsers.append(make_valued_argument_parser(
        "b","bits","Comma-separated operand sizes in bits (positive)",
        [&s](String const& v)->VoidFunction {
            const std::string text(v);
            if(text.empty() || text.back()==',')
                throw std::invalid_argument("invalid size list");
            std::vector<unsigned> sizes;
            std::istringstream input(text);
            std::string token;
            while(std::getline(input,token,',')) {
                const unsigned n=BenchmarkCommandLine::parse_unsigned<unsigned>(String(token));
                if(n==0u || std::find(sizes.begin(),sizes.end(),n)!=sizes.end())
                    throw std::invalid_argument("invalid or duplicate operand size");
                sizes.push_back(n);
            }
            if(sizes.empty()) throw std::invalid_argument("empty size list");
            return [&s,sizes] { s.sizes=sizes; };
        }));
    parsers.append(make_unvalued_argument_parser(
        "","summary-only","Print only aggregate statistics",
        [&s] { s.summary_only=true; }));
    return parsers;
}

template<class F> double measure(F const& f,unsigned n,unsigned samples) {
    using Clock=std::chrono::steady_clock;
    std::vector<double> data;data.reserve(samples);
    for(unsigned sample=0;sample<=samples;++sample) {
        auto start=Clock::now();
        for(unsigned i=0;i<n;++i) { auto value=f(i);consume(value); }
        auto stop=Clock::now();
        if(sample) data.push_back(std::chrono::duration<double,std::nano>(stop-start).count()/n);
    }
    std::sort(data.begin(),data.end());
    return data[samples/2u];
}


struct Pool {
    std::vector<Integer> first,second;
    std::vector<NewDyadic> a,b;
    unsigned exp;
    explicit Pool(unsigned bits):exp(bits/5u+1u) {
        const auto high=shift_left(BigUnsigned(1u),bits-1u);
        for(unsigned i=0;i<32u;++i) {
            first.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(13u+2u*i)),i%4u==0u || i%4u==3u));
            second.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(31u+5u*i)),i%4u==1u || i%4u==3u));
            a.emplace_back(first.back(),exp);
            b.emplace_back(second.back(),exp+1u);
        }
        // Verify exact arithmetic against native Rational outside timing.
        const Integer d1=Integer::from_magnitude(
            shift_left(BigUnsigned(1u),exp),false);
        const Integer d2=Integer::from_magnitude(
            shift_left(BigUnsigned(1u),exp+1u),false);
        for(unsigned i=0u;i<32u;++i) {
            const Rational ra(first[i],d1),rb(second[i],d2);
            if(a[i].to_rational()!=ra || b[i].to_rational()!=rb ||
               (a[i]+b[i]).to_rational()!=ra+rb ||
               (a[i]-b[i]).to_rational()!=ra-rb ||
               (a[i]*b[i]).to_rational()!=ra*rb ||
               pow(a[i],Nat(3u)).to_rational()!=ra*ra*ra ||
               static_cast<bool>(a[i]<b[i])!=static_cast<bool>(ra<rb)) {
                throw std::runtime_error("native Dyadic exact verification mismatch");
            }
        }
    }
};

struct Datum {
    std::string operation;
    unsigned bits;
    double nanoseconds;
};

void run(unsigned bits,unsigned iterations,unsigned samples,
         std::vector<Datum>& data) {
    const Pool x(bits);
    constexpr unsigned mask=31u;
    auto time=[&](char const* name,auto&& operation) {
        const double ns=measure(operation,iterations,samples);
        if(!std::isfinite(ns) || ns<=0.0)
            throw std::runtime_error(std::string("invalid benchmark timing: ")+name);
        data.push_back({name,bits,ns});
    };
    time("construct_integer",[&](unsigned i){return NewDyadic(x.first[i&mask]);});
    time("construct_fraction",[&](unsigned i){return NewDyadic(x.first[i&mask],x.exp);});
    time("copy",[&](unsigned i){return NewDyadic(x.a[i&mask]);});
    time("add",[&](unsigned i){return x.a[i&mask]+x.b[i&mask];});
    time("sub",[&](unsigned i){return x.a[i&mask]-x.b[i&mask];});
    time("mul",[&](unsigned i){return x.a[i&mask]*x.b[i&mask];});
    time("square",[&](unsigned i){return sqr(x.a[i&mask]);});
    time("negate",[&](unsigned i){return -x.a[i&mask];});
    time("abs",[&](unsigned i){return abs(x.a[i&mask]);});
    time("half",[&](unsigned i){return hlf(x.a[i&mask]);});
    time("compare",[&](unsigned i){return cmp(x.a[i&mask],x.b[i&mask]);});
    time("min",[&](unsigned i){return min(x.a[i&mask],x.b[i&mask]);});
    time("max",[&](unsigned i){return max(x.a[i&mask],x.b[i&mask]);});
    time("floor",[&](unsigned i){return floor(x.a[i&mask]);});
    time("ceil",[&](unsigned i){return ceil(x.a[i&mask]);});
    time("round",[&](unsigned i){return round(x.a[i&mask]);});
    time("to_double",[&](unsigned i){return x.a[i&mask].get_d();});
    time("pow3",[&](unsigned i){return pow(x.a[i&mask],Nat(3u));});
}

void print_times(std::vector<Datum> const& data,Selection const& s) {
    std::map<unsigned,double> total;
    double grand_total=0.0;
    if(!s.summary_only) {
        LOGGING_PRINTLN("Native Dyadic operation times | iterations="
                        <<s.iterations<<" | samples="<<s.samples)
        LOGGING_PRINTLN(std::string(55u,'='))
        LOGGING_PRINTLN(std::left<<std::setw(24)<<"Operation"
            <<std::setw(8)<<"Bits"<<std::right<<std::setw(16)<<"Time (ns)")
        LOGGING_PRINTLN(std::string(55u,'-'))
    }
    for(auto const& d:data) {
        total[d.bits]+=d.nanoseconds;
        grand_total+=d.nanoseconds;
        if(!s.summary_only) {
            LOGGING_PRINTLN(std::left<<std::setw(24)<<d.operation
                <<std::setw(8)<<d.bits<<std::right<<std::fixed
                <<std::setprecision(3)<<std::setw(16)<<d.nanoseconds)
        }
    }
    if(!s.summary_only) LOGGING_PRINTLN(std::string(55u,'='))
    LOGGING_PRINTLN("Sum of operation medians (ns per one of each operation):")
    for(auto const& [bits,ns]:total) {
        LOGGING_PRINTLN(std::setw(6)<<bits<<" bits: "<<std::fixed
            <<std::setprecision(3)<<ns)
    }
    LOGGING_PRINTLN("All widths, sum (ns): "<<std::fixed
        <<std::setprecision(3)<<grand_total)
}

} // namespace

int main(int argc,char const* argv[]) {
    BenchmarkCommandLine::initialize_logging();
    try {
        Selection selection;
        CommandLineInterface cli;
        cli.add(logging_command_line());
        cli.add(CommandLineContribution(
            "numeric-dyadic-benchmark",benchmark_parsers(selection)));
        const auto args=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if(!cli.acquire(args))
            return BenchmarkCommandLine::requests_help(args)?0:1;
        LOGGING_PRINTLN("Native Dyadic | iterations="
            <<selection.iterations<<" | samples="<<selection.samples)
        std::vector<Datum> data;
        data.reserve(selection.sizes.size()*18u);
        for(unsigned bits:selection.sizes) {
            LOGGING_PRINTLN("Measuring "<<bits<<" bits ...")
            run(bits,selection.iterations,selection.samples,data);
            LOGGING_PRINTLN("Completed "<<bits<<" bits")
        }
        print_times(data,selection);
        return 0;
    } catch(std::exception const& error) {
        LOGGING_PRINTLN("Dyadic benchmark: "<<error.what())
        return 1;
    }
}
