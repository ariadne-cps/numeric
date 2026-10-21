/***************************************************************************
 * benchmark_compare_dyadic_backends.cpp
 * Paired native NewDyadic / legacy Dyadic comparison (not GMP dyadic).
 ****************************************************************************/
#include "numeric/new_dyadic.hpp"
#include "numeric/dyadic.hpp"
#include "benchmark_command_line.hpp"
#include "logging/logging.hpp"
#include "logging/command_line.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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
    bool details=false;
    double max_slowdown=std::numeric_limits<double>::infinity();
};

// Reuse the same Ariadne CLI and validation as the Integer/Rational
// backends; do not implement another hand-written argv loop.
List<ArgumentParser> comparison_parsers(Selection& s) {
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
    parsers.append(make_unvalued_argument_parser(
        "","details","Print each matched NewDyadic/Dyadic ratio",
        [&s] { s.details=true; }));
    parsers.append(make_valued_argument_parser(
        "","max-slowdown","Fail if worst NewDyadic/Dyadic slowdown exceeds this",
        [&s](String const& v)->VoidFunction {
            const std::string text(v);
            std::size_t used=0u;
            const double maximum=std::stod(text,&used);
            if(used!=text.size() || !std::isfinite(maximum) || maximum<=0.0)
                throw std::invalid_argument("invalid positive slowdown ratio");
            return [&s,maximum] { s.max_slowdown=maximum; };
        }));
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
    std::vector<Dyadic> old_a,old_b;
    std::vector<NewDyadic> now_a,now_b;
    unsigned exp;
    explicit Pool(unsigned bits):exp(bits/5u+1u) {
        const auto high=shift_left(BigUnsigned(1u),bits-1u);
        for(unsigned i=0;i<32u;++i) {
            first.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(13u+2u*i)),i%4u==0u || i%4u==3u));
            second.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(31u+5u*i)),i%4u==1u || i%4u==3u));
            old_a.emplace_back(first.back(),exp);
            old_b.emplace_back(second.back(),exp+1u);
            now_a.emplace_back(first.back(),exp);
            now_b.emplace_back(second.back(),exp+1u);
        }
        for(unsigned i=0;i<32u;++i) {
            if(now_a[i].to_rational()!=Rational(old_a[i]) ||
               now_b[i].to_rational()!=Rational(old_b[i]) ||
               (now_a[i]+now_b[i]).to_rational()!=Rational(old_a[i]+old_b[i]) ||
               (now_a[i]-now_b[i]).to_rational()!=Rational(old_a[i]-old_b[i]) ||
               (now_a[i]*now_b[i]).to_rational()!=Rational(old_a[i]*old_b[i]) ||
               pow(now_a[i],Nat(3u)).to_rational()!=Rational(pow(old_a[i],Nat(3u))) ||
               floor(now_a[i])!=floor(old_a[i]) ||
               ceil(now_a[i])!=ceil(old_a[i]) ||
               round(now_a[i])!=round(old_a[i]) ||
               (now_a[i]<now_b[i])!=(old_a[i]<old_b[i])) {
                throw std::runtime_error("NewDyadic / Dyadic exact verification mismatch");
            }
        }
    }
};
struct Datum {
    std::string operation;
    unsigned bits;
    double reference_ns;
    double candidate_ns;
};

void run(unsigned bits,unsigned iterations,unsigned samples,
         std::vector<Datum>& data) {
    const Pool x(bits);
    constexpr unsigned mask=31u;
    auto pair=[&](char const* name,auto&& old_fn,auto&& new_fn) {
        const double old_ns=measure(old_fn,iterations,samples);
        const double new_ns=measure(new_fn,iterations,samples);
        if(!std::isfinite(old_ns) || !std::isfinite(new_ns) ||
                old_ns<=0.0 || new_ns<=0.0) {
            throw std::runtime_error(std::string("invalid benchmark timing: ")+name);
        }
        data.push_back({name,bits,old_ns,new_ns});
    };
    pair("construct_integer",[&](unsigned i){return Dyadic(x.first[i&mask]);},
         [&](unsigned i){return NewDyadic(x.first[i&mask]);});
    pair("construct_fraction",[&](unsigned i){return Dyadic(x.first[i&mask],x.exp);},
         [&](unsigned i){return NewDyadic(x.first[i&mask],x.exp);});
    pair("copy",[&](unsigned i){return Dyadic(x.old_a[i&mask]);},
         [&](unsigned i){return NewDyadic(x.now_a[i&mask]);});
    pair("add",[&](unsigned i){return x.old_a[i&mask]+x.old_b[i&mask];},
         [&](unsigned i){return x.now_a[i&mask]+x.now_b[i&mask];});
    pair("sub",[&](unsigned i){return x.old_a[i&mask]-x.old_b[i&mask];},
         [&](unsigned i){return x.now_a[i&mask]-x.now_b[i&mask];});
    pair("mul",[&](unsigned i){return x.old_a[i&mask]*x.old_b[i&mask];},
         [&](unsigned i){return x.now_a[i&mask]*x.now_b[i&mask];});
    pair("square",[&](unsigned i){return sqr(x.old_a[i&mask]);},
         [&](unsigned i){return sqr(x.now_a[i&mask]);});
    pair("negate",[&](unsigned i){return -x.old_a[i&mask];},
         [&](unsigned i){return -x.now_a[i&mask];});
    pair("abs",[&](unsigned i){return abs(x.old_a[i&mask]);},
         [&](unsigned i){return abs(x.now_a[i&mask]);});
    pair("half",[&](unsigned i){return hlf(x.old_a[i&mask]);},
         [&](unsigned i){return hlf(x.now_a[i&mask]);});
    pair("compare",[&](unsigned i){return cmp(x.old_a[i&mask],x.old_b[i&mask]);},
         [&](unsigned i){return cmp(x.now_a[i&mask],x.now_b[i&mask]);});
    pair("min",[&](unsigned i){return min(x.old_a[i&mask],x.old_b[i&mask]);},
         [&](unsigned i){return min(x.now_a[i&mask],x.now_b[i&mask]);});
    pair("max",[&](unsigned i){return max(x.old_a[i&mask],x.old_b[i&mask]);},
         [&](unsigned i){return max(x.now_a[i&mask],x.now_b[i&mask]);});
    pair("floor",[&](unsigned i){return floor(x.old_a[i&mask]);},
         [&](unsigned i){return floor(x.now_a[i&mask]);});
    pair("ceil",[&](unsigned i){return ceil(x.old_a[i&mask]);},
         [&](unsigned i){return ceil(x.now_a[i&mask]);});
    pair("round",[&](unsigned i){return round(x.old_a[i&mask]);},
         [&](unsigned i){return round(x.now_a[i&mask]);});
    pair("to_double",[&](unsigned i){return x.old_a[i&mask].get_d();},
         [&](unsigned i){return x.now_a[i&mask].get_d();});
    pair("pow3",[&](unsigned i){return pow(x.old_a[i&mask],Nat(3u));},
         [&](unsigned i){return pow(x.now_a[i&mask],Nat(3u));});
}
struct Aggregate {
    std::size_t count=0u;
    double logarithms=0.0;
    double reference_sum=0.0;
    double candidate_sum=0.0;
    double worst=0.0;
    std::string worst_name;

    void add(Datum const& d) {
        const double ratio=d.candidate_ns/d.reference_ns;
        ++count;
        logarithms+=std::log(ratio);
        reference_sum+=d.reference_ns;
        candidate_sum+=d.candidate_ns;
        if(ratio>worst) {
            worst=ratio;
            worst_name=d.operation+", "+std::to_string(d.bits)+" bits";
        }
    }
    void print(std::string const& label) const {
        if(count==0u || reference_sum<=0.0)
            throw std::runtime_error("empty Dyadic benchmark aggregate");
        LOGGING_PRINTLN(std::left<<std::setw(12)<<label
            <<"  ops="<<std::setw(3)<<count
            <<"  geomean="<<std::fixed<<std::setprecision(3)
            <<std::exp(logarithms/static_cast<double>(count))<<"x"
            <<"  total_time="<<candidate_sum/reference_sum<<"x"
            <<"  worst="<<worst<<"x ("<<worst_name<<")")
    }
};

void print_measurements(std::vector<Datum> const& data,Selection const& s) {
    if(s.summary_only) return;
    LOGGING_PRINTLN("Dyadic / NewDyadic benchmark | iterations="
                    <<s.iterations<<" | samples="<<s.samples)
    LOGGING_PRINTLN(std::string(79u,'='))
    LOGGING_PRINTLN(std::left<<std::setw(24)<<"Operation"
        <<std::setw(8)<<"Bits"<<std::right
        <<std::setw(16)<<"Dyadic (ns)"<<std::setw(16)<<"NewDyadic (ns)")
    LOGGING_PRINTLN(std::string(79u,'-'))
    for(auto const& d:data) {
        LOGGING_PRINTLN(std::left<<std::setw(24)<<d.operation
            <<std::setw(8)<<d.bits<<std::right
            <<std::fixed<<std::setprecision(3)
            <<std::setw(16)<<d.reference_ns<<std::setw(16)<<d.candidate_ns)
    }
    LOGGING_PRINTLN(std::string(79u,'='))
}

void print_comparison(std::vector<Datum> const& data,Selection const& s) {
    std::map<unsigned,Aggregate> widths;
    Aggregate overall;
    if(s.details) {
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("Detailed comparisons (NewDyadic / Dyadic)")
        LOGGING_PRINTLN(std::string(79u,'-'))
        LOGGING_PRINTLN(std::left<<std::setw(24)<<"Operation"
            <<std::setw(8)<<"Bits"<<std::setw(16)<<"Dyadic (ns)"
            <<std::setw(16)<<"NewDyadic (ns)"<<"Slowdown")
    }
    for(auto const& d:data) {
        widths[d.bits].add(d);
        overall.add(d);
        if(s.details) {
            LOGGING_PRINTLN(std::left<<std::setw(24)<<d.operation
                <<std::setw(8)<<d.bits
                <<std::setw(16)<<std::fixed<<std::setprecision(3)
                <<d.reference_ns<<std::setw(16)<<d.candidate_ns
                <<d.candidate_ns/d.reference_ns<<"x")
        }
    }
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN(std::string(79u,'='))
    LOGGING_PRINTLN("NewDyadic / Dyadic slowdown (x; >1 means slower)")
    LOGGING_PRINTLN(std::string(79u,'-'))
    for(auto const& [bits,aggregate]:widths) {
        aggregate.print(std::to_string(bits)+" bits");
    }
    overall.print("ALL");
    LOGGING_PRINTLN(std::string(79u,'='))
    LOGGING_PRINTLN("geomean weights each (operation, size) equally; "
                    "total_time weights by measured cost.")
    if(std::isfinite(s.max_slowdown) && overall.worst>s.max_slowdown) {
        throw std::runtime_error(
            "NewDyadic slowdown exceeds --max-slowdown: "+
            std::to_string(overall.worst)+"x > "+
            std::to_string(s.max_slowdown)+"x");
    }
}

}
int main(int argc,char const* argv[]) {
    BenchmarkCommandLine::initialize_logging();
    try {
        Selection selection;
        CommandLineInterface cli;
        cli.add(logging_command_line());
        cli.add(CommandLineContribution(
            "numeric-dyadic-comparison",comparison_parsers(selection)));
        const auto args=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--max-slowdown","--scheduler","-s","--theme","-t",
             "--verbosity","-v"});
        if(!cli.acquire(args)) {
            return BenchmarkCommandLine::requests_help(args)?0:1;
        }
        LOGGING_PRINTLN("Dyadic vs NewDyadic | iterations="
            <<selection.iterations<<" | samples="<<selection.samples
            <<" | reference=Dyadic (mpf), not GMP")
        std::vector<Datum> data;
        data.reserve(selection.sizes.size()*18u);
        for(unsigned bits:selection.sizes) {
            LOGGING_PRINTLN("Measuring "<<bits<<" bits ...")
            run(bits,selection.iterations,selection.samples,data);
            LOGGING_PRINTLN("Completed "<<bits<<" bits")
        }
        print_measurements(data,selection);
        // Pool checks exact invariants outside the measured samples.
        LOGGING_PRINTLN("Numerical verification PASS: "
            <<selection.sizes.size()*64u<<" paired inputs; "
            <<"signed add/sub/mul, cube, floor/ceil/round and ordering.")
        print_comparison(data,selection);
        return 0;
    } catch(std::exception const& error) {
        LOGGING_PRINTLN("Dyadic comparison benchmark: "<<error.what())
        return 1;
    }
}
