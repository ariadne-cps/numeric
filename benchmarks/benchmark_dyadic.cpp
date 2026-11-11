/***************************************************************************
 *            benchmark_dyadic.cpp
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

#include "numeric/dyadic.hpp"
#include "ariadne-benchmarking.hpp"
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
    Ariadne::Benchmarking::keep_alive(value);
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
            const unsigned n=Benchmarking::parse_unsigned<unsigned>(v);
            if(n==0u) throw std::invalid_argument("iterations must be positive");
            return [&s,n] { s.iterations=n; };
        }));
    parsers.append(make_valued_argument_parser(
        "n","samples","Odd number of samples, at least three",
        [&s](String const& v)->VoidFunction {
            const unsigned n=Benchmarking::parse_unsigned<unsigned>(v);
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
                const unsigned n=Benchmarking::parse_unsigned<unsigned>(String(token));
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

template<class F>
Ariadne::Benchmarking::Statistics measure(
        Ariadne::Benchmarking::Runner const& runner,
        char const* name,F const& f,unsigned n,unsigned samples) {
    auto times=runner.run(name,samples,
        [] {},
        [&] {
            for(unsigned i=0;i<n;++i) { auto value=f(i);consume(value); }
        },
        [] {},
        static_cast<double>(n));
    return Ariadne::Benchmarking::summarize(times);
}


struct Pool {
    std::vector<Integer> first,second;
    std::vector<Dyadic> a,b;
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
    Ariadne::Benchmarking::Statistics statistics;
};

void run(unsigned bits,unsigned iterations,unsigned samples,
         Ariadne::Benchmarking::Runner const& runner,
         std::vector<Datum>& data) {
    const Pool x(bits);
    constexpr unsigned mask=31u;
    auto time=[&](char const* name,auto&& operation) {
        auto stats=measure(runner,name,operation,iterations,samples);
        if(!std::isfinite(stats.median) || stats.median<=0.0)
            throw std::runtime_error(std::string("invalid benchmark timing: ")+name);
        data.push_back({name,bits,stats.median,stats});
    };
    time("construct_integer",[&](unsigned i){return Dyadic(x.first[i&mask]);});
    time("construct_fraction",[&](unsigned i){return Dyadic(x.first[i&mask],x.exp);});
    time("copy",[&](unsigned i){return Dyadic(x.a[i&mask]);});
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

void print_times(std::vector<Datum> const& data,Selection const& s,
                 Ariadne::Benchmarking::SummaryFields const& selected) {
    namespace Bench=Ariadne::Benchmarking;
    Bench::TableOptions style=Bench::TableOptions::ruled();
    style.group_columns={1u};
    Bench::Report report("Native Dyadic operation times",selected);
    report.set_dimension_columns({
        {"Operation"},
        {"Bits",0u,Bench::TableAlignment::Right}
    });
    report.set_table_options(style);
    std::map<unsigned,double> totals;
    double grand_total=0.0;
    for(auto const& datum:data) {
        totals[datum.bits]+=datum.nanoseconds;
        grand_total+=datum.nanoseconds;
        if(!s.summary_only)
            report.add_row({datum.operation,std::to_string(datum.bits)},datum.statistics);
    }
    if(!s.summary_only) {
        LOGGING_PRINTLN("Native Dyadic operation times | iterations="
                        <<s.iterations<<" | samples="<<s.samples)
        report.print();
    }
    Bench::Table summary({
        {"Group"},
        {"Bits",0u,Bench::TableAlignment::Right},
        {"Sum of medians (ns)",0u,Bench::TableAlignment::Right}
    },Bench::TableOptions::ruled());
    for(auto const& [bits,ns]:totals)
        summary.add_row({"TOTAL",std::to_string(bits),Bench::Table::number(ns)});
    summary.add_row({"TOTAL ALL","-",Bench::Table::number(grand_total)});
    LOGGING_PRINTLN(summary.str())
}

} // namespace

int main(int argc,char const* argv[]) {
    try {
        Selection selection;
        Benchmarking::CommandLine cli;
        cli.add(CommandLineContribution(
            "numeric-dyadic-benchmark",benchmark_parsers(selection)));
        if(!cli.acquire(argc,argv)) return 0;
        LOGGING_PRINTLN("Native Dyadic | iterations="
            <<selection.iterations<<" | samples="<<selection.samples)
        Benchmarking::Runner runner(cli,selection.sizes.size());
        std::vector<Datum> data;
        data.reserve(selection.sizes.size()*18u);
        for(unsigned bits:selection.sizes) {
            LOGGING_PRINTLN("Measuring "<<bits<<" bits ...")
            run(bits,selection.iterations,selection.samples,runner,data);
            LOGGING_PRINTLN("Completed "<<bits<<" bits")
        }
        print_times(data,selection,cli.summary_fields());
        return 0;
    } catch(std::exception const& error) {
        LOGGING_PRINTLN("Dyadic benchmark: "<<error.what())
        return 1;
    }
}
