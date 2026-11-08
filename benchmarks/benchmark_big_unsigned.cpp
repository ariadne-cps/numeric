/***************************************************************************
 *            benchmark_big_unsigned.cpp
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

/*! \file benchmark_big_unsigned.cpp
 *  \brief Absolute native BigUnsigned operation timings with exact invariants.
 */

#include "numeric/big_unsigned.hpp"
#include "benchmark_absolute_workload.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace Ariadne;
namespace Bench=Ariadne::Benchmarking;
using AriadneBenchmark::AbsoluteDatum;

bool equal(BigUnsigned const& a,BigUnsigned const& b) {
    return compare(a,b)==0;
}

struct Pool {
    static constexpr unsigned count=32u;
    std::vector<BigUnsigned> a,b;
    std::vector<std::string> decimal_a;
    explicit Pool(unsigned bits) {
        auto const base=shift_left(BigUnsigned(1u),bits-1u);
        for(unsigned i=0u;i<count;++i) {
            a.push_back(add(base,BigUnsigned(101u+3u*i)));
            b.push_back(add(base,BigUnsigned(19u+i)));
            decimal_a.push_back(to_decimal(a.back()));
        }
        // These are correctness checks, not separate timed operations. Verify
        // division, subtraction and square root against exact algebraic laws.
        for(unsigned i=0u;i<count;++i) {
            auto const& x=a[i];
            auto const& y=b[i];
            if(compare(x,y)<0 || !equal(add(sub(x,y),y),x) ||
               !equal(sub(add(x,y),y),x) ||
               !equal(mul(x,y),mul(y,x)) ||
               !equal(shift_right(shift_left(x,65u),65u),x) ||
               !equal(from_decimal(decimal_a[i]),x))
                throw std::runtime_error("BigUnsigned arithmetic invariant failed");
            auto const qr=divrem(x,y);
            if(!equal(add(mul(qr.first,y),qr.second),x) ||
               compare(qr.second,y)>=0 ||
               !equal(quot(x,y),qr.first) ||
               !equal(rem(x,y),qr.second))
                throw std::runtime_error("BigUnsigned Euclidean division failed");
            auto const sr=sqrtrem(x);
            if(!equal(add(square(sr.first),sr.second),x) ||
               compare(sr.second,add(sr.first,sr.first))>0 ||
               !equal(isqrt(x),sr.first))
                throw std::runtime_error("BigUnsigned square-root invariant failed");
        }
    }
};

void measure_width(unsigned bits,
                   AriadneBenchmark::AbsoluteSelection const& selection,
                   Bench::Runner const& runner,
                   std::vector<AbsoluteDatum>& data) {
    Pool const pool(bits);
    constexpr unsigned mask=Pool::count-1u;
    auto time=[&](char const* name,auto&& op,unsigned repetitions) {
        AriadneBenchmark::record_absolute(data,runner,name,bits,repetitions,
                                          selection.samples,std::forward<decltype(op)>(op));
    };
    auto normal=[&](char const* name,auto&& op) {
        time(name,std::forward<decltype(op)>(op),selection.iterations);
    };
    // Keep the historical sqrt workload's reduced repetition count for large
    // operands; each reported value is still normalized ns/operation.
    unsigned const sqrt_iterations=std::max(1u,selection.iterations/
                                          std::max(1u,bits/64u));
    normal("add",[&](unsigned i){return add(pool.a[i&mask],pool.b[i&mask]);});
    normal("sub",[&](unsigned i){return sub(pool.a[i&mask],pool.b[i&mask]);});
    normal("compare",[&](unsigned i){return compare(pool.a[i&mask],pool.b[i&mask]);});
    normal("shl_1",[&](unsigned i){return shift_left(pool.a[i&mask],1u);});
    normal("shr_1",[&](unsigned i){return shift_right(pool.a[i&mask],1u);});
    normal("shl_65",[&](unsigned i){return shift_left(pool.a[i&mask],65u);});
    normal("shr_65",[&](unsigned i){return shift_right(pool.a[i&mask],65u);});
    normal("mul",[&](unsigned i){return mul(pool.a[i&mask],pool.b[i&mask]);});
    normal("square",[&](unsigned i){return square(pool.a[i&mask]);});
    normal("quot",[&](unsigned i){return quot(pool.a[i&mask],pool.b[i&mask]);});
    normal("rem",[&](unsigned i){return rem(pool.a[i&mask],pool.b[i&mask]);});
    normal("divrem",[&](unsigned i){return divrem(pool.a[i&mask],pool.b[i&mask]);});
    normal("gcd",[&](unsigned i){return gcd(pool.a[i&mask],pool.b[i&mask]);});
    time("isqrt",[&](unsigned i){return isqrt(pool.a[i&mask]);},sqrt_iterations);
    time("sqrtrem",[&](unsigned i){return sqrtrem(pool.a[i&mask]);},sqrt_iterations);
    normal("pow3",[&](unsigned i){return pow(pool.a[i&mask],3u);});
    normal("to_decimal",[&](unsigned i){return to_decimal(pool.a[i&mask]);});
    normal("from_decimal",[&](unsigned i){return from_decimal(pool.decimal_a[i&mask]);});
}

} // namespace

int main(int argc,char const* argv[]) {
    Ariadne::
    try {
        AriadneBenchmark::AbsoluteSelection selection;
        Bench::CommandLine cli;
        if(!AriadneBenchmark::acquire_absolute(
              selection,cli,"numeric-big-unsigned-benchmark",argc,argv))
            return 0;
        Bench::Runner runner(cli,selection.bits.size());
        std::vector<AbsoluteDatum> rows;
        rows.reserve(selection.bits.size()*18u);
        for(unsigned bits:selection.bits) {
            LOGGING_PRINTLN("Measuring BigUnsigned "<<bits<<" bits ...")
            measure_width(bits,selection,runner,rows);
        }
        for(unsigned bits:selection.bits) {
            unsigned count=0u;
            for(auto const& row:rows) if(row.bits==bits) ++count;
            if(count!=18u) throw std::runtime_error("BigUnsigned operation count changed");
        }
        AriadneBenchmark::report_absolute("BigUnsigned",selection,
                                           cli.summary_fields(),rows);
        return 0;
    } catch(std::exception const& e) {
        LOGGING_PRINTLN("benchmark_big_unsigned: "<<e.what())
        return 1;
    }
}
