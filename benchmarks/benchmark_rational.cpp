/***************************************************************************
 *            benchmark_rational.cpp
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

/** @file benchmark_rational.cpp
 *  @brief Absolute timing of the native Rational arithmetic operations.
 */

#include "numeric/rational.hpp"
#include "benchmark_absolute_workload.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using namespace Ariadne;
namespace Bench=Ariadne::Benchmarking;
using AriadneBenchmark::AbsoluteDatum;

struct Pool {
    static constexpr unsigned count=32u;
    std::vector<Integer> num_a,den_a,num_b,den_b;
    std::vector<Rational> a,b;

    explicit Pool(unsigned bits) {
        auto const high=shift_left(BigUnsigned(1u),bits-1u);
        for(unsigned i=0u;i<count;++i) {
            auto const sign_a=(i%4u==0u || i%4u==3u);
            auto const sign_b=(i%4u==1u || i%4u==3u);
            num_a.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(13u+2u*i)),sign_a));
            den_a.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(17u+3u*i)),false));
            num_b.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(31u+5u*i)),sign_b));
            den_b.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(19u+7u*i)),false));
            a.emplace_back(num_a.back(),den_a.back());
            b.emplace_back(num_b.back(),den_b.back());
        }
        // Correctness checks are deliberately outside the measured region.
        // Compare canonical identities, not a second copy of the algorithm.
        Rational const one(1u);
        for(unsigned i=0u;i<count;++i) {
            auto const& x=a[i];
            auto const& y=b[i];
            if((x+y)-y!=x || (x-y)+y!=x ||
               (x*y)/y!=x || (x/y)*y!=x ||
               x*rec(x)!=one || x*x!=sqr(x) ||
               !(floor(x)<=ceil(x)) ||
               !(Rational(floor(x))<=x) ||
               !(x<=Rational(ceil(x))) ||
               !((x<y)==(y>x)) ||
               !(x==Rational(x.get_num(),x.get_den())))
                throw std::runtime_error("Rational arithmetic invariant failed");
        }
    }
};

void measure_width(unsigned bits,
                   AriadneBenchmark::AbsoluteSelection const& selection,
                   Bench::Runner const& runner,
                   std::vector<AbsoluteDatum>& data) {
    Pool const pool(bits);
    constexpr unsigned mask=Pool::count-1u;
    auto record=[&](char const* name,auto&& operation) {
        AriadneBenchmark::record_absolute(
            data,runner,name,bits,selection.iterations,selection.samples,
            std::forward<decltype(operation)>(operation));
    };

    // Preserve all 15 historical Rational operation names and operand kinds.
    record("construct_fraction",[&](unsigned i) {
        return Rational(pool.num_a[i&mask],pool.den_a[i&mask]);
    });
    record("copy",[&](unsigned i){return Rational(pool.a[i&mask]);});
    record("compare",[&](unsigned i){
        auto const& a=pool.a[i&mask];
        auto const& b=pool.b[i&mask];
        return a<b ? -1 : b<a ? 1 : 0;
    });
    record("negate",[&](unsigned i){return -pool.a[i&mask];});
    record("add",[&](unsigned i){return pool.a[i&mask]+pool.b[i&mask];});
    record("sub",[&](unsigned i){return pool.a[i&mask]-pool.b[i&mask];});
    record("mul",[&](unsigned i){return pool.a[i&mask]*pool.b[i&mask];});
    record("div",[&](unsigned i){return pool.a[i&mask]/pool.b[i&mask];});
    record("square",[&](unsigned i){return sqr(pool.a[i&mask]);});
    record("reciprocal",[&](unsigned i){return rec(pool.a[i&mask]);});
    record("numerator",[&](unsigned i){return pool.a[i&mask].get_num();});
    record("denominator",[&](unsigned i){return pool.a[i&mask].get_den();});
    record("floor",[&](unsigned i){return floor(pool.a[i&mask]);});
    record("ceil",[&](unsigned i){return ceil(pool.a[i&mask]);});
    record("to_double",[&](unsigned i){return pool.a[i&mask].get_d();});
}

} // namespace

int main(int argc,char const* argv[]) {
    try {
        AriadneBenchmark::AbsoluteSelection selection;
        Bench::CommandLine cli;
        if(!AriadneBenchmark::acquire_absolute(
              selection,cli,"numeric-rational-benchmark",argc,argv))
            return 0;
        Bench::Runner runner(cli,selection.bits.size());
        std::vector<AbsoluteDatum> rows;
        rows.reserve(selection.bits.size()*15u);
        for(unsigned bits:selection.bits) {
            LOGGING_PRINTLN("Measuring Rational "<<bits<<" bits ...")
            measure_width(bits,selection,runner,rows);
        }
        for(unsigned bits:selection.bits) {
            unsigned count=0u;
            for(auto const& row:rows) if(row.bits==bits) ++count;
            if(count!=15u) throw std::runtime_error("Rational operation count changed");
        }
        AriadneBenchmark::report_absolute("Rational",selection,
                                           cli.summary_fields(),rows);
        return 0;
    } catch(std::exception const& e) {
        LOGGING_PRINTLN("benchmark_rational: "<<e.what())
        return 1;
    }
}
