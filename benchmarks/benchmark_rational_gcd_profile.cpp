/***************************************************************************
 *            benchmark_rational_gcd_profile.cpp
 *
 *  Copyright 2026 Ariadne contributors
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

/*! \file benchmark_rational_gcd_profile.cpp
 *  \brief R24 diagnostic for the GCDs used by canonical Rational arithmetic.
 *
 *  Not a replacement for the unchanged official Rational/GMP benchmark.
 *  The binary candidate is local to this executable, never production.
 */

#include "benchmark_exact_arithmetic_workload.hpp"
#include "benchmark_rational_lehmer_candidate.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace Ariadne;
using AriadneBenchmark::RawInteger;
using AriadneBenchmark::RawRational;
using AriadneBenchmark::decimal;
using AriadneBenchmark::exact_arithmetic_median_ns;
using Limb=BigUnsigned::Limb;

constexpr std::size_t binary_capacity=8u;

std::size_t normalized_size(std::array<Limb,binary_capacity> const& data,
                            std::size_t size) {
    while (size!=0u && data[size-1u]==0u) { --size; }
    return size;
}

unsigned trailing_bits(std::array<Limb,binary_capacity> const& data,
                       std::size_t size) {
    std::size_t i=0u;
    while (i<size && data[i]==0u) { ++i; }
    if (i==size) { throw std::logic_error("R24 trailing_bits of zero"); }
    return static_cast<unsigned>(64u*i)+
           static_cast<unsigned>(std::countr_zero(data[i]));
}

void shift_down(std::array<Limb,binary_capacity>& data,std::size_t& size,
                unsigned bits) {
    const std::size_t whole=bits/64u;
    const unsigned part=bits%64u;
    const std::size_t remaining=size-whole;
    for (std::size_t i=0u;i<remaining;++i) {
        const Limb low=data[i+whole]>>part;
        const Limb high=part!=0u && i+whole+1u<size
            ? data[i+whole+1u]<<(64u-part) : Limb{0};
        data[i]=low|high;
    }
    size=normalized_size(data,remaining);
}

int compare_words(std::array<Limb,binary_capacity> const& a,std::size_t na,
                  std::array<Limb,binary_capacity> const& b,std::size_t nb) {
    if (na!=nb) { return na<nb ? -1 : 1; }
    for (std::size_t i=na;i!=0u;--i) {
        if (a[i-1u]!=b[i-1u]) { return a[i-1u]<b[i-1u] ? -1 : 1; }
    }
    return 0;
}

void subtract_words(std::array<Limb,binary_capacity>& a,std::size_t& na,
                    std::array<Limb,binary_capacity> const& b,std::size_t nb) {
    Limb borrow=0u;
    for (std::size_t i=0u;i<na;++i) {
        const Limb other=i<nb ? b[i] : Limb{0};
        const Limb difference=a[i]-other;
        const Limb borrowed=static_cast<Limb>(a[i]<other);
        a[i]=difference-borrow;
        borrow=static_cast<Limb>(borrowed!=0u || difference<borrow);
    }
    if (borrow!=0u) { throw std::logic_error("R24 unsigned underflow"); }
    na=normalized_size(a,na);
}

// Stein algorithm with disjoint, bounded stack magnitudes. This is a
// *diagnostic candidate*, not an architectural dispatch or benchmark oracle.
// Larger inputs use production gcd only for boundary/correctness testing.
BigUnsigned binary_gcd_stack(BigUnsigned const& a,BigUnsigned const& b) {
    if (a.is_zero()) { return b; }
    if (b.is_zero()) { return a; }
    if (a.size()>binary_capacity || b.size()>binary_capacity) {
        return gcd(a,b);
    }
    std::array<Limb,binary_capacity> x,y;
    std::copy(a.limbs().begin(),a.limbs().end(),x.begin());
    std::copy(b.limbs().begin(),b.limbs().end(),y.begin());
    std::size_t nx=a.size(),ny=b.size();
    const unsigned tx=trailing_bits(x,nx);
    const unsigned ty=trailing_bits(y,ny);
    const unsigned common=std::min(tx,ty);
    shift_down(x,nx,tx);
    shift_down(y,ny,ty);
    for (;;) {
        const int order=compare_words(x,nx,y,ny);
        if (order==0) { break; }
        if (order>0) {
            subtract_words(x,nx,y,ny);
            shift_down(x,nx,trailing_bits(x,nx));
        } else {
            subtract_words(y,ny,x,nx);
            shift_down(y,ny,trailing_bits(y,ny));
        }
    }
    return shift_left(BigUnsigned(BigUnsigned::LimbView(x.data(),nx)),common);
}

std::string positive_decimal(mpz_srcptr x) {
    RawInteger copy(x);
    mpz_abs(copy.value,copy.value);
    return decimal(copy.value);
}

struct Case {
    BigUnsigned a,b;
    RawInteger raw_a,raw_b,expected;
};

Case make_case(mpz_srcptr a,mpz_srcptr b) {
    Case result{from_decimal(positive_decimal(a)),from_decimal(positive_decimal(b)),
                RawInteger(a),RawInteger(b),RawInteger()};
    mpz_gcd(result.expected.value,a,b);
    return result;
}

struct Group {
    std::string name;
    std::vector<Case> cases;
};

std::vector<Group> make_groups(unsigned bits) {
    std::vector<Group> groups;
    for (auto const* name : {"construct_a","construct_b","add_den",
                             "add_reduce","sub_reduce","mul_cross_a",
                             "mul_cross_b","div_num","div_den"}) {
        groups.push_back(Group{name,{}});
        groups.back().cases.reserve(32u);
    }
    for (auto const& spec : AriadneBenchmark::rational_specifications(bits)) {
        const RawInteger an(spec.numerator_a),ad(spec.denominator_a);
        const RawInteger bn(spec.numerator_b),bd(spec.denominator_b);
        const RawRational a(an,ad),b(bn,bd);
        groups[0u].cases.push_back(make_case(an.value,ad.value));
        groups[1u].cases.push_back(make_case(bn.value,bd.value));
        groups[2u].cases.push_back(make_case(mpq_denref(a.value),mpq_denref(b.value)));
        groups[5u].cases.push_back(make_case(mpq_numref(a.value),mpq_denref(b.value)));
        groups[6u].cases.push_back(make_case(mpq_numref(b.value),mpq_denref(a.value)));
        groups[7u].cases.push_back(make_case(mpq_numref(a.value),mpq_numref(b.value)));
        groups[8u].cases.push_back(make_case(mpq_denref(a.value),mpq_denref(b.value)));
        RawInteger g,ad_by_g,bd_by_g,scaled_a,scaled_b,sum;
        mpz_gcd(g.value,mpq_denref(a.value),mpq_denref(b.value));
        mpz_divexact(ad_by_g.value,mpq_denref(a.value),g.value);
        mpz_divexact(bd_by_g.value,mpq_denref(b.value),g.value);
        mpz_mul(scaled_a.value,mpq_numref(a.value),bd_by_g.value);
        mpz_mul(scaled_b.value,mpq_numref(b.value),ad_by_g.value);
        // The production add/sub path skips the second GCD when
        // gcd(den_a,den_b)==1, and also when the cross-sum is zero.
        // Time this category only on inputs that actually reach it.
        if (mpz_cmp_ui(g.value,1ul)!=0) {
            mpz_add(sum.value,scaled_a.value,scaled_b.value);
            if (mpz_sgn(sum.value)!=0) {
                groups[3u].cases.push_back(make_case(sum.value,g.value));
            }
            mpz_sub(sum.value,scaled_a.value,scaled_b.value);
            if (mpz_sgn(sum.value)!=0) {
                groups[4u].cases.push_back(make_case(sum.value,g.value));
            }
        }
    }
    return groups;
}

void verify_case(Case const& item) {
    const auto expected=from_decimal(decimal(item.expected.value));
    if (compare(gcd(item.a,item.b),expected)!=0 ||
            compare(binary_gcd_stack(item.a,item.b),expected)!=0 ||
            compare(AriadneBenchmark::LehmerDiagnostic::gcd(item.a,item.b),expected)!=0) {
        throw std::runtime_error("R24 GCD mismatch with independent GMP oracle");
    }
}

void verify_boundary_cases() {
    std::mt19937_64 random(0x52'32'34'20'47'43'44ULL);
    for (unsigned i=0u;i<2048u;++i) {
        const std::size_t na=1u+static_cast<std::size_t>(i%9u);
        const std::size_t nb=1u+static_cast<std::size_t>((i/9u)%9u);
        std::array<Limb,9u> a{},b{};
        for (std::size_t j=0u;j<na;++j) { a[j]=random(); }
        for (std::size_t j=0u;j<nb;++j) { b[j]=random(); }
        a[na-1u]|=1u;
        b[nb-1u]|=1u;
        if (i%17u==0u) { std::fill(a.begin(),a.end(),0u); }
        if (i%19u==0u) { std::fill(b.begin(),b.end(),0u); }
        if (i%23u==0u) { b=a; }
        const BigUnsigned aa(BigUnsigned::LimbView(a.data(),na));
        const BigUnsigned bb(BigUnsigned::LimbView(b.data(),nb));
        RawInteger ga,gb,expected;
        mpz_import(ga.value,aa.size(),-1,sizeof(Limb),0,0,aa.limbs().data());
        mpz_import(gb.value,bb.size(),-1,sizeof(Limb),0,0,bb.limbs().data());
        mpz_gcd(expected.value,ga.value,gb.value);
        verify_case(make_case(ga.value,gb.value));
    }
}

void report_group(Group const& group,unsigned bits,
                  unsigned iterations,unsigned samples) {
    if (group.cases.empty()) { throw std::logic_error("R24 empty group"); }
    for (auto const& item : group.cases) { verify_case(item); }
    const auto index=[&](unsigned i) -> Case const& {
        return group.cases[static_cast<std::size_t>(i)%group.cases.size()];
    };
    const double native=exact_arithmetic_median_ns([&](unsigned i) {
        auto const& item=index(i);
        return gcd(item.a,item.b);
    },iterations,samples);
    const double gmp=exact_arithmetic_median_ns([&](unsigned i) {
        auto const& item=index(i);
        RawInteger result;
        mpz_gcd(result.value,item.raw_a.value,item.raw_b.value);
        return result;
    },iterations,samples);
    std::cout << bits << " " << group.name
              << " active_pairs=" << group.cases.size() << "/32"
              << " native=" << native
              << " gmp=" << gmp << " native/gmp=" << native/gmp;
    const bool bounded=std::all_of(group.cases.begin(),group.cases.end(),
        [](Case const& item) {
            return item.a.size()<=binary_capacity &&
                   item.b.size()<=binary_capacity;
        });
    if (bounded) {
        const double binary=exact_arithmetic_median_ns([&](unsigned i) {
            auto const& item=index(i);
            return binary_gcd_stack(item.a,item.b);
        },iterations,samples);
        std::cout << " binary=" << binary << " binary/native=" << binary/native;
    } else {
        std::cout << " binary=not-applicable";
    }
    if (bounded && AriadneBenchmark::LehmerDiagnostic::available) {
        const double lehmer=exact_arithmetic_median_ns([&](unsigned i) {
            auto const& item=index(i);
            return AriadneBenchmark::LehmerDiagnostic::gcd(item.a,item.b);
        },iterations,samples);
        std::cout << " lehmer=" << lehmer << " lehmer/native=" << lehmer/native;
    } else {
        std::cout << " lehmer=not-applicable";
    }
    std::cout << '\n';
}

} // namespace

int main(int argc,char** argv) {
    try {
        unsigned iterations=3000u,samples=11u;
        if (argc==3) {
            iterations=static_cast<unsigned>(std::stoul(argv[1]));
            samples=static_cast<unsigned>(std::stoul(argv[2]));
        } else if (argc!=1) {
            throw std::invalid_argument("usage: benchmark_rational_gcd_profile [iterations samples]");
        }
        if (iterations==0u || samples==0u || samples%2u==0u) {
            throw std::invalid_argument("iterations>0, odd samples>0 required");
        }
        verify_boundary_cases();
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "R24 GCD diagnostic; iterations=" << iterations
                  << " samples=" << samples
                  << "; all pairs match the unchanged Rational fixture\n";
        for (unsigned bits : {53u,128u,256u,1024u,2048u,4096u}) {
            for (auto const& group : make_groups(bits)) {
                report_group(group,bits,iterations,samples);
            }
        }
        std::cout << "R24 numerical verification PASS\n";
        return 0;
    } catch (std::exception const& e) {
        std::cerr << "R24 diagnostic FAILED: " << e.what() << '\n';
        return 1;
    }
}
