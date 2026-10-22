/***************************************************************************
 *            benchmark_integer_alias_dispatch.cpp
 *
 *  Copyright 2026 Ariadne contributors
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

// B34 diagnostic only: never modifies the official Integer/GMP workload.
// Compile/link with the real NUMERIC_SRC object to measure B33 production
// dispatch on the SAME CPU/compiler as its generic-multiply reference.
#include "numeric/integer.hpp"
#include "numeric/detail/integer_gmp_bridge.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>
#include <gmp.h>

namespace {

using namespace Ariadne;
using Clock=std::chrono::steady_clock;

struct Reference {
    mpz_t value;
    Reference() { mpz_init(value); }
    ~Reference() { mpz_clear(value); }
    Reference(Reference const&)=delete;
    Reference& operator=(Reference const&)=delete;
};

void import_signed(mpz_ptr target,Integer const& input)
{
    const auto limbs=input.magnitude().limbs();
    mpz_import(target,limbs.size(),-1,sizeof(BigUnsigned::Limb),
               0,0,limbs.data());
    if (input.negative()) mpz_neg(target,target);
}

Integer generic_multiply(Integer const& left,Integer const& right)
{
    return Integer::from_magnitude(
        mul(left.magnitude(),right.magnitude()),
        left.negative()!=right.negative());
}

Integer native_magnitude_square(Integer const& input)
{
    return Integer::from_magnitude(square(input.magnitude()),false);
}

void verify(Integer const& left,Integer const& right)
{
    Reference a,b,expected;
    import_signed(a.value,left);
    import_signed(b.value,right);
    mpz_mul(expected.value,a.value,b.value);
    const GmpIntegerBoundary native(left*right);
    const GmpIntegerBoundary generic(generic_multiply(left,right));
    if (mpz_cmp(native.value(),expected.value)!=0 ||
        mpz_cmp(generic.value(),expected.value)!=0) {
        throw std::runtime_error("Integer alias diagnostic: GMP mismatch");
    }
    mpz_mul(expected.value,a.value,a.value);
    const GmpIntegerBoundary aliased(left*left);
    const GmpIntegerBoundary square_value(native_magnitude_square(left));
    const GmpIntegerBoundary direct(generic_multiply(left,left));
    if (mpz_cmp(aliased.value(),expected.value)!=0 ||
        mpz_cmp(square_value.value(),expected.value)!=0 ||
        mpz_cmp(direct.value(),expected.value)!=0) {
        throw std::runtime_error("Integer self-product: GMP mismatch");
    }
}

struct Pools {
    std::vector<Integer> left,right,copies;
};

Pools make_pools(unsigned bits)
{
    constexpr std::size_t count=32u;
    std::mt19937_64 generator(0xB34000ULL+bits);
    Pools pools;
    pools.left.reserve(count);
    pools.right.reserve(count);
    pools.copies.reserve(count);
    const std::size_t words=(bits+63u)/64u;
    const unsigned top_bits=bits%64u;
    for (std::size_t index=0u;index<count;++index) {
        std::vector<std::uint64_t> a(words),b(words);
        for (auto& word:a) word=generator();
        for (auto& word:b) word=generator();
        if (top_bits!=0u) {
            const auto mask=(std::uint64_t{1}<<top_bits)-1u;
            a.back()&=mask;
            b.back()&=mask;
            a.back()|=std::uint64_t{1}<<(top_bits-1u);
            b.back()|=std::uint64_t{1}<<(top_bits-1u);
        } else {
            a.back()|=std::uint64_t{1}<<63u;
            b.back()|=std::uint64_t{1}<<63u;
        }
        pools.left.push_back(Integer::from_magnitude(
            BigUnsigned(std::span(a)),index%2u!=0u));
        pools.right.push_back(Integer::from_magnitude(
            BigUnsigned(std::span(b)),index%3u!=0u));
        pools.copies.push_back(pools.left.back());
        verify(pools.left.back(),pools.right.back());
        verify(pools.left.back(),pools.copies.back());
    }
    return pools;
}

volatile std::uint64_t consumed=0u;

template<class Operation>
double sample_ns(Operation&& operation,Pools const& pool,
                 unsigned iterations,unsigned sample_number)
{
    std::uint64_t checksum=0u;
    const auto start=Clock::now();
    for (unsigned rep=0u;rep<iterations;++rep) {
        const std::size_t index=(static_cast<std::size_t>(rep)+
                                 sample_number*7u)%pool.left.size();
        Integer result=operation(index);
        checksum^=result.magnitude().limb(0u);
    }
    const auto stop=Clock::now();
    consumed=consumed^checksum;
    return std::chrono::duration<double,std::nano>(stop-start).count()/
           static_cast<double>(iterations);
}

enum class Mode {
    native_distinct, generic_distinct, native_copy, generic_copy,
    native_alias, generic_alias, direct_square
};

constexpr std::array<Mode,7> modes{
    Mode::native_distinct,Mode::generic_distinct,
    Mode::native_copy,Mode::generic_copy,
    Mode::native_alias,Mode::generic_alias,Mode::direct_square
};

double measure(Mode mode,Pools const& pool,
               unsigned iterations,unsigned sample_number)
{
    auto const& a=pool.left;
    auto const& b=pool.right;
    auto const& c=pool.copies;
    switch (mode) {
    case Mode::native_distinct:
        return sample_ns([&](std::size_t i) {return a[i]*b[i];},
                         pool,iterations,sample_number);
    case Mode::generic_distinct:
        return sample_ns([&](std::size_t i) {
            return generic_multiply(a[i],b[i]);},
            pool,iterations,sample_number);
    case Mode::native_copy:
        return sample_ns([&](std::size_t i) {return a[i]*c[i];},
                         pool,iterations,sample_number);
    case Mode::generic_copy:
        return sample_ns([&](std::size_t i) {
            return generic_multiply(a[i],c[i]);},
            pool,iterations,sample_number);
    case Mode::native_alias:
        return sample_ns([&](std::size_t i) {return a[i]*a[i];},
                         pool,iterations,sample_number);
    case Mode::generic_alias:
        return sample_ns([&](std::size_t i) {
            return generic_multiply(a[i],a[i]);},
            pool,iterations,sample_number);
    case Mode::direct_square:
        return sample_ns([&](std::size_t i) {
            return native_magnitude_square(a[i]);},
            pool,iterations,sample_number);
    default:
        throw std::logic_error("Unknown diagnostic mode");
    }
}

unsigned parse_positive(const char* arg)
{
    const auto value=std::stoul(arg);
    if (value==0u || value>std::numeric_limits<unsigned>::max()) {
        throw std::invalid_argument("Expected a positive unsigned count");
    }
    return static_cast<unsigned>(value);
}

int run(int argc,char const* const* argv)
{
    unsigned iterations=20000u,samples=31u;
    std::vector<unsigned> bits{53u,128u,256u,512u,1024u,2048u,4096u};
    for (int i=1;i<argc;++i) {
        const std::string option=argv[i];
        if (i+1>=argc) {
            throw std::invalid_argument("Expected value after "+option);
        }
        const std::string argument=argv[++i];
        if (option=="--iterations") iterations=parse_positive(argument.c_str());
        else if (option=="--samples") samples=parse_positive(argument.c_str());
        else if (option=="--bits") {
            bits.clear();
            std::size_t pos=0u;
            while (pos<argument.size()) {
                const std::size_t end=argument.find(',',pos);
                const std::string segment=argument.substr(
                    pos,end==std::string::npos?end:end-pos);
                bits.push_back(parse_positive(segment.c_str()));
                if (end==std::string::npos) break;
                pos=end+1u;
            }
        } else {
            throw std::invalid_argument("Unknown option "+option);
        }
    }
    if (bits.empty() || samples<3u || samples%2u==0u) {
        throw std::invalid_argument("Need sizes and odd samples >=3");
    }

    std::cout<<"B34 alias diagnostic (not the official Integer benchmark)"
             <<" | iterations="<<iterations<<" | samples="<<samples
             <<" | fused="<<BigUnsigned::selected_multiply_uses_fused()
             <<"\n";
    std::cout<<"All measured results are verified against GMP before timing.\n";
    std::cout<<std::fixed<<std::setprecision(3);
    for (unsigned width:bits) {
        const auto inputs=make_pools(width);
        std::array<std::vector<double>,modes.size()> readings;
        for (auto& values:readings) values.reserve(samples);
        // Rotate the measurement order between samples to limit drift.
        for (unsigned sample=0u;sample<samples;++sample) {
            for (std::size_t step=0u;step<modes.size();++step) {
                const std::size_t k=(step+sample)%modes.size();
                readings[k].push_back(
                    measure(modes[k],inputs,iterations,sample));
            }
        }
        std::array<double,modes.size()> medians{};
        for (std::size_t k=0u;k<modes.size();++k) {
            auto& values=readings[k];
            std::sort(values.begin(),values.end());
            medians[k]=values[values.size()/2u];
        }
        std::cout<<"bits="<<width
                 <<" distinct(native/generic)="<<medians[0]<<"/"<<medians[1]
                 <<" equal_copy(native/generic)="<<medians[2]<<"/"<<medians[3]
                 <<" alias(native/generic/square)="<<medians[4]
                 <<"/"<<medians[5]<<"/"<<medians[6]<<" ns\n";
    }
    std::cout<<"GMP verification PASS\n";
    return 0;
}

} // namespace

int main(int argc,char const* argv[])
{
    try { return run(argc,argv); }
    catch (std::exception const& e) {
        std::cerr<<"B34 diagnostic ERROR: "<<e.what()<<"\n";
        return 1;
    }
}
