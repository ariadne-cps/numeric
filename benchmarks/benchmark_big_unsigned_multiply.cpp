/***************************************************************************
 *            benchmark_big_unsigned_multiply.cpp
 *
 *  Copyright  2026  Ariadne contributors
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

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
#include <gmp.h>

namespace {

using Ariadne::BigUnsigned;
using Clock = std::chrono::steady_clock;

struct GmpInteger {
    mpz_t value;
    GmpInteger() { mpz_init(value); }
    ~GmpInteger() { mpz_clear(value); }
    GmpInteger(GmpInteger const&) = delete;
    GmpInteger& operator=(GmpInteger const&) = delete;
};

template<class Function>
double median_ns(Function&& function,
                 unsigned iterations,unsigned samples)
{
    std::vector<double> times;
    times.reserve(samples);
    function();
    for (unsigned sample=0u; sample<samples; ++sample) {
        const auto start = Clock::now();
        for (unsigned i=0u; i<iterations; ++i) {
            function();
        }
        const auto finish = Clock::now();
        const auto elapsed =
            std::chrono::duration<double,std::nano>(finish-start).count();
        times.push_back(elapsed/static_cast<double>(iterations));
    }
    std::sort(times.begin(),times.end());
    return times[times.size()/2u];
}

void import_gmp(mpz_ptr target,BigUnsigned const& x)
{
    const auto digits = x.limbs();
    mpz_import(target,digits.size(),-1,sizeof(BigUnsigned::Limb),
               0,0,digits.data());
}

void verify_gmp(BigUnsigned const& actual,mpz_srcptr expected)
{
    GmpInteger converted;
    import_gmp(converted.value,actual);
    if (mpz_cmp(converted.value,expected)!=0) {
        throw std::runtime_error("BigUnsigned multiply differs from GMP");
    }
}

void measure(unsigned bits,unsigned iterations,unsigned samples)
{
    std::mt19937_64 engine(0xB3B3B3ULL+bits);
    const std::size_t n=(bits+63u)/64u;
    std::vector<BigUnsigned::Limb> aa(n),bb(n);
    for (auto& v:aa) v=engine();
    for (auto& v:bb) v=engine();
    const auto top=static_cast<unsigned>(bits%64u);
    if (top!=0u) {
        aa.back() &= (std::uint64_t{1}<<top)-1u;
        bb.back() &= (std::uint64_t{1}<<top)-1u;
    }
    aa.back() |= std::uint64_t{1}<<((bits-1u)%64u);
    bb.back() |= std::uint64_t{1}<<((bits-1u)%64u);
    const BigUnsigned a(aa),b(bb);
    BigUnsigned::Scratch scratch;
    const auto result_size = a.size()+b.size();
    auto output=scratch.acquire(result_size);
    GmpInteger ga,gb,gr;
    import_gmp(ga.value,a);
    import_gmp(gb.value,b);
    mpz_mul(gr.value,ga.value,gb.value);
    BigUnsigned::mul_into_portable(a.limbs(),b.limbs(),output);
    verify_gmp(BigUnsigned(output),gr.value);
    BigUnsigned::mul_into(a.limbs(),b.limbs(),output);
    verify_gmp(BigUnsigned(output),gr.value);
    verify_gmp(mul(a,b),gr.value);

    const auto portable=median_ns([&] {
        BigUnsigned::mul_into_portable(a.limbs(),b.limbs(),output);
    },iterations,samples);
    const auto native=median_ns([&] {
        BigUnsigned::mul_into(a.limbs(),b.limbs(),output);
    },iterations,samples);
    const auto raw=median_ns([&] {
        mpz_mul(gr.value,ga.value,gb.value);
    },iterations,samples);
    const auto owned=median_ns([&] {
        const BigUnsigned product=mul(a,b);
        if (product.is_zero()) std::abort();
    },iterations,samples);
    // Verify after timing as well, with no exact checks inside timed code.
    verify_gmp(BigUnsigned(output),gr.value);
    verify_gmp(mul(a,b),gr.value);
    std::cout << std::setw(5) << bits
              << std::setw(16) << portable
              << std::setw(16) << native
              << std::setw(16) << raw
              << std::setw(16) << owned
              << std::setw(14) << native/raw
              << '\n';
}

} // namespace

int main(int argc,char** argv)
{
    try {
        const unsigned iterations = argc>1
            ? static_cast<unsigned>(std::stoul(argv[1])) : 20000u;
        const unsigned samples = argc>2
            ? static_cast<unsigned>(std::stoul(argv[2])) : 31u;
        if (iterations==0u || samples==0u) {
            throw std::invalid_argument("iterations and samples must be positive");
        }
        std::cout << "BigUnsigned multiplication; native available="
                  << BigUnsigned::has_native_multiply()
                  << "; iterations=" << iterations
                  << "; samples=" << samples << '\n';
        std::cout << "bits portable_ns selected_ns gmp_mpz_ns owned_ns selected/GMP\n";
        for (unsigned bits: {53u,128u,256u,1024u,2048u}) {
            measure(bits,iterations,samples);
        }
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
