// SPDX-License-Identifier: MIT
// Reproduce A23/A24 test_views_and_scratch exact oracle fixtures.
// Standalone C++20; no Numeric, GMP, Boost or external runtime.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

Nat multiply(Nat left,Nat right) {
    Nat result;
    while (!right.zero()) {
        if ((right.w.front()&1u)!=0u) result=Nat::add(result,left);
        right=right.shr(1u);
        if (!right.zero()) left=left.shl(1u);
    }
    return result;
}

Nat low_limbs(Nat value,std::size_t count) {
    if (value.w.size()>count) value.w.resize(count);
    value.trim();
    return value;
}

std::array<std::string,5u> a23_buffer_arithmetic_fixtures() {
    std::array<std::string,5u> out;
    for (std::size_t part=0u;part<out.size();++part) {
        out[part]="# A23 buffer arithmetic visits "+
            std::to_string(part*1000u)+".."+
            std::to_string(part*1000u+999u)+".\n"
            "# repetition n m low lhs_hex rhs_hex sum_hex difference_hex product_low_hex\n";
    }
    const auto max=std::numeric_limits<std::uint64_t>::max();
    std::mt19937_64 engine(0xB22B22ULL);
    for (unsigned repetition=0u;repetition<5000u;++repetition) {
        const std::size_t n=engine()%28u;
        const std::size_t m=engine()%28u;
        std::vector<std::uint64_t> left(n),right(m);
        for (auto& word:left) word=engine();
        for (auto& word:right) word=engine();
        if (repetition%7u==0u) std::fill(left.begin(),left.end(),max);
        if (repetition%11u==0u) std::fill(right.begin(),right.end(),max);
        const Nat a(left),b(right);
        const std::size_t low=repetition%5u==0u ? 0u :
            repetition%5u==1u ? 1u :
            repetition%5u==2u ? 4u :
            repetition%5u==3u ? 15u : engine()%65u;
        const int comparison=a.compare(b);
        const std::string difference=comparison>=0 ?
            Nat::sub(a,b).hex() : "-";
        out[repetition/1000u]+=row({
            std::to_string(repetition),std::to_string(n),std::to_string(m),
            std::to_string(low),a.hex(),b.hex(),Nat::add(a,b).hex(),
            difference,low_limbs(multiply(a,b),low).hex()});
    }
    return out;
}
}

int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        const auto a23=a23_buffer_arithmetic_fixtures();
        for (std::size_t part=0u;part<a23.size();++part) {
            independent_oracle::output(
                dir,"big_unsigned_a23_views_"+std::to_string(part)+".tsv",
                a23[part],write);
        }
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"BigUnsigned views generator: "<<error.what()<<'\n';
        return 1;
    }
}
