// SPDX-License-Identifier: MIT
// Reproduce B134/B135 test_views_and_scratch exact oracle fixtures.
// Standalone C++20; no Numeric, GMP, Boost or external runtime.
#include <limits>
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

std::array<std::string,5u> b134_buffer_arithmetic_fixtures() {
    std::array<std::string,5u> out;
    for (std::size_t part=0u;part<out.size();++part) {
        out[part]="# B134 buffer arithmetic visits "+
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

std::size_t b135_part(std::size_t ordinal) {
    return ordinal<5775u ? ordinal/825u :
        7u+(ordinal-5775u)/275u;
}

struct A24Fixtures {
    std::string sources;
    std::array<std::string,10u> expected;
};

A24Fixtures b135_shifted_fixtures() {
    A24Fixtures result;
    result.sources=
        "# B135 R38 shifted-operation base sources.\n"
        "# source_ordinal width case_index length_b a_hex b_hex\n";
    for (std::size_t part=0u;part<result.expected.size();++part) {
        const std::size_t first=part<7u ? part*825u :
            5775u+(part-7u)*275u;
        const std::size_t count=part<7u ? 825u : 275u;
        result.expected[part]=
            "# B135 R38 shifted-operation visits "+
            std::to_string(first)+".."+
            std::to_string(first+count-1u)+".\n"
            "# ordinal width case_index sa sb comparison sum_hex difference_hex\n";
    }
    const auto max=std::numeric_limits<std::uint64_t>::max();
    std::mt19937_64 engine(0xB38038ULL);
    std::size_t source_ordinal=0u,ordinal=0u;
    for (const std::size_t width:
            {0u,1u,4u,16u,17u,31u,32u,33u,64u,65u,129u}) {
        for (unsigned case_index=0u;case_index<24u;++case_index) {
            const std::size_t length_b=case_index%5u==0u ? width :
                case_index%5u==1u ? 0u :
                case_index%5u==2u ? 1u :
                case_index%5u==3u ? width/2u : width+1u;
            std::vector<std::uint64_t> left(width),right(length_b);
            for (auto& word:left) word=engine();
            for (auto& word:right) word=engine();
            if (case_index%7u==0u) {
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
            }
            const Nat a(left),b(right);
            result.sources+=row({
                std::to_string(source_ordinal++),std::to_string(width),
                std::to_string(case_index),std::to_string(length_b),
                a.hex(),b.hex()});
            for (const unsigned sa:{0u,1u,7u,31u,63u}) {
                for (const unsigned sb:{0u,1u,7u,31u,63u}) {
                    const Nat shifted_a=a.shl(sa);
                    const Nat shifted_b=b.shl(sb);
                    const int comparison=shifted_a.compare(shifted_b);
                    const Nat difference=comparison>=0 ?
                        Nat::sub(shifted_a,shifted_b) :
                        Nat::sub(shifted_b,shifted_a);
                    // Freeze the segment before row construction increments ordinal.
                    const auto part=b135_part(ordinal);
                    result.expected[part]+=row({
                        std::to_string(ordinal++),std::to_string(width),
                        std::to_string(case_index),std::to_string(sa),
                        std::to_string(sb),std::to_string(comparison),
                        Nat::add(shifted_a,shifted_b).hex(),difference.hex()});
                }
            }
        }
    }
    if (source_ordinal!=264u || ordinal!=6600u)
        throw std::runtime_error("wrong B135 shifted case count");
    return result;
}
}

int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        const auto a23=b134_buffer_arithmetic_fixtures();
        for (std::size_t part=0u;part<a23.size();++part) {
            independent_oracle::output(
                dir,"big_unsigned_b134_views_"+std::to_string(part)+".tsv",
                a23[part],write);
        }
        const auto a24=b135_shifted_fixtures();
        independent_oracle::output(
            dir,"big_unsigned_b135_shifted_sources.tsv",
            a24.sources,write);
        for (std::size_t part=0u;part<a24.expected.size();++part) {
            independent_oracle::output(
                dir,"big_unsigned_b135_shifted_"+std::to_string(part)+".tsv",
                a24.expected[part],write);
        }
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"BigUnsigned views generator: "<<error.what()<<'\n';
        return 1;
    }
}
