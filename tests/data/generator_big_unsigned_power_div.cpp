// SPDX-License-Identifier: MIT
// Reproduce original D104 power and D103 short division frozen test fixtures.
// Standalone C++20; no Numeric, GMP, MPFR, Boost or other oracle dependency.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

// Multiply by binary shift/add rather than Numeric's multiplication kernels.
Nat multiply(Nat left, Nat right) {
    Nat result;
    while (!right.zero()) {
        if ((right.w.front() & 1u)!=0u) result=Nat::add(result,left);
        right=right.shr(1u);
        if (!right.zero()) left=left.shl(1u);
    }
    return result;
}
Nat power(Nat base, unsigned exponent) {
    Nat result(1u);
    while (exponent!=0u) {
        if ((exponent&1u)!=0u) result=multiply(result,base);
        exponent>>=1u;
        if (exponent!=0u) base=multiply(base,base);
    }
    return result;
}
std::size_t bit_length(Nat const& x) {
    return x.zero()?0u:(x.w.size()-1u)*64u +
        64u-static_cast<std::size_t>(std::countl_zero(x.w.back()));
}
// Quotient with aligned-subtraction long division, independent of Numeric.
std::pair<Nat,Nat> divide(Nat remainder, Nat const& divisor) {
    if (divisor.zero()) throw std::invalid_argument("zero divisor");
    Nat quotient;
    while (remainder.compare(divisor)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(divisor);
        Nat subtrahend=divisor.shl(shift);
        if (subtrahend.compare(remainder)>0) subtrahend=divisor.shl(--shift);
        remainder=Nat::sub(remainder,subtrahend);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {quotient,remainder};
}
std::string power_fixture() {
    std::string out=
        "# D104: test_power frozen exact integer references, 1513 original assertions on pow outputs.\n"
        "# kind key1 key2 exponent source_hex result_hex.\n"
        "# R=original limb count and mode, S=sample base, B=boundary bit count.\n"
        "# All sources and outputs canonical lowercase hexadecimal, zero is 0.\n"
        "# Generated from independent std::mt19937_64 exact state and arbitrary-precision integers.\n"
        "# Each result checked both by built-in exponentiation and repeated squaring.\n"
        "# No GMP/mpz/FloatAP/BigUnsigned or runtime transcendental oracle is used.\n";
    auto append=[&](char kind, std::size_t key1, std::size_t key2,
                    unsigned exponent,Nat const& base) {
        out+=row({std::string(1,kind),std::to_string(key1),std::to_string(key2),
                  std::to_string(exponent),base.hex(),power(base,exponent).hex()});
    };
    std::mt19937_64 generator(0xB7400A11ULL);
    for (std::size_t count:{0u,1u,2u,3u,4u,5u,8u,9u,16u,17u}) {
        for (unsigned mode=0u;mode<12u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=generator();
            if (mode%4u==0u) std::fill(words.begin(),words.end(),0u);
            if (mode%4u==1u) std::fill(words.begin(),words.end(),UINT64_MAX);
            if (mode%4u==2u && !words.empty()) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            }
            if (mode!=0u && mode!=1u && mode!=2u && mode!=3u && mode!=7u) continue;
            const Nat base(std::move(words));
            for (unsigned exponent:{0u,1u,2u,8u,16u,17u})
                append('R',count,mode,exponent,base);
        }
    }
    const Nat sample(std::vector<std::uint64_t>{0xfffffffffffffffdULL,17u,1u});
    for (unsigned exponent:{0u,1u,2u,3u,5u,6u,7u,8u,9u,
                            12u,15u,16u,17u,31u,32u,33u})
        append('S',0u,0u,exponent,sample);
    for (unsigned bits:{1u,2u,53u,64u,65u,127u,128u,129u,
                        255u,256u,257u,341u,342u,511u,512u,513u,
                        1023u,1024u,1025u}) {
        std::vector<std::uint64_t> words((bits+63u)/64u,UINT64_MAX);
        if (bits%64u!=0u) words.back()=(std::uint64_t{1}<<(bits%64u))-1u;
        const Nat base(std::move(words));
        const unsigned limit=1024u/bits;
        for (unsigned exponent:{limit==0u?0u:limit-1u,limit,limit+1u})
            append('B',bits,0u,exponent,base);
    }
    return out;
}
std::string short_div_fixture() {
    std::string out=
        "# D103 BigUnsigned short-owned division; frozen independent native-integer oracle.\n"
        "# n m mode original_left_hex original_right_hex quotient_hex remainder_hex\n"
        "# 864 input visits: standardized mt19937_64 seed 0xB74D1A1DE, 6x6x24.\n"
        "# Original source words and all mode%6 transformations are preserved.\n"
        "# Canonical lowercase hexadecimal, zero is 0, division-by-zero results are - -.\n"
        "# Expected values certified by native arbitrary-precision integer quotient/remainder,\n"
        "# and checked with a=b*q+r and 0<=r<b. No GMP, MPFR, FloatAP, BigUnsigned oracle.\n";
    std::mt19937_64 generator(0xB74D1A1DEULL);
    for (std::size_t n=0u;n<=5u;++n) {
        for (std::size_t m=0u;m<=5u;++m) {
            for (unsigned mode=0u;mode<24u;++mode) {
                std::vector<std::uint64_t> left(n),right(m);
                for (auto& word:left) word=generator();
                for (auto& word:right) word=generator();
                if (mode%6u==0u) std::fill(left.begin(),left.end(),0u);
                if (mode%6u==1u) {
                    std::fill(left.begin(),left.end(),UINT64_MAX);
                    std::fill(right.begin(),right.end(),1u);
                }
                if (mode%6u==2u && !right.empty()) {
                    std::fill(right.begin(),right.end(),0u);
                    right.back()=std::uint64_t{1}<<(mode%64u);
                }
                if (mode%6u==3u && n==m) right=left;
                if (mode%6u==4u && n==m && !right.empty()) {
                    right=left;
                    right.front()=0u;
                    left.front()=UINT64_MAX;
                }
                if (mode!=0u && mode!=1u && mode!=2u &&
                    mode!=3u && mode!=4u && mode!=5u &&
                    mode!=11u && mode!=17u && mode!=23u) continue;
                const Nat a(std::move(left)),b(std::move(right));
                std::string q="-",r="-";
                if (!b.zero()) {
                    auto [quotient,remainder]=divide(a,b);
                    // Independently certify the exact Euclidean identity.
                    if (remainder.compare(b)>=0 ||
                        Nat::add(multiply(b,quotient),remainder).compare(a)!=0)
                        throw std::runtime_error("invalid division oracle");
                    q=quotient.hex();
                    r=remainder.hex();
                }
                out+=row({std::to_string(n),std::to_string(m),std::to_string(mode),
                          a.hex(),b.hex(),q,r});
            }
        }
    }
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"big_unsigned_power.tsv",power_fixture(),write);
        independent_oracle::output(dir,"big_unsigned_short_div.tsv",short_div_fixture(),write);
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"Power/division generator: "<<error.what()<<'\n';
        return 1;
    }
}
