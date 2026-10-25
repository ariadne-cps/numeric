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
// A3: exact 128/64 division edge and B26 fixture. Original input streams,
// including all conditional mt19937_64 draws, are preserved byte-for-byte.
std::string div128_fixture() {
    std::string out=
        "# A3 exact 128/64 quotient/remainder oracle, 192 edge + 2048 B26 visits.\n"
        "# kind ordinal divisor_hex low_hex high_hex numerator_hex quotient_hex remainder_hex\n"
        "# E: original 8*4*6 edge cases; B: original mt19937_64(0xB260D1B1) visits.\n"
        "# All original conditional PRNG draws and 128/64 preconditions retained.\n"
        "# Independent aligned-subtraction integer long division, no GMP/MPFR/Boost/Numeric.\n";
    std::size_t rows=0u;
    const auto append=[&](char kind,std::size_t ordinal,
                          std::uint64_t divisor,std::uint64_t low,
                          std::uint64_t high) {
        if (divisor==0u || high>=divisor)
            throw std::runtime_error("invalid A3 128/64 source");
        const Nat source(std::vector<std::uint64_t>{low,high});
        const Nat d(divisor);
        const auto [q,r]=divide(source,d);
        if (r.compare(d)>=0 ||
            Nat::add(multiply(q,d),r).compare(source)!=0)
            throw std::runtime_error("invalid A3 quotient/remainder identity");
        out+=row({std::string(1,kind),std::to_string(ordinal),
                  d.hex(),Nat(low).hex(),Nat(high).hex(),
                  source.hex(),q.hex(),r.hex()});
        ++rows;
    };
    const auto max=UINT64_MAX;
    const std::array<std::uint64_t,8> edge_divisors{
        1u,2u,3u,0xffffffffULL,0x100000000ULL,
        0x8000000000000000ULL,max-1u,max};
    const std::array<std::uint64_t,6> edge_lows{
        0u,1u,0xffffffffULL,0x8000000000000000ULL,max-1u,max};
    std::size_t edge_index=0u;
    for (const auto divisor:edge_divisors) {
        const std::array<std::uint64_t,4> highs{
            0u,divisor>1u?1u:0u,divisor/2u,divisor-1u};
        for (const auto high:highs)
            for (const auto low:edge_lows)
                append('E',edge_index++,divisor,low,high);
    }
    std::mt19937_64 wide_engine(0xB260D1B1ULL);
    for (unsigned repetition=0u;repetition<2048u;++repetition) {
        const std::uint64_t divisor=repetition%11u==0u ? max :
            repetition%11u==1u ? std::uint64_t{1}<<63u :
            repetition%11u==2u ? 3u :
            repetition%11u==3u ? std::uint64_t{1}<<32u :
            wide_engine()|std::uint64_t{1};
        const std::uint64_t high=repetition%5u==0u ?
            divisor-1u : repetition%5u==1u ? 0u :
            wide_engine()%divisor;
        const std::uint64_t low=repetition%7u==0u ?
            max : repetition%7u==1u ? 0u : wide_engine();
        append('B',repetition,divisor,low,high);
    }
    if (rows!=2240u) throw std::runtime_error("wrong A3 case count");
    return out;
}

std::string div_limb_random_fixture() {
    // A4: exact regeneration of all original 6,000 visits and RNG draws.
    std::string out=
        "# A4 exact original random BigUnsigned limb division, 6000 source-keyed visits.\n"
        "# repetition divisor_decimal dividend_hex quotient_hex remainder_hex\n"
        "# Full mt19937_64(0xD1D1D1) draw sequence including overwritten inputs.\n"
        "# Exact Euclidean arithmetic, 0 as canonical zero; no GMP/MPFR/Numeric.\n";
    std::mt19937_64 engine(0xD1D1D1ULL);
    std::size_t rows=0u;
    for (unsigned repetition=0u;repetition<6000u;++repetition) {
        const std::size_t count=repetition%23u==0u ? 0u :
                                repetition%29u==0u ? 64u :
                                repetition%7u+1u;
        std::vector<std::uint64_t> digits(count);
        for (auto& limb:digits) limb=engine();
        if (repetition%13u==0u)
            std::fill(digits.begin(),digits.end(),UINT64_MAX);
        if (repetition%17u==0u)
            std::fill(digits.begin(),digits.end(),0u);
        if (repetition%19u==0u && !digits.empty()) {
            std::fill(digits.begin(),digits.end(),0u);
            digits.back()=1u;
        }
        const std::uint64_t divisor=repetition%11u==0u ? 1u :
                                    repetition%11u==1u ? UINT64_MAX :
                                    repetition%11u==2u ? std::uint64_t{1}<<63u :
                                    repetition%11u==3u ? 2u :
                                    engine()|std::uint64_t{1};
        const Nat n(std::move(digits)),d(divisor);
        const auto [q,r]=divide(n,d);
        if (r.compare(d)>=0 ||
            Nat::add(multiply(q,d),r).compare(n)!=0)
            throw std::runtime_error("invalid A4 Euclidean identity");
        out+=row({std::to_string(repetition),std::to_string(divisor),
                  n.hex(),q.hex(),r.hex()});
        ++rows;
    }
    if (rows!=6000u) throw std::runtime_error("wrong A4 row count");
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"big_unsigned_power.tsv",power_fixture(),write);
        independent_oracle::output(dir,"big_unsigned_short_div.tsv",short_div_fixture(),write);
        independent_oracle::output(dir,"big_unsigned_div128by64.tsv",div128_fixture(),write);
        independent_oracle::output(dir,"big_unsigned_div_limb_random.tsv",div_limb_random_fixture(),write);
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"Power/division generator: "<<error.what()<<'\n';
        return 1;
    }
}
