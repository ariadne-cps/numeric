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

// A5: a separate portable base-2 bitwise division oracle for wide
// dividend / 64-bit divisor. No production divide_wide or GMP involved.
// The 65th remainder bit is tracked as a carry before modulo-2^64
// subtraction, avoiding compiler-specific 128-bit integer types.
std::pair<Nat,std::uint64_t> divide_limb_bits(Nat const& n,
                                               std::uint64_t divisor) {
    if (divisor==0u) throw std::invalid_argument("zero limb divisor");
    std::vector<std::uint64_t> quotient(n.w.size(),0u);
    std::uint64_t remainder=0u;
    for (std::size_t i=n.w.size();i>0u;--i) {
        for (unsigned bit=64u;bit>0u;) {
            --bit;
            const bool carry=(remainder>>63u)!=0u;
            remainder=(remainder<<1u)|((n.w[i-1u]>>bit)&1u);
            if (carry || remainder>=divisor) {
                remainder-=divisor;
                quotient[i-1u]|=std::uint64_t{1}<<bit;
            }
        }
    }
    return {Nat(std::move(quotient)),remainder};
}
std::array<std::string,4> limb_stress_fixtures() {
    std::array<std::string,4> out{
        "# A5 B48 independent exact limb division, all 560 original visits.\n"
        "# count divisor_decimal mode dividend_hex quotient_hex remainder_hex\n"
        "# Original mt19937_64(0xB48032D0) with all overwritten draws preserved.\n"
        "# Exact independent Euclidean division; no Numeric, GMP, MPFR or Boost.\n",
        "# A5 B49 independent exact limb division, all 784 original visits.\n"
        "# count divisor_decimal mode dividend_hex quotient_hex remainder_hex\n"
        "# Original mt19937_64(0xB490C0DE) with all overwritten draws preserved.\n"
        "# Exact independent Euclidean division; no Numeric, GMP, MPFR or Boost.\n",
        "# A5 B57 independent exact limb division, all 1800 original visits.\n"
        "# count divisor_decimal mode dividend_hex quotient_hex remainder_hex\n"
        "# Original mt19937_64(0xB570641) with all overwritten draws preserved.\n"
        "# Exact independent Euclidean division; no Numeric, GMP, MPFR or Boost.\n",
        "# A5 B61 independent exact limb division, all 1560 original visits.\n"
        "# count divisor_decimal mode dividend_hex quotient_hex remainder_hex\n"
        "# Original mt19937_64(0xB610D1A6) with all overwritten draws preserved.\n"
        "# Exact independent Euclidean division; no Numeric, GMP, MPFR or Boost.\n"};
    std::array<std::size_t,4> rows{};
    const auto append=[&](std::size_t group,std::size_t count,
                          std::uint64_t divisor,unsigned mode,
                          std::vector<std::uint64_t> const& words) {
        if (words.size()!=count)
            throw std::runtime_error("wrong A5 source limb count");
        const Nat input(words);
        const auto [quotient,remainder]=divide_limb_bits(input,divisor);
        if (remainder>=divisor ||
            Nat::add(multiply(quotient,Nat(divisor)),Nat(remainder))
                .compare(input)!=0)
            throw std::runtime_error("wrong independent A5 Euclidean identity");
        out[group]+=row({std::to_string(count),std::to_string(divisor),
                         std::to_string(mode),input.hex(),quotient.hex(),
                         Nat(remainder).hex()});
        ++rows[group];
    };
    const auto max=UINT64_MAX;
    std::mt19937_64 narrow_engine(0xB48032D0ULL);
    for (std::size_t count:{1u,2u,3u,4u,5u,8u,16u,32u,64u,96u}) {
        for (std::uint64_t divisor:
                {1ULL,2ULL,3ULL,7ULL,9ULL,10ULL,31ULL,63ULL,
                 255ULL,65535ULL,2147483648ULL,4294967295ULL,
                 4294967296ULL,18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<4u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& word:words) word=narrow_engine();
                if (mode==0u) std::fill(words.begin(),words.end(),max);
                else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==2u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=max;
                    words.back()|=1u;
                }
                append(0u,count,divisor,mode,words);
            }
        }
    }
    std::mt19937_64 reciprocal_engine(0xB490C0DEULL);
    for (std::size_t count:{7u,8u,9u,16u,32u,64u,96u,128u}) {
        for (std::uint64_t divisor:
                {1ULL,2ULL,3ULL,5ULL,7ULL,9ULL,31ULL,
                 63ULL,255ULL,65535ULL,2147483648ULL,
                 4294967295ULL,4294967296ULL,
                 18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<7u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& word:words) word=reciprocal_engine();
                if (mode==0u) std::fill(words.begin(),words.end(),max);
                else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==2u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=max;
                    words.back()=1u;
                } else if (mode==3u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=std::uint64_t{1}<<63u;
                } else if (mode==4u) {
                    std::fill(words.begin(),words.end(),0u);
                } else if (mode==5u) {
                    for (std::size_t i=0u;i<count;i+=2u)
                        words[i]=0u;
                }
                append(1u,count,divisor,mode,words);
            }
        }
    }
    const auto pi1=[&](std::size_t group,std::uint64_t seed,
                       std::initializer_list<std::size_t> counts) {
        std::mt19937_64 engine(seed);
        for (std::size_t count:counts) {
            for (std::uint64_t divisor:
                    {1ULL,2ULL,3ULL,5ULL,7ULL,35ULL,255ULL,
                     65535ULL,4294967295ULL,4294967296ULL,
                     4294967297ULL,9223372036854775807ULL,
                     9223372036854775808ULL,9223372036854775809ULL,
                     18446744073709551615ULL}) {
                for (unsigned mode=0u;mode<8u;++mode) {
                    std::vector<std::uint64_t> words(count);
                    for (auto& word:words) word=engine();
                    if (mode==0u) std::fill(words.begin(),words.end(),max);
                    else if (mode==1u) std::fill(words.begin(),words.end(),0u);
                    else if (mode==2u && !words.empty()) {
                        std::fill(words.begin(),words.end(),0u);
                        words.front()=1u;
                    } else if (mode==3u && !words.empty()) {
                        std::fill(words.begin(),words.end(),0u);
                        words.back()=1u;
                    } else if (mode==4u && !words.empty()) {
                        std::fill(words.begin(),words.end(),0u);
                        words.back()=max;
                    } else if (mode==5u) {
                        for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
                    } else if (mode==6u) {
                        for (std::size_t i=0u;i<count;i+=3u) words[i]=max;
                    }
                    append(group,count,divisor,mode,words);
                }
            }
        }
    };
    pi1(2u,0xB570641ULL,
        {0u,1u,7u,8u,15u,16u,17u,31u,32u,33u,63u,64u,65u,96u,128u});
    pi1(3u,0xB610D1A6ULL,
        {0u,1u,7u,8u,9u,15u,16u,31u,32u,63u,64u,65u,128u});
    if (rows!=std::array<std::size_t,4>{560u,784u,1800u,1560u})
        throw std::runtime_error("wrong A5 fixture row counts");
    return out;
}

// A7: exact-limb dividend sources must be computed outside production
// arithmetic. Keep the complete original PRNG stream and all 6,200
// quotient/divisor visits, including overwritten random input words.
std::array<std::string,4> exact_limb_fixtures() {
    const auto header=[](std::string const& segment) {
        return std::string(
            "# A7 exact single-limb division, original 6,200 visits.\n"
            "# length divisor_index mode divisor_decimal quotient_hex dividend_hex\n"
            "# ")+segment+"\n"
            "# Original mt19937_64(0x22E1AC7) including overwritten word draws.\n"
            "# Products independently obtained by integer multiplication, zero canonical.\n";
    };
    std::array<std::string,4> out{
        header("Segment S: original lengths 0,1,2,4,5,8,9,16; 4,960 visits."),
        header("Segment M: original length 65; 620 visits."),
        header("Segment H0: original length 129, divisor indices 0..77; 312 visits."),
        header("Segment H1: original length 129, divisor indices 78..154; 308 visits.")};
    std::array<std::size_t,4> rows{};
    const std::uint64_t max=UINT64_MAX;
    std::mt19937_64 random(0x22E1AC7ULL);
    std::vector<std::uint64_t> divisors={
        1u,2u,3u,6u,7u,10u,12u,17u,65537u,max,max-1u};
    for (unsigned bit=0u;bit<64u;++bit) {
        divisors.push_back(std::uint64_t{1}<<bit);
        divisors.push_back((random()>>(63u-bit))|(std::uint64_t{1}<<bit));
    }
    for (std::uint64_t odd=1u;odd<32u;odd+=2u)
        divisors.push_back(odd);
    if (divisors.size()!=155u)
        throw std::runtime_error("invalid A7 divisor count");
    for (std::size_t length:{0u,1u,2u,4u,5u,8u,9u,16u,65u,129u}) {
        for (std::size_t index=0u;index<divisors.size();++index) {
            const std::uint64_t divisor=divisors[index];
            for (unsigned mode=0u;mode<4u;++mode) {
                std::vector<std::uint64_t> words(length);
                for (auto& word:words) word=random();
                if (mode==1u) std::fill(words.begin(),words.end(),max);
                if (mode==2u) std::fill(words.begin(),words.end(),0u);
                if (mode==3u) {
                    std::fill(words.begin(),words.end(),0u);
                    if (!words.empty()) {
                        words.back()=1u;
                        words.front()|=1u;
                    }
                }
                const Nat quotient(std::move(words));
                const Nat dividend=quotient.mul_small(divisor);
                // Validate via a separately implemented bitwise
                // quotient/remainder algorithm, not GMP or Numeric.
                const auto [verified,rem]=divide_limb_bits(dividend,divisor);
                if (rem!=0u || verified.compare(quotient)!=0)
                    throw std::runtime_error("A7 exact-product division identity");
                const std::size_t group=length<=16u ? 0u :
                     length==65u ? 1u : index<78u ? 2u : 3u;
                out[group]+=row({std::to_string(length),std::to_string(index),
                                 std::to_string(mode),std::to_string(divisor),
                                 quotient.hex(),dividend.hex()});
                ++rows[group];
            }
        }
    }
    if (rows!=std::array<std::size_t,4>{4960u,620u,312u,308u})
        throw std::runtime_error("wrong A7 fixture row count");
    return out;
}

// A8: B57 multi-limb division by powers of two. The original
// mt19937_64 stream is replayed in full; quotient/remainder are
// computed via standalone Nat shifts, never Numeric's division.
std::array<std::string,3> b57_power_fixtures() {
    const auto header=[](std::size_t visits) {
        return std::string("# A8 B57 division by 2^(64*(m-1)+bit), ")+
            std::to_string(visits)+
            " visits.\n"
            "# nlimbs mlimbs bit mode dividend_hex divisor_hex quotient_hex remainder_hex\n"
            "# Original mt19937_64(0xB570B172), all draws including overwritten words.\n"
            "# Independent shift/mask exact division, canonical lowercase hexadecimal.\n";
    };
    std::array<std::string,3> out{header(1200u),header(480u),header(480u)};
    std::array<std::size_t,3> rows{};
    std::mt19937_64 engine(0xB570B172ULL);
    for (std::size_t n:{1u,2u,8u,16u,32u,64u,65u,96u,128u}) {
        for (std::size_t m:{1u,2u,3u,8u,16u,32u,64u,65u}) {
            for (unsigned bit:{0u,1u,17u,32u,63u}) {
                for (unsigned mode=0u;mode<6u;++mode) {
                    std::vector<std::uint64_t> left(n);
                    for (auto& limb:left) limb=engine();
                    if (mode==0u) std::fill(left.begin(),left.end(),UINT64_MAX);
                    else if (mode==1u) std::fill(left.begin(),left.end(),0u);
                    else if (mode==2u) {
                        std::fill(left.begin(),left.end(),0u);
                        left.front()=1u;
                    } else if (mode==3u) {
                        std::fill(left.begin(),left.end(),0u);
                        left.back()=1u;
                    } else if (mode==4u) {
                        for (std::size_t i=0u;i<n;i+=2u) left[i]=0u;
                    }
                    const Nat dividend(std::move(left));
                    const std::size_t shift=64u*(m-1u)+bit;
                    const Nat divisor=Nat::power2(shift);
                    const Nat quotient=dividend.shr(shift);
                    const Nat remainder=Nat::sub(dividend,quotient.shl(shift));
                    if (remainder.compare(divisor)>=0 ||
                        Nat::add(quotient.shl(shift),remainder).compare(dividend)!=0)
                        throw std::runtime_error("invalid A8 Euclidean identity");
                    const std::size_t group=n<=32u?0u:n<=65u?1u:2u;
                    out[group]+=row({
                        std::to_string(n),std::to_string(m),
                        std::to_string(bit),std::to_string(mode),
                        dividend.hex(),divisor.hex(),quotient.hex(),
                        remainder.hex()});
                    ++rows[group];
                }
            }
        }
    }
    if (rows!=std::array<std::size_t,3>{1200u,480u,480u})
        throw std::runtime_error("wrong A8 visit counts");
    return out;
}

// A9: B55 native divisor-length selection cases. Restore original
// randomized inputs and calculate division by independent Nat long division.
std::array<std::string,3> b55_fixtures() {
    const auto header=[](std::size_t visits) {
        return std::string("# A9 B55 seeded multiprecision division oracle, ")+
            std::to_string(visits)+" original visits.\n"
            "# m extra mode dividend_hex divisor_hex quotient_hex remainder_hex\n"
            "# Original mt19937_64(0xB550DEC1), all conditional overwrites preserved.\n"
            "# Exact independent Euclidean quotient/remainder; no GMP/MPFR/Numeric.\n";
    };
    std::array<std::string,3> out{header(216u),header(144u),header(72u)};
    std::array<std::size_t,3> rows{};
    std::mt19937_64 engine(0xB550DEC1ULL);
    const auto max=UINT64_MAX;
    for (std::size_t m:{7u,8u,9u,16u,32u,64u}) {
        for (std::size_t extra:{2u,9u,32u}) {
            const std::size_t n=m+extra;
            for (unsigned mode=0u;mode<24u;++mode) {
                std::vector<std::uint64_t> left(n),right(m);
                for (auto& word:left) word=engine();
                for (auto& word:right) word=engine();
                switch (mode%8u) {
                case 0u:
                    std::fill(left.begin(),left.end(),max);
                    std::fill(right.begin(),right.end(),max);
                    break;
                case 1u:
                    std::fill(left.begin(),left.end(),0u);
                    left.back()=max;
                    std::fill(right.begin(),right.end(),0u);
                    right.back()=max;
                    break;
                case 2u:
                    std::fill(left.begin(),left.end(),0u);
                    left.front()=max;
                    left.back()=max;
                    std::fill(right.begin(),right.end(),0u);
                    right.front()=3u;
                    right.back()=max;
                    break;
                case 3u:
                    for (std::size_t i=0u;i<n;i+=2u) left[i]=0u;
                    for (std::size_t i=0u;i<m;i+=2u) right[i]=0u;
                    break;
                case 4u:
                    for (std::size_t i=1u;i+1u<n;++i) left[i]=0u;
                    for (std::size_t i=1u;i+1u<m;++i) right[i]=0u;
                    break;
                case 5u:
                    std::fill_n(left.begin(),n/2u,0u);
                    std::fill_n(right.begin(),m/2u,0u);
                    break;
                case 6u:
                    for (std::size_t i=1u;i<n;i+=3u) left[i]=max;
                    for (std::size_t i=1u;i<m;i+=3u) right[i]=max;
                    break;
                default: break;
                }
                left.back()|=std::uint64_t{1}<<63u;
                right.back()|=std::uint64_t{1}<<63u;
                const Nat dividend(std::move(left)),divisor(std::move(right));
                const auto [quotient,remainder]=divide(dividend,divisor);
                if (remainder.compare(divisor)>=0 ||
                    Nat::add(multiply(quotient,divisor),remainder)
                        .compare(dividend)!=0)
                    throw std::runtime_error("invalid A9 exact division identity");
                const std::size_t group=m<=9u?0u:m<=32u?1u:2u;
                out[group]+=row({std::to_string(m),std::to_string(extra),
                                 std::to_string(mode),dividend.hex(),
                                 divisor.hex(),quotient.hex(),remainder.hex()});
                ++rows[group];
            }
        }
    }
    if (rows!=std::array<std::size_t,3>{216u,144u,72u})
        throw std::runtime_error("wrong A9 source visit counts");
    return out;
}

std::array<std::string,2> b27_fixtures() {
    const auto header=[] {
        return std::string(
            "# A10 B27 fused division independent exact oracle, fixed 2 cases plus 2048 random.\n"
            "# kind repetition numerator_hex denominator_hex quotient_hex remainder_hex\n"
            "# F fixed adversarial cases, R original mt19937_64(0xB27071) visits.\n"
            "# Euclidean integer quotient/remainder, no GMP or Numeric.\n");
    };
    std::array<std::string,2> out{header(),header()};
    std::array<std::size_t,2> rows{};
    const auto add=[&](char kind,unsigned index,
                       std::vector<std::uint64_t> left,
                       std::vector<std::uint64_t> right) {
        const Nat a(std::move(left)),b(std::move(right));
        const auto [q,r]=divide(a,b);
        if (r.compare(b)>=0 || Nat::add(multiply(q,b),r).compare(a)!=0)
            throw std::runtime_error("invalid A10 exact division");
        const std::size_t group=kind=='F'||index<1024u ? 0u : 1u;
        out[group]+=row({std::string(1,kind),std::to_string(index),
                         a.hex(),b.hex(),q.hex(),r.hex()});
        ++rows[group];
    };
    const auto max=UINT64_MAX;
    add('F',0u,{max,max,max,max},
        {0xe04e3de411f558a7ULL,0xba43f747b13ac11aULL});
    add('F',1u,{0x65a601688b16b11bULL,0x5b7d0bd5828e0a0cULL,
                0x7ad9be6c73e93fdaULL,0xf370dcf1b5187984ULL,
                0x48293917315898ccULL},
        {0xc28ba8299ce45291ULL,0x91b35c306fef072dULL});
    std::mt19937_64 engine(0xB27071ULL);
    for (unsigned repetition=0u;repetition<2048u;++repetition) {
        const std::size_t n=repetition%31u==0u?96u:2u+engine()%32u;
        const std::size_t m=repetition%17u==0u?2u:2u+engine()%n;
        std::vector<std::uint64_t> left(n),right(m);
        for (auto& word:left) word=engine();
        for (auto& word:right) word=engine();
        if (repetition%13u==0u) std::fill(left.begin(),left.end(),max);
        if (repetition%19u==0u) std::fill(right.begin(),right.end(),max);
        if (repetition%23u==0u) {
            std::fill(right.begin(),right.end(),0u);
            right.back()=1u;
        }
        if (repetition%29u==0u) {
            std::fill(right.begin(),right.end(),0u);
            right.back()=std::uint64_t{1}<<63u;
        }
        add('R',repetition,std::move(left),std::move(right));
    }
    if (rows!=std::array<std::size_t,2>{1026u,1024u})
        throw std::runtime_error("wrong A10 case count");
    return out;
}

std::string b48_multiprecision_fixture() {
    std::string out=
        "# A11 B48 multiprecision division, 336 original visits.\n"
        "# count divisor_decimal mode numerator_hex quotient_hex remainder_hex\n"
        "# mt19937_64(0xB48A11D1), all original random draws and transforms.\n"
        "# Independent Euclidean arithmetic, no GMP/MPFR/Numeric/Boost.\n";
    std::mt19937_64 engine(0xB48A11D1ULL);
    const auto max=UINT64_MAX;
    std::size_t rows=0u;
    for (const std::size_t count:{2u,4u,5u,8u,16u,32u,64u,96u}) {
        for (const std::uint64_t divisor:
                {1ULL,3ULL,5ULL,7ULL,9ULL,16ULL,31ULL,
                 34ULL,63ULL,255ULL,65535ULL,4294967295ULL,
                 4294967296ULL,18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<3u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& limb:words) limb=engine();
                if (mode==0u) std::fill(words.begin(),words.end(),max);
                else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                }
                const Nat n(std::move(words)),d(divisor);
                const auto [q,r]=divide(n,d);
                if (r.compare(d)>=0 ||
                    Nat::add(multiply(q,d),r).compare(n)!=0)
                    throw std::runtime_error("invalid A11 Euclidean identity");
                out+=row({std::to_string(count),std::to_string(divisor),
                          std::to_string(mode),n.hex(),q.hex(),r.hex()});
                ++rows;
            }
        }
    }
    if (rows!=336u) throw std::runtime_error("wrong A11 visit count");
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
        const auto a5=limb_stress_fixtures();
        independent_oracle::output(dir,"big_unsigned_b48.tsv",a5[0u],write);
        independent_oracle::output(dir,"big_unsigned_b49.tsv",a5[1u],write);
        independent_oracle::output(dir,"big_unsigned_b57.tsv",a5[2u],write);
        independent_oracle::output(dir,"big_unsigned_b61.tsv",a5[3u],write);
        const auto a7=exact_limb_fixtures();
        independent_oracle::output(dir,"big_unsigned_exact_limb_small.tsv",a7[0u],write);
        independent_oracle::output(dir,"big_unsigned_exact_limb_65.tsv",a7[1u],write);
        independent_oracle::output(dir,"big_unsigned_exact_limb_129a.tsv",a7[2u],write);
        independent_oracle::output(dir,"big_unsigned_exact_limb_129b.tsv",a7[3u],write);
        const auto a8=b57_power_fixtures();
        independent_oracle::output(dir,"big_unsigned_b57_power_short.tsv",a8[0u],write);
        independent_oracle::output(dir,"big_unsigned_b57_power_medium.tsv",a8[1u],write);
        independent_oracle::output(dir,"big_unsigned_b57_power_large.tsv",a8[2u],write);
        const auto a9=b55_fixtures();
        independent_oracle::output(dir,"big_unsigned_b55_short.tsv",a9[0u],write);
        independent_oracle::output(dir,"big_unsigned_b55_medium.tsv",a9[1u],write);
        independent_oracle::output(dir,"big_unsigned_b55_large.tsv",a9[2u],write);
        const auto a10=b27_fixtures();
        independent_oracle::output(dir,"big_unsigned_b27_first.tsv",a10[0u],write);
        independent_oracle::output(dir,"big_unsigned_b27_second.tsv",a10[1u],write);
        independent_oracle::output(dir,"big_unsigned_b48_multiprecision.tsv",b48_multiprecision_fixture(),write);
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"Power/division generator: "<<error.what()<<'\n';
        return 1;
    }
}
