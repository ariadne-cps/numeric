// SPDX-License-Identifier: MIT
// Reproduce exact source-keyed Rational R7 and R12 IEEE-754 binary64 fixtures.
// Standalone C++20: no Numeric, GMP, MPFR, FloatMP, Boost or host FP rounding.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

std::size_t bit_length(Nat const& x) {
    if (x.zero()) return 0u;
    return (x.w.size()-1u)*64u + 64u -
           static_cast<std::size_t>(std::countl_zero(x.w.back()));
}

// Bit-aligned exact division; the quotient is computed without native FP.
Nat quotient(Nat remainder, Nat const& divisor) {
    if (divisor.zero()) throw std::invalid_argument("zero denominator");
    Nat result;
    while (remainder.compare(divisor)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(divisor);
        Nat shifted=divisor.shl(shift);
        if (shifted.compare(remainder)>0) shifted=divisor.shl(--shift);
        remainder=Nat::sub(remainder,shifted);
        result=Nat::add(result,Nat::power2(shift));
    }
    return result;
}

// GMP mpq_get_d semantics: exact rational -> binary64 rounded toward zero.
// Do not consult the host's floating-point rounding mode.
std::uint64_t binary64_truncated(Nat const& numerator,Nat const& denominator,
                                 bool negative=false) {
    if (denominator.zero()) throw std::invalid_argument("zero denominator");
    const std::uint64_t sign=negative ? (std::uint64_t{1}<<63u) : 0u;
    if (numerator.zero()) return sign;

    std::int64_t exponent=static_cast<std::int64_t>(bit_length(numerator)) -
                          static_cast<std::int64_t>(bit_length(denominator));
    const bool below_power_of_two=exponent>=0
        ? numerator.compare(denominator.shl(static_cast<std::size_t>(exponent)))<0
        : numerator.shl(static_cast<std::size_t>(-exponent)).compare(denominator)<0;
    if (below_power_of_two) --exponent;
    if (exponent>1023) return sign | UINT64_C(0x7ff0000000000000);
    if (exponent< -1074) return sign;

    Nat numerator_scaled=numerator,denominator_scaled=denominator;
    if (exponent< -1022) {
        numerator_scaled=numerator.shl(1074u);
        const Nat significand=quotient(numerator_scaled,denominator_scaled);
        if (significand.w.size()>1u) throw std::runtime_error("subnormal mantissa width");
        return sign | (significand.zero()?0u:significand.w.front());
    }
    const std::int64_t scaling=52-exponent;
    if (scaling>=0)
        numerator_scaled=numerator.shl(static_cast<std::size_t>(scaling));
    else
        denominator_scaled=denominator.shl(static_cast<std::size_t>(-scaling));
    const Nat significand=quotient(numerator_scaled,denominator_scaled);
    if (significand.w.size()!=1u || significand.w.front() < (std::uint64_t{1}<<52u)
        || significand.w.front()>=(std::uint64_t{1}<<53u))
        throw std::runtime_error("normal mantissa width");
    return sign | (static_cast<std::uint64_t>(exponent+1023)<<52u)
           | (significand.w.front()-(std::uint64_t{1}<<52u));
}

struct Ratio { Nat numerator,denominator; bool negative=false; };
std::string rational_r7() {
    std::string output=
        "# R7: independent exact binary64 truncation of original rational sources.\n"
        "# mode_index bits probe_index expected_raw_uint64_bits\n"
        "# 4 host modes x 8 widths x 12 probes; expected from BigInt exact exponent+quotient arithmetic, not Rational or GMP.\n";
    for (unsigned mode=0u;mode<4u;++mode) {
        for (unsigned bits:{53u,128u,256u,448u,512u,1024u,2048u,4096u}) {
            const Nat base=Nat::power2(bits);
            const Nat twice=base.shl(1u);
            const Nat triple=base.mul_small(3u);
            const std::array<Ratio,12u> cases={{{Nat::add(base,Nat(13u)),Nat::add(twice,Nat(17u)),false},
                {Nat::sub(base,Nat(13u)),Nat::add(twice,Nat(17u)),false},
                {Nat::add(base,Nat(31u)),Nat::add(twice,Nat(19u)),false},
                {Nat::add(base,Nat(1u)),Nat::sub(base,Nat(1u)),false},
                {Nat::sub(base,Nat(1u)),Nat::add(base,Nat(1u)),false},
                {Nat::add(twice,Nat(1u)),Nat::add(base,Nat(1u)),false},
                {Nat::add(base,Nat(13u)),Nat::add(twice,Nat(17u)),true},
                {Nat::add(base,Nat(1u)),Nat(1u),false},
                {Nat::add(base,Nat(31u)),Nat::add(triple,Nat(19u)),false},
                {Nat::sub(base,Nat(1u)),Nat::add(triple,Nat(17u)),false},
                {Nat::add(base,Nat(13u)),Nat::add(triple,Nat(19u)),true},
                {Nat::add(triple,Nat(1u)),Nat::add(twice,Nat(17u)),false}}};
            for (std::size_t index=0u;index<cases.size();++index) {
                const auto& r=cases[index];
                output+=row({std::to_string(mode),std::to_string(bits),std::to_string(index),
                    std::to_string(binary64_truncated(r.numerator,r.denominator,r.negative))});
            }
        }
    }
    return output;
}
std::string rational_r12() {
    std::string output=
        "# R12 exact unsigned one-limb Rational::get_d oracle (all seeded sources).\n"
        "# mode_index sample numerator denominator expected_positive_bits expected_negative_bits\n"
        "# Four FE rounding modes, 512 visits each, both signs: 4096 bit-exact conversion outcomes.\n"
        "# Expected binary64 sign/exponent/significand computed independently using exact BigInt quotient and bit lengths, no GMP/MPFR/Numeric.\n";
    for (unsigned mode=0u;mode<4u;++mode) {
        std::mt19937_64 engine(0x1200530012812026ULL);
        for (unsigned sample=0u;sample<512u;++sample) {
            std::uint64_t n=engine() | std::uint64_t{1};
            std::uint64_t d=engine() | std::uint64_t{1};
            if (sample%4u==0u) {
                n=(std::uint64_t{1}<<53u)+engine()%100000u;
                d=2u*n+engine()%100000u;
            } else if (sample%4u==1u) {
                n=std::uint64_t{1}<<(engine()%63u);
                d=std::uint64_t{1}<<(engine()%63u);
            } else if (sample%4u==2u) {
                n=1u+engine()%1024u;
            }
            const std::uint64_t bits=binary64_truncated(Nat(n),Nat(d));
            output+=row({std::to_string(mode),std::to_string(sample),
                std::to_string(n),std::to_string(d),std::to_string(bits),
                std::to_string(bits|(std::uint64_t{1}<<63u))});
        }
    }
    return output;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"rational_r7.tsv",rational_r7(),write);
        independent_oracle::output(dir,"rational_r12.tsv",rational_r12(),write);
        return 0;
    } catch (std::exception const& e) {
        std::cerr<<"R7/R12 binary64 generator: "<<e.what()<<'\n';
        return 1;
    }
}
