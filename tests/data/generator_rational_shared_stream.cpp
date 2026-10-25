// SPDX-License-Identifier: MIT
// Reproduce the R29, R13 and R5 exact binary64 fixtures from their shared source stream.
// Standalone C++20: no Numeric, GMP, independent exact arithmetic, Boost, independent exact arithmetic or native FP arithmetic.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;
using Words=std::vector<std::uint64_t>;

std::size_t bit_length(Nat const& x) {
    if (x.zero()) return 0u;
    return (x.w.size()-1u)*64u+64u-static_cast<std::size_t>(std::countl_zero(x.w.back()));
}
Nat divide(Nat remainder,Nat const& denominator) {
    if (denominator.zero()) throw std::invalid_argument("zero denominator");
    Nat result;
    while (remainder.compare(denominator)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(denominator);
        Nat aligned=denominator.shl(shift);
        if (aligned.compare(remainder)>0) aligned=denominator.shl(--shift);
        remainder=Nat::sub(remainder,aligned);
        result=Nat::add(result,Nat::power2(shift));
    }
    return result;
}
std::uint64_t binary64_truncated(Nat const& numerator,Nat const& denominator,bool negative=false) {
    if (denominator.zero()) throw std::invalid_argument("zero denominator");
    if (numerator.zero()) return 0u;
    const std::uint64_t sign=negative?(std::uint64_t{1}<<63u):0u;
    std::int64_t exponent=static_cast<std::int64_t>(bit_length(numerator))-
                          static_cast<std::int64_t>(bit_length(denominator));
    const bool below=exponent>=0
        ? numerator.compare(denominator.shl(static_cast<std::size_t>(exponent)))<0
        : numerator.shl(static_cast<std::size_t>(-exponent)).compare(denominator)<0;
    if (below) --exponent;
    if (exponent>1023) return sign|UINT64_C(0x7ff0000000000000);
    if (exponent< -1074) return 0u;
    Nat top=numerator,bottom=denominator;
    if (exponent< -1022) {
        const Nat significand=divide(numerator.shl(1074u),denominator);
        if (significand.w.size()>1u) throw std::runtime_error("subnormal range");
        if (significand.zero()) return 0u;
        return sign|significand.w.front();
    }
    const std::int64_t shift=52-exponent;
    if (shift>=0) top=numerator.shl(static_cast<std::size_t>(shift));
    else bottom=denominator.shl(static_cast<std::size_t>(-shift));
    const Nat significand=divide(top,bottom);
    if (significand.w.size()!=1u || significand.w.front()<(std::uint64_t{1}<<52u)
        || significand.w.front()>=(std::uint64_t{1}<<53u))
        throw std::runtime_error("normal significand range");
    return sign|(static_cast<std::uint64_t>(exponent+1023)<<52u)
        |(significand.w.front()-(std::uint64_t{1}<<52u));
}
std::mt19937_64 initialize_shared_source() {
    std::mt19937_64 engine(0x6672616374696f6eULL);
    // R0 precedes R29/R13/R5 in the original test, and consumes
    // all 1200 visits regardless of the 120 saved samples.
    for (unsigned iteration=0u;iteration<1200u;++iteration) {
        const std::size_t words=1u+iteration%17u;
        for (unsigned operand=0u;operand<4u;++operand) {
            for (std::size_t j=0u;j<words;++j) (void)engine();
            (void)engine(); // original sign draw
        }
    }
    return engine;
}
std::string generate_r29_mode(std::mt19937_64& engine,unsigned mode) {
    std::string out=
        "# R29 original exact mt19937_64 sparse prefix source operands and IEEE-754 binary64 results\n"
        "# mode_index bits numerator_high denominator_high source_num source_den positive_uint64_bits negative_uint64_bits\n"
        "# All 4 rounding modes * 6 widths * 3 numerator highs * 3 denominator highs, plus both signs; no sampling\n";
    if (mode!=0u) out.clear();
    for (unsigned bits:{128u,256u,320u,512u,1024u,2048u}) {
            const Nat base=Nat::power2(bits);
            for (std::uint64_t high_n:{1u,3u,7u})
                for (std::uint64_t high_d:{2u,3u,5u}) {
                    const Nat n=Nat::add(base.mul_small(high_n),Nat((engine()&65535u)+1u));
                    const Nat d=Nat::add(base.mul_small(high_d),Nat((engine()&65535u)+1u));
                    const std::uint64_t positive=binary64_truncated(n,d);
                    out+=row({std::to_string(mode),std::to_string(bits),std::to_string(high_n),
                        std::to_string(high_d),n.dec(),d.dec(),std::to_string(positive),
                        std::to_string(binary64_truncated(n,d,true))});
                }
        }
    return out;
}
std::string generate_r13(std::mt19937_64& engine,unsigned mode,unsigned limbs) {
    std::string out=
        "# R13 independent source-keyed native conversion oracle, original seeded random limb streams\n"
        "# mode_index limbs sample raw_signed_numerator raw_denominator expected_binary64_bits\n"
        "# Signed numerator and denominator are exact from C++ std::mt19937_64 limb vectors; no GMP/independent exact arithmetic/Numeric\n";
    for (unsigned sample=0u;sample<96u;++sample) {
        Words nw(limbs),dw(limbs);
        for (unsigned j=0u;j<limbs;++j) {
            nw[j]=engine();dw[j]=engine();
        }
        nw.back()=(engine()&3u)+1u;
        dw.back()=(engine()&3u)+2u;
        const Nat n(std::move(nw)),d(std::move(dw));
        const bool negative=sample%2u!=0u;
        out+=row({std::to_string(mode),std::to_string(limbs),std::to_string(sample),
            (negative?"-":"")+n.dec(),d.dec(),
            std::to_string(binary64_truncated(n,d,negative))});
    }
    return out;
}
void run_all(std::filesystem::path const& dir,bool write) {
    auto engine=initialize_shared_source();
    std::string r29data;
    std::string r5data=
        "# R5: exact frozen IEEE-754 binary64 truncation oracle for original random shift probes.\n"
        "# mode_index shift original_num original_den pos_div_bits neg_div_bits pos_mul_bits neg_mul_bits\n"
        "# 4 modes x 12 shifts x 4 signed quotient cases: all 192 original results, no thinning.\n"
        "# All source mt19937_64 words and expected exact binary64 patterns independently derived, no GMP/independent exact arithmetic/Numeric.\n";
    for (unsigned mode=0u;mode<4u;++mode) {
        r29data+=generate_r29_mode(engine,mode);
        for (unsigned limbs=5u;limbs<=8u;++limbs) {
            const std::string filename="rational_r13_mode"+std::to_string(mode)
                +"_limbs"+std::to_string(limbs)+".tsv";
            independent_oracle::output(dir,filename,generate_r13(engine,mode,limbs),write);
        }
        for (unsigned shift:{0u,1u,52u,53u,1022u,1023u,1024u,1073u,1074u,1075u,1076u,2048u}) {
            const Nat numerator(engine()|std::uint64_t{1});
            const Nat denominator(engine()|std::uint64_t{1});
            const Nat n_div=numerator,d_div=denominator.shl(shift);
            const Nat n_mul=numerator.shl(shift),d_mul=denominator;
            const auto pos_div=binary64_truncated(n_div,d_div);
            const auto pos_mul=binary64_truncated(n_mul,d_mul);
            r5data+=row({std::to_string(mode),std::to_string(shift),numerator.dec(),
                denominator.dec(),std::to_string(pos_div),
                std::to_string(binary64_truncated(n_div,d_div,true)),
                std::to_string(pos_mul),std::to_string(binary64_truncated(n_mul,d_mul,true))});
        }
    }
    independent_oracle::output(dir,"rational_r29.tsv",r29data,write);
    independent_oracle::output(dir,"rational_r5.tsv",r5data,write);
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        run_all(dir,write);
        return 0;
    } catch (std::exception const& e) {
        std::cerr<<"R29/R13/R5 shared-stream generator: "<<e.what()<<'\n';
        return 1;
    }
}
