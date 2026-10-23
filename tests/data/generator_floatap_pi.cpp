// SPDX-License-Identifier: MIT
// Rigorous integer-only Machin-formula pi oracle, independent C++20.
// No GMP, MPFR, FloatAP, Numeric, Boost or native floating point.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

// Binary long division by a small odd positive integer (always <= 57121).
Nat div_small(Nat const& n,std::uint64_t denominator) {
    if (denominator==0u)throw std::invalid_argument("zero divisor");
    std::vector<std::uint64_t> words(n.w.size(),0u);
    std::uint64_t rem=0u;
    for(std::size_t word=n.w.size();word>0;--word) {
        const std::uint64_t value=n.w[word-1u];
        for(unsigned bit=64u;bit>0u;--bit) {
            rem=(rem<<1u)|((value>>(bit-1u))&1u);
            if(rem>=denominator) {
                rem-=denominator;
                words[word-1u]|=std::uint64_t{1}<<(bit-1u);
            }
        }
    }
    return Nat(std::move(words));
}
struct Interval {Nat lower,upper;};
Interval atan_inverse(unsigned denominator) {
    // A=2^1200. Odd denominators ensure every positive nonzero series
    // term has strictly positive fractional part on this dyadic grid.
    Nat quotient=div_small(Nat::power2(1200u),denominator);
    Nat lower,upper;
    for(unsigned k=0u;k<400u;++k) {
        const Nat floor_term=div_small(quotient,2u*k+1u);
        if(k%2u==0u) {
            lower=Nat::add(lower,floor_term);
            upper=Nat::add(upper,Nat::add(floor_term,Nat(1u)));
        } else {
            lower=Nat::sub(lower,Nat::add(floor_term,Nat(1u)));
            upper=Nat::sub(upper,floor_term);
        }
        quotient=div_small(quotient,static_cast<std::uint64_t>(denominator)*denominator);
    }
    // The 400-term alternating sum is a strict LOWER approximation:
    // the 401st term is positive and < 2^-1200 (at these denominators).
    // Expand the upper endpoint by one fixed-point unit as a tail bound.
    return {std::move(lower),Nat::add(upper,Nat(1u))};
}
Nat rounded_quotient(Nat const& value,unsigned shift,bool nearest) {
    Nat truncated=value.shr(shift);
    const Nat remainder=Nat::sub(value,truncated.shl(shift));
    if(nearest) {
        const int comparison=remainder.shl(1u).compare(Nat::power2(shift));
        if(comparison>0 || (comparison==0 && !truncated.zero() &&
                              (truncated.w.front()&1u)!=0u))
            truncated=Nat::add(truncated,Nat(1u));
    } else if(!remainder.zero())truncated=Nat::add(truncated,Nat(1u));
    return truncated;
}
std::pair<std::string,std::string> fraction(Nat numerator,unsigned denominator_bits) {
    if(numerator.zero())return {"0","1"};
    const std::size_t power=std::min(static_cast<std::size_t>(denominator_bits),
                                      numerator.trailing_zero_bits());
    numerator=numerator.shr(power);
    return {numerator.dec(),Nat::power2(denominator_bits-power).dec()};
}
std::string generate() {
    std::string out=
        "# Frozen pi oracle: rigorous 384-bit bracketing and nearest-even 320-bit value.\n"
        "# precision_lower_upper_bits nearest_bits lower_num lower_den upper_num upper_den near_num near_den\n"
        "# Computed by signed-integer fixed-point arithmetic, B=2^1200:\n"
        "# pi=16*atan(1/5)-4*atan(1/239), alternating Taylor partial sums\n"
        "# for k=0..399, per-term outward division and explicit next-term bound.\n"
        "# Independent of FloatAP, FloatMP, MPFR, GMP and floating-point pi.\n"
        "# Both fixed-point interval endpoints certify identical binary directed\n"
        "# 384-bit bounds and identical nearest-even 320-bit reference.\n";
    const Interval a=atan_inverse(5u),b=atan_inverse(239u);
    const Nat low=Nat::sub(a.lower.mul_small(16u),b.upper.mul_small(4u));
    const Nat high=Nat::sub(a.upper.mul_small(16u),b.lower.mul_small(4u));
    constexpr unsigned floor_shift=1200u-384u+2u;
    constexpr unsigned nearest_shift=1200u-320u+2u;
    // Require the ENTIRE certified pi interval to fit into one dyadic cell.
    const Nat lower_cell=low.shr(floor_shift),upper_cell=high.shr(floor_shift);
    if(lower_cell.compare(upper_cell)!=0)
        throw std::runtime_error("pi interval does not certify 384-bit floor");
    const Nat lower_nearest=rounded_quotient(low,nearest_shift,true);
    const Nat upper_nearest=rounded_quotient(high,nearest_shift,true);
    if(lower_nearest.compare(upper_nearest)!=0)
        throw std::runtime_error("pi interval does not certify 320-bit nearest");
    const auto [l_num,l_den]=fraction(lower_cell,382u);
    const auto [u_num,u_den]=fraction(rounded_quotient(high,floor_shift,false),382u);
    const auto [n_num,n_den]=fraction(lower_nearest,318u);
    out+=row({"384","320",l_num,l_den,u_num,u_den,n_num,n_den});
    return out;
}
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatap_pi.tsv",generate(),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"pi generator: "<<e.what()<<'\n';
        return 1;
    }
}
