// SPDX-License-Identifier: MIT
// Independent, reproducible R28, R30 and R22 rational golden data.
// Standalone C++20: standard library only; no Numeric, GMP, MPFR or Boost.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

std::size_t bit_length(Nat const& n) {
    return n.zero() ? 0u : (n.w.size()-1u)*64u +
        64u-static_cast<std::size_t>(std::countl_zero(n.w.back()));
}
std::pair<Nat,Nat> divide(Nat remainder, Nat const& divisor) {
    if(divisor.zero()) throw std::invalid_argument("division by zero");
    Nat quotient;
    while(remainder.compare(divisor)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(divisor);
        Nat aligned=divisor.shl(shift);
        if(aligned.compare(remainder)>0) aligned=divisor.shl(--shift);
        remainder=Nat::sub(remainder,aligned);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {quotient,remainder};
}
Nat multiply(Nat left, Nat right) {
    Nat result;
    while(!right.zero()) {
        if((right.w.front()&1u)!=0u) result=Nat::add(result,left);
        right=right.shr(1u);
        if(!right.zero()) left=left.shl(1u);
    }
    return result;
}
struct Signed {
    Nat mag;
    bool negative=false;
    Signed(Nat x,bool sign=false) : mag(std::move(x)),negative(sign&&!mag.zero()) {}
    Signed opposite() const { return Signed(mag,!negative); }
    std::string dec() const { return (negative?"-":"")+mag.dec(); }
};
Signed sum(Signed const& a,Signed const& b) {
    if(a.negative==b.negative) return Signed(Nat::add(a.mag,b.mag),a.negative);
    const int cmp=a.mag.compare(b.mag);
    if(cmp>=0) return Signed(Nat::sub(a.mag,b.mag),a.negative);
    return Signed(Nat::sub(b.mag,a.mag),b.negative);
}
Signed product(Signed const& a,Signed const& b) {
    return Signed(multiply(a.mag,b.mag),a.negative!=b.negative);
}
struct Fraction {
    Signed num;
    Nat den;
    Fraction(Signed n,Nat d) : num(std::move(n)),den(std::move(d)) {
        if(den.zero()) throw std::invalid_argument("zero fraction denominator");
        const Nat g=Nat::gcd(num.mag,den);
        if(g.compare(Nat(1u))!=0) {
            auto [q,r]=divide(num.mag,g);
            if(!r.zero()) throw std::runtime_error("non-exact numerator reduction");
            num.mag=std::move(q);
            auto [dq,dr]=divide(den,g);
            if(!dr.zero()) throw std::runtime_error("non-exact denominator reduction");
            den=std::move(dq);
        }
        if(num.mag.zero()) num.negative=false;
    }
    std::string str() const { return num.dec()+"/"+den.dec(); }
};
Fraction add(Fraction const& a,Fraction const& b,bool subtract=false) {
    Signed na=product(a.num,Signed(b.den));
    Signed nb=product(b.num,Signed(a.den));
    if(subtract) nb=nb.opposite();
    return Fraction(sum(na,nb),multiply(a.den,b.den));
}
Fraction mul(Fraction const& a,Fraction const& b) {
    return Fraction(product(a.num,b.num),multiply(a.den,b.den));
}
Fraction quotient(Fraction const& a,Fraction const& b) {
    if(b.num.mag.zero()) throw std::invalid_argument("division by zero fraction");
    const bool neg=a.num.negative!=b.num.negative;
    return Fraction(Signed(multiply(a.num.mag,b.den),neg),multiply(a.den,b.num.mag));
}
std::string r28() {
    std::string out=
        "# R28 independent exact signed add/sub rational oracle at original heap thresholds.\n"
        "# bits mode a_num a_den b_num b_den sum_num sum_den diff_num diff_den\n"
        "# 9 original precision boundaries, all 4 sign combinations. Canonical fractions derived from exact signed integers.\n";
    for(unsigned bits:{127u,128u,255u,256u,257u,319u,320u,448u,512u}) {
        const Nat base=Nat::power2(bits);
        for(unsigned mode=0u;mode<4u;++mode) {
            const Nat an=Nat::add(base,Nat(13u)),ad=Nat::add(base.mul_small(2u),Nat(17u));
            const Nat bn=Nat::add(base,Nat(31u)),bd=Nat::add(base.mul_small(3u),Nat(19u));
            const Signed raw_a(an,(mode&1u)!=0u),raw_b(bn,(mode&2u)!=0u);
            const Fraction a(raw_a,ad);
            const Fraction b(raw_b,bd);
            const Fraction s=add(a,b),d=add(a,b,true);
            out+=row({std::to_string(bits),std::to_string(mode),raw_a.dec(),ad.dec(),
                     raw_b.dec(),bd.dec(),s.num.dec(),s.den.dec(),d.num.dec(),d.den.dec()});
        }
    }
    return out;
}
std::string r30() {
    std::string out=
        "# R30 independently derived exact canonical signed fractions, including sparse and dense fallbacks.\n"
        "# bits high_pair low_pair signs expected_num expected_den\n"
        "# 10 widths x 7 ratios x 6 tails x 4 signs plus 9 dense fallbacks.\n";
    constexpr std::array<std::pair<unsigned,unsigned>,7u> high_pairs={{{1u,2u},{1u,3u},
       {2u,4u},{3u,9u},{3u,3u},{2u,3u},{1u,1u}}};
    constexpr std::array<std::pair<std::uint64_t,std::uint64_t>,6u> low_pairs={{{13u,17u},
       {63u,126u},{1u,0u},{7u,14u},{0u,0u},{UINT64_MAX,UINT64_MAX}}};
    for(unsigned bits:{64u,128u,255u,256u,320u,512u,1024u,2048u,4096u,8192u}) {
        const Nat base=Nat::power2(bits);
        for(unsigned h=0u;h<high_pairs.size();++h) {
            auto const& [hn,hd]=high_pairs[h];
            for(unsigned l=0u;l<low_pairs.size();++l) {
                auto const& [ln,ld]=low_pairs[l];
                const Nat n=Nat::add(base.mul_small(hn),Nat(ln));
                const Nat d=Nat::add(base.mul_small(hd),Nat(ld));
                for(unsigned signs=0u;signs<4u;++signs) {
                    const Fraction x(Signed(n,((signs&1u)!=0u)!=((signs&2u)!=0u)),d);
                    out+=row({std::to_string(bits),std::to_string(h),std::to_string(l),
                              std::to_string(signs),x.num.dec(),x.den.dec()});
                }
            }
        }
        if(bits>=128u) {
            const Nat mid=Nat::power2(64u);
            const Nat n=Nat::add(Nat::add(base,mid),Nat(11u));
            const Nat d=Nat::add(Nat::add(base.mul_small(2u),mid),Nat(17u));
            const Fraction dense(Signed(n),d);
            out+=row({std::to_string(bits),"7","0","0",dense.num.dec(),dense.den.dec()});
        }
    }
    return out;
}
std::string r22() {
    std::string out=
        "# Rational R22 original-source exact scalar cancellation oracle.\n"
        "# width_bits factor_index sign_mask canonical_lhs canonical_rhs product lhs_div_rhs\n"
        "# Every original 10 x 4 x 4 case retained: 53/127/128/255/256/447/448/511/512/1024 bits; scalar factors 1,1/15,14,14/15.\n"
        "# Expected canonical rational references computed by independent exact signed integer cross products and Euclidean gcd.\n";
    const std::array<Fraction,4u> factors={
        Fraction(Signed(Nat(1u)),Nat(1u)),Fraction(Signed(Nat(1u)),Nat(15u)),
        Fraction(Signed(Nat(14u)),Nat(1u)),Fraction(Signed(Nat(14u)),Nat(15u))};
    for(unsigned bits:{53u,127u,128u,255u,256u,447u,448u,511u,512u,1024u}) {
        const Nat base=Nat::power2(bits);
        const Nat n=Nat::add(base,Nat(1u)).mul_small(15u);
        const Nat d=base.mul_small(14u);
        for(unsigned index=0u;index<factors.size();++index) {
            for(unsigned mask=0u;mask<4u;++mask) {
                const Fraction a(Signed(n,(mask&1u)!=0u),d);
                const Fraction b(Signed(factors[index].num.mag,(mask&2u)!=0u),factors[index].den);
                const Fraction m=mul(a,b),q=quotient(a,b);
                out+=row({std::to_string(bits),std::to_string(index),std::to_string(mask),
                          a.str(),b.str(),m.str(),q.str()});
            }
        }
    }
    return out;
}

std::string r16() {
    std::string out=
        "# R16 exact independent Rational add/sub and alias results at all 6 widths.\n"
        "# bits pair_index order lhs rhs sum difference twice_lhs zero_lhs\n"
        "# 7 coprime/shared-denominator structural classes and both orders; original 84 cases preserved.\n";
    auto fraction=[](Nat n,Nat d,bool negative=false) {
        return Fraction(Signed(std::move(n),negative),std::move(d));
    };
    for (unsigned bits:{53u,128u,256u,1024u,2048u,4096u}) {
        const Nat base=Nat::power2(bits);
        const Nat next=Nat::add(base,Nat(1u));
        const std::array<std::pair<Fraction,Fraction>,7u> cases={{
            {fraction(Nat(1u),base),fraction(Nat(1u),next)},
            {fraction(Nat(7u),base),fraction(Nat(11u),next,true)},
            {fraction(Nat(1u),base,true),fraction(Nat(1u),next)},
            {fraction(Nat(1u),base),fraction(Nat(1u),base,true)},
            {fraction(Nat(1u),base.mul_small(2u)),fraction(Nat(1u),next.mul_small(2u))},
            {fraction(Nat(1u),base.mul_small(2u),true),fraction(Nat(1u),next.mul_small(2u))},
            {fraction(Nat(1u),Nat(6u)),fraction(Nat(1u),Nat(10u))}
        }};
        for (unsigned pair=0u;pair<cases.size();++pair) {
            for (unsigned order=0u;order<2u;++order) {
                const Fraction& a=order==0u?cases[pair].first:cases[pair].second;
                const Fraction& b=order==0u?cases[pair].second:cases[pair].first;
                out+=row({std::to_string(bits),std::to_string(pair),std::to_string(order),
                    a.str(),b.str(),add(a,b).str(),add(a,b,true).str(),
                    add(a,a).str(),add(a,a,true).str()});
            }
        }
    }
    return out;
}
std::string r19() {
    std::string out=
        "# R19 signed cross-GCD multiplication, division and alias fixture.\n"
        "# width factor_index sign_mask canonical_lhs canonical_rhs product lhs_div_rhs rhs_div_lhs\n"
        "# All original 96 sign/factor cases plus 6 self-alias cases, exact independent integers.\n";
    for (unsigned bits:{53u,128u,256u,1024u,2048u,4096u}) {
        const Nat base=Nat::power2(bits);
        const Nat p3=Nat::add(base,Nat(3u));
        const Nat p5=Nat::add(base,Nat(5u));
        const Nat p7=Nat::add(base,Nat(7u));
        const Nat p9=Nat::add(base,Nat(9u));
        const Fraction original(Signed(p3),p5);
        const std::array<Fraction,5u> factors={{
            Fraction(Signed(p9),p7),Fraction(Signed(p9),p3),
            Fraction(Signed(p5),p9),Fraction(Signed(p5),p3),original
        }};
        for (unsigned factor=0u;factor<factors.size();++factor) {
            for (unsigned mask=0u;mask<(factor==4u?1u:4u);++mask) {
                const Fraction lhs(Signed(original.num.mag,(mask&1u)!=0u),original.den);
                const Fraction rhs(Signed(factors[factor].num.mag,(mask&2u)!=0u),
                                   factors[factor].den);
                out+=row({std::to_string(bits),std::to_string(factor),std::to_string(mask),
                    lhs.str(),rhs.str(),mul(lhs,rhs).str(),
                    quotient(lhs,rhs).str(),quotient(rhs,lhs).str()});
            }
        }
    }
    return out;
}
std::string r20() {
    std::string out=
        "# R20 independently exact signed cross-sum/difference oracle, no GMP at test time.\n"
        "# width pair_index sign_mask canonical_lhs canonical_rhs exact_sum exact_difference\n"
        "# Preserve all 10 width thresholds, 4 structural operand pairs and 4 sign combinations.\n"
        "# Fraction tokens use canonical num/den, calculated by exact integer arithmetic.\n";
    for (unsigned bits:{53u,128u,255u,256u,320u,447u,448u,511u,512u,1024u}) {
        const Nat base=Nat::power2(bits);
        const Fraction first(Signed(Nat::add(base,Nat(13u))),
                             Nat::add(base.mul_small(2u),Nat(17u)));
        const std::array<std::pair<Fraction,Fraction>,4u> pairs={{
            {first,Fraction(Signed(Nat::add(base,Nat(31u))),
                            Nat::add(base.mul_small(3u),Nat(19u)))},
            {Fraction(Signed(Nat::sub(base,Nat(13u))),
                      Nat::add(base.mul_small(3u),Nat(23u))),
             Fraction(Signed(Nat::add(base,Nat(19u))),
                      Nat::add(base.mul_small(2u),Nat(31u)))},
            {Fraction(Signed(Nat::add(base,Nat(1u))),Nat::add(base,Nat(2u))),
             Fraction(Signed(Nat::add(base,Nat(2u))),Nat::add(base,Nat(3u)))},
            {first,first}
        }};
        for (unsigned pair=0u;pair<pairs.size();++pair) {
            for (unsigned mask=0u;mask<4u;++mask) {
                const Fraction& raw_a=pairs[pair].first;
                const Fraction& raw_b=pairs[pair].second;
                const Fraction a(Signed(raw_a.num.mag,(mask&1u)!=0u),raw_a.den);
                const Fraction b(Signed(raw_b.num.mag,(mask&2u)!=0u),raw_b.den);
                out+=row({std::to_string(bits),std::to_string(pair),std::to_string(mask),
                    a.str(),b.str(),add(a,b).str(),add(a,b,true).str()});
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
        independent_oracle::output(dir,"rational_r28.tsv",r28(),write);
        independent_oracle::output(dir,"rational_r30.tsv",r30(),write);
        independent_oracle::output(dir,"rational_r22.tsv",r22(),write);
        independent_oracle::output(dir,"rational_r16.tsv",r16(),write);
        independent_oracle::output(dir,"rational_r19.tsv",r19(),write);
        independent_oracle::output(dir,"rational_r20.tsv",r20(),write);
        return 0;
    } catch(std::exception const& error) {
        std::cerr<<"Rational exact generator: "<<error.what()<<'\n';
        return 1;
    }
}
