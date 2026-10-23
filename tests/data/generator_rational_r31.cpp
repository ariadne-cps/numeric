// SPDX-License-Identifier: MIT
// Source-reproducible R31/R32 exact rational arithmetic fixtures.
// Standalone C++20 using only independent unsigned-limb arithmetic.
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
std::size_t bit_length(Nat const& x) {
    return x.zero()?0u:(x.w.size()-1u)*64u+
        64u-static_cast<std::size_t>(std::countl_zero(x.w.back()));
}
Nat exact_divide(Nat x,Nat const& divisor) {
    if (divisor.zero()) throw std::invalid_argument("zero divisor");
    Nat result;
    while (x.compare(divisor)>=0) {
        std::size_t shift=bit_length(x)-bit_length(divisor);
        Nat selected=divisor.shl(shift);
        if (selected.compare(x)>0) selected=divisor.shl(--shift);
        x=Nat::sub(x,selected);
        result=Nat::add(result,Nat::power2(shift));
    }
    if (!x.zero()) throw std::runtime_error("non-exact rational reduction");
    return result;
}
struct Signed {
    Nat mag;
    bool neg=false;
    Signed(Nat magnitude,bool negative=false):mag(std::move(magnitude)),neg(negative&&!mag.zero()) {}
    Signed opposite() const { return Signed(mag,!neg); }
    std::string dec() const { return (neg?"-":"")+mag.dec(); }
};
Signed sum(Signed const& a,Signed const& b) {
    if (a.neg==b.neg) return Signed(Nat::add(a.mag,b.mag),a.neg);
    const int cmp=a.mag.compare(b.mag);
    if (cmp>=0) return Signed(Nat::sub(a.mag,b.mag),a.neg);
    return Signed(Nat::sub(b.mag,a.mag),b.neg);
}
Signed product(Signed const& a,Signed const& b) {
    return Signed(multiply(a.mag,b.mag),a.neg!=b.neg);
}
struct Fraction {
    Signed num;
    Nat den;
    Fraction(Signed numerator,Nat denominator):num(std::move(numerator)),den(std::move(denominator)) {
        if (den.zero()) throw std::invalid_argument("zero denominator");
        const Nat divisor=Nat::gcd(num.mag,den);
        if (divisor.compare(Nat(1u))!=0) {
            num.mag=exact_divide(num.mag,divisor);
            den=exact_divide(den,divisor);
        }
    }
};
Fraction add(Fraction const& a,Fraction const& b,bool subtract=false) {
    Signed left=product(a.num,Signed(b.den));
    Signed right=product(b.num,Signed(a.den));
    if (subtract) right=right.opposite();
    return Fraction(sum(left,right),multiply(a.den,b.den));
}
Fraction mul(Fraction const& a,Fraction const& b) {
    return Fraction(product(a.num,b.num),multiply(a.den,b.den));
}
Fraction div(Fraction const& a,Fraction const& b) {
    if (b.num.mag.zero()) throw std::invalid_argument("divide by zero fraction");
    return Fraction(Signed(multiply(a.num.mag,b.den),a.num.neg!=b.num.neg),
                    multiply(a.den,b.num.mag));
}
void append_fraction(std::vector<std::string>& fields,Fraction const& fraction) {
    fields.push_back(fraction.num.dec());
    fields.push_back(fraction.den.dec());
}
void generate(std::filesystem::path const& dir,bool write) {
    constexpr auto header=
        "# R31/R32 exact canonical source and results, independently derived using signed BigInt and Euclid.\n"
        "# bits source_index sign_mask left_num left_den right_num right_den sum_num sum_den diff_num diff_den mul_num mul_den div_num div_den\n"
        "# original index ";
    for (unsigned bits:{128u,192u,255u,256u,320u,448u,512u,1024u,2048u,4096u}) {
        const Nat base=Nat::power2(bits);
        unsigned current_chunk=~0u;
        std::string output;
        for (unsigned index=0u;index<24u;++index) {
            const unsigned chunk=bits==4096u?index/6u:bits==2048u?index/12u:0u;
            if (chunk!=current_chunk) {
                if (current_chunk!=~0u) {
                    independent_oracle::output(dir,"rational_r31_"+std::to_string(bits)+"_"+
                        std::to_string(current_chunk)+".tsv",output,write);
                }
                const unsigned start=bits==4096u?chunk*6u:bits==2048u?chunk*12u:0u;
                const unsigned stop=bits==4096u?start+5u:bits==2048u?start+11u:23u;
                output=std::string(header)+std::to_string(start)+".."+std::to_string(stop)+
                    ", all four signs, no omitted cases.\n";
                current_chunk=chunk;
            }
            const Nat an=Nat::add(base,Nat(13u+2u*index));
            const Nat ad=Nat::add(base.mul_small(2u),Nat(17u+3u*index));
            const Nat bn=Nat::add(base,Nat(31u+5u*index));
            const Nat bd=Nat::add(base.mul_small(3u),Nat(19u+7u*index));
            for (unsigned signs=0u;signs<4u;++signs) {
                const Fraction lhs(Signed(an,(signs&1u)!=0u),ad);
                const Fraction rhs(Signed(bn,(signs&2u)!=0u),bd);
                std::vector<std::string> fields={std::to_string(bits),std::to_string(index),std::to_string(signs)};
                append_fraction(fields,lhs);
                append_fraction(fields,rhs);
                append_fraction(fields,add(lhs,rhs));
                append_fraction(fields,add(lhs,rhs,true));
                append_fraction(fields,mul(lhs,rhs));
                append_fraction(fields,div(lhs,rhs));
                output+=row(fields);
            }
        }
        independent_oracle::output(dir,"rational_r31_"+std::to_string(bits)+"_"+
            std::to_string(current_chunk)+".tsv",output,write);
    }
}
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        generate(dir,write);
        return 0;
    } catch (std::exception const& e) {
        std::cerr<<"R31/R32 generator: "<<e.what()<<'\n';
        return 1;
    }
}
