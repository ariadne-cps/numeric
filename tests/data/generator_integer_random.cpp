// SPDX-License-Identifier: MIT
// Reproduce the original B17 and B18 source-keyed signed integer test data.
// Standalone C++20; no Numeric, GMP, independent exact arithmetic, Boost or independent exact arithmetic.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

Nat parse_natural(std::string const& value) {
    Nat answer;
    for(char c:value) {
        if(c<'0' || c>'9') throw std::invalid_argument("invalid decimal digit");
        answer=Nat::add(answer.mul_small(10u),Nat(static_cast<std::uint64_t>(c-'0')));
    }
    return answer;
}
struct Signed {
    Nat magnitude;
    bool negative=false;
    explicit Signed(std::string const& decimal) {
        negative=!decimal.empty() && decimal.front()=='-';
        magnitude=parse_natural(decimal.substr(negative?1u:0u));
        if(magnitude.zero()) negative=false;
    }
    Signed(Nat value,bool sign):magnitude(std::move(value)),negative(sign&&!magnitude.zero()) {}
    std::string dec() const {return (negative?"-":"")+magnitude.dec();}
};
Signed add(Signed const& a,Signed const& b) {
    if(a.negative==b.negative)
        return Signed(Nat::add(a.magnitude,b.magnitude),a.negative);
    const int cmp=a.magnitude.compare(b.magnitude);
    if(cmp>=0)return Signed(Nat::sub(a.magnitude,b.magnitude),a.negative);
    return Signed(Nat::sub(b.magnitude,a.magnitude),b.negative);
}
Signed subtract(Signed const& a,Signed const& b) {
    return add(a,Signed(b.magnitude,!b.negative));
}
Nat multiply(Nat lhs,Nat rhs) {
    Nat result;
    while(!rhs.zero()) {
        if((rhs.w.front()&1u)!=0u) result=Nat::add(result,lhs);
        rhs=rhs.shr(1u);
        if(!rhs.zero()) lhs=lhs.shl(1u);
    }
    return result;
}
Signed product(Signed const& a,Signed const& b) {
    return Signed(multiply(a.magnitude,b.magnitude),a.negative!=b.negative);
}
std::size_t bit_length(Nat const& x) {
    return x.zero()?0u:(x.w.size()-1u)*64u +
        64u-static_cast<std::size_t>(std::countl_zero(x.w.back()));
}
std::pair<Nat,Nat> unsigned_quotrem(Nat remainder,Nat const& divisor) {
    if(divisor.zero())throw std::invalid_argument("division by zero");
    Nat quotient;
    while(remainder.compare(divisor)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(divisor);
        Nat aligned=divisor.shl(shift);
        if(aligned.compare(remainder)>0)aligned=divisor.shl(--shift);
        remainder=Nat::sub(remainder,aligned);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {quotient,remainder};
}
std::pair<Signed,Signed> quotrem(Signed const& a,Signed const& b) {
    auto [q,r]=unsigned_quotrem(a.magnitude,b.magnitude);
    return {Signed(std::move(q),a.negative!=b.negative),
            Signed(std::move(r),a.negative)};
}
std::string b17() {
    std::string out=
        "# B17 independent signed arithmetic/shift oracle using exact arbitrary-precision integer calculations.\n"
        "# rep lhs rhs sum difference product trunc_quotient trunc_remainder trunc_right_shift signed_left_shift (decimal)\n"
        "# Original C++ mt19937_64(0x51C0FFEE), generated in full then sampled 1/10: rep = 10*k+(7*k+3)%10.\n"
        "# bits=rep%130; every shift count 0..129 survives, including 31/32/33, 63/64/65 and 127/128/129.\n"
        "# Negative right shift truncates toward zero; left shift preserves sign.\n";
    std::mt19937_64 engine(0x51C0FFEEULL);
    for(unsigned rep=0u;rep<2500u;++rep) {
        std::string lhs=std::to_string(engine());
        std::string rhs=std::to_string(engine()|std::uint64_t{1});
        for(unsigned j=0u;j<rep%8u;++j) {
            lhs+=std::to_string(engine());
            rhs+=std::to_string(engine());
        }
        if(rep%2u!=0u)lhs.insert(lhs.begin(),'-');
        if(rep%3u!=0u)rhs.insert(rhs.begin(),'-');
        if(rep%17u==0u)lhs="0";
        const unsigned sample=rep/10u;
        if(rep!=10u*sample+(7u*sample+3u)%10u)continue;
        const Signed a(lhs),b(rhs);
        const auto [q,r]=quotrem(a,b);
        const std::size_t bits=rep%130u;
        out+=row({std::to_string(rep),lhs,rhs,add(a,b).dec(),subtract(a,b).dec(),
            product(a,b).dec(),q.dec(),r.dec(),
            Signed(a.magnitude.shr(bits),a.negative).dec(),
            Signed(a.magnitude.shl(bits),a.negative).dec()});
    }
    return out;
}
std::string b18() {
    std::string out=
        "# B18 independent arithmetic oracle using exact arbitrary-precision signed integer calculations.\n"
        "# rep lhs rhs sum difference product trunc_quotient trunc_remainder (all decimal)\n"
        "# Original C++ mt19937_64(0x1A7E6EED), generated in full then sampled 1/10: rep = 10*k+(7*k+3)%10.\n"
        "# Includes all 7 operand-length modes, 6 sign combinations and original zero-dividend cases.\n";
    std::mt19937_64 engine(0x1A7E6EEDULL);
    for(unsigned rep=0u;rep<1700u;++rep) {
        std::string lhs=std::to_string(engine());
        std::string rhs=std::to_string(engine()|std::uint64_t{1});
        for(unsigned j=0u;j<rep%7u;++j) {
            lhs+=std::to_string(engine());
            rhs+=std::to_string(engine());
        }
        if(rep%2u!=0u)lhs.insert(lhs.begin(),'-');
        if(rep%3u!=0u)rhs.insert(rhs.begin(),'-');
        if(rep%19u==0u)lhs="0";
        const unsigned sample=rep/10u;
        if(rep!=10u*sample+(7u*sample+3u)%10u)continue;
        const Signed a(lhs),b(rhs);
        const auto [q,r]=quotrem(a,b);
        out+=row({std::to_string(rep),lhs,rhs,add(a,b).dec(),subtract(a,b).dec(),
            product(a,b).dec(),q.dec(),r.dec()});
    }
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"integer_b17.tsv",b17(),write);
        independent_oracle::output(dir,"integer_b18.tsv",b18(),write);
        return 0;
    } catch(std::exception const& error) {
        std::cerr<<"Integer B17/B18 generator: "<<error.what()<<'\n';
        return 1;
    }
}
