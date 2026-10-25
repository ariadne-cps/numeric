// SPDX-License-Identifier: MIT
// Reproducible original R6/R8/R10/R11/R21/R26/R27 GCD golden fixtures.
// Standalone portable C++20, exact unsigned limbs; no Numeric/GMP/independent exact arithmetic/Boost.
#include "oracle_exact_integer.hpp"
using independent_oracle::Nat;
using independent_oracle::row;
namespace {
constexpr std::uint64_t MAX=UINT64_MAX;
std::string s(unsigned long long n) {return std::to_string(n);}
std::string h(Nat const& n) {return n.hex();}
std::string g(Nat const& x,Nat const& y) {return Nat::gcd(x,y).hex();}
std::string r6() {
    std::string result="# R6: independent BigInt Euclid GCD, original 7 widths x 4 offsets x two operand orders\n"
        "# bits offset order lhs_hex rhs_hex expected_gcd_hex\n";
    for(unsigned bits:{128u,192u,256u,257u,320u,512u,1024u})
        for(std::uint64_t offset:{1u,3u,17u,65537u}) {
            const Nat base=Nat::power2(bits);
            const Nat left=Nat::add(base,Nat(offset*3u));
            const Nat right=Nat::add(base,Nat(offset*5u));
            result+=row({s(bits),s(offset),"0",h(left),h(right),g(left,right)});
            result+=row({s(bits),s(offset),"1",h(right),h(left),g(right,left)});
        }
    return result;
}
std::string r8() {
    std::string result="# R8: original seeded 640 single-word GCD cases and 36 wide/narrow cases\n"
        "# kind index width_or_divisor lhs_hex rhs_hex gcd_hex ; S: index=i,width_or_divisor=0; W: index=bits,width_or_divisor=d\n";
    std::mt19937_64 random(0x8890ABCDEF234ULL);
    for(unsigned i=0;i<640;++i) {
        const std::uint64_t x=i%11u==0u?0u:random();
        const std::uint64_t y=i%17u==0u?0u:random();
        const Nat a(x),b(y);
        result+=row({"S",s(i),"0",h(a),h(b),g(a,b)});
    }
    for(unsigned bits:{128u,255u,256u,257u,1024u,4096u}) {
        const Nat wide=Nat::add(Nat::power2(bits),Nat(77));
        for(std::uint64_t d:std::array<std::uint64_t,6>{1u,2u,3u,17u,65537u,UINT64_MAX}) {
            const Nat small(d);
            result+=row({"W",s(bits),s(d),h(wide),h(small),g(wide,small)});
        }
    }
    return result;
}
std::string r10() {
    std::string result="# R10: 9 widths x 5 ordered operand pairs plus the original self-alias probe, independent BigInt GCD.\n"
        "# bits case_index lhs_hex rhs_hex expected_gcd_hex; case 5 = gcd(close,close)\n";
    for(unsigned bits:{53u,128u,255u,256u,257u,512u,1024u,2048u,4096u}) {
        const Nat base=Nat::power2(bits);
        const Nat close=Nat::add(base,Nat(41));
        const Nat apart=Nat::add(base.shl(3),Nat(93));
        const Nat multiple=close.mul_small(7);
        const std::vector<std::pair<Nat,Nat>> pairs={
            {close,apart},{apart,close},{close,multiple},{apart,Nat(65537)},
            {multiple,Nat()},{close,close}
        };
        for(std::size_t i=0;i<pairs.size();++i) {
            auto const& [a,b]=pairs[i];
            result+=row({s(bits),s(i),h(a),h(b),g(a,b)});
        }
    }
    return result;
}
std::string r11() {
    std::string result="# R11: independent Euclidean gcd at q=1/2/3 transitions and 4/5-limb boundary.\n"
        "# bits dividend_index lhs_hex rhs_hex expected_gcd_hex\n";
    for(unsigned bits:{53u,128u,192u,255u,256u,257u,320u,512u,1024u,2048u,4096u}) {
        const Nat unit=Nat::power2(bits);
        const Nat divisor=Nat::add(unit,Nat(5));
        const Nat twice=unit.mul_small(2),thrice=unit.mul_small(3);
        const std::vector<Nat> inputs={
            Nat::add(twice,Nat(7)),Nat::add(twice,Nat(10)),
            Nat::add(twice,Nat(17)),Nat::add(thrice,Nat(23)),divisor
        };
        for(std::size_t i=0;i<inputs.size();++i)
            result+=row({s(bits),s(i),h(inputs[i]),h(divisor),g(inputs[i],divisor)});
    }
    return result;
}
std::string r26(std::size_t count,std::mt19937_64& random) {
    std::string result="# R26 40 original trial inputs per width; 9 original q/variant GCD results per trial.\n"
        "# count trial source_low_hex source_b_hex q2v0 q2v1 q2v2 q3v0 q3v1 q3v2 q4v0 q4v1 q4v2\n"
        "# Frozen results from independent exact unsigned BigInt Euclid, including 127/128/129-limb boundary.\n";
    const Nat top=Nat::power2(64*(count-1));
    for(unsigned trial=0;trial<40;++trial) {
        const Nat low(random()|1ULL);
        const Nat b=Nat::add(top,low);
        std::vector<std::string> fields={s(count),s(trial),h(low),h(b)};
        for(std::uint64_t q:{2ULL,3ULL,4ULL})for(unsigned variant=0;variant<3;++variant) {
            const Nat p=b.mul_small(q);
            const Nat a=variant==0?p:variant==1?Nat::add(p,Nat(1)):
                Nat::sub(Nat::add(p,b),Nat(1));
            fields.push_back(g(a,b));
        }
        result+=row(fields);
    }
    return result;
}
std::string r27() {
    std::string result="# R27 source-keyed original 7 * 96 shared-prefix plus 8 * 96 scalar-tail GCD cases.\n"
        "# group count repetition lhs_hex rhs_hex expected_gcd_hex (P=prefix,S=scalar)\n"
        "# Full original std::mt19937_64 stream, frozen source operands and exact independent Euclidean GCD.\n";
    std::mt19937_64 random(0x5277705245464958ULL);
    for(std::size_t count:{2u,3u,4u,5u,7u,8u,9u})for(unsigned rep=0;rep<96;++rep) {
        std::vector<std::uint64_t> left(count),right(count);
        for(auto& word:left) word=random();
        left.back()|=1ULL<<63;
        right=left;
        const std::uint64_t diff=rep%5==0?0ULL:rep%5==1?1ULL:
            rep%5==2?6ULL:rep%5==3?65537ULL:1ULL<<63;
        left.front()=MAX;
        right.front()=MAX-diff;
        const Nat a(std::move(left)),b(std::move(right));
        result+=row({"P",s(count),s(rep),h(a),h(b),g(a,b)});
    }
    for(std::size_t count:{1u,2u,3u,4u,5u,7u,8u,9u})for(unsigned rep=0;rep<96;++rep) {
        std::vector<std::uint64_t> limbs(count);
        for(auto& word:limbs) word=random();
        limbs.back()|=1ULL<<63;
        const std::uint64_t d=rep%6==0?1ULL:rep%6==1?2ULL:
            rep%6==2?65537ULL:rep%6==3?1ULL<<63:
            rep%6==4?MAX:random()|1ULL;
        const Nat a(std::move(limbs)),b(d);
        result+=row({"S",s(count),s(rep),h(a),h(b),g(a,b)});
    }
    return result;
}
std::string r21() {
    std::string result="# R21 independent 40 exact dense GCD oracles, source-keyed by width and trial.\n"
        "# count trial lhs_hex rhs_hex gcd_hex ; use original mt19937_64 stream\n"
        "# Fixed gcd factor cases and their aliasing tests remain unchanged.\n";
    std::mt19937_64 random(0x21BCD128ULL);
    for(std::size_t count:{2u,4u,5u,8u,16u,64u,127u,128u,129u,160u})for(unsigned trial=0;trial<4;++trial){
        std::vector<std::uint64_t> x(count),y(count);
        for(auto& word:x) word=random();
        for(auto& word:y) word=random();
        x.back()|=1ULL<<63;
        y.back()|=1ULL<<63;
        const Nat a(std::move(x)),b(std::move(y));
        result+=row({s(count),s(trial),h(a),h(b),g(a,b)});
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"big_unsigned_r6.tsv",r6(),write);
        independent_oracle::output(dir,"big_unsigned_r8.tsv",r8(),write);
        independent_oracle::output(dir,"big_unsigned_r10.tsv",r10(),write);
        independent_oracle::output(dir,"big_unsigned_r11.tsv",r11(),write);
        std::mt19937_64 random(0x52632633474344ULL);
        for(std::size_t count:{2u,3u,4u,5u,8u,16u,127u,128u,129u})
            independent_oracle::output(dir,"big_unsigned_r26_"+s(count)+".tsv",r26(count,random),write);
        independent_oracle::output(dir,"big_unsigned_r27.tsv",r27(),write);
        independent_oracle::output(dir,"big_unsigned_r21.tsv",r21(),write);
        return 0;
    } catch(std::exception const& e) {std::cerr<<"GCD generator: "<<e.what()<<'\n';return 1;}
}
