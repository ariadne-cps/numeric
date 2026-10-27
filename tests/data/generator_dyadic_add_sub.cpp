// SPDX-License-Identifier: MIT
// Independent exact signed dyadic addition/subtraction, no Numeric/GMP/Boost.
#include "oracle_exact_integer.hpp"
namespace {
using independent_oracle::Nat;
struct Value { Nat magnitude; bool negative; std::size_t exponent; };
Value normalized(Value x) {
    if (x.magnitude.zero()) return {Nat(),false,0u};
    const auto remove=std::min(x.exponent,x.magnitude.trailing_zero_bits());
    x.magnitude=x.magnitude.shr(remove); x.exponent-=remove; return x;
}
std::string coefficient(Value const& x) { return (x.negative?"-":"")+x.magnitude.dec(); }
Value sum(Value const& a,Value const& b,bool subtract) {
    const auto exponent=std::max(a.exponent,b.exponent);
    const auto left=a.magnitude.shl(exponent-a.exponent);
    const auto right=b.magnitude.shl(exponent-b.exponent);
    const bool negative_b=b.negative!=subtract;
    if (a.negative==negative_b) return normalized({Nat::add(left,right),a.negative,exponent});
    const int order=left.compare(right);
    return normalized({order>=0?Nat::sub(left,right):Nat::sub(right,left),
                       order>=0?a.negative:negative_b,exponent});
}
std::string fixture(unsigned limbs) {
    std::string result="# Independent exact dyadic signed add/sub at dispatch boundaries.\n"
        "# id lhs_coefficient lhs_exponent rhs_coefficient rhs_exponent sum_coefficient sum_exponent difference_coefficient difference_exponent\n";
    const Nat a=Nat::sub(Nat::power2(64u*limbs),Nat(1));
    unsigned id=0u;
    for (unsigned base:{0u,64u*limbs+3u})
    for (unsigned shift:{0u,1u,63u,64u})
    for (bool reversed:{false,true})
    for (unsigned shape=0u;shape<4u;++shape)
    for (bool negative_a:{false,true}) for (bool negative_b:{false,true}) {
        const Nat b=shape==0u?a:shape==1u?Nat::sub(a,Nat(2)):
                    shape==2u?Nat::add(a,Nat(2)):Nat();
        const auto left=normalized({a,negative_a,base+(reversed?shift:0u)});
        const auto right=normalized({b,negative_b,base+(reversed?0u:shift)});
        const auto added=sum(left,right,false),subtracted=sum(left,right,true);
        result+=independent_oracle::row({std::to_string(id++),coefficient(left),
            std::to_string(left.exponent),coefficient(right),std::to_string(right.exponent),
            coefficient(added),std::to_string(added.exponent),coefficient(subtracted),
            std::to_string(subtracted.exponent)});
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        for (unsigned limbs:{1u,5u,16u,17u,31u,32u,64u,65u})
            independent_oracle::output(dir,"dyadic_add_sub_"+std::to_string(limbs)+".tsv",fixture(limbs),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
