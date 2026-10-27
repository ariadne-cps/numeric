// SPDX-License-Identifier: MIT
// Independent small exact-rational oracle for rounded minimum/maximum.
#include "oracle_exact_integer.hpp"
#include <numeric>
namespace {
using I=std::int64_t;
struct Value { char kind='F'; I n=0,d=1; bool negative=false; };
Value finite(I n,I d=1,bool negative_zero=false) {
    const I divisor=std::gcd(n,d);
    return {'F',n/divisor,d/divisor,n<0||(n==0&&negative_zero)};
}
Value source(std::string const& text) {
    if (text=="nan") return {'N',0,1,false};
    if (text=="inf" || text=="-inf") return {'I',0,1,text[0]=='-'};
    if (text=="+0" || text=="-0") return finite(0,1,text[0]=='-');
    const auto slash=text.find('/');
    return finite(std::stoll(text.substr(0,slash)),std::stoll(text.substr(slash+1)));
}
Value rounded(Value x,unsigned bits,char mode) {
    if (x.kind!='F' || x.n==0) return x;
    const bool negative=x.n<0;
    const I magnitude=negative?-x.n:x.n;
    int exponent=0;
    if (magnitude>=x.d) { while (magnitude>=x.d*(I{1}<<(exponent+1))) ++exponent; }
    else { while (magnitude*(I{1}<<(-exponent))<x.d) --exponent; }
    const int scale=static_cast<int>(bits)-1-exponent;
    const I n=scale>=0?magnitude*(I{1}<<scale):magnitude;
    const I d=scale<0?x.d*(I{1}<<(-scale)):x.d;
    I quotient=n/d;
    const I remainder=n%d;
    bool increment=false;
    if (remainder!=0) {
        if (mode=='N') increment=2*remainder>d || (2*remainder==d && quotient%2!=0);
        else if (mode=='D') increment=negative;
        else if (mode=='U') increment=!negative;
    }
    quotient+=increment?1:0;
    if (negative) quotient=-quotient;
    return scale>=0?finite(quotient,I{1}<<scale):finite(quotient*(I{1}<<(-scale)));
}
bool less(Value a,Value b) {
    if (a.kind=='I') return a.negative && !(b.kind=='I'&&b.negative);
    if (b.kind=='I') return !b.negative;
    return a.n*b.d<b.n*a.d;
}
Value selected(Value a,Value b,bool maximum) {
    if (a.kind=='N') return b;
    if (b.kind=='N') return a;
    if (a.kind=='F' && b.kind=='F' && a.n==0 && b.n==0)
        return finite(0,1,maximum?a.negative&&b.negative:a.negative||b.negative);
    const bool order=less(a,b);
    return maximum?(order?b:a):(order?a:b);
}
std::string fixture() {
    // Inputs have |n|<=5 and d<=4; all scaled cross-products fit int64_t.
    const std::vector<std::pair<std::string,std::string>> pairs={
        {"nan","nan"},{"nan","5/4"},{"5/4","nan"},{"nan","inf"},{"-inf","nan"},
        {"+0","+0"},{"+0","-0"},{"-0","+0"},{"-0","-0"},
        {"+0","1/4"},{"-1/4","-0"},{"1/4","+0"},{"-0","-1/4"},
        {"1/4","5/4"},{"5/4","1/4"},{"-1/4","-5/4"},{"-5/4","-1/4"},
        {"-5/4","5/4"},{"5/4","-5/4"},{"inf","-inf"},{"-inf","inf"},
        {"inf","5/4"},{"-5/4","-inf"},{"5/4","3/2"},{"3/2","5/4"},
        {"-3/2","-5/4"},{"-5/4","-3/2"}};
    std::string result="# Independent exact rounded FloatMP min/max including NaN, infinity and signed zero.\n"
        "# id operation mode lhs_bits lhs_source rhs_bits rhs_source result_bits kind numerator denominator negative\n";
    unsigned id=0u;
    for (auto const& [a,b]:pairs) for (auto [ap,bp]:std::array<std::pair<unsigned,unsigned>,3>{{{2,7},{7,2},{7,7}}})
    for (char mode:{'N','D','U','Z'}) for (bool maximum:{false,true}) {
        const unsigned p=std::min(ap,bp);
        const auto value=rounded(selected(rounded(source(a),ap,'N'),rounded(source(b),bp,'N'),maximum),p,mode);
        result+=independent_oracle::row({std::to_string(id++),maximum?"max":"min",std::string(1,mode),
            std::to_string(ap),a,std::to_string(bp),b,std::to_string(p),std::string(1,value.kind),
            std::to_string(value.n),std::to_string(value.d),value.negative?"1":"0"});
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatmp_minmax.tsv",fixture(),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
