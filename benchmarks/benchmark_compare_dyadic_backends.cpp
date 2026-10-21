/***************************************************************************
 * benchmark_compare_dyadic_backends.cpp
 * Paired native NewDyadic / legacy Dyadic comparison (not GMP dyadic).
 ****************************************************************************/
#include "numeric/new_dyadic.hpp"
#include "numeric/dyadic.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace Ariadne;

template<class T> void consume(T const& value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&value) : "memory");
#else
    static void const* volatile sink=nullptr;
    sink=&value;
#endif
}
unsigned positive(std::string const& s) {
    std::size_t used=0u;
    const unsigned long v=std::stoul(s,&used);
    if(used!=s.size() || !v || v>std::numeric_limits<unsigned>::max())
        throw std::invalid_argument("invalid positive unsigned");
    return static_cast<unsigned>(v);
}
std::vector<unsigned> bits_arg(std::string const& s) {
    std::vector<unsigned> v; std::istringstream in(s);std::string piece;
    while(std::getline(in,piece,',')) {
        auto const n=positive(piece);
        if(std::find(v.begin(),v.end(),n)!=v.end()) throw std::invalid_argument("duplicate bit width");
        v.push_back(n);
    }
    if(v.empty() || s.back()==',')throw std::invalid_argument("invalid bits");
    return v;
}
template<class F> double measure(F const& f,unsigned n,unsigned samples) {
    using Clock=std::chrono::steady_clock;
    std::vector<double> data;data.reserve(samples);
    for(unsigned sample=0;sample<=samples;++sample) {
        auto start=Clock::now();
        for(unsigned i=0;i<n;++i) { auto value=f(i);consume(value); }
        auto stop=Clock::now();
        if(sample) data.push_back(std::chrono::duration<double,std::nano>(stop-start).count()/n);
    }
    std::sort(data.begin(),data.end());
    return data[samples/2u];
}

struct Pool {
    std::vector<Integer> first,second;
    std::vector<Dyadic> old_a,old_b;
    std::vector<NewDyadic> now_a,now_b;
    unsigned exp;
    explicit Pool(unsigned bits):exp(bits/5u+1u) {
        const auto high=shift_left(BigUnsigned(1u),bits-1u);
        for(unsigned i=0;i<32u;++i) {
            first.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(13u+2u*i)),i%4u==0u || i%4u==3u));
            second.push_back(Integer::from_magnitude(
                add(high,BigUnsigned(31u+5u*i)),i%4u==1u || i%4u==3u));
            old_a.emplace_back(first.back(),exp);
            old_b.emplace_back(second.back(),exp+1u);
            now_a.emplace_back(first.back(),exp);
            now_b.emplace_back(second.back(),exp+1u);
        }
        for(unsigned i=0;i<32u;++i) {
            if(now_a[i].to_rational()!=Rational(old_a[i]) ||
               now_b[i].to_rational()!=Rational(old_b[i]) ||
               (now_a[i]+now_b[i]).to_rational()!=Rational(old_a[i]+old_b[i]) ||
               (now_a[i]-now_b[i]).to_rational()!=Rational(old_a[i]-old_b[i]) ||
               (now_a[i]*now_b[i]).to_rational()!=Rational(old_a[i]*old_b[i]) ||
               pow(now_a[i],Nat(3u)).to_rational()!=Rational(pow(old_a[i],Nat(3u))) ||
               floor(now_a[i])!=floor(old_a[i]) ||
               ceil(now_a[i])!=ceil(old_a[i]) ||
               round(now_a[i])!=round(old_a[i]) ||
               (now_a[i]<now_b[i])!=(old_a[i]<old_b[i])) {
                throw std::runtime_error("NewDyadic / Dyadic exact verification mismatch");
            }
        }
    }
};
void run(unsigned bits,unsigned iterations,unsigned samples) {
    const Pool x(bits);
    constexpr unsigned mask=31u;
    auto pair=[&](char const* name,auto&& old_fn,auto&& new_fn) {
        const double old_ns=measure(old_fn,iterations,samples);
        const double new_ns=measure(new_fn,iterations,samples);
        std::cout<<std::left<<std::setw(23)<<name<<std::setw(7)<<bits
                 <<std::right<<std::fixed<<std::setprecision(3)
                 <<std::setw(16)<<old_ns<<std::setw(16)<<new_ns
                 <<std::setw(12)<<new_ns/old_ns<<'\n';
    };
    pair("construct_integer",[&](unsigned i){return Dyadic(x.first[i&mask]);},
         [&](unsigned i){return NewDyadic(x.first[i&mask]);});
    pair("construct_fraction",[&](unsigned i){return Dyadic(x.first[i&mask],x.exp);},
         [&](unsigned i){return NewDyadic(x.first[i&mask],x.exp);});
    pair("copy",[&](unsigned i){return Dyadic(x.old_a[i&mask]);},
         [&](unsigned i){return NewDyadic(x.now_a[i&mask]);});
    pair("add",[&](unsigned i){return x.old_a[i&mask]+x.old_b[i&mask];},
         [&](unsigned i){return x.now_a[i&mask]+x.now_b[i&mask];});
    pair("sub",[&](unsigned i){return x.old_a[i&mask]-x.old_b[i&mask];},
         [&](unsigned i){return x.now_a[i&mask]-x.now_b[i&mask];});
    pair("mul",[&](unsigned i){return x.old_a[i&mask]*x.old_b[i&mask];},
         [&](unsigned i){return x.now_a[i&mask]*x.now_b[i&mask];});
    pair("square",[&](unsigned i){return sqr(x.old_a[i&mask]);},
         [&](unsigned i){return sqr(x.now_a[i&mask]);});
    pair("negate",[&](unsigned i){return -x.old_a[i&mask];},
         [&](unsigned i){return -x.now_a[i&mask];});
    pair("abs",[&](unsigned i){return abs(x.old_a[i&mask]);},
         [&](unsigned i){return abs(x.now_a[i&mask]);});
    pair("half",[&](unsigned i){return hlf(x.old_a[i&mask]);},
         [&](unsigned i){return hlf(x.now_a[i&mask]);});
    pair("compare",[&](unsigned i){return cmp(x.old_a[i&mask],x.old_b[i&mask]);},
         [&](unsigned i){return cmp(x.now_a[i&mask],x.now_b[i&mask]);});
    pair("min",[&](unsigned i){return min(x.old_a[i&mask],x.old_b[i&mask]);},
         [&](unsigned i){return min(x.now_a[i&mask],x.now_b[i&mask]);});
    pair("max",[&](unsigned i){return max(x.old_a[i&mask],x.old_b[i&mask]);},
         [&](unsigned i){return max(x.now_a[i&mask],x.now_b[i&mask]);});
    pair("floor",[&](unsigned i){return floor(x.old_a[i&mask]);},
         [&](unsigned i){return floor(x.now_a[i&mask]);});
    pair("ceil",[&](unsigned i){return ceil(x.old_a[i&mask]);},
         [&](unsigned i){return ceil(x.now_a[i&mask]);});
    pair("round",[&](unsigned i){return round(x.old_a[i&mask]);},
         [&](unsigned i){return round(x.now_a[i&mask]);});
    pair("to_double",[&](unsigned i){return x.old_a[i&mask].get_d();},
         [&](unsigned i){return x.now_a[i&mask].get_d();});
    pair("pow3",[&](unsigned i){return pow(x.old_a[i&mask],Nat(3u));},
         [&](unsigned i){return pow(x.now_a[i&mask],Nat(3u));});
}
}
int main(int argc,char const* argv[]) {
    unsigned iterations=20000u,samples=31u;
    std::vector<unsigned> widths{53u,128u,256u,1024u,2048u,4096u};
    try {
        for(int i=1;i<argc;++i) {
            std::string key=argv[i];
            if(++i>=argc)throw std::invalid_argument("missing argument "+key);
            std::string value=argv[i];
            if(key=="--bits") widths=bits_arg(value);
            else if(key=="--iterations") iterations=positive(value);
            else if(key=="--samples") samples=positive(value);
            else throw std::invalid_argument("unknown option "+key);
        }
        if(samples<3u || samples%2u==0u)
            throw std::invalid_argument("samples must be odd and >=3");
        std::cout<<"Dyadic vs NewDyadic | iterations="<<iterations
                 <<" samples="<<samples<<" | reference=legacy Dyadic (mpf)\n";
        std::cout<<std::left<<std::setw(23)<<"Operation"<<std::setw(7)<<"Bits"
                 <<std::right<<std::setw(16)<<"Dyadic(ns)"
                 <<std::setw(16)<<"NewDyadic(ns)"<<std::setw(12)<<"New/Old\n";
        for(unsigned bits:widths)run(bits,iterations,samples);
        std::cout<<"Exact verification PASS: native Rational vs legacy Dyadic\n";
    } catch(std::exception const& e) {
        std::cerr<<"Dyadic comparison failed: "<<e.what()<<'\n';
        return 1;
    }
    return 0;
}
