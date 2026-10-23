// SPDX-License-Identifier: MIT
// Independent test-data generation only. No Numeric, GMP, MPFR or Boost.
#ifndef ARIADNE_TEST_ORACLE_EXACT_INTEGER_HPP
#define ARIADNE_TEST_ORACLE_EXACT_INTEGER_HPP
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace independent_oracle {
struct Nat {
    std::vector<std::uint64_t> w; // little endian
    Nat() = default;
    explicit Nat(std::uint64_t value) { if (value) w.push_back(value); }
    explicit Nat(std::vector<std::uint64_t> data) : w(std::move(data)) { trim(); }
    void trim() { while (!w.empty() && w.back() == 0) w.pop_back(); }
    bool zero() const { return w.empty(); }
    static Nat power2(std::size_t bits) {
        std::vector<std::uint64_t> words(bits/64+1,0);
        words[bits/64] = std::uint64_t{1} << (bits%64);
        return Nat(std::move(words));
    }
    int compare(Nat const& other) const {
        if (w.size() != other.w.size()) return w.size()<other.w.size() ? -1 : 1;
        for (std::size_t i=w.size();i>0;--i) if (w[i-1]!=other.w[i-1]) return w[i-1]<other.w[i-1] ? -1 : 1;
        return 0;
    }
    static Nat add(Nat const& a, Nat const& b) {
        std::vector<std::uint64_t> out(std::max(a.w.size(),b.w.size())+1,0);
        std::uint64_t carry=0;
        for (std::size_t i=0;i+1<out.size();++i) {
            const std::uint64_t x=i<a.w.size()?a.w[i]:0;
            const std::uint64_t y=i<b.w.size()?b.w[i]:0;
            std::uint64_t z=x+y;
            const bool carry1=z<x;
            const std::uint64_t result=z+carry;
            const bool carry2=result<z;
            out[i]=result;
            carry=carry1||carry2;
        }
        out.back()=carry;
        return Nat(std::move(out));
    }
    static Nat sub(Nat const& a, Nat const& b) {
        if (a.compare(b)<0) throw std::runtime_error("Nat subtraction underflow");
        std::vector<std::uint64_t> out=a.w;
        std::uint64_t borrow=0;
        for (std::size_t i=0;i<out.size();++i) {
            const std::uint64_t y=i<b.w.size()?b.w[i]:0;
            const std::uint64_t z=out[i]-y;
            const bool borrow1=out[i]<y;
            const std::uint64_t result=z-borrow;
            const bool borrow2=z<borrow;
            out[i]=result;
            borrow=borrow1||borrow2;
        }
        if (borrow) throw std::runtime_error("Nat subtraction final borrow");
        return Nat(std::move(out));
    }
    Nat shl(std::size_t shift) const {
        if (zero() || shift==0) return *this;
        const std::size_t offset=shift/64;
        const unsigned bits=static_cast<unsigned>(shift%64);
        std::vector<std::uint64_t> out(w.size()+offset+1,0);
        for (std::size_t i=0;i<w.size();++i) {
            out[i+offset] |= w[i]<<bits;
            if (bits) out[i+offset+1] |= w[i]>>(64-bits);
        }
        return Nat(std::move(out));
    }
    Nat shr(std::size_t shift) const {
        if (zero() || shift==0) return *this;
        const std::size_t offset=shift/64;
        if (offset>=w.size()) return Nat();
        const unsigned bits=static_cast<unsigned>(shift%64);
        std::vector<std::uint64_t> out(w.size()-offset,0);
        for(std::size_t i=0;i<out.size();++i) {
            out[i]=w[i+offset]>>bits;
            if(bits && i+offset+1<w.size()) out[i]|=w[i+offset+1]<<(64-bits);
        }
        return Nat(std::move(out));
    }
    std::size_t trailing_zero_bits() const {
        if (zero()) throw std::runtime_error("trailing_zero_bits(0)");
        std::size_t i=0;
        while (w[i]==0) ++i;
        return i*64u+static_cast<std::size_t>(std::countr_zero(w[i]));
    }
    static Nat gcd(Nat a,Nat b) {
        if (a.zero()) return b;
        if (b.zero()) return a;
        const std::size_t common=std::min(a.trailing_zero_bits(),b.trailing_zero_bits());
        a=a.shr(a.trailing_zero_bits());
        b=b.shr(b.trailing_zero_bits());
        while (a.compare(b)!=0) {
            if(a.compare(b)>0) { a=sub(a,b); a=a.shr(a.trailing_zero_bits()); }
            else { b=sub(b,a); b=b.shr(b.trailing_zero_bits()); }
        }
        return a.shl(common);
    }
    Nat mul_small(std::uint64_t multiplier) const {
        Nat x,term=*this;
        for (;multiplier;multiplier>>=1) {
            if (multiplier&1) x=add(x,term);
            if (multiplier>1) term=add(term,term);
        }
        return x;
    }
    std::string hex() const {
        if (zero()) return "0";
        std::ostringstream out;
        out<<std::hex<<std::nouppercase<<w.back();
        for(std::size_t i=w.size()-1;i>0;--i) out<<std::setw(16)<<std::setfill('0')<<w[i-1];
        return out.str();
    }
    std::uint32_t div_base_1e9() {
        constexpr std::uint64_t d=1000000000;
        std::uint64_t rem=0;
        for(std::size_t i=w.size();i>0;--i) {
            const std::uint64_t value=w[i-1];
            const std::uint64_t top=(rem<<32)|(value>>32);
            const std::uint64_t qtop=top/d;
            rem=top%d;
            const std::uint64_t low=(rem<<32)|(value&0xffffffffULL);
            const std::uint64_t qlow=low/d;
            rem=low%d;
            w[i-1]=(qtop<<32)|qlow;
        }
        trim();
        return static_cast<std::uint32_t>(rem);
    }
    std::string dec() const {
        if (zero()) return "0";
        Nat current=*this;
        std::vector<std::uint32_t> blocks;
        while (!current.zero()) blocks.push_back(current.div_base_1e9());
        std::ostringstream out;
        out<<blocks.back();
        for(std::size_t i=blocks.size()-1;i>0;--i)
            out<<std::setw(9)<<std::setfill('0')<<blocks[i-1];
        return out.str();
    }
};
inline std::string row(std::vector<std::string> const& fields) {
    std::ostringstream out;
    for(std::size_t i=0;i<fields.size();++i) { if(i) out<<'\t'; out<<fields[i]; }
    return out.str()+"\n";
}
inline std::filesystem::path data_directory(int argc,char** argv,bool& write) {
    write=false;
    std::filesystem::path path="tests/data";
    for(int i=1;i<argc;++i) {
        std::string a=argv[i];
        if(a=="--write") write=true;
        else if(a=="--verify") write=false;
        else if(a=="--data-dir" && i+1<argc) path=argv[++i];
        else throw std::invalid_argument("Usage: generator [--verify|--write] [--data-dir DIR]");
    }
    return path;
}
inline void output(std::filesystem::path const& dir,std::string const& name,std::string const& expected,bool write) {
    const auto path=dir/name;
    if(write) {
        std::ofstream out(path,std::ios::binary);
        if (!out) throw std::runtime_error("Cannot open output: "+path.string());
        out<<expected;
        if (!out) throw std::runtime_error("Cannot write output: "+path.string());
    } else {
        std::ifstream in(path,std::ios::binary);
        if(!in) throw std::runtime_error("Missing TSV: "+path.string());
        const std::string found{std::istreambuf_iterator<char>{in},std::istreambuf_iterator<char>{}};
        if(found!=expected) {
            std::size_t pos=0;
            while(pos<found.size()&&pos<expected.size()&&found[pos]==expected[pos]) ++pos;
            throw std::runtime_error("TSV mismatch: "+path.string()+" byte "+std::to_string(pos)+" found_length="+std::to_string(found.size())+" regenerated_length="+std::to_string(expected.size()));
        }
    }
    std::cout<<(write?"WROTE ":"VERIFIED ")<<name<<" ("<<expected.size()<<" bytes)\n";
}
}
#endif
