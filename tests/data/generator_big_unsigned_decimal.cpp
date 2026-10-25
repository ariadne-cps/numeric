// SPDX-License-Identifier: MIT
// Reproduce A20/A21 decimal-parsing oracle visits.
// Standalone C++20; no Numeric, GMP, Boost or external runtime.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

Nat parse_decimal(std::string const& digits) {
    Nat value;
    for (const char digit:digits) {
        value=value.mul_small(10u);
        value=Nat::add(value,Nat(static_cast<std::uint64_t>(digit-'0')));
    }
    return value;
}

void replay_b76(std::mt19937_64& engine) {
    for (unsigned repetition=0u;repetition<2800u;++repetition) {
        const std::size_t n=repetition%113u==0u ? 96u :
            static_cast<std::size_t>(engine()%20u);
        const std::size_t m=repetition%71u==0u ? 32u :
            static_cast<std::size_t>(engine()%12u);
        for (std::size_t i=0u;i<n;++i) (void)engine();
        for (std::size_t i=0u;i<m;++i) (void)engine();
    }
}

std::string a20_decimal_random_fixture() {
    std::string out=
        "# A20: original 1800 randomized decimal parsing visits from test_gcd_decimal().\n"
        "# repetition length source_decimal expected_hex canonical_decimal\n"
        "# Exact source stream after replaying all 2800 B76 mt19937_64(0xB16B00B5) visits.\n"
        "# Expected values computed with independent multiple-precision integer parsing; no GMP/Numeric/Boost.\n";
    std::mt19937_64 engine(0xB16B00B5ULL);
    replay_b76(engine);
    for (unsigned repetition=0u;repetition<1800u;++repetition) {
        const std::size_t length=repetition%149u==0u ? 1900u :
            1u+static_cast<std::size_t>(engine()%370u);
        std::string digits(length,'0');
        for (auto& digit:digits) {
            digit=static_cast<char>('0'+static_cast<int>(engine()%10u));
        }
        const Nat expected=parse_decimal(digits);
        out+=row({std::to_string(repetition),std::to_string(length),
                  digits,expected.hex(),expected.dec()});
    }
    return out;
}

void replay_a20(std::mt19937_64& engine) {
    for (unsigned repetition=0u;repetition<1800u;++repetition) {
        const std::size_t length=repetition%149u==0u ? 1900u :
            1u+static_cast<std::size_t>(engine()%370u);
        for (std::size_t i=0u;i<length;++i) (void)engine();
    }
}

std::string a21_decimal_structured_fixture() {
    std::string out=
        "# A21: all 484 structured decimal parsing visits from test_gcd_decimal().\n"
        "# ordinal source_decimal expected_hex canonical_decimal\n"
        "# Preserves every original check_long_parse input and the shared mt19937_64 stream.\n"
        "# Exact expectations computed with independent multiple-precision integers; no GMP/Numeric/Boost.\n";
    std::mt19937_64 engine(0xB16B00B5ULL);
    replay_b76(engine);
    replay_a20(engine);
    std::size_t ordinal=0u;
    const auto append=[&](std::string const& digits) {
        const Nat expected=parse_decimal(digits);
        out+=row({std::to_string(ordinal++),digits,
                  expected.hex(),expected.dec()});
    };
    for (const std::size_t length:{1u,19u,31u,32u,76u,77u,78u,95u}) {
        append(std::string(length,'0'));
        append(std::string(length-1u,'0')+"7");
        append(std::string(length,'9'));
    }
    for (const std::size_t length:
            {308u,309u,310u,311u,312u,616u,617u,1234u,1900u}) {
        append(std::string(length,'9'));
        append("1"+std::string(length-2u,'0')+"1");
        append(std::string(length-1u,'0')+"1");
        std::string mixed(length,'0');
        mixed.front()='8';
        mixed[18u]='7';
        mixed.back()='6';
        append(mixed);
        std::string alternating(length,'9');
        for (std::size_t i=1u;i<length;i+=2u) alternating[i]='0';
        append(alternating);
    }
    for (const std::size_t length:
            {16u,17u,18u,19u,31u,32u,33u,35u,36u,37u,38u,
             39u,40u,53u,54u,55u,76u,77u,78u,79u,
             307u,308u,309u,310u,311u,617u}) {
        append(std::string(length,'9'));
        append(std::string(length,'0'));
        append("1"+std::string(length-2u,'0')+"1");
        append(std::string(length-1u,'0')+"7");
        std::string alternating(length,'0');
        for (std::size_t i=0u;i<length;i+=2u) alternating[i]='9';
        append(alternating);
        std::string random_digits(length,'0');
        for (auto& digit:random_digits)
            digit=static_cast<char>('0'+static_cast<int>(engine()%10u));
        append(random_digits);
    }
    for (const std::size_t length:
            {18u,19u,20u,31u,32u,33u,36u,37u,38u,39u,40u,
             56u,57u,58u,75u,76u,77u,78u,94u,95u,96u,
             113u,114u,115u,151u,152u,153u,307u,308u,309u,
             310u,311u,616u,617u,618u,619u,1233u,1234u,1235u}) {
        append(std::string(length,'9'));
        append(std::string(length,'0'));
        append("1"+std::string(length-2u,'0')+"1");
        append(std::string(length-1u,'0')+"7");
        std::string alternating(length,'0');
        for (std::size_t i=0u;i<length;i+=2u) alternating[i]='9';
        append(alternating);
        std::string random_digits(length,'0');
        for (auto& digit:random_digits)
            digit=static_cast<char>('0'+static_cast<int>(engine()%10u));
        append(random_digits);
    }
    for (const unsigned bits:{128u,256u,1024u,2048u,4096u}) {
        const Nat base=Nat::power2(bits-1u);
        for (const unsigned long offset:{0ul,13ul,31ul,75ul,187ul})
            append(Nat::add(base,Nat(offset)).dec());
    }
    if (ordinal!=484u)
        throw std::runtime_error("wrong A21 structured decimal case count");
    return out;
}
}

int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(
            dir,"big_unsigned_a20_decimal_random.tsv",
            a20_decimal_random_fixture(),write);
        independent_oracle::output(
            dir,"big_unsigned_a21_decimal_structured.tsv",
            a21_decimal_structured_fixture(),write);
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"BigUnsigned decimal generator: "<<error.what()<<'\n';
        return 1;
    }
}
