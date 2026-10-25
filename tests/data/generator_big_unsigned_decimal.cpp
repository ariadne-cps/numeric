// SPDX-License-Identifier: MIT
// Reproduce A20's 1,800 randomized decimal-parsing oracle visits.
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
}

int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(
            dir,"big_unsigned_a20_decimal_random.tsv",
            a20_decimal_random_fixture(),write);
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"BigUnsigned decimal generator: "<<error.what()<<'\n';
        return 1;
    }
}
