// SPDX-License-Identifier: MIT
// Standalone oracle for representability in signed/unsigned 64-bit integers.
#include "oracle_exact_integer.hpp"
namespace {
std::string fixture() {
    using independent_oracle::Nat;
    const Nat signed_limit=Nat::power2(63);
    const Nat unsigned_limit=Nat::power2(64);
    std::string result="# Independent integer conversion boundaries.\n# source signed64_ok unsigned64_ok signed32_ok unsigned32_ok\n";
    const std::array<Nat,14> values={Nat(),Nat(1),
        Nat::sub(Nat::power2(31),Nat(1)),Nat::power2(31),Nat::add(Nat::power2(31),Nat(1)),
        Nat::sub(Nat::power2(32),Nat(1)),Nat::power2(32),Nat::add(Nat::power2(32),Nat(1)),Nat::sub(signed_limit,Nat(1)),
        signed_limit,Nat::add(signed_limit,Nat(1)),Nat::sub(unsigned_limit,Nat(1)),
        unsigned_limit,Nat::add(unsigned_limit,Nat(1))};
    for (auto const& value:values) {
        for (bool negative:{false,true}) {
            if (negative && value.zero()) continue;
            const bool signed_ok=negative?value.compare(signed_limit)<=0:value.compare(signed_limit)<0;
            const bool unsigned_ok=!negative && value.compare(unsigned_limit)<0;
            result+=independent_oracle::row({(negative?"-":"")+value.dec(),
                signed_ok?"1":"0",unsigned_ok?"1":"0",
                value.compare(Nat::power2(31))<(negative?1:0)?"1":"0",
                !negative && value.compare(Nat::power2(32))<0?"1":"0"});
        }
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"integer_conversion.tsv",fixture(),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
