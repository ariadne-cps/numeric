// SPDX-License-Identifier: MIT
// Reproduce B131/B132 decimal-parsing oracle visits.
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

std::string b131_decimal_random_fixture() {
    std::string out=
        "# B131: original 1800 randomized decimal parsing visits from test_gcd_decimal().\n"
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

std::string b132_decimal_structured_fixture() {
    std::string out=
        "# B132: all 484 structured decimal parsing visits from test_gcd_decimal().\n"
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
        throw std::runtime_error("wrong B132 structured decimal case count");
    return out;
}

std::string b133_decimal_format_fixture() {
    std::string out=
        "# B133: all 1509 decimal-formatting oracle visits from test_gcd_decimal().\n"
        "# ordinal source_hex expected_decimal\n"
        "# Preserves original source streams and boundary constructions; no GMP/Numeric/Boost.\n";
    std::size_t ordinal=0u;
    const auto append=[&](Nat const& source) {
        out+=row({std::to_string(ordinal++),source.hex(),source.dec()});
    };
    const auto max=std::numeric_limits<std::uint64_t>::max();

    std::mt19937_64 leaf_engine(0xB250C0DEULL);
    for (const std::size_t count:
            {7u,8u,9u,10u,15u,16u,17u,32u,64u,128u}) {
        for (unsigned mode=0u;mode<5u;++mode) {
            std::vector<std::uint64_t> words(count,0u);
            if (mode==0u) {
                words.back()=1u;
            } else if (mode==1u) {
                std::fill(words.begin(),words.end(),max);
            } else if (mode==2u) {
                words.front()=1u;
                words.back()=1u;
            } else {
                for (auto& word:words) word=leaf_engine();
                if (mode==4u)
                    for (std::size_t i=1u;i+1u<count;++i) words[i]=0u;
                words.back()|=std::uint64_t{1}<<63u;
            }
            append(Nat(std::move(words)));
        }
    }

    for (const unsigned long groups:{8ul,9ul,10ul,16ul,17ul,18ul}) {
        const Nat power=parse_decimal(
            "1"+std::string(18ul*groups,'0'));
        append(power);
        append(Nat::sub(power,Nat(1u)));
        append(Nat::add(power,Nat(1u)));
    }

    std::mt19937_64 divq_decimal_engine(0xB300DEC1ULL);
    for (const std::size_t count:{8u,9u,10u,12u,16u,17u,32u,64u}) {
        for (unsigned mode=0u;mode<16u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=divq_decimal_engine();
            if (mode%8u==0u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u; words.back()=1u;
            } else if (mode%8u==1u) {
                std::fill(words.begin(),words.end(),max);
            } else if (mode%8u==2u) {
                std::fill(words.begin(),words.end(),0u); words.back()=1u;
            } else if (mode%8u==3u) {
                std::fill(words.begin(),words.end(),0u); words.back()=max;
            } else if (mode%8u==4u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u; words.back()=std::uint64_t{1}<<63u;
            } else if (mode%8u==5u) {
                for (std::size_t i=1u;i+1u<count;++i) words[i]=0u;
            }
            words.back()|=1u;
            append(Nat(std::move(words)));
        }
    }

    std::mt19937_64 b42_engine(0xB42DEC1A7ULL);
    for (const std::size_t count:
            {8u,9u,10u,15u,16u,17u,18u,32u,33u,64u,65u,128u,256u}) {
        for (unsigned mode=0u;mode<16u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b42_engine();
            switch (mode%8u) {
            case 0u: std::fill(words.begin(),words.end(),0u); break;
            case 1u: std::fill(words.begin(),words.end(),max); break;
            case 2u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u; break;
            case 3u:
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u; words.back()=std::uint64_t{1}<<63u; break;
            case 4u:
                for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
                break;
            case 5u:
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL}); break;
            case 6u: words.back()=1u; break;
            default: break;
            }
            append(Nat(std::move(words)));
        }
    }

    std::mt19937_64 b46_engine(0xB460DEC1A7ULL);
    for (const std::size_t count:
            {8u,9u,10u,11u,12u,13u,14u,15u,16u,17u,32u}) {
        for (unsigned mode=0u;mode<32u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b46_engine();
            if (mode%8u==0u) std::fill(words.begin(),words.end(),max);
            else if (mode%8u==1u) {
                std::fill(words.begin(),words.end(),0u); words.back()=1u;
            } else if (mode%8u==2u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u; words.back()=1u;
            } else if (mode%8u==3u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            } else if (mode%8u==4u) {
                std::fill_n(words.begin(),count/2u,0u);
            } else if (mode%8u==5u) {
                for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
            } else if (mode%8u==6u) {
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL});
            }
            words.back()|=1u;
            append(Nat(std::move(words)));
        }
    }

    std::mt19937_64 b51_engine(0xB510D3C1ULL);
    for (const std::size_t count:
            {7u,8u,9u,10u,12u,15u,16u,17u,18u,32u,33u,64u}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b51_engine();
            switch (mode%8u) {
            case 0u: std::fill(words.begin(),words.end(),max); break;
            case 1u:
                std::fill(words.begin(),words.end(),0u); words.back()=1u; break;
            case 2u:
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u; words.back()=1u; break;
            case 3u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u; break;
            case 4u: std::fill_n(words.begin(),count/2u,0u); break;
            case 5u:
                for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
                break;
            case 6u:
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL}); break;
            default: break;
            }
            words.back()|=1u;
            append(Nat(std::move(words)));
        }
    }

    std::mt19937_64 b59_engine(0xB590DEC1ULL);
    for (const std::size_t count:
            {15u,16u,17u,31u,32u,33u,63u,64u,65u}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b59_engine();
            switch (mode%8u) {
            case 0u:
                std::fill(words.begin(),words.end(),0u); words.back()=1u; break;
            case 1u: std::fill(words.begin(),words.end(),max); break;
            case 2u:
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u; words.back()=max; break;
            case 3u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u; break;
            case 4u:
                for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
                break;
            case 5u:
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL}); break;
            case 6u:
                for (std::size_t i=1u;i+1u<count;++i) words[i]=0u;
                break;
            default: break;
            }
            words.back()|=std::uint64_t{1}<<63u;
            append(Nat(std::move(words)));
        }
    }

    for (const unsigned long digits:
            {144ul,288ul,576ul,1152ul,2304ul}) {
        const Nat power=parse_decimal("1"+std::string(digits,'0'));
        for (const unsigned long factor:{1ul,7ul,1234567ul,9999999ul}) {
            const Nat base=power.mul_small(factor);
            for (const unsigned long offset:{0ul,1ul,11ul,99ul})
                append(Nat::add(base,Nat(offset)));
        }
    }

    for (const unsigned long digits:
            {9ul,18ul,27ul,72ul,144ul,288ul,576ul}) {
        const Nat power=parse_decimal("1"+std::string(digits,'0'));
        for (const unsigned long offset:
                {0ul,1ul,9ul,10ul,11ul,99ul,100ul,101ul,999ul}) {
            append(Nat::add(power,Nat(offset)));
            if (offset!=0ul) append(Nat::sub(power,Nat(offset)));
        }
    }

    for (unsigned k=0u;k<=9u;++k) {
        const unsigned long digits=9ul*(1ul<<k);
        const Nat power=parse_decimal("1"+std::string(digits,'0'));
        append(power);
        append(Nat::sub(power,Nat(1u)));
        append(Nat::add(power,Nat(1u)));
        append(Nat::add(power,Nat(9u)));
        append(parse_decimal(
            std::string("123456789")+std::string(digits,'0')+
            std::string("000000001")));
    }
    if (ordinal!=1509u)
        throw std::runtime_error("wrong B133 decimal formatting case count");
    return out;
}

}

int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(
            dir,"big_unsigned_b131_decimal_random.tsv",
            b131_decimal_random_fixture(),write);
        independent_oracle::output(
            dir,"big_unsigned_b132_decimal_structured.tsv",
            b132_decimal_structured_fixture(),write);
        independent_oracle::output(
            dir,"big_unsigned_b133_decimal_format.tsv",
            b133_decimal_format_fixture(),write);
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"BigUnsigned decimal generator: "<<error.what()<<'\n';
        return 1;
    }
}
