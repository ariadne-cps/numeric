/***************************************************************************
 *            benchmark_integer_small_power_profile.cpp
 *
 *  Copyright 2026 Ariadne contributors
 ***************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file benchmark_integer_small_power_profile.cpp
 *  \brief B64 independent GMP-oracle power/square ownership diagnostics.
 *
 *  Production pow(), square(), and the official 84-case Integer/GMP
 *  benchmark are not modified. All signed inputs are exactly the
 *  historical GMP-generated 32-input Integer corpus. Alternative
 *  computations return owning Integer results with the same semantics.
 *  Independent timings are not additive.
 */

#include "benchmark_exact_arithmetic_report.hpp"

#include <bit>
#include <cstddef>
#include <random>
#include <vector>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace Ariadne {
// B65 internal comparison entrypoint, not part of the public API.
BigUnsigned square_short_symmetric_diagnostic(BigUnsigned const& input);
}

namespace {

using namespace Ariadne;
using namespace AriadneBenchmark;

// B64 control: generic exponentiation by squaring with deferred ownership
// of the first power. This is a whole-algorithm, not exponent-3-only,
// alternative; it never edits its input and works for exponent zero.
Integer pow_lazy_diagnostic(Integer const& input,unsigned exponent)
{
    if (exponent==0u) {
        return Integer::from_magnitude(BigUnsigned(std::uint64_t{1}),false);
    }
    const bool negative=input.negative() && (exponent&1u)!=0u;
    BigUnsigned accumulator;
    BigUnsigned squared;
    BigUnsigned const* current=&input.magnitude();
    bool has_factor=false;
    while (exponent!=0u) {
        if ((exponent&1u)!=0u) {
            if (!has_factor) {
                accumulator=*current;
                has_factor=true;
            } else {
                accumulator=mul(accumulator,*current);
            }
        }
        exponent>>=1u;
        if (exponent!=0u) {
            squared=square(*current);
            current=&squared;
        }
    }
    return Integer::from_magnitude(std::move(accumulator),negative);
}

// B65: exact generic left-to-right binary exponentiation for all
// unsigned powers. The original magnitude is borrowed until the first
// square; no initial BigUnsigned copy is materialized for exponent >=2.
// Exponent 3 *naturally* follows square(x) then mul(x^2,x), unlike
// a benchmark-specific cube fast path.
Integer pow_msb_diagnostic(Integer const& input,unsigned exponent)
{
    if (exponent==0u) {
        return Integer::from_magnitude(BigUnsigned(std::uint64_t{1}),false);
    }
    const bool negative=input.negative() && (exponent&1u)!=0u;
    const BigUnsigned& base=input.magnitude();
    const BigUnsigned* current=&base;
    BigUnsigned result;
    unsigned bit=std::bit_floor(exponent)>>1u;
    while (bit!=0u) {
        result=square(*current);
        current=&result;
        if ((exponent&bit)!=0u) {
            result=mul(result,base);
        }
        bit>>=1u;
    }
    if (current==&base) {
        result=base; // Exponent 1 still requires a distinct owned result.
    }
    return Integer::from_magnitude(std::move(result),negative);
}

// Lower-bound diagnostic for cube ownership overhead, NOT a production
// special-case recommendation. Both operands remain the same exact input.
Integer pow3_direct_chain_diagnostic(Integer const& input)
{
    const auto& magnitude=input.magnitude();
    return Integer::from_magnitude(
        mul(square(magnitude),magnitude),input.negative());
}

Integer square_magnitude_diagnostic(Integer const& input)
{
    return Integer::from_magnitude(square(input.magnitude()),false);
}

Integer square_short_symmetric_owning_diagnostic(Integer const& input)
{
    return Integer::from_magnitude(
        square_short_symmetric_diagnostic(input.magnitude()),false);
}

void require_equal_to_gmp(Integer const& actual,mpz_srcptr expected,
                          std::string const& context)
{
    RawInteger got;
    const auto view=actual.magnitude().limbs();
    mpz_import(got.value,view.size(),-1,sizeof(BigUnsigned::Limb),
               0,0,view.data());
    if (actual.negative()) {
        mpz_neg(got.value,got.value);
    }
    if (mpz_cmp(got.value,expected)!=0) {
        throw std::runtime_error("B64 GMP mismatch: "+context);
    }
}

void verify_input(Integer const& native,RawInteger const& gmp,
                  unsigned bits,unsigned index)
{
    const std::string context=std::to_string(bits)+" bits case "+
                              std::to_string(index);
    require_equal_to_gmp(native,gmp.value,context+" input");
    RawInteger expected;
    for (const unsigned exponent:{0u,1u,2u,3u,4u,5u,7u,8u,9u,12u}) {
        mpz_pow_ui(expected.value,gmp.value,exponent);
        require_equal_to_gmp(pow(native,Nat(exponent)),expected.value,
                             context+" production exponent "+std::to_string(exponent));
        require_equal_to_gmp(pow_lazy_diagnostic(native,exponent),expected.value,
                             context+" lazy exponent "+std::to_string(exponent));
        require_equal_to_gmp(pow_msb_diagnostic(native,exponent),expected.value,
                             context+" msb exponent "+std::to_string(exponent));
        if (exponent==3u) {
            require_equal_to_gmp(pow3_direct_chain_diagnostic(native),
                                 expected.value,context+" direct cube");
        }
    }
    mpz_mul(expected.value,gmp.value,gmp.value);
    require_equal_to_gmp(native*native,expected.value,context+" square product");
    require_equal_to_gmp(square_magnitude_diagnostic(native),
                         expected.value,context+" magnitude square");
    require_equal_to_gmp(square_short_symmetric_owning_diagnostic(native),
                         expected.value,context+" short symmetric square");
}

// B65 additional GMP oracle for the 0..8-limb short-square dispatch
// boundary, including zero, saturated, alternating sparse, and random
// magnitudes. These are extra diagnostic tests, not benchmark inputs.
void verify_short_squares()
{
    std::mt19937_64 generator(0xB650A11ULL);
    RawInteger expected,input_gmp;
    for (const std::size_t count:{0u,1u,2u,3u,4u,5u,6u,7u,8u}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            std::vector<std::uint64_t> digits(count);
            for (auto& word:digits) word=generator();
            if (mode%6u==0u) {
                std::fill(digits.begin(),digits.end(),std::uint64_t{0});
            } else if (mode%6u==1u) {
                std::fill(digits.begin(),digits.end(),
                          std::numeric_limits<std::uint64_t>::max());
            } else if (mode%6u==2u && !digits.empty()) {
                std::fill(digits.begin(),digits.end(),std::uint64_t{0});
                digits.back()=std::uint64_t{1}<<63u;
            } else if (mode%6u==3u) {
                for (std::size_t i=0u;i<count;i+=2u) digits[i]=0u;
            }
            const BigUnsigned magnitude(digits);
            const Integer signed_input=Integer::from_magnitude(
                magnitude,(mode&1u)!=0u);
            const auto limbs=magnitude.limbs();
            mpz_import(input_gmp.value,limbs.size(),-1,
                       sizeof(BigUnsigned::Limb),0,0,limbs.data());
            mpz_mul(expected.value,input_gmp.value,input_gmp.value);
            require_equal_to_gmp(
                square_short_symmetric_owning_diagnostic(signed_input),
                expected.value,"B65 forced short square "+std::to_string(count)+
                               " limbs mode "+std::to_string(mode));
        }
    }
}

void verify(unsigned bits,
            IntegerPool<IntegerOperations> const& native,
            IntegerPool<GmpIntegerOperations> const& gmp)
{
    for (unsigned i=0u;i<IntegerPool<IntegerOperations>::pool_size;++i) {
        verify_input(native.a[i],gmp.a[i],bits,i);
    }
    // Extra identities cover the zero and unit cases which the official
    // positive-magnitude near-power-of-two generator cannot produce.
    for (const std::string literal:{"0","1","-1","2","-2"}) {
        const Integer a{String(literal)};
        const RawInteger b(literal);
        verify_input(a,b,bits,32u);
    }
}

void report(char const* name,double native,double gmp,double baseline)
{
    std::ostringstream output;
    output << std::left << std::setw(27) << name
           << std::right << std::fixed << std::setprecision(3)
           << std::setw(11) << native << " ns  "
           << std::setw(11) << native/gmp << "x GMP  "
           << std::setw(11) << native/baseline << "x production";
    LOGGING_PRINTLN(output.str())
}

void run_size(unsigned bits,unsigned iterations,unsigned samples)
{
    IntegerPool<IntegerOperations> native(bits);
    IntegerPool<GmpIntegerOperations> gmp(bits);
    verify(bits,native,gmp);
    const auto index=[](unsigned i) { return i%32u; };

    const double pow_gmp=exact_arithmetic_median_ns([&](unsigned i) {
        return GmpIntegerOperations::power3(gmp.a[index(i)]);
    },iterations,samples);
    const double pow_native=exact_arithmetic_median_ns([&](unsigned i) {
        return pow(native.a[index(i)],Nat(3u));
    },iterations,samples);
    const double pow_lazy=exact_arithmetic_median_ns([&](unsigned i) {
        return pow_lazy_diagnostic(native.a[index(i)],3u);
    },iterations,samples);
    const double pow_msb=exact_arithmetic_median_ns([&](unsigned i) {
        return pow_msb_diagnostic(native.a[index(i)],3u);
    },iterations,samples);
    const double pow_direct=exact_arithmetic_median_ns([&](unsigned i) {
        return pow3_direct_chain_diagnostic(native.a[index(i)]);
    },iterations,samples);
    const double square_gmp=exact_arithmetic_median_ns([&](unsigned i) {
        return GmpIntegerOperations::square(gmp.a[index(i)]);
    },iterations,samples);
    const double square_native=exact_arithmetic_median_ns([&](unsigned i) {
        return native.a[index(i)]*native.a[index(i)];
    },iterations,samples);
    const double square_direct=exact_arithmetic_median_ns([&](unsigned i) {
        return square_magnitude_diagnostic(native.a[index(i)]);
    },iterations,samples);
    const double square_short=exact_arithmetic_median_ns([&](unsigned i) {
        return square_short_symmetric_owning_diagnostic(native.a[index(i)]);
    },iterations,samples);

    // Guard against any result corruption or retained mutable input state.
    verify(bits,native,gmp);
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN(bits << " bits | 32 official GMP-verified signed inputs"
                    << " | iterations=" << iterations
                    << " | samples=" << samples)
    LOGGING_PRINTLN("  GMP pow3: " << std::fixed
                    << std::setprecision(3) << pow_gmp << " ns")
    report("production pow(z,3)",pow_native,pow_gmp,pow_native);
    report("generic lazy power",pow_lazy,pow_gmp,pow_native);
    report("generic MSB-first power",pow_msb,pow_gmp,pow_native);
    report("direct cube control",pow_direct,pow_gmp,pow_native);
    LOGGING_PRINTLN("  GMP square: " << std::fixed
                    << std::setprecision(3) << square_gmp << " ns")
    report("production z*z",square_native,square_gmp,square_native);
    report("direct magnitude square",square_direct,square_gmp,square_native);
    if (bits<=256u) {
        report("forced short symmetric",square_short,square_gmp,square_native);
    }
    LOGGING_PRINTLN("Diagnostic rows are independently timed owning results; "
                    "no stage percentages are additive.")
}

} // namespace

int main(int argc,const char* argv[])
{
    try {
        BenchmarkCommandLine::initialize_logging();
        ExactArithmeticSelection selection;
        selection.sizes={53u,128u,256u,1024u,2048u,4096u};
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution(
            "numeric-integer-small-power-profile",
            exact_arithmetic_parsers(selection,false)));
        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }
        verify_short_squares();
        for (unsigned bits:selection.sizes) {
            run_size(bits,selection.iterations,selection.samples);
        }
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("Numerical verification PASS (independent GMP "
                        "signed powers 0/1/2/3/4/5/7/8/9/12, square, "
                        "extra zero and signed-unit identities, "
                        "forced 0..8-limb symmetric squares).")
        return 0;
    } catch (std::exception const& error) {
        std::cerr << "B64 small power diagnostic FAILED: "
                  << error.what() << '\n';
        return 1;
    }
}
