/***************************************************************************
 *            benchmark_integer_small_operations_profile.cpp
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

/*! \file benchmark_integer_small_operations_profile.cpp
 *  \brief B67 owning Integer versus owning magnitude primitive diagnostics.
 *
 *  Diagnostic only; the original 84-operation official benchmark remains
 *  completely unchanged. Signed Integer/GMP rows have identical semantics.
 *  "abs core" rows use the corresponding absolute magnitudes, omit the
 *  signed wrapper, and are NOT an additive decomposition of either row.
 *  Unlike a benchmark-specific shortcut, all paths run production APIs.
 */

#include "benchmark_exact_arithmetic_report.hpp"

#include <cstddef>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {

using namespace Ariadne;
using namespace AriadneBenchmark;

void check_magnitude(BigUnsigned const& native,mpz_srcptr expected,
                     std::string const& context)
{
    RawInteger actual;
    const auto limbs=native.limbs();
    mpz_import(actual.value,limbs.size(),-1,
               sizeof(BigUnsigned::Limb),0,0,limbs.data());
    if (mpz_cmp(actual.value,expected)!=0) {
        throw std::runtime_error("B67 GMP unsigned oracle mismatch: "+context);
    }
}

std::string_view abs_digits(std::string const& signed_decimal)
{
    std::string_view view(signed_decimal);
    if (!view.empty() && (view.front()=='-' || view.front()=='+')) {
        view.remove_prefix(1u);
    }
    return view;
}

// B68 candidates preserve the COMPLETE signed value/result contract.
// Only String adapter and sign-handling steps differ from production.
Integer parse_direct_signed_diagnostic(std::string const& source)
{
    std::string_view digits(source);
    if (digits.empty()) {
        throw std::invalid_argument("Integer empty decimal literal");
    }
    const bool negative=digits.front()=='-';
    if (negative || digits.front()=='+') {
        digits.remove_prefix(1u);
    }
    return Integer::from_magnitude(from_decimal(digits),negative);
}

std::string format_direct_signed_diagnostic(Integer const& value)
{
    std::string digits=to_decimal(value.magnitude());
    if (value.negative()) {
        digits.insert(digits.begin(),'-');
    }
    return digits;
}

// This is semantically identical to Integer::mul for nonaliased inputs;
// unlike unsigned B67 core rows, the result retains signed semantics.
Integer multiply_direct_signed_diagnostic(Integer const& first,
                                          Integer const& second)
{
    return Integer::from_magnitude(
        mul(first.magnitude(),second.magnitude()),
        first.negative()!=second.negative());
}

Integer square_direct_signed_diagnostic(Integer const& input)
{
    return Integer::from_magnitude(square(input.magnitude()),false);
}

void verify_signed_diagnostics(unsigned bits,
                               IntegerPool<IntegerOperations> const& native,
                               IntegerPool<GmpIntegerOperations> const& gmp)
{
    RawInteger expected,got;
    auto check=[&](Integer const& value,mpz_srcptr oracle,
                   std::string const& context) {
        const auto words=value.magnitude().limbs();
        mpz_import(got.value,words.size(),-1,
                   sizeof(BigUnsigned::Limb),0,0,words.data());
        if (value.negative()) {
            mpz_neg(got.value,got.value);
        }
        if (mpz_cmp(got.value,oracle)!=0) {
            throw std::runtime_error("B68 signed GMP mismatch: "+context);
        }
    };

    for (unsigned i=0u;i<32u;++i) {
        const std::string context=std::to_string(bits)+" bits #"+
                                  std::to_string(i);
        check(parse_direct_signed_diagnostic(native.spec[i].a),
              gmp.a[i].value,context+" parse a");
        check(parse_direct_signed_diagnostic(native.spec[i].b),
              gmp.b[i].value,context+" parse b");
        if (format_direct_signed_diagnostic(native.a[i])!=
                decimal(gmp.a[i].value) ||
                format_direct_signed_diagnostic(native.b[i])!=
                decimal(gmp.b[i].value)) {
            throw std::runtime_error("B68 signed decimal output differs: "+
                                     context);
        }

        mpz_mul(expected.value,gmp.a[i].value,gmp.b[i].value);
        check(multiply_direct_signed_diagnostic(native.a[i],native.b[i]),
              expected.value,context+" multiply");
        mpz_mul(expected.value,gmp.a[i].value,gmp.a[i].value);
        check(square_direct_signed_diagnostic(native.a[i]),
              expected.value,context+" square");
    }

    // The original pool contains negative inputs but not negative zero.
    // Preserve normalization and handling of leading signs explicitly.
    for (const std::string text:{"0","-0","+0","1","-1",
                                  "18446744073709551616",
                                  "-18446744073709551616"}) {
        const Integer ref{String(text)};
        const Integer candidate=parse_direct_signed_diagnostic(text);
        if (candidate.negative()!=ref.negative() ||
                format_direct_signed_diagnostic(candidate)!=
                std::string(ref.literal())) {
            throw std::runtime_error("B68 sign/zero normalization mismatch");
        }
    }
}

template<class Reference,class Production,class Candidate>
void report_signed_candidate(char const* name,
                             Reference&& reference,Production&& production,
                             Candidate&& candidate,
                             unsigned iterations,unsigned samples)
{
    const double g=exact_arithmetic_median_ns(
        std::forward<Reference>(reference),iterations,samples);
    const double n=exact_arithmetic_median_ns(
        std::forward<Production>(production),iterations,samples);
    const double c=exact_arithmetic_median_ns(
        std::forward<Candidate>(candidate),iterations,samples);
    std::ostringstream line;
    line << std::left << std::setw(22) << name
         << std::fixed << std::setprecision(3)
         << std::right << std::setw(12) << g
         << std::setw(14) << n
         << std::setw(14) << c
         << std::setw(14) << c/n << "x"
         << std::setw(14) << c/g << "x";
    LOGGING_PRINTLN(line.str())
}

void verify(unsigned bits,IntegerPool<IntegerOperations> const& native,
            IntegerPool<GmpIntegerOperations> const& gmp)
{
    RawInteger first,second,expected;
    for (unsigned i=0u;i<32u;++i) {
        const auto& a=native.a[i].magnitude();
        const auto& b=native.b[i].magnitude();
        mpz_abs(first.value,gmp.a[i].value);
        mpz_abs(second.value,gmp.b[i].value);
        const std::string label=std::to_string(bits)+" bits #"+
                                std::to_string(i);
        check_magnitude(a,first.value,label+" input a");
        check_magnitude(b,second.value,label+" input b");
        check_magnitude(from_decimal(abs_digits(native.spec[i].a)),
                        first.value,label+" decimal construction");
        if (to_decimal(a)!=decimal(first.value)) {
            throw std::runtime_error("B67 decimal magnitude mismatch: "+label);
        }

        mpz_mul(expected.value,first.value,second.value);
        check_magnitude(mul(a,b),expected.value,label+" mul");
        mpz_mul(expected.value,first.value,first.value);
        check_magnitude(square(a),expected.value,label+" square");
        mpz_pow_ui(expected.value,first.value,3u);
        check_magnitude(mul(square(a),a),expected.value,label+" cube");

        const auto& numerator=native.division_a[i].magnitude();
        const auto& divisor=native.division_b[i].magnitude();
        mpz_abs(first.value,gmp.division_a[i].value);
        mpz_abs(second.value,gmp.division_b[i].value);
        check_magnitude(numerator,first.value,label+" division a");
        check_magnitude(divisor,second.value,label+" division b");
        if (mpz_sgn(second.value)==0) {
            throw std::runtime_error("B67 unexpected zero divisor");
        }
        mpz_tdiv_q(expected.value,first.value,second.value);
        check_magnitude(quot(numerator,divisor),
                        expected.value,label+" quotient");
        mpz_tdiv_r(expected.value,first.value,second.value);
        check_magnitude(rem(numerator,divisor),
                        expected.value,label+" remainder");
    }
}

template<class GMP,class Native,class Core>
void report_row(char const* name,GMP&& gmp,Native&& native,
                Core&& core,unsigned iterations,unsigned samples)
{
    const double g=exact_arithmetic_median_ns(
        std::forward<GMP>(gmp),iterations,samples);
    const double n=exact_arithmetic_median_ns(
        std::forward<Native>(native),iterations,samples);
    const double c=exact_arithmetic_median_ns(
        std::forward<Core>(core),iterations,samples);
    std::ostringstream line;
    line << std::left << std::setw(21) << name
         << std::fixed << std::setprecision(3)
         << std::right << std::setw(11) << g
         << std::setw(11) << n
         << std::setw(11) << c
         << std::setw(11) << n/g << "x"
         << std::setw(11) << c/g << "x"
         << std::setw(11) << n/c << "x";
    LOGGING_PRINTLN(line.str())
}

void measure(unsigned bits,unsigned iterations,unsigned samples)
{
    IntegerPool<IntegerOperations> native(bits);
    IntegerPool<GmpIntegerOperations> gmp(bits);
    verify(bits,native,gmp);
    const auto ix=[](unsigned i)->unsigned { return i%32u; };
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN(bits<<" bits | 32 GMP-verified OFFICIAL integer/division input pairs"
                   <<" | iterations="<<iterations<<" | samples="<<samples)
    LOGGING_PRINTLN(std::left<<std::setw(21)<<"operation"
                   <<std::right<<std::setw(11)<<"GMP ns"
                   <<std::setw(11)<<"Integer ns"
                   <<std::setw(11)<<"abs core ns"
                   <<std::setw(12)<<"Integer/GMP"
                   <<std::setw(12)<<"core/GMP"
                   <<std::setw(12)<<"Integer/core")

    report_row("construct_decimal",
        [&](unsigned i) {
            return GmpIntegerOperations::parse(gmp.spec[ix(i)].a);
        },[&](unsigned i) {
            return IntegerOperations::parse(native.spec[ix(i)].a);
        },[&](unsigned i) {
            return from_decimal(abs_digits(native.spec[ix(i)].a));
        },iterations,samples);

    report_row("copy",
        [&](unsigned i) { return GmpIntegerOperations::copy(gmp.a[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::copy(native.a[ix(i)]); },
        [&](unsigned i) { return BigUnsigned(native.a[ix(i)].magnitude()); },
        iterations,samples);

    report_row("compare",
        [&](unsigned i) { return GmpIntegerOperations::compare(gmp.a[ix(i)],gmp.b[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::compare(native.a[ix(i)],native.b[ix(i)]); },
        [&](unsigned i) { return compare(native.a[ix(i)].magnitude(),native.b[ix(i)].magnitude()); },
        iterations,samples);

    report_row("mul",
        [&](unsigned i) { return GmpIntegerOperations::mul(gmp.a[ix(i)],gmp.b[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::mul(native.a[ix(i)],native.b[ix(i)]); },
        [&](unsigned i) { return mul(native.a[ix(i)].magnitude(),native.b[ix(i)].magnitude()); },
        iterations,samples);

    report_row("square",
        [&](unsigned i) { return GmpIntegerOperations::square(gmp.a[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::square(native.a[ix(i)]); },
        [&](unsigned i) { return square(native.a[ix(i)].magnitude()); },
        iterations,samples);

    report_row("quotient",
        [&](unsigned i) { return GmpIntegerOperations::quotient(gmp.division_a[ix(i)],gmp.division_b[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::quotient(native.division_a[ix(i)],native.division_b[ix(i)]); },
        [&](unsigned i) { return quot(native.division_a[ix(i)].magnitude(),native.division_b[ix(i)].magnitude()); },
        iterations,samples);

    report_row("remainder",
        [&](unsigned i) { return GmpIntegerOperations::remainder(gmp.division_a[ix(i)],gmp.division_b[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::remainder(native.division_a[ix(i)],native.division_b[ix(i)]); },
        [&](unsigned i) { return rem(native.division_a[ix(i)].magnitude(),native.division_b[ix(i)].magnitude()); },
        iterations,samples);

    report_row("pow3",
        [&](unsigned i) { return GmpIntegerOperations::power3(gmp.a[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::power3(native.a[ix(i)]); },
        [&](unsigned i) {
            auto const& a=native.a[ix(i)].magnitude();
            return mul(square(a),a);
        },iterations,samples);

    report_row("to_decimal",
        [&](unsigned i) { return GmpIntegerOperations::to_decimal(gmp.a[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::to_decimal(native.a[ix(i)]); },
        [&](unsigned i) { return to_decimal(native.a[ix(i)].magnitude()); },
        iterations,samples);

    verify(bits,native,gmp);
    verify_signed_diagnostics(bits,native,gmp);
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN("B68 SAME SIGNED SEMANTICS: independent owning A/B"
                    " comparisons | "+std::to_string(bits)+" bits")
    LOGGING_PRINTLN(std::left << std::setw(22) << "operation"
                    << std::right << std::setw(12) << "GMP ns"
                    << std::setw(14) << "production ns"
                    << std::setw(14) << "candidate ns"
                    << std::setw(15) << "candidate/prod"
                    << std::setw(16) << "candidate/GMP")
    report_signed_candidate("construct_signed",
        [&](unsigned i) {
            return GmpIntegerOperations::parse(gmp.spec[ix(i)].a);
        },[&](unsigned i) {
            return IntegerOperations::parse(native.spec[ix(i)].a);
        },[&](unsigned i) {
            return parse_direct_signed_diagnostic(native.spec[ix(i)].a);
        },iterations,samples);
    report_signed_candidate("to_decimal_signed",
        [&](unsigned i) {
            return GmpIntegerOperations::to_decimal(gmp.a[ix(i)]);
        },[&](unsigned i) {
            return IntegerOperations::to_decimal(native.a[ix(i)]);
        },[&](unsigned i) {
            return format_direct_signed_diagnostic(native.a[ix(i)]);
        },iterations,samples);
    report_signed_candidate("mul_signed",
        [&](unsigned i) {
            return GmpIntegerOperations::mul(gmp.a[ix(i)],gmp.b[ix(i)]);
        },[&](unsigned i) {
            return IntegerOperations::mul(native.a[ix(i)],native.b[ix(i)]);
        },[&](unsigned i) {
            return multiply_direct_signed_diagnostic(
                native.a[ix(i)],native.b[ix(i)]);
        },iterations,samples);
    report_signed_candidate("square_signed",
        [&](unsigned i) {
            return GmpIntegerOperations::square(gmp.a[ix(i)]);
        },[&](unsigned i) {
            return IntegerOperations::square(native.a[ix(i)]);
        },[&](unsigned i) {
            return square_direct_signed_diagnostic(native.a[ix(i)]);
        },iterations,samples);
    verify_signed_diagnostics(bits,native,gmp);
    LOGGING_PRINTLN("B68: every candidate preserves the signed GMP result. "
                    "Rows are timed independently, not additive. "
                    "A/B timings can still contain host drift.")
    LOGGING_PRINTLN("abs core = unsigned production BigUnsigned operations "
                    "on the SAME magnitudes. Signed overhead cannot "
                    "be inferred as an exact difference because "
                    "signed controls can take different branches.")
    LOGGING_PRINTLN("All rows timed independently; do not sum components.")
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
            "numeric-integer-small-operations-profile",
            exact_arithmetic_parsers(selection,false)));
        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }
        for (const unsigned bits:selection.sizes) {
            measure(bits,selection.iterations,selection.samples);
        }
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("Numerical verification PASS: all magnitude core "
                        "operations AND B68 signed candidates "
                        "versus independent GMP before/after timing.")
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"B67 owning-vs-core diagnostic FAILED: "
                 <<error.what()<<'\n';
        return 1;
    }
}
