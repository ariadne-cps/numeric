/***************************************************************************
 *            test_float_error.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

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

#include <sstream>

#include "numeric/float_bounds.hpp"
#include "numeric/float_error.hpp"
#include "numeric/float_error.tpl.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/positive.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

template<class F>
constexpr bool check_concept() {
    using PR = typename F::PrecisionType;
    using E = Error<F>;
    return requires(PR pr, Nat m, F f, E x, E y) {
        E(pr);
        E(f);
        E(m, pr);
        E(ExactDouble(1), pr);
        E(TwoExp(0), pr);
        x = m;
        x + y;
        x * y;
        x += y;
        x *= y;
        +x;
        -x;
        pm(x);
        nul(x);
        pos(x);
        neg(x);
        add(x, y);
        mul(x, y);
        sqr(x);
        pow(x, m);
        fma(x, y, x);
        sqrt(x);
        exp(x);
        log(x);
        log2(x);
        max(x, y);
        min(x, y);
        abs(x);
        mag(x);
        same(x, y);
        refines(x, y);
        refinement(x, y);
        cast_positive(x);
        x.precision();
        x.characteristics();
        x.raw();
        x.generic();
    };
}

static_assert(check_concept<FloatDP>());
static_assert(check_concept<FloatMP>());

template<class F>
class TestFloatError : public Ariadne::Testing::ContextFixture {
    typename F::PrecisionType precision;

  public:
    TestFloatError(Ariadne::Testing::Context& ctx, typename F::PrecisionType pr)
        : ContextFixture(ctx), precision(pr) {}
    void test_error();
};

template<class F>
void TestFloatError<F>::test_error() {
    const auto& pr = precision;
    using E = Error<F>;

    Nat saved_output_places = E::output_places;
    E::set_output_places(4u);
    ARIADNE_TESTING_EXPECT_EQUAL(E::output_places, Nat(4u));
    E::set_output_places(saved_output_places);

    E zero(pr);
    E one(1u, pr);
    E two_error(2u, pr);

    ARIADNE_TESTING_EXPECT_EQUAL(zero.raw(), F(0, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(one.raw(), F(1, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(two_error.precision(), pr);
    ARIADNE_TESTING_EXPECT_EQUAL(two_error.characteristics(), pr);

    Positive<F> positive_two(F(2, pr));
    E from_positive(positive_two);
    ARIADNE_TESTING_EXPECT_EQUAL(from_positive.raw(), F(2, pr));

    UpperBound<F> upper_two(F(2, pr));
    E from_upper(upper_two);
    ARIADNE_TESTING_EXPECT_EQUAL(from_upper.raw(), F(2, pr));

    ARIADNE_TESTING_EXPECT_EXCEPTION(E(F(-1, pr)));
    UpperBound<F> negative_upper(F(-1, pr));
    ARIADNE_TESTING_EXPECT_EXCEPTION(E{negative_upper});

    PositiveBounds<F> positive_bounds(F(1, pr), F(2, pr));
    E from_bounds(positive_bounds);
    ARIADNE_TESTING_EXPECT_EQUAL(from_bounds.raw(), F(2, pr));

    ValidatedUpperNumber generic_upper = upper_two;
    E from_generic_upper(generic_upper, pr);
    ARIADNE_TESTING_EXPECT_EQUAL(from_generic_upper.raw(), F(2, pr));

    UpperBound<F> converted_upper = from_upper;
    ARIADNE_TESTING_EXPECT_EQUAL(converted_upper.raw(), F(2, pr));

    E assigned_natural(pr);
    assigned_natural = Nat(2u);
    ARIADNE_TESTING_EXPECT_EQUAL(assigned_natural.raw(), F(2, pr));

    E from_exact(ExactDouble(2), pr);
    E from_twoexp(TwoExp(1), pr);
    ARIADNE_TESTING_EXPECT_EQUAL(from_exact.raw(), F(2, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(from_twoexp.raw(), F(2, pr));

    ARIADNE_TESTING_EXPECT_EQUAL((one + two_error).raw(), F(3, pr));
    ARIADNE_TESTING_EXPECT_EQUAL((two_error * two_error).raw(), F(4, pr));

    E accumulated = one;
    accumulated += two_error;
    ARIADNE_TESTING_EXPECT_EQUAL(accumulated.raw(), F(3, pr));
    accumulated *= two_error;
    ARIADNE_TESTING_EXPECT_EQUAL(accumulated.raw(), F(6, pr));

    ARIADNE_TESTING_EXPECT_EQUAL((+two_error).raw(), F(2, pr));
    ARIADNE_TESTING_EXPECT_EQUAL((-two_error).raw(), F(-2, pr));

    auto symmetric = pm(two_error);
    ARIADNE_TESTING_EXPECT_EQUAL(symmetric.lower_raw(), F(-2, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(symmetric.upper_raw(), F(2, pr));

    ARIADNE_TESTING_EXPECT_EQUAL(nul(two_error).raw(), F(0, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(pos(two_error).raw(), F(2, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(neg(two_error).raw(), F(-2, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(add(one, two_error).raw(), F(3, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(mul(two_error, two_error).raw(), F(4, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(sqr(two_error).raw(), F(4, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(pow(two_error, Nat(3u)).raw(), F(8, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(fma(two_error, two_error, one).raw(), F(5, pr));

    ARIADNE_TESTING_EXPECT_EQUAL(sqrt(E(4u, pr)).raw(), F(2, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(exp(zero).raw(), F(1, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(log(one).raw(), F(0, pr));

    ARIADNE_TESTING_EXPECT_EQUAL(max(one, two_error).raw(), F(2, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(min(one, two_error).raw(), F(1, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(abs(two_error).raw(), F(2, pr));
    ARIADNE_TESTING_EXPECT_EQUAL(mag(two_error).raw(), F(2, pr));

    ARIADNE_TESTING_EXPECT_EQUAL((two_error + Nat(3u)).raw(), F(5, pr));
    ARIADNE_TESTING_EXPECT_EQUAL((Nat(3u) + two_error).raw(), F(5, pr));
    ARIADNE_TESTING_EXPECT_EQUAL((two_error * Nat(3u)).raw(), F(6, pr));
    ARIADNE_TESTING_EXPECT_EQUAL((Nat(3u) * two_error).raw(), F(6, pr));
    ARIADNE_TESTING_EXPECT_EQUAL((two_error / Nat(2u)).raw(), F(1, pr));

    ARIADNE_TESTING_PROVE(same(one, E(1u, pr)));
    ARIADNE_TESTING_PROVE(not same(one, two_error));
    ARIADNE_TESTING_PROVE(refines(one, two_error));
    ARIADNE_TESTING_EXPECT_EQUAL(refinement(one, two_error).raw(), F(1, pr));

    ValidatedErrorNumber generic = two_error.generic();
    E assigned(pr);
    assigned = generic;
    ARIADNE_TESTING_EXPECT_EQUAL(assigned.raw(), F(2, pr));
    if constexpr (Same<F, FloatMP>) {
        // Native error conversion must retain the FloatMP backend.
        ARIADNE_TESTING_EXPECT_EQUAL(generic.class_name(), String("FloatMP"));
    }

    StringStream stream;
    Operations<E>::_write(stream, two_error);
    ARIADNE_TESTING_PROVE(not stream.str().empty());

    StringStream direct_output;
    direct_output << two_error;
    ARIADNE_TESTING_PROVE(not direct_output.str().empty());

    StringStream positive_input("2");
    E parsed(pr);
    Operations<E>::_read(positive_input, parsed);
    ARIADNE_TESTING_EXPECT_EQUAL(parsed.raw(), F(2, pr));

    StringStream direct_input("3");
    direct_input >> parsed;
    ARIADNE_TESTING_EXPECT_EQUAL(parsed.raw(), F(3, pr));

    StringStream negative_input("-1");
    ARIADNE_TESTING_EXPECT_EXCEPTION(Operations<E>::_read(negative_input, parsed));

    StringStream direct_negative_input("-1");
    ARIADNE_TESTING_EXPECT_EXCEPTION(direct_negative_input >> parsed);

    StringStream failed_input;
    failed_input.setstate(std::ios::badbit);
    Operations<E>::_read(failed_input, parsed);
    ARIADNE_TESTING_PROVE(failed_input.bad());

    StringStream direct_failed_input;
    direct_failed_input.setstate(std::ios::badbit);
    direct_failed_input >> parsed;
    ARIADNE_TESTING_PROVE(direct_failed_input.bad());

    if constexpr (Same<F, FloatDP>) {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<E>(), String("FloatDPError"));
    } else {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<E>(), String("FloatMPError"));
    }
}

} // namespace

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("float_error");
    suite.add("FloatDP", [](Ariadne::Testing::Context& ctx) {
        TestFloatError<FloatDP>(ctx, dp).test_error();
    });
    suite.add("FloatMP-128", [](Ariadne::Testing::Context& ctx) {
        TestFloatError<FloatMP>(ctx, MultiplePrecision(128_bits)).test_error();
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
