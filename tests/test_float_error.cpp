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

#include "numeric/float_error.hpp"
#include "numeric/float_error.tpl.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/positive.hpp"

#include "utility/test.hpp"

using namespace Ariadne;

namespace {

template<class F>
constexpr bool check_concept()
{
    using PR=typename F::PrecisionType;
    using E=Error<F>;
    return requires(PR pr, Nat m, F f, E x, E y) {
        E(pr); E(f); E(m,pr); E(ExactDouble(1),pr); E(TwoExp(0),pr);
        x=m;
        x+y; x*y; x+=y; x*=y;
        +x; -x; pm(x);
        nul(x); pos(x); neg(x); add(x,y); mul(x,y); sqr(x); pow(x,m); fma(x,y,x);
        sqrt(x); exp(x); log(x); log2(x);
        max(x,y); min(x,y); abs(x); mag(x);
        same(x,y); refines(x,y); refinement(x,y); cast_positive(x);
        x.precision(); x.characteristics(); x.raw(); x.generic();
    };
}

static_assert(check_concept<FloatDP>());
static_assert(check_concept<FloatMP>());

template<class F>
void test_float_error(typename F::PrecisionType pr)
{
    using E=Error<F>;

    E zero(pr);
    E one(1u,pr);
    E two_error(2u,pr);

    ARIADNE_TEST_EQUALS(zero.raw(),F(0,pr));
    ARIADNE_TEST_EQUALS(one.raw(),F(1,pr));
    ARIADNE_TEST_EQUALS(two_error.precision(),pr);
    ARIADNE_TEST_EQUALS(two_error.characteristics(),pr);

    Positive<F> positive_two(F(2,pr));
    E from_positive(positive_two);
    ARIADNE_TEST_EQUALS(from_positive.raw(),F(2,pr));

    UpperBound<F> upper_two(F(2,pr));
    E from_upper(upper_two);
    ARIADNE_TEST_EQUALS(from_upper.raw(),F(2,pr));

    PositiveBounds<F> positive_bounds(F(1,pr),F(2,pr));
    E from_bounds(positive_bounds);
    ARIADNE_TEST_EQUALS(from_bounds.raw(),F(2,pr));

    ValidatedUpperNumber generic_upper=upper_two;
    E from_generic_upper(generic_upper,pr);
    ARIADNE_TEST_EQUALS(from_generic_upper.raw(),F(2,pr));

    PositiveUpperBound<F>& mutable_upper=from_upper;
    ARIADNE_TEST_EQUALS(mutable_upper.raw(),F(2,pr));
    E const& const_from_upper=from_upper;
    PositiveUpperBound<F> const& const_upper=const_from_upper;
    ARIADNE_TEST_EQUALS(const_upper.raw(),F(2,pr));

    E from_exact(ExactDouble(2),pr);
    E from_twoexp(TwoExp(1),pr);
    ARIADNE_TEST_EQUALS(from_exact.raw(),F(2,pr));
    ARIADNE_TEST_EQUALS(from_twoexp.raw(),F(2,pr));

    ARIADNE_TEST_EQUALS((one+two_error).raw(),F(3,pr));
    ARIADNE_TEST_EQUALS((two_error*two_error).raw(),F(4,pr));

    E accumulated=one;
    accumulated+=two_error;
    ARIADNE_TEST_EQUALS(accumulated.raw(),F(3,pr));
    accumulated*=two_error;
    ARIADNE_TEST_EQUALS(accumulated.raw(),F(6,pr));

    ARIADNE_TEST_EQUALS((+two_error).raw(),F(2,pr));
    ARIADNE_TEST_EQUALS((-two_error).raw(),F(-2,pr));

    auto symmetric=pm(two_error);
    ARIADNE_TEST_EQUALS(symmetric.lower_raw(),F(-2,pr));
    ARIADNE_TEST_EQUALS(symmetric.upper_raw(),F(2,pr));

    ARIADNE_TEST_EQUALS(nul(two_error).raw(),F(0,pr));
    ARIADNE_TEST_EQUALS(pos(two_error).raw(),F(2,pr));
    ARIADNE_TEST_EQUALS(neg(two_error).raw(),F(-2,pr));
    ARIADNE_TEST_EQUALS(add(one,two_error).raw(),F(3,pr));
    ARIADNE_TEST_EQUALS(mul(two_error,two_error).raw(),F(4,pr));
    ARIADNE_TEST_EQUALS(sqr(two_error).raw(),F(4,pr));
    ARIADNE_TEST_EQUALS(pow(two_error,Nat(3u)).raw(),F(8,pr));
    ARIADNE_TEST_EQUALS(fma(two_error,two_error,one).raw(),F(5,pr));

    ARIADNE_TEST_EQUALS(sqrt(E(4u,pr)).raw(),F(2,pr));
    ARIADNE_TEST_EQUALS(exp(zero).raw(),F(1,pr));
    ARIADNE_TEST_EQUALS(log(one).raw(),F(0,pr));

    ARIADNE_TEST_EQUALS(max(one,two_error).raw(),F(2,pr));
    ARIADNE_TEST_EQUALS(min(one,two_error).raw(),F(1,pr));
    ARIADNE_TEST_EQUALS(abs(two_error).raw(),F(2,pr));
    ARIADNE_TEST_EQUALS(mag(two_error).raw(),F(2,pr));

    ARIADNE_TEST_EQUALS((two_error+Nat(3u)).raw(),F(5,pr));
    ARIADNE_TEST_EQUALS((Nat(3u)+two_error).raw(),F(5,pr));
    ARIADNE_TEST_EQUALS((two_error*Nat(3u)).raw(),F(6,pr));
    ARIADNE_TEST_EQUALS((Nat(3u)*two_error).raw(),F(6,pr));
    ARIADNE_TEST_EQUALS((two_error/Nat(2u)).raw(),F(1,pr));

    ARIADNE_TEST_ASSERT(same(one,E(1u,pr)));
    ARIADNE_TEST_ASSERT(not same(one,two_error));
    ARIADNE_TEST_ASSERT(refines(one,two_error));
    ARIADNE_TEST_EQUALS(refinement(one,two_error).raw(),F(1,pr));

    ValidatedErrorNumber generic=two_error.generic();
    E assigned(pr);
    assigned=generic;
    ARIADNE_TEST_EQUALS(assigned.raw(),F(2,pr));

    StringStream stream;
    Operations<E>::_write(stream,two_error);
    ARIADNE_TEST_ASSERT(not stream.str().empty());

    StringStream direct_output;
    direct_output << two_error;
    ARIADNE_TEST_ASSERT(not direct_output.str().empty());

    StringStream positive_input("2");
    E parsed(pr);
    Operations<E>::_read(positive_input,parsed);
    ARIADNE_TEST_EQUALS(parsed.raw(),F(2,pr));

    StringStream direct_input("3");
    direct_input >> parsed;
    ARIADNE_TEST_EQUALS(parsed.raw(),F(3,pr));

    StringStream negative_input("-1");
    ARIADNE_TEST_FAIL(Operations<E>::_read(negative_input,parsed));

    StringStream direct_negative_input("-1");
    ARIADNE_TEST_FAIL(direct_negative_input >> parsed);

    StringStream failed_input;
    failed_input.setstate(std::ios::badbit);
    Operations<E>::_read(failed_input,parsed);
    ARIADNE_TEST_ASSERT(failed_input.bad());

    StringStream direct_failed_input;
    direct_failed_input.setstate(std::ios::badbit);
    direct_failed_input >> parsed;
    ARIADNE_TEST_ASSERT(direct_failed_input.bad());

    if constexpr (Same<F,FloatDP>) {
        ARIADNE_TEST_EQUALS(class_name<E>(),String("FloatDPError"));
    } else {
        ARIADNE_TEST_EQUALS(class_name<E>(),String("FloatMPError"));
    }
}

} // namespace

int main()
{
    test_float_error<FloatDP>(dp);
    test_float_error<FloatMP>(MultiplePrecision(128_bits));
    return ARIADNE_TEST_FAILURES;
}
