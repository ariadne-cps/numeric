/***************************************************************************
 *            test_logical_sequence.cpp
 *
 *  Copyright  2026  Ariadne contributors
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric. If not, see <https://www.gnu.org/licenses/>.
 */

#include "foundation/logical.hpp"
#include "numeric/integer.hpp"
#include "numeric/sequence.hpp"
#include "numeric/logical_sequence.hpp"
#include "ariadne-testing.hpp"
#include <iostream>
#include <sstream>
using namespace Ariadne;
class TestLogicalSequence : public Ariadne::Testing::ContextFixture { public: using ContextFixture::ContextFixture; Void test(); };
int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("logical_sequence");
    suite.add("disjunction and conjunction", [](Ariadne::Testing::Context& context) {
        TestLogicalSequence(context).test();
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
Void TestLogicalSequence::test() {
    Sequence<LowerKleenean> seq([](Natural n){return n==2 ? LowerKleenean(true) : LowerKleenean(indeterminate);});
    LowerKleenean some=disjunction(seq);
    ARIADNE_TESTING_EXPECT_POSSIBLY( not some.check(2_eff));
    ARIADNE_TESTING_EXPECT_DEFINITELY( some.check(3_eff));
    ARIADNE_TESTING_EXPECT_DEFINITELY( some.check(4_eff));
    LogicalInterface* some_copy=some.repr().pointer()->_copy();
    ARIADNE_TESTING_PROVE( some_copy->_check(3_eff)==LogicalValue::TRUE);
    delete some_copy;
    Sequence<UpperKleenean> upper_seq([](Natural n){return n==2 ? UpperKleenean(false) : UpperKleenean(indeterminate);});
    UpperKleenean all=conjunction(upper_seq);
    ARIADNE_TESTING_EXPECT_POSSIBLY( all.check(2_eff));
    ARIADNE_TESTING_EXPECT_NOT_POSSIBLY( all.check(3_eff));
    ARIADNE_TESTING_EXPECT_DEFINITELY( not all.check(4_eff));
    LogicalInterface* all_copy=all.repr().pointer()->_copy();
    ARIADNE_TESTING_PROVE( all_copy->_check(3_eff)==LogicalValue::FALSE);
    delete all_copy;

    // Verify the observable representation of the two lazy sequences.
    // This also exercises the polymorphic _write implementations.
    std::ostringstream disjunction_output;
    some.repr().pointer()->_write(disjunction_output);
    std::ostringstream expected_disjunction;
    expected_disjunction << "disjunction(" << seq[0u] << "," << seq[1u]
                         << "," << seq[2u] << ",...)";
    ARIADNE_TESTING_EXPECT_EQUAL(disjunction_output.str(),expected_disjunction.str());

    std::ostringstream conjunction_output;
    all.repr().pointer()->_write(conjunction_output);
    std::ostringstream expected_conjunction;
    expected_conjunction << "conjunction(" << upper_seq[0u] << ","
                         << upper_seq[1u] << "," << upper_seq[2u] << ",...)";
    ARIADNE_TESTING_EXPECT_EQUAL(conjunction_output.str(),expected_conjunction.str());
}
