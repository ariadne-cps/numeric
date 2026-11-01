/***************************************************************************
 *            test_logical_sequence.cpp
 *
 *  Copyright  2026  Ariadne contributors
 *
 ****************************************************************************/
#include "foundation/logical.hpp"
#include "numeric/integer.hpp"
#include "numeric/sequence.hpp"
#include "numeric/logical_sequence.hpp"
#include "ariadne-testing.hpp"
#include <iostream>
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
}
