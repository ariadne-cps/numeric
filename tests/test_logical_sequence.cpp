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
#include "utility/test.hpp"
using namespace Ariadne;
class TestLogicalSequence { public: Void test(); };
Int main() {
    ARIADNE_TEST_CLASS(TestLogicalSequence,TestLogicalSequence());
    return ARIADNE_TEST_FAILURES;
}
Void TestLogicalSequence::test() {
    Sequence<LowerKleenean> seq([](Natural n){return n==2 ? LowerKleenean(true) : LowerKleenean(indeterminate);});
    ARIADNE_TEST_ASSIGN_CONSTRUCT(LowerKleenean, some, disjunction(seq));
    ARIADNE_TEST_ASSERT(possibly(not some.check(2_eff)));
    ARIADNE_TEST_ASSERT(definitely(some.check(3_eff)));
    ARIADNE_TEST_ASSERT(definitely(some.check(4_eff)));
    LogicalInterface* some_copy=some.repr().pointer()->_copy();
    ARIADNE_TEST_ASSERT(some_copy->_check(3_eff)==LogicalValue::TRUE);
    delete some_copy;
    Sequence<UpperKleenean> upper_seq([](Natural n){return n==2 ? UpperKleenean(false) : UpperKleenean(indeterminate);});
    ARIADNE_TEST_ASSIGN_CONSTRUCT(UpperKleenean, all, conjunction(upper_seq));
    ARIADNE_TEST_ASSERT(possibly(all.check(2_eff)));
    ARIADNE_TEST_ASSERT(not possibly(all.check(3_eff)));
    ARIADNE_TEST_ASSERT(definitely(not all.check(4_eff)));
    LogicalInterface* all_copy=all.repr().pointer()->_copy();
    ARIADNE_TEST_ASSERT(all_copy->_check(3_eff)==LogicalValue::FALSE);
    delete all_copy;
}
