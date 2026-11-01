# Migration of Numeric tests to Ariadne Testing

2026-11-01: `test_operators` and `test_complex` migrated on branch `testing`.

- Both test binaries now parse `--concise` using `Ariadne::Testing::run_suite(suite, argc, argv)`.
- `test_operators` retains its five functional groups as independent cases; its metadata helper records two assertions per invocation. Execution-only operations become explicit no-throw checks, and the invalid operator-construction case still requires an exception.
- `test_complex` retains the six original fixture groups plus Fourier computation, with fresh fixtures per group and the original FloatDPBounds precision. The code continues to avoid enabling commented-out tests.
- Existing runtime assertions are migrated to evidence-aware checks; any inconclusive evidence is now explicitly classified rather than coerced into a boolean. This semantic difference is intentional and must be checked by CI.
- The new suite macros record failed checks without disabling subsequent checks. Errors thrown from one group no longer suppress unrelated groups.

CMake links the two migrated executables with `ariadne-testing`; all other targets retain their prior configuration and CTest arguments. Legacy test binaries may still ignore `--concise` until migrated.

**Validation status:** The source transformation and CMake registration have been checked, but a C++ build or CTest run was not possible in this environment because GitHub cloning/network access is unavailable from the build container. Treat compilation and the Numeric CI matrix as a mandatory acceptance gate; do not interpret the migration as verified until those checks pass.

## 2026-11-01: Float ball and directed floating-point bounds

- Migrated `test_float_ball` (3 precision configurations x 4 groups = 12 cases) and `test_float_lower_upper_bound` (2 precision configurations x 6 groups = 12 cases) to Ariadne Testing's Context-based suites.
- Each original method is registered in its original order, with fresh precision-specific fixtures. Existing assertions, precision choices, exception expectations and numeric predicates remain present. Binary predicates use explicit evidence-preserving macro checks; relational comparisons map to their matching relation checks.
- Direct invocation defaults to full console traces. CTest's pre-existing `--concise` argument now takes effect for these two targets. Both link against `ariadne-testing` rather than solely including the legacy macro headers.
- Validation pending on the new commit's CI matrix; the previous green CI results cover `test_operators` and `test_complex` only.

## 2026-11-01: Temporary crash diagnostics for directed bounds

The first CI of `bf0c061` passes Release but reports SIGSEGV for
`test_float_lower_upper_bound` in Debug. Concise intentionally emits
nothing until `suite_finished`, which a SIGSEGV bypasses. Temporarily
register ONLY this executable with `--full` in CTest so that the
flushed `[RUN]` and `[CHECK]` records identify the failing case and
last reached assertion in CI's `--output-on-failure` transcript. Revert
to `--concise` once the root cause is fixed and tested across platforms.
This is a diagnostic change, not a relaxation or deletion of assertions.

## 2026-11-01: Restore concise after crash breadcrumb support

The temporary `--full` CTest exception for `test_float_lower_upper_bound`
is removed. The Testing submodule now emits fatal-signal-only current-case
and last-CHECK breadcrumbs in `--concise`; all normal runs retain the
summary-only contract. The previously observed Debug SIGSEGV is not itself
fixed by this change and remains a separate regression to investigate.

## 2026-11-01: Fix MSVC generic fixture method registration

MSVC rejected the explicitly spelled dependent pointer-to-member parameter
`void (Fixture::*method)()` in the two newly migrated precision-parametric
suites (`test_float_ball` and `test_float_lower_upper_bound`) with C2653 /
C2091 and cascading C2664 failures. Both local add-case lambdas now use
`auto method`, allowing the compiler to deduce the exact member pointer
from `&Fixture::test_...`. Invocation through `(fixture.*method)()`,
registration order, precision variants, test bodies and assertions are
unchanged. The generic-lambda shape was syntax-checked with GCC and Clang
under C++20 -Wall -Wextra -Wpedantic -Werror. MSVC acceptance and downstream
runtime verification remain CI gates.

## 2026-11-01: FloatBounds precision suites

Migrated `test_float_bounds` to Ariadne Testing: 16 existing functional groups
registered in their original order for FloatDP, FloatMP at 53 bits and FloatMP
at 128 bits (48 independent cases). The fixture preserves the precision choice,
all numeric oracles, compile-time contracts and all active assertions. Existing
`COMPARE` assertions now use the corresponding relation-specific checks; this
requires the short-lived operand snapshot support already pinned in `ss/tst`.
The previous `SAME` predicates remain identity checks and are not weakened
to equality. Eight print-only diagnostics are retained as typed inspection
events; three warnings retain their warning severity. One-line checks missing
trailing semicolons under the old block macros now have explicit terminators.
The migrated executable accepts `--concise` and is linked to `ariadne-testing`.
Windows member-method registration uses generic-lambda deduction to avoid the
MSVC pointer-to-member parsing failure previously seen in the ball/directed
suites. All three original precision instantiations remain enabled.

This is a source-level migration, not a completed build validation: verify the
numeric Ubuntu, Windows, macOS, Sanitizers and Coverage CI matrix before
accepting it, and classify any failure rather than modifying numerical oracles.

## 2026-11-01: Dyadic suite migration

Migrated `test_dyadic` to Ariadne Testing: the original eight fixture
methods and five standalone TSV/contract oracles remain registered in exactly
their original order (13 cases). The original mathematical checks (159
equality comparisons, 99 propositions, 9 binary predicates and 15 expected
exceptions), compile-time concept check and file row-count/completeness
conditions are retained. A nonasserting constructor smoke check still builds
`Dyadic q(3.25_q2)` before the value checks.

Oracle row identities have been moved from direct stdout lines to scoped
Testing diagnostic contexts attached to every row assertion, providing useful
row information in concise failure summaries without printing thousands of
progress lines during CTest. The original 25 execution-only operations retain
their no-throw checks; 15 ostream streaming operations use a local
`std::ostringstream` instead of writing directly to stdout, with each
rendered output captured through an inspection event (the seven original
PRINT diagnostics are also retained). Thus default/full console output remains
diagnostic while `--concise` stays silent on success.

The test maintains the original default-writer state transitions and test
group ordering. Existing oracle values, numerical expectations, required
exceptions, and precision choices have not been relaxed. The executable
uses `run_suite(suite, argc, argv)` and links with `ariadne-testing`.
CI validation on all five platforms/workflows remains the acceptance gate.

## 2026-11-01: Repair explicit contexts in free-standing Dyadic oracle cases

The initial Dyadic migration mistakenly left the fixture-only assertion
shorthands in five standalone functions that accept `Testing::Context& context`.
Those shorthands expand to `this->ctx()` and cannot compile outside a member
function. Change every such assertion to its explicit `_CTX(context, ...)`
form: 41 equalities, 17 propositions, and 8 expected-exception checks.
All original predicates, expected values, exception requirements, case order,
and TSV completeness oracles are unchanged. Fixture methods continue using
the appropriate fixture shorthands. Acceptance requires rerunning the full
Numeric CI matrix on the new commit; a source audit alone is insufficient.

## 2026-11-01: Integer suite migration

Migrated `test_integer` from the legacy macro runner to Ariadne Testing on
branch `testing`, preserving the ten method groups and their original order:
constructors, comparisons, literals, arithmetic, small types, combinatorics,
native operations, signed regression, owned results, and integer square root.
The standalone test now accepts `--concise` through
`run_suite(suite, argc, argv)` and links with `ariadne-testing`.

All original computational expressions and oracle expected values are
unchanged. Thirteen legacy constructor macros now declare the same named
`Integer` variables in the original scope and emit typed construction
diagnostics. Twelve legacy relational comparisons use the exact matching
relation macro (including strict vs non-strict ordering and inequality).
Existing assertions and expected-exception checks remain active, and the
integer-conversion TSV retains its row completeness checks; each row now
adds an assertion scope for diagnosing failures without a successful-run log
flood. Compile-time concept checks remain unchanged.

Static inspection verified the absence of active legacy macros and checked
terminating punctuation for new assertions. No C++ binary was built locally:
acceptance requires the full Numeric CI matrix on the new commit. If any
test fails, classify the regression as test migration, test expectation,
or Numeric behavior rather than weakening or deleting the oracle.

## 2026-11-01: Rounded float suite migration

Migrated `test_rounded_float` on branch `testing`, retaining its original
sequence: numeric cast checks, then nine rounded-float method groups for
FloatDP, then the same nine groups for FloatMP at 128-bit precision (19
Ariadne Testing cases). The fixture method declarations, rounding-mode
transitions, compile-time concepts, certified oracle TSV data and exact
directed comparison thresholds are retained. Each case receives the same
precision as before. The test suite executes cases serially and in the
original order; existing shared floating-point rounding-mode transitions are
not artificially reset between method groups.

All active legacy assertions are converted to matching Ariadne Testing
checks: 86 relational checks preserve their exact relation, 71 propositions
remain propositions, the 25 EQUALS and 20 EQUAL comparisons remain equality
checks, and four `models` predicates retain their predicate. The 40 legacy
execution-only operations become no-throw checks, 12 explicit constructor
declarations remain in their original scope, and three print inspections
become typed diagnostics. The 20 direct progress `cout` messages now emit
reporter notes rather than bypassing `--concise`; the two existing warnings
remain warnings. The former anonymous-namespace numeric-cast function
accepts an explicit `Testing::Context&`, and only its assertions use
`_CTX` macros.

CTest already supplies `--concise`; CMake now links the target with
`ariadne-testing`, and direct invocation defaults to full diagnostics.
Source-level macro/terminator audits were performed, but compilation and
the full Numeric CI matrix are still required for acceptance. Do not alter
oracle expectations to make the converted tests pass.

## 2026-11-01: Resolve systematic inconclusive Rounded tests

The first CI of `cfb7fc7` reported 18/19 cases as Inconclusive because
`Rounded<FloatDP>` and `Rounded<FloatMP>` comparisons return
`ApproximateKleenean`. Its `definitely` is never true and `possibly` is
never false, so a generic evidence-preserving assertion cannot classify
either a likely-true or an unlikely-false comparison. This classification
is deliberate in `ariadne-testing`; the generic framework is unchanged.

This test file now makes its historical *decision-based* semantics explicit:
legacy EQUAL/EQUALS and COMPARE assertions use local, lifetime-safe decided
relation checks, invoking ADL `decide` only when the predicate type is
`ApproximateKleenean`. Other logical types continue through the framework's
normal evidence classification. The local relation wrapper evaluates each
operand once and snapshots both operands' display representations before
temporaries expire. The original relation operators and numeric thresholds
are unchanged. The original ASSERT checks use a locally named approximate
decision check (and do not label a merely likely outcome as a mathematical
proof); unlike the previous legacy warning-only handling for approximate
ASSERTs, a false approximate decision now fails visibly, increasing test
sensitivity rather than silently accepting a dubious result.

The numeric-cast case now also checks that `ApproximateKleenean(true)` and
`ApproximateKleenean(false)` make opposite decisions while both retain
`Evidence::Indeterminate` under generic evidence classification. This
guards against incorrectly treating all indeterminate evidence as passing.
The suite still has 19 cases and all prior numerical oracles. No other
consumer or test framework source is changed.

Compilation of the test-local helper syntax under GCC/Clang C++20 was checked
using a minimal standalone mock; acceptance still requires full Numeric CI
for this exact commit on Ubuntu, macOS, Windows, Sanitizers, and Coverage.

## 2026-11-01: Rational test suite migration

Migrated `test_rational` to Ariadne Testing on `numeric/testing`, in the
original test order. Three preliminary cases (`integral_binary64_scale`,
`distinct_prefixes`, `contract_edges`) are followed by the independent
27-row scalar comparison TSV case, then the ten existing `TestRational`
method groups (14 cases total). The 10 fixture methods now inherit
`ContextFixture`; all four free functions accept explicit
`Testing::Context&` and use only `_CTX(context, ...)` assertions. The
original numeric assertions, operator predicates, exception checks and
expectations are unchanged, including comparisons involving infinities,
NaNs, binary64 directed rounding and frozen TSV oracles.

There are 230 equality expectations, 67 propositional assertions,
26 binary-predicate checks, and 12 exception expectations, all retained.
The single construction macro now creates the same named local `Rational`
and emits an inspection. Scalar comparison row identities are attached as
scoped diagnostics rather than printed unconditionally. Diagnostics for
R5 double-conversion mismatches are routed through `Context::warning`
instead of writing straight to stderr; the original conversion_matches
flag and final assertion remain unchanged. One existing source expression
lacked a trailing semicolon because the legacy macro was statement-shaped:
an explicit semicolon was inserted for the new macro.

CTest already passes `--concise`. CMake links this target to
`ariadne-testing`, and its direct invocation defaults to full reporting.
Source inspection confirmed that no active legacy assertion macro remains,
that no fixture-only assertion shorthand is used in the free functions,
and that all converted assertion invocations terminate correctly.
Full CI of this commit remains mandatory; do not relax oracle coverage
or numerical expectations to make the regression green.

## 2026-11-01: Unify Dyadic fixtures and repair Rational lambda captures

All five formerly free-standing Dyadic oracle/contract functions are now
`TestDyadic` fixture methods. This keeps the original order of all 13 cases
and their isolated fixture instances, while registering every method with
`ARIADNE_TESTING_ADD_TEST_TO(suite, TestDyadic, ...)`. Within these five
methods, all 66 checks (41 equalities, 17 propositions, 8 exception checks)
use the regular fixture form rather than `_CTX(context, ...)`; three
diagnostic TSV-row scopes now use `this->ctx().scope(...)`. No oracle row,
expected value, completeness condition, or computation was changed.
There is no need for free functions in this suite; free-standing tests are
still supported by the framework via explicit `_CTX` macros in general.

The preceding Rational migration (commit `c59dde8`) could not compile:
four local `[]` lambdas in `TestRational` contained checks that expand to
`this->ctx()`. The lambdas `check_result`, `verify_fraction`,
`check_bounds`, and `check_dyadic` now capture `[this]` explicitly.
The independent `parse_r16` lambda remains capture-free. The eight
affected expectation statements are unchanged; this is a capture-only
compilation fix, not a changed numeric oracle. CMake and framework sources
are unchanged.

Acceptance requires the complete Numeric CI matrix on this exact commit;
prior failed builds of `test_rational` do not qualify as validated.

## 2026-11-01: Consolidate runtime test bodies inside fixtures

As a coding convention, every runtime test body in the application repositories
belongs to a ContextFixture-derived class. Standalone functions are reserved
for true compile-time checks and oracle generator executables; the Testing
framework itself deliberately exercises its explicit-Context API separately.

This change moves the four standalone Rational cases into TestRational while
preserving their order ahead of the ten existing cases, complete TSV coverage
and all original comparisons. It introduces TestDecimal and TestExtended
fixtures for their two and three existing cases, respectively, and a templated
TestFloatError fixture for the original FloatDP and FloatMP parameterizations.
The numeric-casts case becomes TestNumericCasts; the Fourier case becomes a
method of TestComplex, retaining its expected-no-throw check and constructor
configuration. The templated operator metadata helper is moved into
TestOperators, where it uses fixture assertion macros.

Test-only registration helper functions remain outside the fixtures if they
only build parameterized Suite entries, not run assertions. Mathematical and
compile-time oracle values remain unchanged. All these edits require CI before
acceptance; preserving source expressions is not the same as verified runtime
equivalence.

## 2026-11-01: Legacy Numeric runtime helpers moved into test classes

All 32 standalone runtime test functions from `test_float` are now methods
of the internal `TestFloatStandalone` class; `main` invokes the same methods
at the same points, interleaved with the same FloatDP/FloatMP parametrized
fixture invocations. `test_big_unsigned_contracts_and_controls` is now
a `TestBigUnsigned` method, called first in the original test sequence.
Its eight oracle/helper functions and the `LimbDivisionGolden` reader type
are private members/nested helpers of that same class; the immutable oracle
path constants stay namespace-scope configuration.

In `test_number`, the four local model predicates are members of
`TestNumber<Y>`, and `test_number_get` with its ExactNumber and
ApproximateNumber specializations is now a class member. The original
template specialization behavior and call sites remain. `test_real`
moves the excessive-precision runtime preflight out of `main` and into
`TestReal` before all original test phases; the certified pi TSV reader
and its result struct are private class members. The `Ariadne` namespace's
literal operator remains a non-member because C++ requires user-defined
literal operators to be declared at namespace scope.

This commit is a structural refactor of the four yet-unmigrated executables,
not a conversion of their legacy runner/macro API. Original check macro counts
and numeric oracles were compared before and after the transformation:
BigUnsigned 265 ASSERT / 30 EQUALS / 16 FAIL, Real 39 ASSERT / 62 EQUALS
/ 3 FAIL / 7 WITHIN, Number 43 ASSERT / 28 EQUALS / 5 FAIL
/ 49 BINARY_PREDICATE / 4 WITHIN, and Float 289 ASSERT / 389 EQUALS
/ 39 FAIL / 28 BINARY_PREDICATE, all unchanged. Full CI is required for
acceptance; only the classes and call qualification have changed.

## 2026-11-01: Parameterized case registration belongs to fixture types

Removed the remaining free-standing run-time case registration helpers in
`test_rounded_float`, `test_float_lower_upper_bound`,
`test_float_bounds`, and `test_float_ball`. Each corresponding fixture
now declares a public static `add_cases(Suite&, variant, precision)`
method. The implementation preserves all exact original case labels,
precision variants, method pointers, construction semantics, and order.
`main` calls `Fixture<...>::add_cases(...)` rather than a free function.

Parameterized registration callbacks remain local to these methods, not
free-standing runtime test functions. Compile-time concept checks, frozen
oracle fixtures and generator programs are unaffected. CI is required for
acceptance; no numerical oracle or expected result was altered.

## 2026-11-01: Real test conversion to Ariadne Testing

Migrated `test_real` from the legacy runner to Ariadne Testing. The
`TestReal` fixture now derives from `ContextFixture`; its existing ten
functional methods and the earlier excessive-precision preflight are
registered as eleven separate cases in the previous execution order. Global
FloatDP approximation output precision is still configured at the beginning
of the expression-evaluation phase, immediately after preflight. The original
certified-pi TSV loader stays a private fixture helper.

All active expected numeric values, exception checks, executable smoke
operations and oracle inputs are retained. The conversion covers 39 assertions,
61 equality expectations, 3 exception expectations, 78 execute/no-throw
operations, 7 within-tolerance checks, 4 unary-predicate checks,
11 ordering comparisons, 12 explicit object constructions, and 36
active print/inspection events. The `std::boolalpha` stream manipulator
remains a process stream state change rather than being passed as an
unresolved overloaded function to the typed inspection API.
Two old statement-shaped equality macros were missing semicolons;
these were added when converting to the new macros. One previously
commented-out PRINT expression remains commented out.

CTest runs `test_real --concise`, whereas direct invocation defaults to
full reporter output. No independent oracle was rewritten and no failure
expectation was relaxed. Compile and runtime CI on every platform are required
before accepting the conversion; historical decision semantics for Real
predicates must be examined if generic evidence makes any case inconclusive.

## 2026-11-01: Real incompleteness and deprecated FloatDP construction

The initial `test_real` migration, commit `ba2ee99`, built and executed
but CTest reported two INCONCLUSIVE cases (arithmetic, transcendental):
`Real` equality is a non-decidable `NegatedSierpinskian` predicate.
Rendering equal values such as `x == 2.5_dy` or
`sqrt(Real(4)) == 2.0_dy` does not make generic evidence conclusive.
Neither blanket `decide()` nor treating unknown as pass is justified.

Twenty arithmetic and three transcendental equality checks with a `Real`
operand were replaced with explicit independently certified value oracles.
Each obtains `actual.compute(Accuracy(64_bits)).get()` as a rigorously
validated `DyadicBounds`, verifies that the original dyadic-exact expected
value belongs to the interval, and proves that the interval width is at
most `Dyadic(1,40u)` (2^-40). These are strong but *finite-resolution*
numerical checks, not claims of mathematically decidable equality.
The original expressions, expected constants, operation sequence and
other checks are unchanged. Every oracle carries its source expression
as a named diagnostic scope so failures remain attributable.

The two active `Real{deprecated_float}` constructor smoke tests are
replaced by `Real(Dyadic(binary_float))` and
`Real(ExactDouble(binary_float.get_d()))`, both preserving the exact
binary64 value while invoking supported constructors. The public
`Real(FloatDP)` API declaration and implementation remain for consumer
compatibility; this change removes its deprecated *use* within Numeric.
The separate `numeric_cast<Real>(FloatDP)` test is deliberately preserved:
`numeric_cast` is already specialized to construct
`Real(ExactDouble(a.get_d()))`, so it does not call the deprecated
constructor and still provides valuable public API coverage.

Eight remaining unconditional stdout debugging/progress diagnostics in
`test_real` are routed through the test reporter (notes/inspection) to
keep successful `--concise` execution quiet. The process's formatting
manipulator remains untouched.

This patch must be tested on Numeric's complete Ubuntu, macOS, Windows,
Sanitizers, and Coverage CI matrix; CI success is not inferred from the
static source audit. In particular, the bound-width expectation is a
numerical property that must be checked at runtime, not weakened merely
to silence an inconclusive outcome.

## 2026-11-01: Remove deprecated Real(FloatDP) and certify sqrt tolerance

Deleted the public `Real(FloatDP)` constructor declaration from
`include/numeric/real.hpp` and its implementation from `src/real.cpp`.
This intentionally removes a deprecated API; callers should use
`Real(Dyadic(float_value))` or `Real(ExactDouble(float_value.get_d()))`
when preserving the exact binary64 value. Existing
`numeric_cast<Real>(FloatDP)` remains supported by its already explicit
ExactDouble specialization, rather than via the removed constructor.

The new certified-value tests fixed the nondecidable-equality inconclusive
cases, but `test_real / transcendental` remained inconclusive on the
original `sqrt(x)` versus `sqrt(ax)` tolerance check. The general
`EXPECT_WITHIN` observes a Real-valued distance and may return unknown
when its interval overlaps the exact boundary. The sqrt case now evaluates
the FloatDPApproximation `sqrt(ax)` into its *exact* dyadic raw value and
computes certified 128-bit-accuracy dyadic bounds for `sqrt(x)`. It checks
both enclosure endpoints are within the **same original eps** of that raw
reference; no tolerance increase, ignored unknown, or numeric-oracle change.
The other six tolerance assertions are unchanged.

The whole Numeric CI matrix must pass before accepting this patch.

## 2026-11-01: Remove final indirect FloatDP literal construction

The public `Real(FloatDP)` declaration/definition had been removed, but the
`test_real` constructor case still called it indirectly through
`Real xe(1.5_exact)`: the `_exact` suffix produces `FloatDP`.
Use `Real xe{Dyadic(1.5_exact)}`, preserving the original exact dyadic
binary64 value and subsequent certified-bound test. This is the same
conversion as the former deprecated implementation, now explicit.
A compile-time contract in `tests/check_numeric.cpp` requires that
`Real` not be constructible from `FloatDP`, and checks that the
recommended `Dyadic` and `ExactDouble` alternatives are constructible.
An audit of numeric test files found no other `Real` constructions
from the `_exact` literal. The complete CI matrix is still required.

## 2026-11-01: Certify all seven Real transcendental tolerance tests

The previous migration fixed the `sqrt(x)` tolerance assertion by comparing
certified dyadic bounds against the exact dyadic value of the FloatDP
approximation. CI then revealed the same inconclusive phenomenon in
`exp(x) within 8*eps of exp(ax)`: `EXPECT_WITHIN` applied to a computable
`Real` difference cannot, in general, produce a conclusive comparison
at the prescribed tolerance. A displayed approximate error of `~0.~` is
not a proof of the requested numerical bound.

Rather than change the generic Testing framework, every transcendental
comparison (`sqrt`, `exp`, `log`, `sin`, `cos`, `tan`, `atan`) now
uses the fixture-local `expect_certified_within` helper. It evaluates the
same Real expression to a certified DyadicBounds enclosure at Accuracy(128_bits),
converts the same FloatDPApproximation expression's raw float to its exact
Dyadic value, and checks the lower and upper endpoints against that exact
reference minus/plus the **original tolerance**. All seven names and
multipliers are unchanged: eps, 8*eps, eps, 2*eps, 2*eps, 3*eps, eps.
The interval containment rule certifies that every possible value enclosed
by the algorithm is within tolerance; it does not infer mathematical
certainty from floating-point display or swallow inconclusive evidence.

This strengthens the interpretation relative to treating unknown as a pass,
and may legitimately fail if an actual numerical bound is too weak or the
original tolerance is not met. The seven case checks remain part of the
single original transcendental case. No changes to numeric implementation,
API, or generic Testing semantics are made by this patch. Acceptance
requires the complete CI matrix.

## 2026-11-01: Number suite conversion to Ariadne Testing

Migrated `test_number` on branch `testing` to ContextFixture-based
Ariadne Testing cases. Its original sequence is now 20 cases: the four
TestNumbers groups; the original class/get/comparisons groups for each of
ApproximateNumber, ValidatedNumber, EffectiveNumber and ExactNumber;
and the original operations group for each of the four directed types.
Empty specializations and previously unused methods remain inactive.

All 43 ASSERT checks, 31 EQUAL/EQUALS checks, 29 typed exception checks,
5 untyped exception checks, 4 tolerance checks and 49 binary predicate
checks retain their expectations and original expressions. The 9
`models` predicate calls are now fixture-local proposition checks, since
the testing framework's generic binary-predicate macro deliberately uses
a noncapturing lambda. Of the other binary predicates, 10 retain
explicit negation via EXPECT_NOT_BINARY_PREDICATE. The 3 runtime
same-type checks become compile-time contracts. All 29 smoke operations,
4 assignments, 22 construction declarations and 23 explicit print
inspections remain. The frozen 27-row Integer exponent oracle is intact.
Notifications and progress diagnostics are reporter notes/inspections.

The old unconditional `INCOMPLETE` stderr marker is not emitted during
a successful concise run; the prior incompleteness of the suite remains
visible in this journal, including the previously uninvoked methods.
CTest now links the target to `ariadne-testing` and passes the already
registered `--concise`; direct invocation uses full reporting.
The source and assertion counts have been audited. No local binary
was built because the build environment cannot access GitHub;
the Numeric CI matrix remains the acceptance gate, and any failing
numeric predicate must be classified before changing an oracle.

## 2026-11-01: Sync dependencies and certify Number approximation tolerances

Bumped `ss/fou` to `8c907955aaebbf766ce4edef96fd9d6548c49afb` and
`ss/log` to `e31295ba0d10a85b8bd07f00b0f438e0827f61d7`.
Both now pin the same Configuration (`d57304ec`), Testing (`cf1b9ad`),
and Utility (`38e5013`) commits. Added the direct `testing|configuration`
chain to the Configuration consistency check and added a Testing check
across the direct `testing` gitlink and the `foundation|testing` and
`logging|testing` chains. Other dependency gitlinks are unchanged.

The migrated `test_number` previously had three inconclusive `get`
cases: approximate output differences were small enough but the
`EXPECT_WITHIN` relation on `FloatDPApproximation` or
`FloatMPApproximation` could not yield definite evidence. All four
original approximate-get tolerance checks now convert the *raw binary
floating-point output* to an exact Dyadic and then Rational; they check
that the resulting exact rational absolute error is <= the original
tolerance using a conclusive ordering assertion. The original `2^-52`,
`2^-128`, and `3e-16` tolerances are unchanged. In particular
`Dyadic(3e-16)` preserves the literal's exact binary64 value rather
than substituting a rounded decimal or a larger margin. Diagnostic
scopes still identify the original checked expression, and both
precision paths and the 3/5 oracle remain active.

Three expected `ParadigmError` regression checks deliberately invoke
unsupported conversions from `FloatDPApproximation` to `Rational`,
`FloatDPBounds` and `FloatDPBall`. The three
`NumberGetterMixin::_get_as` branches now throw the *same exception*
without writing unsolicited `Warning: Cannot convert` lines directly
to stderr. All successful conversion code paths remain unchanged;
the failing conversion paths remain observable to callers via the
exception and to the test framework via `EXPECT_THROWS`.

This is source-level validation only. The full Numeric CI matrix
(Ubuntu, macOS, Windows, Sanitizers, Coverage) must pass on this commit
before treating the change as runtime-verified. Do not weaken numerical
expectations or silence unclassified failures.

## 2026-11-01: Restrict dependency checks to shared Configuration and Utility

Follow-up to the dependency synchronization above: `Testing` itself is
allowed to be pinned to different commits by Numeric, Foundation, and
Logging. Removed the cross-chain `require_same_dependency_commit(Testing)`
check; equal Testing gitlinks were a coincidental property of the
currently selected dependency revisions, not an architectural requirement.

The `Configuration` consistency check **still** includes
`testing|configuration` as requested, so the direct Testing submodule's
Configuration gitlink must agree with the other Configuration chains.
The pre-existing `Utility` consistency check remains unchanged.
No gitlinks, numerical checks, or test bodies were altered.
