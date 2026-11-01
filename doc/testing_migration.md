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
