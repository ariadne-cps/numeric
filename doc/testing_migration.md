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
