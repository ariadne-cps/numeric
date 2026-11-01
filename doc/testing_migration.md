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
