# Migration of Numeric tests to Ariadne Testing

2026-11-01: `test_operators` and `test_complex` migrated on branch `testing`.

- Both test binaries now parse `--concise` using `Ariadne::Testing::run_suite(suite, argc, argv)`.
- `test_operators` retains its five functional groups as independent cases; its metadata helper records two assertions per invocation. Execution-only operations become explicit no-throw checks, and the invalid operator-construction case still requires an exception.
- `test_complex` retains the six original fixture groups plus Fourier computation, with fresh fixtures per group and the original FloatDPBounds precision. The code continues to avoid enabling commented-out tests.
- Existing runtime assertions are migrated to evidence-aware checks; any inconclusive evidence is now explicitly classified rather than coerced into a boolean. This semantic difference is intentional and must be checked by CI.
- The new suite macros record failed checks without disabling subsequent checks. Errors thrown from one group no longer suppress unrelated groups.

CMake links the two migrated executables with `ariadne-testing`; all other targets retain their prior configuration and CTest arguments. Legacy test binaries may still ignore `--concise` until migrated.

**Validation status:** The source transformation and CMake registration have been checked, but a C++ build or CTest run was not possible in this environment because GitHub cloning/network access is unavailable from the build container. Treat compilation and the Numeric CI matrix as a mandatory acceptance gate; do not interpret the migration as verified until those checks pass.
