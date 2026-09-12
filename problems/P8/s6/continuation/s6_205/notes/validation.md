# Reproducibility and verification boundary

The read-only report has eighteen explicitly enumerated sources and twenty top-level fields. The scientific sources include the complete degree partition, radial moments, original-cutoff tail and exact restoration identity, with seven written notes and two test files.

The audit contains 29 named exact checks, 37 scalar entries, 39 gates, nine controls and 127 rejected scope inputs. Independent scientific tests check actual mixed source-time/spatial coefficients at three separated momentum scales, every finite cell and time jet, all physical pairs, and two spatial Cauchy resolutions. They also check all four radial moments, the original removed-union geometry, half-band power tails, zero-transfer preservation and the borderline positive-majorant limitation.

Exact integer-order validation occurs before memoization. Tests prewarm valid cache entries and then reject booleans, floating-point lookalikes, symbolic/nonfinite values and out-of-range inputs. APIs for finite-cell bounds reject all fifteen candidate UV cells instead of silently assigning convergent moments to them.

The integrated scientific file contains 390 tests. Test counts alone are not a proof of continuous estimates: the joint-domain Cauchy argument, dominated convergence and regulator removal are written mathematics, not FORMALIZED proofs.

Before freezing, formatting/lint and a fresh private scientific process must pass. Frozen bytes are copied exactly and the repository scientific replay is rerun. The report is independently rebuilt through the read-only native process, preserving original SymPy. Fresh ordinary and CLI checks also use original SymPy. Only the complete P8 regression uses the separately audited exact GCD adapter and its original-function tuple comparisons. Final run outcomes and hashes belong in the publication audit; a running process is not a passed check.
