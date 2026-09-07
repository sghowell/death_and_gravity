# Opt-in exact P8 regression runner

The 2026-09-06 continuation encountered very slow univariate polynomial
cancellation in the unchanged SymPy 1.14.0 Gaussian-integer PRS path.
An isolated legacy symbolic-time quadratic bridge still passed under the
ordinary pytest command in 82.77 seconds, but a combined run spent more
than 33 minutes in that routine. A second ordinary run hit the same
routine in another legacy test. These runs were stopped, not counted as
passed. No calculation source, dependency or certificate was altered.

The [optional runner](../scripts/p8_exact_regression.py) changes only the
in-process GCD dispatch during pytest, with an explicit version guard:

```sh
.venv/bin/python scripts/p8_exact_regression.py --self-check
.venv/bin/python -m pytest scripts/tests/test_p8_exact_regression.py -q
PYTHONHASHSEED=0 .venv/bin/python scripts/p8_exact_regression.py problems/P8 -q
```

This is an **adapted exact-arithmetic regression**, not a successful run
of the unmodified full-suite command. Targeted child tests and standalone
certificate CLIs also run with the ordinary, unadapted interpreter.

## Why the descent is exact

Only nonzero univariate polynomials over Q(i) or Z[i] are admitted, and
each operand must separately be real or purely imaginary. Remove its
unit phase 1 or i and calculate the real-domain GCD. For real rational
polynomials the monic GCD over Q(i) is fixed by complex conjugation, so
belongs to Q[x] and equals the monic rational GCD. Extension of a field
does not introduce a new common divisor of relatively prime univariate
polynomials, by the Bezout identity. Over Z[i], the same primitive-part
argument and the integer contents give the normalized real integer GCD;
the removed phases are units. Restore the original Gaussian coefficient
domain and put each phase back in its corresponding cofactor.

Every accelerated call checks both exact polynomial identities
`gcd*cofactor_left=left` and `gcd*cofactor_right=right`, raising an error
on disagreement even under Python optimization. No floating-point values
or approximate zeros are introduced. Mixed complex coefficients, multiple
generators, zero operands and other domains use the original algorithm.
The original method is restored on exit, including exceptions. Unexpected
versions or a preexisting replacement are rejected.

Before running pytest, 128 fixtures compare the entire normalized triple
against the original Gaussian implementation, including both coefficient
domains, all four unit phases and nonintegral rational contents. The
[separate tests](../scripts/tests/test_p8_exact_regression.py) audit the
fallback, restoration, version and deliberately corrupted-cofactor gates.
These finite comparisons supplement the algebraic justification; they do
not replace it or constitute a theorem-prover formalization.

The installed SymPy implementation was inspected directly: real Q/Z
domains use their dedicated GCD methods, while the Gaussian domains take
the generic dense PRS route. The distinction is also discussed in
[SymPy issue 23131](https://github.com/sympy/sympy/issues/23131). The runner
is deliberately local and opt-in, not a global dependency modification.
