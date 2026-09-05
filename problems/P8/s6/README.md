# P8 S6 — Adopted conditional UV framework and first result

The user delegated the UV-assumptions choice on 2026-09-04. We adopt a
vacuum-based, Lorentz-invariant, local/causal/unitary UV target with explicit
EFT matching to the bounce, and retain gravitational dispersion corrections.
This resolves the policy-choice blocker; it does not establish UV viability.

Read the [new contract](FORMULATION.md), [primary-source review](notes/literature.md)
and [S6.1 proof and limitations](notes/vacuum-extension.md).

## First result: vacuum information is independent of the smooth clock tube

An explicit smooth D-only family equals the old action on an open set
containing the **whole** clock tube, but near X=0 becomes Einstein gravity
plus a healthy massive scalar with arbitrary coefficient lambda of
`((partial psi)^2)^2`. The independently expanded vacuum amplitude has
`b2=4*lambda`. Thus opposite leading vacuum positivity diagnostics can be
attached without changing any existing on-tube result. This is an algebraic
freedom statement, **not a UV-compatible bouncing model**.

Requiring connected global real analyticity and exact open-tube equality
instead forces the original F, which has no finite constant-clock Minkowski
vacuum. That stronger conditional obstruction excludes neither approximate
matching, other witnesses, nor all UV completions.

The [report](certificates/vacuum-extension.json) pins the previous checkpoint,
all new source/test/proof inputs, the vacuum algebra and independent
heavy-mediator matching control. The full UV matching gate remains open;
neither pole subtraction nor the earlier M*tau bound supplies it. Bound
bookkeeping refuses a verdict when gravitational/truncation bounds are unknown.

## Reproduce

```sh
PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/control/src:problems/P8/s5/scattering/src:problems/P8/s6/src .venv/bin/python -m p8_uv.verify --check
.venv/bin/python -m pytest problems/P8/s6/tests -q
.venv/bin/ruff check problems/P8/s6
git diff --check
```

Without `--check` the verifier prints a candidate report; it never overwrites
one. Smoothness/identity-theorem arguments are written proofs, not formalized
by the finite machine replay. Earlier certificates remain unchanged.

Next: S6.2 tests a common vacuum-to-tube EFT matching construction for the
positive-lambda member. The splice transition's health, heavy gap and
quantum origin must be established or the ansatz rejected. M1 interactions,
all orders, nonlinear stability and P8(a) remain separate unfinished work.
