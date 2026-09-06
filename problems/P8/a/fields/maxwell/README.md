# A.16 — photon-field conformal QSEI

A new free-Maxwell field branch, with an exact all-Hadamard absolute QSEI
on smooth spatially flat conformal strips. The conformal reference anomaly
and the independent symbolic finite-curvature parameter beta_M are retained.
The coefficient is hbar/(8pi^2) times the squared proper operator

    L_H f=f''-2Hf'+(3H^2/4-3Hdot/2)f.

[FORMULATION.md](FORMULATION.md) states the field, local state domain,
prescription family, H2 boundary traces and conditional focusing match.
[The proof](notes/proof.md) gives the local positive-type argument, full
source/sign dictionary and exact radiation IBP check.
[Sources](notes/sources.md) distinguish new Maxwell physics from the
immutable scalar/general-clock inputs.

The rational curvature envelopes do not assert that any supplied metric
satisfies them. All classical-source and cosmological signs stay visible.
A negative control proves that the naive cap Hmax*tau<=1 cannot satisfy the
frozen A.1 zeta=0 contraction test even when quantum constants are tiny.
No exact new SEE solution, cosmological focusing or P8 closure is claimed.

## Read-only replay

From the repository root, run:

```sh
.venv/bin/pytest -q problems/P8/a/fields/maxwell/tests
.venv/bin/ruff check problems/P8/a/fields/maxwell
.venv/bin/python -c 'from pathlib import Path; import sys; root=Path("problems/P8/a/fields/maxwell"); old=root.parents[1]/"applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei"; sys.path[:0]=[str(root/"src"), str(old.parents[2]/"qsei"/"src"), *(str(p/"src") for p in (old,*old.parents))]; from p8a_maxwell.verify import main; main()' --check
```

The verifier only reads: it checks pinned A.1/A.12 reports and replays their
lineage, exact symbolic identities, independent Fraction reconstruction and
all source hashes against [the certificate](certificates/maxwell-qsei.json).
Without `--check` it prints a candidate report to stdout; it does not write
or repair a certificate. Changed claims, constants, sources or parent pins
fail validation. Independent coordinate tests include the parent-owned
`test_maxwell_covariant_audit.py` unchanged.

## Public calculation interfaces

- `stress.reference_jets` / `reference_stress`: every physical component,
  trace and EED, with the finite beta_M term explicit.
- `functional.proper_operator` / `ibp_coefficients` / `radiation_control`:
  exact-clock and endpoint-sensitive sampler identities.
- `bounds.envelope` / `qsei_coefficients`: exact rational cap inputs and a
  mandatory `beta_m` argument; generic symbolic functions are separate.
- `focusing.geometric_constants`: conditional Ricci constants with explicit
  `cosmological_constant` and signed `other_eed_lower`.

Numerical interfaces reject booleans, binary floats, nonfinite values and
unresolved symbolic inputs. Neither the code nor the report certifies an
uninspected metric, arbitrary endpoint data or an unspecified source.
