# S6.41: mass retuning preserves matching but widens the matter cone

The separately named source-centered trace-mass update exactly retains
the original auxiliary CD/M1 action at zero curl. Its full quotient is
uniformly invertible and its source preserves the original primary
lapse degeneracy on the clock tube. The complete nonlinear secondary
constraint system is not claimed solved.

For small positive curl, the complete rolling scalar constraints give
positive kinetic and principal gradient matrices. Nevertheless one
physical high-frequency scalar speed exceeds the matter light cone at
every finite noncentral time. This fails the original witness acceptance
condition: positive energy is not enough. No below-cutoff exclusion,
general affine no-go, UV completion or original P8 closure is inferred.

See the [literal action and scope](FORMULATION.md),
[matching and primary-constraint proof](notes/geometry.md), and
[physical constraints and cone obstruction](notes/dynamics.md).

## Read-only replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/affine/kinetic/retuned/tests", "-q", "-p", "no:faulthandler"]))'
```

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'import sys, runpy; from pathlib import Path; from sympy.core.random import seed; seed(0); sys.path[:0]=[str(p.resolve()) for p in sorted(Path("problems/P8").rglob("src"))]; sys.argv=["p8_affine_retuned.verify", "--check"]; runpy.run_module("p8_affine_retuned.verify", run_name="__main__")'
```

Both commands seed Python hashing and SymPy's arithmetic RNG. They do
not use the broad-regression GCD adapter. The diagnosed host
faulthandler timer plugin remains disabled. Without `--check`, the CLI
prints the report and does not write files.
