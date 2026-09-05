# CD/M1 coupled free propagation

S5.7.CD completes the coupled quadratic normalization and a conservative
local free-energy estimate for the fixed matter-backed bounce. The
[formulation](FORMULATION.md) freezes the scope and the
[derivation](notes/normalization.md) supplies the analytic proof.

The full finite-q scalar mixing and all time-dependent generators survive
the calculation. On overlapping charts, q>=10^20 is sufficient for
positive coupled oscillator potential. If q0>=2*10^20 at the center of
abs(t-t0)<=ell0/100, every exact free scalar/tensor solution has energy
ratio between 7/11 and 11/7 throughout the window. This very conservative
threshold is **not an interacting cutoff** or an optimized scale.

At the bounce, the independent exact q=8 control has positive scalar
kinetic matrix but indefinite oscillator potential. This diagnoses the
failure of that positive-energy test, not a dynamical instability.

## Replay

From the repository root, using the existing environment:

```bash
export PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src:problems/P8/s5/matter/physical/control/src
.venv/bin/python -m p8_m1_control.verify --check
.venv/bin/python -m pytest -q problems/P8/s5/matter/physical/control/tests
```

The verifier is read-only. Without `--check` it prints the deterministic
report to stdout. The report pins its source files and the prior S5.6.CD
certificate/source chain; the complete P8 suite also replays earlier gates.

## What follows

The next M1 gate must substitute the exact inverse **phase-space** map

    Q=a^-3/2*T^-1*Y, p=a^-3/2*T^T*(P-S*Y)

into S5.6's cubic/quartic Hamiltonians and bound the resulting kernels on
an explicit nonexceptional momentum domain. S generally depends on q;
using only the kinetic square root or substituting free velocities into
the old interaction kernels would omit terms. Free energy estimates do
not establish weak interaction strength, radiative protection, or UV
admissibility. All of those remain separate open gates.
