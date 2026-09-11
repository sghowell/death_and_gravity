"""Off-shell first metric/field variation of the frozen clock-transparent map."""

from functools import cache

import sympy as s
from p8_vacuum_clock_transparent_map import gate


@cache
def data():
    q = gate.data()["gate"]
    X = gate.X
    e, xi, r0, r1, Ephi = s.symbols("epsilon delta_X R0 delta_R E_Phi", real=True)
    mapped = q.subs(X, 1 + e * xi) * (r0 + e * r1)
    bad = (1 - (1 + e * xi)) * (r0 + e * r1)
    first = s.diff(mapped, e).subs(e, 0)
    checks = {f"same_clock_gate_jet_{j}": s.diff(q, X, j).subs(X, 1) for j in range(8)}
    checks.update(
        {
            "first_general_clock_variation": s.factor(first),
            "second_general_clock_variation": s.factor(s.diff(mapped, e, 2).subs(e, 0)),
            "off_shell_metric_chain_correction": s.factor(Ephi * first),
            "nonzero_zeroth_only_map_control": s.factor(
                s.diff(bad, e).subs(e, 0) + xi * r0
            ),
            "unchanged_first_nonidentity_clock_variation": s.factor(
                s.diff(q, X, 8).subs(X, 1) - 6435 * s.factorial(8)
            ),
        }
    )
    return {
        "same_map": "Fhat(Psi,g)=Psi+q8(X)R_old[Psi,g], X=(partial Psi)^2/kappa, with the exact S6.164 polynomial and its fixed parent coefficients.",
        "variational_statement": "On the entire unit clock X=1, Fhat=Psi and its first variation is deltaPsi: q8(1)=q8'(1)=0 kills delta(q8 R_old) for arbitrary compact metric and field variations, even if R_old and the parent scalar Euler derivative are nonzero. The equality is an identity along the clock, so spacetime derivatives of its zero first variation also vanish.",
        "classical_stress_consequence": "By the variational chain rule, the action pullback has the same classical metric Euler derivative on that clock. The extra off-shell term is the parent scalar Euler derivative paired with delta_g Fhat=0. No parent field equation or inverse map is assumed. The map therefore cannot turn the nonzero minimal Einstein null residual into a solution.",
        "negative_control": "Equality of Fhat and Psi in value alone is insufficient: replacing q8 by1-X leaves a nonzero first variation -R_old deltaX and hence an off-shell metric chain term. The actual eightfold zero is essential.",
        "scope": "This is the classical action/first-variation statement. It is not a bound on the full scalar quantum Jacobian, a global inverse, or a new interacting stress estimate. Other leading curvature and derivative operators can change the null equation; they must be matched explicitly.",
        "checks": checks,
        "gates": {
            "value_and_first_jet_both_zero": bool(
                q.subs(X, 1) == 0 and s.diff(q, X).subs(X, 1) == 0
            ),
            "zeroth_only_control_changes_metric_variation": bool(
                s.diff(bad, e).subs(e, 0) != 0
            ),
            "off_shell_no_background_equation_assumed": True,
        },
    }
