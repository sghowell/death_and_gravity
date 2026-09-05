"""Necessary disformal-to-quartic-Horndeski map for the pinned open D tube.

Derive from ADM kinetic mixing and the covariant Horndeski tensor identity.
No assertion that tensor normalization alone supplies an Einstein-scalar seed.
"""

from functools import cache

import sympy as sp
from p8 import adm
from p8_s5 import nonlinear_d as model


@cache
def identities():
    X = sp.Symbol("X", positive=True)
    C, W, M = sp.symbols("C W M", positive=True)
    CX, WX, BX, A3 = sp.symbols("C_X W_X B_X A3", real=True)
    K, KK, V, temporal_shift = sp.symbols("K KK V temporal_shift", real=True)
    shift = sp.sqrt(X)*CX/C*V+temporal_shift
    # In a Horndeski seed K_E has no lapse-velocity term in its ADM action.
    # Pullback: K_E^i_j=(K^i_j+shift*delta^i_j)/sqrt(W).
    # The mapped K^2-Kij^2 coefficient is fixed to -M^2/2 by target G_T.
    trace, norm = K+3*shift, KK+2*K*shift+3*shift**2
    kinetic = -M**2*(trace**2-norm)/2
    KV = sp.expand(kinetic).coeff(K, 1).coeff(V, 1)
    solved_CX = sp.solve(KV-X**sp.Rational(3, 2)*A3, CX)[0]
    J = W-X*WX
    # Horndeski identity in the seed: G_E=F_E-2*X_E*dF_E/dX_E.
    # Inverse tensor map is used before imposing that identity.
    seed_F = M**2/sp.sqrt(C*W)
    seed_G = M**2*sp.sqrt(W)/C**sp.Rational(3, 2)
    seed_F_X = sp.diff(seed_F, C)*CX+sp.diff(seed_F, W)*WX
    seed_X_X = J/W**2
    horndeski_residual = seed_G-seed_F+2*(X/W)*seed_F_X/seed_X_X
    reduced_condition = sp.simplify(horndeski_residual*C*J/(seed_F*W))
    # W=C+B*X, W_X=C_X+B+X*B_X. Eliminate B=(W-C)/X.
    in_BX = sp.factor(reduced_condition.subs(WX, CX+(W-C)/X+X*BX))
    solved_BX = sp.solve(in_BX, BX)[0]
    required_J = sp.factor(J.subs(WX, CX+(W-C)/X+X*solved_BX).subs(CX, solved_CX))
    coefficients = adm.coefficients()
    h_locus = {coefficients["A1"]: 2*coefficients["F2X"], coefficients["A3"]: 0}
    data = model.functions()
    d, u = model.d, model.u
    a3 = data["a3"]
    Cstar = sp.exp((1-X**2)/d**3)
    y = sp.Symbol("y", positive=True)
    Bstar = 4*sp.Integral(sp.exp((1-y**2)/d**3), (y, 1, X))/d**3
    Wstar = Cstar+X*Bstar
    Jstar = sp.simplify(Wstar-X*sp.diff(Wstar, X))
    factor = 1-2*X**2/d**3
    clock_factor = sp.factor(factor.subs(X, 1))
    root = sp.sqrt(sp.real_root(2, 3)-1)
    target = {coefficients["X"]: X, coefficients["F2"]: -M**2/2,
              coefficients["F2X"]: 0, coefficients["A1"]: 0, coefficients["A3"]: M**2*a3}
    target_mix = coefficients["C"].subs(target)
    residuals = {"KV_pullback": sp.simplify(KV+2*M**2*sp.sqrt(X)*CX/C),
                 "target_ADM_KV": sp.simplify(target_mix-M**2*X**sp.Rational(3, 2)*a3),
                 "Horndeski_tensor_identity": sp.simplify((coefficients["GT"]-coefficients["Lambda"]).subs(h_locus)),
                 "Horndeski_lapse_mixing_zero": sp.simplify(coefficients["C"].subs(h_locus)),
                 "necessary_W_condition": sp.simplify(reduced_condition-(W-C-X*WX-X*CX)),
                 "necessary_B_X": sp.simplify(solved_BX+2*CX/X),
                 "necessary_J": sp.simplify(required_J-C*(1-X**2*A3/(2*M**2))),
                 "C_primitive": sp.simplify(sp.diff(Cstar, X)/Cstar+X*a3/2),
                 "B_primitive": sp.simplify(sp.diff(Bstar, X)+2*sp.diff(Cstar, X)/X),
                 "map_jacobian": sp.simplify(Jstar-Cstar*factor),
                 "prior_spatial_map": sp.simplify(Cstar-sp.exp(-2*data["omega"].subs(model.N, 1/sp.sqrt(X)))),
                 "prior_Lambda": sp.simplify(clock_factor-data["Lambda"]),
                 "positive_root": sp.simplify(clock_factor.subs(u, root)),
                 "negative_root": sp.simplify(clock_factor.subs(u, -root)),
                 "positive_root_polynomial": sp.simplify(((1+u**2)**3-2).subs(u, root)),
                 "root_polynomial_monotonicity": sp.diff((1+u**2)**3-2, u)-6*u*(1+u**2)**2,
                 "normalized_C": Cstar.subs(X, 1)-1,
                 "normalized_W": sp.simplify(Wstar.subs(X, 1)-1)}
    return {"required_C_X_over_C": sp.factor(solved_CX/C),
            "required_B_X": solved_BX, "required_J": required_J,
            "normalized_C_example": Cstar, "normalized_B_example": Bstar,
            "clock_J_over_C": clock_factor, "positive_crossing_time_over_tau": root,
            "negative_crossing_time_over_tau": -root,
            "root_bracket_signs": [sp.expand(((1+u**2)**3-2).subs(u, v)) for v in (0, 1)],
            "residuals": residuals}


def derive():
    result = identities()
    failed = {key: sp.simplify(value) for key, value in result["residuals"].items()
              if sp.simplify(value) != 0}
    if failed:
        raise ValueError(f"Horndeski-map identity failed: {failed}")
    return {**{key: str(value) for key, value in result.items() if key != "residuals"},
            "exact_residuals": {key: "0" for key in result["residuals"]},
            "result": "No smooth everywhere-invertible one-clock C(phi,X),B(phi,X) disformal map of the exact D tube to a quartic Horndeski action",
            "metric_vs_map_control": "C=W=1 at X=1, even at both zero-J times: regular metric determinants do not imply invertible field transformation",
            "Einstein_seed_subcase": "Exact open-tube G_T=F_T=M^2 forces C=W=1 for an EH seed; C_X=0 contradicts the nonzero D lapse mixing",
            "scope": "Exact tube map only; failed coordinates are not a singularity of the original witness, a general parent exclusion, or a bound for arbitrary tube-deforming corrections"}
