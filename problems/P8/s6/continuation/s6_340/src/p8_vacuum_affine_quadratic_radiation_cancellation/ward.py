"""Independent inverse Ward identity and exact real-TT boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree


@cache
def data():
    checks = {}
    eta = s.diag(1, -1, -1, -1)
    pv = s.Matrix(s.symbols("p0:4"))
    w = s.Symbol("omega", nonzero=True)
    kv = s.Matrix([w, 0, 0, w])
    rv = pv + kv
    p2 = (pv.T * eta * pv)[0]
    r2 = (rv.T * eta * rv)[0]
    d = s.factor(r2 - p2)
    Kp, Kr = s.symbols("Kp Kr")
    Gamma = (pv * rv.T + rv * pv.T) * (Kr - Kp) / d - eta * (Kr + Kp) / 2
    checks["arbitrary_inverse_full_offshell_Ward"] = (
        kv.T * eta * Gamma - (pv * Kr - rv * Kp).T
    ).applyfunc(s.factor)
    checks["free_inverse_literal_minimal_tensor"] = (
        Gamma.subs({Kp: p2 - 1, Kr: r2 - 1}) - tree.tensor(pv, rv, 1, eta)
    ).applyfunc(s.factor)

    # Null-transverse two-direction homogeneous tensors are physical-TT invisible.
    b, c = s.symbols("homogeneous_b homogeneous_c")
    a = (pv.T * eta * kv)[0]
    hom = b * (pv * kv.T + kv * pv.T - a * eta) + c * kv * kv.T
    checks["two_direction_homogeneous_Ward"] = (kv.T * eta * hom).applyfunc(s.expand)
    ep = s.diag(0, 1, -1, 0)
    cross = s.zeros(4)
    cross[1, 2] = cross[2, 1] = 1
    for label, eps in (("plus", ep), ("cross", cross)):
        checks["two_direction_homogeneous_" + label] = s.expand(
            sum(hom[i, j] * eps[i, j] for i in range(4) for j in range(4))
        )
    A, B, D = s.symbols("A B D")
    pk = s.Symbol("pk", nonzero=True)
    # Coefficients of the independent p and k in k_mu H^{mu nu}.
    solution = s.solve((A * pk, B * pk + D), (A, D), dict=True)
    assert solution == [{A: 0, D: -B * pk}]
    return {
        "checks": checks,
        "gates": {
            "literal_minimal_scalar_vertex_normalization": True,
            "full_arbitrary_endpoint_inverse_Ward_not_only_soft": True,
            "nonzero_transfer_only_assumption_for_homogeneous_coefficient": True,
            "possible_unprojected_kk_form_factor_not_declared_absent": True,
            "real_null_two_direction_scope_not_four_hard_Weyl_coordinate": True,
        },
        "whole_generic_inverse_representative": Gamma,
        "whole_null_transverse_homogeneous_tensor": hom,
        "whole_homogeneous_solution": {str(k): v for k, v in solution[0].items()},
        "whole_physical_boundary": "For r=p+k and null k, Gamma=(pr+rp)DD_K-eta(K(r^2)+K(p^2))/2 obeys k.Gamma=p K(r^2)-r K(p^2). The free inverse gives exactly the frozen minimal tensor. With only p,k, a homogeneous symmetric transverse tensor has no pp coefficient for p.k!=0, so its physical TT projection vanishes. Analyticity extends through the soft limit; no unprojected kk coefficient is fixed.",
        "whole_not_inferred": "This two-direction real-null statement does not remove four-hard-direction curvature operators, off-shell graviton response, massless-graviton singularities or independent higher-order physical matching.",
    }
