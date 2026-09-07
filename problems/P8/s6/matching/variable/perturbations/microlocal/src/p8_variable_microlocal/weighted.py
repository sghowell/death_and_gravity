"""Exact quadratic completion underlying the semiclassical auxiliary weights.

Write k=sqrt(K), q=(pi,xi_g,chi_B). Curvature remainders and the xi_f
remainder have weight k^-1, while the two lapse remainders have weight
one. The resulting leading auxiliary block has no frequency denominator.
The full physical evolution proof must also exclude a nonuniform exact
Cauchy-map denominator; the point fixtures alone do not do so.
"""

from functools import cache

import sympy as sp

from . import leading
from . import stueckelberg as st

RG, RF, LG, LF, T = sp.symbols("r_g r_f lambda_g lambda_f t_f", real=True)
VRG, VRF = sp.symbols("r_g_prime r_f_prime", real=True)
REMAINDERS = (RG, RF, LG, LF, T)


@cache
def derive():
    d = leading.derive()
    u, a, b, c = (d[key] for key in ("u", "a", "b", "c"))
    shifts = dict(zip(st.LEADING_AUXILIARIES, REMAINDERS, strict=True))
    substitution = {q: d["stationary_fields"][q]+shifts[q] for q in st.LEADING_AUXILIARIES}
    substitution.update({
        st.VPSIG: -d["U"]*st.VPI-sp.diff(d["U"], u)*st.PI+VRG,
        st.VPSIF: d["V"]*st.VPI+sp.diff(d["V"], u)*st.PI+VRF,
    })
    completed_L1 = d["boundary_subtracted_L1"].subs(substitution, simultaneous=True)
    completed_L0 = d["boundary_subtracted_L0"].subs(substitution, simultaneous=True)
    auxiliary_hessian = sp.hessian(d["boundary_subtracted_L1"], st.LEADING_AUXILIARIES).applyfunc(sp.factor)
    expected = sp.Matrix([
        [2*a, 0, 2*a, 0, 0], [0, 2*c*b, 0, 2*c*b, 0],
        [2*a, 0, 0, 0, 0], [0, 2*c*b, 0, 0, 0],
        [0, 0, 0, 0, a*c**2*d["P"]/(c+d["y"])],
    ])
    # Rows are (E_r_g/k²,E_r_f/k²,E_lambda_g/k,E_lambda_f/k,E_t_f/k).
    # Columns are (lambda_g,lambda_f,k*r_g,k*r_f,k*t_f).
    leading_auxiliary = sp.diag(2*a, 2*c*b, 2*a, 2*c*b, d["xf_hessian"])
    return {**d, "remainders": REMAINDERS,
            "remainder_substitution": substitution,
            "completed_L1": completed_L1, "completed_L0": completed_L0,
            "auxiliary_hessian": auxiliary_hessian,
            "expected_auxiliary_hessian": expected,
            "leading_weighted_auxiliary": leading_auxiliary,
            "leading_weighted_auxiliary_inverse": leading_auxiliary.inv(),
            "weights": {"r_g": -1, "r_f": -1, "lambda_g": 0, "lambda_f": 0, "t_f": -1},
            "physical_equation_weight": -2}


@cache
def checks():
    d = derive()
    result = {}
    for i, entry in enumerate(d["auxiliary_hessian"]-d["expected_auxiliary_hessian"]):
        result[f"leading_auxiliary_literal_hessian_{i}"] = sp.factor(entry)
    # Check the completed coefficients instead of expanding one enormous
    # expression. All expressions are quadratic perturbation polynomials.
    r = sp.Matrix(REMAINDERS)
    difference = d["completed_L1"]-d["stationary_L1"]-(r.T*d["auxiliary_hessian"]*r)[0]/2
    variables = (*st.PHYSICAL, *st.PHYSICAL_VELOCITIES, *REMAINDERS)
    for i, x in enumerate(variables):
        for j, y in enumerate(variables[i:], i):
            result[f"exact_stationary_completion_{i}_{j}"] = sp.factor(sp.diff(difference, x, y))
    result["no_remaining_xif_velocity"] = sp.factor(sp.diff(d["completed_L0"], st.VXF))
    recovered = sp.hessian(d["completed_L0"], st.PHYSICAL_VELOCITIES).subs(
        dict.fromkeys((*REMAINDERS, VRG, VRF), 0))
    for i, value in enumerate(recovered-d["kinetic"]):
        result[f"principal_kinetic_not_changed_by_completion_{i}"] = sp.factor(value)
    for i, value in enumerate(d["leading_weighted_auxiliary"]*d["leading_weighted_auxiliary_inverse"]-sp.eye(5)):
        result[f"weighted_inverse_{i}"] = sp.factor(value)
    return result
