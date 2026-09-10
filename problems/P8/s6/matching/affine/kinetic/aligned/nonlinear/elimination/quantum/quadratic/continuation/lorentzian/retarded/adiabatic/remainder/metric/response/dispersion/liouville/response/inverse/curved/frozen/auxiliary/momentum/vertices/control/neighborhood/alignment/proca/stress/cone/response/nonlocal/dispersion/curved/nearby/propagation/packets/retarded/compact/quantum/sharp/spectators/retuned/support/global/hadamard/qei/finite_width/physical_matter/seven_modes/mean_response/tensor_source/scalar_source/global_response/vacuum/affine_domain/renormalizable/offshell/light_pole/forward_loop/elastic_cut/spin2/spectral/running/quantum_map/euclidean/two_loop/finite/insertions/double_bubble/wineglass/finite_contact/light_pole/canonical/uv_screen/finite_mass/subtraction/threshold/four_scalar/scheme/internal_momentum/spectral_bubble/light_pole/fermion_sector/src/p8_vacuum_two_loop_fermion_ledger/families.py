"""Noncommuting background-degree expansion and a complete nine-row ledger."""

from functools import cache

import sympy as sp


def catalog():
    result = []
    for k in range(3):
        for j in range(k + 1):
            n = k - j
            result.append(
                {
                    "id": f"scalar_Phi{2 * k}_W{j}_F{2 * n}",
                    "external_background_degree": 2 * k,
                    "contracted_boson": "Phi",
                    "tree_W_insertions": j,
                    "fermion_Hessian_background_degree": 2 * n,
                    "scalar_Y_power": n + 1,
                    "gauge_a_power": 0,
                    "trace_coefficient": sp.Rational((-1) ** j, 2),
                    "estimate_status": (
                        "Complete paired covariance-insertion quadratic sector bounded by S6.136"
                        if (k, j) == (1, 1)
                        else "Local-quartic W W subset with inner/outer reference pairing bounded by S6.135; full nonlocal W sector remains"
                        if (k, j) == (2, 2)
                        else "Not evaluated or bounded by this ledger"
                    ),
                }
            )
    for k in range(3):
        result.append(
            {
                "id": f"gauge_Phi{2 * k}",
                "external_background_degree": 2 * k,
                "contracted_boson": "gauge",
                "tree_W_insertions": 0,
                "fermion_Hessian_background_degree": 2 * k,
                "scalar_Y_power": k,
                "gauge_a_power": 1,
                "trace_coefficient": sp.Rational(1, 2),
                "estimate_status": "Not evaluated or bounded by this ledger",
            }
        )
    return result


@cache
def data():
    z = sp.symbols("background_degree_marker")
    D, W, F0, F2, F4 = sp.symbols("D W F0 F2 F4", commutative=False)
    inverse = D - z * z * D * W * D + z**4 * D * W * D * W * D
    Hessian = F0 + z * z * F2 + z**4 * F4
    product = sp.expand(inverse * Hessian / 2)
    expected = {
        0: D * F0 / 2,
        2: (D * F2 - D * W * D * F0) / 2,
        4: (D * F4 - D * W * D * F2 + D * W * D * W * D * F0) / 2,
    }
    # An independent finite matrix inverse recursion, with noncommuting factors.
    d = sp.Matrix([[2, sp.Rational(1, 3)], [sp.Rational(1, 3), 1]])
    w = sp.Matrix([[1, 2], [2, 3]])
    expansion = d - z * z * d * w * d + z**4 * d * w * d * w * d
    residual = (d.inv() + z * z * w) * expansion - sp.eye(2)
    h = sp.symbols("covariance_insertion_marker")
    f = sp.Matrix([[3, 1], [1, 2]])
    P = -d * f * d
    quartic_loop = -sp.trace((d + h * P) * w * (d + h * P) * w) / 4
    quadratic_loop = sp.trace((d + h * P) * w) / 2
    checks = {
        f"noncommuting_scalar_background_degree_{degree}": sp.expand(
            product.coeff(z, degree) - expr
        )
        for degree, expr in expected.items()
    }
    for i in range(2):
        for j in range(2):
            for degree in (0, 2, 4):
                checks[f"inverse_left_identity_{i}_{j}_{degree}"] = sp.expand(
                    residual[i, j]
                ).coeff(z, degree)
    checks.update(
        {
            "quartic_covariance_sector_matches_master_trace": sp.factor(
                sp.diff(quartic_loop, h).subs(h, 0)
                - sp.trace(d * w * d * w * d * f) / 2
            ),
            "quadratic_covariance_sector_matches_master_trace": sp.factor(
                sp.diff(quadratic_loop, h).subs(h, 0) + sp.trace(d * w * d * f) / 2
            ),
            "six_scalar_rows": sum(r["contracted_boson"] == "Phi" for r in catalog())
            - 6,
            "three_gauge_rows": sum(r["contracted_boson"] == "gauge" for r in catalog())
            - 3,
            "nine_total_rows": len(catalog()) - 9,
            "four_four_point_rows": sum(
                r["external_background_degree"] == 4 for r in catalog()
            )
            - 4,
            "three_two_point_rows": sum(
                r["external_background_degree"] == 2 for r in catalog()
            )
            - 3,
            "two_vacuum_constant_rows": sum(
                r["external_background_degree"] == 0 for r in catalog()
            )
            - 2,
        }
    )
    return {
        "scalar_covariance_formal_background_expansion": inverse,
        "fermion_scalar_Hessian_background_expansion": Hessian,
        "scalar_trace_background_degree_coefficients": {
            str(k): v for k, v in expected.items()
        },
        "complete_degree_zero_two_four_catalog": catalog(),
        "quartic_master_identity": "Gamma_2,F[Phi^4]=1/2 Tr[D F4-D W2 D F2+D W2 D W2 D F0]+1/2 Tr[D_A F_AA,4], before counterterm insertion and canonical conversion.",
        "W_dictionary": "W2 is the complete quadratic-in-Phi Hessian of the frozen exactly H-eliminated quartic action, including all nonlocal heavy contractions, not just L Phi^2/2.",
        "F_dictionary": "F0,F2,F4 are the Phi-background degree 0,2,4 coefficients of the scalar Hessian of the entire one-loop fermion determinant; no extra factorials are absorbed after this definition.",
        "scope": "Complete primitive two-loop fermion-sector ownership through four external Phi derivatives, not completed integrals, counterterm coefficients, scheme conversion or an error budget.",
        "checks": checks,
    }
