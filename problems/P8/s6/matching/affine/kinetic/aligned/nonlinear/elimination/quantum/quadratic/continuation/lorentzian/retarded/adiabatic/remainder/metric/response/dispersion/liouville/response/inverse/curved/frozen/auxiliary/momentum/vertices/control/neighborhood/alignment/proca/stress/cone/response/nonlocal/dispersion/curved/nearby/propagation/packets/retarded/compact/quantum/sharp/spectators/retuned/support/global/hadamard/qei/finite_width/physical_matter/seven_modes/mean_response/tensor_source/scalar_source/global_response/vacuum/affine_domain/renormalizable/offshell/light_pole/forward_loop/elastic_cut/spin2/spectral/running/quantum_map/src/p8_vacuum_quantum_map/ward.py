"""Explicit change-of-variable controls and one-loop external-pole bookkeeping."""

from functools import cache

import sympy as sp


@cache
def data():
    y, r = sp.symbols("Gaussian_field cubic_map_parameter", real=True)
    hbar, K = sp.symbols("loop_order positive_Gaussian_kernel", positive=True)
    F = y + r * y**3
    jac = sp.diff(F, y)
    weight = (
        1
        - K * r * y**4 / hbar
        + r**2 * (K**2 * y**8 / (2 * hbar**2) - K * y**6 / (2 * hbar))
    )

    def average(expr):
        result = sp.Integer(0)
        for (power,), coefficient in sp.Poly(sp.expand(expr), y).terms():
            if power % 2 == 0:
                result += (
                    coefficient * sp.factorial2(power - 1) * (hbar / K) ** (power // 2)
                )
        return sp.expand(result)

    checks = {}
    moment_rows = {}
    for n in (0, 2, 4, 6):
        corrected = average(sp.series(jac * weight * F**n, r, 0, 3).removeO())
        baseline = average(y**n)
        moment_rows[n] = corrected
        for order in (0, 1, 2):
            checks[f"Gaussian_source_moment_{n}_map_order_{order}"] = sp.expand(
                corrected.coeff(r, order) - (baseline if order == 0 else 0)
            )
    no_jac = average(weight).coeff(r, 1)
    no_source = average(jac * weight * y**2).coeff(r, 1)
    no_sextic = average(jac * (weight + r**2 * K * y**6 / (2 * hbar))).coeff(r, 2)
    checks.update(
        {
            "omitted_Jacobian_nonzero_control": sp.expand(no_jac + 3 * hbar / K),
            "omitted_source_nonzero_control": sp.expand(no_source + 6 * hbar**2 / K**2),
            "omitted_sextic_nonzero_control": sp.expand(
                no_sextic - 15 * hbar**2 / (2 * K**2)
            ),
        }
    )
    a0, a1, w = sp.symbols(
        "physical_tree_amplitude physical_one_loop_amplitude composite_one_loop_overlap"
    )
    overlap = 1 - hbar * w
    amplitude = a0 + hbar * a1
    residue = sp.series(overlap**2, hbar, 0, 2).removeO()
    pole_numerator = sp.series(overlap**4 * amplitude, hbar, 0, 2).removeO()
    amputated = sp.series(pole_numerator / residue**4, hbar, 0, 2).removeO()
    lsZ = sp.series(residue**2 * amputated, hbar, 0, 2).removeO()
    checks.update(
        {
            "one_particle_residue": sp.expand(residue - (1 - 2 * hbar * w)),
            "four_external_pole_numerator": sp.expand(
                pole_numerator - a0 - hbar * (a1 - 4 * w * a0)
            ),
            "amputated_vertex_keeps_composite_overlap": sp.expand(
                amputated - a0 - hbar * (a1 + 4 * w * a0)
            ),
            "LSZ_cancels_all_four_external_overlap_factors": sp.expand(lsZ - amplitude),
        }
    )
    L = sp.Symbol("diagnostic_quartic")
    S = K * y * y / 2 + L * y**4 / 24
    transformed = S.subs(y, F)
    one_loop_transformed = sp.diff(transformed, y, 2)
    # The additive log(K)/2 is irrelevant; use ratios with nonzero constant one.
    gamma1 = sp.series(
        sp.log(one_loop_transformed / K) / 2 - sp.log(jac), y, 0, 5
    ).removeO()
    gamma1_of_F = sp.series(sp.log((K + L * F**2 / 2) / K) / 2, y, 0, 5).removeO()
    difference4 = sp.expand(gamma1 - gamma1_of_F).coeff(y, 4)
    checks["off_shell_effective_action_not_scalar_control"] = sp.expand(
        difference4 + 24 * r * r + L * r / K
    )
    # Literal scalar chain rule, including the term that vanishes only at a classical extremum.
    sx = (K * y + L * y**3 / 6).subs(y, F)
    sxx = K + L * F**2 / 2
    checks["Hessian_chain_rule_keeps_equation_term"] = sp.expand(
        one_loop_transformed - jac**2 * sxx - sx * sp.diff(F, y, 2)
    )
    return {
        "finite_dimensional_formal_map": F,
        "full_Gaussian_weight_through_map_order_two": weight,
        "full_source_moments": moment_rows,
        "omitting_Jacobian_changes_partition_function": no_jac,
        "omitting_transformed_source_changes_two_point_function": no_source,
        "omitting_generated_sextic_changes_partition_function": no_sextic,
        "one_particle_overlap": overlap,
        "one_particle_residue": residue,
        "four_external_pole_numerator": pole_numerator,
        "ordinary_amputated_four_point_vertex": amputated,
        "proper_LSZ_amplitude": lsZ,
        "off_shell_effective_action_difference_quartic": difference4,
        "scope": "Gaussian moments are independent finite-dimensional diagnostic controls, not a substitute for the actual derivative-map Wick calculation. The continuous formal identity uses transformed action, Jacobian and sources. Pole factorization applies to the stable light particle and the full transformed action only; ordinary off-shell effective actions are not scalars under nonlinear maps.",
        "checks": checks,
    }
