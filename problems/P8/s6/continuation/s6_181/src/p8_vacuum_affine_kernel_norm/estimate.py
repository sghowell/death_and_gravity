"""Full-band summation, same-kernel identification and current inverse transfer."""

from functools import cache

import sympy as s
from p8_proca_rank_one_inverse import kernel
from p8_vacuum_affine_coupled_response import inverse as previous_inverse

from . import tail, threshold

K_UPPER = s.Integer(300000000000)


@cache
def data():
    small_sum = 24 * threshold.C_SMALL
    tail_sum = 48 * tail.C_TAIL
    direct = 35 * (small_sum + tail_sum)
    C, gamma = s.symbols("C gamma", positive=True)
    a0 = s.Rational(15625, 6144)
    beta = C * (a0 + K_UPPER / gamma)
    weight = (4 * beta + 1) ** 2
    m, t, omega, p = s.symbols("m t omega p", positive=True)
    original = previous_inverse.data()
    params = original["weak_log_majorant_parameters"]
    mapped = {params["C"]: C, params["K_L1"]: K_UPPER, params["gamma"]: gamma}
    checks = {
        "small_band_total": small_sum - 4800000000,
        "tail_band_total": tail_sum - 960000000,
        "all_bands_numeric_norm": direct - 201600000000,
        "same_parent_inverse_weight_substitution": s.factor(
            original["explicit_exponential_weight"].subs(mapped) - weight
        ),
        "same_parent_inverse_beta_substitution": s.factor(
            original["integrable_composition_constant"].subs(mapped) - beta
        ),
        "mass_fraction_rescaling": s.factor(
            1 - 4 * m * m / (m * omega) ** 2 - (1 - 4 / omega**2)
        ),
        "mass_kernel_integrated_scaling": m / m - 1,
        "band_Laplace_normalization": s.factor(
            s.integrate(-2 * s.exp(-p * t) * s.sin(omega * t), (t, 0, s.oo))
            + 2 * omega / (p * p + omega * omega)
        ),
        "same_source_zero_instantaneous": kernel.data()["instantaneous_inverse"],
        "same_source_positive_static_moment": kernel.data()["spectral_static_moment"]
        - s.Rational(1, 4),
    }
    return {
        "threshold_band_sum_upper": small_sum,
        "tail_band_sum_upper": tail_sum,
        "unrounded_full_half_line_kernel_L1_upper": direct,
        "rounded_full_half_line_kernel_L1_upper": K_UPPER,
        "mass_scope": "For every fixed m>0: K_m(t)=m K_1(mt), so the complete half-line L1 norm is mass independent; this includes the same current mass1000",
        "identification": "Absolute L1 summation of all dyadic sine kernels; their Laplace transforms sum by the nonnegative spectral density and partition to the source-pinned R(p^2)=1/F_m(p^2). Uniqueness of Laplace transforms identifies the same retained inverse, with no pole, tail or instantaneous term deleted.",
        "current_inverse_beta_upper": beta,
        "current_inverse_weight_choice": weight,
        "current_adapted_inverse_C0_formula": 2
        * s.exp(weight)
        * s.Max(a0, K_UPPER / gamma),
        "remaining_unevaluated_constant": "The full current curved weak-log remainder majorant C, including its matching and fixed local contacts, is still not numerically bounded here. This scalar bound does not establish a small or stable full inverse.",
        "local_lapse_recovery": "The already quantified actual C1<15000 is unchanged, not replaced by the scalar kernel norm.",
        "checks": checks,
        "gates": {
            "positive_numeric_kernel_upper": direct > 0,
            "rounded_bound_strictly_larger": direct < K_UPPER,
            "same_old_kernel_L1_object": kernel.data()["L1_half_line_kernel"] is True,
        },
    }
