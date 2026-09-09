"""Independent Fraction reconstruction; no SymPy or production imports."""

from fractions import Fraction as F


def replay():
    r = F(1, 100)
    d = (F(21, 10), F(9), F(70), F(800))
    c = (F(4), F(128), F(16384), F(1048576))

    def pieces(v):
        h, hd, hdd, h3 = v
        return (
            3 * (h**4 + 2 * h * h * hd),
            18 * h3 + 90 * hd * hd + 90 * h * hdd + 108 * h * h * hd,
        )

    qp = F(3, 4) * d[0] ** 2 + F(3, 2) * d[1]
    qf = F(3, 4) * c[0] ** 2 + F(3, 2) * c[1]
    p = (F(7, 2) + F(11, 5) * d[0] * r + F(2, 3) * qp * r * r) ** 2 / r**3
    f = (F(7, 2) * (1 + 2 * c[0]) + F(25, 12) * qf) ** 2
    vp, vf = pieces(d), pieces(c)
    C0 = F(7, 18) * (p + f) + F(13, 12600) * r * vp[0] + F(13, 1080) * vf[0]
    Cb = F(13, 12600) * r * vp[1] + F(13, 1080) * vf[1]
    Wp = (F(11, 10) + d[0] * r) ** 2 / r
    Wf = (F(11, 10) + F(27, 20) * c[0]) ** 2
    costs = {
        "scalar_C0": C0,
        "scalar_Cbeta": Cb,
        "past_Wick_square_cost": Wp,
        "future_Wick_square_cost": Wf,
        "Wick_square_cost": Wp + Wf,
        "past_caps": d,
        "future_caps": c,
        "actual_past_ratio": r,
        "positive_C0_margin": 5000000 - C0,
        "positive_Cbeta_margin": 5000000 - Cb,
    }
    delta, zeta, sigma = F(1, 10**8), F(1, 5000), F(5)
    gain = 3 * F(19, 10) + F(39, 35) * F(19, 10) ** 2 * r
    margin = (
        gain
        - F(18, 5)
        - 5000000 * delta
        - (Wp + Wf) * zeta
        - F(13, 35) * (1 + r) * sigma
    )
    gate = {
        "weighted_delta": delta,
        "zeta_upper_kappa_Phi_square_over_three": zeta,
        "sigma": sigma,
        "strict_focusing_margin": margin,
        "margin_above_nine_fiftieths": margin - F(9, 50),
        "future_geometry_and_state_amplitude_are_conditional_hypotheses_not_verified_from_past": True,
    }
    lam = delta / 360
    response = (F(16, 25), F(1728, 25), F(155136, 25), F(14770176, 25))
    errors = [v * lam for v in response]
    caps = (F(1, 100), F(1, 2), F(3), F(16))
    zeta2 = F(20, 9) * delta * (1 - 4 * lam)
    past = {
        "delta": delta,
        "lambda": lam,
        "past_x": [-r, F(0)],
        "backward_y_bounds": [F(50, 27), F(2)],
        "a_past_lower": F(1),
        "C3_errors": errors,
        "strict_C3_margins": [v - e for v, e in zip(caps, errors, strict=True)],
        "anchor_H0_times_tau": F(2),
        "initial_and_past_zeta_squared_upper": zeta2,
        "strict_squared_field_gate_margin": zeta * zeta - zeta2,
        "sigma": F(0),
        "beta_S": F(0),
        "actual_scalar_state_and_SEE_past": True,
        "future_state_cap_proved_on_every_shorter_segment": False,
    }
    thermal = {
        "actual_thermal_past_at_delta_max": past,
        "lambda_to_delta": F(1, 360),
        "generic_ODE_comparison_coupling_margin": F(1, 16) - lam,
        "state_Wick_to_energy_coefficient_without_hbar_over_pi_squared": F(5, 24),
        "endpoint_zeta_squared": F(50),
        "endpoint_x_strict_bounds": [F(1, 16), F(1, 4)],
        "thermal_branch_is_not_a_new_QEI_only_endpoint_argument": True,
        "future_Wick_bound_is_a_conditional_theorem_hypothesis": True,
    }
    assert C0 < 5000000 and Cb < 5000000
    assert margin == F(63325013, 350000000) > F(9, 50)
    assert min(past["strict_C3_margins"]) > 0
    assert past["strict_squared_field_gate_margin"] > 0
    return {"costs": costs, "gate": gate, "thermal": thermal}


def spectral_moment(xi):
    """Integrate k^3 and k(theta-k)^2 by their monomial coefficients."""
    xi = F(xi)
    radial = (1 - 2 * xi) / 4 + 2 * xi * (F(1, 2) - F(2, 3) + F(1, 4))
    return radial / 4
