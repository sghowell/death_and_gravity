"""Literal all-time action normalization and center fourteen-variable CCR bridge."""

from functools import cache

import sympy as sp
from p8_coupled_momentum import constraints, quadratic
from p8_proca_physical_matter import bridge as physical_bridge

from . import model


@cache
def data():
    D = sp.diag(
        physical_bridge.data()["old_density_to_fixed_spatial_phase"], sp.eye(10)
    )
    OmegaQ = sp.diag(
        model.J2, model.J2, model.J2, model.J2, model.data()["Proca_symplectic_form"]
    )
    k = constraints.wave_scale
    cases = {
        ("tensor", "tensor"): k * k / 2,
        ("tensor_p", "tensor_p"): 2,
        ("tensor2", "tensor2"): k * k / 2,
        ("tensor2_p", "tensor2_p"): 2,
        ("Px", "Px"): 10**6 + k * k,
        ("Py", "Py"): sp.Integer(10**6),
        ("Pz", "Pz"): sp.Integer(10**6),
        ("Wx", "Wx"): sp.Integer(1),
        ("Wy", "Wy"): 1 + k * k / 10**6,
        ("Wz", "Wz"): 1 + k * k / 10**6,
        ("tensor", "tensor_p"): sp.Integer(0),
        ("tensor", "matter"): sp.Integer(0),
        ("tensor_p", "Py"): sp.Integer(0),
        ("Wx", "Py"): sp.Integer(0),
    }
    checks = {
        "all_fourteen_physical_phase_CCR_normalizations": model.clean(
            D * model.data()["actual_full_symplectic_form"] * D.T - OmegaQ
        )
    }
    for pair, value in cases.items():
        checks["literal_all_time_quadratic_" + pair[0] + "_" + pair[1]] = sp.factor(
            quadratic.pair(*pair)["actual"] - value
        )
    # Restore physical orthonormal fields at a fixed time, retaining
    # coordinate volume a^3 and density momenta.
    a, m = model.a, model.m
    T, P, W, Pi, kc = sp.symbols("T P W Pi nonzero_comoving_k", real=True)
    tensor = a**3 * ((P / a**3) ** 2 + (kc / a) ** 2 * T * T / 4)
    long = a**3 * ((m * m + (kc / a) ** 2) * (Pi / a**2) ** 2 + (W / a) ** 2) / 2
    trans = (
        a**3
        * (m * m * (Pi / a**2) ** 2 + (1 + (kc / a) ** 2 / m**2) * (W / a) ** 2)
        / 2
    )
    checks.update(
        {
            "actual_tensor_density_and_wave_normalization": sp.factor(
                tensor - P * P / a**3 - a * kc * kc * T * T / 4
            ),
            "actual_longitudinal_Proca_density_normalization": sp.factor(
                long - (m * m / a + kc * kc / a**3) * Pi * Pi / 2 - a * W * W / 2
            ),
            "actual_transverse_Proca_density_normalization": sp.factor(
                trans
                - m * m * Pi * Pi / (2 * a)
                - (a + kc * kc / (a * m * m)) * W * W / 2
            ),
        }
    )
    return {
        "full_density_to_S6_102_spatial_phase": D,
        "physical_stress_bridge_time": sp.Integer(0),
        "full_fixed_spatial_phase_symplectic_form": OmegaQ,
        "actual_sector_order": "old scalar density(4), tensor T/PT(2), tensor S/PS(2), Cartesian W/Pi(6)",
        "checks": checks,
    }
