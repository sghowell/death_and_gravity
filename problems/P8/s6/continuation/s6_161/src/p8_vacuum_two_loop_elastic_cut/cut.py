"""Physical-channel ownership and two-loop optical-theorem normalization."""

from functools import cache

import sympy as s


def rows():
    return [
        {
            "state": "two_Phi",
            "first_loop_order": 1,
            "through_two_loops": "TREE_SQUARE_AND_TREE_ONE_LOOP_INTERFERENCE",
        },
        {
            "state": "three_Phi",
            "first_loop_order": None,
            "through_two_loops": "FORBIDDEN_BY_EXACT_PHI_PARITY",
        },
        {
            "state": "heavy_or_fermionic_states",
            "first_loop_order": None,
            "through_two_loops": "KINEMATICALLY_ABOVE_THE_WINDOW",
        },
        {
            "state": "two_gauge_bosons",
            "first_loop_order": 3,
            "through_two_loops": "ONE_FERMION_LOOP_ON_EACH_SIDE_FIRST_CONTRIBUTES_AT_THREE",
        },
        {
            "state": "four_or_more_Phi",
            "first_loop_order": 3,
            "through_two_loops": "AT_LEAST_THREE_PHASE_SPACE_LOOP_WEIGHTS",
        },
        {
            "state": "free_decoupled_spectators",
            "first_loop_order": None,
            "through_two_loops": "NO_SCALAR_SCATTERING_COUPLING_IN_THIS_VACUUM",
        },
    ]


def validate_rows(value):
    if value != rows():
        raise ValueError("The fixed low-energy intermediate-state ownership changed")
    return True


@cache
def data():
    h, A0, R, I, beta, pi, lam, B = s.symbols(
        "h A0 ReA1 ImA1 beta pi lambda B1", real=True
    )
    modulus = s.expand((A0 + h * (R + s.I * I)) * (A0 + h * (R - s.I * I)))
    rho = s.expand(beta * modulus / (32 * pi))
    z = s.Symbol("s", positive=True)
    kernel = s.integrate((z - 2) ** -3, (z, 4, 6))
    rho2_upper = s.Rational(3, 5) * 73 * lam * B / (16 * 3)
    cut2_upper = s.Rational(2, 3) * kernel * rho2_upper
    return {
        "intermediate_state_ownership": rows(),
        "two_particle_loop_degree": "L_total=L_left+L_right+n_intermediate-1. Through L_total=2, only the two-Phi tree/one-loop interference is present on this window.",
        "density_one": "rho1(s)=beta/(32 pi) Integral_0^1 A0(s,z)^2 dz",
        "density_two": "rho2(s)=beta/(16 pi) Integral_0^1 A0(s,z) Re A1(s,z) dz",
        "two_loop_density_absolute_upper": rho2_upper,
        "two_loop_integrated_cut_absolute_upper": cut2_upper,
        "formal_cut_definition": "I(h)=h I1+h^2 I2, Ij=(2/pi) Integral_4^6 rhoj(s)/(s-2)^3 ds. Do not assert rho2>=0 or replace the full higher-order spectrum by this truncated polynomial.",
        "checks": {
            "tree_one_loop_interference_factor_two": modulus.coeff(h, 1) - 2 * A0 * R,
            "first_cut_identical_particle_normalization": rho.coeff(h, 0)
            - beta * A0 * A0 / (32 * pi),
            "second_cut_identical_particle_normalization": rho.coeff(h, 1)
            - beta * A0 * R / (16 * pi),
            "one_loop_square_is_next_cut_order": modulus.coeff(h, 2) - R * R - I * I,
            "selected_crossing_kernel_weight": kernel - s.Rational(3, 32),
            "uniform_second_density_bound": rho2_upper - s.Rational(73, 80) * lam * B,
            "integrated_second_cut_bound": cut2_upper - s.Rational(73, 1280) * lam * B,
            "two_loop_cut_relative_to_tree": cut2_upper / (4 * lam)
            - s.Rational(73, 5120) * B,
            "two_loop_two_particle_cut_order": 0 + 1 + 2 - 1 - 2,
            "first_gauge_pair_cut_order": 1 + 1 + 2 - 1 - 3,
            "first_four_Phi_cut_order": 0 + 0 + 4 - 1 - 3,
        },
    }
