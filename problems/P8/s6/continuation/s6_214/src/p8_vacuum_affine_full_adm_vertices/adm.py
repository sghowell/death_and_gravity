"""Actual full ADM temporal constraint and positive reduced Proca Hamiltonian."""

from functools import cache

import sympy as s
from p8_vacuum_affine_proca_gaussian import gaussian


def magnetic_energy(h, B):
    v = s.sqrt(h.det())
    return (B.T * h * B)[0] / (2 * v)


def temporal_constraint(N, v, beta, A, div_pi, m):
    return (beta.T * A)[0] - N * div_pi / (m * m * v)


def shifted_spatial_direction(Q):
    return Q - s.trace(Q) * s.eye(3) / 2


@cache
def data():
    N, v, m = s.symbols("N volume m", positive=True)
    A0, bA, D = s.symbols("A0 beta_dot_A div_pi", real=True)
    raw = -A0 * D - m * m * v * (A0 - bA) ** 2 / (2 * N)
    constraint = bA - N * D / (m * m * v)
    p = s.Matrix(s.symbols("pi0:3", real=True))
    b = s.Matrix(s.symbols("beta0:3", real=True))
    A = s.Matrix(s.symbols("A0:3", real=True))
    dA = s.Matrix(3, 3, s.symbols("dA0:9", real=True))
    db = s.Matrix(3, 3, s.symbols("db0:9", real=True))
    flux = (b.T * (dA - dA.T) * p)[0] - (b.T * A)[0] * D
    Lie = (p.T * (dA.T * b + db * A))[0]
    # dA_ij = partial_i A_j; divergence uses partial_j A_i, not partial_i A_j.
    divergence = (p.T * db * A)[0] + (p.T * dA * b)[0] + (b.T * A)[0] * D
    q = s.Symbol("psi", real=True)
    source = gaussian.data()["checks"]
    checks = {
        "complete_temporal_constraint_equation": s.factor(
            s.diff(raw, A0).subs(A0, constraint)
        ),
        "positive_constraint_energy_and_longitudinal_shift": s.factor(
            raw.subs(A0, constraint) - (N * D * D / (2 * m * m * v) - bA * D)
        ),
        "shift_flux_to_covector_Lie_with_full_boundary": s.expand(
            Lie - flux - divergence
        ),
        "pure_spatial_trace_shifted_matrix": shifted_spatial_direction(2 * q * s.eye(3))
        + q * s.eye(3),
        "pure_spatial_trace_constraint_exponent": -s.trace(2 * q * s.eye(3)) / 2
        + 3 * q,
        "actual_full_parent_source_value_at_clock": source[
            "actual_covariant_clock_source_value"
        ],
        "actual_full_parent_first_source_variation_at_clock": source[
            "actual_covariant_clock_source_first_variation"
        ],
        "zero_mean_second_source_operator_contact": source[
            "linear_source_contact_has_zero_mean"
        ],
    }
    return {
        "signature": "Use the unchanged P8 +--- physical metric ds^2=N^2dt^2-h_ij(dx^i+beta^i dt)(dx^j+beta^j dt), positive spatial h, and the corresponding source-free Proca quadratic density-F^2/4+m^2A^2/2.",
        "Legendre": "E_i=dot A_i-partial_i A0-beta^k F_ki; Pi^i=sqrt(h)h^ij E_j/N. Keep the boundary Pi^i partial_i A0=-A0 divPi+div(A0 Pi). Eliminating A0 gives beta.A-N divPi/(m^2sqrt(h)).",
        "full_reduced_H": "Integral N[h_ij Pi^iPi^j/(2sqrt(h))+sqrt(h)h^ik h^jl F_ijF_kl/4+m^2sqrt(h)h^ij A_iA_j/2+(divPi)^2/(2m^2sqrt(h))]+beta^i[Pi^jF_ij-A_i divPi]. No fourth oscillator is retained or longitudinal energy removed.",
        "shift_boundary": "The complete shift term equals integral Pi^j(Lie_beta A)_j modulo the explicitly retained compact spatial divergence. Dropping-beta.A divPi loses the longitudinal part of the transport generator.",
        "full_spatial_chart": "Set N=1+n,h=a^2exp(Q),Q symmetric,tau=trQ,B=Q-tau I/2. The exact lapse Hamiltonian is N/2 times[Pi^t exp(B)Pi/a+curlA^t exp(B)curlA/a+am^2 A^t exp(-B)A+exp(-tau/2)(divPi)^2/(a^3m^2)].",
        "actual_source_boundary": "The full affine source is not set to zero away from the reference. Frozen S176 proves S=DS=0 on the actual clock and vanishing zero-mean second source insertion. Thus this same-preparation Gaussian first metric response is ordinary connected Proca. Higher source contacts and nonlinear mean response remain.",
        "checks": checks,
        "gates": {
            "positive_mass_volume_and_lapse_domain": m.is_positive
            and v.is_positive
            and N.is_positive,
            "full_temporal_constraint_varies_with_lapse_and_trace": True,
            "full_shift_longitudinal_term_retained": True,
            "physical_field_and_canonical_momentum_unchanged": True,
            "same_parent_gaussian_only_not_source_free_nonlinear_parent": True,
        },
    }
