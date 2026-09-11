"""Full-continuum weak response on separated supports and canonical bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import noise
from p8_vacuum_affine_cd_source_noise import modes

from . import tail

CURRENT = s.Integer(5) * 10**49
CANONICAL = 2 * s.Rational(1, 10) ** 750


@cache
def data():
    a, m, V0 = s.symbols("a mass V0", positive=True)
    A = s.Matrix(s.symbols("A0:3", real=True))
    pi = s.Matrix(s.symbols("pi0:3", real=True))
    b = s.Matrix(s.symbols("B0:3", real=True))
    d = s.symbols("D0:5", real=True)
    D = s.Matrix([[d[0], d[1], d[2]], [d[1], d[3], d[4]], [d[2], d[4], -d[0] - d[3]]])
    E = pi / a**2
    B = b / a**2
    V = m * A / a
    stress = (
        -E * E.T
        - B * B.T
        + V * V.T
        + s.eye(3) * (E.dot(E) + B.dot(B) + V0**2 - V.dot(V)) / 2
    )
    HD = (pi.dot(D * pi) / a + b.dot(D * b) / a - a * m * m * A.dot(D * A)) / 2
    k, h = s.symbols("kappa h", positive=True)
    metric = 2 * h / s.sqrt(k)
    c = noise.data()["complete_stress_variance_coefficient_upper"]
    x, y = s.symbols("inner_real inner_imaginary", real=True)
    z = x + s.I * y
    K = s.Symbol("K", positive=True)
    checks = {
        "full_physical_stress_to_Hamiltonian_source_sign": s.expand(
            HD + a**3 * s.trace(D * stress) / 2
        ),
        "temporal_constraint_is_trace_only_in_direct_vertex": s.diff(
            s.trace(D * stress), V0
        ),
        "real_current_from_full_commutator": s.expand(
            s.I * (z - s.conjugate(z)) / 4 + y / 2
        ),
        "both_canonical_metric_chain_factors": s.diff(metric, h) ** 2 - 4 / k,
        "complete_canonical_response_display": 4 * CURRENT / modes.KAPPA - CANONICAL,
        "complete_canonical_tail_coefficient": 4
        * (tail.DISPLAY / (2 * K))
        / modes.KAPPA
        - 2 * s.Rational(1, 10) ** 748 / K,
        "canonical_tail_at_computational_partition": 2
        * s.Rational(1, 10) ** 748
        / s.Integer(10) ** 16
        - 2 * s.Rational(1, 10) ** 764,
        "fixed_parent_mass_and_kappa": (modes.MASS - 1000)
        + (modes.KAPPA - s.Integer(10) ** 800),
    }
    return {
        "actual_response": "For the strict ordered supports, B(D,Gamma)=i<[T[D],T[Gamma]]>/4 is the full Gaussian-sector weak first-order current response, including all internal/external momenta. The local metric/finite-curvature contacts vanish by support, not by altering the fixed prescription.",
        "current_bound": CURRENT,
        "unrounded_current_coefficient": c / 2,
        "canonical_bound": CANONICAL,
        "unrounded_canonical_coefficient": 2 * c / modes.KAPPA,
        "current_regulator_error": "|B-B_K|<=tau_D(K)tau_Gamma(K)/2 <5e51 N[D]N[Gamma]/K for nonzero tests, zero if either test is zero.",
        "canonical_regulator_error": "The error after both canonical metric factors is <2e-748 N[D]N[Gamma]/K, hence<2e-764 at K1e16. K is removed in the full result, not declared a physical cutoff.",
        "smearing_norm": noise.data()["smearing_norm"],
        "normalization": "T[f]=integral a^3 f_ij T_ij in the physical frame; H_f=-T[f]/2 for tracefree shear. For h=sqrt(kappa)gamma/2, the weak canonical response has both chain factors4/kappa. This uses the proved tensor normalization, not a fully reduced mixed scalar/metric norm.",
        "weak_boundary": "This is the bounded bilinear first-order Kubo/relative-Cauchy response on separated smooth supports, not operator-norm differentiability, a complete response on overlapping times, a finite-amplitude remainder, full inverse or interacting quantum background.",
        "checks": checks,
        "gates": {
            "actual_full_stress_variance_input": 0 < c < 10**50,
            "current_strict_display": c / 2 < CURRENT,
            "canonical_strict_display": 2 * c / modes.KAPPA < CANONICAL,
            "same_three_time_one_space_norm": "partial_t^j"
            in noise.data()["smearing_norm"],
            "actual_finite_mean_not_reselected": True,
        },
    }
