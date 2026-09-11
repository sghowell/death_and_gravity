"""Actual canonical measure, arbitrary-lapse modes and physical metric readouts."""

from functools import cache

import sympy as s
from p8_constant_proca import quantum as ordinary
from p8_vector_state import wkb
from p8_vector_variation import modes

ZETA = s.Rational(1, 10**6)
KAPPA = s.Integer(10) ** 800
MASS = s.Integer(1000)


@cache
def data():
    kappa, zeta = s.symbols("kappa zeta", positive=True)
    aa, ss, ff = s.symbols("canonical_A source_S canonical_F_squared", real=True)
    W = aa / s.sqrt(kappa * zeta)
    J = s.sqrt(kappa / zeta) * ss
    original = -kappa * zeta * ff / (4 * kappa * zeta) + kappa * (W - ss) ** 2 / 2
    canonical = -ff / 4 + aa**2 / (2 * zeta) - J * aa + kappa * ss**2 / 2
    n, a, m, k = s.symbols("N a m k", positive=True)
    velocity, A, A0, pi = s.symbols("A_dot A A0 pi_A", real=True)
    q = k * k / (a * a)
    g2 = {"transverse": a, "longitudinal": a * m * m / (q + m * m)}
    longitudinal = (
        a * (velocity - k * A0) ** 2 / (2 * n)
        + a**3 * m * m * A0 * A0 / (2 * n)
        - n * a * m * m * A * A / 2
    )
    solved = s.solve(s.diff(longitudinal, A0), A0)[0]
    reduced = s.factor(longitudinal.subs(A0, solved))
    checks = {
        "canonical_source_and_contact_normalization": s.simplify(original - canonical),
        "chosen_mass_matches_earlier_ordinary_sector": 1 / ZETA - MASS**2,
        "chosen_zeta_within_frozen_parent_range": ZETA - s.Rational(1, 1000000),
        "literal_temporal_constraint": s.factor(
            solved - k * velocity / (k * k + a * a * m * m)
        ),
        "literal_longitudinal_reduced_action": s.factor(
            reduced
            - (g2["longitudinal"] * velocity**2 / (2 * n) - n * a * m * m * A * A / 2)
        ),
    }
    packets = {}
    v, p = s.symbols("v p", real=True)
    for kind, g in g2.items():
        ham = n * (pi * pi / (2 * g) + g * (q + m * m) * A * A / 2)
        replacements = {A: v / s.sqrt(g), pi: s.sqrt(g) * p}
        energy = s.factor(s.diff(ham, n).subs(replacements) / a**3)
        pressure = s.factor((-s.diff(ham, a) / (3 * n * a * a)).subs(replacements))
        z = q / (q + m * m)
        expected_p = (
            (p * p + (q - m * m) * v * v)
            if kind == "transverse"
            else ((1 + 2 * z) * p * p - (q + m * m) * v * v)
        ) / (6 * a**3)
        checks[kind + "_metric_energy_before_canonical_substitution"] = s.factor(
            energy - (p * p + (q + m * m) * v * v) / (2 * a**3)
        )
        checks[kind + "_metric_pressure_before_canonical_substitution"] = s.factor(
            pressure - expected_p
        )
        # The transform is symplectic at fixed metric, not just a mode scaling.
        C = s.diag(1 / s.sqrt(g), s.sqrt(g))
        Omega = s.Matrix([[0, 1], [-1, 0]])
        checks[kind + "_canonical_symplectic_map"] = s.simplify(C.T * Omega * C - Omega)
        Hsym = s.symbols("H", real=True)
        rate = s.factor(s.diff(s.log(g), a) * a * Hsym / 2)
        desired = Hsym * (s.Rational(1, 2) + (z if kind == "longitudinal" else 0))
        checks[kind + "_full_fixed_comoving_canonical_rate"] = s.factor(rate - desired)
        actual_matrix = s.hessian(pressure, (v, p))
        mapping = {
            a: modes.scale,
            m: modes.mass2 ** s.Rational(1, 2),
            k: modes.scale * s.sqrt(modes.q),
        }
        checks[kind + "_old_ordinary_pressure_readout_bridge"] = (
            actual_matrix.subs(mapping)
            - ordinary.matrices()[kind]["new_canonical_pressure"]
        ).applyfunc(s.simplify)
        checks[kind + "_old_ordinary_energy_readout_bridge"] = (
            s.hessian(energy, (v, p)).subs(mapping)
            - ordinary.matrices()[kind]["new_canonical_energy"]
        ).applyfunc(s.simplify)
        packets[kind] = {
            "normalizer_squared": g,
            "physical_Hamiltonian": ham,
            "energy": energy,
            "pressure": pressure,
            "canonical_rate": rate,
        }
    return {
        "zeta": ZETA,
        "kappa": KAPPA,
        "canonical_mass": MASS,
        "canonical_field": "A_mu=sqrt(kappa*zeta) W_mu; the source is J_mu=sqrt(kappa/zeta) S_mu",
        "quantization_choice": "Canonical three-polarization Proca phase measure and the same covariant dimensional finite prescription as S6.82, mu=m=1000. This defines only the retained conditional Gaussian vector sector.",
        "signature_bridge": "g_old=-g_P8, with the same lower A components, connection, volume, positive physical mode energy and canonical phase covariance; the ordinary actions coincide.",
        "complete_canonical_source": J,
        "modes": packets,
        "checks": checks,
    }


@cache
def clock():
    u = wkb.u
    a = (1 + u * u) ** 2
    H = s.diff(a, u) / a
    k, m = s.symbols("k m", positive=True)
    q = k * k / a**2
    z = q / (q + m * m)
    checks = {}
    for kind, g2 in (("transverse", a), ("longitudinal", a * m * m / (q + m * m))):
        d = s.diff(s.log(g2), u) / 2
        U = s.factor(s.diff(d, u) + d * d)
        checks[kind + "_actual_clock_operator_equals_frozen_state_operator"] = s.factor(
            U - wkb.frequency(kind)["U"].subs(wkb.z, z)
        )
    checks["actual_clock_metric_and_old_Hubble"] = s.factor(H - wkb.background()["H"])
    return {
        "a": a,
        "H": H,
        "mass_time": MASS,
        "tau": s.S.One,
        "M_tau": s.sqrt(KAPPA),
        "checks": checks,
    }
