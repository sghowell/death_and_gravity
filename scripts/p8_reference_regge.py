"""PRESCRIPTION-2 reduction and a finite-window/Regge gravity interface.

Read-only exact diagnostic. The written proofs, not finite fixtures, carry
the functional jet and all-energy claims. No common quantum reference,
Regge envelope, physical error, or replacement model is silently supplied.
"""

import hashlib
import json
from fractions import Fraction as F
from math import factorial

import sympy as s
from p8_match1_rate_input import REPO, exact
from p8_prescription2_gravity import (
    DETECTOR,
    kernel_error_bound,
    radial_moment,
)
from p8_prescription2_gravity import (
    inputs as previous_inputs,
)
from p8_rate4_candidate import original_borns

EXTRA_SOURCES = {
    "scripts/p8_prescription2_gravity.py": "5bff460ef3cc0f3e4fb293bdb54bfee76668be72dc8c47074d9adc7edf77ad1d",
    "docs/assessment-2026-09-22-p8-prescription2-gravity.md": "cadecb406ae0c5ede7b58839ae19f31c474647b55c644fd45d7b841579dc0d3f",
}
AFFINE_REPORT = "problems/P8/s6/continuation/s6_174/certificates/polynomial-vacuum-analytic-affine-parent.json"
AFFINE_SHA = "e0cf546e77cf12ff4e0f2e6bef706ff6de1eddc40cacff3b0b151fcb7009a4f4"


def inputs():
    manifest = previous_inputs()
    for name, expected in {**EXTRA_SOURCES, AFFINE_REPORT: AFFINE_SHA}.items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected
        assert name not in manifest or manifest[name] == expected
        manifest[name] = expected
    root = (REPO / AFFINE_REPORT).parent.parent
    for name, expected in json.loads((REPO / AFFINE_REPORT).read_text())[
        "source_sha256"
    ].items():
        path = (root / name).resolve()
        assert path.is_relative_to(root)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
        relative = str(path.relative_to(REPO))
        assert relative not in manifest or manifest[relative] == expected
        manifest[relative] = expected
    return manifest


def log_inverse_upper(detector):
    """Rational bound log(1/E) <= (7/10)*ceil(log2(1/E)), no float log."""
    detector = exact(detector)
    if not 0 < detector < 1:
        raise ValueError("Require 0<E<1")
    scaled, count = detector, 0
    while scaled < 1:
        scaled *= 2
        count += 1
    return F(7 * count, 10)


def regge_tail_bound(*, detector, split, residue_over_slope, remainder_arc_error):
    """Conditional normalized tail bound for a supplied finite-angle envelope.

    For BOTH cuts at s>=T and E<=q<=1 require
    |rho(s,-q^2)| <= C*s^2*(s/T)^(-alpha*q^2) + residual,
    with alpha>0. residue_over_slope is C/alpha, NOT inferred from kappa.
    remainder_arc_error bounds the complete normalized integral of residual.
    Does not bound the large complex contour or observable conversion.
    """
    detector, split, ratio, remainder = map(
        exact, (detector, split, residue_over_slope, remainder_arc_error)
    )
    if not 0 < detector < 1 or split < 32 or ratio < 0 or remainder < 0:
        raise ValueError("Require 0<E<1, T>=32 and explicit nonnegative bounds")
    # P_E = integral_E^1 phi(q)/q dq. The log majorant is proved separately.
    pole_upper = (
        log_inverse_upper(detector) - F(4, 3) + 3 * detector / 2 - detector**3 / 6
    )
    assert pole_upper > 0
    # 2/pi < 2/3. Both unequal crossed denominators are retained in the proof.
    return (
        F(2, 3)
        * (split / (split - 3)) ** 3
        * ratio
        * pole_upper
        / radial_moment(0, lower=detector)
        + remainder
    )


def finite_window_error(
    *,
    detector,
    split,
    finite_forward_moment,
    residue_over_slope,
    regge_remainder_error,
    contour_error,
    observable_error,
):
    """Conditional a_E >= -returned_error; every physical input is required.

    finite_forward_moment = (2/pi)*int_[32,T) rho_positive(s,1)/(s-3)^3 ds.
    All errors use the normalized arc. observable_error must include the
    positive-reference mismatch and low-cut/IR conversions without overlaps.
    This interface checks numbers, not their physical provenance.
    """
    detector, split, moment, ratio, remainder, contour, observable = map(
        exact,
        (
            detector,
            split,
            finite_forward_moment,
            residue_over_slope,
            regge_remainder_error,
            contour_error,
            observable_error,
        ),
    )
    if split < 32 or min(moment, ratio, remainder, contour, observable) < 0:
        raise ValueError("Require T>=32 and explicit nonnegative complete inputs")
    if split == 32 and moment != 0:
        raise ValueError("An empty finite window must have zero moment")
    return (
        kernel_error_bound(
            detector=detector,
            uv_s_min=32,
            positive_forward_moment=moment,
            other_error=0,
        )
        + regge_tail_bound(
            detector=detector,
            split=split,
            residue_over_slope=ratio,
            remainder_arc_error=remainder,
        )
        + contour
        + observable
    )


def audit():
    before = inputs()
    # Inherited actual rank-60 quotient, not a claim to rederive its 60x60 matrix.
    p = s.Symbol("p", positive=True)
    tau, sigma = 3 * (2 * p**3 - 1) / p, (8 * p + 5) / (8 * p**2)
    D, eta = s.diag(tau, sigma, sigma, sigma), s.diag(1, -1, -1, -1)
    Xi = eta - D.inv()
    assert s.simplify((s.eye(4) + Xi * D).det() + tau * sigma**3) == 0
    assert (
        4**3 - 4 - 4 == 56
    )  # Gauge directions and RETAINED trace, not 60 auxiliaries.
    assert 2 * F(3, 5) ** 3 < 1  # tau<0 throughout 1/3<p<3/5.

    # Exact block Gaussian elimination, including indefinite algebraic blocks.
    gaussian_checks = 0
    for C, B, A in (
        (s.diag(2, 3), s.Matrix([[1, 2], [3, -1]]), s.Matrix([[11, 2], [2, 17]])),
        (s.diag(-2, 5), s.Matrix([[2, -3], [1, 4]]), s.Matrix([[13, 1], [1, 19]])),
        (s.Matrix([[3, 1], [1, 2]]), s.Matrix([[1], [2]]), s.Matrix([[7]])),
    ):
        n, m = C.rows, A.rows
        full = C.row_join(B).col_join(B.T.row_join(A))
        transform = (
            s.eye(n).row_join(-C.inv() * B).col_join(s.zeros(m, n).row_join(s.eye(m)))
        )
        reduced = A - B.T * C.inv() * B
        assert transform.det() == 1
        assert transform.T * full * transform == s.diag(C, reduced)
        assert full.det() == C.det() * reduced.det()
        gaussian_checks += 1
    k, c, mixing, jac = s.symbols("k c mixing jac", nonzero=True)
    assert s.Matrix([[jac, mixing * k], [0, 1]]).det() == jac
    # Integrating an algebraic block cannot discard the retained kinetic block.
    full = s.Matrix([[c, mixing], [mixing, k**2 + 1]])
    assert s.expand(full.det() / c - (k**2 + 1 - mixing**2 / c)) == 0
    assert s.diff(full.det() / c, k) == 2 * k
    # Explicit canonical algebraic toy: sqrt(det constraint bracket)=|C|.
    C = s.Matrix([[3, 1], [1, 2]])
    brackets = s.zeros(2).row_join(-C).col_join(C.row_join(s.zeros(2)))
    assert brackets.det() == C.det() ** 2

    # Heavy Gaussian Schur formula with field-dependent operator and source.
    x, H = s.symbols("x H")
    K, J, base = 2 + x + x**2, x**2 + 2 * x**3, x**2 / 2 + x**4
    action = base - K * H**2 / 2 + H * J
    mean = J / K
    reduced = base + J**2 / (2 * K)
    assert s.cancel(action.subs(H, mean) - reduced) == 0
    schur = (
        s.diff(action, x, 2) - s.diff(action, x, H) ** 2 / s.diff(action, H, 2)
    ).subs(H, mean)
    covariant_formula = (
        s.diff(base, x, 2)
        + mean * s.diff(J, x, 2)
        - mean**2 * s.diff(K, x, 2) / 2
        + (s.diff(J, x) - s.diff(K, x) * mean) ** 2 / K
    )
    assert s.cancel(schur - s.diff(reduced, x, 2)) == 0
    assert s.cancel(schur - covariant_formula) == 0
    assert s.cancel(covariant_formula - s.diff(base, x, 2)) != 0

    # Finite controls of the written all-variation clock-degree argument.
    # The physical h has order EIGHT, not the proposed counterterm's 1024.
    jet_checks = 0
    t = s.Symbol("t")
    for order in (2, 4, 8):
        source = t**order * (1 + t + t**2)
        exchange = source**2 / (2 * (1 + t + 2 * t**2))
        leading = s.series(exchange, t, 0, 2 * order + 1).removeO().expand()
        assert leading == t ** (2 * order) / 2
        hessian_leading = (
            s.series(s.diff(exchange, t, 2), t, 0, 2 * order - 1).removeO().expand()
        )
        assert hessian_leading == order * (2 * order - 1) * t ** (2 * order - 2)
        jet_checks += 1

    # The actual S174 source instead has a nonzero quadratic clock jet.
    lapse, trace, hubble, h = s.symbols("lapse trace hubble h", nonzero=True)
    source = (
        ((1 - t * lapse) ** 2 - 1)
        / h
        * (3 * hubble + t * trace - 3 * hubble * (1 - t * lapse))
    )
    source2 = s.expand(source).coeff(t, 2)
    assert s.factor(source2 + 2 * lapse * (trace + 3 * hubble * lapse) / h) == 0
    assert source2.subs({hubble: 0, h: 1}) == -2 * lapse * trace
    # Flat Euclidean symbol control only; NOT the curved physical Green function.
    w, k, zeta = s.symbols("w k zeta", positive=True)
    momentum = s.Matrix([w, k])
    square = w**2 + k**2
    P = (1 + zeta * square) * s.eye(2) - zeta * momentum * momentum.T
    inverse = (s.eye(2) + zeta * momentum * momentum.T) / (1 + zeta * square)
    assert s.simplify(P * inverse - s.eye(2)) == s.zeros(2)
    residual = s.simplify((s.eye(2) - inverse)[0, 0])
    assert s.factor(residual - zeta * k**2 / (1 + zeta * square)) == 0
    assert residual.subs(k, 0) == 0 and residual.subs({k: 1, w: 1, zeta: 1}) == F(1, 3)
    assert s.diff(source2**2 / 2, lapse, 2) != 0

    # A regulator/evanescent conversion changes a finite reference; transport it.
    eps, pole, finite, evanescent = s.symbols("eps pole finite evanescent")
    raw = (pole / eps + finite) * (1 + eps * evanescent)
    assert s.expand(raw).coeff(eps, 0) == finite + pole * evanescent
    boundary = s.Symbol("boundary")
    assert (
        s.expand((finite + pole * evanescent) + (boundary - pole * evanescent))
        == finite + boundary
    )

    # Exact integration and all-energy majorant for a finite-angle Regge tail.
    q, u, alpha = s.symbols("q u alpha", positive=True)
    phi = 1 - 3 * q / 2 + q**3 / 2
    primitive = s.log(q) - 3 * q / 2 + q**3 / 6
    assert s.simplify(s.diff(primitive, q) - phi / q) == 0
    exponent = alpha * q**2
    primitive_energy = -(u ** (-exponent)) / exponent
    assert s.simplify(s.diff(primitive_energy, u) - u ** (-1 - exponent)) == 0
    assert s.limit(primitive_energy, u, s.oo) == 0
    assert sum(F(7, 10) ** j / factorial(j) for j in range(5)) > 2
    log_checks = 0
    for detector in (DETECTOR, F(1, 2), F(1, 3), F(3, 4), F(15, 16)):
        # Independent dyadic bracket plus the exp(7/10)>2 proof above.
        upper = log_inverse_upper(detector)
        power = int(upper / F(7, 10))
        assert 2 ** (power - 1) * detector < 1 <= 2**power * detector
        log_checks += 1
    energy = s.Symbol("energy", positive=True)
    assert s.cancel(s.diff(energy / (energy - 3), energy) + 3 / (energy - 3) ** 2) == 0
    # Difference of the J integrand and 1/s is positive for every s>=32.
    z = s.Symbol("z", nonnegative=True)
    numerator = s.fraction(s.factor(energy**2 / (energy - 3) ** 3 - 1 / energy))[0]
    assert all(c > 0 for c in s.Poly(numerator.subs(energy, 32 + z), z).all_coeffs())

    # Rodrigues proves all-spin positivity of the Regge comparison fixture.
    # These identities only cross-check its finite polynomial integration steps.
    angular = s.Symbol("angular", real=True)
    rodrigues_checks = 0
    for ell in range(9):
        rodrigues = s.diff((angular**2 - 1) ** ell, angular, ell) / (
            2**ell * factorial(ell)
        )
        assert s.expand(rodrigues - s.legendre(ell, angular)) == 0
        assert s.legendre(ell, -angular) == (-1) ** ell * s.legendre(ell, angular)
        rodrigues_checks += 1
    # The reflected even-spin exponential is no larger on the forward window.
    assert 32 - 4 - 2 > 0

    coefficient = regge_tail_bound(
        detector=DETECTOR, split=32, residue_over_slope=1, remainder_arc_error=0
    )
    assert coefficient < 60
    for T in (32, 64, 256, 1024):
        assert (
            regge_tail_bound(
                detector=DETECTOR, split=T, residue_over_slope=1, remainder_arc_error=0
            )
            <= coefficient
        )
    _, _, kappa, lam, _ = original_borns()
    assert 60 / kappa < lam / 10**198
    # Budget sensitivity only: gamma is NOT supplied physical Regge data.
    assert 60 * F(10**196) / kappa < lam / 100
    supplied = {
        "detector": DETECTOR,
        "split": F(64),
        "finite_forward_moment": F(3, 7),
        "residue_over_slope": F(2, 11),
        "regge_remainder_error": F(1, 13),
        "contour_error": F(1, 17),
        "observable_error": F(1, 19),
    }
    base_error = finite_window_error(**supplied)
    for name in ("regge_remainder_error", "contour_error", "observable_error"):
        assert finite_window_error(
            **{**supplied, name: supplied[name] + F(1, 23)}
        ) == base_error + F(1, 23)
    assert (
        finite_window_error(**{**supplied, "finite_forward_moment": F(6, 7)})
        > base_error
    )
    assert (
        finite_window_error(**{**supplied, "residue_over_slope": F(4, 11)}) > base_error
    )
    rejected = 0
    for name in supplied:
        incomplete = {key: value for key, value in supplied.items() if key != name}
        try:
            finite_window_error(**incomplete)
        except TypeError:
            rejected += 1
        else:
            raise AssertionError("Missing physical input accepted")
    bad = (
        {"detector": 0},
        {"detector": 1},
        {"detector": 0.5},
        {"split": 31},
        {"split": 32},
        {"split": True},
        {"finite_forward_moment": -1},
        {"residue_over_slope": -1},
        {"regge_remainder_error": -1},
        {"contour_error": -1},
        {"observable_error": -1},
    )
    for update in bad:
        try:
            finite_window_error(**{**supplied, **update})
        except (TypeError, ValueError):
            rejected += 1
        else:
            raise AssertionError("Invalid input accepted")
    assert inputs() == before
    return {
        "milestone": "P8.PRESCRIPTION-2.REFERENCE-REDUCTION + GRAVITY.FINITE-WINDOW",
        "outcome": "SCOPED_REFERENCE_REDUCTION_AND_FORWARD_MOMENT_FREE_TAIL_INTERFACE",
        "physical_model_adopted": False,
        "tracks_complete": False,
        "protected_input_files": len(before),
        "algebraic_connection_modes": 56,
        "retained_vector_components": 4,
        "gaussian_block_fixtures": gaussian_checks,
        "measure_scope": "ultralocal determinants vanish only in the stated perturbative dimensional measure; constrained physical determinant not supplied",
        "heavy_source_clock_jet_controls": jet_checks,
        "heavy_source_result": "tree source exchange begins at field degree16; source-dependent formal one-loop term begins at degree14; free heavy metric determinant remains",
        "proca_source_result": "nonzero degree2 clock source gives degree4 exchange and possible degree2 one-loop response; flat Euclidean symbol is a nonzero control, not a curved physical bound",
        "reference_conversion": "evanescent finite shifts require opposite boundary transport, not reset-to-zero",
        "dyadic_log_controls": log_checks,
        "rodrigues_controls": rodrigues_checks,
        "positive_spectral_fixture": "written all-even-spin Regge fixture has divergent full J but a finite finite-angle arc; not a full unitary/crossing UV model",
        "normalized_tail_coefficient_upper": str(coefficient),
        "conditional_tail_bound": "<60*C/alpha + supplied residual error at original detector; contour and observable errors remain separate",
        "conditional_weak_gravity_sensitivity": "IF C/alpha<=gamma/kappa and gamma<=1e196, leading Regge tail <lambda/100; no gamma bound established",
        "invalid_or_missing_inputs_rejected": rejected,
        "remaining": "dimensional/gauge/constraint/state reference and Proca-current response; actual finite-window moment, Regge envelope and all physical errors; M/V/G/B/R/P8 remain open",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
