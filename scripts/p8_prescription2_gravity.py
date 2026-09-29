"""PRESCRIPTION-2 initial packet and a conditional massive gravity-arc kernel.

No quantum model is adopted and no physical spectral moment is filled in.
The all-spin comparison uses a written ball-autocorrelation proof. Exact
finite checks verify its algebra, not the omitted physical hypotheses.
"""

import hashlib
import json
from fractions import Fraction as F
from math import comb

import sympy as s
from p8_match1_audit import POINTS, solve_fraction
from p8_match1_rate_input import REPO, exact
from p8_prescription1 import inputs as previous_inputs
from p8_rate4_candidate import dot, original_borns, rational
from p8_rate4_d8 import matrix, rate4_d8_coefficients

EXTRA_SOURCES = {
    "scripts/p8_prescription1.py": "b5b53643fa9afd4c803d5d1942457a87080eb0810f09a2872a39ebf64b04e2f6",
    "problems/P8/s6/FORMULATION.md": "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901",
    "docs/assessment-2026-09-21-p8-rate4-d8-candidate.md": "7dc612d7b3def5b97232c96d798292b64d62f0ef71b56eded20b4584748d8d1a",
}
DETECTOR = F(1, 256)


def inputs():
    manifest = previous_inputs()
    for name, expected in EXTRA_SOURCES.items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected
        if name in manifest:
            assert manifest[name] == expected
        manifest[name] = expected
    return manifest


def radial_moment(power, *, lower):
    """Integral_lower^1 q^(2*power+1) phi(q) dq, phi=1-3q/2+q^3/2."""
    lower = exact(lower)
    if type(power) is not int or power < 0 or not 0 <= lower < 1:
        raise ValueError("Require nonnegative integer power and 0<=lower<1")
    return sum(
        c * (1 - lower ** (2 * power + j + 2)) / (2 * power + j + 2)
        for j, c in ((0, F(1)), (1, F(-3, 2)), (3, F(1, 2)))
    )


def overlap_partial_wave(ell, *, energy_squared, lower):
    """Exact comparison-kernel integral, NOT the full gravity spectral kernel."""
    energy = exact(energy_squared)
    if type(ell) is not int or ell < 0 or energy < 5:
        raise ValueError("Require nonnegative integer spin and s>=5")
    return sum(
        F((-1) ** k * comb(ell, k) * comb(ell + k, k))
        * radial_moment(k, lower=lower)
        / (energy - 4) ** k
        for k in range(ell + 1)
    )


def arc_local_vector(*, detector):
    """Normalized [v^2] at s=2+v,t=-q^2,u=2+q^2-v, not the forward L."""
    normalizer = radial_moment(0, lower=detector)
    m2 = radial_moment(1, lower=detector) / normalizer
    m4 = radial_moment(2, lower=detector) / normalizer
    return F(0), F(2), m2, 32 + 16 * m2 + 12 * m4


def arc_weights(*, detector):
    return tuple(
        solve_fraction(
            list(zip(*matrix(), strict=True)), arc_local_vector(detector=detector)
        )
    )


def kernel_error_bound(*, detector, uv_s_min, positive_forward_moment, other_error):
    """Conditional normalized arc error; no spectral/IR/contour input inferred.

    J=(2/pi) int_S^infinity rho(s,1)/(s-3)^3 ds must be finite and supplied
    for a positive partial-wave reference. other_error covers the COMPLETE
    normalized difference from the physical observable (unitarity mismatch,
    contour boundary and other unaccounted terms), without double counting.
    """
    detector, threshold, moment, error = map(
        exact, (detector, uv_s_min, positive_forward_moment, other_error)
    )
    if not 0 < detector < 1 or threshold < 32 or moment < 0 or error < 0:
        raise ValueError("Require 0<E<1, S>=32 and explicit nonnegative errors")
    normalizer = radial_moment(0, lower=detector)
    angular = detector**2 / 2
    crossing = 3 * radial_moment(1, lower=0) / (2 * (threshold - 3))
    return (angular + crossing) * moment / normalizer + error


def scalar_forward_moment(*, masses_squared, residues):
    """J for supplied leading scalar exchanges only; not the full gravity J."""
    masses, weights = tuple(map(exact, masses_squared)), tuple(map(exact, residues))
    if not masses or len(masses) != len(weights):
        raise ValueError("Require equally sized nonempty spectral data")
    if min(masses) < 32 or min(weights) < 0:
        raise ValueError("Require z>=32 and nonnegative residues")
    return 2 * sum(r / (z - 3) ** 3 for z, r in zip(masses, weights, strict=True))


def source_coefficients(*, masses_squared, scalar_sources, curvature_sources):
    """Leading local J (z+Box)^-1 J /2 coefficients, not a full parent."""
    z, h, a = (
        tuple(map(exact, data))
        for data in (masses_squared, scalar_sources, curvature_sources)
    )
    if not z or len(z) != len(h) or len(z) != len(a) or min(z) <= 0:
        raise ValueError("Require equally sized nonempty data with z>0")
    return (
        sum(x * x / (8 * mass) for mass, x in zip(z, h, strict=True)),
        sum(x * y / (2 * mass) for mass, x, y in zip(z, h, a, strict=True)),
        sum(y * y / (2 * mass) for mass, y in zip(z, a, strict=True)),
    )


def audit():
    before = inputs()
    q, v, x, r = s.symbols("q v x r", real=True)
    phi = 1 - 3 * q / 2 + q**3 / 2
    half = s.Rational(1, 2)
    cap_height = (1 - q) / 2
    overlap = 2 * s.pi * cap_height**2 * (half - cap_height / 3)
    assert s.expand(overlap / (4 * s.pi * half**3 / 3) - phi) == 0
    assert s.factor(phi) == (q - 1) ** 2 * (q + 2) / 2

    moments_checked = 0
    for lower in (F(0), DETECTOR, F(1, 2)):
        for power in range(7):
            literal = s.integrate(q ** (2 * power + 1) * phi, (q, rational(lower), 1))
            assert literal == rational(radial_moment(power, lower=lower))
            moments_checked += 1
    assert radial_moment(0, lower=0) == F(1, 10)
    assert radial_moment(1, lower=0) == F(3, 140)
    W = radial_moment(0, lower=DETECTOR)
    assert F(1, 11) < W < F(1, 10)

    wave_checks = 0
    for energy in (F(8), F(32), F(64)):
        for ell in range(9):
            for lower in (F(0), DETECTOR):
                literal = s.integrate(
                    s.expand(q * phi * s.legendre(ell, 1 - 2 * q**2 / (energy - 4))),
                    (q, rational(lower), 1),
                )
                assert literal == rational(
                    overlap_partial_wave(ell, energy_squared=energy, lower=lower)
                )
                wave_checks += 1
    negatives = {}
    for ell in (128, 256, 512):
        full = overlap_partial_wave(ell, energy_squared=32, lower=0)
        cut = overlap_partial_wave(ell, energy_squared=32, lower=DETECTOR)
        assert full > 0 and abs(cut - full) <= DETECTOR**2 / 2
        if ell >= 256:
            assert cut < 0
            negatives[str(ell)] = "strictly negative, exact Fraction arithmetic"
        else:
            assert cut > 0
    # Abel-Poisson identity used in the written finite-gap impossibility proof.
    generating = 1 / s.sqrt(1 - 2 * r * x + r * r)
    poisson = (1 - r * r) / (1 - 2 * r * x + r * r) ** s.Rational(3, 2)
    assert s.simplify(2 * r * s.diff(generating, r) + generating - poisson) == 0

    # Exact ingredients of the written crossing and retained-pole estimates.
    energy, y = s.symbols("energy y", positive=True)
    assert s.diff((energy - 2 - y) ** -3, y) == 3 / (energy - 2 - y) ** 4
    t_channel = (2 + 2 * q**2 - (2 + v) * (2 + q**2 - v)) / (-(q**2))
    assert s.expand(-t_channel).coeff(v, 2) == 1 / q**2
    pole_primitive = s.log(q) - 3 * q / 2 + q**3 / 6
    assert s.simplify(s.diff(pole_primitive, q) - phi / q) == 0

    S = (2 + v) ** 2 + q**4 + (2 + q**2 - v) ** 2
    U = -(q**2) * (2 + v) * (2 + q**2 - v)
    local = arc_local_vector(detector=DETECTOR)
    weights = arc_weights(detector=DETECTOR)
    for actual, word in zip(local, (s.Integer(1), S, U, S**2), strict=True):
        coeff = s.expand(word).coeff(v, 2)
        assert s.integrate(
            s.expand(q * phi * coeff), (q, rational(DETECTOR), 1)
        ) / rational(W) == rational(actual)
    symbolic_weights = s.Matrix([list(map(rational, local))]) * s.Matrix(matrix()).inv()
    assert list(symbolic_weights) == list(map(rational, weights))
    assert sum(map(abs, weights)) < 1
    tail_checks = []
    for i in range(6):
        for j in range(4):
            degree = 2 * i + 3 * j
            if not 5 <= degree <= 11:
                continue
            coefficient = s.expand(S**i * U**j).coeff(v, 2)
            projected = s.integrate(
                s.expand(q * phi * coefficient), (q, rational(DETECTOR), 1)
            ) / rational(W)
            samples = [
                sum(a * a for a in point) ** i * s.prod(point) ** j for point in POINTS
            ]
            difference = F(projected) - dot(weights, samples)
            tail_checks.append((abs(difference) / 64**degree, i, j))
    maximum, i_max, j_max = max(tail_checks)
    assert (i_max, j_max) == (3, 0)
    # Written Cauchy + sample bound for EVERY word of weighted degree>=12.
    assert F(328, 64**2) < F(2, 7) ** 2
    assert F(576, 64**3) < F(2, 7) ** 3
    assert F(26, 64**2) < F(2, 7) ** 2
    assert F(12, 64**3) < F(2, 7) ** 3
    infinite_bound = (1 + sum(map(abs, weights))) * F(2, 7) ** 12
    assert infinite_bound < maximum < F(1, 800000)

    n, g, kappa, lam, born = original_borns()
    moment = scalar_forward_moment(masses_squared=(n,), residues=(g * g,))
    assert moment < 5 * lam
    unit_bound = kernel_error_bound(
        detector=DETECTOR, uv_s_min=32, positive_forward_moment=1, other_error=0
    )
    assert unit_bound < F(1, 80)
    heavy_error = kernel_error_bound(
        detector=DETECTOR, uv_s_min=32, positive_forward_moment=moment, other_error=0
    )
    assert heavy_error < lam / 16
    assert 2000 * lam * maximum < lam / 400
    assert 66 / kappa < lam / 10**198  # Written log(256)<6 and W>1/11.

    # Existing four physical targets still fit the new arc, without a fifth rate.
    hard, conversion = (
        (F(1, 7), F(-2, 9), F(3, 11), F(-5, 13)),
        (F(1, 101), F(-1, 103), F(2, 107), F(-3, 109)),
    )
    c = rate4_d8_coefficients(
        born=born, other_hard_real=hard, rate_conversion=conversion
    )
    residuals = tuple(
        h + a * d / 2 for h, a, d in zip(hard, born, conversion, strict=True)
    )
    assert dot(local, c) == -dot(weights, residuals)

    # A small sixth moment alone does NOT control the gravity third moment.
    moment_controls = 0
    for z in (32, 64, 128, 256):
        residue = lam * z**6
        assert residue / z**6 == lam
        assert (
            scalar_forward_moment(masses_squared=(z,), residues=(residue,))
            > 2 * lam * z**3
        )
        moment_controls += 1

    # Explicit local scalar/curvature source correlations and clock-localized lift.
    source_controls = 0
    for masses, h, a in (
        ((32,), (2,), (3,)),
        ((32, 64), (2, 1), (3, -1)),
        ((n,), (g,), (0,)),
    ):
        flat, mixed, curved = source_coefficients(
            masses_squared=masses, scalar_sources=h, curvature_sources=a
        )
        assert mixed * mixed <= 4 * flat * curved
        if len(masses) == 1:
            assert mixed * mixed == 4 * flat * curved
        source_controls += 1
    H, h, a, z, R, field = s.symbols("H h a z R field")
    J = h * field**2 / 2 + a * R
    assert s.expand((-z * H**2 / 2 + H * J).subs(H, J / z) - J**2 / (2 * z)) == 0
    jet_controls = 0
    t = s.Symbol("t")
    for order in (4, 8, 16):
        X = 1 + 2 * t + 3 * t * t
        switch = (1 - X) ** order / (X**order + (1 - X) ** order)
        truncated = (
            s.series(switch * (1 + t + t**3), t, 0, order + 1).removeO().expand()
        )
        assert truncated == 2**order * t**order
        jet_controls += 1
    assert (x**1024 + (1 - x) ** 1024).subs(x, 1) == 1

    rejected = 0
    bad_calls = (
        lambda: radial_moment(True, lower=0),
        lambda: radial_moment(-1, lower=0),
        lambda: radial_moment(0, lower=1),
        lambda: overlap_partial_wave(-1, energy_squared=32, lower=0),
        lambda: overlap_partial_wave(0, energy_squared=4, lower=0),
        lambda: kernel_error_bound(
            detector=0, uv_s_min=32, positive_forward_moment=1, other_error=0
        ),
        lambda: kernel_error_bound(
            detector=DETECTOR, uv_s_min=31, positive_forward_moment=1, other_error=0
        ),
        lambda: kernel_error_bound(
            detector=DETECTOR, uv_s_min=32, positive_forward_moment=-1, other_error=0
        ),
        lambda: kernel_error_bound(
            detector=DETECTOR, uv_s_min=32, positive_forward_moment=1, other_error=-1
        ),
        lambda: kernel_error_bound(
            detector=DETECTOR, uv_s_min=32, positive_forward_moment=1
        ),
        lambda: scalar_forward_moment(masses_squared=(), residues=()),
        lambda: scalar_forward_moment(masses_squared=(32,), residues=(-1,)),
        lambda: scalar_forward_moment(masses_squared=(16,), residues=(1,)),
        lambda: source_coefficients(
            masses_squared=(0,), scalar_sources=(1,), curvature_sources=(1,)
        ),
        lambda: source_coefficients(
            masses_squared=(32,), scalar_sources=(1, 2), curvature_sources=(1,)
        ),
        lambda: arc_local_vector(detector=0.5),
    )
    for call in bad_calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
        else:
            raise AssertionError("Incomplete or invalid input accepted")
    assert inputs() == before
    return {
        "milestone": "P8.PRESCRIPTION-2 + GRAVITY.APPLICABILITY, initial packet",
        "outcome": "BOUNDARY_PROPOSAL_AND_CONDITIONAL_ALL_SPIN_COMPARISON_NOT_PHYSICAL_MATCHING",
        "physical_model_adopted": False,
        "protected_input_files": len(before),
        "independent_radial_moments": moments_checked,
        "independent_low_spin_integrals": wave_checks,
        "finite_detector_exact_negative_controls": negatives,
        "all_spin_comparison": "written positive ball-autocorrelation kernel proof, every s>=32 and every spin",
        "finite_gap_obstruction": "written Abel-Poisson proof; not a no-go for approximate finite-resolution bounds",
        "arc_weights": list(map(str, weights)),
        "arc_weight_norm": str(sum(map(abs, weights))),
        "tail_words_checked": len(tail_checks),
        "tail_operator_norm": str(maximum),
        "tail_maximizing_word": "S^3",
        "D8_tail_hypothesis_consequence": "<lambda/400 for N64<=2000lambda, NOT a proof of the hypothesis",
        "conditional_kernel_error_per_forward_moment": str(unit_bound),
        "kernel_error_upper_bound": "<=J/80 plus separately supplied complete physical errors; strict for J>0",
        "original_heavy_tree_forward_moment": "<5lambda; selected leading scalar exchange only",
        "original_heavy_tree_angular_crossing_error": "<lambda/16; not the full gravitational allowance",
        "explicit_t_channel_pole_arc": "positive and <66/kappa <1e-198lambda; kept, not discarded",
        "sixth_vs_third_moment_controls": moment_controls,
        "correlated_source_controls": source_controls,
        "clock_switch_jet_controls": jet_controls,
        "clock_jet_statement": "proposed finite local lift has zero reference variations through degree1023; formal first-order counterterm effect only",
        "invalid_inputs_rejected": rejected,
        "remaining": "complete renormalized covariant reference, physical complementary boundary review, full positive forward moment or sharper spectral estimate, finite-coupling/contour/low-cut errors, interacting state and global bounce; M/V/G/B/R/P8 remain open",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
