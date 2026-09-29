"""Read-only PRESCRIPTION-1 feasibility checks, not a physical parent certificate.

Tests explicit cutoff ansatz obstructions, a selected one-loop UV-shell bound,
a restricted positive-spectral matching mechanism and the unchanged S240
heavy-state high-momentum readout. No new action, cutoff or state is adopted.
"""

import hashlib
import json
from fractions import Fraction as F
from math import factorial

import sympy as s
from p8_match1_audit import POINTS, solve_fraction
from p8_match1_rate_input import REPO, exact
from p8_rate4_candidate import dot, original_borns, rational
from p8_rate4_d8 import matrix
from p8_rate4_d8_witness import inputs as witness_inputs
from p8_rate4_parent import local_case, positive_after

EXTRA_SOURCES = {
    "scripts/p8_rate4_d8_witness.py": "4e917ab534680acbf920e1a31e4edde2180ff9e9014712c76c4696d239f53c37",
    "problems/P8/s6/continuation/s6_235/notes/diagrams.md": "b0567c3784a8a90a24b7021aeb4cda4a915bc2831af1d4d422565dd8ec162694",
    "problems/P8/s6/continuation/s6_263/FORMULATION.md": "deee3c325ff63b0fde44ca2db94d960524160ba472d16b1748a31174dbd0dce7",
    "problems/P8/s6/continuation/s6_264/FORMULATION.md": "38e1cbf81e1b24a0a3d4f25b933aeb365f8476a05ce78d075441b0f7c704572f",
    "problems/P8/s6/continuation/s6_265/notes/integrability.md": "e9d7468d8e5f5877b8c937e49495213f5dfe33f0514ade16dcbc05a8ac0af4c3",
    "problems/P8/s6/continuation/s6_278/notes/dispersion.md": "95e4a6fdabb87aa0a72aa287e12dfa57b6b8912cbbe0fefa1e17444e292ec609",
    "problems/P8/s6/continuation/s6_277/notes/scope.md": "d2efbfeab21d6c44fe97fd649011a3c6f40af312779fc7b83ae0dc5e10c98966",
}
WEIGHTS = local_case(4)[1]
AMPLIFICATION = 1 + sum(map(abs, WEIGHTS))


def inputs():
    manifest = witness_inputs()
    for name, expected in EXTRA_SOURCES.items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected
        if name in manifest:
            assert manifest[name] == expected
        manifest[name] = expected
    assert "problems/P8/s6/continuation/s6_240/notes/bounds.md" in manifest
    assert "problems/P8/s6/continuation/s6_240/notes/state.md" in manifest
    return manifest


def channel_power_projection(degree):
    """L P4 sum(channel**degree); mass-one conventions, not a positive UV arc."""
    if type(degree) is not int or degree < 0:
        raise ValueError("Require a nonnegative integer channel degree")
    forward = F(degree * (degree - 1) * 2 ** (degree - 2)) if degree >= 2 else F(0)
    return forward - dot(WEIGHTS, [sum(a**degree for a in q) for q in POINTS])


def shell_bound(*, cutoff, prefactor_abs):
    """Entire short-proper-time bubble contribution AFTER the four-rate fit.

    The supplied prefactor is an actual sector coefficient, not a bound on
    all sectors. e**(4/cutoff**2)<3 for cutoff>=2. The mass can be any
    nonnegative mass squared; retaining e**(-mass_squared*tau) improves it.
    """
    cutoff, prefactor = map(exact, (cutoff, prefactor_abs))
    if cutoff < 2 or prefactor < 0:
        raise ValueError("Require cutoff>=2 and a nonnegative sector prefactor")
    return prefactor * 9 * AMPLIFICATION * 4**5 / (600 * cutoff**10)


def shell_coefficient(*, tau_power, mass_squared):
    """Exact coefficient of cutoff**(-2*tau_power) in the matched UV shell."""
    if type(tau_power) is not int or tau_power < 1:
        raise ValueError("Require a positive integer tau power")
    mass = exact(mass_squared)
    if mass < 0:
        raise ValueError("Require nonnegative mass squared")
    return sum(
        channel_power_projection(k)
        * F(factorial(k), factorial(2 * k + 1))
        * (-mass) ** (tau_power - k)
        / (factorial(tau_power - k) * tau_power)
        for k in range(1, tau_power + 1)
    )


def shell_series_error(*, cutoff, mass_squared, through):
    """Written exponential remainder bound for a finite exact shell series."""
    cutoff, mass = map(exact, (cutoff, mass_squared))
    if type(through) is not int or through < 0:
        raise ValueError("Require a nonnegative integer truncation order")
    if mass < 0 or cutoff <= 0 or cutoff**2 < mass + 4:
        raise ValueError("Require cutoff**2>=mass_squared+4 and cutoff>0")
    z = (mass + 4) / cutoff**2
    return (
        18
        * AMPLIFICATION
        * z ** (through + 1)
        / ((through + 1) * factorial(through + 1))
    )


def spectral_budget(*, masses_squared, residues):
    """One sufficient moment in the RESTRICTED positive scalar-exchange class.

    Not a representation theorem for generic UV/gravity matching. Every
    residue is supplied. Contact, other-spin and curved data are excluded.
    """
    masses = tuple(map(exact, masses_squared))
    weights = tuple(map(exact, residues))
    if not masses or len(masses) != len(weights):
        raise ValueError("Require equally sized nonempty spectral data")
    if min(masses) < 32 or min(weights) < 0:
        raise ValueError("Require masses squared>=32 and nonnegative residues")
    moment = sum(r / mass**6 for mass, r in zip(masses, weights, strict=True))
    return 6 * AMPLIFICATION * 16**5 * moment


def heavy_state_tail(*, comoving_split, kappa, mode_majorant):
    """S240 renormalized nonlocal readout tail on [-1,1], through five jets.

    Fixed comoving split, no state truncation/reset. Local heat-counteraction
    stays intact. This is not an interacting response or light/metric bound.
    """
    split, coupling, majorant = map(exact, (comoving_split, kappa, mode_majorant))
    if split <= 0 or coupling <= 0 or majorant < 0:
        raise ValueError("Require positive split/kappa and nonnegative majorant")
    return 29 * majorant / (coupling * split**2)


def audit():
    before = inputs()
    n0, g, kappa0, lam, _ = original_borns()

    # Independent symbolic and Fraction low-degree fits and projections.
    independent = solve_fraction(list(zip(*matrix(), strict=True)), [0, 2, 0, 32])
    assert tuple(independent) == WEIGHTS
    v, x, a = s.symbols("v x a", real=True)
    project_checks = 0
    for degree in range(13):
        forward = s.expand((2 + v) ** degree + (2 - v) ** degree).coeff(v, 2)
        symbolic = forward - sum(
            rational(w) * sum(s.Integer(c) ** degree for c in q)
            for w, q in zip(WEIGHTS, POINTS, strict=True)
        )
        assert symbolic == rational(channel_power_projection(degree))
        project_checks += 1
    assert all(channel_power_projection(k) == 0 for k in range(5))
    assert channel_power_projection(5) == -F(2255640, 3259)

    # Direct parameter-integral calculation of sixteen exact series coefficients.
    series_checks = 0
    for mass in (1, 3):
        for order in range(1, 9):
            integrand = s.expand((a * x * (1 - x) - mass) ** order - (-mass) ** order)
            polynomial = s.Poly(s.integrate(integrand, (x, 0, 1)), a)
            literal = sum(
                coefficient * rational(channel_power_projection(power[0]))
                for power, coefficient in polynomial.terms()
            ) / (factorial(order) * order)
            assert literal == rational(
                shell_coefficient(tau_power=order, mass_squared=mass)
            )
            series_checks += 1
    leading = shell_coefficient(tau_power=5, mass_squared=1)
    assert leading == -F(18797, 45169740)
    assert 9 * AMPLIFICATION * 4**5 / 600 == F(2612352, 81475) < 33
    shell_controls = 0
    for cutoff in (3, 10, 1000, 10**100):
        partial = sum(
            shell_coefficient(tau_power=j, mass_squared=1) / cutoff ** (2 * j)
            for j in range(1, 13)
        )
        error = shell_series_error(cutoff=cutoff, mass_squared=1, through=12)
        assert abs(partial) + error < shell_bound(cutoff=cutoff, prefactor_abs=1)
        shell_controls += 1
    # All cutoff>=1000: the proven order-six error is smaller than the leading term.
    error_ratio = shell_series_error(cutoff=1000, mass_squared=1, through=5) / (
        abs(leading) / 1000**10
    )
    assert error_ratio == F(1571493, 4812032) < 1
    assert shell_bound(cutoff=10**100, prefactor_abs=1) < lam / 10**398
    C = -g * g * (3 / (n0 - 2) - 2 / (n0 - 2) ** 2)
    assert 0 < C * C / 288 < 1  # C^2/(32*pi^2) is the selected S235 bubble prefactor.

    # RG/scheme compensation does not fix the physical integration constant.
    boundary_controls = 0
    for physical_projection in (-8 * lam, F(0), 8 * lam):
        for scale in (2, 10):
            calculated_piece = leading / scale**10
            supplied_boundary = physical_projection - calculated_piece
            assert calculated_piece + supplied_boundary == physical_projection
            boundary_controls += 1

    # Exact finite-cutoff local-field tests, not tests of the original P8 parent.
    p2, mass, cutoff2 = s.symbols("p2 mass cutoff2", positive=True)
    regulated = s.exp(-(p2 + mass) / cutoff2) / (p2 + mass)
    derivative = s.simplify(s.diff(p2 * regulated, p2) * s.exp((p2 + mass) / cutoff2))
    assert (
        s.simplify(derivative.subs(p2, cutoff2) + cutoff2 / (cutoff2 + mass) ** 2) == 0
    )
    assert s.simplify(s.diff(p2 / (p2 + mass), p2) - mass / (p2 + mass) ** 2) == 0
    t, c0, c2, c4 = s.symbols("t c0 c2 c4")
    time_kernel = lambda z: c0 + c2 * z**2 / 2 + c4 * z**4 / 24
    determinant = s.expand(
        time_kernel(2 * t) * time_kernel(4 * t) - time_kernel(3 * t) ** 2
    )
    assert determinant.coeff(t, 2) == c0 * c2
    # For the proper-time kernel C(0)>0 and C''(0)<0 by their positive integrals.
    momentum, radius, cutoff = s.symbols("momentum radius cutoff", positive=True)
    band_delta = (s.sin(cutoff * radius) - cutoff * radius * s.cos(cutoff * radius)) / (
        2 * s.pi**2 * radius**3
    )
    literal_delta = s.integrate(
        momentum * s.sin(momentum * radius), (momentum, 0, cutoff)
    ) / (2 * s.pi**2 * radius)
    assert s.simplify(literal_delta - band_delta) == 0
    assert s.simplify(band_delta.subs(radius, s.pi / cutoff)) == cutoff**3 / (
        2 * s.pi**4
    )

    # Positive scalar exchanges supply a moment, but P4 itself is not a positive arc.
    n = s.Symbol("n", positive=True)
    KQ = s.factor(
        2 / (n - 2) ** 3
        - sum(
            rational(w) * sum(1 / (n - channel) for channel in q)
            for w, q in zip(WEIGHTS, POINTS, strict=True)
        )
    )
    positive_after(-KQ, n, 32)
    spectral_checks = 0
    for masses, residues in (
        ((F(32),), (F(1),)),
        ((F(32), F(64)), (lam, 2 * lam)),
        ((F(64), F(1024), n0), (F(1, 10), F(3, 10), g * g)),
    ):
        exact_projection = sum(
            rational(r) * KQ.subs(n, rational(m))
            for m, r in zip(masses, residues, strict=True)
        )
        bound = spectral_budget(masses_squared=masses, residues=residues)
        assert exact_projection < 0 and abs(exact_projection) < rational(bound)
        spectral_checks += 1
    assert spectral_budget(masses_squared=(32, 64), residues=(lam, 2 * lam)) < lam / 4
    phi, heavy, h, z, positive_quartic, mu = s.symbols(
        "phi heavy h z positive_quartic mu"
    )
    potential = (
        mu * phi**2 / 2
        + z * heavy**2 / 2
        - h * heavy * phi**2 / 2
        + (3 * h**2 / z + positive_quartic) * phi**4 / 24
    )
    completed = (
        mu * phi**2 / 2
        + z * (heavy - h * phi**2 / (2 * z)) ** 2 / 2
        + positive_quartic * phi**4 / 24
    )
    assert s.expand(potential - completed) == 0

    # A real existing-state readout bound: retain every local heat contribution.
    p, n = s.symbols("p n", positive=True)
    primitive = 64 * p**3 / (3 * n * (p**2 + 16 * n) ** s.Rational(3, 2))
    assert (
        s.simplify(s.diff(primitive, p) - p**2 / (n + p**2 / 16) ** s.Rational(5, 2))
        == 0
    )
    assert F(256, 9) < 29
    tail = heavy_state_tail(comoving_split=10**100, kappa=kappa0, mode_majorant=10**301)
    assert tail < F(1, 10**697)
    tail_controls = 0
    for split in (1, 32, 10**64, 10**100):
        assert heavy_state_tail(
            comoving_split=2 * split, kappa=kappa0, mode_majorant=10**301
        ) * 4 == heavy_state_tail(
            comoving_split=split, kappa=kappa0, mode_majorant=10**301
        )
        tail_controls += 1
    assert n0 < (10**100) ** 2 and 2 * 10**64 < 10**100

    rejected = 0
    bad_calls = (
        lambda: shell_bound(cutoff=1, prefactor_abs=1),
        lambda: shell_bound(cutoff=2, prefactor_abs=-1),
        lambda: shell_bound(cutoff=2.0, prefactor_abs=1),
        lambda: shell_bound(cutoff=2),
        lambda: shell_coefficient(tau_power=True, mass_squared=1),
        lambda: shell_coefficient(tau_power=0, mass_squared=1),
        lambda: shell_coefficient(tau_power=5, mass_squared=-1),
        lambda: shell_series_error(cutoff=2, mass_squared=1, through=5),
        lambda: shell_series_error(cutoff=3, mass_squared=1, through=-1),
        lambda: spectral_budget(masses_squared=(), residues=()),
        lambda: spectral_budget(masses_squared=(32,), residues=(1, 2)),
        lambda: spectral_budget(masses_squared=(16,), residues=(1,)),
        lambda: spectral_budget(masses_squared=(32,), residues=(-1,)),
        lambda: heavy_state_tail(comoving_split=0, kappa=1, mode_majorant=1),
        lambda: heavy_state_tail(comoving_split=1, kappa=0, mode_majorant=1),
        lambda: heavy_state_tail(comoving_split=1, kappa=1, mode_majorant=-1),
        lambda: channel_power_projection(-1),
    )
    for call in bad_calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
        else:
            raise AssertionError("Invalid feasibility input accepted")
    assert inputs() == before
    return {
        "milestone": "P8.PRESCRIPTION-1",
        "outcome": "FEASIBILITY_SCREEN_COMPLETE_NO_COMMON_QUANTUM_PARENT",
        "physical_parent_adopted": False,
        "protected_input_files": len(before),
        "proper_time_physical_local_covariance": "EXCLUDED_ANSATZ; not an exclusion of a computational regulator or original P8",
        "spatial_cutoff_original_local_field": "NOT_LOCAL; coarse-grained observables remain possible",
        "cutoff_obstruction_identities": 5,
        "channel_power_projection_checks": project_checks,
        "independent_shell_series_coefficients": series_checks,
        "shell_series_enclosure_controls": shell_controls,
        "selected_shell_bound": "<=33*abs(prefactor)/Lambda^10 for Lambda>=2, strict for nonzero prefactor; not the complete matching error",
        "selected_shell_leading_coefficient": str(leading),
        "selected_shell_sign": "negative for positive prefactor, mass_squared=1, Lambda>=1000",
        "selected_shell_at_test_scale_relative_to_lambda": "<1e-398 for Lambda=1e100, prefactor<=1; no physical cutoff adopted",
        "RG_boundary_compensation_controls": boundary_controls,
        "positive_spectral_moment_constant": str(6 * AMPLIFICATION * 16**5),
        "spectral_projection_controls": spectral_checks,
        "spectral_halfline_sign_proofs": 1,
        "stable_scalar_sector_square_completion": True,
        "heavy_state_tail_antiderivative_identity": True,
        "heavy_state_tail_scaling_controls": tail_controls,
        "heavy_reference_C5_normalized_tail": "<1e-697 above fixed comoving P=1e100, on [-1,1]; local terms retained",
        "invalid_inputs_rejected": rejected,
        "scope": "Written cutoff-ansatz obstructions and scoped bound mechanisms only. No all-sector matching, physical cutoff, all-spin gravity functional, interacting-state/bounce theorem, D8 witness or P8 closure.",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
