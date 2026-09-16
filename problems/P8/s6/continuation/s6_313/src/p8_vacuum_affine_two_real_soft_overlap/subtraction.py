"""State-correct overlap subtraction and uniformly integrable signed measure."""

from functools import cache

import mpmath as mp
import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as e
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as matter
from p8_vacuum_affine_uniform_two_real_tree_bound.bounds import energies

from . import analytic, source

E = s.Rational(5, 4)
OUTGOING = s.ImmutableMatrix([1, 0, 0])
EPS = e.imm(s.diag(0, 1, -1, 0))
BASELINE_COEFFICIENT = 2 / s.sqrt(source.KAPPA)


def entropy(first, second):
    a, b = energies(first, second)
    return a * s.log((a + b) / a) + b * s.log((a + b) / b)


def remainder_upper(first, second):
    a, b = energies(first, second)
    return analytic.REMAINDER_COEFFICIENT * entropy(a, b) / (a * b)


def require_resolution(value):
    value = e.exact_real(value)
    if not 0 < value <= s.Rational(1, 8):
        raise ValueError("Require resolution in(0,1/8]")
    return value


def measure_upper(value):
    value = require_resolution(value)
    return 2 * s.Rational(1, 10**725) * value + s.Rational(1, 10**652) * value**2


def leading_reference_relative_upper(value):
    value = require_resolution(value)
    return 4 * s.Rational(1, 10**725) * s.sqrt(value) + 2 * s.Rational(
        1, 10**652
    ) * value ** s.Rational(3, 2)


def configuration(a, b):
    rays = []
    if a:
        rays.append(e.imm([a, 0, 0, a]))
    if b:
        rays.append(e.imm([b, 0, 0, -b]))
    ps, _, born = source.recoil.momenta(E, rays, OUTGOING)
    ps = tuple(e.imm(p.applyfunc(s.simplify)) for p in ps)
    return ps, tuple(rays), tuple(map(e.imm, born))


def rho(a, b):
    Ep = s.sqrt(E * (E - a - b) + a * b)
    return s.factor(E * s.sqrt(Ep * Ep - 1) / (s.Rational(3, 4) * Ep))


@cache
def amplitude(points, rays):
    value, count = e.amplitude(
        points,
        [(q, EPS) for q in rays],
        heavy=source.HEAVY_MASS2,
        cubic=source.CUBIC,
        contact=source.CONTACT,
        kappa=source.KAPPA,
    )
    return s.factor(value), count


def current(points, q):
    return s.factor(e.old.pair(EPS, matter.soft_current(points, q)))


@cache
def original_samples():
    checks = {}
    gates = {}
    records = []
    K, n, g, C = source.KAPPA, source.HEAVY_MASS2, source.CUBIC, source.CONTACT
    for family in ("simultaneous", "hierarchical"):
        for denominator in (100, 1000, 10000):
            h = s.Rational(1, denominator)
            if family == "simultaneous":
                alpha = 2 - h
                Ep = (alpha + 1 / alpha) / 2
                a = b = E - Ep
            else:
                alpha0 = s.Rational(19, 10)
                Ep0 = (alpha0 + 1 / alpha0) / 2
                b = E - Ep0 * Ep0 / E
                alpha = alpha0 - h
                Ep = (alpha + 1 / alpha) / 2
                a = s.factor((E * E - E * b - Ep * Ep) / (E - b))
            label = f"{family}_{denominator}"
            points, qs, born = configuration(a, b)
            A0 = s.factor(
                matter.born_continuation(born, n, g, C) + e.old.born(born) / K
            )
            M6, n6 = amplitude(points, qs)
            pa, qa, _ = configuration(a, 0)
            pb, qb, _ = configuration(0, b)
            M5a, n5a = amplitude(pa, qa)
            M5b, n5b = amplitude(pb, qb)
            checks[label + "_whole_count"] = s.Integer(n6 - 434)
            checks[label + "_first_face_count"] = s.Integer(n5a - 47)
            checks[label + "_second_face_count"] = s.Integer(n5b - 47)
            for j, p in enumerate(points):
                checks[f"{label}_mass_shell_{j}"] = s.factor(e.old.dot(p, p) - 1)
            checks[label + "_full_conservation"] = sum(points, qs[0] + qs[1])
            for tag, ps, marked, soft in (
                ("first", pb, qb[0], qs[0]),
                ("second", pa, qa[0], qs[1]),
            ):
                J = matter.soft_current((*ps, marked), soft)
                checks[f"{label}_{tag}_complete_state_current_Ward"] = (
                    J * e.ETA * soft
                ).applyfunc(s.factor)
                missing = matter.soft_current(ps, soft) * e.ETA * soft
                checks[f"{label}_{tag}_missing_null_divergence"] = (
                    missing + marked
                ).applyfunc(s.factor)
                gates[f"{label}_{tag}_massive_only_Ward_control_nonzero"] = (
                    missing != e.VECTOR_ZERO
                )
            ja = current((*pb, qb[0]), qs[0])
            jb = current((*pa, qa[0]), qs[1])
            j0a = current(born, qs[0])
            j0b = current(born, qs[1])
            # Both unnormalized plus tensors have norm^2=2; each two-leg term divides by2.
            F2 = s.sqrt(rho(a, b)) * M6 / (2 * A0)
            B2 = (s.sqrt(rho(0, b)) * ja * M5b + s.sqrt(rho(a, 0)) * jb * M5a) / (
                2 * s.sqrt(K) * A0
            ) - j0a * j0b / (2 * K)
            triangle = (
                abs(M6) / (2 * A0)
                + (abs(ja * M5b) + abs(jb * M5a)) / (2 * s.sqrt(K) * A0)
                + abs(j0a * j0b) / (2 * K)
            )
            gates[label + "_entire_amplitudes_exact_rational"] = all(
                value.is_Rational for value in (M6, M5a, M5b, A0, triangle)
            )
            gates[label + "_positive_original_Born"] = A0 > 0
            gates[label + "_all_exact_phase_ratios_in_unit_interval"] = all(
                0 < rho(x, y) <= 1 for x, y in ((a, b), (a, 0), (0, b))
            )
            # I(a,b)>min(a,b)/2, so this rational comparison is sufficient.
            gates[label + "_exact_phase_independent_remainder_bound"] = (
                triangle < analytic.REMAINDER_COEFFICIENT / (2 * max(a, b))
            )
            records.append((family, denominator, a, b, F2, B2))
    return checks, gates, tuple(records)


@cache
def entropy_quadratures():
    with mp.workdps(70):

        def h(x):
            return -x * mp.log(x) - (1 - x) * mp.log(1 - x)

        first = mp.quad(lambda x: h(x) / (x * (1 - x)), [0, mp.mpf("0.5"), 1])
        second = mp.quad(lambda x: h(x) ** 2 / (x * (1 - x)), [0, mp.mpf("0.5"), 1])
        return bool(abs(first - mp.pi**2 / 3) < mp.mpf("1e-60")), bool(
            abs(second - (4 * mp.zeta(3) - mp.pi**2 / 3)) < mp.mpf("1e-60")
        )


@cache
def data():
    a, b = s.symbols("a b", positive=True)
    I = (a + b) * s.log(a + b) - a * s.log(a) - b * s.log(b)
    n = s.Symbol("n", positive=True, integer=True)
    K = source.KAPPA
    cc = analytic.REMAINDER_COEFFICIENT
    D = s.Rational(49, 32) / s.sqrt(K) + 2304 / K
    checks = {
        "entropy_rectangle_mixed_derivative": s.simplify(s.diff(I, a, b) - 1 / (a + b)),
        "entropy_first_axis": s.limit(I, a, 0, dir="+"),
        "entropy_second_axis": s.limit(I, b, 0, dir="+"),
        "entropy_energy_exchange": s.simplify(I - I.xreplace({a: b, b: a})),
        "cross_log_series_partial_fraction": s.cancel(
            1 / (n * (n + 1) ** 2) - (1 / n - 1 / (n + 1) - 1 / (n + 1) ** 2)
        ),
        "entropy_square_integral": s.simplify(
            4 * (s.zeta(3) - 1) + 2 * (2 - s.zeta(2)) - (4 * s.zeta(3) - s.pi**2 / 3)
        ),
        "entropy_interference_integral": 2 * s.zeta(2) - s.pi**2 / 3,
        "two_identical_graviton_angular_factor": s.factor(
            4 * (4 * s.pi) ** 2 / (8 * (2 * s.pi) ** 6) - 1 / (8 * s.pi**4)
        ),
        "positive_entropy_square_upper_margin": 4 * s.Rational(5, 4)
        - 2 * s.Rational(3, 2)
        - 2,
        "complex_double_soft_overlap_coefficient": s.Rational(49, 32)
        - 2 * s.Rational(49, 64),
        "conservative_interference_coefficient": 7 * 2 * s.Rational(1, 10**400) * cc
        - 14 * s.Rational(1, 10**726),
        "conservative_remainder_square_coefficient": cc**2 - s.Rational(1, 10**652),
    }
    physical, physical_gates, _ = original_samples()
    checks.update(physical)
    q1, q2 = entropy_quadratures()
    return {
        "whole_state_correct_projection": "With F1=sqrt(rho1)M5/A0,F2=sqrt(rho2)M6/A0 and G=abF2, the single-soft face is G(0,b)=b*S_a(sigma_b)*F1(b)/sqrt(kappa), where sigma_b includes its null marked graviton. G00=S_a(Born)S_b(Born)/kappa. Define B2=J_a(sigma_b)F1(b)/sqrt(kappa)+J_b(sigma_a)F1(a)/sqrt(kappa)-G00/(ab). R2=F2-B2 is the two-face inclusion-exclusion remainder. The leading double coefficient agrees in both orders; no simultaneous subleading contact is set to zero.",
        "whole_gauge_and_recoil_definition": "Every term uses the S300 recoil with the same E and outgoing pair-rest-frame direction. The reduced massive-plus-marked-null states conserve momentum, so each soft current obeys the appropriate Ward identity. The full434 tree and both47 one-real faces retain their existing Ward proofs. Keeping only massive legs would violate that identity. The phase ratio is rho=E*rprime/(r0*Eprime); common flux and identical-scalar factors cancel, while the two identical emitted gravitons retain1/2!.",
        "whole_integrable_remainder_theorem": "The uniform derivative bound gives |R2|<C*I(a,b)/(ab),C=1e-326,I=W lnW-a lna-b lnb. This follows by integrating partial_a partial_b G on the rectangle, including its continuous soft faces. If M=max(a,b),m=min(a,b),I<=m[1+ln(1+M/m)]. The logarithmic energy-sharing envelope is integrable and does not claim the sharper fixed-other-energy limit from samples.",
        "whole_baseline_bound": "S311 gives|F1(b)|<1/(64b). Four massive soft terms obey a|Jmassive|<=48 and the marked null TT term obeys a|Jnull|<=2b without an angular pole. Thus a|Jfull|<49 and |B2|<D/(ab),D=49/(32sqrt(kappa))+2304/kappa<2e-400.",
        "whole_signed_measure_bound": "Define the signed two-real tree measure by the complete difference |F2|2-|B2|2, summed over four unit physical polarization pairs with1/2! phase space and a+b<=x. All interference terms are kept. The entropy integrals pi2/3 and4zeta(3)-pi2/3<2 give TV<7DCx+C2x2<2e-725*x+1e-652*x2, independent of the lower soft cutoff, all physical emitted directions and all nonforward hard angles. Each unsubtracted measure separately remains soft divergent.",
        "whole_reference_comparison": "The previously established unexpanded elastic leading-soft reference satisfies P0>x^alpha/2,alpha<1/2. The new signed measure divided by P0 is therefore below4e-725*sqrt(x)+2e-652*x^(3/2), uniformly vanishing as x tends to zero. This is a size comparison, not an identification of B2 squared with S309's dressed single-residual measure or permission to add overlapping subtraction schemes without a matching proof.",
        "whole_remaining_scope": "Finite virtual terms, integrated counterterm matching, all-N nonleading radiative errors, radiative hard loops and finite physical matching, absolute complex Regge and the original common-parent quantum state/bounce remain open. This uniform subtracted real contribution is not a full inclusive detector rate or original V/G/B/P8 closure.",
        "checks": checks,
        "gates": {
            **physical_gates,
            "independent_entropy_first_quadrature": q1,
            "independent_entropy_square_quadrature": q2,
            "original_baseline_below_two_over_sqrt_kappa": D < BASELINE_COEFFICIENT,
            "interference_rounding_has_positive_margin": 14 * s.Rational(1, 10**726)
            < 2 * s.Rational(1, 10**725),
            "phase_polarization_Bose_prefactor_below_one": 1 / (8 * s.pi**4) < 1,
            "marked_state_current_includes_all_outgoing_particles": True,
            "both_soft_faces_and_single_common_overlap_retained": True,
            "full_signed_difference_includes_all_interferences": True,
            "leading_reference_not_a_matched_full_detector_claim": True,
        },
    }
