"""Complete ordered Ward cancellation, scalar projectors and clock-chart contact."""

from functools import cache

import sympy as s
from p8_vacuum_affine_full_adm_vertices import chart
from p8_vacuum_affine_prepared_ward_reconstruction import kinematics as k
from p8_vacuum_affine_prepared_ward_reconstruction import ward

from . import scalar

t, x, y, z = k.COORDS
a, H = k.a, k.H
p = s.Symbol("positive_transfer", positive=True)
Pi = s.diag(1, 0, 0)


def projected(zeta, eta, c):
    return 2 * (zeta - H * eta) * s.eye(3) - 2 * c * Pi


def ordered_local(nD, zD, bD, etaD, cD, nG, zG, etaG, cG, A, B):
    source = (
        2 * nD * (s.diff(A, t) * etaG + A * cG)
        - 2 * (a * a * B + A) * bD * etaG
        - 6 * a * a * zD * ((s.diff(B, t) + 2 * H * B) * etaG + B * nG + B * cG)
    )
    traceG = 6 * (zG - H * etaG) - 2 * cG
    detector = a * a * B * (etaD * s.diff(traceG, t) - cD * traceG)
    return source, detector


@cache
def ward_data():
    A, B = s.Function("A")(t), s.Function("B")(t)
    density = s.diag(A, B, B, B)
    eta = s.Function("eta")(*k.COORDS)
    chi = s.Matrix([s.Function("chi" + str(j))(*k.COORDS) for j in range(3)])
    xi = s.Matrix([eta, *chi])
    n = s.Function("n")(*k.COORDS)
    beta = s.Matrix([s.Function("beta" + str(j))(*k.COORDS) for j in range(3)])
    Q = s.Matrix(
        3, 3, lambda i, j: s.Function("Q" + str(min(i, j)) + str(max(i, j)))(*k.COORDS)
    )
    divchi = s.trace(k.spatial_jac(chi))
    trace = s.trace(Q)
    expected_source = (
        2 * n * (s.diff(A, t) * eta + A * divchi)
        + 2 * (a * a * B + A) * beta.dot(k.grad(eta))
        - a
        * a
        * trace
        * ((s.diff(B, t) + 2 * H * B) * eta + B * s.diff(eta, t) + B * divchi)
    )
    expected_detector = a * a * B * (eta * s.diff(trace, t) + chi.dot(k.grad(trace)))
    checks = {
        "full_general_source_density_and_ADM_contact": s.factor(
            ward.source_integrand((n, beta, Q), xi, density) - expected_source
        ),
        "full_synchronous_detector_Lie_and_ADM_contact": s.factor(
            ward.detector_integrand(xi, (s.S.Zero, s.zeros(3, 1), Q), density)
            - expected_detector
        ),
    }
    nD, zD, bD, etaD, cD = [
        s.Function(name)(t) for name in ("nD", "zD", "bD", "etaD", "cD")
    ]
    nG, zG, bG, etaG, cG = [
        s.Function(name)(t) for name in ("nG", "zG", "bG", "etaG", "cG")
    ]
    eD, eG = s.exp(-s.I * p * x), s.exp(s.I * p * x)
    D = (nD * eD, s.Matrix([s.I * bD * eD / p, 0, 0]), 2 * zD * eD * s.eye(3))
    xiG = s.Matrix([etaG * eG, -s.I * cG * eG / p, 0, 0])
    xiD = s.Matrix([etaD * eD, s.I * cD * eD / p, 0, 0])
    QG = projected(zG, etaG, cG) * eG
    sw, dw = ordered_local(nD, zD, bD, etaD, cD, nG, zG, etaG, cG, A, B)
    actual_s = s.factor(
        ward.source_integrand(D, xiG, density).subs(s.diff(etaG, t), nG)
    )
    actual_d = s.factor(
        ward.detector_integrand(xiD, (s.S.Zero, s.zeros(3, 1), QG), density)
    )
    checks.update(
        {
            "opposite_Fourier_full_ordered_source": s.factor(actual_s - sw),
            "opposite_Fourier_synchronous_detector": s.factor(actual_d - dw),
            "local_expression_no_inverse_momentum": s.factor(
                s.diff(actual_s + actual_d, p)
            ),
        }
    )
    for sign, et, cc, nn, bb, ee in (
        (1, etaG, cG, nG, bG, eG),
        (-1, etaD, cD, nD, bD, eD),
    ):
        xx = s.Matrix([et * ee, -s.I * sign * cc * ee / p, 0, 0])
        ng, bg, Qg = k.gauge(xx[0], xx[1:, 0])
        substitutions = {s.diff(et, t): nn, s.diff(cc, t): bb - p * p * et / a**2}
        checks[f"Fourier_sign_{sign}_lapse"] = s.expand(
            ng.subs(substitutions) - nn * ee
        )
        checks[f"Fourier_sign_{sign}_shift"] = s.simplify(
            bg.subs(substitutions) - s.Matrix([-s.I * sign * bb * ee / p, 0, 0])
        )
        checks[f"Fourier_sign_{sign}_spatial"] = s.expand(
            Qg - (2 * H * et * s.eye(3) + 2 * cc * Pi) * ee
        )
    pole = s.factor(
        s.expand(ward.pair(ward.lie_density(density, xiG), chart.first(*D, a))).coeff(
            p, -2
        )
    )
    checks["nonzero_deleted_chart_IR_pole_value"] = s.factor(
        pole - 2 * a * a * A * bD * s.diff(cG, t)
    )
    return {
        "general_ordered_source": expected_source,
        "synchronous_ordered_detector": expected_detector,
        "scalar_ordered_source": sw,
        "scalar_ordered_detector": dw,
        "missing_chart_inverse_square_control": pole,
        "prepared_projector": "For each side I=Iminus or Iplus, eta=I n, c=I b-|P|^2 I(a^-2 eta), chi=-i P c/|P|^2. Qsyn=2(zeta-H eta)I-2 c Phat Phat^T, with zeta=v+delta n. Both primitive derivatives are positive; source and detector use different endpoints.",
        "infrared": "The projector has norm1 away from P=0. Every shift inverse momentum cancels in the complete Ward expression: source beta_D appears only as beta_D dot grad eta_G, and synchronous detector chi_D only as spatial transport. A deleted source chart leaves the displayed nonzero inverse-square pole.",
        "boundary": "Use the original initial source germ and final detector germ. Cut the advanced synchronous detector inside the initial zero-output neighborhood as in S215; no detector time derivative enters M. The identity is for the covariant renormalized Gaussian current, not a sharp finite-band Ward identity.",
        "checks": checks,
        "gates": {
            "full_source_and_detector_identities_separate": True,
            "nonzero_deleted_chart_inverse_square_pole": pole != 0,
            "both_opposite_Fourier_leg_signs_checked": True,
            "no_inverse_transfer_in_final_scalar_forms": True,
            "no_arbitrary_zero_mean_stress": True,
            "no_finite_sharp_band_Ward_claim": True,
            "no_two_retarded_cross_kernels_identified": True,
        },
    }


@cache
def clock_data():
    N, D, nD, nG, pressure, rho = s.symbols("N delta nD nG pressure rho", real=True)
    R = 1 + 2 * D * (N**-2 - 1)
    Qclock = -s.log(R) / 2
    mixed = (8 * D**2 - 6 * D) * nD * nG * s.eye(3)
    density = s.diag(
        -(a**3) * rho / 2, -a * pressure / 2, -a * pressure / 2, -a * pressure / 2
    )
    contact = ward.pair(density, chart.first(0, s.zeros(3, 1), mixed, a))
    e = s.Symbol("epsilon", real=True)
    direction = s.Matrix(s.symbols("P0:3", real=True))
    projector = direction * direction.T / (direction.dot(direction))
    phase_map = s.Matrix(
        [
            scalar.lapse_numerator() / (2 * (scalar.J0 + scalar.dJ)),
            scalar.v
            + scalar.delta * scalar.lapse_numerator() / (2 * (scalar.J0 + scalar.dJ)),
            scalar.pv / 2,
        ]
    )
    checks = {
        "physical_clock_first_log_metric_jet": s.factor(
            s.diff(Qclock, N).subs(N, 1) - 2 * D
        ),
        "physical_clock_second_log_metric_jet": s.factor(
            s.diff(Qclock, N, 2).subs(N, 1) - (8 * D**2 - 6 * D)
        ),
        "full_extra_one_current_contact": s.factor(
            contact - 3 * a**3 * pressure * (8 * D**2 - 6 * D) * nD * nG / 2
        ),
        "projector_idempotent": s.simplify(projector * projector - projector),
        "projector_unit_trace": s.factor(s.trace(projector) - 1),
        "projector_Frobenius_one": s.factor(s.trace(projector.T * projector) - 1),
        "bounded_projector_along_any_ray": s.simplify(
            projector.subs({v: e * v for v in direction}, simultaneous=True) - projector
        ),
        "linear_phase_map_recovered": phase_map.jacobian(scalar.Z) * scalar.Z
        - phase_map,
    }
    return {
        "clock_extra_contact": contact,
        "phase_to_scalar_metric_comparison_map": phase_map,
        "clock_map": "The actual physical metric is g=C ghat+D du du with C=R^-1/2. N and contravariant shift are unchanged. In scalar logarithmic spatial coordinates Qphysical=2v I-(log R)/2 I, not just its first jet. Its second clock jet supplies the extra displayed one-current term beyond the full ADM chart contact.",
        "phase_scope": "Insert the current QG1 coefficient-sector linear lapse map only as a comparison pullback. The full Gaussian response makes the coupled lapse equation nonlocal. Neither this map nor its linear-family Hessian performs the full quantum Schur elimination.",
        "zero_transfer": "The nonzero-P scalar constraint has no literal global P0 shift counterpart. For R3 Sobolev fields the bounded projector can be assigned any value at the measure-zero origin; this does not identify that assignment with a physical homogeneous solution or prove direction-independent continuity there.",
        "checks": checks,
        "gates": {
            "nonzero_clock_second_chart_contact": contact != 0,
            "both_ADM_and_clock_chart_contacts_retained": True,
            "directional_order_zero_projector_not_continuous_at_origin": True,
            "global_homogeneous_constraint_kept_separate": True,
            "phase_comparison_not_full_quantum_reduction": True,
        },
    }
