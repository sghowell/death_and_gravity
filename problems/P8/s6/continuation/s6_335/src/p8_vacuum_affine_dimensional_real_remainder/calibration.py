"""Exact original-parameter recoil calibrations, not proofs by finite sampling."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor, vertices
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import recoil
from p8_vacuum_affine_one_newton_inclusive_assembly import forward

from . import source


@cache
def data():
    checks = {}
    gates = {}

    def exact(name, value):
        checks[name] = s.factor(s.together(value))

    def positive(name, value):
        gates[name] = s.factor(value).is_positive is True

    fn = lambda A: s.trace(s.conjugate(A).T * A)
    B = s.Integer(4000000000)
    eta = tensor.ETA
    cases = 0
    for t0 in (s.Integer(2), s.Integer(3)):
        E0 = (t0 + 1 / t0) / 2
        tprime = t0 - s.Rational(1, 4096)
        Ep = (tprime + 1 / tprime) / 2
        omega = E0 - Ep * Ep / E0
        for n, u in (
            (
                (s.Rational(3, 5), s.Rational(4, 5), 0),
                (0, s.Rational(3, 5), s.Rational(4, 5)),
            ),
            ((0, 0, 1), (s.Rational(3, 5), s.Rational(4, 5), 0)),
        ):
            ps, k, pborn = recoil.momenta(E0, omega, n, u)
            ps, k = tensor.kinematics(ps, k)
            basis = tensor.frame(k)
            transfers = [
                ((pborn[0] + pborn[j]).T * eta * (pborn[0] + pborn[j]))[0]
                for j in (1, 2, 3)
            ]
            tau = min(-transfers[1], -transfers[2])
            dd = min(s.S.One, tau)
            assert 0 < omega < dd / 192
            J0 = sum(
                (
                    s.Matrix([(p.T * eta * f)[0] for f in basis])
                    * s.Matrix([(p.T * eta * f)[0] for f in basis]).T
                    / (p.T * eta * k)[0]
                    for p in pborn
                ),
                s.zeros(2),
            )
            S = tensor.matter_core(ps, k, basis, 0, 1, 0)
            AS = -sum(1 / z for z in transfers)
            GS = tensor.gravity_core(ps, k, basis)
            G0 = forward.gravity_born(4 * E0 * E0, u[2], 1)
            rS = S - AS * J0
            positive(
                "physical_auxiliary_remainder_" + str(cases),
                (45216 / dd**2) ** 2 - fn(rS),
            )
            exact(
                "physical_Born_gravity_" + str(cases),
                vertices.component_engine(4)["born"](pborn) - G0,
            )
            extend = lambda p: s.Matrix([*p, 0, 0])
            exact(
                "physical_D6_Born_" + str(cases),
                vertices.component_engine(6)["born"]([extend(p) for p in pborn])
                - G0
                - AS,
            )
            positive("physical_auxiliary_Born_" + str(cases), AS)
            positive("physical_gravity_Born_gap_" + str(cases), G0 - 8 / dd)
            AM = forward.matter_born(
                4 * E0 * E0, u[2], source.HEAVY_MASS2, source.CUBIC
            )
            exact(
                "original_matter_Born_" + str(cases),
                source.CONTACT
                + source.CUBIC**2 * sum(1 / (source.HEAVY_MASS2 - z) for z in transfers)
                - AM,
            )
            M = tensor.matter_core(
                ps, k, basis, source.HEAVY_MASS2, source.CUBIC, source.CONTACT
            )
            A0 = AM + G0 / source.KAPPA
            for ep in (s.S.Zero, s.Rational(1, 32), s.Rational(1, 8)):
                tt = ep / (1 + ep)
                GD = GS + 2 * tt * S
                G0D = G0 + 2 * tt * AS
                rg = GD - G0D * J0
                positive(
                    "physical_canonical_gravity_remainder_"
                    + str(cases)
                    + "_"
                    + str(ep),
                    (s.Integer(30000000000) / dd**2) ** 2 - fn(rg),
                )
                AD = AM + G0D / source.KAPPA
                combined = (M + GD / source.KAPPA) / AD - J0
                positive(
                    "original_full47_normalized_remainder_"
                    + str(cases)
                    + "_"
                    + str(ep),
                    (B / dd) ** 2 - fn(combined),
                )
                positive(
                    "original_Born_ratio_upper_" + str(cases) + "_" + str(ep),
                    s.Rational(19, 18) - AD / A0,
                )
            cases += 1
    low_count = cases
    cases = 0
    for t0 in (s.Integer(2), s.Integer(3)):
        E = (t0 + 1 / t0) / 2
        tprime = t0 - s.Rational(1, 12)
        Ep = (tprime + 1 / tprime) / 2
        omega = E - Ep**2 / E
        r = s.Rational(3, 100)
        for n, u in (
            (
                (s.Rational(3, 5), s.Rational(4, 5), 0),
                (0, s.Rational(3, 5), s.Rational(4, 5)),
            ),
            ((0, 0, 1), (2 * r / (1 + r * r), 0, (1 - r * r) / (1 + r * r))),
            ((0, 0, -1), (2 * r / (1 + r * r), 0, -(1 - r * r) / (1 + r * r))),
        ):
            ps, k, born = recoil.momenta(E, omega, n, u)
            ps, k = tensor.kinematics(ps, k)
            basis = tensor.frame(k)
            q = [
                ((born[0] + born[j]).T * eta * (born[0] + born[j]))[0]
                for j in (1, 2, 3)
            ]
            dd = min(s.S.One, -q[1], -q[2])
            m = dd / 30000
            assert dd / 192 < omega <= s.Rational(1, 8)
            for left, right in tensor.PARTS:
                for part in (left, right):
                    P = sum((ps[i] for i in part), s.zeros(4, 1))
                    positive(
                        "high_actual_physical_pair_gap_" + str(cases) + "_" + str(part),
                        abs((P.T * eta * P)[0]) - m,
                    )
            S = tensor.matter_core(ps, k, basis, 0, 1, 0)
            G = tensor.gravity_core(ps, k, basis)
            M = tensor.matter_core(
                ps, k, basis, source.HEAVY_MASS2, source.CUBIC, source.CONTACT
            )
            G0 = forward.gravity_born(4 * E * E, u[2], 1)
            AS = -sum(1 / z for z in q)
            AM = forward.matter_born(4 * E * E, u[2], source.HEAVY_MASS2, source.CUBIC)
            for e in (s.S.Zero, s.Rational(1, 32), s.Rational(1, 8)):
                theta = e / (1 + e)
                GD = G + 2 * theta * S
                positive(
                    "high_actual_canonical_GD_" + str(cases) + "_" + str(e),
                    (s.Integer(9000000000) / m**2) ** 2 - fn(GD),
                )
                AD = AM + (G0 + 2 * theta * AS) / source.KAPPA
                normalized = (M + GD / source.KAPPA) / AD
                positive(
                    "high_original_all47_normalized_" + str(cases) + "_" + str(e),
                    (s.Integer(2000000000000000000) / dd) ** 2 - fn(normalized),
                )
            cases += 1
    return {
        "checks": checks,
        "gates": gates,
        "whole_calibration_domain": {
            "low_recoil_states": low_count,
            "high_recoil_states": cases,
            "epsilon_values": (s.S.Zero, s.Rational(1, 32), s.Rational(1, 8)),
            "all_original_source_parameters_retained": True,
        },
        "whole_verification_boundary": "Exact rational original-source recoils with low and high energies, noncoplanar directions and near-forward/backward Born transfers. Full canonical trace response and all47 interferences are retained. The finite calibrations test formulas and margins; the uniform conclusions require the general symbolic soft identities and written D6 component/physical-gap bounds.",
    }
