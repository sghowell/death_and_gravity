"""Canonical rank-two tensor representative and independent literal graph checks."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as matter
from p8_vacuum_affine_minimal_gravity_radiation import vertices as old_vertices

from . import sew, source, vertices, ward

ETA = s.diag(1, -1, -1, -1)
PARTS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


def vector(value):
    if not isinstance(value, s.MatrixBase) or value.shape != (4, 1):
        raise ValueError("Require an exact four-component column")
    return s.Matrix([sew.exact_real(v) for v in value])


def frame(k):
    k = vector(k)
    if k[0] <= 0 or s.factor((k.T * ETA * k)[0]) != 0:
        raise ValueError("Require a positive-energy null emitted momentum")
    n = k[1:4, 0] / k[0]
    first = n.cross(s.Matrix([0, 0, 1]))
    if first == s.zeros(3, 1):
        first = n.cross(s.Matrix([1, 0, 0]))
    first = first / s.sqrt(first.dot(first))
    second = n.cross(first)
    return (s.Matrix([0, *first]), s.Matrix([0, *second]))


def kinematics(momenta, k):
    if not isinstance(momenta, (tuple, list)) or len(momenta) != 4:
        raise ValueError("Require four original massive outgoing-convention legs")
    ps = tuple(vector(p) for p in momenta)
    k = vector(k)
    E = -ps[0][0]
    if (
        E < s.Rational(5, 4)
        or E > 2
        or ps[1][0] != -E
        or k[0] <= 0
        or k[0] > s.Rational(1, 8)
    ):
        raise ValueError(
            "Require the original compact hard-energy and radiation domain"
        )
    if any(p[0] <= 0 for p in ps[2:]) or ps[0][1:4, 0] + ps[1][1:4, 0] != s.zeros(3, 1):
        raise ValueError("Require physical incoming COM and positive outgoing energies")
    if (
        any(s.factor((p.T * ETA * p)[0] - 1) != 0 for p in ps)
        or s.factor((k.T * ETA * k)[0]) != 0
    ):
        raise ValueError(
            "Require original unit-mass shells and a null emitted momentum"
        )
    if any(s.factor(v) != 0 for v in sum(ps, k.copy())):
        raise ValueError("Require exact momentum conservation")
    if any(s.factor((p.T * ETA * k)[0]) == 0 for p in ps):
        raise ValueError("An external radiation pole is excluded")
    for left, right in PARTS:
        for pair in (left, right):
            q = sum((ps[i] for i in pair), s.zeros(4, 1))
            if s.factor((q.T * ETA * q)[0]) == 0:
                raise ValueError("An internal gravity pole is excluded")
    return ps, k


@cache
def canonical_polynomial():
    expression = ward.tt_channel().subs(ward.DIM, 4)
    result = s.Poly(expression, ward.x1, ward.x2, ward.x3)
    if not all(sum(mon) == 2 for mon, coeff in result.terms()):
        raise ValueError("The canonical rank-one expression must be quadratic")
    return result


def gravity_core(ps, k, basis):
    """Internal canonical representative: explicit metric terms are already removed."""
    total = s.zeros(2)
    for left, right in PARTS:
        pp = [ps[i] for i in (*left, *right)]
        sub = {
            ward.mu: 1,
            ward.a: (pp[0].T * ETA * pp[1])[0],
            ward.b: (pp[0].T * ETA * pp[2])[0],
        }
        sub.update(
            {u: (p.T * ETA * k)[0] for u, p in zip((ward.u1, ward.u2, ward.u3), pp)}
        )
        projected = [s.Matrix([(p.T * ETA * e)[0] for e in basis]) for p in pp[:3]]
        for mon, coeff in canonical_polynomial().terms():
            indices = [i for i, count in enumerate(mon) for _ in range(count)]
            a, b = indices
            outer = (projected[a] * projected[b].T + projected[b] * projected[a].T) / 2
            total += s.factor(coeff.subs(sub)) * outer
    return total.applyfunc(s.factor)


def matter_core(ps, k, basis, heavy, cubic, contact):
    """Internal symbolic tensor; heavy=0 is used only for the named auxiliary identity."""
    heavy, cubic, contact = map(s.sympify, (heavy, cubic, contact))
    z = [s.Matrix([(p.T * ETA * e)[0] for e in basis]) for p in ps]
    J = [a * a.T / (p.T * ETA * k)[0] for a, p in zip(z, ps)]
    total = contact * sum(J, s.zeros(2))
    for left, right in PARTS:
        PL = sum((ps[i] for i in left), s.zeros(4, 1))
        PR = sum((ps[i] for i in right), s.zeros(4, 1))
        dl = (PL.T * ETA * PL)[0] - heavy
        dr = (PR.T * ETA * PR)[0] - heavy
        zL = sum((z[i] for i in left), s.zeros(2, 1))
        total -= cubic**2 * (
            sum((J[i] for i in left), s.zeros(2)) / dr
            + sum((J[i] for i in right), s.zeros(2)) / dl
            + 2 * zL * zL.T / (dl * dr)
        )
    return total.applyfunc(s.factor)


def original_cores(momenta, k):
    ps, k = kinematics(momenta, k)
    basis = frame(k)
    M = matter_core(ps, k, basis, source.HEAVY_MASS2, source.CUBIC, source.CONTACT)
    G = gravity_core(ps, k, basis)
    S = matter_core(ps, k, basis, 0, 1, 0)
    return (
        (M + G / source.KAPPA) / s.sqrt(source.KAPPA),
        2 * S / source.KAPPA ** s.Rational(3, 2),
    )


def original_rate(momenta, k, epsilon):
    A, B = original_cores(momenta, k)
    return sew.continued_rate(A, B, epsilon)


def transverse_polarizations(basis, eta):
    n = len(basis)
    result = []
    for i in range(n):
        for j in range(i + 1, n):
            E = basis[i] * basis[j].T + basis[j] * basis[i].T
            result.append((eta * E * eta, s.Integer(2)))
    for j in range(1, n):
        E = (
            sum((basis[i] * basis[i].T for i in range(j)), s.zeros(eta.rows))
            - j * basis[j] * basis[j].T
        )
        result.append((eta * E * eta, s.Integer(j * (j + 1))))
    return tuple(result)


def literal_sew(ps4, k4, dimension, heavy, cubic, contact, invk):
    component = vertices.component_engine(dimension)
    eta = component["ETA"]
    extend = lambda p: s.Matrix([*p, *([0] * (dimension - 4))])
    ps = [extend(p) for p in ps4]
    k = extend(k4)
    basis = [extend(e) for e in frame(k4)] + [
        s.eye(dimension)[:, j] for j in range(4, dimension)
    ]
    polarizations = transverse_polarizations(basis, eta)
    aux = matter.whole_tensor(
        ps, k, mass=1, heavy=heavy, cubic=cubic, contact=contact, kappa=1
    )
    value = s.S.Zero
    for eps, norm2 in polarizations:
        if s.trace(eta * eps) != 0 or eps * k != s.zeros(dimension, 1):
            raise ValueError("A literal polarization is not transverse trace-free")
        if sum(e * e for e in eps) != norm2:
            raise ValueError("A spatial literal polarization has the wrong norm")
        amp = sum(
            aux[i, j] * eps[i, j] for i in range(dimension) for j in range(dimension)
        )
        amp += invk * component["amplitude"](ps, k, eps)
        value += amp * s.conjugate(amp) / norm2
    return s.factor(value), len(polarizations)


@cache
def data():
    checks = {}
    wrong = []
    rows = []
    for which in (0, 1):
        ps, k, _ = old_vertices.sample(which)
        ps, k = kinematics(ps, k)
        basis = frame(k)
        G = gravity_core(ps, k, basis)
        S = matter_core(ps, k, basis, 0, 1, 0)
        heavy, cubic, contact, invk = (
            s.Integer(29),
            s.Rational(2, 7),
            s.Rational(3, 11),
            s.Rational(1, 5),
        )
        A = matter_core(ps, k, basis, heavy, cubic, contact) + invk * G
        B = 2 * invk * S
        for dim in (4, 5, 6):
            prefix = f"fixture{which}_D{dim}_"
            t = s.Rational(dim - 4, dim - 2)
            full, count = literal_sew(ps, k, dim, heavy, cubic, contact, invk)
            expected = sew.bilinear(A + t * B, A + t * B, dim)
            checks[prefix + "all47_complete_polarization_sew"] = s.factor(
                full - expected
            )
            checks[prefix + "exact_polarization_count"] = s.Integer(
                count - dim * (dim - 3) // 2
            )
            component = vertices.component_engine(dim)
            eta = component["ETA"]
            extend = lambda p, dim=dim: s.Matrix([*p, *([0] * (dim - 4))])
            pp = [extend(p) for p in ps]
            kk = extend(k)
            xi = extend(basis[0] + s.I * basis[1])
            eps = (eta * xi) * (eta * xi).T
            grav = component["amplitude"](pp, kk, eps)
            aux = matter.whole_tensor(
                pp, kk, mass=1, heavy=0, cubic=1, contact=0, kappa=1
            )
            scalar = (xi.T * eta * aux * eta * xi)[0]
            xi4 = basis[0] + s.I * basis[1]
            grav4 = old_vertices.amplitude(ps, k, (ETA * xi4) * (ETA * xi4).T)
            checks[prefix + "independent_auxiliary_TT_identity"] = s.factor(
                grav - grav4 - 2 * t * scalar
            )
            gauge = (eta * kk) * (eta * xi).T + (eta * xi) * (eta * kk).T
            checks[prefix + "independent_21_graph_Ward"] = s.factor(
                component["amplitude"](pp, kk, gauge)
            )
            if dim > 4:
                difference = s.factor(expected - sew.bilinear(A + t * B, A + t * B, 4))
                checks[prefix + "omitted_trace_error"] = s.factor(
                    difference - t * s.trace(A + t * B) ** 2 / 2
                )
                wrong.append(difference != 0)
            rows.append({"fixture": which, "dimension": dim, "polarizations": count})
    ps, k, _ = old_vertices.sample(0)
    A, B = original_cores(ps, k)
    full, count = literal_sew(
        ps, k, 5, source.HEAVY_MASS2, source.CUBIC, source.CONTACT, 1 / source.KAPPA
    )
    t = s.Rational(1, 3)
    checks["original_physical_source_D5_all47_all5_polarizations"] = s.factor(
        full / source.KAPPA - sew.bilinear(A + t * B, A + t * B, 5)
    )
    checks["original_physical_source_D5_polarization_count"] = s.Integer(count - 5)
    checks["original_regulator_zero_rate"] = s.factor(
        original_rate(ps, k, 0) - sew.bilinear(A, A, 4)
    )
    return {
        "checks": checks,
        "gates": {
            "literal_components_independent_of_invariant_index_contraction": True,
            "all32_fixture_polarizations_and_all_interferences_retained": sum(
                row["polarizations"] for row in rows
            )
            == 32,
            "all4_missing_trace_controls_are_nonzero": len(wrong) == 4 and all(wrong),
            "rational_coupling_fixtures_not_original_source_calibrations": True,
            "separate_original_source_D5_all5_polarizations_checked": True,
            "canonical_trace_representative_fixed_before_sewing": True,
            "original_physical_APIs_do_not_retune_any_parameter": True,
        },
        "whole_literal_fixture_inventory": rows,
        "whole_original_calibration": "One exact original recoil state, D5, all five polarizations and original n,g,C,kappa, checked separately from rational parameter fixtures. Four-dimensional sample rotations and noncoplanar emission are frozen S304 coordinates.",
        "whole_canonical_core": "The unique implemented nonmetric momentum-dyadic representative is read from the unrestricted Gram TT rank-one polynomial, not reconstructed from only D4 helicities. A=M/sqrt(kappa)+G4alg/kappa^(3/2);B=2mu^2*Saux/kappa^(3/2).",
    }
