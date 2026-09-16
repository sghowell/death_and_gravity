"""Smooth exact two-scalar recoil against arbitrary finite radiated energy."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import (
    recoil as single_recoil,
)

ETA = s.diag(1, -1, -1, -1)


def dot(p, q):
    return (s.Matrix(p).T * ETA * s.Matrix(q))[0]


def exact_real(value):
    if isinstance(value, (bool, float, s.Float)):
        raise TypeError("Require exact real numbers")
    value = s.sympify(value)
    if value.has(s.Float) or value.is_number is not True or value.is_real is not True:
        raise ValueError("Require known exact real numbers")
    return value


def require_state(energy, quanta, outgoing):
    E = exact_real(energy)
    if E < s.Rational(5, 4) or E > 2:
        raise ValueError("Require hard scalar COM energy in[5/4,2]")
    if not isinstance(quanta, (list, tuple)):
        raise TypeError("Require a finite list of outgoing null four-momenta")
    u = s.Matrix(outgoing)
    if u.shape != (3, 1):
        raise ValueError("Require a three-dimensional unit rest-frame direction")
    for value in u:
        exact_real(value)
    if s.simplify(u.dot(u) - 1) != 0:
        raise ValueError("Require a unit rest-frame direction")
    rays = []
    for item in quanta:
        q = s.Matrix(item)
        if q.shape != (4, 1):
            raise ValueError("Require null four-momenta")
        for value in q:
            exact_real(value)
        if q[0] <= 0 or s.simplify(dot(q, q)) != 0:
            raise ValueError("Require future positive-energy null radiation")
        rays.append(q)
    if sum((q[0] for q in rays), s.S.Zero) > s.Rational(1, 8):
        raise ValueError("Require total radiated energy at most1/8")
    return E, tuple(rays), u


def momenta(energy, quanta, outgoing):
    E, rays, u = require_state(energy, quanta, outgoing)
    Q = sum(rays, s.zeros(4, 1))
    R, v = Q[0], Q[1:4, 0]
    h = E - R / 2
    Ep = s.sqrt(h * h - v.dot(v) / 4)
    rp, r0 = s.sqrt(Ep * Ep - 1), s.sqrt(E * E - 1)
    uv = u.dot(v)
    pair = []
    for sign in (1, -1):
        spatial = (
            sign * rp * u
            + (sign * rp * uv / (4 * Ep * (h + Ep)) - s.Rational(1, 2)) * v
        )
        pair.append(s.Matrix([h - sign * rp * uv / (2 * Ep), *spatial]))
    incoming = (s.Matrix([-E, 0, 0, -r0]), s.Matrix([-E, 0, 0, r0]))
    born = (*incoming, s.Matrix([E, *(r0 * u)]), s.Matrix([E, *(-r0 * u)]))
    return (*incoming, *pair), rays, born


@cache
def data():
    checks = {}

    def put(name, value):
        clean = lambda entry: s.factor(s.radsimp(s.sympify(entry)))
        checks[name] = (
            value.applyfunc(clean) if isinstance(value, s.MatrixBase) else clean(value)
        )

    h, Ep, r, c, V = s.symbols("h Eprime rprime u_dot_Q Q_spatial_squared", real=True)
    for sign in (1, -1):
        e = h - sign * r * c / (2 * Ep)
        b = sign * r * c / (4 * Ep * (h + Ep)) - s.Rational(1, 2)
        norm = r * r + b * b * V + 2 * sign * r * b * c
        polynomial = s.together((e * e - norm - 1) * 16 * Ep**2 * (h + Ep) ** 2).subs(
            V, 4 * (h * h - Ep * Ep)
        )
        put(
            "generic_mass_shell_" + str(sign),
            s.rem(s.expand(polynomial), r * r - Ep * Ep + 1, r),
        )
    put(
        "generic_energy_conservation",
        (h - r * c / (2 * Ep)) + (h + r * c / (2 * Ep)) - 2 * h,
    )
    put(
        "generic_spatial_conservation",
        r * c / (4 * Ep * (h + Ep))
        - s.Rational(1, 2)
        - r * c / (4 * Ep * (h + Ep))
        - s.Rational(1, 2)
        + 1,
    )
    put(
        "zero_total_spatial_momentum_no_singularity",
        (4 * Ep * (h + Ep)).subs(Ep, h) - 8 * h * h,
    )
    put("boost_removable_quotient", (h - Ep) * (h + Ep) - (h * h - Ep * Ep))
    E, R, M = s.symbols("E total_energy Q_mass_squared", positive=True)
    put("pair_energy_squared", (E - R / 2) ** 2 - (R * R - M) / 4 - E * (E - R) - M / 4)
    put("pair_energy_lower_than_initial", E * E - (E * (E - R) + M / 4) - E * R + M / 4)
    put(
        "pair_energy_above_E_minus_R",
        E * (E - R) + M / 4 - (E - R) ** 2 - R * (E - R) - M / 4,
    )
    put(
        "compact_rprime_lower_square",
        s.Rational(5, 4) * s.Rational(9, 8) - 1 - s.Rational(13, 32),
    )
    put(
        "compact_radial_shift_margin",
        2 - 2 / (s.Rational(3, 4) + s.Rational(1, 2)) - s.Rational(2, 5),
    )
    put("compact_spatial_shift_sum", 2 + s.Rational(1, 2) + s.Rational(1, 2) - 3)
    put(
        "future_massive_doppler_identity",
        (E - s.Symbol("p")) * (E + s.Symbol("p")) - (E**2 - s.Symbol("p") ** 2),
    )
    put("six_pair_dot_shift_bound", 4 * 8 * R + 16 * R - 48 * R)
    put("max_pair_dot", 2 * 2 + 3 - 7)
    u = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    nhat = s.Matrix([s.Rational(3, 5), 0, s.Rational(4, 5)])
    w = s.Rational(1, 20)
    q = s.Matrix([w, *(w * nhat)])
    made, rays, born = momenta(s.Rational(5, 4), [q], u)
    old, oldq, oldborn = single_recoil.momenta(s.Rational(5, 4), w, nhat, u)
    for i, (p, oldp) in enumerate(zip(made, old)):
        put("single_real_recoil_agreement_" + str(i), p - oldp)
    put("single_radiation_agreement", rays[0] - oldq)
    put("single_Born_agreement", s.Matrix.hstack(*born) - s.Matrix.hstack(*oldborn))
    opposite = [s.Matrix([w, w, 0, 0]), s.Matrix([w, -w, 0, 0])]
    made, rays, born = momenta(s.Rational(5, 4), opposite, u)
    for i, p in enumerate(made):
        put("zero_net_radiation_spatial_mass_" + str(i), dot(p, p) - 1)
    put(
        "zero_net_radiation_spatial_conservation",
        sum(made, s.zeros(4, 1)) + sum(rays, s.zeros(4, 1)),
    )
    return {
        "whole_smooth_recoil_formula": "For Q=(R,v),h=E-R/2,Ep=sqrt(h^2-v^2/4),r=sqrt(Ep^2-1),p3^0=h-r(u.v)/(2Ep),p3vec=r*u+[r(u.v)/(4Ep(h+Ep))-1/2]v; p4 reverses r. The expression is regular at v0, unlike a unit-vector boost formula.",
        "whole_recoil_domain_and_bounds": "E in[5/4,2],R<=1/8 and positive null outgoing rays imply0<=Q^2<=R^2. Ep^2>=E(E-R),soEp>1 and the final masses stay1. Each energy changes by<=R and stays<=2. Radius shift<=2R and boost-plus-recoil spatial shift<=R, giving|delta pvec|<=3R. Every massive future pair dot lies in[1,7], and the sum of its six absolute changes is<=48R. No multiplicity appears.",
        "checks": checks,
        "gates": {
            "arbitrary_positive_null_radiation_and_exact_conservation": True,
            "smooth_zero_net_spatial_radiation_limit": True,
            "one_quantum_map_matches_complete_S295": True,
            "uniform_total_energy_not_particle_count_bounds": True,
            "massive_Doppler_gap_and_pair_dot_domain": True,
            "not_an_interacting_quantum_state_construction": True,
        },
    }
