"""Full physical quadratic normal matter observable in the fixed spatial gauge."""

from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model as old
from p8_coupled_momentum import constraints, quadratic
from p8_coupled_vertices import physical
from p8_physical import geometry
from p8_physical import jets as j

from . import coefficients

CHANNELS = quadratic.CHANNELS


@cache
def at_time(point):
    point = constraints.exact(point)
    if abs(point) > sp.Rational(1, 100):
        raise ValueError("Require an exact time in the declared bounce strip")
    d = coefficients.data()
    return {
        "rows": {
            key: tuple(sp.factor(v.subs(old.u, point)) for v in row)
            for key, row in d["rows"].items()
        },
        "background": {
            key: sp.factor(v.subs(old.u, point))
            for key, v in old.coefficients()["background"].items()
        },
        "Jnew": sp.factor(d["Jnew"].subs(old.u, point)),
    }


def construct(ctx, f, t, tp, W, Pi, point=0):
    if not isinstance(ctx, j.Context) or ctx.n not in (2, 3, 4):
        raise TypeError("Require a two-through-four-leg labelled spatial context")
    if type(f) is not dict or set(f) != {"zeta", "scalar_p", "chi", "chi_p"}:
        raise ValueError("Require all four named scalar phase fields")
    if any(
        len(matrix) != 3 or any(len(row) != 3 for row in matrix) for matrix in (t, tp)
    ):
        raise ValueError("Require three by three tensor fields and momenta")
    if len(W) != 3 or len(Pi) != 3:
        raise ValueError("Require three vector fields and three vector momenta")
    fields = (
        tuple(f.values())
        + tuple(x for matrix in (t, tp) for row in matrix for x in row)
        + tuple(W)
        + tuple(Pi)
    )
    if any(
        not isinstance(value, j.Jet)
        or value.context is not ctx
        or not value.homogeneous(0).is_zero()
        or not (value - value.homogeneous(1)).is_zero()
        for value in fields
    ):
        raise ValueError(
            "Require zero-background linear fields in one labelled context"
        )
    for matrix in (t, tp):
        if not j.trace(matrix).is_zero() or any(
            not (matrix[a][b] - matrix[b][a]).is_zero()
            for a in range(3)
            for b in range(3)
        ):
            raise ValueError("Require symmetric tracefree tensor fields and momenta")
        if any(
            not sum(matrix[a][b].derivative(b) for b in range(3)).is_zero()
            for a in range(3)
        ):
            raise ValueError("Require transverse tensor fields and momenta")
    data = at_time(point)
    bg = data["background"]
    geo = geometry.derive(ctx, f["zeta"], t)
    momentum = constraints.derive(
        ctx,
        geo,
        f["zeta"],
        t,
        f["scalar_p"],
        tp,
        f["chi"],
        f["chi_p"],
        W,
        Pi,
        bg["H"],
        bg["ell"],
        order=min(2, ctx.n - 1),
    )
    values = physical.invariants(ctx, geo, momentum["momentum"], f, t, tp, W, Pi, bg)
    zero = (0,) * len(old.VARIABLES)
    A = data["rows"][zero]
    force = ctx.jet()
    linear_derivative = ctx.jet()
    for powers, row in data["rows"].items():
        if not any(powers):
            continue
        term = ctx.jet(1)
        for val, power in zip(values, powers, strict=True):
            term *= val**power
        force += row[1] * term
        linear_derivative += row[2] * term.homogeneous(1)
    n1 = -force.homogeneous(1) / A[2]
    n2 = -(force.homogeneous(2) + A[3] * n1 * n1 / 2 + linear_derivative * n1) / A[2]
    n = n1 + n2
    force_residual = A[2] * n + A[3] * n1 * n1 / 2 + force + linear_derivative * n1
    ell, hh = bg["ell"], bg["h"]
    dc, gradient = values[1], values[7]
    common = (
        ell * ell / 2
        + ell * dc
        + dc * dc / 2
        - 3 * ell * ell * n / (2 * hh)
        - 3 * ell * dc * n1 / hh
        + sp.Rational(3, 4) * ell * ell * (1 + 3 * hh) * n1 * n1 / hh**2
    )
    rho = common + gradient / 2
    pressure = common - gradient / 6
    N = ctx.jet(1) + n
    ratio = (hh - 1 + N.power(-2)) / hh
    direct = (ell + dc) ** 2 * ratio.power(
        sp.Rational(3, 2)
    ) / 2 + gradient * ratio.power(sp.Rational(1, 2)) / 2
    direct_p = (ell + dc) ** 2 * ratio.power(
        sp.Rational(3, 2)
    ) / 2 - gradient * ratio.power(sp.Rational(1, 2)) / 6
    truncate = lambda value: sum((value.homogeneous(d) for d in range(3)), ctx.jet())
    rho, pressure, direct, direct_p = map(truncate, (rho, pressure, direct, direct_p))
    checks = {
        "linear_stationary_lapse_force": force_residual.homogeneous(1).is_zero(),
        "quadratic_stationary_lapse_force": force_residual.homogeneous(2).is_zero(),
        "direct_physical_density_after_complete_lapse": (rho - direct).is_zero(),
        "direct_physical_pressure_after_complete_lapse": (
            pressure - direct_p
        ).is_zero(),
    }
    if not all(checks.values()):
        raise ValueError("A full physical matter reconstruction failed")
    return {
        "context": ctx,
        "geometry": geo,
        "momentum": momentum,
        "values": values,
        "n1": n1,
        "n2": n2,
        "rho": rho,
        "pressure": pressure,
        "checks": checks,
        "missing_if_lapse_truncated": -3 * ell * ell * n2 / (2 * hh),
    }


def pair(left, right, point=0):
    if (
        type(left) is not str
        or type(right) is not str
        or left not in CHANNELS
        or right not in CHANNELS
    ):
        raise ValueError("Require declared physical quadratic phase channels")
    point = constraints.exact(point)
    at_time(point)
    return _pair(left, right, point)


@cache
def _pair(left, right, point):
    ctx = constraints.context(((1, 0, 0), (-1, 0, 0)), symbolic_scale=True)
    result = construct(ctx, *quadratic.fields(ctx, left, right), point=point)
    return {
        **{
            key: sp.factor(result[key].coefficient(ctx.full))
            for key in ("rho", "pressure", "n2", "missing_if_lapse_truncated")
        },
        "checks": result["checks"],
    }


@cache
def matrices(point=0):
    point = constraints.exact(point)
    rho = sp.zeros(14)
    pressure = sp.zeros(14)
    missing = sp.zeros(14)
    checks = {}
    for i, left in enumerate(CHANNELS):
        for q, right in enumerate(CHANNELS[i:], start=i):
            data = pair(left, right, point)
            for key, M in (
                ("rho", rho),
                ("pressure", pressure),
                ("missing_if_lapse_truncated", missing),
            ):
                M[i, q] = data[key]
                if q != i:
                    M[q, i] = data[key].subs(
                        constraints.wave_scale, -constraints.wave_scale
                    )
            for key, value in data["checks"].items():
                checks[left + "_" + right + "_" + key] = value
    return {
        "rho": rho,
        "pressure": pressure,
        "missing_if_lapse_truncated": missing,
        "checks": checks,
    }


@cache
def nonzero_output_control(point=0):
    ctx = constraints.context(((1, 0, 0), (0, 1, 0), (-1, -1, 0)), symbolic_scale=True)
    result = construct(ctx, *physical.mixed_fields(ctx), point=point)
    return {
        "checks": result["checks"],
        "nonzero_output_coefficients": {
            mask: {
                key: sp.factor(result[key].coefficient(mask))
                for key in ("rho", "pressure", "n2", "missing_if_lapse_truncated")
            }
            for mask in (3, 5, 6)
        },
    }
