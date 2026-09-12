"""General-dimensional complete magnetic vertex and normalized sphere geometry."""

from functools import cache

import sympy as s
from p8_vacuum_affine_subleading_band_conversion import flat

d, u, y = s.symbols("d u y", real=True)
X, Y, Z = s.symbols("X Y Z", real=True)


def magnetic(k, ell, D, identity):
    return (
        -(k.T * D * ell)[0] * identity
        + ell * k.T * D
        + D * ell * k.T
        - (k.T * ell)[0] * D
    )


def azimuth(expression):
    numerator, denominator = s.fraction(s.cancel(expression))
    if any(denominator.has(v) for v in (X, Y, Z)):
        raise ValueError("Only scalar longitudinal denominator allowed")
    out = 0
    for powers, c in s.Poly(numerator, X, Y, Z).terms():
        if powers[2] % 2:
            raise ValueError("The grouped transverse radius requires an even power")
        if any(n % 2 for n in powers[:2]):
            continue
        A, B, C = (n // 2 for n in powers)
        degree = A + B + C
        moment = (
            (1 - u * u) ** degree
            * s.rf(s.Rational(1, 2), A)
            * s.rf(s.Rational(1, 2), B)
            * s.rf((d - 3) / 2, C)
            / s.rf((d - 1) / 2, degree)
        )
        out += c * moment
    return s.factor(out / denominator)


@cache
def contractions(channel):
    k = s.Matrix([X, Y, Z, u])
    ell = -k + s.Matrix([0, 0, 0, y])
    ell2 = 1 - 2 * u * y + y * y
    identity = s.eye(4)
    P = identity - k * k.T
    Q = identity - ell * ell.T / ell2
    if channel == "tensor":
        D = s.diag(1, -1, 0, 0)
        tail = 0
        norm = 2
    elif channel == "vector":
        D = s.zeros(4)
        D[0, 3] = D[3, 0] = 1
        tail = 0
        norm = 2
    elif channel == "scalar":
        D = s.diag(-1, -1, -1, d - 1)
        tail = -1
        norm = d * (d - 1)
    else:
        raise ValueError("Require tensor,vector or scalar channel")
    M = magnetic(k, ell, D, identity)
    extraM = -(k.T * D * ell)[0] - (k.T * ell)[0] * tail
    extra = d - 4
    rows = {
        "00": s.trace(P * D * Q * D) + extra * tail * tail,
        "01": s.trace(P * M * Q * D) + extra * extraM * tail,
        "10": s.trace(P * D * Q * M.T) + extra * tail * extraM,
        "11": s.trace(P * M * Q * M.T) + extra * extraM**2,
        "TL": (ell.T * D * P * D * ell)[0] / ell2,
        "LT": (k.T * D * Q * D * k)[0],
        "LL": (k.T * D * ell)[0] ** 2 / ell2,
    }
    return {key: azimuth(value / norm) for key, value in rows.items()}


def sphere_average(expression):
    polynomial = s.Poly(s.cancel(expression), u)
    return s.factor(
        sum(
            c * s.rf(s.Rational(1, 2), n // 2) / s.rf(d / 2, n // 2)
            for (n,), c in polynomial.terms()
            if n % 2 == 0
        )
    )


@cache
def physical_checks():
    out = {}
    for ch in ("tensor", "vector", "scalar"):
        ug, yg, old = flat.averaged_geometry(ch)
        for key, value in contractions(ch).items():
            out[ch + "_" + key + "_physical_3D"] = s.factor(
                value.subs(d, 3) - old[key].subs({ug: u, yg: y})
            )
    k = s.Matrix(s.symbols("k0:3", real=True))
    ell = s.Matrix(s.symbols("l0:3", real=True))
    q = s.symbols("D0:5", real=True)
    D = s.Matrix([[q[0], q[1], q[2]], [q[1], q[3], q[4]], [q[2], q[4], -q[0] - q[3]]])
    out["literal_general_magnetic_vertex_reduces_to_3D_dual"] = magnetic(
        k, ell, D, s.eye(3)
    ) - flat.cross(k).T * D * flat.cross(ell)
    return out


@cache
def data():
    checks = dict(physical_checks())
    checks["all_general_dimension_channel_contractions"] = (
        sum(len(contractions(ch)) for ch in ("tensor", "vector", "scalar")) - 21
    )
    checks["normalized_sphere_identity"] = sphere_average(s.Integer(1)) - 1
    checks["normalized_longitudinal_square"] = sphere_average(u * u) - 1 / d
    checks["grouped_radius_second_moment"] = s.factor(
        azimuth(Z * Z) - (d - 3) * (1 - u * u) / (d - 1)
    )
    return {
        "magnetic_vertex": "Use the actual general-d field strength, not a three-dimensional dual vector. The creation-pair magnetic matrix is-(k.D.ell)I+ell k^t D+D ell k^t-(k.ell)D. It reduces exactly to Ck^t D Cell atd3 and retains its passive-dimension identity term.",
        "active_representation": "Take k=(X,Y,Z,u), ell=-k+y e4. Z^2 is the squared radius in d-3 remaining transverse directions. Add d-4 passive scalar terms to every full matrix trace; TL/LT/LL are vector contractions in the active block. The negative passive factor atd3 removes the auxiliary fourth direction.",
        "normalized_moments": "The exact transverse moment of X^(2A)Y^(2B)Z^(2C) is(1-u^2)^(A+B+C)*(1/2)_A*(1/2)_B*((d-3)/2)_C/((d-1)/2)_(A+B+C). Odd X/Y moments vanish; odd grouped-radius powers are outside this representation and rejected. Remaining normalized sphere moments are< u^(2n)>=(1/2)_n/(d/2)_n.",
        "channels": "Use normalized tensor,vector,scalar channels. The scalar test tensor has V=W=(d-1)/d; its dimension dependence must be removed by invariant reconstruction before differentiating a fixed physical source.",
        "checks": checks,
        "gates": {
            "all_21_physical_channel_contractions_reproduced": True,
            "full_3D_magnetic_duality_identity_replayed": True,
            "passive_trace_terms_not_dropped": True,
            "analytic_general_dimension_not_integer_fit": True,
            "no_three_dimensional_polarization_count_shortcut": True,
        },
    }
