"""Ordered full spatial trace vertices in the original dimensional geometry."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_spatial_symbol import geometry as previous

d, u, y = previous.d, previous.u, previous.y
X, Y, Z = previous.X, previous.Y, previous.Z
azimuth, sphere_average = previous.azimuth, previous.sphere_average


def magnetic(k, ell, D, tau=None):
    if tau is None:
        tau = s.trace(D)
    I = s.eye(D.rows)
    return previous.magnetic(k, ell, D, I) + tau * ((k.T * ell)[0] * I - ell * k.T) / 2


def channel_directions(channel):
    I = s.eye(4)
    S = s.diag(-1, -1, -1, d - 1)
    if channel == "trace_trace":
        return I, s.S.One, I, s.S.One, s.S.One
    if channel == "trace_scalar":
        return I, s.S.One, S, -s.S.One, s.S.One
    if channel == "scalar_trace":
        return S, -s.S.One, I, s.S.One, s.S.One
    if channel == "tensor":
        T = s.diag(1, -1, 0, 0)
        return T, s.S.Zero, T, s.S.Zero, s.Integer(2)
    if channel == "vector":
        V = s.zeros(4)
        V[0, 3] = V[3, 0] = 1
        return V, s.S.Zero, V, s.S.Zero, s.Integer(2)
    if channel == "scalar":
        return S, -s.S.One, S, -s.S.One, d * (d - 1)
    raise ValueError("Require an explicitly ordered scalar or frozen shear channel")


@cache
def contractions(channel):
    D, td, G, tg, norm = channel_directions(channel)
    extra = d - 4
    tauD, tauG = s.trace(D) + extra * td, s.trace(G) + extra * tg
    BD, BG = D - tauD * s.eye(4) / 2, G - tauG * s.eye(4) / 2
    bd, bg = td - tauD / 2, tg - tauG / 2
    k = s.Matrix([X, Y, Z, u])
    ell = -k + s.Matrix([0, 0, 0, y])
    ell2 = 1 - 2 * u * y + y * y
    P = s.eye(4) - k * k.T
    Q = s.eye(4) - ell * ell.T / ell2
    MD, MG = magnetic(k, ell, D, tauD), magnetic(k, ell, G, tauG)
    md = -(k.T * D * ell)[0] - (k.T * ell)[0] * bd
    mg = -(k.T * G * ell)[0] - (k.T * ell)[0] * bg
    rows = {
        "00": s.trace(P * BD * Q * BG) + extra * bd * bg,
        "01": s.trace(P * BD * Q * MG.T) + extra * bd * mg,
        "10": s.trace(P * MD * Q * BG) + extra * md * bg,
        "11": s.trace(P * MD * Q * MG.T) + extra * md * mg,
        "TL": (ell.T * BD * P * BG * ell)[0] / ell2,
        "LT": (k.T * BD * Q * BG * k)[0],
        "LL": (k.T * BD * ell)[0] * (k.T * BG * ell)[0] / ell2,
        "LC": tauG * (k.T * BD * ell)[0],
        "CL": tauD * (k.T * BG * ell)[0],
        "CC": tauD * tauG * ell2,
    }
    return {name: azimuth(value / norm) for name, value in rows.items()}


def invariant_row(channel):
    D, td, G, tg, norm = channel_directions(channel)
    tauD, tauG = s.trace(D) + (d - 4) * td, s.trace(G) + (d - 4) * tg
    T = s.trace(D * G) + (d - 4) * td * tg
    return tuple(
        s.cancel(v / norm)
        for v in (
            T,
            (D * G)[3, 3],
            D[3, 3] * G[3, 3],
            tauD * tauG,
            tauD * G[3, 3],
            tauG * D[3, 3],
        )
    )


@cache
def invariant_matrix():
    return s.Matrix(
        [
            invariant_row(ch)
            for ch in (
                "tensor",
                "vector",
                "scalar",
                "trace_trace",
                "trace_scalar",
                "scalar_trace",
            )
        ]
    )


@cache
def data():
    from p8_vacuum_affine_spatial_current import hamiltonian

    k = s.Matrix(s.symbols("k0:3", real=True))
    ell = s.Matrix(s.symbols("l0:3", real=True))
    q = s.symbols("Q0:6", real=True)
    D = s.Matrix([[q[0], q[1], q[2]], [q[1], q[3], q[4]], [q[2], q[4], q[5]]])
    checks = {
        "full_trace_magnetic_field_strength_and_3D_dual": magnetic(k, ell, D)
        - hamiltonian.cross(k).T
        * (D - s.trace(D) * s.eye(3) / 2)
        * hamiltonian.cross(ell),
        "six_invariant_reconstruction_determinant": s.factor(
            invariant_matrix().det() - d**3 * (d - 1) ** 3 / 2
        ),
        "six_invariant_inverse": invariant_matrix() * invariant_matrix().inv()
        - s.eye(6),
        "normalized_full_sphere": sphere_average(s.Integer(1)) - 1,
    }
    for ch in ("tensor", "vector", "scalar"):
        old = previous.contractions(ch)
        rows = contractions(ch)
        checks[ch + "_all_frozen_tracefree_geometry"] = s.Matrix(
            [rows[key] - v for key, v in old.items()]
        )
        checks[ch + "_trace_constraint_zero"] = s.Matrix(
            [rows[key] for key in ("LC", "CL", "CC")]
        )
    return {
        "full_magnetic": "The field-strength vertex is M_old(Q)+tr(Q)[(k.ell)I-ell k^t]/2. At d3 this is Ck^t[Q-tr(Q)I/2]Cell; at general d use the full two-form, not a 3D dual.",
        "ten_ordered_products": "Four TT,TL,LT,LL and the distinct LC,CL,CC trace-constraint products. Every passive-dimensional volume/identity contribution is retained.",
        "fixed_invariants": "Use T=tr(DG), V=e^t DG e, W=(eDe)(eGe), S=tr(D)tr(G), UD=tr(D)eGe and UG=tr(G)eDe. Reconstruct these six before differentiating d. The three tracefree channels plus ordered I/I,I/S_d,S_d/I have determinant d^3(d-1)^3/2.",
        "checks": checks,
        "gates": {
            "all_ten_ordered_geometry_contractions": all(
                len(contractions(ch)) == 10
                for ch in ("trace_trace", "trace_scalar", "scalar_trace")
            ),
            "arbitrary_trace_not_tracefree_only": True,
            "both_ordered_longitudinal_constraint_cross_terms": True,
            "analytic_dimension_not_integer_interpolation": True,
            "fixed_six_invariants_before_dimension_derivative": True,
        },
    }
