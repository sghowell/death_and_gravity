"""Literal trace-curvature spectator and named trace-form kinetic screen.

Use the frozen source signature (-+++), derivative index last, and an
unrestricted connection. A sign of an AUXILIARY mass matrix is not a
kinetic verdict; the two physical time-kinetic sectors are built here.
"""
from fractions import Fraction
from functools import cache
from itertools import permutations

import sympy as sp
from p8_affine import connection as parent


def nonzero_rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Use a literal exact rational kinetic coefficient")
    result = sp.Rational(value)
    if result == 0:
        raise ValueError("Zero is the auxiliary, not the kinetic, rank chart")
    return result


@cache
def homothetic():
    """Trace commutator cancellation, not a Levi-Civita Ricci substitution."""
    left = sp.Matrix(4, 4, sp.symbols("left_0:16", real=True))
    right = sp.Matrix(4, 4, sp.symbols("right_0:16", real=True))
    trace_commutator = sp.expand(sp.trace(left*right-right*left))
    a_mu_nu, a_nu_mu, mixed_log, zeta = sp.symbols(
        "d_mu_A_nu d_nu_A_mu mixed_log zeta", real=True)
    trace_curvature = (mixed_log+4*a_mu_nu)-(mixed_log+4*a_nu_mu)+trace_commutator
    field_strength = a_mu_nu-a_nu_mu
    f_squared = sp.Symbol("F_squared", real=True)
    return {
        "trace_commutator": trace_commutator,
        "homothetic_curl": sp.expand(trace_curvature-4*field_strength),
        "Maxwell_normalization": sp.expand(-zeta*(16*f_squared)/64+zeta*f_squared/4),
    }


@cache
def sectors():
    """Two actual 64-component vectors in unique homogeneous SO(3) sectors.

    The harmonic polynomial x^3-3*x*y^2 supplies an STF spin-3 component;
    epsilon_ijk supplies the axial pseudoscalar. All temporal components
    vanish. The full time-kinetic matrix includes the trace-free projector:
    d²[sum_i (tr(dot C_i²)-(tr dot C_i)²/4)/2]/d dot C².
    """
    axial, stf = sp.zeros(64, 1), sp.zeros(64, 1)
    for triple in permutations(range(1, 4)):
        inversions = sum(triple[i] > triple[j] for i in range(3) for j in range(i+1, 3))
        axial[parent.index(*triple)] = (-1)**inversions
    stf[parent.index(1, 1, 1)] = 1
    for triple in set(permutations((1, 2, 2))):
        stf[parent.index(*triple)] = -1
    kinetic = sp.zeros(64)
    for a in range(4):
        for b in range(4):
            for i in range(1, 4):
                kinetic[parent.index(a, b, i), parent.index(b, a, i)] = 1
                kinetic[parent.index(a, a, i), parent.index(b, b, i)] -= sp.Rational(1, 4)
    full = parent.quadratic()
    return {
        "axial": sp.ImmutableMatrix(axial), "stf": sp.ImmutableMatrix(stf),
        "kinetic": sp.ImmutableMatrix(kinetic), "mass": full["hessian"],
        "source": full["source"], "projective": full["gauge"],
    }


@cache
def checks():
    data = sectors()
    out = homothetic().copy()
    out["full_tracefree_kinetic_projective_kernel"] = data["kinetic"]*data["projective"]
    for name, sign, norm in (("axial", -1, 6), ("stf", 1, 4)):
        vector = data[name]
        out[name+"_kinetic_eigen_equation"] = data["kinetic"]*vector-sign*vector
        out[name+"_full_mass_eigen_equation"] = (data["mass"]*vector+4*parent.P**2*vector).applyfunc(sp.factor)
        out[name+"_no_scalar_Hessian_source"] = sp.factor((data["source"].T*vector)[0])
        out[name+"_projective_trace"] = data["projective"].T*vector
        out[name+"_norm"] = (vector.T*vector)[0]-norm
        out[name+"_temporal_multiplier_components"] = sp.Matrix([
            vector[parent.index(a, b, c)] for a in range(4) for b in range(4)
            for c in range(4) if 0 in (a, b, c)])
    out["orthogonal_sectors"] = (data["stf"].T*data["axial"])[0]
    # STF trace vanishes on any pair, not only the projective trace.
    out["stf_all_pair_traces"] = sp.Matrix([
        sum(data["stf"][parent.index(i, i, k)] for i in range(1, 4))
        for k in range(1, 4)])
    return out


def sign_control(zeta):
    value = nonzero_rational(zeta)
    return {"zeta": value, "stf_velocity_squared_coefficient": 2*value,
            "axial_velocity_squared_coefficient": -3*value,
            "one_strictly_negative": bool(2*value < 0 or -3*value < 0),
            "coefficient_product": -6*value**2}
