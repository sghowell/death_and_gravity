"""Literal 4.8 auxiliary elimination, physical clock, and rolling normalization."""

from functools import cache

import sympy as sp

from . import model as m


@cache
def spring():
    aa, bb, cc = sp.symbols("A B C", positive=True)
    x, z, u, v = sp.symbols("Edot Sdot shift_g shift_f", real=True)
    lagrangian = aa*(x-u)**2+bb*(z-v)**2+cc*(u-v)**2
    solution = sp.solve([sp.diff(lagrangian, u), sp.diff(lagrangian, v)], (u, v))
    reduced = m.cancel(lagrangian.subs(solution))
    harmonic = aa*bb*cc/(aa*bb+aa*cc+bb*cc)
    deviations = sp.Matrix([u-solution[u], v-solution[v]])
    matrix = sp.Matrix([[aa+cc, -cc], [-cc, bb+cc]])
    return {"A": aa, "B": bb, "C": cc, "x": x, "z": z, "u": u, "v": v,
            "L": lagrangian, "solution": solution, "reduced": reduced,
            "harmonic": harmonic, "matrix": matrix,
            "completed_square": (deviations.T*matrix*deviations)[0]}


@cache
def derive():
    k, xi = sp.symbols("k Xi", positive=True)
    mu = sp.Symbol("mu", real=True)
    d = spring()
    assignment = {d["A"]: k**2/2, d["B"]: m.F*m.y**3*k**2/(2*m.G*m.c),
                  d["C"]: m.a**2*xi/(2*m.G)}
    inverse_c = 1+m.G*m.c/(m.F*m.y**3)+m.G*k**2/(m.a**2*xi)
    cal_c = 1/inverse_c
    K = m.G*m.a**3*cal_c/(16*m.Ng)
    U = m.Ng*m.a**3*mu/16
    # dT=Ne dt. The physical-clock conversion must precede normalization.
    kt, ut = K*m.Ne, U/m.Ne
    bare_frequency = m.cancel(ut/kt)
    speed = m.cancel(sp.diff(bare_frequency, k, 2)*m.Ae**2/2)
    return {"k": k, "Xi": xi, "mu": mu, "assignment": assignment,
            "C": cal_c, "K": K, "U": U, "K_T": kt, "U_T": ut,
            "bare_frequency_squared": bare_frequency, "physical_principal_speed_squared": speed,
            "literal_K": m.cancel(m.G*m.a**3*d["harmonic"].subs(assignment)/(8*m.Ng*k**2)),
            "literal_U": m.Ng*m.a**3*mu/16}


def canonical_frequency(K, U, clock):
    """For real positive K: Q=sqrt(2K)q, Omega²=U/K-(sqrt K)''/sqrt K.

    This is a differential identity, not a lower bound or an eigenmode claim.
    Input K and U must be written in the same clock as `clock`.
    """
    root = sp.sqrt(K)
    return U/K-sp.diff(root, clock, 2)/root


def checks():
    d, q = spring(), derive()
    time = sp.Symbol("T", real=True)
    root = sp.Function("root_K", positive=True)(time)
    norm, normdot, stiff = sp.symbols("Q Qdot U", real=True)
    connection = sp.diff(root, time)/root
    unintegrated = ((normdot-connection*norm)**2-stiff*norm**2/root**2)/2
    frequency = canonical_frequency(root**2, stiff, time)
    integrated = (normdot**2-frequency*norm**2)/2
    boundary_derivative = -sp.diff(connection, time)*norm**2/2-connection*norm*normdot
    return {"stationary_action": m.cancel(d["reduced"]-d["harmonic"]*(d["x"]-d["z"])**2),
            "full_positive_completed_square": m.cancel(d["L"]-d["reduced"]-d["completed_square"]),
            "positive_constraint_determinant": m.cancel(d["matrix"].det()-(d["A"]*d["B"]+d["A"]*d["C"]+d["B"]*d["C"])),
            "literal_vector_inertia": m.cancel(q["literal_K"]-q["K"]),
            "literal_vector_stiffness": m.cancel(q["literal_U"]-q["U"]),
            "physical_clock_inertia": m.cancel(q["K_T"]-m.G*m.a**3*m.s*q["C"]/16),
            "physical_clock_stiffness": m.cancel(q["U_T"]-m.a**3*q["mu"]/(16*m.s)),
            "physical_principal_speed": m.cancel(q["physical_principal_speed_squared"]-(m.r/m.s)**2*q["mu"]/q["Xi"]),
            "canonical_time_boundary": sp.simplify(unintegrated-integrated-boundary_derivative)}


def source_normalization_checks():
    """Literal source 4.8/4.7/4.11 only; no unprovided rescaling of q."""
    q = derive()
    # The printed 4.12 coefficient in its stated q=k(E-S) convention.
    source_K = m.G*m.a**3*q["C"]/(8*m.Ng)
    source_U = m.Ng*m.a**3*q["mu"]*(m.F*m.y**2/(m.G+m.F*m.y**2))/8
    return {"4_11_C_matches_literal": m.cancel(q["C"]-1/(m.G*q["k"]**2/(m.a**2*q["Xi"])
                                                                        +m.G*m.c/(m.F*m.y**3)+1)),
            "printed_4_12_K_over_literal": m.cancel(source_K/q["K"]),
            "printed_4_12_U_over_literal": m.cancel(source_U/q["U"]),
            "printed_frequency_over_literal": m.cancel((source_U/source_K)/(q["U"]/q["K"]))}


def negative_controls():
    d = derive()
    time = sp.Symbol("T", real=True)
    rate = sp.Symbol("rate", positive=True)
    arbitrary_field_factor = m.F*m.y**2/(m.G+m.F*m.y**2)
    return {"omit_composite_clock_and_ruler": m.cancel(d["mu"]/d["Xi"]-d["physical_principal_speed_squared"]),
            "extra_source_mass_factor_changes_speed": m.cancel((arbitrary_field_factor-1)*d["physical_principal_speed_squared"]),
            "zero_mu_does_not_fix_canonical_frequency": sp.simplify(canonical_frequency(sp.exp(2*rate*time), 0, time))}
