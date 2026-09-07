"""Fixed-physical-frame dictionary for a formal stationary-f truncation.

The curvature-square term is retained.  These identities do not assert a
controlled heavy-mode inverse, an actual low-energy gap, or a UV verdict.
"""

from fractions import Fraction
from functools import lru_cache

import sympy as sp


def _rational(value, name):
    if isinstance(value, (bool, float, sp.Float)):
        raise TypeError(f"{name} must be an exact rational")
    if not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError(f"{name} must be an exact rational")
    return sp.Rational(value)


def _positive(value, name):
    result = _rational(value, name)
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _relative_error(value, name):
    result = _rational(value, name)
    if not 0 <= result < 1:
        raise ValueError(f"{name} must lie in [0,1)")
    return result


@lru_cache(maxsize=1)
def coefficients():
    """Return the symbolic IBP normal form, with q=(log r)_,phi.

    ``F_X2_coefficient`` multiplies X**2; it is not F_XX (twice it).
    The display retains kappa*(Ricci**2-R**2/3) as a separate operator.
    """
    phi, x = sp.symbols("phi X", real=True)
    m2 = sp.Symbol("M_squared", positive=True)
    rr = sp.Function("r")(phi)
    qq = sp.Function("q")(phi)
    kk = sp.Function("kappa")(phi)
    b0 = sp.Function("beta0")(phi)
    b1 = sp.Function("beta1")(phi)
    q1, q2 = sp.diff(qq, phi), sp.diff(qq, phi, 2)
    k1, k2, k3 = (sp.diff(kk, phi, n) for n in (1, 2, 3))
    f2x = -2 * qq * k1
    a1 = -4 * qq * k1
    kx = 6 * ((qq**2 + q1) * k1 + qq * k2)
    fx2 = 2 * (k1 * (2 * qq**3 + 2 * qq * q1 + q2)
               + k2 * (qq**2 + 2 * q1) + qq * k3)
    return {
        "symbols": {"phi": phi, "X": x, "M_squared": m2,
                    "r": rr, "q": qq, "kappa": kk,
                    "beta0": b0, "beta1": b1},
        "q_definition": sp.diff(rr, phi) / rr,
        "kappa_definition": m2**2 * rr**3 / (4 * b1),
        "curvature_squared_coefficient": kk,
        "F2": -m2 * (1 + rr**2) / 2 + f2x * x,
        "F2_X": f2x,
        "A1": a1,
        "A2": -a1,
        "A3": sp.S.Zero,
        "A4": sp.S.Zero,
        "A5": sp.S.Zero,
        "K": kx * x,
        "K_X": kx,
        "F": (-2 * (b0 + 3 * rr * b1)
              + (sp.Rational(1, 2) - 3 * m2 * sp.diff(rr, phi)**2) * x
              + fx2 * x**2),
        "F_X2_coefficient": fx2,
        "matter": sp.Symbol("Y", real=True) / 2,
    }


def reparametrization():
    """Exact scalar-only clock weights; g and chi are not changed."""
    jj, xx = sp.symbols("J X_theta", positive=True)
    f2, f2x, a1, a3 = sp.symbols("F2 F2_X A1 A3", real=True)
    xp = jj**2 * xx
    gt = -2 * f2 + 2 * xp * a1
    xi = a1 - 2 * f2x
    return {
        "symbols": {"J": jj, "X_theta": xx, "F2": f2,
                    "F2_X": f2x, "A1": a1, "A3": a3},
        "X_parent": xp,
        "F2_X_theta": jj**2 * f2x,
        "A1_theta": jj**2 * a1,
        "A3_theta": jj**4 * a3,
        "Xi_parent": xi,
        "Xi_theta": jj**2 * xi,
        "GT_parent": gt,
        "GT_theta": -2 * f2 + 2 * xx * jj**2 * a1,
        "I_parent": xp * xi / gt,
        "I_theta": xx * jj**2 * xi / (-2 * f2 + 2 * xx * jj**2 * a1),
    }


def center(c=3, x=1, parent_planck_squared=1, time_scale=1):
    """Actual coefficient-label center, not a new constant-clock vacuum.

    M_star**2=5*M**2 grants a best fit of the leading tensor value at X=1.
    It does not change the parent's fixed free-chi background normalization.
    ``formal.GT`` is the scalar-tensor coefficient with curvature squared
    retained separately, not a kinetic response of the full higher-order ODE.
    """
    cc = _rational(c, "c")
    xx = _rational(x, "X_theta")
    m2 = _positive(parent_planck_squared, "parent_planck_squared")
    tau = _positive(time_scale, "time_scale")
    if not 2 < cc <= 4:
        raise ValueError("c must lie in (2,4]")
    if not sp.Rational(9, 10) < xx < sp.Rational(11, 10):
        raise ValueError("X_theta must lie in the open clock tube")
    kbar = 64 / cc - sp.Rational(801, 100)
    beta1 = m2 * 32 / (tau**2 * cc * (cc - 2))
    kap = m2 * tau**2 * cc * (cc - 2) / 16
    target_m2 = 5 * m2
    return {
        "c": cc, "delta": cc - 2, "X_theta": xx,
        "r": sp.Integer(2), "r_phi": sp.S.Zero,
        "kappa_phi": sp.S.Zero, "kbar": kbar,
        "clock_J_squared": m2 * kbar / tau**2,
        "beta1": beta1, "kappa": kap,
        "M_star_squared": target_m2,
        "formal": {"F2": -target_m2 / 2, "F2_X": sp.S.Zero,
                   "A1": sp.S.Zero, "A3": sp.S.Zero,
                   "GT": target_m2, "Xi": sp.S.Zero, "I": sp.S.Zero},
        "CD": {"F2": -target_m2 * xx / 2, "F2_X": -target_m2 / 2,
               "A1": sp.S.Zero, "A3": target_m2 / xx,
               "GT": target_m2 * xx, "Xi": target_m2, "I": sp.S.One},
        "normalized_F2_defect": (xx - 1) / 2,
        "normalized_F2_X_defect": sp.Rational(1, 2),
        "normalized_Xi_defect": -sp.S.One,
        "I_defect": -sp.S.One,
        "curvature_square_over_Mstar_tau2": cc * (cc - 2) / 80,
        "weyl_square_over_Mstar_tau2": cc * (cc - 2) / 160,
        "controlled_reduction_claimed": False,
        "matter_background_retuned": False,
    }


def remainder_floor(relative_xi_error=0, target_planck_squared=1):
    """Necessary weighted C1 coefficient remainder at theta=0, X_theta=1.

    For |Xi_eff/M_star**2-1|<=eta, the weighted l1 remainder
    |Delta A1|+2|Delta F2_X| is at least the returned floor.
    """
    eta = _relative_error(relative_xi_error, "relative_xi_error")
    m2 = _positive(target_planck_squared, "target_planck_squared")
    return {"weighted_C1_floor": m2 * (1 - eta),
            "normalized_floor": 1 - eta,
            "coefficient_norm": "abs(Delta_A1)+2*abs(Delta_F2_X)",
            "error_domain": "normalized Xi coefficient at theta=0, X_theta=1"}


def invariant_remainder_floor(kinetic_error=0, invariant_error=0,
                              target_planck_squared=1):
    """Necessary floor when both GT and I=X*Xi/GT are approximately matched."""
    et = _relative_error(kinetic_error, "kinetic_error")
    ei = _relative_error(invariant_error, "invariant_error")
    m2 = _positive(target_planck_squared, "target_planck_squared")
    factor = (1 - et) * (1 - ei)
    return {"weighted_C1_floor": m2 * factor, "normalized_floor": factor,
            "coefficient_norm": "abs(Delta_A1)+2*abs(Delta_F2_X)",
            "error_domain": "abs(GT/M_star**2-1)<=et, abs(I-1)<=ei, X_theta=1"}


def contacts():
    """Inverse-metric convention and the induced free-chi operators."""
    kap, planck = sp.symbols("kappa Q", positive=True)
    ric_vv, ric_scalar, yy = sp.symbols("Ric_chichi R Y", real=True)
    return {
        "symbols": {"kappa": kap, "Q": planck, "Ric_chichi": ric_vv,
                    "R": ric_scalar, "Y": yy},
        "inverse_metric_factor": 2 * kap / planck,
        "chi_curvature_contact": kap * (ric_vv - ric_scalar * yy / 3) / planck,
        "pure_chi_Y2_coefficient_mod_leading_metric_equations": 2 * kap / (3 * planck**2),
        "mixed_clock_contacts_retained": True,
        "same_physical_metric_after_redefinition": False,
    }


def constant_r_schur():
    """A source-preserving flat quadratic normalization regression only."""
    gg, ff, nu = sp.symbols("G F nu", positive=True)
    dd = sp.Symbol("D", real=True)
    total = gg + ff
    mass2 = nu * total / (gg * ff)
    kernel = gg * dd + nu - nu**2 / (ff * dd + nu)
    return {
        "symbols": {"G": gg, "F": ff, "nu": nu, "D": dd},
        "kernel": kernel,
        "response": 1 / (total * dd) + ff / (gg * total * (dd + mass2)),
        "D2_coefficient": -ff**2 / nu,
        "kappa": ff**2 / (2 * nu),
        "weyl_coefficient": ff**2 / (4 * nu),
        "source_tensor_contact_coefficient": ff**2 / (2 * nu * total**2),
        "pure_chi_Y2_coefficient": ff**2 / (3 * nu * total**2),
        "full_physical_observable_retained": True,
        "actual_VARIABLE_vacuum_claimed": False,
    }


def checks():
    """Exact algebraic replay; the general IBP proof is in notes/dictionary.md."""
    data = coefficients()
    z = data["symbols"]
    phi, q, kap = z["phi"], z["q"], z["kappa"]
    qp = sp.diff(q, phi)
    a, b, c, d = -2*q, 2*(q**2-qp), 2*q, 2*qp+q**2
    j = 4*q*(sp.diff(kap, phi)+kap*q)
    residues = {
        "B_square_Hessian": a*a-4*q*q,
        "B_square_box": 2*a*c+4*c*c-(a+4*c)**2/3+4*q*q,
        "B_square_vHv": 2*a*b+8*q*(q*q-qp),
        "B_square_Xbox": (2*a*d+2*b*c+8*c*d
                          -sp.Rational(2, 3)*(a+4*c)*(b+4*d)+4*q*(q*q+2*qp)),
        "B_square_X2": b*b+2*b*d+4*d*d-(b+4*d)**2/3+12*q*q*qp,
        "A1_after_weighted_Ricci_IBP": 4*kap*q*q-j-data["A1"],
        "quartic_Horndeski_locus": data["A1"]-2*data["F2_X"],
        "K_after_IBP": -12*kap*q*qp+sp.Rational(3, 2)*sp.diff(j, phi)-data["K_X"],
        "F_X2_after_IBP": (4*sp.diff(kap, phi)*q*(q*q-qp)
                           -4*kap*(qp**2+q*sp.diff(qp, phi))
                           +sp.diff(j, phi, 2)/2-data["F_X2_coefficient"]),
    }
    clock = reparametrization()
    residues["clock_invariant_I"] = clock["I_theta"]-clock["I_parent"]
    residues["clock_tensor_coefficient"] = clock["GT_theta"]-clock["GT_parent"]
    schur = constant_r_schur()
    ss = schur["symbols"]
    residues["source_schur_inverse"] = schur["kernel"]*schur["response"]-1
    residues["source_schur_D2"] = (sp.diff(schur["kernel"], ss["D"], 2).subs(ss["D"], 0)/2
                                   -schur["D2_coefficient"])
    residues["source_contact_Y2"] = (sp.Rational(2, 3)*schur["kappa"]/(ss["G"]+ss["F"])**2
                                     -schur["pure_chi_Y2_coefficient"])
    # All ten Ricci components and all four scalar-gradient components.
    eta = sp.diag(1, -1, -1, -1)
    vals = sp.symbols("R00 R01 R02 R03 R11 R12 R13 R22 R23 R33", real=True)
    ric = sp.zeros(4)
    positions = tuple((row, col) for row in range(4) for col in range(row, 4))
    for value, (i, j) in zip(vals, positions, strict=True):
        ric[i, j] = ric[j, i] = value
    vv = sp.Matrix(sp.symbols("v0:4", real=True))
    scalar = sp.trace(eta*ric)
    yy = (vv.T*eta*vv)[0]
    tcov = vv*vv.T-eta*yy/2
    kap0, planck = sp.symbols("kap0 Planck", positive=True)
    delta_inverse = 2*kap0/planck*(eta*ric*eta-scalar*eta/6)
    contraction = lambda left, right: sum(left[i, j]*right[i, j]
                                          for i in range(4) for j in range(4))
    residues["inverse_metric_chi_sign"] = (
        contraction(tcov, delta_inverse)/2
        -kap0/planck*((vv.T*eta*ric*eta*vv)[0]-scalar*yy/3))
    residues["inverse_metric_curvature_cancel"] = (
        -planck*contraction(ric-eta*scalar/2, delta_inverse)/2
        +kap0*(sp.trace(eta*ric*eta*ric)-scalar**2/3))
    residues["canonical_chi_stress_square"] = sp.trace(eta*tcov*eta*tcov)-yy**2
    residues["canonical_chi_stress_trace"] = sp.trace(eta*tcov)+yy
    for cc, xx in ((3, 1), (sp.Rational(5, 2), sp.Rational(19, 20)),
                   (4, sp.Rational(21, 20))):
        fixture = center(cc, xx)
        tag = f"{cc}_{xx}"
        formal, cd = fixture["formal"], fixture["CD"]
        residues[f"center_Xi_{tag}"] = cd["A1"]-2*cd["F2_X"]-cd["Xi"]
        residues[f"center_I_{tag}"] = xx*cd["Xi"]/cd["GT"]-cd["I"]
        residues[f"center_defect_{tag}"] = ((formal["Xi"]-cd["Xi"])
                                             /fixture["M_star_squared"]+1)
    return {key: sp.factor(sp.expand(value)) for key, value in residues.items()}


def calibration():
    return {
        "center": center(),
        "exact_remainder": remainder_floor(),
        "ten_percent_invariant_errors": invariant_remainder_floor(
            sp.Rational(1, 10), sp.Rational(1, 10)),
        "whole_tube_normalized_I_defect": -sp.S.One,
        "center_normalized_F2_X_remainder_required": sp.Rational(1, 2),
        "center_A3_remainder_infimum_over_open_tube": sp.Rational(10, 11),
        "curvature_square_retained": True,
        "GT_is_scalar_tensor_part_only": True,
        "scalar_only_clock_change": True,
        "order_six_can_contribute": True,
        "controlled_heavy_EFT_claimed": False,
        "UV_or_cutoff_claimed": False,
    }
