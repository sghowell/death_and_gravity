"""The unchanged analytic target's exact clock jets and leading Euler sources."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family
from p8_exceptional_vacuum import analytic

from . import homogeneous


@cache
def data():
    old = family.data()
    u, X = family.u, family.X
    tree = family.original.data()["original_retuned_tree_scalar"]
    F0 = s.factor(old["F"].subs(X, 1))
    FX = s.factor(s.diff(old["F"], X).subs(X, 1))
    D = 1 + u * u
    h = D**3
    H = 4 * u / D
    Hd = s.diff(H, u)
    chidot = 1 / (10 * D**6)
    Y = chidot * chidot
    rhoR = 24 * (1 + 7 * u * u) / D**5
    rhoA = -12 * (1 + 5 * u * u) / D**5
    d = homogeneous.data()
    z = d["symbols"]
    t = z["t"]
    replace = {
        z["H"]: H,
        z["Hdot"]: Hd,
        z["h"]: h,
        z["F0"]: F0,
        z["FX"]: FX,
        z["A3"]: 1 / h,
        s.diff(z["A3"], t): s.diff(1 / h, u),
        z["chidot"]: chidot,
    }
    pieces = {
        name: {key: s.factor(value.xreplace(replace)) for key, value in row.items()}
        for name, row in d["complete_piecewise_first_metric_variations"].items()
    }
    checks = {
        "full_analytic_scalar_clock_value": s.factor(F0 - tree.subs(X, 1)),
        "full_analytic_scalar_first_clock_jet": s.factor(
            FX - s.diff(tree, X).subs(X, 1)
        ),
        "full_analytic_curvature_clock_value": s.factor(
            old["F2"].subs(X, 1) + s.Rational(1, 2)
        ),
        "full_analytic_curvature_first_clock_jet": s.factor(
            s.diff(old["F2"], X).subs(X, 1) + 1 / (2 * h)
        ),
        "full_analytic_A3_clock_value": s.factor(old["A3"].subs(X, 1) - 1 / h),
        "full_analytic_A4_clock_value": s.factor(
            old["A4"].subs(X, 1) + 1 / h + 7 / (4 * h * h)
        ),
        "full_analytic_A5_clock_value": s.factor(old["A5"].subs(X, 1) - 1 / (h * h)),
        "exact_CD_Hubble": s.factor(s.diff(D * D, u) / (D * D) - H),
        "original_M1_conserved_momentum": s.diff(D**6 * chidot, u),
        "original_M1_momentum_value": D**6 * chidot - s.Rational(1, 10),
        "target_scalar_value_from_pressure_equation": s.factor(
            F0 + 2 * Hd + 3 * H * H + Y / 2
        ),
        "target_scalar_first_jet_from_null_equation": s.factor(
            2 * FX + rhoR + rhoA + Y + 2 * Hd
        ),
        "independent_curvature_source": s.factor(
            pieces["curvature_nonEinstein"]["rho_Euler"] - rhoR
        ),
        "independent_A3_source": s.factor(pieces["A3"]["rho_Euler"] - rhoA),
        "entire_time_energy_equation": s.factor(
            sum(row["rho_Euler"] for row in pieces.values())
        ),
        "entire_time_pressure_equation": s.factor(
            sum(row["P_Euler"] for row in pieces.values())
        ),
        "entire_time_null_equation": s.factor(
            sum(row["null_Euler"] for row in pieces.values())
        ),
        "bounce_F0": F0.subs(u, 0) + s.Rational(1601, 200),
        "bounce_FX": FX.subs(u, 0) + s.Rational(2001, 200),
        "bounce_geometric_null": rhoR.subs(u, 0) + rhoA.subs(u, 0) - 12,
    }
    table = {name: dict(row) for name, row in pieces.items() if name != "scalar_F"}
    table["canonical_clock"] = {
        "rho_Euler": s.Rational(1, 2),
        "P_Euler": s.Rational(1, 2),
        "null_Euler": s.Integer(1),
    }
    table["scalar_F_minus_canonical"] = {
        key: s.factor(pieces["scalar_F"][key] - table["canonical_clock"][key])
        for key in table["canonical_clock"]
    }
    bounce = {
        name: {key: value.subs(u, 0) for key, value in row.items()}
        for name, row in table.items()
    }
    return {
        "normalized_by_physical_kappa": analytic.KAPPA,
        "actual_CD_scale_factor": D * D,
        "actual_CD_Hubble": H,
        "actual_original_M1_velocity": chidot,
        "full_analytic_target_clock_F0": F0,
        "full_analytic_target_clock_FX": FX,
        "complete_first_metric_variation_table": table,
        "actual_bounce_table": bounce,
        "source_ancestry": "Read the unchanged S6.109 full analytic family. Its original scalar is the actual CD_matter witness plus the frozen retuned margin times(X-1)^2/h^2. The rational/Gaussian switch and that margin leave the needed first clock jets unchanged; no constant-branch replacement is asserted globally.",
        "interpretation": "The extra curvature value and all A_i action values vanish on the unit clock, yet the curvature X derivative and the lapse derivative in A3 contribute to the first metric equations. The original rolling M1 scalar is retained. The canonical-clock reference is the kinetic X/2 term; the entire lower potential remains in F minus that reference. Values are dimensionless Euler sources; multiply by kappa for the physical action. Individual pieces are not independently conserved physical matter tensors.",
        "classical_scope": "This rechecks the leading target background equations and identifies their required parent matching data. It does not construct a healthy propagating UV parent, quantum-corrected bounce or response bound. The small S6.171 free stress remains an added source needing state-aware background control.",
        "checks": checks,
    }
