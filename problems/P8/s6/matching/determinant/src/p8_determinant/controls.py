"""Actual proportional solutions, a non-solution branch control, and auxiliary map."""

from functools import cache

import sympy as sp

from . import potential


@cache
def proportional_stiff():
    t = sp.Symbol("T", positive=True)
    ratios, gs, beta = (sp.Integer(1), sp.Integer(2), sp.Integer(3)), (1, 2, 3), (1, -1, 1)
    scales = tuple(y*t**sp.Rational(1, 3) for y in ratios)
    lapses = ratios
    d = potential.derive_three()
    point = {d["lambda"]: 1} | dict(zip(d["a"], scales, strict=True))
    point |= dict(zip(d["n"], lapses, strict=True)) | dict(zip(d["beta"], beta, strict=True))
    interaction_rho = tuple(sp.simplify(r.subs(point)) for r in d["rho"])
    interaction_p = tuple(sp.simplify(p.subs(point)) for p in d["pressure"])
    cosmo = tuple(-rho/g for rho, g in zip(interaction_rho, gs, strict=True))
    fields = tuple(sp.sqrt(sp.Rational(2*g, 3))*sp.log(t) for g in gs)
    hs = tuple(sp.diff(a, t)/(n*a) for a, n in zip(scales, lapses, strict=True))
    matter_rho = tuple(sp.diff(phi, t)**2/(2*n**2) for phi, n in zip(fields, lapses, strict=True))
    residuals = {}
    for i, (a, n, g, h, cc, rho, ri, pi, phi) in enumerate(zip(
            scales, lapses, gs, hs, cosmo, matter_rho, interaction_rho, interaction_p, fields, strict=True)):
        residuals[f"Einstein_00_{i}"] = sp.simplify(g*(3*h**2-cc)-rho-ri)
        residuals[f"Einstein_space_{i}"] = sp.simplify(g*(-2*sp.diff(h, t)/n-3*h**2+cc)-rho-pi)
        residuals[f"canonical_scalar_{i}"] = sp.simplify(sp.diff(a**3*sp.diff(phi, t)/n, t))
    k = sum(g*y**2 for g, y in zip(gs, ratios, strict=True))
    return {"T": t, "ratios": ratios, "Gs": gs, "beta": beta, "lambda": 1,
            "a": scales, "n": lapses, "fields": fields, "Hs": hs, "matter_rho": matter_rho,
            "matter_pressure": matter_rho, "interaction_rho": interaction_rho,
            "interaction_pressure": interaction_p, "cosmological_constants": cosmo,
            "S_y": sum(b*y for b, y in zip(beta, ratios, strict=True)), "K": k,
            "scaled_H_prime": sp.diff(hs[0]/sp.sqrt(k), t), "residuals": residuals}


@cache
def proportional_de_sitter():
    t = sp.Symbol("T", real=True)
    h0 = sp.Symbol("H0", real=True)
    base = proportional_stiff()
    ys, gs = base["ratios"], base["Gs"]
    scales = tuple(y*sp.exp(h0*t) for y in ys)
    hs = tuple(sp.diff(a, t)/(n*a) for a, n in zip(scales, ys, strict=True))
    cc = tuple(3*h**2-rho/g for h, rho, g in zip(hs, base["interaction_rho"], gs, strict=True))
    residuals = {}
    for i, (g, h, cosmological, ri, pi) in enumerate(zip(
            gs, hs, cc, base["interaction_rho"], base["interaction_pressure"], strict=True)):
        residuals[f"vacuum_00_{i}"] = sp.simplify(g*(3*h**2-cosmological)-ri)
        residuals[f"vacuum_space_{i}"] = sp.simplify(g*(-3*h**2+cosmological)-pi)
    return {"T": t, "H0": h0, "a": scales, "n": ys, "Hs": hs,
            "cosmological_constants": cc, "matter_nulls": (0, 0, 0),
            "residuals": residuals, "scaled_H_prime": sp.S.Zero}


def singular_sum_bianchi_control():
    """Bianchi-only point: NOT an actual solution or a viable singular branch."""
    beta, aa, nn, vv = (1, 1, -1), (1, 2, 1), (1, 1, 2), (1, -1, 0)
    d = potential.derive_three()
    point = {d["lambda"]: 1} | dict(zip(d["a"], aa, strict=True))
    point |= dict(zip(d["n"], nn, strict=True)) | dict(zip(d["beta"], beta, strict=True))
    point |= dict(zip(d["velocities"], vv, strict=True))
    return {"beta": beta, "a": aa, "n": nn, "velocities": vv,
            "S_a": d["S_a"].subs(point), "S_N": d["S_N"].subs(point),
            "S_a_dot": d["S_a_dot"].subs(point),
            "Q_i": tuple(sp.Rational(v, n) for v, n in zip(vv, nn, strict=True)),
            "bianchi": tuple(sp.factor(value.subs(point)) for value in d["bianchi"]),
            "actual_solution": False, "branch_health_claim": False}


@cache
def auxiliary_map():
    lam = sp.Symbol("lambda", positive=True)
    ww = sp.symbols("w0:4", nonzero=True)
    uu = sp.symbols("U0:4", nonzero=True)
    lag = lam*sp.prod(ww)*(3-sum(u/w for u, w in zip(uu, ww, strict=True)))
    stationary = dict(zip(ww, uu, strict=True))
    return {"lambda": lam, "w": ww, "U": uu, "L_aux": lag,
            "B": -3*lam/2, "p_over_beta": lam/2,
            "stationary_Euler": tuple(sp.factor(sp.diff(lag, w).subs(stationary)) for w in ww),
            "reduced_action": sp.factor(lag.subs(stationary)),
            "expected_reduced_action": -lam*sp.prod(uu),
            "matter_stays_on_original_EH_leaves": True,
            "matter_on_auxiliary_w": False}


def source_relocation_control():
    """Putting a canonical source on w changes its actual stationarity equation."""
    d = auxiliary_map()
    ww = d["w"]
    # Unit scalar coordinate velocity, V=0, at a unit auxiliary coframe.
    matter = sp.prod(ww)/(2*ww[0]**2)
    identity = dict.fromkeys(ww, sp.S.One)
    return tuple(sp.diff(matter, w).subs(identity) for w in ww)


def checks():
    stiff, vacuum, aux, singular = proportional_stiff(), proportional_de_sitter(), auxiliary_map(), singular_sum_bianchi_control()
    return dict(stiff["residuals"]) | dict(vacuum["residuals"]) | {
        "stiff_K36": stiff["K"]-36,
        "stiff_scaled_H_prime": sp.simplify(stiff["scaled_H_prime"]+1/(18*stiff["T"]**2)),
        "auxiliary_stationary_value": sp.expand(aux["reduced_action"]-aux["expected_reduced_action"]),
        "moved_canonical_source_changes_auxiliary_lapse": source_relocation_control()[0]+sp.Rational(1, 2),
        **{f"auxiliary_diagonal_Euler_{i}": value for i, value in enumerate(aux["stationary_Euler"])},
        **{f"singular_sum_Bianchi_only_{i}": value for i, value in enumerate(singular["bianchi"])},
    }
