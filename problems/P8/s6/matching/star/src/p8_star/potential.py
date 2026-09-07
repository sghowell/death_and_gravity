"""Literal lapse/scale variations of -2 sqrt|g_i| sum beta_n e_n(S_i)."""

from functools import cache

import sympy as sp

from .exact import coefficients, positive


@cache
def derive():
    r, n, c = sp.symbols("R N c", positive=True)
    ni, nu, ai, au = sp.symbols("n_i n_u a_i a_u", positive=True)
    beta = sp.symbols("beta0:5", real=True)
    hi, hu = sp.symbols("H_i H_u", real=True)
    a = beta[0]+3*beta[1]*r+3*beta[2]*r**2+beta[3]*r**3
    b = beta[1]+3*beta[2]*r+3*beta[3]*r**2+beta[4]*r**3
    j = beta[1]+2*beta[2]*r+beta[3]*r**2
    # Positive FLRW square-root eigenvalues (N,R,R,R), in leaf-first orientation.
    elementary = (1, n+3*r, 3*n*r+3*r**2, 3*n*r**2+r**3, n*r**3)
    eigenvalue_polynomial = sum(bn*en for bn, en in zip(beta, elementary, strict=True))
    lag = -2*ai**3*(ni*a.subs(r, au/ai)+nu*b.subs(r, au/ai))
    point = {au: r*ai, nu: n*ni}
    rho_i = sp.factor((-sp.diff(lag, ni)/ai**3).subs(point))
    pressure_i = sp.factor((sp.diff(lag, ai)/(3*ni*ai**2)).subs(point))
    rho_u = sp.factor((-sp.diff(lag, nu)/au**3).subs(point))
    pressure_u = sp.factor((sp.diff(lag, au)/(3*nu*au**2)).subs(point))
    null_i = sp.factor(rho_i+pressure_i)
    null_u = sp.factor(rho_u+pressure_u)
    r_dot_tau_i = r*(n*hu-hi)
    bianchi = sp.factor(sp.diff(rho_i, r)*r_dot_tau_i+3*hi*null_i)
    return {"R": r, "N": n, "c": c, "beta": beta,
            "n_i": ni, "n_u": nu, "a_i": ai, "a_u": au,
            "H_i": hi, "H_u": hu, "A": a, "B": b, "J": j,
            "elementary": elementary, "eigenvalue_polynomial": eigenvalue_polynomial,
            "L": lag, "rho_i": rho_i, "pressure_i": pressure_i,
            "rho_u": rho_u, "pressure_u": pressure_u,
            "null_i": null_i, "null_u": null_u,
            "R_dot_tau_i": r_dot_tau_i, "bianchi": bianchi}


def evaluate(beta, R, N):
    """Evaluate one link; N=n_u/n_i, not the inverse lapse ratio or speed c."""
    beta = coefficients(beta)
    r, n = positive(R, "R"), positive(N, "N")
    d = derive()
    point = dict(zip(d["beta"], beta, strict=True)) | {d["R"]: r, d["N"]: n}
    return {key: sp.factor(d[key].subs(point)) for key in
            ("A", "B", "J", "rho_i", "pressure_i", "rho_u", "pressure_u", "null_i", "null_u")}


def checks():
    d = derive()
    r, n, j, a, b = (d[key] for key in ("R", "N", "J", "A", "B"))
    return {
        "literal_eigenvalue_polynomial": sp.expand(d["eigenvalue_polynomial"]-a-n*b),
        "A_derivative": sp.expand(sp.diff(a, r)-3*j),
        "reciprocal_derivative": sp.expand(3*b-r*sp.diff(b, r)-3*j),
        "leaf_density": sp.expand(d["rho_i"]-2*a),
        "leaf_pressure": sp.expand(d["pressure_i"]+2*(a+(n-r)*j)),
        "central_density": sp.cancel(d["rho_u"]-2*b/r**3),
        "central_pressure": sp.cancel(d["pressure_u"]+2*j/(n*r**2)+2*sp.diff(b, r)/(3*r**2)),
        "leaf_null": sp.expand(d["null_i"]-2*(r-n)*j),
        "central_null": sp.cancel(d["null_u"]-2*j*(1-r/n)/r**3),
        "unfactored_bianchi": sp.factor(d["bianchi"]-6*n*j*(r*d["H_u"]-d["H_i"])),
        "paired_null_cancellation": sp.cancel(d["null_u"]+d["null_i"]/(n*r**3)),
    }
