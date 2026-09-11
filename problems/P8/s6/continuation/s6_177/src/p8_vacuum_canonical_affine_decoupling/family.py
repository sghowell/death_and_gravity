"""Complete fixed canonical functions and an exactly anchored affine family."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family as base
from p8_exceptional_vacuum import analytic
from p8_vacuum_flat_dirac_hadamard.symbols import rational

K0 = analytic.KAPPA
N = base.N
LAMBDA = analytic.VACUUM_LAMBDA_BAR
GAMMA = analytic.FIXED_GAMMA
ZETA = s.Rational(1, 10**6)


def require_scope(kappa, normalized_Y, zeta=ZETA):
    k, x, z = map(rational, (kappa, normalized_Y, zeta))
    if k < K0:
        raise ValueError("Require the anchored family kappa>=kappa0")
    if not -s.Rational(1, 4 * N) < x < s.Rational(6, 5):
        raise ValueError("Outside the fixed canonical gradient strip")
    if z != ZETA:
        raise ValueError("This family fixes zeta and the canonical vector mass")
    return k, x, z


def at_base_coordinates(v, x, ratio):
    """u=v/sqrt(ratio), X=x/ratio; v,x are FIXED canonical/base coordinates."""
    d = base.data()
    Rb = d["R"].subs({base.u: v, base.X: x}, simultaneous=True)
    Fb = d["F"].subs({base.u: v, base.X: x}, simultaneous=True)
    return 1 + (Rb - 1) / ratio, Fb / ratio


@cache
def data():
    u, X, phi, Y = s.symbols("u X Phi Y", real=True)
    K = s.Symbol("kappa", positive=True)
    f = s.Function("fixed_f")(phi, Y)
    r = s.Function("fixed_r")(phi, Y)
    pull = {phi: s.sqrt(K) * u, Y: K * X}
    RK = 1 + r.subs(pull, simultaneous=True) / K
    FK = f.subs(pull, simultaneous=True) / K
    inv = {u: phi / s.sqrt(K), X: Y / K}
    back = lambda e: e.subs(inv, simultaneous=True).doit()
    rx, rp = s.diff(r, Y), s.diff(r, phi)
    A3 = s.diff(RK, X) / X
    A4 = -A3 - s.Rational(7, 4) * s.diff(RK, X) ** 2 / RK
    A5 = s.diff(RK, X) ** 2 / (RK * X)
    canonical = {
        "F": f,
        "curvature": -(K + r) / 2,
        "A3": rx / Y,
        "A4": -rx / Y - 7 * rx**2 / (4 * (K + r)),
        "A5": rx**2 / ((K + r) * Y),
    }
    checks = {
        "full_independent_F_fixed": s.simplify(K * back(FK) - f),
        "full_nonminimal_r_fixed": s.simplify(K * (back(RK) - 1) - r),
        "physical_R_X_chain": s.simplify(back(s.diff(RK, X)) - rx),
        "physical_R_u_chain": s.simplify(back(s.diff(RK, u)) - rp / s.sqrt(K)),
        "canonical_A3_four_field_scaling": s.simplify(back(A3) / K - canonical["A3"]),
        "canonical_A4_four_field_scaling": s.simplify(back(A4) / K - canonical["A4"]),
        "canonical_A5_six_field_scaling": s.simplify(back(A5) / K**2 - canonical["A5"]),
        "full_fixed_f_K_derivative": s.diff(f, K),
        "full_fixed_r_K_derivative": s.diff(r, K),
        "full_fixed_a3_K_derivative": s.diff(rx / Y, K),
        "dependent_A4_not_falsely_fixed": s.simplify(
            s.diff(canonical["A4"], K) - 7 * rx**2 / (4 * (K + r) ** 2)
        ),
        "dependent_A5_not_falsely_fixed": s.simplify(
            s.diff(canonical["A5"], K) + rx**2 / ((K + r) ** 2 * Y)
        ),
    }
    literalR, literalF = at_base_coordinates(base.u, base.X, s.S.One)
    checks["actual_full_base_R_anchor"] = literalR - base.data()["R"]
    checks["actual_full_base_F_anchor"] = literalF - base.data()["F"]
    # The fixed canonical strip uses the base variables, NOT u,X held fixed.
    rb = s.Symbol("base_R", real=True)
    t = s.Symbol("ratio_kappa_to_kappa0", positive=True)
    convex = 1 + (rb - 1) / t
    checks["domain_R_convex_combination"] = s.expand(convex - ((1 - 1 / t) + rb / t))
    checks["domain_lower_margin"] = s.factor(
        convex
        - (s.Rational(1, 2) + s.Rational(1, 8 * N)) / t
        - (1 - 1 / t)
        - (rb - s.Rational(1, 2) - s.Rational(1, 8 * N)) / t
    )
    # The lower dictionary is defined by an analytic parameter integral:
    # q=(3X/4) R(X)^(3/4) int_0^1 sqrt(t) RX(u,tX) Ru(u,tX) R(u,tX)^(-7/4) dt.
    # Its local coefficients solve the full ODE, including the singular 1/(2X).
    a, b, c = s.symbols("a b c", real=True)
    ap, bp, cp = s.symbols("a_phi b_phi c_phi", real=True)
    Rjet = 1 + a * X**2 + b * X**3 + c * X**4
    Ujet = ap * X**2 + bp * X**3 + cp * X**4
    q4, q5, q6 = s.symbols("q4 q5 q6", real=True)
    qjet = q4 * X**4 + q5 * X**5 + q6 * X**6
    residual = (
        s.series(
            s.diff(qjet, X)
            + (1 / (2 * X) - 3 * s.diff(Rjet, X) / (4 * Rjet)) * qjet
            - 3 * s.diff(Rjet, X) * Ujet / (4 * Rjet),
            X,
            0,
            6,
        )
        .removeO()
        .expand()
    )
    sol = s.solve([residual.coeff(X, j) for j in (3, 4, 5)], (q4, q5, q6))
    for j in (3, 4, 5):
        checks["vacuum_regular_q_ODE_order_" + str(j)] = s.factor(
            residual.subs(sol).coeff(X, j)
        )
    return {
        "kappa0": K0,
        "fixed_switch_order": N,
        "fixed_lambda": LAMBDA,
        "fixed_gamma": GAMMA,
        "fixed_canonical_definitions": "f(Phi,Y)=kappa0 F_base(Phi/sqrt(kappa0),Y/kappa0); r(Phi,Y)=kappa0 [R_base(Phi/sqrt(kappa0),Y/kappa0)-1]; a3=r_Y/Y. These are the complete S6.109 functions, not their Taylor germs.",
        "exact_physical_family": {"R": RK, "F": FK},
        "exact_canonical_coefficients": canonical,
        "fixed_domain": "Every real Phi and -kappa0/4096<Y<6kappa0/5 for every kappa>=kappa0. The physical u,X strip is rescaled, not incorrectly held fixed.",
        "uniform_R_enclosure": [
            s.Rational(1, 2) + s.Rational(1, 8 * N),
            s.Rational(6, 5),
        ],
        "regular_lower_q_integral": "q(u,X)=(3X/4) R(u,X)^(3/4) integral_0^1 sqrt(t) R_X(u,tX) R_u(u,tX) R(u,tX)^(-7/4) dt; analytic at X=0 from either sign.",
        "regular_q_leading_coefficients": sol,
        "fixed_affine_construction": "At each kappa use the generic-R S6.174 regular metric chart, source-centered isotropic full-quotient mass, lower q and shift B; zeta=1e-6. Href(Phi)=H_base(Phi/sqrt(kappa0)) is fixed.",
        "anchor_boundary": "At kappa0 the entire action and maps equal S6.174/S6.176. The same bounce is NOT asserted to solve the off-base family.",
        "checks": checks,
    }


@cache
def germs():
    u, X = base.u, base.X
    d = base.data()
    # Literal derivatives of the complete rational-switch/exponential functions.
    at0 = {u: 0, X: 0}
    F = d["F"]
    targets = {
        (0, 0): 0,
        (1, 0): 0,
        (2, 0): -1,
        (0, 1): s.Rational(1, 2),
        (1, 1): 0,
        (3, 0): 0,
        (4, 0): -8 * N,
        (2, 1): 0,
        (0, 2): 2 * LAMBDA * K0,
    }
    checks = {}
    for (i, j), wanted in targets.items():
        actual = s.diff(F, u, i, X, j).subs(at0)
        checks["literal_full_F_derivative_" + str(i) + "_" + str(j)] = s.simplify(
            actual - wanted
        )
    checks["full_F_even_in_Phi"] = s.factor(F.subs(u, -u) - F)
    checks["full_R_even_in_Phi"] = s.factor(d["R"].subs(u, -u) - d["R"])
    checks["literal_full_R_YY_leading"] = s.simplify(
        s.diff(d["R"], X, 2).subs(at0) + 2 * N
    )
    checks["canonical_quartic_DHOST_fixed"] = s.cancel(-2 * N / K0 + 2 * GAMMA)
    checks["canonical_quartic_potential_fixed"] = s.cancel(-N / (3 * K0) + GAMMA / 3)
    Fxxx = s.simplify(s.diff(F, X, 3).subs(at0))
    Rxxx = s.simplify(s.diff(d["R"], X, 3).subs(at0))
    Ruuxx = s.simplify(s.diff(d["R"], u, 2, X, 2).subs(at0))
    checks["full_fixed_F_Y_cubed_nonzero_coefficient"] = Fxxx + s.Rational(655104, 25)
    checks["full_fixed_a3_Y_nonzero_coefficient"] = Rxxx - 6 * N
    checks["full_fixed_a3_Phi_squared_nonzero_coefficient"] = Ruuxx - 12 * N
    return {
        "mass_squared": s.S.One,
        "canonical_quadratic": "(Y-Phi^2)/2",
        "canonical_quartic": "lambda Y^2-gamma Phi^4/3+2gamma (L4-L3)",
        "full_flat_limit": "f(Phi,Y)+a3(Phi,Y)(L3-L4), with ALL higher independent interactions retained",
        "no_cubic_scalar_vertex": True,
        "fixed_higher_coefficients": {
            "f_Y_cubed": Fxxx / (6 * K0 * K0),
            "a3_Y": Rxxx / (2 * K0 * K0),
            "a3_Phi_squared": Ruuxx / (2 * K0 * K0),
        },
        "checks": checks,
    }
