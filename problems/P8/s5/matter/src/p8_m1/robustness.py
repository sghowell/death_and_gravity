"""Exact canonical-matter principal-cone normal form; no loop computation.

All conclusions require a two-derivative principal action in a regular
chart, K positive definite, and the frozen equal matter entries K22=G22.
No eigenvalue of a singular velocity chart is used at gamma crossing.
"""

from functools import cache

import sympy as sp

from . import nonlinear as model


def normal_form(kinetic, gradient):
    if kinetic.shape != (2, 2) or gradient.shape != (2, 2):
        raise ValueError("Expected two 2x2 matrices")
    if kinetic != kinetic.T or gradient != gradient.T:
        raise ValueError("Principal matrices must be symmetric")
    if sp.cancel(kinetic[1, 1]-gradient[1, 1]) != 0:
        raise ValueError("Equal canonical-matter principal entries are required")
    A, B, Q = kinetic[0, 0], kinetic[0, 1], kinetic[1, 1]
    C, D = gradient[0, 0], gradient[0, 1]
    k = sp.cancel(A-B**2/Q)
    if Q.is_zero or k.is_zero:
        raise ValueError("Degenerate normalization; require Q>0 and k>0")
    a, b = sp.cancel(A-C), sp.cancel(B-D)
    return {"k": k, "Q": Q, "a": a, "b": b,
            "sigma": sp.cancel((a-2*B*b/Q)/k), "r2": sp.cancel(b**2/(k*Q))}


@cache
def algebra_checks():
    A, B, C, D, Q, z = sp.symbols("A B C D Q z", real=True)
    K = sp.Matrix([[A, B], [B, Q]])
    G = sp.Matrix([[C, D], [D, Q]])
    nf = normal_form(K, G)
    k, s, r2 = nf["k"], nf["sigma"], nf["r2"]
    shear = sp.Matrix([[1, 0], [-B/Q, 1]])
    target = sp.Matrix([[k*(1-s), -nf["b"]], [-nf["b"], Q]])
    return {
        **{f"kinetic_shear_{i}{j}": sp.cancel((shear.T*K*shear-sp.diag(k, Q))[i, j])
           for i in range(2) for j in range(2)},
        **{f"gradient_shear_{i}{j}": sp.cancel((shear.T*G*shear-target)[i, j])
           for i in range(2) for j in range(2)},
        "characteristic_polynomial": sp.cancel((G-z*K).det()-k*Q*((1-z)*(1-s-z)-r2)),
        "gradient_determinant": sp.cancel(G.det()-k*Q*(1-s-r2)),
        "causal_difference_determinant": sp.cancel((K-G).det()+nf["b"]**2),
    }


@cache
def chart_checks():
    """Check both derived principal charts against crossing-safe invariants.

    Small-symbol formulas are the outputs of the pinned metric-derived
    coupled.py and gamma.py actions. The old full replay remains required.
    """
    T, theta, lam, delta, l, J, Pgamma = sp.symbols("T theta Lambda delta l J Pgamma", real=True)
    w = l*(3*delta-1)
    J0 = J+w**2/2
    E = lam-T*(1-3*delta)
    Ku = sp.Matrix([[T**2*J0/theta**2, T*w/(2*theta)], [T*w/(2*theta), sp.Rational(1, 2)]])
    Gu = sp.Matrix([[Pgamma/theta**2, -l*lam/(2*theta)], [-l*lam/(2*theta), sp.Rational(1, 2)]])
    Kg = sp.Matrix([[T**2*J0/lam**2, -T*w/(2*lam)], [-T*w/(2*lam), sp.Rational(1, 2)]])
    Gg = sp.Matrix([[Pgamma/lam**2, l/2], [l/2, sp.Rational(1, 2)]])
    nf_u, nf_g = normal_form(Ku, Gu), normal_form(Kg, Gg)
    target_r2 = l**2*E**2/(2*T**2*J)
    target_sigma = (T**2*J0-Pgamma-T*w*l*E)/(T**2*J)
    checks = {}
    for name, nf in (("unitary", nf_u), ("gamma", nf_g)):
        checks[f"{name}_r2"] = sp.cancel(nf["r2"]-target_r2)
        checks[f"{name}_sigma"] = sp.cancel(nf["sigma"]-target_sigma)
    # Independent elimination of the shift/lapse principal velocity blocks.
    vd, sd, ss, n, b, sig = sp.symbols("vd sd s n b Sigma", real=True)
    L = -3*T*vd**2+sig*n**2+6*theta*n*vd+sd**2/2+w*n*sd+b*(2*theta*n-2*T*vd-l*ss)
    sol = sp.solve([sp.diff(L, n), sp.diff(L, b)], (n, b))
    direct_K = (sp.hessian(L.subs(sol), (vd, sd))/2).subs(sig, J0-3*theta**2/T)
    checks.update({f"direct_constraint_K_{i}{j}": sp.cancel((direct_K-Ku)[i, j])
                   for i in range(2) for j in range(2)})
    f2, fx, a3, X = sp.symbols("F2 F2X A3 X", real=True)
    gt = -2*f2
    dn = X*(-4*fx-X*a3)/(2*gt)
    lambda_n = gt+4*X*fx+gt*dn
    checks["covariant_exceptional_relation"] = sp.cancel(lambda_n-gt*(1-3*dn)+2*X*(2*fx+X*a3))
    return checks


@cache
def controls():
    bounce = sp.Matrix([[6, sp.Rational(1, 20)], [sp.Rational(1, 20), sp.Rational(1, 2)]])
    mismatch = sp.diag(sp.Rational(599, 100), sp.Rational(1, 2))
    z, epsilon = sp.symbols("z epsilon", real=True)
    diag_only = bounce-sp.diag(epsilon, 0)
    outside = sp.Matrix([[sp.Rational(9, 10), sp.Rational(1, 20)],
                         [sp.Rational(1, 20), sp.Rational(9, 10)]])
    return {
        "bounce": bounce, "mismatch_gradient": mismatch,
        "mismatch_speeds": (1-1/sp.sqrt(1199), 1+1/sp.sqrt(1199)),
        "residuals": {
            "mismatch_characteristic": sp.expand((mismatch-z*bounce).det()
                                                   -(1199*z**2-2398*z+1198)/400),
            "diagonal_characteristic": sp.expand((diag_only-z*bounce).det()
                                                   -sp.Rational(1199, 400)*(1-z)*(1-200*epsilon/1199-z)),
            "outside_M1_characteristic": sp.expand((outside-z*sp.eye(2)).det()
                                                    -(sp.Rational(19, 20)-z)*(sp.Rational(17, 20)-z)),
        },
    }


@cache
def deformation():
    """Actual F deformation retaining the entire pinned background and tube.

    epsilon is fixed and finite, not a derived loop coefficient. For epsilon
    >=0 the family is principal-healthy/causal; -1<epsilon<0 is healthy but
    superluminal; epsilon=-1 degenerates at the bounce. The action's Ia
    completion is unmodified because only F changes.
    """
    epsilon, z = sp.symbols("epsilon z", real=True)
    f, X, N, u, d = model.functions(), model.X, model.N, model.u, model.d
    delta_F = epsilon*f["J"]*(X-1)**2/(4*d**2)
    delta_S = sp.factor(sp.diff(N*delta_F.subs(X, N**-2), N, 2).subs(N, 1)/2)
    checks = {
        "background_F": delta_F.subs(X, 1),
        "background_FX": sp.diff(delta_F, X).subs(X, 1),
        "background_Fphi": sp.diff(delta_F, u).subs(X, 1),
        "lapse_sigma_shift": sp.cancel(delta_S-epsilon*f["J"]/d**2),
        "canonical_tail_J": sp.limit(d*f["J"], u, sp.oo)-4,
    }
    for label, denom in (("unitary", f["Theta"]), ("gamma", f["Lambda"])):
        base_K = f["K"] if label == "unitary" else f["K_gamma"]
        base_G = f["G"] if label == "unitary" else f["G_gamma"]
        new_K = base_K+sp.diag(delta_S/denom**2, 0)
        detK = sp.factor(base_K.det())
        # Use the rank-one determinant identity, avoiding expansion of the
        # already-pinned high-degree rational witness.
        coefficient = sp.cancel(delta_S*base_K[1, 1]/denom**2/detK)
        checks[f"{label}_rank_one_ratio"] = sp.cancel(coefficient-epsilon/d**2)
        checks.update({f"{label}_base_luminal_{i}{j}": sp.cancel((base_K-base_G)[i, j])
                       for i in range(2) for j in range(2)})
        # det(G-z(K+deltaK))/detK = (1-z)^2-z(1-z)*coefficient.
        checks[f"{label}_speed_factorization"] = sp.expand((1-z)**2-z*(1-z)*epsilon/d**2
                                                           -(1-z)*(1-(1+epsilon/d**2)*z))
        checks[f"{label}_new_determinant"] = sp.cancel(new_K.det()/detK-1-epsilon/d**2)
    # With p=4 the old normalized clock field has X=d*Z/(2p)=d*Z/8.
    # Coefficients of 1,Z,Z^2 are the additional canonical tail terms.
    canonical_Z = sp.Symbol("canonical_Z", real=True)
    tail_coefficients = (epsilon*f["J"]/(4*d**2), -epsilon*f["J"]/(16*d), epsilon*f["J"]/256)
    checks["canonical_tail_polynomial"] = sp.cancel(delta_F.subs(X, d*canonical_Z/8)
                                                   -sum(c*canonical_Z**k for k, c in enumerate(tail_coefficients)))
    for k, coeff in enumerate(tail_coefficients):
        checks[f"canonical_tail_limit_{k}"] = sp.limit(coeff, u, sp.oo)
    return {"epsilon": epsilon, "delta_F": delta_F, "delta_S": delta_S,
            "tail_coefficients": tail_coefficients, "residuals": checks}


@cache
def deformation_tail_bounds():
    """Quantitative old-contract remainder and first-two derivative bounds.

    P=d*J_star, Z=1-x^2=1/d, and d/dpsi=Z*d/dx/sqrt(8).
    First derivatives are squared to keep every reported majorant rational.
    The three added tail coefficients have no remaining X dependence.
    """
    from .series import compact_even, polynomial_bound, x

    P = compact_even(model.d*model.functions()["J"])
    Z = 1-x**2
    epsilon = deformation()["epsilon"]
    rows, checks = [], {}
    for k, (coefficient, n, actual) in enumerate(zip(
            (sp.Rational(1, 4), -sp.Rational(1, 16), sp.Rational(1, 256)), (3, 2, 1),
            deformation()["tail_coefficients"], strict=True)):
        P1 = sp.expand(Z*sp.diff(P, x)-2*n*x*P)
        P2 = sp.expand(Z*sp.diff(P1, x)-2*n*x*P1)
        compact = coefficient*epsilon*Z**n*P
        checks[f"tail_coefficient_{k}"] = sp.cancel(compact.subs(x, model.u/sp.sqrt(model.d))-actual)
        checks[f"tail_Dpsi_{k}"] = sp.expand(Z*sp.diff(Z**n*P, x)-Z**n*P1)
        checks[f"tail_Dpsi2_{k}"] = sp.expand(Z*sp.diff(Z**n*P1, x)-Z**n*P2)
        checks[f"tail_X_first_{k}"] = sp.diff(actual, model.X)
        checks[f"tail_X_second_{k}"] = sp.diff(actual, model.X, 2)
        checks[f"tail_X_mixed_{k}"] = sp.diff(actual, model.X, model.u)
        rows.append({"canonical_polynomial_degree": k, "decay_power_n": n,
                     "compact_P": str(P), "compact_first_factor": str(P1), "compact_second_factor": str(P2),
                     "coefficient_bound_per_abs_epsilon": str(8*abs(coefficient)),
                     "first_derivative_squared_bound_per_epsilon_squared": str(coefficient**2*polynomial_bound(P1)**2/8),
                     "first_derivative_squared_decay_power": 2*n,
                     "second_derivative_bound_per_abs_epsilon": str(abs(coefficient)*polynomial_bound(P2)/8),
                     "X_and_mixed_derivatives": "0"})
    return {"method": "P=dJ<8 for values; exact compact polynomial l1 for first squared and second canonical derivatives",
            "domain": "all real phi, fixed finite epsilon; divide each bound by d to its recorded decay power",
            "rows": rows, "residuals": checks,
            "new_canonical_kinetic_ge_one_quarter_sufficient_d_for_epsilon_ge_zero": "10481/1680+2*epsilon"}


@cache
def audit_checks():
    f = model.functions()
    X = model.X
    return {**algebra_checks(), **chart_checks(), **controls()["residuals"], **deformation()["residuals"],
            **deformation_tail_bounds()["residuals"],
            "pinned_exceptional_tube": sp.cancel(f["A3"]+2*sp.diff(f["F2"], X)/X),
            "pinned_bounce_J": f["J"].subs(model.u, 0)-sp.Rational(1199, 800),
            **{f"pinned_bounce_{i}{j}": (f["K_gamma"].subs(model.u, 0)-controls()["bounce"])[i, j]
               for i in range(2) for j in range(2)}}
