"""Two unrestricted distortion traces and their rolling quadratic Schur action.

Only frozen S6.37/S6.38 scientific inputs are imported.  The zero-Schur
identity is a QUADRATIC rolling identity, not a nonlinear inverse on a tube.
Source signature (-+++); the Schur Lorentz matrix below is its negative.
"""

from functools import cache

import sympy as sp
from p8_affine import connection as old
from p8_affine_kinetic import vector as previous

A, B = sp.symbols("A B", real=True)
P, X, H, HP, Q = old.P, previous.X, previous.H, previous.HP, previous.Q
ETA = sp.ImmutableMatrix(sp.diag(1, -1, -1, -1))


def _clean(value):
    value = sp.factor(value)
    return sp.S.Zero if value.is_zero is True else value


def _matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(_clean))


def _null_remainder(value):
    return _clean(sp.rem(sp.Poly(sp.expand(value), A),
                         sp.Poly(A**2-14*A*B+B**2, A)).as_expr())


@cache
def traces():
    """Construct V and U on all 64 components; impose no torsion symmetry."""
    v, u = sp.zeros(4, 64), sp.zeros(4, 64)
    for mu in range(4):
        for a in range(4):
            v[mu, old.index(a, mu, a)] += 1
            v[mu, old.index(a, a, mu)] -= sp.Rational(1, 4)
            u[mu, old.index(mu, a, a)] += old.SIGN[mu]*old.SIGN[a]
            u[mu, old.index(a, a, mu)] -= sp.Rational(1, 4)
    full = v.col_join(u)
    quotient = old.quotient()
    return {"V": sp.ImmutableMatrix(v), "U": sp.ImmutableMatrix(u),
            "full": sp.ImmutableMatrix(full),
            "quotient": _matrix(full*quotient["embedding"]),
            "projective_residual": full*old.quadratic()["gauge"],
            "source_metric": -ETA, "Schur_Lorentz_matrix": ETA}


@cache
def schur():
    """Solve all quotient blocks for both traces before restricting p=1/2."""
    quotient = old.quotient()
    n = traces()["quotient"]
    inverse_n = sp.zeros(60, 8)
    for indices, block in zip(quotient["blocks"], quotient["block_matrices"], strict=True):
        solution = block.inv(method="DM")*n.extract(range(8), indices).T
        for i, component in enumerate(indices):
            for j in range(8):
                inverse_n[component, j] = _clean(solution[i, j])
    d = _matrix(n*inverse_n)
    clock = {P: sp.Rational(1, 2)}
    expected = sp.Rational(3, 8)*sp.kronecker_product(sp.Matrix([[1, -7], [-7, 1]]), ETA)
    return {"M": quotient["hessian"], "M_clock": _matrix(quotient["hessian"].subs(clock)),
            "N": n, "inverse_N_transpose": sp.ImmutableMatrix(inverse_n),
            "inverse_N_transpose_clock": _matrix(inverse_n.subs(clock)),
            "D": d, "D_clock": _matrix(d.subs(clock)),
            "D_clock_expected": sp.ImmutableMatrix(expected),
            "embedding": quotient["embedding"],
            "stationary": _matrix(traces()["full"]*old.eliminate()["solution"])}


@cache
def covariant():
    """Derive all Ustar coefficients from the full 64-source stationary field.

    Ustar=a*v+b*H.v+c*v*Box(phi)+d*v*(v.H.v).  The parent affine
    coefficient c_X is distinct from the vector coefficient named c.
    """
    actual = schur()["stationary"][4:8, :]
    poly = sp.Poly(actual[0], *old.H_SYMBOLS)
    a = _clean(poly.coeff_monomial(1)/old.S)
    b = _clean(-actual[1]/(old.S*old.H[0, 1]))
    c = _clean(poly.coeff_monomial(old.H[1, 1])/old.S)
    d = _clean((poly.coeff_monomial(old.H[0, 0])/old.S+b+c)/old.S**2)
    gradient = sp.Matrix([old.S, 0, 0, 0])
    rebuilt = (a*gradient-b*old.S*old.H[:, 0]
               +c*gradient*old.literal()["trace_H"]+d*gradient*old.literal()["vHv"])
    generic = {name: _clean(value.subs(old.S**2, -X))
               for name, value in zip(("a", "b", "c", "d"), (a, b, c, d), strict=True)}
    fphi = HP*(1+X)/(2*H**2)
    px = -1/(8*H*P)
    cx = 2*(1-4*P)*px/X-2*(P-2*P**2)/X**2
    substitutions = {old.PX: px, old.PP: fphi/(4*P), old.CX: cx,
                     old.F3: -(Q+fphi)/(4*P*X)}
    coefficients = {name: _clean(value.subs(substitutions)) for name, value in generic.items()}
    expected = {
        "a": 3*(P-1)*Q/(4*P)+3*(8*P**2-8*P-1)*fphi/(16*P**2),
        "b": (2*P-1)/X+(3*P+1)/(4*H*P**2),
        "c": (2*P-1)**2/(2*X),
        "d": -(4*P**2-1)/(2*X**2)-(12*P**2+1)/(16*H*P**2*X),
    }
    clock = {P: sp.Rational(1, 2), X: -1, Q: 0}
    # The frozen lower ODE gives q_X=0 on its q(u,-1)=0 boundary.
    derivatives = {name: _clean((sp.diff(value, X)+px*sp.diff(value, P)).subs(clock))
                   for name, value in coefficients.items()}
    return {"generic": generic, "coefficients": coefficients, "expected": expected,
            "reconstruction_residual": _matrix(actual-rebuilt),
            "coefficient_residuals": {name: _clean(coefficients[name]-expected[name]) for name in expected},
            "clock_values": {name: _clean(value.subs(clock)) for name, value in coefficients.items()},
            "clock_X_derivatives": derivatives}


@cache
def rolling():
    """Use the full frozen ADM Hessian/inverse, including shift and lapse jets."""
    adm = previous.adm_linearization()
    symbols = adm["symbols"]
    eps, n, hubble = symbols["eps"], symbols["n"], symbols["Hubble"]
    data = covariant()
    coefficients = {name: data["clock_values"][name]+2*eps*n*data["clock_X_derivatives"][name]
                    for name in ("a", "b", "c", "d")}
    gradient = sp.Matrix([1, 0, 0, 0])
    full = (coefficients["a"]*gradient+coefficients["b"]*adm["Hv"]
            +coefficients["c"]*gradient*adm["Box"]+coefficients["d"]*gradient*adm["vHv"])
    first = _matrix(full.applyfunc(lambda value: sp.expand(value).coeff(eps, 1)))
    first_cd = _matrix(first.subs(HP, 3*hubble*H/2))
    expected = sp.Matrix([3*symbols["n_dot"]/(2*H)-27*hubble*n/(8*H),
                          *(5*entry/(2*H) for entry in symbols["n_gradient"])])
    combined = _matrix(A*adm["first_variation_CD"]+B*first_cd)
    alpha = 3*(A-B)/(2*H)
    lapse = -3*(A+9*B)*hubble/(8*H)
    alpha_prime = -9*(A-B)*hubble/(4*H)
    d = _clean(-alpha_prime-lapse)
    e = (A-B)/H
    curl_value = _clean(15*(A-B)*hubble/(4*H)-lapse)
    return {"symbols": symbols, "U_first_variation_before_clock": first,
            "U_first_variation": first_cd, "U_expected": sp.ImmutableMatrix(expected),
            "T_first_variation": combined, "alpha": alpha, "alpha_prime": alpha_prime,
            "T_lapse_value": lapse, "d": d, "e": e,
            "curl_lapse_dot": -e, "curl_lapse_value": curl_value,
            "background_source": previous.background()["source"],
            "background_stationary": previous.background()["stationary"]}


@cache
def rank_one():
    """No inverse of gamma is taken: this API includes the null-Schur cases."""
    data = schur()
    combine = A*sp.eye(4).row_join(sp.zeros(4))+B*sp.zeros(4).row_join(sp.eye(4))
    n = _matrix(combine*data["N"])
    w = _matrix(data["inverse_N_transpose_clock"]*combine.T)
    full = _matrix(combine*traces()["full"])
    gamma = 3*(A**2-14*A*B+B**2)/8
    kept = old.quotient()["kept"]
    v_columns = [kept.index(old.index((mu+1) % 4, mu, (mu+1) % 4)) for mu in range(4)]
    u_columns = [kept.index(old.index(mu, (mu+1) % 4, (mu+1) % 4)) for mu in range(4)]
    u_signs = sp.diag(*(old.SIGN[mu]*old.SIGN[(mu+1) % 4] for mu in range(4)))
    return {"N": n, "W": w, "full": full, "gamma": gamma,
            "D": _matrix(n*w), "M": data["M_clock"],
            "V_minor": n.extract(range(4), v_columns),
            "U_minor": n.extract(range(4), u_columns), "U_minor_signs": u_signs,
            "null_roots_A_over_B": (7-4*sp.sqrt(3), 7+4*sp.sqrt(3)),
            "full_vector_rank_for_nonzero_AB": 4}


@cache
def action_identity():
    """Multiplier reduction retains ALL 64 connection Euler equations.

    Write y=k-kstar and add ell.(T-Tstar-Ny).  Eliminating y gives
    S_K[T]+ell.(T-Tstar)-ell.D.ell/2, valid including gamma=0.
    EL(S_K) is the volume-normalized functional derivative; no
    frozen-frequency approximation of the Maxwell operator is used.
    """
    data = rank_one()
    ell = sp.Matrix(sp.symbols("ell_0:4", real=True))
    delta = sp.Matrix(sp.symbols("Delta_0:4", real=True))
    kinetic_euler = sp.Matrix(sp.symbols("EulerK_0:4", real=True))
    y = data["W"]*ell
    algebraic = _clean((y.T*data["M"]*y)[0]/2+(ell.T*(delta-data["N"]*y))[0])
    reduced = _clean((ell.T*delta)[0]-(ell.T*data["D"]*ell)[0]/2)
    full_m = old.quadratic()["hessian"].subs(P, sp.Rational(1, 2))
    e = schur()["embedding"]
    return {"multiplier": ell, "Delta": delta, "kinetic_Euler": kinetic_euler,
            "algebraic_substitution": algebraic, "reduced_multiplier_action": reduced,
            "substitution_residual": _clean(algebraic-reduced),
            "quotient_Euler_residual": _matrix(data["M"]*y-data["N"].T*ell),
            "full_Euler_residual": _matrix(full_m*e*y-data["full"].T*ell),
            "multiplier_Euler_residual": _matrix(sp.Matrix([sp.diff(reduced, component) for component in ell])
                                                   -delta+data["D"]*ell),
            "stationary_y": _matrix(-data["W"]*kinetic_euler),
            "stationary_full_Euler_residual": _matrix(-full_m*e*data["W"]*kinetic_euler
                                                       +data["full"].T*kinetic_euler)}


@cache
def null_schur():
    """Quadratic rolling null identity; not a nonlinear fixed-T inverse."""
    data, action = rank_one(), action_identity()
    remainder = lambda matrix: sp.ImmutableMatrix(matrix.applyfunc(_null_remainder))
    energy = (action["stationary_y"].T*data["M"]*action["stationary_y"])[0]/2
    roll = rolling()
    n, ndot = roll["symbols"]["n"], roll["symbols"]["n_dot"]
    q = sp.Symbol("q_spatial_normalized", positive=True)
    kinetic = sp.Symbol("lambda_normalized", real=True, nonzero=True)
    scalar_density = kinetic*q*(roll["curl_lapse_dot"]*ndot+roll["curl_lapse_value"]*n)**2/2
    return {"relation": A**2-14*A*B+B**2,
            "N_W_remainder": remainder(data["N"]*data["W"]),
            "null_Gram_remainder": remainder(data["W"].T*data["M"]*data["W"]),
            "T_minus_Tstar_remainder": remainder(data["N"]*action["stationary_y"]),
            "stationary_algebraic_energy_remainder": _null_remainder(energy),
            "multiplier_action_remainder": _null_remainder(action["reduced_multiplier_action"]
                                                            -(action["multiplier"].T*action["Delta"])[0]),
            "kernel_dimension": 56, "kernel_radical_dimension": 4,
            "remaining_nondegenerate_kernel_quotient_dimension": 52,
            "q_spatial": q, "kinetic_coefficient": kinetic,
            "Maxwell_scalar_density_divided_by_a_cubed": scalar_density,
            "lapse_dot_squared_coefficient": _clean(sp.diff(scalar_density, ndot, 2)/2),
            "nonlinear_open_tube_inverse_claim": False}


def _exact(value, name):
    if isinstance(value, (bool, float, str)):
        raise TypeError(f"{name} must be an exact real constant")
    value = sp.sympify(value)
    if not isinstance(value, sp.Expr) or value.has(sp.Float) or value.free_symbols:
        raise TypeError(f"{name} must be an exact real constant")
    if value.is_real is not True or value.is_finite is not True:
        raise ValueError(f"{name} must be finite and real")
    return value


def require_parameters(a, b, *, kinetic_coefficient=1, mass_squared=1, time_scale=1):
    """All nonzero real curl coefficients; no gamma=0 or negative sign dropped."""
    a, b, zeta, mass, tau = (_exact(value, name) for value, name in
                            ((a, "A"), (b, "B"), (kinetic_coefficient, "kinetic_coefficient"),
                             (mass_squared, "mass_squared"), (time_scale, "time_scale")))
    if sp.simplify(a**2+b**2).is_positive is not True:
        raise ValueError("Require (A,B) != (0,0)")
    if zeta.is_zero is not False:
        raise ValueError("The kinetic coefficient must be nonzero")
    if mass.is_positive is not True or tau.is_positive is not True:
        raise ValueError("Physical M² and tau must be positive")
    return {"A": a, "B": b, "zeta_physical": zeta, "M_squared": mass, "tau": tau,
            "gamma": sp.simplify(3*(a**2-14*a*b+b**2)/8)}


def units(a, b, *, kinetic_coefficient=1, mass_squared=1, time_scale=1):
    data = require_parameters(a, b, kinetic_coefficient=kinetic_coefficient,
                              mass_squared=mass_squared, time_scale=time_scale)
    return data | {"lambda_normalized": data["zeta_physical"]/(data["M_squared"]*data["tau"]**2),
                   "D_physical_factor": 1/data["M_squared"], "trace_normalized_factor": data["tau"]}


@cache
def checks():
    s, r, c, roll, action, null = schur(), rank_one(), covariant(), rolling(), action_identity(), null_schur()
    hcal = roll["symbols"]["Hubble"]
    result = {
        "projective_invariance": traces()["projective_residual"],
        "frozen_V_bridge": traces()["V"]-previous.vector_map()["full"],
        "all_block_inverse": _matrix(s["M"]*s["inverse_N_transpose"]-s["N"].T),
        "on_clock_eight_Schur": _matrix(s["D_clock"]-s["D_clock_expected"]),
        "U_covariant_reconstruction": c["reconstruction_residual"],
        "U_rolling_first_variation": _matrix(roll["U_first_variation"]-roll["U_expected"]),
        "background_all_sources": roll["background_source"],
        "background_all_connection": roll["background_stationary"],
        "combination_Schur": _matrix(r["D"]-r["gamma"]*ETA),
        "rank_A_minor": _matrix(r["V_minor"]-A*sp.eye(4)),
        "rank_B_minor": _matrix(r["U_minor"]-B*r["U_minor_signs"]),
        "rolling_d": _clean(roll["d"]-3*(7*A+3*B)*hcal/(8*H)),
        "rolling_e": _clean(roll["e"]-(A-B)/H),
        "rolling_curl_value": _clean(roll["curl_lapse_value"]-3*(11*A-B)*hcal/(8*H)),
        "null_Maxwell_dot_coefficient": _clean(null["lapse_dot_squared_coefficient"]
                                                -null["kinetic_coefficient"]*null["q_spatial"]*(A-B)**2/(2*H**2)),
    }
    result.update({f"U_coefficient_{name}": value for name, value in c["coefficient_residuals"].items()})
    result.update({f"action_{name}": value for name, value in action.items() if name.endswith("residual")})
    result.update({f"null_{name}": value for name, value in null.items() if name.endswith("remainder")})
    return result


@cache
def calibration():
    return {"source_signature": "-+++", "full_connection_dimension": 64,
            "projective_dimension": 4, "quotient_dimension": 60,
            "two_trace_rank": 8, "nonzero_combination_rank": 4,
            "D_clock": schur()["D_clock"], "gamma": rank_one()["gamma"],
            "d": rolling()["d"], "e": rolling()["e"],
            "zero_Schur_is_quadratic_rolling_only": True,
            "nonlinear_null_inverse_or_full_UV_claim": False}
