"""Selective quotient-vector elimination and actual CD rolling linearization.

Only frozen S6.37 source modules are imported.  A Maxwell-like kinetic term
for the displayed quotient vector is a NEW action, not a completion claim.
All formulas initially use M=tau=1 and the source signature (-+++).
"""

from functools import cache

import sympy as sp
from p8_affine import connection as old
from p8_affine import dictionary, lower

P = old.P
X = sp.Symbol("x", negative=True)
H = sp.Symbol("h", positive=True)
HP = sp.Symbol("h_phi", real=True)
Q = sp.Symbol("q", real=True)


def _clean(value):
    answer = sp.factor(value)
    return sp.S.Zero if answer.is_zero is True else answer


def _matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(_clean))


@cache
def vector_map():
    """V_mu=kappa^a_{mu a}-kappa^a_{a mu}/4, without a torsion restriction."""
    tensor = sp.zeros(4, 64)
    for mu in range(4):
        for a in range(4):
            tensor[mu, old.index(a, mu, a)] += 1
            tensor[mu, old.index(a, a, mu)] -= sp.Rational(1, 4)
    quotient = old.quotient()
    return {"full": sp.ImmutableMatrix(tensor),
            "quotient": _matrix(tensor*quotient["embedding"]),
            "projective_residual": tensor*old.quadratic()["gauge"]}


@cache
def schur():
    """Derive the four-vector Schur matrix using every frozen quotient block.

    The complement has dimension 56.  No 56-component ansatz or discarded
    Euler equation is used.  D's sign is NOT a full coupled spectrum test.
    """
    previous = old.quotient()
    n = vector_map()["quotient"]
    inverse_n = sp.zeros(60, 4)
    for indices, block in zip(previous["blocks"], previous["block_matrices"], strict=True):
        answer = block.inv(method="DM")*n.extract(range(4), indices).T
        for local, component in enumerate(indices):
            for j in range(4):
                inverse_n[component, j] = _clean(answer[local, j])
    d = _matrix(n*inverse_n)
    dinv = _matrix(d.inv())
    constrained_lift = _matrix(inverse_n*dinv)
    projector = _matrix(sp.eye(60)-constrained_lift*n)
    target = _matrix(vector_map()["full"]*old.eliminate()["solution"])
    return {"M": previous["hessian"], "N": n, "D": d, "D_inverse": dinv,
            "inverse_N_transpose": sp.ImmutableMatrix(inverse_n),
            "lift": constrained_lift, "complement_projector": projector,
            "Vstar": target, "full_stationary": old.eliminate()["solution"],
            "quotient_stationary": old.eliminate()["quotient_solution"],
            "complement_dimension": 56}


@cache
def rest_coefficients():
    """Read Vstar=a*v+b*H.v+c*v*Box(phi)+d*v*(v.H.v) from all 64 sources."""
    result = schur()["Vstar"]
    zero = sp.Poly(result[0], *old.H_SYMBOLS)
    a = _clean(zero.coeff_monomial(1)/old.S)
    b = _clean(-result[1]/(old.S*old.H[0, 1]))
    c = _clean(zero.coeff_monomial(old.H[1, 1])/old.S)
    d = _clean((zero.coeff_monomial(old.H[0, 0])/old.S+b+c)/old.S**2)
    v = sp.Matrix([old.S, 0, 0, 0])
    hv = -old.S*old.H[:, 0]
    reconstructed = (a*v+b*hv+c*v*old.literal()["trace_H"]
                     +d*v*old.literal()["vHv"])
    return {"a": _clean(a.subs(old.S**2, -X)),
            "b": _clean(b.subs(old.S**2, -X)),
            "c": _clean(c.subs(old.S**2, -X)),
            "d": _clean(d.subs(old.S**2, -X)),
            "reconstruction_residual": _matrix(result-reconstructed)}


@cache
def cd_coefficients():
    """Actual coefficient jets of the CD lift, with its smooth q kept explicit.

    p²=(h-1-x)/(4h), p_x=-1/(8hp), f_phi=h_phi(1+x)/(2h²).
    q is the unique frozen lower.ode() coefficient with q(u,-1)=0.
    """
    fphi = HP*(1+X)/(2*H**2)
    px = -1/(8*H*P)
    cx = 2*(1-4*P)*px/X-2*(P-2*P**2)/X**2
    substitutions = {old.PX: px, old.PP: fphi/(4*P), old.CX: cx,
                     old.F3: -(Q+fphi)/(4*P*X)}
    actual = {name: _clean(rest_coefficients()[name].subs(substitutions))
              for name in ("a", "b", "c", "d")}
    expected = {
        "a": 3*(P+1)*Q/(4*P)+3*(8*P**2+8*P-1)*fphi/(16*P**2),
        "b": -(2*P-1)/X-(3*P+1)/(4*H*P**2),
        "c": (2*P-1)*(2*P+3)/(2*X),
        "d": -(4*P**2-1)/(2*X**2)-(12*P**2-7)/(16*H*P**2*X),
    }
    clock = {P: sp.Rational(1, 2), X: -1, Q: 0}
    qx = sp.Symbol("q_x", real=True)
    derivatives = {name: _clean((sp.diff(value, X)+px*sp.diff(value, P)
                                +qx*sp.diff(value, Q)).subs(clock).subs(qx, 0))
                   for name, value in actual.items()}
    return {"coefficients": actual, "expected": expected, "f_phi": fphi,
            "clock_values": {name: _clean(value.subs(clock))
                             for name, value in actual.items()},
            "clock_X_derivatives": derivatives,
            "residuals": {name: _clean(actual[name]-expected[name]) for name in actual}}


def _exact_constant(value, name):
    if isinstance(value, (bool, float, str)):
        raise TypeError(f"{name} requires an exact scalar, not a boolean, float, or string")
    result = sp.sympify(value)
    if not isinstance(result, sp.Expr):
        raise TypeError(f"{name} requires an exact scalar")
    if result.has(sp.Float) or result.free_symbols:
        raise ValueError(f"{name} must be an exact constant")
    if result.is_real is not True or result.is_finite is not True:
        raise ValueError(f"{name} must be real and finite")
    return result


def require_domain(p, *, mass_squared=1, time_scale=1):
    """Closed CD tube and positive physical units, without float sampling."""
    p, mass_squared, time_scale = (_exact_constant(value, name) for value, name in
                                  ((p, "p"), (mass_squared, "mass_squared"),
                                   (time_scale, "time_scale")))
    if p.is_positive is not True:
        raise ValueError("The positive p branch is required")
    if ((p**2-sp.Rational(9, 40)).is_nonnegative is not True
            or (sp.Rational(11, 40)-p**2).is_nonnegative is not True):
        raise ValueError("Require 9/40 <= p² <= 11/40")
    if mass_squared.is_positive is not True or time_scale.is_positive is not True:
        raise ValueError("Mass squared and time scale must be positive")
    return p, mass_squared, time_scale


def units(*, mass_squared=1, time_scale=1, kinetic_coefficient=1):
    """Restore units for the named positive curl coefficient, not a gap proof.

    V_norm=tau*V_phys; after extracting M²*tau² from the action,
    lambda_norm=zeta_phys/(M²*tau²), even though zeta_phys is dimensionless.
    """
    mass_squared, time_scale, kinetic_coefficient = (
        _exact_constant(value, name) for value, name in
        ((mass_squared, "mass_squared"), (time_scale, "time_scale"),
         (kinetic_coefficient, "kinetic_coefficient")))
    if any(value.is_positive is not True for value in
           (mass_squared, time_scale, kinetic_coefficient)):
        raise ValueError("All three physical normalization parameters must be positive")
    return {"M_squared": mass_squared, "tau": time_scale,
            "zeta_physical": kinetic_coefficient,
            "lambda_normalized": kinetic_coefficient/(mass_squared*time_scale**2),
            "D_physical_factor": 1/mass_squared, "V_normalized_factor": time_scale,
            "isolated_center_mass_squared_physical": 8*mass_squared/(3*kinetic_coefficient),
            "isolated_center_mass_squared_normalized": 8*mass_squared*time_scale**2/(3*kinetic_coefficient)}


@cache
def domain_bounds():
    """Continuous rational bounds; the narrower lower p bracket is essential."""
    lo, hi = sp.Rational(47, 100), sp.Rational(11, 20)
    polynomial = P**3+2*P**2+P-1
    fmin = polynomial.subs(P, lo)
    d0_lower = 3*fmin/(2*hi)
    minus_di_lower = (5-4*hi)*(4*lo-1)/(32*hi**2)
    return {"p_lower": lo, "p_upper": hi,
            "lower_p_squared_margin": sp.Rational(9, 40)-lo**2,
            "upper_p_squared_margin": hi**2-sp.Rational(11, 40),
            "polynomial": polynomial, "polynomial_lower": fmin,
            "polynomial_derivative": sp.diff(polynomial, P),
            "D00_lower": d0_lower, "minus_Dii_lower": minus_di_lower,
            "D_inverse_norm_upper": max(1/d0_lower, 1/minus_di_lower),
            "D_inverse_norm_margin_below_24": 24-max(1/d0_lower, 1/minus_di_lower)}


@cache
def background():
    """Every affine source vanishes on the actual constant-X rolling trajectory."""
    values = dictionary.lift()
    clock = {dictionary.x: -1}
    substitutions = {old.S: 1, old.P: sp.Rational(1, 2), old.F3: 0,
                     old.PP: values["p_phi"].subs(clock),
                     old.PX: values["px"].subs(clock),
                     old.CP: sp.diff(values["c"], dictionary.u).subs(clock),
                     old.CX: values["cx"].subs(clock)}
    substitutions.update({old.H[0, i]: 0 for i in range(4)})
    raw_source = _matrix(old.quadratic()["source"].subs(substitutions))
    stationary = _matrix(old.eliminate()["solution"].subs(substitutions))
    u = dictionary.u
    scale = (1+u**2)**2
    return {"u": u, "a": scale, "Hubble": 4*u/(1+u**2),
            "h": values["h"], "source": raw_source, "stationary": stationary,
            "q_clock": sp.S.Zero,
            "qx_clock": _clean(lower.ode()["forcing"].subs(clock)),
            "p_clock": _clean(values["p"].subs(clock)),
            "c_clock": _clean(values["c"].subs(clock)),
            "quartic_clock": _clean(values["quartic"].subs(clock)),
            "p_phi_clock": _clean(values["p_phi"].subs(clock)),
            "f_phi_clock": _clean(values["f_phi"].subs(clock))}


@cache
def adm_linearization():
    """Differentiate the full ADM metric, inverse, and clock Hessian to first order.

    N=1+eps*n, N_i=eps*shift_i, gamma_ij=a²(1+2eps*zeta)delta_ij.
    Scalar shift is a specialization; cancellation holds for arbitrary shift_i.
    """
    eps = sp.Symbol("eps", real=True)
    a = sp.Symbol("a", positive=True)
    hubble, n, ndot, zeta, zetadot = sp.symbols("Hubble n n_dot zeta zeta_dot", real=True)
    ni = sp.symbols("n_1:4", real=True)
    shift = sp.symbols("shift_1:4", real=True)
    shift_dot = sp.symbols("shift_dot_1:4", real=True)
    zeta_i = sp.symbols("zeta_1:4", real=True)
    shift_gradient = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"shift_{i+1}_{j+1}", real=True))

    def truncate(value):
        expanded = sp.expand(value)
        return expanded.coeff(eps, 0)+eps*expanded.coeff(eps, 1)

    metric = sp.diag(-1-2*eps*n, *(a**2*(1+2*eps*zeta) for _ in range(3)))
    inverse = sp.diag(-1+2*eps*n, *((1-2*eps*zeta)/a**2 for _ in range(3)))
    for i in range(3):
        metric[0, i+1] = metric[i+1, 0] = eps*shift[i]
        inverse[0, i+1] = inverse[i+1, 0] = eps*shift[i]/a**2
    derivatives = [sp.zeros(4) for _ in range(4)]
    derivatives[0][0, 0] = -2*eps*ndot
    for i in range(3):
        derivatives[0][0, i+1] = derivatives[0][i+1, 0] = eps*shift_dot[i]
        derivatives[i+1][0, 0] = -2*eps*ni[i]
        for j in range(3):
            derivatives[0][i+1, j+1] = (2*a**2*(hubble+eps*(2*hubble*zeta+zetadot))
                                       *int(i == j))
            derivatives[j+1][0, i+1] = derivatives[j+1][i+1, 0] = eps*shift_gradient[i, j]
            for k in range(3):
                derivatives[k+1][i+1, j+1] = 2*eps*a**2*zeta_i[k]*int(i == j)
    hessian = sp.Matrix(4, 4, lambda mu, nu: truncate(-sum(
        inverse[0, lam]*(derivatives[mu][lam, nu]+derivatives[nu][lam, mu]
                         -derivatives[lam][mu, nu])/2 for lam in range(4))))
    v = sp.Matrix([1, 0, 0, 0])
    raised = inverse*v
    hv = (hessian*raised).applyfunc(truncate)
    vhv = truncate((raised.T*hessian*raised)[0])
    box = truncate(sum(inverse[i, j]*hessian[i, j] for i in range(4) for j in range(4)))
    jets = cd_coefficients()
    coefficients = {name: jets["clock_values"][name]
                    +2*eps*n*jets["clock_X_derivatives"][name] for name in ("a", "b", "c", "d")}
    vector = (coefficients["a"]*v+coefficients["b"]*hv
              +coefficients["c"]*v*box+coefficients["d"]*v*vhv).applyfunc(truncate)
    first = _matrix(vector.diff(eps))
    first_cd = _matrix(first.subs(HP, 3*hubble*H/2))
    expected = sp.Matrix([-3*ndot/(2*H)-3*hubble*n/(8*H),
                          *(-5*value/(2*H) for value in ni)])
    return {"symbols": {"eps": eps, "a": a, "Hubble": hubble, "n": n,
                         "n_dot": ndot, "zeta": zeta, "zeta_dot": zetadot,
                         "n_gradient": ni, "shift": shift, "shift_dot": shift_dot,
                         "shift_gradient": shift_gradient, "zeta_gradient": zeta_i},
            "metric": sp.ImmutableMatrix(metric), "inverse": sp.ImmutableMatrix(inverse),
            "metric_inverse_residual": _matrix((metric*inverse-sp.eye(4)).applyfunc(truncate)),
            "Hessian": sp.ImmutableMatrix(hessian), "Hv": sp.ImmutableMatrix(hv),
            "vHv": vhv, "Box": box, "Vstar": sp.ImmutableMatrix(vector),
            "first_variation": first, "first_variation_CD": first_cd,
            "expected_first_variation": sp.ImmutableMatrix(expected),
            "first_variation_residual": _matrix(first_cd-expected)}


@cache
def rolling():
    """The exact-gradient linear change W=V+d[3n/(2h)] keeps F unchanged."""
    previous = background()
    u, h, hubble = previous["u"], previous["h"], previous["Hubble"]
    amplitude = 3/(2*h)
    lapse_term = -3*hubble/(8*h)
    mass_time_shift = _clean(-sp.diff(amplitude, u)-lapse_term)
    electric_gradient = _clean(sp.diff(-5/(2*h), u)-lapse_term)
    return {"u": u, "h": h, "Hubble": hubble, "a": previous["a"],
            "gradient_amplitude": amplitude, "Vstar_lapse_value": lapse_term,
            "mass_time_shift": mass_time_shift, "mass_spatial_shift": 1/h,
            "Vstar_curl_lapse_dot": -1/h, "Vstar_curl_lapse_value": electric_gradient,
            "physical_center_mass_squared_times_zeta": sp.Rational(8, 3),
            "full_spectrum_or_degeneracy_claim": False}


@cache
def checks():
    """Exact action, projection, background, and full first-variation identities."""
    data = schur()
    m, n, lift = data["M"], data["N"], data["lift"]
    projector, dinv = data["complement_projector"], data["D_inverse"]
    d0 = 3*(P**3+2*P**2+P-1)/(2*P)
    di = (4*P-5)*(4*P-1)/(32*P**2)
    actual = rolling()
    result = {
        "projective_invariance": vector_map()["projective_residual"],
        "literal_inverse_N": _matrix(m*data["inverse_N_transpose"]-n.T),
        "D_formula": _matrix(data["D"]-sp.diag(d0, di, di, di)),
        "D_inverse": _matrix(data["D"]*dinv-sp.eye(4)),
        "fixed_vector_right_inverse": _matrix(n*lift-sp.eye(4)),
        "constrained_Euler": _matrix(m*lift-n.T*dinv),
        "complement_kernel": _matrix(n*projector),
        "complement_idempotent": _matrix(projector*projector-projector),
        "complement_rank": _clean(sp.trace(projector)-56),
        "Schur_action": _matrix(lift.T*m*lift-dinv),
        "action_cross_term": _matrix(projector.T*m*lift),
        "covariant_reconstruction": rest_coefficients()["reconstruction_residual"],
        "background_source": background()["source"],
        "background_connection": background()["stationary"],
        "background_qx": background()["qx_clock"],
        "background_p": background()["p_clock"]-sp.Rational(1, 2),
        "background_c": background()["c_clock"],
        "background_quartic": background()["quartic_clock"],
        "background_p_phi": background()["p_phi_clock"],
        "background_f_phi": background()["f_phi_clock"],
        "ADM_inverse": adm_linearization()["metric_inverse_residual"],
        "rolling_first_variation": adm_linearization()["first_variation_residual"],
        "rolling_clock_relation": _clean(sp.diff(actual["h"], actual["u"])
                                          -3*actual["Hubble"]*actual["h"]/2),
        "time_mass_shift": _clean(actual["mass_time_shift"]-21*actual["Hubble"]/(8*actual["h"])),
        "curl_value_shift": _clean(actual["Vstar_curl_lapse_value"]-33*actual["Hubble"]/(8*actual["h"])),
        "center_D": _matrix(data["D"].subs(P, sp.Rational(1, 2))
                             -sp.diag(3, -3, -3, -3)/8),
    }
    result.update({f"CD_coefficient_{name}": value
                   for name, value in cd_coefficients()["residuals"].items()})
    return result


@cache
def calibration():
    """Exact bounds and scope flags for a separate read-only report wrapper."""
    return {"quotient_dimension": 60, "vector_dimension": 4, "remaining_auxiliary_dimension": 56,
            "D": schur()["D"], "D_inverse": schur()["D_inverse"],
            "bounds": domain_bounds(), "first_variation": adm_linearization()["first_variation_CD"],
            "center_mass_squared_times_zeta_over_M_squared": sp.Rational(8, 3),
            "rolling_background_preserved": True, "physical_metric_and_chi_unchanged": True,
            "point_transformation_removes_lapse_velocity": False,
            "exact_gradient_linear_transformation_available": True,
            "full_coupled_health_gap_or_UV_claim": False}
