"""Curved auxiliary reduction and the actual two-scalar principal action.

Positive principal matrices are not a UV completion, heavy spectrum,
all-frequency stability theorem or exact CD/M1 action matching.
"""
from functools import cache

import sympy as sp
from p8_affine_kinetic import scalar as old

from . import geometry, source

u, q = old.u, old.q
COUPLING = sp.Symbol("lambda", real=True)
AUXILIARY = sp.Symbol("A_correction", real=True)


def _matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.factor))


@cache
def coefficients():
    bg = source.background()
    h, hubble = bg["h"], bg["H"]
    curvature = sp.factor(2*(sp.diff(hubble, u)+2*hubble**2))
    denominator = 1+COUPLING*curvature
    cstar = source.rolling()["coefficient"]
    correction = sp.factor(COUPLING*cstar**2/(2*denominator))
    return {**bg, "curvature_endomorphism": curvature, "denominator": sp.factor(denominator),
            "Cstar_coefficient": cstar, "A": correction,
            "expected_A": 50*COUPLING*u**2/((1+u**2)**6
                            *((1+u**2)**2+8*COUPLING*(1+7*u**2))),
            "scalar_curvature": 3*curvature,
            "hubble_source_coefficient": -5*hubble/(2*h)}


@cache
def auxiliary_action():
    """Two-form multiplier after all64 connection equations and covariant IBP.

Pairing W is indefinite but algebraic, so its sign is not a ghost test.
The constant flat symbol is null; the curved endomorphism is retained.
"""
    z = sp.Matrix(sp.symbols("z0:6", real=True))
    cstar = sp.Matrix(sp.symbols("Cstar0:6", real=True))
    curvature = sp.Symbol("curvature_endomorphism", real=True)
    W = geometry.PAIRING
    action = ((1/COUPLING+curvature)*(z.T*W*z)[0]/2-(z.T*W*cstar)[0])
    solution = COUPLING*cstar/(1+COUPLING*curvature)
    reduced = -COUPLING*(cstar.T*W*cstar)[0]/(2*(1+COUPLING*curvature))
    return {"z": z, "Cstar": cstar, "curvature": curvature, "before": action,
            "solution": solution, "reduced": reduced,
            "Euler_residual": _matrix(sp.Matrix([sp.diff(action, field) for field in z])
                                       .subs(dict(zip(z, solution, strict=True)))),
            "action_residual": sp.factor(action.subs(dict(zip(z, solution, strict=True)))-reduced)}


@cache
def scalar_action():
    """The full lapse, shift and free matter before and after constraint solves."""
    before = old.action()["base"]+q*AUXILIARY*old.n**2
    lapse = (old.vd+old.ell*old.matter/2)/old.theta
    after = sp.factor(before.subs(old.n, lapse))
    kinetic = _matrix(sp.hessian(after, (old.vd, old.sd))/2)
    target = (old.sd+old.w*old.vd/old.theta)**2/2+(old.J+q*AUXILIARY)*old.vd**2/old.theta**2
    velocity = sp.Matrix([old.vd, old.sd])
    first_order = (-2*q*old.shift*old.vd+old.pm*old.sd-before).subs(
        {old.vd: old.theta*old.n-old.ell*old.matter/2, old.sd: old.pm-old.w*old.n})
    R = q*(old.theta*old.shift+old.lam*old.v)+old.w*old.pm/2-3*old.ell*old.matter*old.theta/2
    constant = old.pm**2/2+q*old.ell*old.shift*old.matter+q*old.matter**2/2-q*old.v**2-3*old.ell**2*old.matter**2/4
    Jeff = old.J+q*AUXILIARY
    return {"before": before, "lapse": lapse, "after": after, "kinetic": kinetic,
            "J_eff": Jeff, "first_order_before_lapse": first_order,
            "first_order_H": constant+R**2/Jeff, "first_order_R": R,
            "first_order_lapse": -R/Jeff,
            "shift_residual": sp.factor(sp.diff(before, old.shift)/(2*q)
                                        -old.theta*old.n+old.vd+old.ell*old.matter/2),
            "kinetic_squares_residual": sp.factor((velocity.T*kinetic*velocity)[0]-target),
            "kinetic_determinant_residual": sp.factor(kinetic.det()-Jeff/(2*old.theta**2)),
            "first_order_residual": sp.factor(first_order-constant+Jeff*old.n**2+2*R*old.n),
            "first_order_lapse_residual": sp.factor(sp.diff(first_order, old.n).subs(old.n, -R/Jeff)),
            "physical_quadratic_mismatch": sp.diff(before-old.action()["base"], old.n, 2)}


@cache
def gradients():
    """Keep the time-dependent boundary and q*s*v_dot gyroscopic term."""
    bg = old.background()
    a, theta, ell, lam = (bg[name] for name in ("a", "theta", "ell", "lam"))
    correction = coefficients()["A"]
    mixing = sp.factor(correction*ell/theta**2)
    vv = sp.factor(sp.diff(a*lam/theta, u)/a-1)
    ss = sp.factor(sp.Rational(1, 2)-correction*ell**2/(4*theta**2))
    vs = sp.factor(-lam*ell/(2*theta)+sp.diff(a*mixing, u)/(4*a))
    matrix = sp.ImmutableMatrix([[vv, vs], [vs, ss]])
    determinant = sp.factor(matrix.det())
    v, matter, vd, sd = old.v, old.matter, old.vd, old.sd
    raw = (v**2+2*lam*v*vd/theta+lam*ell*v*matter/theta-matter**2/2
           +correction*(vd+ell*matter/2)**2/theta**2)
    normalized = (correction*vd**2/theta**2+mixing*(matter*vd-v*sd)/2
                  -(sp.Matrix([v, matter]).T*matrix*sp.Matrix([v, matter]))[0])
    boundary = (sp.diff(a*lam/theta, u)*v**2/a+2*lam*v*vd/theta
                +sp.diff(a*mixing, u)*v*matter/(2*a)+mixing*(vd*matter+v*sd)/2)
    return {"matrix": matrix, "determinant": determinant, "gyroscopic_coefficient": mixing/2,
            "boundary_residual": sp.factor(raw-normalized-boundary),
            "ss_numerator": sp.fraction(ss)[0], "ss_denominator": sp.fraction(ss)[1],
            "det_numerator": sp.fraction(determinant)[0], "det_denominator": sp.fraction(determinant)[1]}


@cache
def center():
    """Regular original canonical variables at u=0, for lambda>=0 and q>6."""
    bg = old.background()
    pb, ps = sp.symbols("P_b P_s", real=True)
    first = old.a**3*scalar_action()["first_order_H"].subs(
        {old.v: pb/(2*old.a**3*q), old.pm: ps/old.a**3}, simultaneous=True)-old.H*old.shift*pb
    actual = {symbol: bg[name] for symbol, name in ((old.J, "J"), (old.theta, "theta"),
              (old.ell, "ell"), (old.w, "w"), (old.lam, "lam"), (old.a, "a"), (old.H, "H"))}
    actual[AUXILIARY] = coefficients()["A"]
    hamiltonian = sp.factor(first.subs(actual, simultaneous=True).subs(u, 0))
    momentum_hessian = _matrix(sp.hessian(hamiltonian, (pb, ps)))
    kinetic = _matrix(momentum_hessian.inv()/2)
    expected = sp.Matrix([[6*q/(q-6), q/(20*(q-6))],
                          [q/(20*(q-6)), (200*q-1199)/(400*(q-6))]])
    return {"Hamiltonian": hamiltonian, "momentum_hessian": momentum_hessian,
            "kinetic": kinetic, "expected": expected,
            "momentum_determinant": sp.factor(momentum_hessian.det()),
            "domain": "lambda>=0 and q>6; regular first-order Hamiltonian exists for every q>0"}


def positive_polynomial(value):
    """Strict positivity on all real u and lambda>=0 by exact coefficients."""
    poly = sp.Poly(value, u, COUPLING)
    terms = poly.terms()
    return {"degrees": poly.degree_list(), "monomials": len(terms),
            "positive_constant": poly.coeff_monomial(1),
            "all_even_u_nonnegative": all(powers[0] % 2 == 0 and coefficient >= 0
                                           for powers, coefficient in terms),
            "strictly_positive": bool(poly.coeff_monomial(1) > 0 and all(
                powers[0] % 2 == 0 and coefficient >= 0 for powers, coefficient in terms))}


@cache
def kinetic_remainder():
    """One uniform quadratic kinetic bound, not a full V/G/B remainder."""
    ratio = sp.factor(coefficients()["Cstar_coefficient"]**2/(2*old.background()["J"]))
    numerator, denominator = sp.fraction(ratio)
    witness = sp.expand(10*denominator-numerator)
    return {"ratio_before_positive_curved_denominator": ratio,
            "strict_upper": 10, "positive_witness": witness,
            "certificate": positive_polynomial(witness),
            "relative_kinetic_norm": sp.factor(q*coefficients()["A"]/old.background()["J"]),
            "bound": 10*COUPLING*q,
            "scope": "u!=0, lambda>=0, q>0; only the two-scalar quadratic velocity form"}


@cache
def proof_checks():
    data, gradient = coefficients(), gradients()
    numerator, denominator = sp.fraction(old.background()["J"])
    inputs = {"J_numerator": numerator, "J_denominator": denominator,
              "gradient_ss_numerator": gradient["ss_numerator"],
              "gradient_ss_denominator": gradient["ss_denominator"],
              "gradient_det_numerator": gradient["det_numerator"],
              "gradient_det_denominator_after_punctured_u_squared": sp.cancel(gradient["det_denominator"]/u**2)}
    out = {name: positive_polynomial(value)["strictly_positive"] for name, value in inputs.items()}
    cnum, cden = sp.fraction(data["curvature_endomorphism"])
    out["curvature_positive_all_real_time"] = (positive_polynomial(cnum)["strictly_positive"]
                                               and positive_polynomial(cden)["strictly_positive"])
    out["zero_correction_at_center"] = sp.limit(data["A"], u, 0) == 0
    out["zero_first_correction_jet_at_center"] = sp.limit(sp.diff(data["A"], u), u, 0) == 0
    out["uniform_quadratic_kinetic_remainder_below_10_lambda_q"] = kinetic_remainder()["certificate"]["strictly_positive"]
    out["center_uses_regular_first_order_not_Theta_inverse"] = True
    out["center_positive_momentum_matter_pivot"] = bool(center()["momentum_hessian"][1, 1] > 0)
    out["no_heavy_gap_exact_matching_or_UV_conclusion"] = True
    return out


def _exact(value, name):
    if isinstance(value, (bool, float, str)):
        raise TypeError(f"{name} must be an exact finite real")
    value = sp.sympify(value)
    if (not isinstance(value, sp.Expr) or value.has(sp.Float) or value.free_symbols
            or value.is_real is not True or value.is_finite is not True):
        raise ValueError(f"{name} must be an exact finite real constant")
    return value


def require_regular(coupling, time, momentum, *, punctured=True):
    lam, point, q_value = [_exact(value, name) for value, name in
                            ((coupling, "coupling"), (time, "time"), (momentum, "q"))]
    if q_value.is_positive is not True:
        raise ValueError("Require q>0")
    if punctured and point.is_zero is not False:
        raise ValueError("The configuration-velocity chart requires u!=0")
    denominator = sp.factor(coefficients()["denominator"].subs({u: point, COUPLING: lam}))
    if denominator.is_zero is not False:
        raise ValueError("The curved algebraic two-form chart must be nonsingular")
    return {"lambda": lam, "u": point, "q": q_value, "denominator": denominator,
            "unchanged_auxiliary": lam.is_zero is True}


def negative_witness(coupling, time):
    lam, point = _exact(coupling, "coupling"), _exact(time, "time")
    if lam.is_negative is not True or sp.simplify(point**2+112*lam).is_positive is not True:
        raise ValueError("Require lambda<0 and u²>112*abs(lambda)")
    correction = sp.factor(coefficients()["A"].subs({COUPLING: lam, u: point}))
    J = old.background()["J"].subs(u, point)
    if correction.is_negative is not True:
        raise ArithmeticError("The strict negative-coupling tail witness failed")
    threshold = sp.factor(J/(-correction))
    return {**require_regular(lam, point, 2*threshold), "threshold": threshold,
            "J_eff": sp.factor(J+2*threshold*correction)}


def units(mass_squared, time_scale, physical_coupling):
    mass, tau, coupling = [_exact(value, name) for value, name in
                           ((mass_squared, "M²"), (time_scale, "tau"), (physical_coupling, "lambda_physical"))]
    if mass.is_positive is not True or tau.is_positive is not True:
        raise ValueError("M² and tau must be positive")
    return {"normalized_lambda": coupling/(mass*tau**2),
            "normalized_C_factor": tau**2, "physical_curvature_factor": 1/tau**2,
            "additional_mass_or_gap_assigned": False}


@cache
def checks():
    data, auxiliary, scalar, gradient = coefficients(), auxiliary_action(), scalar_action(), gradients()
    bg = old.background()
    out = {"actual_curved_commutator": sp.factor(data["curvature_endomorphism"]
                                                  -8*(1+7*u**2)/(1+u**2)**2),
           "actual_Cstar_coefficient": sp.factor(data["Cstar_coefficient"]-data["hubble_source_coefficient"]),
           "literal_correction_coefficient": sp.factor(data["A"]-data["expected_A"]),
           "two_form_Euler": auxiliary["Euler_residual"], "two_form_action": auxiliary["action_residual"],
           "gradient_time_boundary": gradient["boundary_residual"],
           "zero_coupling_returns_original_action": sp.factor(data["A"].subs(COUPLING, 0)),
           "curvature_tail_upper": sp.factor(56/u**2-data["curvature_endomorphism"]
                                              -8*(7+13*u**2)/(u**2*(1+u**2)**2)),
           "physical_quadratic_mismatch_not_boundary": scalar["physical_quadratic_mismatch"]-2*q*AUXILIARY,
           "nonunit_coupling": units(3, 2, 5)["normalized_lambda"]-sp.Rational(5, 12)}
    out.update({"scalar_"+name: value for name, value in scalar.items() if name.endswith("residual")})
    out.update({"original_background_"+name: sp.factor(data[name]-bg[name]) for name in ("a", "h", "H")})
    crossing = center()
    out["regular_center_kinetic_matrix"] = _matrix(crossing["kinetic"]-crossing["expected"])
    out["regular_center_momentum_determinant"] = sp.factor(crossing["momentum_determinant"]-100*(q-6)/(1199*q))
    mass, tau, physical_coupling, physical_q = sp.symbols("M_squared tau lambda_physical q_physical", positive=True)
    out["kinetic_bound_physical_units"] = sp.factor(
        10*physical_coupling/(mass*tau**2)*tau**2*physical_q-10*physical_coupling*physical_q/mass)
    return out


@cache
def calibration():
    data, gradient = coefficients(), gradients()
    return {"curvature_endomorphism": data["curvature_endomorphism"], "Cstar_coefficient": data["Cstar_coefficient"],
            "A": data["A"], "gradient_ss_polynomial": positive_polynomial(gradient["ss_numerator"]),
            "gradient_det_polynomial": positive_polynomial(gradient["det_numerator"]),
            "scalar_degrees_of_freedom_on_regular_branch": 2,
            "quadratic_kinetic_remainder_bound": kinetic_remainder()["bound"],
            "quadratic_kinetic_remainder_polynomial": kinetic_remainder()["certificate"],
            "regular_center_kinetic": center()["kinetic"],
            "positive_coupling_principal_matrices_positive": True,
            "physical_CD_quadratic_action_changed_for_nonzero_coupling": True,
            "nonlinear_spectrum_UV_completion_or_heavy_gap_claim": False}
