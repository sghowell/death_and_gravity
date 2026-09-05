"""Exact CD curvature/variation and conditional local-coefficient estimates.

P8 uses +--- and R=-6*(Hdot+2*H**2).  All metric variations are performed
before restricting to FLRW.  Components below are lower orthonormal ones,
not mixed components.  The source denominator M**2/ell**2 never vanishes.
"""

from functools import cache

import sympy as sp

from . import heat_kernel as hk

U = sp.Symbol("u", real=True)
TAU = sp.Symbol("tau", positive=True)
X = sp.Symbol("x", real=True)  # Compact time, NOT the clock kinetic invariant.
V = sp.Symbol("v", real=True)  # x**2, in [0,1].
MTAU_S58 = sp.Integer(10)**324


def compact_derivative(expression, dimension):
    """ell**(w+1)*d_t[ell**(-w)*f(x)], ell=tau*sqrt(1+u**2)."""
    return sp.expand((1-X**2)*sp.diff(expression, X)-dimension*X*expression)


def flrw(hubble, dt):
    """Covariant curvature contractions and the already-varied R2 tensor.

    I_ab=2R R_ab-g_ab R2/2+2(g_ab Box-nabla_a nabla_b)R.
    Weyl vanishes as a tensor on flat FLRW, so its full first variation does
    too.  That statement says nothing about its second variation.
    """
    hd = dt(hubble)
    ricci = -6*(hd+2*hubble**2)
    ric2 = 12*(hd**2+3*hubble**2*hd+3*hubble**4)
    riem2 = 12*((hd+hubble**2)**2+hubble**4)
    box_r = dt(dt(ricci))+3*hubble*dt(ricci)
    i00 = -6*ricci*(hd+hubble**2)-ricci**2/2+6*hubble*dt(ricci)
    ip = (2*ricci*(hd+3*hubble**2)+ricci**2/2
          - 2*dt(dt(ricci))-4*hubble*dt(ricci))
    values = {"R": ricci, "Ric2": ric2, "Riem2": riem2,
              "E4": riem2-4*ric2+ricci**2, "C2": riem2-2*ric2+ricci**2/3,
              "Rdot": dt(ricci), "Rddot": dt(dt(ricci)), "BoxR": box_r,
              "I_R2_00": i00, "I_R2_pressure": ip,
              "I_A2_00": i00/72, "I_A2_pressure": ip/72}
    values["A2_chosen_representative"] = (riem2-ric2)/180+ricci**2/72+box_r/30
    return {key: sp.factor(value) for key, value in values.items()}


@cache
def cosmic():
    return flrw(4*U/(TAU*(1+U**2)), lambda expression: sp.diff(expression, U)/TAU)


def compact():
    """Polynomials f(x) in the dimensionful identity invariant=ell**(-w)*f."""
    return {"R": -24*(1+6*X**2),
            "Ric2": 192*(1+8*X**2+28*X**4),
            "Riem2": 192*(1+4*X**2+20*X**4),
            "E4": 1536*X**2*(1+2*X**2), "C2": sp.Integer(0),
            "Rdot": -240*X+576*X**3,
            "Rddot": -240+2688*X**2-3456*X**4,
            "BoxR": -240-192*X**2+3456*X**4,
            "I_R2_00": 288*(1-16*X**2+36*X**4),
            "I_R2_pressure": 576*(1-2*X**2-6*X**4),
            "I_A2_00": 4*(1-16*X**2+36*X**4),
            "I_A2_pressure": 8*(1-2*X**2-6*X**4),
            "A2_chosen_representative": sp.Rational(32, 3)*X**2*(8+37*X**2)}


def dimension(key):
    return {"R": 2, "Rdot": 3}.get(key, 4)


@cache
def algebra_checks():
    direct, polynomials = cosmic(), compact()
    ell = TAU*sp.sqrt(1+U**2)
    checks = {"compact_"+key: sp.simplify(value.subs(X, U/sp.sqrt(1+U**2))/ell**dimension(key)-direct[key])
              for key, value in polynomials.items()}
    checks.update({
        "compact_first_curvature_derivative": sp.expand(compact_derivative(polynomials["R"], 2)-polynomials["Rdot"]),
        "compact_second_curvature_derivative": sp.expand(compact_derivative(polynomials["Rdot"], 3)-polynomials["Rddot"]),
        "variation_trace": sp.expand(polynomials["I_R2_00"]-3*polynomials["I_R2_pressure"]-6*polynomials["BoxR"]),
        "variation_conservation": sp.expand(compact_derivative(polynomials["I_R2_00"], 4)
                                              + 12*X*(polynomials["I_R2_00"]+polynomials["I_R2_pressure"])),
    })
    witness = hk.nonclosure_witness()
    a_cd = (1+U**2)**2
    substitutions = {witness[name]: sp.diff(a_cd, U, j)/TAU**j
                     for j, name in enumerate(("a", "adot", "addot", "a3", "a4"))}
    checks["minisuperspace_scale_variation"] = sp.factor(
        witness["euler_lagrange"].subs(substitutions, simultaneous=True)
        - a_cd**2*direct["I_R2_pressure"]/12)
    mass = sp.Symbol("M", positive=True)
    chi_dot = mass/(10*TAU*(1+U**2)**6)
    checks["free_M1_background_equation"] = sp.diff(a_cd**3*chi_dot, U)/TAU
    return checks


def polynomial_range(expression):
    """Exact quadratic extrema on compact v=x**2 in [0,1], no sampling."""
    polynomial = sp.Poly(sp.expand(expression), X)
    if any(monomial[0] % 2 for monomial, _ in polynomial.terms()):
        raise ValueError("Expected an even compact polynomial")
    in_v = sum(coefficient*V**(monomial[0]//2) for monomial, coefficient in polynomial.terms())
    if sp.degree(in_v, V) > 2:
        raise ValueError("This exact extrema checker is limited to quadratics in x^2")
    locations = {sp.Integer(0), sp.Integer(1)}
    derivative = sp.diff(in_v, V)
    if derivative != 0:
        for point in sp.solve(derivative, V):
            if point.is_Rational is not True:
                raise ValueError("Expected rational stationary locations")
            if 0 <= point <= 1:
                locations.add(point)
    values = [in_v.subs(V, point) for point in sorted(locations)]
    return {"min": min(values), "max": max(values),
            "absolute_sup": max(abs(value) for value in values),
            "candidate_v": sorted(locations)}


def bounds():
    keys = ("R", "Ric2", "Riem2", "E4", "C2", "BoxR", "I_R2_00", "I_R2_pressure",
            "I_A2_00", "I_A2_pressure", "A2_chosen_representative")
    return {key: {"dimension": dimension(key), **polynomial_range(compact()[key])} for key in keys}


def exact_finite_parameter(value):
    """No rounded Float, infinity or unproved-finite symbol enters an enclosure."""
    value = sp.sympify(value)
    if value.has(sp.Float) or value.is_finite is not True:
        raise ValueError("Bounds require exact finite parameters; rounded Float inputs are not enclosures")
    return value


def constant_log_source_bound(mtau, log_abs, hbar=1, rational=True):
    """Conditional |2 Delta cR I_R2|/(M2/ell2), all center times.

    log_abs is an externally supplied upper bound on the CONSTANT scale
    comparison |ln(mu/mu0)|.  It is not mu(t), a fitted finite coefficient,
    or a bound on the complete renormalized effective action/stress tensor.

    Exact rational/algebraic or finite symbolic inputs with the required sign
    assumptions are accepted; machine/decimal Float and nonfinite values are
    rejected because they are not certified outward enclosures.
    """
    mtau, log_abs, hbar = map(exact_finite_parameter, (mtau, log_abs, hbar))
    if mtau.is_positive is not True or log_abs.is_nonnegative is not True or hbar.is_nonnegative is not True:
        raise ValueError("Need positive M*tau and nonnegative log/hbar bounds")
    coefficient = sp.Rational(7, 6) if rational else 168/(4*sp.pi)**2
    return coefficient*hbar*log_abs/mtau**2


def finite_local_source_bound(mtau, c_r_abs):
    """Conditional bound requiring independent matching input |c_R(mu0)|."""
    mtau, c_r_abs = map(exact_finite_parameter, (mtau, c_r_abs))
    if mtau.is_positive is not True or c_r_abs.is_nonnegative is not True:
        raise ValueError("Need positive M*tau and a supplied nonnegative finite-coefficient bound")
    return 12096*c_r_abs/mtau**2


def tensor_principal_checks():
    """A flat TT highest-derivative coefficient test, not a mode/cone bound.

    C2=E4+2 Ric2-2 R2/3, and Ric_lin=-Box h/2 for TT h.  Each real
    polarization has E_ij E_ij=1 and h=2Y/M at a=1.  Only the indicated
    principal coefficients are compared; the Einstein on-shell kernel is
    never a denominator.  Curved scalar/matter/constraint reduction is open.
    """
    omega, momentum, amplitude, mass, field, c_c, c_r = sp.symbols("omega k h M Y cC cR")
    polarization = sp.diag(1, -1, 0)/sp.sqrt(2)
    box_h = -(omega**2-momentum**2)*amplitude*polarization
    ricci_linear = -box_h/2
    c2_quadratic = 2*sp.trace(ricci_linear.T*ricci_linear)
    expected = (omega**2-momentum**2)**2*amplitude**2/2
    canonical = (c_c*c2_quadratic).subs(amplitude, 2*field/mass)
    f = -mass**2*hk.R/2+c_r*hk.R**2
    effective_planck = -2*sp.diff(f, hk.R)
    return {"TT_norm": sp.simplify(sp.trace(polarization.T*polarization)-1),
            "TT_Weyl_quadratic": sp.expand(c2_quadratic-expected),
            "TT_Weyl_canonical_coefficient": sp.expand(canonical-2*c_c*(omega**2-momentum**2)**2*field**2/mass**2),
            "R_squared_two_derivative_tensor_coefficient": sp.expand(effective_planck-mass**2+4*c_r*hk.R)}


def controls():
    t = sp.Symbol("t", positive=True)
    radiation = flrw(1/(2*t), lambda expression: sp.diff(expression, t))
    c = compact()
    wrong_density = (-6*c["R"]*(4+8*X**2)-c["R"]**2/2)
    finite = sp.Symbol("cR_finite", real=True)
    return {"R_flat_radiation_R": radiation["R"],
            "R_flat_radiation_I_R2_00": radiation["I_R2_00"],
            "R_flat_radiation_I_R2_pressure": radiation["I_R2_pressure"],
            "drop_curvature_derivatives_I00_error_at_x_half": sp.expand(wrong_density-c["I_R2_00"]).subs(X, sp.Rational(1, 2)),
            "wrong_R_sign_curvature_error_at_bounce": -2*c["R"].subs(X, 0),
            "drop_ell_derivative_Rdot_error": sp.expand((1-X**2)*sp.diff(c["R"], X)-c["Rdot"]),
            "full_A2_density_bounce": c["A2_chosen_representative"].subs(X, 0),
            "bulk_A2_variation_bounce_not_zero": c["I_A2_00"].subs(X, 0),
            "finite_local_coefficient_remains_free": 2*finite*c["I_R2_00"].subs(X, 0)}


def report():
    result = {}
    for key, record in bounds().items():
        result[key] = {name: list(map(str, value)) if isinstance(value, list) else str(value)
                       for name, value in record.items()}
        result[key]["compact_polynomial"] = str(sp.expand(compact()[key]))
        result[key]["cosmic_expression"] = str(cosmic()[key])
    return {"exact_all_time_invariants_and_variations": result,
            "local_finite_R2_source_ratio": "12096*abs(cR(mu0))/(M*ell)^2; finite cR is matching input, not bounded here",
            "isolated_constant_log_source_ratio": "168*hbar*abs(log(mu/mu0))/((4*pi)^2*(M*ell)^2)",
            "rational_log_source_enclosure": "(7/6)*hbar*L/(M*tau)^2 if abs(log(mu/mu0))<=L; uses pi>3 and ell>=tau",
            "S5_8_Mtau": str(MTAU_S58),
            "conditional_unit_log_hbar_one_source_ratio_upper": str(constant_log_source_bound(MTAU_S58, 1)),
            "tensor_R2_two_derivative_coefficient_ratio": "abs(delta M_T^2)/M^2<=672*abs(cR)/(M*ell)^2 for this added local term",
            "tensor_flat_Weyl_highest_derivative_coefficient": "L2_cC=2*cC*(Box Y)^2/M^2, TT norm one; not a curved finite-band operator bound"}
