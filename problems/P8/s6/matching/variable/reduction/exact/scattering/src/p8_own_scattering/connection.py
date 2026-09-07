"""Own-f coefficient-16 connection and a conditional finite-window transfer.

This classical homogeneous ODE calculation uses the actual canonical tensor
variable, not the coefficient-80 coupled-relative mode.  The artificial free
exterior defines endpoint coordinates; it is not a physical vacuum extension.
The fixed-window conclusion requires the separately certified |R_delta| <= 44.
"""

from fractions import Fraction
from functools import cache

import sympy as sp

RHO = sp.sqrt(7) / 2
PHASE_SCALE = 4 * sp.sqrt(2)
HALF_WIDTH = Fraction(1, 100)
DELTA_MAX = Fraction(1, 10**6)
REMAINDER_BOUND = Fraction(44)


def _exact_real(value, name, *, positive=False):
    """Validate before any cached use; binary floats and booleans are excluded."""
    if isinstance(value, bool) or value is sp.true or value is sp.false:
        raise TypeError(f"{name} must be exact, not bool")
    value = sp.sympify(value)
    if not isinstance(value, sp.Expr):
        raise TypeError(f"{name} must be a scalar expression")
    if value.has(sp.Float) or value.is_finite is not True or value.is_real is not True:
        raise ValueError(f"{name} must be an exact finite real expression")
    if positive and value.is_positive is not True:
        raise ValueError(f"{name} must be strictly positive")
    return value


@cache
def coefficients():
    """Right Jost solution and its LEFT time-frequency coefficients.

    psi_R ~ exp(i*rho*t) at +infinity, and
    psi_R ~ A*exp(i*rho*t)+B*exp(-i*rho*t) at -infinity.
    'Reflected flux' is the corresponding mathematical scattering convention,
    not a quantum-production or physical spatial-scattering assertion.
    """
    t = sp.symbols("t", real=True)
    w = (1 - sp.tanh(t)) / 2
    aa = (
        sp.gamma(1 - sp.I * RHO) * sp.gamma(-sp.I * RHO)
        / (sp.gamma(sp.Rational(3, 2) - sp.I * RHO)
           * sp.gamma(-sp.Rational(1, 2) - sp.I * RHO))
    )
    bb = sp.I / sp.sinh(sp.pi * RHO)
    return {
        "rho": RHO,
        "t": t,
        "right_Jost": sp.exp(sp.I * RHO * t)
        * sp.hyper((-sp.Rational(1, 2), sp.Rational(3, 2)), (1 - sp.I * RHO,), w),
        "A": aa,
        "B": bb,
        "A_abs_squared": sp.coth(sp.pi * RHO)**2,
        "B_abs_squared": sp.csch(sp.pi * RHO)**2,
        "transmitted_flux": sp.tanh(sp.pi * RHO)**2,
        "reflected_flux": sp.sech(sp.pi * RHO)**2,
        "radial_singular_max": sp.coth(sp.pi * RHO / 2),
        "radial_singular_min": sp.tanh(sp.pi * RHO / 2),
    }


@cache
def even_odd():
    """Ordinary Gauss 2F1 solutions, with W_x(even, odd)=1.

    E_plus and O_plus multiply x**(1/2+i*rho) at x>0.  Conjugates
    multiply the negative log-frequency power.  At x<0 parity leaves
    E_plus unchanged and negates O_plus.
    """
    x = sp.symbols("x", real=True)
    aa = -sp.Rational(1, 4) + sp.I * RHO / 2
    bb = sp.conjugate(aa)
    ee = (
        sp.sqrt(sp.pi) * sp.gamma(sp.I * RHO)
        * 8**(sp.Rational(1, 4) + sp.I * RHO / 2)
        / (sp.gamma(aa) * sp.gamma(sp.Rational(3, 4) + sp.I * RHO / 2))
    )
    oo = (
        sp.gamma(sp.Rational(3, 2)) * sp.gamma(sp.I * RHO)
        * 8**(-sp.Rational(1, 4) + sp.I * RHO / 2)
        / (sp.gamma(sp.Rational(1, 4) + sp.I * RHO / 2)
           * sp.gamma(sp.Rational(5, 4) + sp.I * RHO / 2))
    )
    return {
        "x": x, "a": aa, "b": bb,
        "even": sp.hyper((aa, bb), (sp.Rational(1, 2),), -8 * x**2),
        "odd": x * sp.hyper((aa + sp.Rational(1, 2), bb + sp.Rational(1, 2)),
                            (sp.Rational(3, 2),), -8 * x**2),
        "E_plus": ee, "O_plus": oo,
        "p_plus": sp.Rational(1, 2) + sp.I * RHO,
        "p_minus": sp.Rational(1, 2) - sp.I * RHO,
        "Wronskian_x": sp.Integer(1),
    }


def time_transfer():
    """LEFT -> RIGHT coefficients of exp(+/-i*rho*t); determinant +1."""
    data = coefficients()
    aa, bb = data["A"], data["B"]
    return sp.Matrix([[sp.conjugate(aa), -sp.conjugate(bb)], [-bb, aa]])


def power_transfer(delta):
    """LEFT -> RIGHT |u|**(1/2 +/- i*rho) coefficients; determinant -1.

    Here x=u/sqrt(delta).  The minus determinant comes solely from opposite
    left/right radial orientations; it is not a canonical Cauchy determinant.
    """
    delta = _exact_real(delta, "delta", positive=True)
    data = coefficients()
    phase = sp.exp(2 * sp.I * RHO * sp.log(PHASE_SCALE / sp.sqrt(delta)))
    return sp.Matrix([
        [-sp.conjugate(data["B"]), sp.conjugate(data["A"]) * phase],
        [data["A"] / phase, -data["B"]],
    ])


def unitary_power_coordinates():
    """Real cosine/sine amplitudes -> normalized conjugate power amplitudes."""
    return sp.Matrix([[1, -sp.I], [1, sp.I]]) / sp.sqrt(2)


def real_transfer(delta):
    unitary = unitary_power_coordinates()
    return unitary.conjugate().T * power_transfer(delta) * unitary


def plane_wave_frame(t):
    """Unitary map from normalized wave coefficients to (psi, psi_t/rho)."""
    t = _exact_real(t, "t")
    plus, minus = sp.exp(sp.I * RHO * t), sp.exp(-sp.I * RHO * t)
    return sp.Matrix([[plus, minus], [sp.I * plus, -sp.I * minus]]) / sp.sqrt(2)


def canonical_endpoint_map(u, delta):
    """(Y,Y_u) -> (psi,psi_t/rho), Y=sqrt(r)*psi, r=sqrt(u^2+delta/8)."""
    u = _exact_real(u, "u")
    delta = _exact_real(delta, "delta", positive=True)
    radius = sp.sqrt(u**2 + delta / 8)
    return sp.Matrix([
        [1 / sp.sqrt(radius), 0],
        [-u / (2 * RHO * radius**sp.Rational(3, 2)), sp.sqrt(radius) / RHO],
    ])


def endpoint_map(u, delta, k, k_u):
    """Actual (Q,Q_u) -> (psi,psi_t/rho), with Y=sqrt(k)*Q.

    The input derivative is Q_u=tau*Q_T, not Q_T.  The actual kinetic
    coefficient is k=b^3/N>0; k_u MUST be retained.  Overall constant M
    normalization cancels from this classical transfer.
    """
    k = _exact_real(k, "k", positive=True)
    k_u = _exact_real(k_u, "k_u")
    whitening = sp.Matrix([[sp.sqrt(k), 0], [k_u / (2 * sp.sqrt(k)), sp.sqrt(k)]])
    return canonical_endpoint_map(u, delta) * whitening


@cache
def identities():
    """Exact pullbacks and algebra used by the stated analytic proof.

    Hypergeometric connection identities themselves are source-audited in
    notes/sources.md; these residuals check their parameter dictionary and
    normalization.  No sampled ODE solution is used as a proof of a bound.
    """
    output = {}

    def add_matrix(name, matrix):
        for row in range(matrix.rows):
            for col in range(matrix.cols):
                output[f"{name}_{row}{col}"] = sp.simplify(matrix[row, col])

    t = sp.symbols("t", real=True)
    psi = sp.Function("psi")(t)
    yy = sp.sqrt(sp.cosh(t)) * psi
    pulled = (sp.diff(yy, t, 2) - sp.tanh(t) * sp.diff(yy, t) + 2 * yy)
    expected = sp.sqrt(sp.cosh(t)) * (
        sp.diff(psi, t, 2) + (sp.Rational(7, 4) + sp.Rational(3, 4) / sp.cosh(t)**2) * psi
    )
    output["coefficient16_Liouville"] = sp.simplify(sp.expand(pulled - expected))
    output["rho_squared"] = RHO**2 - sp.Rational(7, 4)

    w, mu = sp.symbols("w mu", positive=True)
    f, fp, fpp = sp.symbols("F Fp Fpp")
    wp = -2 * w * (1 - w)
    wpp = sp.diff(wp, w) * wp
    transformed = wp**2 * fpp + (wpp + 2 * sp.I * mu * wp) * fp + 3 * w * (1 - w) * f
    hypergeometric = w * (1 - w) * fpp + (1 - sp.I * mu - 2 * w) * fp + sp.Rational(3, 4) * f
    output["right_Jost_hypergeometric"] = sp.expand(transformed - 4 * w * (1 - w) * hypergeometric)
    aa, bb = even_odd()["a"], even_odd()["b"]
    output["even_parameter_sum"] = sp.simplify(aa + bb + sp.Rational(1, 2))
    output["even_parameter_product"] = sp.simplify(aa * bb - sp.Rational(1, 2))
    output["odd_parameter_sum"] = sp.simplify(aa + bb + 1 - sp.Rational(1, 2))
    output["odd_parameter_product"] = sp.simplify((aa + sp.Rational(1, 2)) * (bb + sp.Rational(1, 2)) - sp.Rational(1, 2))
    output["even_ODE_pullback"] = sp.cancel(
        (-16 * fp - 32 * w * fpp + 16 * f / (1 - w)) * (1 - w)
        + 32 * (w * (1 - w) * fpp + (1 - w) * fp / 2 - f / 2)
    )
    output["odd_ODE_pullback"] = sp.cancel(
        (-48 * fp - 32 * w * fpp + 16 * f / (1 - w)) * (1 - w)
        + 32 * (w * (1 - w) * fpp + 3 * (1 - w) * fp / 2 - f / 2)
    )
    output["gamma_B_reflection"] = sp.simplify(
        sp.pi / sp.sin(sp.pi * sp.I * mu) / (-sp.pi) - sp.I / sp.sinh(sp.pi * mu)
    )
    gamma_i_sq = sp.pi / (mu * sp.sinh(sp.pi * mu))
    gamma_half_sq = sp.pi / sp.cosh(sp.pi * mu)
    denominator = (mu**2 + sp.Rational(1, 4)) * gamma_half_sq**2 / (mu**2 + sp.Rational(1, 4))
    output["gamma_A_modulus_recurrence"] = sp.simplify(mu**2 * gamma_i_sq**2 / denominator - sp.coth(sp.pi * mu)**2)
    output["Wronskian_modulus"] = sp.simplify(sp.coth(sp.pi * mu)**2 - sp.csch(sp.pi * mu)**2 - 1)
    output["reflected_flux"] = sp.simplify(sp.csch(sp.pi * mu)**2 / sp.coth(sp.pi * mu)**2 - sp.sech(sp.pi * mu)**2)

    a, ac, b, bc, phase = sp.symbols("A Abar B Bbar phase", nonzero=True)
    forward = sp.Matrix([[ac, -bc], [-b, a]])
    backward = sp.Matrix([[a, bc], [b, ac]])
    add_matrix("time_inverse", forward * backward - (a * ac - b * bc) * sp.eye(2))
    radial = sp.Matrix([[-bc, ac * phase], [a / phase, -b]])
    output["time_determinant"] = sp.expand(forward.det() - a * ac + b * bc)
    output["radial_determinant"] = sp.expand(radial.det() + a * ac - b * bc)
    # The even-potential involution uses the already derived pure-imaginary B.
    add_matrix("radial_involution", (radial**2 - (a * ac - b * bc) * sp.eye(2)).subs(bc, -b))
    ar, ai, beta = sp.symbols("ar ai beta", real=True)
    unitary = unitary_power_coordinates()
    radial_real = sp.Matrix([[ar, ai + beta], [ai - beta, -ar]])
    radial_form = sp.Matrix([[sp.I * beta, ar - sp.I * ai], [ar + sp.I * ai, -sp.I * beta]])
    add_matrix("real_radial_formula", unitary.conjugate().T * radial_form * unitary - radial_real)
    add_matrix("power_unitarity", unitary.conjugate().T * unitary - sp.eye(2))
    waves = plane_wave_frame(t)
    add_matrix("wave_unitarity", waves.conjugate().T * waves - sp.eye(2))
    output["wave_frame_determinant"] = sp.simplify(waves.det() + sp.I)

    u = sp.symbols("u", real=True)
    delta, kinetic = sp.symbols("delta k", positive=True)
    kinetic_u = sp.symbols("k_u", real=True)
    radius = sp.sqrt(u**2 + delta / 8)
    output["asinh_clock_derivative"] = sp.simplify(sp.diff(sp.asinh(sp.sqrt(8) * u / sp.sqrt(delta)), u) - 1 / radius)
    output["canonical_map_determinant"] = sp.simplify(canonical_endpoint_map(u, delta).det() - 1 / RHO)
    output["physical_map_determinant"] = sp.simplify(endpoint_map(u, delta, kinetic, kinetic_u).det() - kinetic / RHO)
    eps, scale = sp.symbols("epsilon scale", positive=True)
    output["physical_log_phase"] = sp.simplify(
        sp.exp(2 * sp.I * mu * (sp.log(scale) - sp.log(eps)))
        - sp.exp(2 * sp.I * mu * sp.expand_log(sp.log(scale / eps), force=True))
    )

    z1, z2, potential = sp.symbols("z1 z2 potential", real=True)
    generator = sp.Matrix([[0, mu], [-mu - potential / mu, 0]])
    state = sp.Matrix([z1, z2])
    output["energy_symmetric_factor"] = sp.expand((state.T * (generator + generator.T) * state)[0] + 2 * potential * z1 * z2 / mu)
    length, rr = sp.symbols("L radius", positive=True)
    output["tail_rationalization"] = sp.factor(
        (1 - length / rr) - (rr**2 - length**2) / (rr * (rr + length))
    )
    output["tail_difference_numerator"] = sp.expand(
        (rr**2 - length**2) - (rr - length) * (rr + length)
    )
    x = sp.symbols("x", real=True)
    quotient = x**6 - 4 * x**5 + 5 * x**4 - 4 * x**2 + 4
    output["pi_upper_positive_integral_division"] = sp.expand(x**4 * (1 - x)**4 - (1 + x**2) * quotient + 4)
    output["pi_upper_polynomial_integral"] = sp.integrate(quotient, (x, 0, 1)) - sp.Rational(22, 7)
    return output


def calibration():
    """Exact sufficient constants; |R_delta|<=44 is an external checked premise."""
    central = REMAINDER_BOUND * (HALF_WIDTH**2 + 2 * HALF_WIDTH * Fraction(1, 2000))
    tail = 3 * DELTA_MAX / (32 * HALF_WIDTH**2)
    budget = (central + tail) * Fraction(4, 5)
    return {
        "half_width": HALF_WIDTH,
        "delta_max": DELTA_MAX,
        "remainder_bound": REMAINDER_BOUND,
        "small_radius_upper": Fraction(1, 2000),
        "central_integral_upper": central,
        "tail_integral_upper": tail,
        "perturbation_budget_upper": budget,
        "budget_gate": Fraction(1, 200),
        "reference_flow_upper": Fraction(2),
        "transfer_error_upper": Fraction(4, 399),
        "transfer_error_gate": Fraction(1, 80),
        "reference_B_lower": Fraction(1, 40),
        "effective_B_lower": Fraction(1, 40) - Fraction(4, 399),
        "effective_B_gate": Fraction(1, 80),
        "exp_one_upper": Fraction(49, 18),
        "exp_three_fifths_upper": Fraction(1549, 850),
        "exp_seventeen_fourths_upper": Fraction(14641, 192),
        "conditional_premise": "actual canonical remainder |R_delta(u)|<=44 on the full fixed window",
        "endpoint_input": "(Q,Q_u)=(Q,tau*Q_T), with actual k=b^3/N and k_u",
        "exterior_scope": "fictitious free exterior only; no physical vacuum extension",
        "data_scope": "homogeneous input transfer, not a zero-data matter/source response",
    }


def checks():
    """Rational strict margins, separate from symbolic residual identities."""
    data = calibration()
    f = Fraction
    return {
        "rho_above_5_over_4": f(7, 4) > f(25, 16),
        "rho_below_4_over_3": f(7, 4) < f(16, 9),
        "pi_rho_below_17_over_4": f(22, 7) * f(4, 3) < f(17, 4),
        "exp_one_below_11_over_4": data["exp_one_upper"] < f(11, 4),
        "exp_product_bound": f(11, 4)**4 * f(4, 3) == data["exp_seventeen_fourths_upper"],
        "exp_product_below_80": data["exp_seventeen_fourths_upper"] < 80,
        "reference_energy_factor_below_two": data["exp_three_fifths_upper"] < 2,
        "small_radius_strict": DELTA_MAX / 8 < data["small_radius_upper"]**2,
        "central_integral_value": data["central_integral_upper"] == f(121, 25000),
        "tail_integral_value": data["tail_integral_upper"] == f(3, 3200),
        "budget_value": data["perturbation_budget_upper"] == f(2311, 500000),
        "budget_below_gate": data["perturbation_budget_upper"] < data["budget_gate"],
        "exp_geometric_error": 4 * (1 / (1 - data["budget_gate"] / 2) - 1) == data["transfer_error_upper"],
        "error_below_gate": data["transfer_error_upper"] < data["transfer_error_gate"],
        "effective_B_above_gate": data["effective_B_lower"] > data["effective_B_gate"],
        "actual_delta_subset": DELTA_MAX < f(1, 100),
    }
