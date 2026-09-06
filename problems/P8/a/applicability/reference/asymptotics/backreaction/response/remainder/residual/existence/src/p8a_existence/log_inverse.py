"""Exact bounds for the causal logarithmic inverse, including its positive pole.

The operator is dimensionless and isolated. Its parameter beta is not a
permission to change the physical finite stress prescription. A small bound
for this linear inverse alone is not a semiclassical existence certificate.
"""

import sympy as sp


def exact_nonnegative(value, name="value"):
    if isinstance(value, (bool, float)):
        raise TypeError(name+" must be an exact nonnegative rational")
    value = sp.sympify(value)
    if value.is_Rational is not True or value.is_finite is not True or value < 0:
        raise ValueError(name+" must be an exact nonnegative rational")
    return value


def exact_integer(value, name, minimum=0):
    value = exact_nonnegative(value, name)
    if value.is_Integer is not True or value < minimum:
        raise ValueError(name+" must be an exact integer in its stated domain")
    return int(value)


def dyadic_norm_bound(beta_floor_size, n):
    """C0 norm bound for beta>=-B on L=2**(-2*n), with integer B>=0.

    Uses exp(B)<=3**B, log(2)>1/2, pi**2<10 and a split at 2**n.
    The returned numbers are upper bounds, not evaluations of the kernel.
    """
    b = exact_integer(beta_floor_size, "B")
    n = exact_integer(n, "n", 1)
    denominator = sp.Rational(n, 2)-b
    if denominator <= 0:
        raise ValueError("The logarithmic tail denominator must be positive")
    length = sp.Rational(1, 2**(2*n))
    pole_argument = 3**b*length
    if pole_argument >= 1:
        raise ValueError("The exponential rational majorant requires pole_argument<1")
    pole = 80*pole_argument/(1-pole_argument)
    low_cut = sp.Rational(8, 2**n)
    high_cut = 80/denominator
    return {"B": sp.Integer(b), "n": sp.Integer(n), "length": length,
            "split": sp.Integer(2**n), "pole_argument_upper": pole_argument,
            "pole_norm_upper": pole, "low_cut_norm_upper": low_cut,
            "high_cut_norm_upper": high_cut,
            "norm_upper": pole+low_cut+high_cut}


def contraction_conditions(norm_upper, lipschitz, residual, radius):
    """Check Banach/a-posteriori inequalities given independently proved inputs.

    This does not check that a proposed nonlinear functional has those bounds
    or even belongs to the stated Banach space. That is a separate theorem.
    """
    norm_upper = exact_nonnegative(norm_upper, "operator norm")
    lipschitz = exact_nonnegative(lipschitz, "Lipschitz bound")
    residual = exact_nonnegative(residual, "fixed-point residual")
    radius = exact_nonnegative(radius, "ball radius")
    if radius == 0:
        raise ValueError("The closed ball must have positive radius")
    q = norm_upper*lipschitz
    if q >= 1:
        raise ValueError("A strict contraction factor q<1 is required")
    margin = radius-residual-q*radius
    if margin < 0:
        raise ValueError("The fixed-point map is not certified to preserve the ball")
    return {"q_upper": q, "self_map_margin": margin,
            "distance_upper": residual/(1-q), "radius": radius}


def identities():
    t, r, s, p = sp.symbols("t r s p", positive=True)
    beta, u = sp.symbols("beta u", real=True)
    den = (sp.log(r)+beta)**2+sp.pi**2
    result = {
        "positive_cut_density_jump": sp.simplify(
            1/(sp.log(r)+beta-sp.I*sp.pi)
            -1/(sp.log(r)+beta+sp.I*sp.pi)-2*sp.I*sp.pi/den),
        "simple_positive_pole_residue": sp.limit((s-p)/sp.log(s/p), s, p)-p,
        "integrated_pole": sp.integrate(p*sp.exp(p*u), (u, 0, t))-(sp.exp(p*t)-1),
        "integrated_cut": sp.integrate(sp.exp(-r*u), (u, 0, t))-(1-sp.exp(-r*t))/r,
        "log_tail_antiderivative": sp.diff(-1/(sp.log(r)+beta), r)
            -1/(r*(sp.log(r)+beta)**2),
        "classical_pi_upper_identity": sp.integrate(u**4*(1-u)**4/(1+u**2), (u, 0, 1))
            -(sp.Rational(22, 7)-sp.pi),
        "pi_squared_upper_margin": 10-sp.Rational(22, 7)**2-sp.Rational(6, 49),
    }
    # Direct convolution on polynomial data; Laplace transforms use the
    # exact harmonic-number derivative of Gamma, retaining Euler's constant.
    for degree in range(1, 5):
        convolution = sp.integrate(degree*u**(degree-1)*sp.log(t-u), (u, 0, t))
        result[f"log_convolution_degree_{degree}"] = sp.simplify(
            convolution-t**degree*(sp.log(t)-sp.harmonic(degree)))
        log_moment = sp.harmonic(degree)-sp.EulerGamma-sp.log(s)
        result[f"Laplace_multiplier_degree_{degree}"] = sp.simplify(
            beta-sp.EulerGamma+sp.harmonic(degree)-log_moment-(beta+sp.log(s)))
    return result


def controls():
    s, p = sp.symbols("s p", positive=True)
    return {
        "omitted_pole_Laplace_defect": p/(s-p),
        "wrong_signed_cut_jump": 4*sp.I*sp.pi/(1+sp.pi**2),
        "forward_t_data_derivative_diverges": "(beta-EulerGamma-log(t))/(8*pi^2) as t->0+",
        "inverse_constant_data_not_C1": "the positive cut gives J[1](t)/t unbounded as t->0+",
        "positive_pole_is_not_uniform_in_duration": "norm >= 8*pi^2*(exp(exp(-beta)*L)-1)",
    }


def calibration():
    bound = dyadic_norm_bound(1, 512)
    if bound["norm_upper"] >= sp.Rational(1, 3):
        raise ValueError("The isolated dyadic inverse norm exceeds its rounded cap")
    abstract = contraction_conditions(sp.Rational(1, 3), 1, sp.Rational(1, 10), 1)
    return {"isolated_inverse": bound, "rounded_norm_cap": sp.Rational(1, 3),
            "abstract_conditional_example": abstract,
            "example_is_a_SEE_solution": False}
