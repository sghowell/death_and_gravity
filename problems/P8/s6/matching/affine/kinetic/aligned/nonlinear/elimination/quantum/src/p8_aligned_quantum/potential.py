"""Dimensional-regularization pole and finite MSbar local vector potential."""
from functools import cache

import sympy as sp

from . import kernel

a, b = kernel.a, kernel.b
Lscale = sp.Symbol("log_m0_squared_over_mu_squared", real=True)
epsilon = sp.Symbol("epsilon_DR", positive=True)


@cache
def coefficients():
    longitudinal_weight = (a/b)**sp.Rational(3, 2)
    C = 2*b**2+a**sp.Rational(3, 2)*b**sp.Rational(1, 2)
    B = (2*b**2*(Lscale+sp.log(b)-sp.Rational(1, 2))
         +a**sp.Rational(3, 2)*b**sp.Rational(1, 2)*(Lscale+sp.log(a)-sp.Rational(3, 2)))
    # In d=4-2epsilon, there are d-2 transverse modes. The spatial
    # rescaling of the longitudinal determinant supplies c^(-3+2epsilon).
    weight = 2-2*epsilon+longitudinal_weight*sp.exp(epsilon*sp.log(b/a))
    scalar = -1/epsilon+Lscale+sp.log(b)-sp.Rational(3, 2)
    # The standard Euclidean momentum integral is
    # -m^d Gamma(-d/2)/(2(4pi)^(d/2)), with the MSbar scale included.
    # Divide by m^4/(64pi^2) and use Gamma(e-2)=Gamma(e)/((e-1)(e-2)).
    scalar_integral = (-2*sp.exp(epsilon*(sp.EulerGamma-Lscale-sp.log(b)))*sp.gamma(epsilon)
                       /((epsilon-1)*(epsilon-2)))
    scalar_expansion = sp.series(scalar_integral, epsilon, 0, 1).removeO()
    expansion = sp.series(b**2*weight*scalar, epsilon, 0, 1).removeO()
    residual = sp.factor(sp.expand_log(sp.expand(expansion+C/epsilon-B)))
    return {"pole_weight": C, "finite_weight": B,
            "longitudinal_weight": longitudinal_weight,
            "dimensional_weight": weight,
            "actual_scalar_integral_MSbar_finite_identity": sp.simplify(scalar_expansion-scalar),
            "dimensional_pole_and_finite_identity": residual,
            "scalar_Gamma_pole_residue": sp.limit(epsilon*sp.gamma(-2+epsilon), epsilon, 0)-sp.Rational(1, 2),
            "physical_clock_pole_three_modes": sp.factor(C.subs({a: 1, b: 1})-3),
            "physical_clock_MSbar_finite_vector_constant": sp.factor(B.subs({a: 1, b: 1})-3*Lscale+sp.Rational(5, 2)),
            "renormalization_log_coefficient": sp.factor(sp.diff(B, Lscale)-C)}


@cache
def clock_jets():
    p = kernel.geometry.P
    h = sp.Symbol("h", positive=True)
    p_N = -1/(2*h)
    p_NN = 3/(2*h)-1/(2*h**2)
    substitution = {a: kernel.masses()["a"], b: kernel.masses()["b"]}
    out = {"h": h}
    for name in ("pole_weight", "finite_weight"):
        value = coefficients()[name].subs(Lscale, 0).subs(substitution, simultaneous=True)
        p_jets = [sp.factor(sp.diff(value, p, j).subs(p, sp.Rational(1, 2))) for j in range(3)]
        out[name] = {"value": p_jets[0], "N_first": sp.factor(p_jets[1]*p_N),
                     "N_second": sp.factor(p_jets[2]*p_N**2+p_jets[1]*p_NN)}
    return out


@cache
def checks():
    data = coefficients()
    return {name: data[name] for name in ("dimensional_pole_and_finite_identity", "actual_scalar_integral_MSbar_finite_identity",
                                         "scalar_Gamma_pole_residue",
                                         "physical_clock_pole_three_modes", "physical_clock_MSbar_finite_vector_constant",
                                         "renormalization_log_coefficient")}


@cache
def independent_sign_controls():
    mass, radial = sp.symbols("mass radial_momentum_squared", positive=True)
    mass2 = mass**2
    primitive = sp.log(radial+mass2)+mass2/(radial+mass2)
    # Differentiating 1/2 Tr log(p²+m²) twice in m² gives
    # -1/(32pi²) times this positive radial integrand in four dimensions.
    radial_check = sp.factor(sp.diff(primitive, radial)-radial/(radial+mass2)**2)
    three_dimensional = -mass**3*sp.gamma(-sp.Rational(3, 2))/(2*(4*sp.pi)**sp.Rational(3, 2))
    return {"independent_cutoff_logarithm_primitive": radial_check,
            "independent_three_dimensional_scalar_integral": sp.simplify(three_dimensional+mass**3/(12*sp.pi))}
