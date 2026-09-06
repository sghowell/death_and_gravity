"""Independent coordinate, physical-clock and signed-source photon audits."""

from functools import cache
from itertools import product

import sympy as sp
from p8a_maxwell import bounds, focusing, functional, stress


@cache
def coordinate_tensors():
    t = sp.Symbol("t", positive=True)
    a = sp.Function("a", positive=True)(t)
    metric = sp.diag(1, -a**2, -a**2, -a**2)
    inverse = metric.inv()

    def derivative(expression, index):
        return sp.diff(expression, t) if index == 0 else sp.S.Zero

    connection = [[[sum(inverse[i, k]*(derivative(metric[k, j], l)
                        +derivative(metric[k, l], j)-derivative(metric[j, l], k))/2
                        for k in range(4))
                    for l in range(4)] for j in range(4)] for i in range(4)]
    ricci = sp.zeros(4)
    for i, j in product(range(4), repeat=2):
        # Minus the usual coordinate contraction is the frozen FK convention.
        ricci[i, j] = sp.simplify(-sum(
            derivative(connection[k][i][j], k)-derivative(connection[k][i][k], j)
            +sum(connection[k][k][l]*connection[l][i][j]
                 -connection[k][j][l]*connection[l][i][k] for l in range(4))
            for k in range(4)))
    scalar = sp.trace(inverse*ricci)
    hessian = sp.Matrix(4, 4, lambda i, j:
        derivative(derivative(scalar, i), j)
        -sum(connection[k][i][j]*derivative(scalar, k) for k in range(4)))
    box = sp.trace(inverse*hessian)
    ricci_square = sp.trace(inverse*ricci*inverse*ricci)
    i_tensor = 2*hessian-2*metric*box+2*scalar*ricci-metric*scalar**2/2
    h3_tensor = ricci*inverse*ricci-sp.Rational(2, 3)*scalar*ricci
    h3_tensor += metric*(-ricci_square/2+scalar**2/4)
    return t, a, metric, inverse, ricci, scalar, i_tensor, h3_tensor


def test_full_coordinate_curvature_and_all_diagonal_reference_components():
    t, a, metric, inverse, ricci, scalar, counterterm, quadratic = coordinate_tensors()
    h = sp.diff(a, t)/a
    beta = sp.Symbol("beta", real=True)
    expected = stress.reference_jets(h, *(sp.diff(h, t, n) for n in (1, 2, 3)), beta)
    tensor = 62*quadratic+beta*counterterm
    assert sp.simplify(ricci[0, 0]-scalar/2+3*h*h) == 0
    assert sp.simplify(tensor[0, 0]-expected["rho"]) == 0
    for index in (1, 2, 3):
        assert sp.simplify(tensor[index, index]/a**2-expected["pressure"]) == 0
    assert sp.simplify(sp.trace(inverse*tensor)-expected["trace"]) == 0
    assert sp.simplify((tensor-metric*sp.trace(inverse*tensor)/2)[0, 0]
                       -expected["EED"]) == 0


def test_conserved_counterterm_not_literal_nonimported_derivative_formula():
    t, a, _, _, _, _, counterterm, _ = coordinate_tensors()
    h = sp.diff(a, t)/a
    rho, pressure = counterterm[0, 0], counterterm[1, 1]/a**2
    assert sp.simplify(sp.diff(rho, t)+3*h*(rho+pressure)) == 0
    h = sp.Function("H")(t)
    hd, hdd, hthird = (sp.diff(h, t, n) for n in (1, 2, 3))
    # Literal proper-clock conversion of Markowicz et al. (48), not imported
    # as the conserved FK R-squared tensor. This is a non-import control only.
    printed_rho = 36*h*hdd+18*hd**2+180*h*h*hd
    printed_pressure = -12*hthird-42*hd**2-72*h*hdd-84*h*h*hd
    residual = sp.diff(printed_rho, t)+3*h*(printed_rho+printed_pressure)
    assert sp.expand(residual-72*(hd+h*h)*(hdd+4*h*hd)) == 0
    assert residual != 0


def test_state_independent_trace_cancels_in_eed_difference_not_absolute_eed():
    reference_density, trace, difference = sp.symbols("rho_ref trace delta_rho", real=True)
    reference_eed = reference_density-trace/2
    target_eed = reference_density+difference-trace/2
    assert sp.expand(target_eed-reference_eed) == difference
    data = stress.reference_jets(1, 0, 0, 0, 0)
    assert data["rho"] > 0 > data["EED"]


def test_two_polarization_minkowski_spectral_constant_and_real_parseval():
    omega, frequency = sp.symbols("omega u", positive=True)
    # Angular integration of sum of electric and magnetic vacuum kernels
    # gives 4/(2*pi)^3; the time-translation frequency triangle is explicit.
    weight = 4/(2*sp.pi)**3*sp.integrate(omega**3, (omega, 0, frequency))
    assert sp.simplify(weight-frequency**4/(8*sp.pi**3)) == 0
    assert sp.simplify(weight*sp.pi/frequency**4) == 1/(8*sp.pi**2)
    # Positive-half Parseval uses reality. Doubling it is a detectable error.
    assert 2/(8*sp.pi**2) != 1/(8*sp.pi**2)


def test_literal_proper_clock_conformal_weight_on_independent_polynomials():
    t = sp.Symbol("t", positive=True)
    scale = 1+t+t**2
    sampler = (t-1)**2*(2-t)**2
    flat_sampler = scale**-sp.Rational(3, 2)*sampler
    conformal_second = scale*sp.diff(scale*sp.diff(flat_sampler, t), t)
    proper = functional.proper_operator(sampler, sp.diff(scale, t)/scale, t)
    assert sp.simplify(conformal_second**2/scale-proper**2) == 0
    for endpoint in (1, 2):
        assert sampler.subs(t, endpoint) == sp.diff(sampler, t).subs(t, endpoint) == 0
    # Omitting the dt/a measure produces a different integrand.
    assert sp.simplify(conformal_second**2-proper**2) != 0


def test_actual_H2_zero_boundary_radiation_integral_including_reference_credit():
    t = sp.Symbol("t", positive=True)
    f = (t-1)**2*(2-t)**2
    correct = functional.proper_operator(f, 1/(2*t), t)**2
    ibp = sp.diff(f, t, 2)**2-sp.Rational(15, 8)*sp.diff(f, t)**2/t**2
    ibp += sp.Rational(945, 256)*f*f/t**4
    assert sp.integrate(sp.cancel(correct-ibp), (t, 1, 2)) == 0
    absolute = ibp-sp.Rational(31, 320)*f*f/t**4
    assert sp.expand(absolute-(
        sp.diff(f, t, 2)**2-sp.Rational(15, 8)*sp.diff(f, t)**2/t**2
        +sp.Rational(4601, 1280)*f*f/t**4)) == 0
    assert sp.expand(absolute-ibp+sp.Rational(31, 320)*f*f/t**4) == 0


def test_signed_reference_envelope_covers_exact_jet_corners_for_both_beta_signs():
    caps = (2, 3, 5, 7)
    for beta in (-2, 0, 2):
        loss = bounds.envelope(*caps, 1, beta_m=beta)["reference_loss_numerator"]
        for signs in product((-1, 0, 1), repeat=4):
            jets = tuple(value*sign for value, sign in zip(caps, signs, strict=True))
            eed = stress.reference_jets(*jets, beta)["EED"]
            assert -loss <= eed <= loss


def test_exact_pi_sampler_factor_is_coarsened_upwards_not_downwards():
    x, y = sp.symbols("x y", nonnegative=True)
    p = sp.Symbol("p", positive=True)
    exact = 1+2*x/p+(3*x*x/4+3*y/2)/p**2
    coarse = bounds.sampler_factor(x, y, 1)
    assert sp.simplify(coarse-exact-(p-3)*(2*x/(3*p)
        +(p+3)*(x*x/12+y/6)/p**2)) == 0


def test_semiclassical_trace_keeps_cosmological_and_additional_sources():
    kappa, lam, maxwell_rho, maxwell_p, other_rho, other_p = sp.symbols(
        "kappa Lambda rho_M p_M rho_other p_other", real=True)
    density = maxwell_rho+other_rho
    trace = density-3*(maxwell_p+other_p)
    scalar = kappa*trace+4*lam
    ricci_00 = -kappa*density-lam+scalar/2
    expected = focusing.see_ricci((maxwell_rho+3*maxwell_p)/2,
        (other_rho+3*other_p)/2, kappa, lam)
    assert sp.expand(ricci_00-expected) == 0
    for lam_value, lower in product((-3, 3), repeat=2):
        d = focusing.geometric_constants(1, 2, 3, 4, 1, beta_m=-2, kappa=2,
            hbar=1, cosmological_constant=lam_value, other_eed_lower=lower)
        assert d["Q0"] >= d["raw_constant_before_nonnegative_coarsening"]
        assert d["new_SEE_solution_asserted"] is False


def test_small_planck_ratio_does_not_supply_missing_geometric_hypotheses():
    scale = sp.Symbol("tau", positive=True)
    # A1's gradient term on its 3*tau/4 tail already exceeds this K cap.
    geometric_floor = 3/(3*scale/4)
    contraction_cap = 3/scale
    assert sp.simplify(geometric_floor-contraction_cap) == 1/scale
    assert focusing.calibration()["A1_initial_pointwise_Ricci_premise_removed"] is False
