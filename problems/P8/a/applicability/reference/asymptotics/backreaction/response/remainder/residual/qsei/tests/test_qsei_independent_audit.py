"""Independent phase, Fourier, proper-clock and Fraction-only QSEI audits.

Only this file is owned by the A-track reviewer. The polynomial/Gaussian
controls are algebraic tests, not claimed physical prepared quantum states.
"""

import json
from fractions import Fraction as F
from math import comb

import sympy as sp
from p8a_qsei import mode_bounds, sampling
from p8a_remainder import verify as mode_prior
from p8a_residual import bounds as residual_bounds


def _fraction_constants():
    report = json.loads(mode_prior.REPORT.read_text())
    data = report["derived_finite_geometry"]
    b = list(map(F, data["U_derivative_bounds_through_5_per_delta"][:4]))
    cap, span = F(data["delta_bar"]), F(3)
    naive = {1: 2*span*b[0]}
    naive[2] = span*(2*b[1]+cap*b[0]*naive[1])
    naive[3] = span*(2*b[2]+cap*(2*b[1]*naive[1]+b[0]*naive[2]))
    products = [2*b[j]+cap*sum((comb(j, r)*b[r]*naive[j-r]
                               for r in range(j)), F(0)) for j in range(4)]
    inverse = [2*span*b[0]]+[(products[j-1]+span*products[j])/2 for j in range(1, 4)]
    h0, h1 = F(2), 4+cap*b[0]
    h2 = cap*b[1]+2*h0*h1
    error = [(1+h0)*inverse[0]+inverse[1],
             (1+h0)*inverse[1]+inverse[2]+h1*inverse[0],
             (1+h0)*inverse[2]+inverse[3]+2*h1*inverse[1]+h2*inverse[0]]
    f0 = F(25, 72)
    f1 = F(5, 6)+F(12, 31)*F(25, 36)
    f2 = 1+F(16, 31)*F(5, 6)+F(246, 961)*F(25, 36)
    weighted = F(8, 31)*F(5, 6)+F(164, 961)*F(25, 36)
    auxiliary = f2+F(3, 2)*weighted
    infrared = 33*b[0]
    error_root = 2*(error[0]*f2+2*error[1]*f1+error[2]*f0)+2*infrared*f0
    return {"b": b, "cap": cap, "naive": [naive[j] for j in range(1, 4)],
            "products": products, "inverse": inverse, "h": [h0, h1, h2],
            "error": error, "infrared": infrared, "f": [f0, f1, f2],
            "weighted": weighted, "auxiliary": auxiliary, "error_root": error_root}


def test_qsei_audit_exact_flat_past_transfer_and_nonflat_endpoint_control():
    x, past = sp.symbols("x past", real=True)
    k = sp.Symbol("k", positive=True)
    rate = 2*sp.I*k

    def integral_polynomial(source):
        if source == 0:
            return sp.Integer(0)
        # Independently solve (d_x-2ik)I=source, I(past)=0 using the
        # terminating polynomial particular solution and its homogeneous term.
        return sum((-sp.diff(source, x, j)+sp.exp(rate*(x-past))*sp.diff(source, x, j).subs(x, past))
                   /rate**(j+1) for j in range(sp.degree(source, x)+1))

    source = (x-past)**4
    exact = integral_polynomial(source)
    assert sp.simplify(sp.diff(exact, x)-rate*exact-source) == 0
    assert sp.simplify(exact.subs(x, past)) == 0
    for order in range(4):
        assert sp.simplify(sp.diff(exact, x, order)-integral_polynomial(sp.diff(source, x, order))) == 0
    # Omitting the past endpoint is false for a nonflat source.
    nonflat = integral_polynomial(sp.Integer(1))
    assert sp.simplify(sp.diff(nonflat, x)-sp.exp(rate*(x-past))) == 0
    assert sp.simplify(sp.diff(nonflat, x)) != 0


def test_qsei_audit_phase_volterra_kernel_and_first_ibp_normalization():
    duration, k = sp.symbols("duration k", positive=True)
    kernel = (sp.exp(2*sp.I*k*duration)-1)/(2*sp.I*k)
    assert kernel.subs(duration, 0) == 0
    assert sp.simplify(sp.diff(kernel, duration)-sp.exp(2*sp.I*k*duration)) == 0
    assert sp.simplify(sp.diff(kernel, duration, 2)-2*sp.I*k*sp.diff(kernel, duration)) == 0
    # |exp(i theta)-1|<=2 gives 1/k, not 1/(2k), for this kernel.
    theta = sp.Symbol("theta", real=True)
    modulus_squared = sp.expand_complex((sp.exp(sp.I*theta)-1)*(sp.exp(-sp.I*theta)-1))
    assert sp.trigsimp(modulus_squared-4*sp.sin(theta/2)**2) == 0
    assert modulus_squared.subs(theta, sp.pi) == 4
    independent = _fraction_constants()
    actual = mode_bounds.calibration()
    assert list(map(F, map(str, actual["b_inverse_frequency_constants"]))) == independent["inverse"]
    assert independent["inverse"][0] == 6*independent["b"][0]


def test_qsei_audit_all_mode_derivative_constants_with_fraction_arithmetic():
    separate = _fraction_constants()
    actual = mode_bounds.calibration()
    mapping = {"naive_b_derivatives_per_delta": "naive", "product_Ub_derivatives_per_delta": "products",
               "conformal_Hc_derivative_caps": "h", "ultraviolet_error_derivatives_per_delta": "error"}
    for exported, independent in mapping.items():
        assert list(map(F, map(str, actual[exported]))) == separate[independent]
    assert F(str(actual["infrared_error_per_delta"])) == separate["infrared"]
    assert 9*separate["cap"]*separate["b"][0] <= 1


def test_qsei_audit_error_product_derivatives_need_only_two_sampler_derivatives():
    x = sp.Symbol("x", real=True)
    k = sp.Symbol("k", positive=True)
    phase, hc = sp.Function("b")(x), sp.Function("hc")(x)
    sampler = sp.Function("F")(x)
    error = -sp.I*k*(phase-1)+sp.diff(phase, x)-hc*(phase-1)
    expected_second = (-sp.I*k*sp.diff(phase, x, 2)+sp.diff(phase, x, 3)
                       -sp.diff(hc, x, 2)*(phase-1)-2*sp.diff(hc, x)*sp.diff(phase, x)
                       -hc*sp.diff(phase, x, 2))
    assert sp.expand(sp.diff(error, x, 2)-expected_second) == 0
    assert sp.expand(sp.diff(sampler*error, x, 2)-sp.diff(sampler, x, 2)*error
                     -2*sp.diff(sampler, x)*sp.diff(error, x)-sampler*sp.diff(error, x, 2)) == 0


def test_qsei_audit_infrared_full_complex_parseval_not_half_parseval():
    # A complex modulated Gaussian is an explicit Fourier-only exclusion
    # control: its positive-frequency mass can exceed half its total mass.
    u = sp.Symbol("u", real=True)
    positive_mass = 2*sp.pi*sp.integrate(sp.exp(-(u-2)**2), (u, 1, sp.oo))
    full_mass = 2*sp.pi*sp.sqrt(sp.pi)
    assert sp.simplify((positive_mass-full_mass/2).rewrite(sp.erf)) == sp.pi**sp.Rational(3, 2)*sp.erf(1)
    assert sp.erf(1).is_positive is True
    # Full Parseval is 2pi and integral_0^1 k dk is 1/2. Relative to
    # hbar/(16pi^2), the resulting squared IR norm coefficient is four.
    coefficient = (sp.Rational(1, 4)/sp.pi**3)*(2*sp.pi)*sp.Rational(1, 2)/(sp.Rational(1, 16)/sp.pi**2)
    assert sp.simplify(coefficient) == 4


def test_qsei_audit_two_ibp_uv_moment_and_one_ibp_exclusion():
    alpha, k = sp.symbols("alpha k", positive=True)
    alpha_integral = sp.integrate((alpha+k)**-4, (alpha, 0, sp.oo))
    assert alpha_integral == 1/(3*k**3)
    assert sp.integrate(k*alpha_integral, (k, 1, sp.oo)) == sp.Rational(1, 3)
    first_only = sp.integrate((alpha+k)**-2, (alpha, 0, sp.oo))
    assert sp.integrate(k*first_only, (k, 1, sp.oo)) == sp.oo
    # With support length<=3, the normalized UV square coefficient is
    # 4*3/(3*pi)=4/pi; its square root is safely enlarged to two.
    assert sp.simplify((sp.Rational(1, 4)/sp.pi**3)*3*sp.Rational(1, 3)
                       /(sp.Rational(1, 16)/sp.pi**2)) == 4/sp.pi
    assert 4/sp.pi < 4


def test_qsei_audit_flat_spectral_cross_sign_and_eight_thirds_coefficient():
    k = sp.Symbol("k", positive=True)
    fr, fi, br, bi = sp.symbols("fr fi br bi", real=True)
    fourier, weighted = fr+sp.I*fi, br+sp.I*bi
    amplitude = -sp.I*k*fourier-weighted
    expected = k**2*(fr**2+fi**2)+br**2+bi**2-2*k*sp.im(fourier*sp.conjugate(weighted))
    assert sp.expand_complex(amplitude*sp.conjugate(amplitude)-expected).expand() == 0
    # The k integral contributes -2/3. For real F,hF, half-frequency
    # Parseval contributes -pi times integral F''(hF)', giving +8/3.
    cross = (sp.Rational(1, 4)/sp.pi**3)*(-sp.Rational(2, 3))*(-sp.pi)
    assert sp.simplify(cross/(sp.Rational(1, 16)/sp.pi**2)) == sp.Rational(8, 3)


def test_qsei_audit_exact_proper_clock_norms_and_plateau_caps():
    s = sp.Symbol("s", positive=True)
    a, psi = sp.Function("a", positive=True)(s), sp.Function("psi")(s)
    dx = lambda value: a*sp.diff(value, s)
    f = a**-sp.Rational(3, 2)*psi
    hubble = sp.diff(a, s)/a
    hc = sp.diff(a, s)
    assert sp.simplify(dx(dx(f))/sp.sqrt(a)-sp.diff(psi, s, 2)+2*hubble*sp.diff(psi, s)
                       -(3*hubble**2/4-3*sp.diff(hubble, s)/2)*psi) == 0
    assert sp.simplify(dx(hc*f)/sp.sqrt(a)-hubble*sp.diff(psi, s)
                       -(sp.diff(hubble, s)-hubble**2/2)*psi) == 0
    delta = sp.Symbol("delta", nonnegative=True)
    exact_hubble = 2*s/(4*s*s-delta)
    assert exact_hubble.subs({s: 2, delta: sp.Rational(1, 2)}) == sp.Rational(8, 31)
    assert (-sp.diff(exact_hubble, s)).subs({s: 2, delta: sp.Rational(1, 2)}) == sp.Rational(132, 961)
    assert sp.simplify(sp.diff(-sp.diff(exact_hubble, s), s)
                       +16*s*(4*s*s+3*delta)/(4*s*s-delta)**3) == 0
    separate = _fraction_constants()
    data = sampling.calibration()
    assert list(map(F, map(str, data["sampler_norms_per_proper_second_derivative"]))) == separate["f"]
    assert F(str(data["weighted_first_sampler_norm"])) == separate["weighted"]


def test_qsei_audit_physical_dimensions_and_exact_final_coefficient():
    normalization, eta_star = sp.symbols("A eta_star", positive=True)
    # d_eta=eta_star dx, k_physical=k/eta_star; the amplitude scales as
    # (A eta_star)^(-3/2), and k dk d_alpha supplies eta_star^(-3).
    spectral_scaling = (normalization*eta_star)**-3*eta_star**-3
    proper_second_derivative_norm_scaling = (normalization*eta_star**2)**-3
    assert sp.simplify(spectral_scaling-proper_second_derivative_norm_scaling) == 0
    separate = _fraction_constants()
    actual = sampling.calibration()
    assert separate["auxiliary"] == F(2026, 961)
    assert F(str(actual["auxiliary_flat_kernel_root_coefficient"])) == separate["auxiliary"]
    assert F(str(actual["mode_error_root_coefficient_per_delta"])) == separate["error_root"]
    upper = (separate["auxiliary"]+F(1, 10**14)*separate["error_root"])**2
    assert F(str(actual["maximum_QSEI_coefficient"])) == upper < 5
    assert all(value < sp.Rational(1, 100)
               for value in residual_bounds.reference_accuracy(sp.Rational(1, 10**14)).values())
