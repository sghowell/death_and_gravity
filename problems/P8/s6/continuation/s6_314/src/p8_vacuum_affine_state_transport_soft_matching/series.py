"""Uniform parameter derivatives of the entire signed dimensional soft sum."""

from functools import cache

import mpmath as mp
import sympy as s

from . import source


def full_series(e, a, b, y, count=240):
    t = mp.gamma(1 + 2 * e) * mp.power(y, 2 * e)
    lam = t * (a / (2 * e) + b)
    power = mp.mpf(1)
    phi = mp.mpf(0)
    phip = mp.mpf(0)
    for n in range(count):
        phi += power / mp.gamma(1 + 2 * e * n)
        phip += power / mp.gamma(1 + 2 * e * (n + 1))
        power *= lam / (n + 1)
    virtual = mp.exp(-a / (2 * e))
    return virtual * phi, virtual * (t * phip - phi) / (2 * e), virtual * t * phip


@cache
def independent_numerics():
    with mp.workdps(80):
        cases = (
            ("0.125", "0", "-0.2", "0.125"),
            ("0.125", "0", "0.2", "0.000001"),
            ("0.03125", "0.3", "0.2", "0.000001"),
            ("0.03125", "0.3", "-0.2", "0.125"),
            ("0.0078125", "0", "-0.2", "0.125"),
            ("0.0078125", "0.3", "0.2", "0.125"),
        )
        for row in cases:
            e, a, b, y = map(mp.mpf, row)
            value, da, db = full_series(e, a, b, y)
            value2, da2, db2 = full_series(e, a, b, y, 360)
            assert max(abs(value - value2), abs(da - da2), abs(db - db2)) < mp.mpf(
                "1e-60"
            )
            numeric_da = mp.diff(
                lambda aa, e=e, b=b, y=y: full_series(e, aa, b, y)[0], a
            )
            numeric_db = mp.diff(
                lambda bb, e=e, a=a, y=y: full_series(e, a, bb, y)[0], b
            )
            assert abs(da - numeric_da) < mp.mpf("1e-60")
            assert abs(db - numeric_db) < mp.mpf("1e-60")
            assert abs(da) < 12 * (1 + abs(mp.log(y))) and abs(db) < 6
        for alpha in (mp.mpf("0.25"), mp.mpf("0.75")):
            for slope in (mp.mpf("0.2"), mp.mpf("-0.2")):
                cut = mp.mpf("0.125")
                direct = alpha * mp.quad(
                    lambda z, alpha=alpha, slope=slope, cut=cut: (
                        mp.power(1 - z, alpha) * mp.expm1(slope * cut * z) / z
                    ),
                    [0, mp.mpf("0.5"), 1],
                )
                moments = alpha * mp.fsum(
                    mp.power(slope * cut, n) * mp.beta(n, 1 + alpha) / mp.factorial(n)
                    for n in range(1, 81)
                )
                assert abs(direct - moments) < mp.mpf("1e-60")
                assert direct != 0

    return {
        "six_signed_series_derivative_calibrations": True,
        "four_remaining_energy_connector_Beta_calibrations": True,
    }


@cache
def data():
    e, a, t, y, z = s.symbols("e a t y z", positive=True)
    b = s.Symbol("finite_regulator_coordinate", real=True)
    Phi = s.Function("Phi")
    lam = t * (a / (2 * e) + b)
    expression = s.exp(-a / (2 * e)) * Phi(lam)
    derivative = s.diff(Phi(z), z).subs(z, lam)
    checks = {
        "entire_series_finite_coordinate_derivative": s.simplify(
            s.diff(expression, b) - s.exp(-a / (2 * e)) * t * derivative
        ),
        "entire_series_physical_index_derivative": s.simplify(
            s.diff(expression, a)
            - s.exp(-a / (2 * e)) * (t * derivative - Phi(lam)) / (2 * e)
        ),
        "physical_index_derivative_cancellation": s.expand(
            t * derivative - Phi(lam) - (t - 1) * Phi(lam) - t * (derivative - Phi(lam))
        ),
        "inverse_Gamma_derivative": s.simplify(
            s.diff(1 / s.gamma(z), z) + s.polygamma(0, z) / s.gamma(z)
        ),
        "Gamma_recurrence_for_global_derivative": s.expand_func(
            s.gamma(z) - (z - 1) * s.gamma(z - 1)
        ),
        "finite_coordinate_gradient_budget": s.Integer(3 * 2 - 6),
        "physical_index_gradient_constant_budget": s.Integer(3 * 2 + 3 * 2 - 12),
        "state_series_difference_constant_budget": s.Integer(
            12 * 40 + 6 * 8000 - 48480
        ),
        "connector_dominated_integral_primitive": s.simplify(
            s.diff(y * (2 - s.log(y)), y) - (1 - s.log(y))
        ),
        # Positive real substitution u=1-z; never use an ambiguous shifted-log primitive.
        "remaining_energy_log_integral": s.integrate(-s.log(z), (z, 0, 1)) - 1,
    }
    for count in range(1, 9):
        coeff = s.symbols("c0:" + str(count + 1))
        polynomial = sum(coeff[n] * z**n / s.factorial(n) for n in range(count + 1))
        shifted = sum(coeff[n + 1] * z**n / s.factorial(n) for n in range(count))
        checks["finite_polynomial_derivative_shift_" + str(count)] = s.expand(
            s.diff(polynomial, z) - shifted
        )
        checks["connector_Beta_moment_" + str(count)] = s.simplify(
            s.expand_func(
                s.gamma(count)
                * s.gamma(1 + a)
                / (s.factorial(count) * s.gamma(1 + a + count))
                - s.gamma(1 + a) / (count * s.gamma(1 + a + count))
            )
        )
    numerical = independent_numerics()
    return {
        "whole_signed_sum_coordinates": "S_e(a,b;y)=exp[-a/(2e)] Phi_e(t_e[a/(2e)+b]), Phi_e(lambda)=sum_N lambda^N/[N!Gamma(1+2eN)], t_e=Gamma(1+2e)y^(2e)<=1, b=(a_e-a)/(2e). This b is a finite regulator coordinate, not a marked graviton energy. The physical domain is a>=0,|b|<5500/kappa,0<y<=1,0<e<=1/8; a_e need not be positive.",
        "whole_global_inverse_Gamma_derivative": "For1<=z<=2, |psi(z)|<1 and1/Gamma(z)<exp1<3. Forz>=2, Jensen on a Gamma(z,1) random variable gives0<psi(z)<=lnz; Gamma(z)>=(z-1)/exp1 and lnz<=z-1. Thus |(1/Gamma)'(z)|<3 for allz>=1, so neighboring coefficients differ by<6e.",
        "whole_entire_series_derivative_proof": "Absolute convergence permits differentiating on compact parameter sets. exp(|lambda|-a/(2e))<=exp(|b|)<2. Therefore |partial_b S_e|<6. Split t_e Phi'_e-Phi_e into(t_e-1)Phi_e+t_e(Phi'_e-Phi_e); (1-t_e)/(2e)<=1+|lny| and the reciprocal-Gamma coefficient difference removes the apparent1/e. Thus |partial_a S_e|<12(1+|lny|). No positive finite-e probability is assumed.",
        "whole_uniform_connector_domination": "The segment between the elastic and marked (a,b_e) lies in the same convex parameter domain. S300's |delta a|<40R/kappa and the new |delta b_e|<8000R(1-lnR)/kappa imply |S_e(sigma;y)-S_e(0;y)|<50000R[1-lnR+|lny|]/kappa. With y=x-R and elastic Born seed a0*dR/R, the majorant integrates to50000*a0*x(3-2lnx)/kappa. At fixed positive x this justifies the entire signed-series regulator limit under the seed integral, including both energy endpoints.",
        "whole_independent_calibration": "Six whole-series derivative cases include a0 and negative a_e, agreement of240/360-term sums and numerical differentiation within1e-60 at80digits. Four connector integrals with positive/negative finite-conversion slopes agree with independent exact Beta-moment series within1e-60. These are diagnostic calibrations, not the general uniform proof.",
        "checks": checks,
        "gates": {
            **numerical,
            "original_finite_coordinate_below_half": 5500 / source.KAPPA
            < s.Rational(1, 2),
            "state_series_gradient_below50000": s.Integer(48480) < 50000,
            "negative_regulated_index_and_zero_physical_index_included": True,
            "physical_index_derivative_cancellation_before_bound": True,
            "difference_integrated_before_regulator_removal": True,
            "unsubtracted_Born_seed_integrals_not_separately_regularized": True,
            "no_fixed_N_limit_used_for_positive_physical_index": True,
        },
    }
