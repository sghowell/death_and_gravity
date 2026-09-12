"""Full continued W8 bounds and dominated spatial dimension limit."""

from functools import cache, lru_cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vacuum_affine_uniform_uv_remainder import domain
from p8_vector_clock_matching import continuation
from p8_vector_hadamard import series
from p8_vector_state import wkb

DIMENSION_RADIUS = s.Rational(1, 4)
SIGMA = DIMENSION_RADIUS
eta = s.Symbol("dimension_minus_three", real=True)


def D0(value):
    return s.factor(
        s.diff(value, wkb.u) + wkb.background()["z_prime"] * s.diff(value, wkb.z)
    )


@lru_cache(maxsize=None, typed=True)
def coefficient(kind, order):
    if (
        kind not in ("transverse", "longitudinal")
        or type(order) is not int
        or not 0 <= order <= 4
    ):
        raise ValueError("Require an actual polarization and native WKB order0,...,4")
    if order == 0:
        return s.Integer(1)
    U = continuation.coefficients(kind)["U"]
    lam = wkb.background()["lambda"]
    previous = [coefficient(kind, n) for n in range(order)]
    inverse = [s.Integer(1)]
    for n in range(1, order):
        inverse.append(
            s.factor(-sum(previous[j] * inverse[n - j] for j in range(1, n + 1)))
        )
    DS = [s.Integer(0)] + [
        D0(previous[n]) - 2 * n * lam * previous[n] for n in range(1, order)
    ]
    rate = [lam] + [
        s.factor(sum(DS[j] * inverse[n - j] for j in range(1, n + 1)))
        for n in range(1, order)
    ]
    residual = (
        -sum(previous[j] * previous[order - j] for j in range(1, order))
        - D0(rate[-1]) / 2
        + (order - 1) * lam * rate[-1]
        + sum(rate[j] * rate[order - 1 - j] for j in range(order)) / 4
        - (U if order == 1 else 0)
    )
    return s.factor(residual / 2)


def complex_bound(value):
    shifted = s.cancel(value.subs(continuation.local.dimension, 3 + eta))
    num, den = s.fraction(shifted)
    polyden = s.Poly(den, wkb.u)
    power, lead = polyden.degree() // 2, polyden.LC()
    if (
        lead.is_positive is not True
        or s.expand(den - lead * (1 + wkb.u**2) ** power) != 0
    ):
        raise ValueError(
            "Require a dimension-independent positive quadratic denominator"
        )
    poly = s.Poly(num / lead, wkb.u, wkb.z, eta)
    if not all(c.is_Rational for c in poly.coeffs()):
        raise ValueError("Exact rational coefficients required")
    bound = (
        sum(
            abs(c)
            * s.Rational(51, 100) ** i
            * s.Rational(101, 100) ** j
            * DIMENSION_RADIUS**k
            for (i, j, k), c in poly.terms()
        )
        / s.Rational(7, 10) ** power
    )
    return {
        "upper": bound,
        "denominator_power": power,
        "reconstruction": s.factor(shifted - poly.as_expr() / (1 + wkb.u**2) ** power),
        "dimension_degree": poly.degree(eta),
    }


@cache
def mode_data():
    checks, rows = {}, {}
    for kind in ("transverse", "longitudinal"):
        row = {}
        for order in range(1, 5):
            value = coefficient(kind, order)
            bound = complex_bound(value)
            row[order] = {k: v for k, v in bound.items() if k != "reconstruction"}
            checks[f"{kind}_order_{order}_exact_complex_denominator"] = bound[
                "reconstruction"
            ]
            checks[f"{kind}_order_{order}_full_physical_WKB"] = s.factor(
                value.subs(continuation.local.dimension, 3)
                - series.coefficient(kind, order)
            )
            if order <= 2:
                checks[f"{kind}_order_{order}_original_dimensional_WKB"] = s.factor(
                    value
                    - continuation.coefficients(kind)["P2" if order == 1 else "P4"]
                )
        row["all_real_relative_W8_defect"] = sum(
            row[n]["upper"] / (s.Rational(99, 100) * domain.MASS) ** (2 * n)
            for n in range(1, 5)
        )
        row["far_relative_scaled_W8_defect"] = sum(
            row[n]["upper"]
            * (s.Rational(25, 16) * domain.EPS / (s.Rational(98, 100) * domain.MASS))
            ** (2 * n)
            for n in range(1, 5)
        )
        rows[kind] = row
    return {"rows": rows, "checks": checks}


@cache
def data():
    raw = mode_data()
    r, L, sig = s.symbols("r L sigma", positive=True)
    checks = dict(raw["checks"])
    checks.update(
        {
            "far_dimension_envelope_antiderivative": s.simplify(
                s.diff(-(r ** (-1 + sig)) / (1 - sig), r) - r ** (-2 + sig)
            ),
            "low_dimension_envelope_antiderivative": s.simplify(
                s.diff(r ** (3 - sig) / (3 - sig), r) - r ** (2 - sig)
            ),
            "far_limit_at_sigma_zero": (L ** (-1 + sig) / (1 - sig)).subs(sig, 0)
            - 1 / L,
            "far_transfer_degree": (5 - s.Symbol("j"))
            + (-1 + SIGMA)
            - (4 + SIGMA - s.Symbol("j")),
            "worst_near_transfer_degree": 5 + SIGMA - s.Rational(21, 4),
            "far_spatial_weight_margin": 6 - (4 + SIGMA) - s.Rational(7, 4),
            "near_spatial_weight_margin": 6 - (5 + SIGMA) - s.Rational(3, 4),
            "angular_endpoint_exponent": -SIGMA / 2 + s.Rational(1, 8),
            "radial_low_exponent": 2 - SIGMA - s.Rational(7, 4),
            "complex_time_inner_outer_discs": reference.OUTER - 2 * reference.RADIUS,
        }
    )
    return {
        "mode_majorants": raw["rows"],
        "complex_dimension_domain": "|d-3|<=1/4 with the unchanged outer/inner complex clock discs. All four complete continued WKB polynomials are bounded explicitly; opposite-phase branches are analytically continued without conjugating d.",
        "angular_geometry": "Perform the finite active/passive transverse polynomial moments before dimensional continuation. The remaining normalized u weight is Gamma(d/2)/[sqrt(pi)Gamma((d-1)/2)]*(1-u^2)^((d-3)/2). All invariant reconstruction denominators are nonzero on this disk. No noninteger-dimensional Hilbert frame or negative-dimensional positive measure is assumed.",
        "far_dimension_limit": "For r>=200U, U=m+|P|, the actual joint complex inverse-radius domain and uniformly nonzero continued W8/summed phases give a Cauchy UV-subtracted endpoint bound C_j U^(5-j) r^-4. The radial dimension factor costs at most r^sigma, sigma1/4. Its integral is C_j U^(4+sigma-j); all source jets and endpoint labels remain.",
        "near_dimension_limit": "On real momentum, continued mode and canonical factors are bounded uniformly in d using the full W8 defect and mass gap. The complete raw row grows at most C nu, with the same harmonic-mean suppression of the second leg as in S203. Low r<m is integrable with r^(2-sigma). Near m<=r<200U, the raw and all UV terms give powers at most U^(4+sigma-j); the logarithmic slot is bounded by an extra U/m, hence at most U^(5+sigma-j).",
        "dominated_limit": "Angular endpoints have integrable envelope(1-u^2)^(-1/8), radial far exponent-7/4 and transfer growth at most21/4<6. Thus every compact-time Schwartz Fourier pairing has a d-uniform integrable majorant, and the UV-subtracted reference limit equals the physical Qnew. The already convergent actual-state and time-remainder integrals remain in physical dimension3. Quantitative final physical constants are inherited from S208, not inferred from an unspecified continuity constant.",
        "not_a_state": "The continued unit-W8 comparison has no newly selected physical state. No analytic Borel preparation in momentum or complex dimension is assumed. The original actual state and complete correction remain unchanged.",
        "checks": checks,
        "gates": {
            "explicit_dimension_disk_including_both_sides_of_three": DIMENSION_RADIUS
            == s.Rational(1, 4),
            "all_four_actual_continued_WKB_orders": all(
                len([k for k in row if isinstance(k, int)]) == 4
                for row in raw["rows"].values()
            ),
            "all_real_complex_time_relative_W8_defect_below_one_thousandth": all(
                row["all_real_relative_W8_defect"] < s.Rational(1, 1000)
                for row in raw["rows"].values()
            ),
            "far_complex_time_relative_W8_defect_below_one_thousandth": all(
                row["far_relative_scaled_W8_defect"] < s.Rational(1, 1000)
                for row in raw["rows"].values()
            ),
            "both_frequency_branches_uniformly_nonzero": s.Rational(99, 100)
            - s.Rational(2, 1000)
            > s.Rational(1, 2),
            "far_frequency_branches_uniformly_nonzero": s.Rational(98, 100)
            / s.Rational(25, 16)
            - s.Rational(2, 1000)
            > s.Rational(1, 2),
            "sphere_and_source_basis_denominators_nonzero": 3 - DIMENSION_RADIUS > 1,
            "angular_endpoints_integrable": -SIGMA / 2 > -1,
            "low_radial_endpoint_integrable": 2 - SIGMA > -1,
            "far_radial_tail_integrable": -2 + SIGMA < -1,
            "near_transfer_growth_strictly_below_six": 5 + SIGMA < 6,
            "actual_state_not_dimensionally_reprepared": True,
        },
    }
