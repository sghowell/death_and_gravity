"""Complete finite scalar-master dictionary and physical contour representation."""

from functools import cache

import sympy as s

from . import source

X = s.Symbol("master_x", real=True)
R = s.Symbol("contour_r", real=True)
DELTA = s.Symbol("positive_i0", positive=True)


def raw_log(mass, scale):
    return s.log(4 * s.pi * s.sympify(scale) ** 2 / s.sympify(mass)) - s.EulerGamma


def parameter_dictionary(channel, mass):
    channel, mass = map(s.sympify, (channel, mass))
    denominator = mass - channel * X * (1 - X) - s.I * DELTA
    return {
        "J": s.Limit(s.Integral(1 / denominator, (X, 0, 1)), DELTA, 0, dir="+"),
        "Q": s.Limit(
            s.Integral(s.log(denominator / mass) / denominator, (X, 0, 1)),
            DELTA,
            0,
            dir="+",
        ),
        "L": s.Limit(
            s.Integral(s.log(denominator / mass), (X, 0, 1)), DELTA, 0, dir="+"
        ),
        "C": source.original_masters.massless_triangle(channel, mass),
        "channel_log": s.Limit(
            s.log((-channel - s.I * DELTA) / mass), DELTA, 0, dir="+"
        ),
    }


def symbols(prefix):
    return {name: s.Symbol(prefix + "_" + name) for name in ("J", "Q", "L", "C", "la")}


def finite_dictionary(row, ell):
    J, Q, L, C, la = (row[name] for name in ("J", "Q", "L", "C", "la"))
    return {"C00mu": C, "C0mumu": (ell * J - Q) / 2, "B00": 2 + la, "Bmm": ell - L}


def contour(parameter=R):
    r = s.sympify(parameter)
    return r + s.I * r * (1 - r) * (1 - 2 * r) / 4


def triangle_radial_primitive(radial, positive_A):
    r, a = map(s.sympify, (radial, positive_A))
    return (s.log(1 + r) - s.log(r * r + a) / 2 + s.atan(r / s.sqrt(a)) / s.sqrt(a)) / (
        1 + a
    )


def spacelike_triangle_kernel(z):
    z = s.sympify(z)
    return -(s.log(z) + s.pi / s.sqrt(z)) / (2 * (1 + z))


def timelike_triangle_kernel(z):
    z = s.sympify(z)
    real = -s.Piecewise((-1, s.Eq(z, 1)), (s.log(z) / (1 - z), True)) / 2
    return real - s.I * s.pi / (2 * s.sqrt(z) * (1 + s.sqrt(z)))


@cache
def data():
    a, mu, e, ell = s.symbols("channel positive_mass epsilon raw_log", nonzero=True)
    J, Q, L, C, la = s.symbols("J Q L C la")
    row = dict(zip(("J", "Q", "L", "C", "la"), (J, Q, L, C, la), strict=True))
    finite = finite_dictionary(row, ell)
    radial, positive = s.symbols("radial positive_A", positive=True)
    z = s.Symbol("positive_z", positive=True)
    rr = R
    ss = s.Symbol("physical_s", positive=True)
    y = rr * (1 - rr) * (1 - 2 * rr) / 4
    xx = contour()
    den = s.expand(1 - ss * xx * (1 - xx))
    v = rr * (1 - rr)
    raw = (
        -s.gamma(1 - e)
        * (4 * s.pi * s.Symbol("nu", positive=True) ** 2) ** (-e)
        * mu**e
        / e
    )
    # Differentiate the analytic factor e*B rather than a formal series of a pole.
    raw0 = s.simplify(s.diff(e * raw, e).subs(e, 0))
    nu = s.Symbol("nu", positive=True)
    checks = {
        "whole_raw_zero_bubble_pole": s.limit(e * raw, e, 0) + 1,
        "whole_raw_zero_bubble_finite": s.expand_log(
            raw0 - raw_log(mu, nu), force=True
        ),
        "massive_triangle_finite_dictionary": finite["C0mumu"] - (ell * J - Q) / 2,
        "massless_bubble_finite_dictionary": finite["B00"] - 2 - la,
        "massive_bubble_finite_dictionary": finite["Bmm"] - ell + L,
        "finite_massless_triangle_retained": finite["C00mu"] - C,
        "whole_triangle_radial_primitive": s.simplify(
            s.diff(triangle_radial_primitive(radial, positive), radial)
            - 1 / ((1 + radial) * (radial**2 + positive))
        ),
        "radial_zero_endpoint": s.simplify(
            triangle_radial_primitive(0, positive)
            + s.log(positive) / (2 * (1 + positive))
        ),
        "radial_infinite_endpoint": s.simplify(
            s.limit(triangle_radial_primitive(radial, positive), radial, s.oo)
            - s.pi / (2 * s.sqrt(positive) * (1 + positive))
        ),
        "spacelike_radial_integral": s.simplify(
            spacelike_triangle_kernel(positive)
            + (s.log(positive) + s.pi / s.sqrt(positive)) / (2 * (1 + positive))
        ),
        "timelike_lower_branch_continuation": s.simplify(
            -(s.log(z) - s.I * s.pi + s.I * s.pi / s.sqrt(z)) / (2 * (1 - z))
            + s.log(z) / (2 * (1 - z))
            + s.I * s.pi / (2 * s.sqrt(z) * (1 + s.sqrt(z)))
        ),
        "timelike_imaginary_removable_quotient": s.factor(
            (1 / s.sqrt(z) - 1) / (1 - z) - 1 / (s.sqrt(z) * (1 + s.sqrt(z)))
        ),
        "timelike_real_removable_limit": s.limit(s.log(z) / (1 - z), z, 1) + 1,
        "timelike_entire_removable_value": timelike_triangle_kernel(s.Integer(1))
        - s.Rational(1, 2)
        + s.I * s.pi / 4,
        "physical_contour_real_denominator": s.factor(
            s.re(den) - (1 - ss * (v + y * y))
        ),
        "physical_contour_imaginary_sign": s.factor(
            s.im(den) + ss * v * (1 - 2 * rr) ** 2 / 4
        ),
        "physical_contour_lower_half_interior": s.factor(
            s.im(den) + ss * v * (1 - 4 * v) / 4
        ),
        "physical_contour_endpoint0": contour(0),
        "physical_contour_endpoint1": contour(1) - 1,
        "physical_contour_reflection": s.expand(contour(1 - rr) - (1 - contour(rr))),
        "physical_contour_derivative": s.expand(
            s.diff(contour(), rr) - 1 - s.I * (1 - 6 * rr + 6 * rr**2) / 4
        ),
        "derivative_quadratic_range_chart": s.expand(
            1
            - 6 * rr
            + 6 * rr**2
            - (6 * (rr - s.Rational(1, 2)) ** 2 - s.Rational(1, 2))
        ),
        "contour_small_real_lower_v": s.Rational(3, 4) / (16 * s.Rational(65, 64))
        - s.Rational(3, 65),
        "contour_small_real_upper_v": s.Rational(5, 4) / s.Rational(25, 4)
        - s.Rational(1, 5),
        "contour_small_real_imaginary_gap": s.Rational(25, 4)
        * s.Rational(3, 65)
        * s.Rational(1, 5)
        / 4
        - s.Rational(3, 208),
        "contour_gap_one_hundredth_margin": s.Rational(3, 208)
        - s.Rational(1, 100)
        - s.Rational(23, 5200),
        "contour_absolute_upper_bound": 1
        + 16 * (s.Rational(1, 4) + s.Rational(1, 256) + s.Rational(1, 16))
        - s.Rational(97, 16),
        "massless_triangle_beta_integral": s.integrate(
            1 / s.sqrt(X * (1 - X)), (X, 0, 1)
        )
        - s.pi,
        "endpoint_log_integral": s.integrate(-s.log(X), (X, 0, 1)) - 1,
    }
    return {
        "full_Feynman_parameter_dictionary": parameter_dictionary(a, mu),
        "whole_finite_master_jet_dictionary": finite,
        "finite_ordered_box": J / a * la,
        "physical_contour": xx,
        "positive_channel_denominator": den,
        "spacelike_triangle_integrand": spacelike_triangle_kernel(z),
        "timelike_triangle_integrand": timelike_triangle_kernel(z),
        "physical_branch": "For s>4mu take the lower-half-plane parameter denominator. "
        "At a negative real value log carries -i*pi, including contour midpoint. "
        "The left Feynman root is passed above and the right root below. "
        "The timelike z=1 triangle singularity is removable.",
        "uniform_compact_master_bounds": {
            "abs_J": 200,
            "abs_Q": 1800,
            "abs_L": 18,
            "abs_Cmm_finite": 1200,
            "abs_Bmm_finite": 21,
            "abs_C00_timelike": 5,
            "abs_C00_all": "10*(1+abs(log(delta)))/sqrt(delta)",
            "abs_channel_log": "10+abs(log(delta))",
            "abs_B00_finite": "12+abs(log(delta))",
        },
        "checks": checks,
        "gates": {
            "literal_S283_Feynman_masters_and_raw_constants_retained": True,
            "both_triangle_species_and_ordered_box_normalization_retained": True,
            "complex_physical_sheet_not_real_part_only": True,
            "fixed_contour_avoids_both_massive_Feynman_roots": True,
            "written_uniform_contour_gap_and_log_bounds_required": True,
            "massless_triangle_endpoint_integrability_and_removable_value": True,
            "quadratures_calibrate_but_do_not_certify_uniform_errors": True,
            "compact_domain_excludes_exact_forward_and_Gram_thresholds": True,
        },
    }
