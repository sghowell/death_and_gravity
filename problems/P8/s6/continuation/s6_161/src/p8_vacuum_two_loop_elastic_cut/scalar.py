"""Full scalar parameter polynomials and physical-cut boundary estimates."""

from functools import cache

import sympy as s


@cache
def data():
    H, x, y, M, v, z, delta = s.symbols(
        "H x y M channel complementary_channel delta", real=True
    )
    L = 1 - H
    a, b, c, d = L * x, L * (1 - x), H * y, H * (1 - y)
    triangle = delta * L**2 + M * H
    box = triangle - z * H**2 * y * (1 - y)
    raw = L + M * H - v * a * b - L * H - z * c * d
    A = s.Symbol("quadratic_coefficient", real=True)
    general = delta + (M - 2 * delta) * H + A * H * H
    derivative = s.diff(general, H)
    f = H * (1 - H)
    w = f / derivative
    wp = s.diff(w, H) / derivative
    wpp = s.diff(wp, H) / derivative
    expected = s.diff(f, H) / derivative**2 - f * s.diff(derivative, H) / derivative**3
    expected2 = (
        s.diff(f, H, 2) / derivative**3
        - 3 * s.diff(f, H) * s.diff(derivative, H) / derivative**4
        + 3 * f * s.diff(derivative, H) ** 2 / derivative**5
    )
    r = s.Symbol("r", positive=True)
    tri_r = s.factor(triangle.subs(H, r / (1 + r)) * (1 + r) ** 2)
    box_r = s.factor(box.subs(H, r / (1 + r)) * (1 + r) ** 2)
    beta = s.Symbol("beta", positive=True)
    root = (1 + beta) / 2
    checks = {
        "literal_four_denominator_parameter_polynomial": s.expand(
            raw - box.subs(delta, 1 - v * x * (1 - x))
        ),
        "triangle_is_heavy_parameter_face": s.expand(box.subs(z, 0) - triangle),
        "general_quadratic_matches_box": s.expand(
            general.subs(A, delta - z * y * (1 - y)) - box
        ),
        "first_box_weight_derivative": s.factor(wp - expected),
        "second_box_weight_derivative": s.factor(wpp - expected2),
        "triangle_radial_polynomial": s.expand(tri_r - delta - M * r - M * r * r),
        "box_radial_polynomial": s.expand(
            box_r - delta - M * r - (M - z * y * (1 - y)) * r * r
        ),
        "box_measure_cancels_radial_denominator_scale": s.factor(
            H * (1 - H) * s.diff(r / (1 + r), r)
        ).subs(H, r / (1 + r))
        * (1 + r) ** 4
        - r,
        "box_weight_zero_at_first_endpoint": w.subs(H, 0),
        "box_weight_zero_at_second_endpoint": w.subs(H, 1),
        "light_cut_parameter_root": s.factor(
            (1 - v * x * (1 - x)).subs({x: root, v: 4 / (1 - beta * beta)})
        ),
        "light_cut_root_width": root - (1 - beta) / 2 - beta,
        "threshold_log_integral": s.integrate(-2 * s.log(r), (r, 0, 1)) - 2,
        "triangle_derivative_bound_at_M32": s.Rational(16, 9)
        + s.Rational(384, 27 * 32)
        - s.Rational(20, 9),
        "box_first_derivative_bound_at_M32": s.Rational(16, 3)
        + s.Rational(512, 27 * 32)
        - s.Rational(160, 27),
        "box_second_derivative_bound_at_M32": s.Rational(128, 27)
        + s.Rational(9216, 81 * 32)
        + s.Rational(98304, 243 * 32**2)
        - s.Rational(704, 81),
    }
    return {
        "triangle_parameter_polynomial": triangle,
        "box_parameter_polynomial": box,
        "box_boundary_weight": w,
        "box_first_boundary_weight_derivative": expected,
        "box_second_boundary_weight_derivative": expected2,
        "uniform_parameter_domain": "M>=32; channel and complementary channel in [-2,6]; delta in [-1/2,3/2]; A=delta-z y(1-y) in [-2,2]. On the extended H interval [-1,1], Delta' >=3M/4 and the unique Delta=0 root lies inside.",
        "light_log_integral_bound": "Integral_0^1 |log|1-v x(1-x)|| dx <=2 for v in [-2,6], including the continuous threshold boundary. The imaginary logarithm is bounded by pi<4.",
        "boundary_value_enclosures": {
            "bubble_MS": "ell+6",
            "triangle": "(3 log(M)+22)/M",
            "box": "(6 log(M)+52)/M^2",
        },
        "scalar_graph_ownership": "Three channels, each with one C^2 light bubble, four Cg full heavy-light triangles and four g^2 full two-heavy boxes, multiplied by 1/(2Q). This is the direct MS interaction reference; no fixed old sigma or additional field factor is hidden in it.",
        "checks": checks,
        "scope": "Bounds are taken after the one-dimensional Feynman boundary value and its endpoint subtraction. An integral of |Delta|^-2 through a real zero is not used. Full heavy propagators are retained.",
    }
