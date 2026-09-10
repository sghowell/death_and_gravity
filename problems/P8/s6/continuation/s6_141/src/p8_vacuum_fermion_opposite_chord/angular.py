"""Exact four-dimensional massless chord angular average."""

from functools import cache

import sympy as sp


@cache
def data():
    a, b, c = sp.symbols("a b c")
    A, B, t, x, y, z = sp.symbols("A B t x y z", positive=True)
    primitive = 2 * sp.atan(t * sp.sqrt(B / A)) / sp.sqrt(A * B)
    tangent_integral = sp.limit(primitive, t, sp.oo) - primitive.subs(t, 0)
    mean = 2 * (x + y - z) / (4 * x * y)
    checks = {
        "angular_polynomial_division": sp.factor(
            (1 - c * c) / (a - b * c)
            - (c / b + a / b**2 + (1 - a * a / b**2) / (a - b * c))
        ),
        "tangent_half_angle_denominator": sp.factor(
            (a - b * (1 - t * t) / (1 + t * t)) * (1 + t * t)
            - ((a - b) + (a + b) * t * t)
        ),
        "tangent_half_angle_integral": tangent_integral - sp.pi / sp.sqrt(A * B),
        "tangent_half_angle_antiderivative": sp.simplify(
            sp.diff(primitive, t) - 2 / (A + B * t * t)
        ),
        "ordered_radii_discriminant": sp.expand(
            (x + y) ** 2 - 4 * x * y - (x - y) ** 2
        ),
        "ordered_x_larger_shell_average": sp.factor(mean.subs(z, x - y) - 1 / x),
        "ordered_y_larger_shell_average": sp.factor(mean.subs(z, y - x) - 1 / y),
        "equal_radii_integrable_limit": sp.factor(mean.subs({z: 0, y: x}) - 1 / x),
        "normalized_S3_angle_measure": sp.integrate(sp.sin(t) ** 2, (t, 0, sp.pi))
        * 2
        / sp.pi
        - 1,
        "zero_small_radius_shell_limit": sp.limit(1 / x, y, 0) - 1 / x,
    }
    return {
        "normalized_angular_average": "(2/pi) integral_0^pi sin(theta)^2/[x+y-2sqrt(xy)cos(theta)] dtheta =1/max(x,y)",
        "x_y": "x=q^2,y=l^2; angular averages are over real Euclidean four-vectors",
        "equal_radius_limit": "At x=y>0 the theta=0 singularity cancels against sin(theta)^2; the average is 1/x.",
        "IR_domain": "The diagonal q=l is locally integrable in four relative dimensions. The radial positive majorant also proves integrability at q=l=0.",
        "positive_scalar_mass": "Adding b=1 only decreases the unshifted positive real chord; b=0 is the common upper bound.",
        "checks": checks,
    }
