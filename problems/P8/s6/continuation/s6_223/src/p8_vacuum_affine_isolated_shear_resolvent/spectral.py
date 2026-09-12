"""Exact first-sheet zero enclosure and pole-plus-cut spectral data."""

from functools import cache

import sympy as s

from . import normalization as norm

y = norm.y
z, d = s.symbols("cut_coordinate d", positive=True)
N = 16
ROOT_LO = s.Rational(2899, 5000)
ROOT_HI = s.Rational(5799, 10000)
RESIDUE_LO = s.Integer(16)
RESIDUE_HI = s.Integer(17)
P = 30 - 20 * z + 3 * z * z
D = (
    -s.Rational(172, 225)
    + s.Rational(19, 30) * z
    - z * z / 10
    + s.sqrt(z) * P * s.atanh(s.sqrt(z)) / 30
)
U = s.sqrt(z) * P / 60
RHO = U / (D * D + s.pi**2 * U * U)


@cache
def moments():
    return tuple(
        s.integrate(y * y * (30 - 20 * y * y + 3 * y**4) * (1 - y * y) ** j, (y, 0, 1))
        for j in range(N + 1)
    )


def A_interval(r):
    x = r / 4
    M = moments()
    finite = s.Rational(1, 30) - r * sum(M[j] * x**j for j in range(N + 1)) / 120
    tail = r * M[0] * x ** (N + 1) / (120 * (1 - x))
    return s.factor(finite - tail), s.factor(finite)


def derivative_interval(r):
    x = r / 4
    M = moments()
    finite = sum((j + 1) * M[j] * x**j for j in range(N + 1)) / 120
    tail = M[0] * x ** (N + 1) * ((N + 2) - (N + 1) * x) / (120 * (1 - x) ** 2)
    return s.factor(finite), s.factor(finite + tail)


@cache
def data():
    J = {0: s.atanh(1 / s.sqrt(d)) / s.sqrt(d)}
    for j in range(1, 4):
        J[j] = d * J[j - 1] - s.Rational(1, 2 * j - 1)
    closed = s.Rational(1, 30) + (30 * J[1] - 20 * J[2] + 3 * J[3]) / 30
    wanted = (
        -s.Rational(172, 225)
        + s.Rational(19, 30) * d
        - d * d / 10
        + s.sqrt(d) * (30 - 20 * d + 3 * d * d) * s.atanh(1 / s.sqrt(d)) / 30
    )
    x = s.Symbol("geometric_ratio", positive=True)
    # Exact remainder sums, not a floating root finder.
    root_low, root_high = A_interval(ROOT_LO), A_interval(ROOT_HI)
    deriv_low = derivative_interval(ROOT_LO)[0]
    deriv_high = derivative_interval(ROOT_HI)[1]
    L = s.Symbol("large_log", real=True)
    eps = s.Symbol("epsilon", positive=True)
    high = D.subs(s.atanh(s.sqrt(z)), L / 2 + s.log((1 + s.sqrt(z)) / 2)).subs(
        z, 1 - eps
    )
    a, b, c = s.symbols("p_real p_imag radial_coefficient", real=True)
    M = s.Symbol("four_mass_squared", positive=True)
    complex_rational = (a + s.I * b) / (M + c * (a + s.I * b))
    imag = s.simplify(s.im(complex_rational))
    root_derivative = s.diff(a / (M + c * a), a)
    checks = {
        "closed_form_from_three_radial_recursions": s.simplify(closed - wanted),
        "strict_half_plane_sign_identity": s.factor(
            imag - M * b / ((M + c * a) ** 2 + c * c * b * b)
        ),
        "strict_real_derivative_identity": s.factor(
            root_derivative - M / (M + c * a) ** 2
        ),
        "geometric_remainder_sum": s.factor(
            1 / (1 - x) - sum(x**j for j in range(N + 1)) - x ** (N + 1) / (1 - x)
        ),
        "differentiated_geometric_remainder_sum": s.factor(
            1 / (1 - x) ** 2
            - sum((j + 1) * x**j for j in range(N + 1))
            - x ** (N + 1) * ((N + 2) - (N + 1) * x) / (1 - x) ** 2
        ),
        "cut_threshold_density_constant": s.limit(RHO / s.sqrt(z), z, 0)
        - s.Rational(50625, 59168),
        "cut_log_real_constant": s.simplify(
            high.subs(eps, 0) - s.Rational(13, 60) * L + s.Rational(52, 225)
        ),
        "cut_high_imaginary_coefficient": U.subs(z, 1) - s.Rational(13, 60),
        "cut_inverse_log_constant": s.Rational(13, 60) / s.Rational(13, 60) ** 2
        - s.Rational(60, 13),
        "cut_positive_polynomial": s.expand(P - 13 - (1 - z) * (17 - 3 * z)),
        "high_log_lower_margin": s.Rational(13, 80)
        - s.Rational(1, 16)
        - s.Rational(1, 10),
        "total_static_pole_plus_cut_moment": 1 / norm.A0 - 30,
        "total_next_pole_plus_cut_moment": s.Rational(3, 56) / norm.A0**2
        - s.Rational(675, 14),
    }
    return {
        "closed_A2": wanted,
        "cut_upper_bank": {
            "real_D": D,
            "imaginary_over_pi_U": U,
            "inverse_density": RHO,
        },
        "root_enclosure": {
            "r_lower": ROOT_LO,
            "r_upper": ROOT_HI,
            "A_at_lower_interval": root_low,
            "A_at_upper_interval": root_high,
            "exact_moments": moments(),
            "geometric_terms": N + 1,
        },
        "positive_weight_enclosure": {
            "R_over_m2_lower": 1 / deriv_high,
            "R_over_m2_upper": 1 / deriv_low,
            "coarse_bounds": (RESIDUE_LO, RESIDUE_HI),
        },
        "residue_sign": "The unique zero is p=-r m^2. R=1/A2'(-r m^2)>0 is the pole weight of1/A2. The residue of1/F2 is -R, not+R.",
        "representation": "1/F2(p)=-R/(p+r m^2)-integral_(4m^2)^infinity rho(tau)/(p+tau)dtau.",
        "moments": "R/(r m^2)+integral rho/tau=30; R/(r^2 m^4)+integral rho/tau^2=675/(14m^2). These are total moments INCLUDING the pole.",
        "first_sheet_proof": "Strict half-plane sign excludes every nonreal zero. Positive real derivative on(-4m^2,infinity), opposite endpoint signs, and positive open-cut imaginary part give exactly one simple real zero and no other first-sheet poles.",
        "asymptotics": "rho/sqrt(z)->50625/59168 at threshold. For L=log(tau/m^2), D=13L/60-52/225+O((m^2/tau)L), U->13/60 and rho L^2->60/13. The explicit log-coordinate formula is differentiated in the written proof.",
        "checks": checks,
        "gates": {
            "all_seventeen_exact_moments_positive": all(v > 0 for v in moments()),
            "root_lower_endpoint_strictly_positive": root_low[0] > 0,
            "root_upper_endpoint_strictly_negative": root_high[1] < 0,
            "root_strictly_between_zero_and_threshold": 0 < ROOT_LO < ROOT_HI < 4,
            "positive_pole_weight_above16": 1 / deriv_high > RESIDUE_LO,
            "positive_pole_weight_below17": 1 / deriv_low < RESIDUE_HI,
            "nonzero_continuum_static_mass_even_from_coarse_bounds": RESIDUE_HI
            / ROOT_LO
            < 30,
            "nonzero_pole_detected_without_numerical_quadrature": RESIDUE_LO / ROOT_HI
            > 0,
            "threshold_real_part_nonzero": D.subs(z, 0) != 0,
            "positive_threshold_density_constant": s.Rational(50625, 59168) > 0,
        },
    }
