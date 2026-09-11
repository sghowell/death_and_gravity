"""Strict rational enclosures for the complete heavy source."""

import sympy as s


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("An exact real rational is required")
    return s.Rational(value)


def scalar_source_enclosure(m, M, L, G, Q, U, B):
    m, M, L, G, Q, U, B = map(rational, (m, M, L, G, Q, U, B))
    if not (
        m >= 2 and M >= 2 and L >= 0 and G >= 0 and 0 < Q <= 144 and U >= 1 and B >= 1
    ):
        raise ValueError(
            "Need physical masses, nonnegative couplings and valid Q/power caps"
        )
    if not (m <= U**8 and 3 * M <= B**16):
        raise ValueError("The exact source power caps do not enclose the masses")
    g = G * G
    groups = {
        "differentiated_sunset": 7000 * G * g * U**2 * B**3 / Q**2,
        "mixed_bubble_times_light_double_propagator": 2048 * G * g * U**2 * B / Q**2,
        "correlated_OS_slope_times_tadpole": 70 * G * g * U**2 / (3 * Q**2),
        "proper_cubic_MS_counterterm": 280 * L * G * U / Q**2,
    }
    return {
        "power_caps": {"U": U, "B": B},
        "group_uppers": groups,
        "complete_scalar_and_cubic_source_absolute_upper": sum(groups.values()),
    }


def field_reexpression_upper(G, k0, k1, ell, Q):
    G, k0, k1, ell, Q = map(rational, (G, k0, k1, ell, Q))
    if not (min(G, k0, k1, ell) >= 0 and 0 < Q <= 144):
        raise ValueError(
            "Need nonnegative coefficient/log bounds and valid Q lower bound"
        )
    return G * ((ell + 1) * k0 + k1) / (2 * Q)


def data():
    G, g, Q, U, B = s.symbols("G g Q U B", positive=True)
    return {
        "source_scalar_bound_rule": "Half the G-weighted full scalar OS tadpole plus the proper cubic MS source insertion, all before finite-part extraction.",
        "field_map_rule": "The positive epsilon first field coefficient multiplies the pole of the entire first source. Bound its finite coefficient after full re-expression, not by truncating either factor first.",
        "checks": {
            "sunset_source_factor": G * g * s.Integer(28000) * U**2 * B**3 / (4 * Q**2)
            - 7000 * G * g * U**2 * B**3 / Q**2,
            "bubble_source_factor": G * g * (64 * U / Q) * (64 * U * B / Q) / 2
            - 2048 * G * g * U**2 * B / Q**2,
            "slope_source_factor": G * (s.Rational(2, 3) * g * U / Q) * (70 * U / Q) / 2
            - 70 * G * g * U**2 / (3 * Q**2),
            "cubic_source_factor": s.Symbol("L") * G * (16 / Q) * (70 * U / Q) / 4
            - 280 * s.Symbol("L") * G * U / Q**2,
        },
    }
