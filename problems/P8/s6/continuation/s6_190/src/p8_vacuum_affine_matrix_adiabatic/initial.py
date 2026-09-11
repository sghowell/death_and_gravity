"""Quantitative matching to the unchanged all-order Cauchy state and UV tail."""

from functools import cache

import sympy as s
from p8_vector_hadamard import cutoffs, series
from p8_vector_state import wkb

from . import jets, riccati

PARTITION = s.Integer(10) ** 16
GRAPH_ERROR = s.Integer(10) ** 30
COVARIANCE_ERROR = 2 * s.Integer(10) ** 31
AMAX = s.Rational(25, 16)
RADIAL_SCALE = s.Rational(99, 100) / AMAX**2


@cache
def data():
    checks, rows = {}, {}
    H = wkb.background()["H"]
    for kind in ("transverse", "longitudinal"):
        d = H / 2 if kind == "transverse" else H * (s.Rational(1, 2) + wkb.z)
        checks[kind + "_actual_frame_second_order_frequency"] = s.cancel(
            series.D0(d) + d * d - wkb.frequency(kind)["U"]
        )
        coefficients = {}
        for n in range(1, 6):
            b = series.coefficient_bounds(kind, n)
            coefficients[n] = {
                "value": b["coefficient_upper"],
                "slope": b["frequency_coefficient_slope_upper"],
            }
            checks[f"{kind}_coefficient_reconstruction_{n}"] = b[
                "coefficient_reconstruction"
            ]
            checks[f"{kind}_slope_reconstruction_{n}"] = b["slope_reconstruction"]
        deviation = sum(
            row["value"] / jets.MASS ** (2 * n) for n, row in coefficients.items()
        )
        slope = sum(
            (row["slope"] + 2 * row["value"]) / jets.MASS ** (2 * n)
            for n, row in coefficients.items()
        )
        beta = (6 + slope / (2 * (1 - deviation))) / jets.MASS
        rows[kind] = {
            "actual_fifth_WKB_cutoff": cutoffs.threshold(kind, 5, jets.MASS),
            "exact_full_first_five_coefficient_bounds": coefficients,
            "complex_relative_frequency_deviation": deviation,
            "complex_normalized_graph_rate": beta,
            "complex_graph_modulus": (deviation + beta) / (2 - deviation - beta),
        }
    omega, wd, W, Wd, Wdd, d, dd = s.symbols(
        "omega omega_dot W W_dot W_ddot d d_dot", real=True, nonzero=True
    )
    v = Wd / (2 * W) + d
    graph = (omega - W + s.I * v) / (omega + W - s.I * v)
    derivative = sum(
        s.diff(graph, x) * dx for x, dx in ((omega, wd), (W, Wd), (Wd, Wdd), (d, dd))
    )
    squeeze = wd / (2 * omega) + d
    F = derivative - 2 * s.I * omega * graph - squeeze + graph * squeeze * graph
    residual = omega**2 - dd - d**2 - W**2 - Wdd / (2 * W) + 3 * Wd**2 / (4 * W**2)
    checks["exact_WKB_to_Riccati_residual_identity"] = s.cancel(
        F + s.I * residual * (1 + graph) ** 2 / (2 * omega)
    )
    lower = s.Symbol("lower_frequency", positive=True)
    cutoff = s.Symbol("analysis_partition", positive=True)
    checks["infinite_radial_tail_integral"] = s.integrate(
        lower ** (-7), (lower, cutoff, s.oo)
    ) - 1 / (6 * cutoff**6)
    c = riccati.constants()
    comparison = 10**6 * (
        jets.MASS**11 / PARTITION + c["residual_over_inverse_frequency_tenth"]
    )
    integral = 3 * COVARIANCE_ERROR * 4 / (12 * s.Integer(9)) * PARTITION ** (-6)
    return {
        "unchanged_initial_state": "The S6.55 locally finite all-order Cauchy frequency and slope remain. For initial nu above twice both fifth-order cutoffs, all first five WKB terms are fully on. The entire remaining cutoff-weighted series has both frequency and slope tail below nu^-10/32; no term or state is dropped.",
        "actual_Cauchy_matching": rows,
        "initial_graph_remainder": "On each scalar T/T/L initial block the W10 graph is analytic in the inverse-frequency disc |x|<=1/1000 with modulus<1/10. Formal WKB/Riccati uniqueness matches its first ten Taylor coefficients to the matrix reference. Thus ||r_actual(t0)-r_reference(t0)||<1e33 nu^-11 in the high band.",
        "frequency_envelope": "nu_minus=sqrt(m^2+(99/100)|k|^2/(25/16)^2), and nu_minus<=omega(t,k)<=3 nu_minus throughout every admitted shear history. This is not a polarization gap.",
        "analysis_partition": PARTITION,
        "actual_uniform_graph_error_coefficient": comparison,
        "graph_error_display": GRAPH_ERROR,
        "covariance_error_display": COVARIANCE_ERROR,
        "actual_high_band_result": "For nu_minus>=1e16, ||r_actual-r_reference||<1e30 nu_minus^-10 and the full real symmetrized balanced six-quadrature covariance differs by less than2e31 nu_minus^-10, uniformly on the actual CD slab.",
        "integrable_energy_weighted_covariance_tail_upper": integral,
        "tail_measure": "The integral over nu_minus>=1e16 of omega ||Sigma_actual-Sigma_reference||op d^3k/(2pi)^3 is below1e-65. The finite reference itself still needs the correct UV subtraction; this is only its actual-state integrable remainder.",
        "not_a_response_remainder": "No amplitude-parameter derivatives, covariant subtraction/contact matching, finite-coupling Taylor remainder, feedback inverse or physical cutoff follows from this covariance tail.",
        "checks": checks,
        "gates": {
            "actual_initial_canonical_rate_below_five": s.Rational(3, 2)
            * s.Rational(8, 5)
            < 5,
            "actual_initial_log_frequency_rate_below_two": s.Rational(8, 5) < 2,
            "actual_W10_real_slope_below_three_frequencies": all(
                2
                + sum(
                    v["slope"] / jets.MASS ** (2 * n)
                    for n, v in row["exact_full_first_five_coefficient_bounds"].items()
                )
                < 3
                for row in rows.values()
            ),
            "actual_high_band_frequency_stays_above_half": all(
                1
                - row["complex_relative_frequency_deviation"]
                - 1 / (32 * jets.MASS**11)
                > s.Rational(1, 2)
                for row in rows.values()
            ),
            "actual_half_log_rate_tail_below_one": s.Rational(7, 32) < 1,
            "same_all_order_cutoffs_fully_on_in_high_band": all(
                2 * row["actual_fifth_WKB_cutoff"] < PARTITION for row in rows.values()
            ),
            "complex_frequency_deviation_below_one_ten_thousandth": all(
                row["complex_relative_frequency_deviation"] < s.Rational(1, 10000)
                for row in rows.values()
            ),
            "complex_graph_rate_below_one_hundredth": all(
                row["complex_normalized_graph_rate"] < s.Rational(1, 100)
                for row in rows.values()
            ),
            "complex_W10_graph_below_one_tenth": all(
                row["complex_graph_modulus"] < s.Rational(1, 10)
                for row in rows.values()
            ),
            "higher_cutoff_geometric_tail": s.Rational(1, 64)
            / (1 - s.Rational(1, 2) / jets.MASS**2)
            < s.Rational(1, 32),
            "actual_initial_graph_vs_W10_tail": s.Rational(1, 32) + 1 / jets.MASS < 1,
            "Cauchy_and_actual_tail_display": jets.MASS**11 / 5 + 1 < jets.MASS**11,
            "uniform_actual_graph_error_display": comparison < GRAPH_ERROR,
            "actual_graph_and_reference_inside_covariance_ball": s.Rational(1, 100)
            + GRAPH_ERROR / PARTITION**10
            < s.Rational(1, 10),
            "full_covariance_error_display": 16 * GRAPH_ERROR < COVARIANCE_ERROR,
            "real_metric_exponential_lower_bound": 1 - jets.DELTA
            == s.Rational(99, 100),
            "uniform_frequency_upper_factor_three": 2 < 9 * RADIAL_SCALE,
            "radial_jacobian_constant_below_four": RADIAL_SCALE**3 > s.Rational(1, 16),
            "complete_infinite_energy_weighted_tail_display": integral
            < s.Rational(1, 10) ** 65,
        },
    }
