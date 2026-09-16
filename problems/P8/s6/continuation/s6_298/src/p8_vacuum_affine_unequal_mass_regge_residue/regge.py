"""Controlled spin-order approximation and conditional finite-transfer budget."""

from functools import cache

import sympy as s

from . import source, waves

MU, N, G = source.MU, source.N, source.G
T = s.Symbol("timelike_transfer", positive=True)
CAP, SCALE, EPS, ARC = s.symbols(
    "transfer_cap trajectory_scale_squared residue_ratio_error arc_supremum",
    positive=True,
)
X = s.Symbol("Q_probability_u", positive=True)


def require_budget(cap, scale, error, arc, mass, heavy):
    cap, scale, mass, heavy = map(source.require_mass, (cap, scale, mass, heavy))
    for value in (error, arc):
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, s.Integer, s.Rational))
            or value < 0
        ):
            raise ValueError("Require exact nonnegative declared errors")
    if heavy <= 4 * mass or cap <= 4 * heavy:
        raise ValueError("Require cap above both thresholds and separated masses")
    return cap, scale, s.sympify(error), s.sympify(arc), mass, heavy


def order_loss(species, energy=T, mass=MU, heavy=N):
    a, b, active = waves.kinematics(species, energy, mass, heavy)
    return (
        s.acosh(a / b)
        + s.Rational(1, 4)
        + s.log((energy - 4 * mass) / (energy - 4 * active)) / 2
    )


def closed_budget(cap=CAP, scale=SCALE, error=EPS, arc=ARC, mass=MU, heavy=N, cubic=G):
    cap, scale, error, arc, mass, heavy, cubic = map(
        s.sympify, (cap, scale, error, arc, mass, heavy, cubic)
    )
    q = cap - 4 * mass
    x = cap - 4 * heavy
    A = 2 * cap + 4 * heavy
    trajectory = (
        cubic**2
        * (
            q * q * (s.log(A / q) + s.Rational(3, 4))
            + x * x * (s.log(A / x) + s.log(cap / x) / 2 + 1)
        )
        / (480 * s.pi**2 * heavy**3 * scale)
    )
    return {
        "unknown_transfer_arc": arc / cap,
        "known_full_graph_tail": cubic**2 / (16 * s.pi**2 * heavy * cap),
        "declared_residue_ratio_error": error * cubic**2 / (144 * s.pi**2 * heavy**2),
        "declared_positive_trajectory_error": trajectory,
        "total": arc / cap
        + cubic**2 / (16 * s.pi**2 * heavy * cap)
        + error * cubic**2 / (144 * s.pi**2 * heavy**2)
        + trajectory,
    }


@cache
def data():
    u = X
    d, order = s.symbols("reduced_threshold order", positive=True)
    B = s.Symbol("scaled_Q_denominator", positive=True)
    checks = {
        "whole_Q_probability_denominator_Jacobian": s.factor(
            ((B - 1) * (B - d)).subs(B, 1 / u) * u * u - (1 - u) * (1 - d * u)
        ),
        "whole_Q_probability_power_Jacobian": s.powsimp(
            (B ** (-order - 1)).subs(B, 1 / u) / u - u**order
        ),
        "Q_reference_normalization": s.integrate(u * u / s.sqrt(1 - u), (u, 0, 1))
        - s.Rational(16, 15),
        "Q_reference_mean_excess": s.integrate(u * s.sqrt(1 - u), (u, 0, 1))
        - s.Rational(4, 15),
        "Q_reference_mean_excess_ratio": s.Rational(4, 15) / s.Rational(16, 15)
        - s.Rational(1, 4),
        "exponential_defect_initial": (s.exp(-u) - 1 + u).subs(u, 0),
        "exponential_defect_derivative": s.diff(s.exp(-u) - 1 + u, u) - (1 - s.exp(-u)),
        "two_positive_kernel_product_defect": s.expand(
            1 - (1 - u) * (1 - d) - (u + d - u * d)
        ),
    }
    q, A, U = s.symbols("threshold_offset logarithm_anchor cap", positive=True)
    primitive = q * q * s.log(A / q) / 2 + q * q / 4
    checks["whole_logarithmic_threshold_primitive"] = s.factor(
        s.diff(primitive, q) - q * s.log(A / q)
    )
    light = q * q * (s.log(A / q) + s.Rational(3, 4)) / 2
    heavy = q * q * (s.log(A / q) + s.log(U / q) / 2 + 1) / 2
    checks["whole_light_trajectory_integral"] = s.factor(
        s.diff(light, q) - q * (s.log(A / q) + s.Rational(1, 4))
    )
    checks["whole_heavy_trajectory_integral"] = s.factor(
        s.diff(heavy, q) - q * (s.log(A / q) + s.log(U / q) / 2 + s.Rational(1, 4))
    )
    checks["light_log_endpoint_vanishes"] = s.limit(light, q, 0, dir="+")
    checks["heavy_log_endpoint_vanishes"] = s.limit(heavy, q, 0, dir="+")
    qi, qa = s.symbols("qi qa", positive=True)
    checks["light_density_bound_normalization"] = s.factor(
        s.Rational(1, 2) * 2 / s.pi / (240 * s.pi) - 1 / (240 * s.pi**2)
    )
    checks["heavy_density_kinematic_cancellation"] = s.factor(
        qa / qi * qa * qi - qa * qa
    )
    checks["two_external_slope_tail_normalization"] = 2 * s.Rational(
        1, 32
    ) - s.Rational(1, 16)
    checks["residue_error_total_slope_normalization"] = 2 * s.Rational(
        1, 288
    ) - s.Rational(1, 144)
    checks["trajectory_integrated_normalization"] = s.Rational(1, 240) / 2 - s.Rational(
        1, 480
    )
    pieces = closed_budget()
    checks["complete_conditional_error_assembly"] = s.expand(
        sum(v for k, v in pieces.items() if k != "total") - pieces["total"]
    )
    slope, ell = s.symbols("F1prime normalized_Regge_log_slope", real=True)
    checks["selected_forward_mismatch_sign"] = (
        (-ell / source.K) - (-2 * slope / source.K) + (ell - 2 * slope) / source.K
    )
    return {
        "whole_light_order_loss_weight": order_loss("light"),
        "whole_heavy_order_loss_weight": order_loss("heavy"),
        "whole_conditional_closed_error_budget": pieces,
        "positive_Q_order_proof": "For v=cosh(eta)>1 and delta>=0, write Q_(2+delta)(v)/Q2(v) as E[exp(-delta*(eta-log u))] with positive normalized weight u^2/sqrt((1-u)(1-exp(-2eta)*u)),0<u<1. Its mean of -log u is <=mean((1-u)/u)<=1/4. The last bound follows from the beta(3,1/2) value and the opposite monotonicity covariance with (1-d*u)^(-1/2). Therefore0<=1-Q_(2+delta)/Q2<=delta*(acosh(v)+1/4), for the whole Legendre function, not its leading large-v term.",
        "conditional_cut_error": "Assume one isolated factorized even pole, its continued two-channel unitarity law, alpha(T)=2+delta(T) with0<=delta(T)<=T/L, and the normalized nonkinematic heavy/light residue ratio zeta within epsilon of1. Kinematic residue ratio is[(T-4a)/(T-4mu)]^(alpha/2). Relative to the exact spin2 cut, the full factor is zeta*[(T-4a)/(T-4mu)]^(delta/2)*Q_(2+delta)(v_a)/Q2(v_a). Its absolute defect is <=epsilon+delta*order_loss_a(T). For the light channel zeta=1 identically. Threshold logarithms are integrable; no uniformly small relative error at the threshold is assumed.",
        "finite_transfer_contour": "Let ell(t) be the selected one-loop correction to log[f(t)/f(0)], not an unnormalized residue. Assume its stated slit disk toU>4n has no other omitted cuts or singularities and an arc norm <=B. Cauchy's derivative formula gives ell'(0)=integral rho_ell/(pi*T^2)+arc. The closed bound is |ell'(0)-2F1prime(0)| <=B/U+known_tail+declared_residue_error+declared_trajectory_error. Subtracting the normalized value at zero is a definition, not a new finite matching choice. Divide the mismatch bykappa for this selected contribution to the finite v^2 contour coefficient.",
        "bound_derivation": "Both pair weights satisfy2ImF1/(pi*T^2)<=g^2*(T-4a)^2/(240pi^2*n^3*T^2), using the exact Bose P2 coefficient majorant and both denominator gaps. For T<=U, acosh(v_light)<=log[(2U+4n)/(T-4mu)] and acosh(v_heavy)<=log[(2U+4n)/(T-4n)]. The heavy kinematic factor adds half log[U/(T-4n)]. Multiply delta<=T/L, use(T-4a)^2/T<=T-4a, and integrate the full logarithmic primitives to obtain the two displayed terms. The known tail and total moment use the entire previous massive graphs.",
        "not_established": "No model values for U,L,epsilon or B, no absence of other Regge singularities, no actual complex-J unitary completion, no analytic trajectory term, and no high-energy complex-nu remainder or all-order massless infrared limit is established. Integer-spin unitarity and an arbitrarily small graviton coupling cannot fix these inputs. This is a proved conditional error map, not a finite-gravity P8 verdict.",
        "checks": checks,
        "gates": {
            "whole_Q_function_order_error_not_only_asymptotics": True,
            "threshold_nonuniformity_integrated_not_suppressed": True,
            "two_mass_kinematic_residue_ratio_retained": True,
            "normalized_log_residue_not_unsubtracted_f_slope": True,
            "unknown_transfer_arc_explicit": True,
            "no_numerical_Regge_matching_inputs_invented": True,
            "complex_nu_and_trajectory_and_other_cuts_still_open": True,
        },
    }
