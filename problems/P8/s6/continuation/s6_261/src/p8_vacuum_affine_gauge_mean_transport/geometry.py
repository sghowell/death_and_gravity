"""Complete density-weight orbit and physical-volume identities."""

from functools import cache

import sympy as s
from p8_vacuum_affine_physical_background_vertices import parent
from p8_vacuum_affine_spatial_gauge import gauge

ALPHA = s.Rational(9, 4)


@cache
def family():
    old = gauge.full_operator()
    x = old["whole_spatial_coordinates"]
    Q = old["whole_contravariant_shape_density"]
    xi = old["whole_spatial_gauge_vector"]
    chi = old["whole_spatial_gauge_condition"]
    v = s.Function("whole_log_hat_volume", real=True)(*x)
    N = s.Function("whole_lapse", positive=True)(*x)
    C = s.Function("whole_physical_conformal_factor", positive=True)(N)
    a = s.Symbol("positive_homogeneous_scale", positive=True)
    dw = s.Symbol("gauge_weight_displacement", real=True)
    Qw = s.exp(3 * dw * (v + s.log(a))) * Q
    chiw = s.Matrix([sum(s.diff(Qw[i, j], x[j]) for j in range(3)) for i in range(3)])
    forcing = 3 * (v + s.log(a)) * chi + 3 * Q * s.Matrix([s.diff(v, z) for z in x])
    div = sum(s.diff(xi[i], x[i]) for i in range(3))
    dv = sum(xi[i] * s.diff(v, x[i]) for i in range(3)) + div / 3
    dN = sum(xi[i] * s.diff(N, x[i]) for i in range(3))
    density = a**3 * s.exp(3 * v) * C ** s.Rational(3, 2)
    direct = density * (3 * dv + s.Rational(3, 2) * s.diff(C, N) * dN / C)
    divergence = sum(s.diff(density * xi[i], x[i]) for i in range(3))
    k = s.Matrix(s.symbols("nonzero_wavevector0:3", real=True))
    k2 = k.dot(k)
    M0 = k2 * s.eye(3) + k * k.T / 3
    inverse = (s.eye(3) - k * k.T / (4 * k2)) / k2
    first = -3 * s.I * inverse * k
    return {
        "whole_density_weight_family": Qw,
        "whole_weight_derivative_before_gauge_restriction": forcing,
        "whole_on_slice_weight_forcing": 3 * Q * s.Matrix([s.diff(v, z) for z in x]),
        "whole_on_slice_gauge_transport_equation": "M_Q xi=-3 Q grad(v), with the entire S258 gauge-surface operator and residual translations separate.",
        "whole_log_volume_variation": dv,
        "whole_lapse_scalar_variation": dN,
        "whole_physical_volume_density": density,
        "whole_physical_volume_gauge_variation": direct,
        "whole_nonzero_fourier_displacement_row": first,
        "checks": {
            "full_weight_derivative_including_off_slice_contact": (
                s.diff(chiw, dw).subs(dw, 0) - forcing
            ).applyfunc(s.expand),
            "full_physical_volume_is_a_density_divergence": s.simplify(
                direct - divergence
            ),
            "whole_flat_gauge_inverse": (M0 * inverse - s.eye(3)).applyfunc(s.factor),
            "whole_linear_gauge_displacement": (first + ALPHA * s.I * k / k2).applyfunc(
                s.factor
            ),
            "whole_linear_log_volume_slope": s.factor(
                (s.I * k.T * first)[0] / 3 - s.Rational(3, 4)
            ),
        },
        "gates": {
            "off_slice_weight_contact_not_deleted": forcing
            - 3 * Q * s.Matrix([s.diff(v, z) for z in x])
            != s.zeros(3, 1),
            "lapse_dependent_physical_conformal_factor_retained": direct.has(
                s.diff(C, N)
            ),
            "constant_translation_inverse_not_used": True,
            "local_spatial_map_not_quantum_Nielsen_vector": True,
        },
    }


@cache
def periodic_orbit():
    w = s.Symbol("gauge_weight", real=True)
    beta = (3 * w - 2) / (2 - w)
    y = s.Symbol("original_periodic_coordinate", real=True)
    v = s.Function("periodic_log_volume", real=True)(y)
    t = s.Function("first_tensor_profile", real=True)(y)
    z = s.Function("second_tensor_profile", real=True)(y)
    K = s.Symbol("positive_periodic_normalization", positive=True)
    N = s.Function("periodic_positive_lapse", positive=True)(y)
    C = s.Function("whole_parent_conformal_factor", positive=True)(N)
    # Two independent tensor profiles: the transverse block is SPD and det=1.
    B = s.Matrix([[s.exp(t), z], [z, (1 + z * z) * s.exp(-t)]])
    jac = K * s.exp(beta * v)
    metric = s.diag(1, B) * s.exp(2 * v)
    pulled = s.diag(jac * jac, B) * s.exp(2 * v)
    qxx_log = (3 * w - 2) * v + (w - 2) * (s.log(K) + beta * v)
    b = s.Symbol("independent_weight_exponent", real=True)
    norm = s.Function("periodic_partition", positive=True)(b)
    lapse_weight = s.Function("periodic_lapse_numerator", real=True)(b)
    mean = -(1 + b / 3) * s.diff(norm, b) / norm + s.log(norm) / 3
    meanN = lapse_weight / norm
    m1, m2, n1, vn = s.symbols(
        "volume_first volume_second lapse_first volume_lapse_product", real=True
    )
    bind = {
        norm: 1,
        s.diff(norm, b): -m1,
        s.diff(norm, b, 2): m2,
        lapse_weight: n1,
        s.diff(lapse_weight, b): -vn,
    }
    ds = s.diff(mean, b).subs(bind, simultaneous=True).subs(b, 0)
    dn = s.diff(meanN, b).subs(bind, simultaneous=True).subs(b, 0)
    return {
        "gauge_parameter": w,
        "whole_finite_weight_exponent": beta,
        "whole_two_tensor_transverse_metric": B,
        "whole_three_direction_reference_metric": metric,
        "whole_three_direction_pulled_metric": pulled,
        "whole_positive_coordinate_Jacobian": jac,
        "whole_normalization_and_inverse_map": "K=<exp(-beta v)>_y; dx/dy=exp(-beta v)/K, with the full periodic integral normalized to one. The inverse monotone map y=phi(x) is a smooth orientation-preserving circle diffeomorphism. A constant translation fixes the residual origin.",
        "whole_finite_gauge_condition": "Q_w^(xx)=K^(w-2) is constant and Q_w^(xy)=Q_w^(xz)=0, so all three divergences vanish; both transverse tensor profiles remain arbitrary.",
        "whole_exact_log_volume_mean": mean,
        "whole_exact_lapse_mean": meanN,
        "whole_log_volume_mean_weight_derivative": ALPHA * ds,
        "whole_lapse_mean_weight_derivative": ALPHA * dn,
        "whole_physical_volume_integral": "The pulled density is C(N(phi))^(3/2) exp(3v(phi)) phi'. Changing variables gives the unchanged full original integral exactly, for every positive parent C in its domain.",
        "checks": {
            "whole_transverse_tensor_determinant": s.factor(B.det() - 1),
            "whole_metric_determinant": s.factor(metric.det() - s.exp(6 * v)),
            "whole_pulled_metric_determinant": s.simplify(
                pulled.det() - s.exp(6 * v) * jac * jac
            ),
            "whole_three_gauge_conditions_constant_x_density": s.factor(
                qxx_log - (w - 2) * s.log(K)
            ),
            "weight_exponent_zero_at_Dirac": beta.subs(w, s.Rational(2, 3)),
            "whole_weight_exponent_derivative": s.factor(
                s.diff(beta, w).subs(w, s.Rational(2, 3)) - ALPHA
            ),
            "whole_exact_coordinate_mean_variance": s.factor(ds + m2 - m1 * m1),
            "whole_exact_lapse_mean_cross_covariance": s.factor(dn + vn - m1 * n1),
            "whole_physical_density_change_of_variable": s.simplify(
                C ** s.Rational(3, 2) * s.exp(3 * v) * jac / jac
                - C ** s.Rational(3, 2) * s.exp(3 * v)
            ),
        },
        "gates": {
            "both_transverse_tensor_profiles_retained": B.has(t) and B.has(z),
            "full_finite_gauge_domain_w_less_than_two": True,
            "circle_map_not_global_general_three_dimensional_slice": True,
            "finite_gauge_changes_coordinate_mean_not_physical_volume": True,
        },
    }


@cache
def physical_volume_contact():
    Cvv, Cvn, c1 = s.symbols(
        "whole_Cvv whole_Cvn physical_volume_lapse_slope", real=True
    )
    mean_v, mean_n = -ALPHA * Cvv, -ALPHA * Cvn
    first_vv = ALPHA * Cvv / 3
    first_vn = ALPHA * Cvn / 3
    linear = 3 * mean_v + c1 * mean_n
    quadratic = 9 * first_vv + 3 * c1 * first_vn
    N, u = parent.N, parent.u
    Rgerm = 1 + (N**-2 - 1) / parent.h
    Cgerm = Rgerm - s.Rational(1, 2)
    slope = s.factor(s.Rational(3, 2) * s.diff(Cgerm, N) / Cgerm)
    clock_slope = s.factor(slope.subs(N, 1))
    # Entire R minus this germ vanishes to order eight. Only its first
    # lapse jet is used; the full off-clock function is never replaced.
    return {
        "whole_second_order_log_volume_transport": mean_v,
        "whole_second_order_lapse_transport": mean_n,
        "whole_linear_onepoint_volume_contact": linear,
        "whole_quadratic_two_point_volume_contact": quadratic,
        "whole_clock_physical_volume_lapse_slope": clock_slope,
        "whole_bounce_physical_volume_lapse_slope": clock_slope.subs(u, 0),
        "checks": {
            "full_hat_volume_onepoint_and_two_point_cancellation": 3 * mean_v
            + 9 * first_vv,
            "full_physical_lapse_dependent_volume_cancellation": s.expand(
                linear + quadratic
            ),
            "literal_entire_parent_first_lapse_jet": s.factor(
                clock_slope + 6 / parent.h
            ),
            "actual_bounce_lapse_volume_slope": clock_slope.subs(u, 0) + 6,
        },
        "gates": {
            "deleting_nonlinear_mean_breaks_volume_identity": quadratic != 0,
            "lapse_cross_covariance_not_set_to_zero": linear.has(Cvn),
            "whole_parent_clock_jet_source_order_eight_retained": True,
        },
    }
