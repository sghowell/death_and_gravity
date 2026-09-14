"""All mixed implicit-root contacts and the original physical-volume chain."""

from functools import cache

import sympy as s

from . import source as q


def first_second(d, gradient, mixed, hessian, d2):
    gradient = s.Matrix(gradient)
    mixed = s.Matrix(mixed)
    first = -gradient / d
    second = (
        -(s.Matrix(hessian) + mixed * first.T + first * mixed.T + d2 * first * first.T)
        / d
    )
    return first, second


def third_direction(d, c, b, a, e, f, g, k, r=0):
    first = -c / d
    second = -(a + 2 * b * first + e * first**2) / d
    third = (
        -(
            r
            + 3 * k * first
            + 3 * g * first**2
            + f * first**3
            + 3 * b * second
            + 3 * e * first * second
        )
        / d
    )
    return first, second, third


def pullback(first, second, jacobian, coordinate_hessians):
    first = s.Matrix(first)
    jacobian = s.Matrix(jacobian)
    gradient = jacobian.T * first
    hessian = jacobian.T * s.Matrix(second) * jacobian
    for coefficient, contact in zip(first, coordinate_hessians, strict=True):
        hessian += coefficient * s.Matrix(contact)
    return gradient, hessian


@cache
def jets():
    d = s.Symbol("nonzero_full_C_N", real=True, nonzero=True)
    e = s.Symbol("whole_C_NN", real=True)
    gradient = s.Matrix(s.symbols("whole_C_z0:12", real=True))
    mixed = s.Matrix(s.symbols("whole_C_Nz0:12", real=True))
    hessian = s.zeros(12)
    for i in range(12):
        for j in range(i, 12):
            hessian[i, j] = hessian[j, i] = s.Symbol(
                "whole_C_z_" + str(i) + "_" + str(j), real=True
            )
    first, second = first_second(d, gradient, mixed, hessian, e)
    c, b, a, f, g, k, r = s.symbols(
        "C_direction C_N_direction C_direction2 C_NNN C_NN_direction C_N_direction2 C_direction3",
        real=True,
    )
    n1, n2, n3 = third_direction(d, c, b, a, e, f, g, k, r)
    t, x = s.symbols("direction_parameter lapse_displacement", real=True)
    taylor = (
        d * x
        + c * t
        + (e * x * x + 2 * b * x * t + a * t * t) / 2
        + (f * x**3 + 3 * g * x * x * t + 3 * k * x * t * t + r * t**3) / 6
    )
    ansatz = n1 * t + n2 * t * t / 2 + n3 * t**3 / 6
    expanded = s.Poly(s.expand(taylor.subs(x, ansatz)), t)
    return {
        "whole_all_twelve_first_implicit_derivatives": first,
        "whole_all_144_mixed_second_implicit_derivatives": second,
        "whole_general_first_second_third_directional_derivatives": [n1, n2, n3],
        "whole_directional_third_derivative_convention": "c=C_z[v], b=C_Nz[v], a=C_zz[v,v], e=C_NN, f=C_NNN, g=C_NNz[v], k=C_Nzz[v,v], r=C_zzz[v,v,v]. On the twelve stated invariant coordinates the literal full C is quadratic, so r=0. In nonlinear canonical field coordinates r need not vanish. The full symmetric third tensor is recovered by polarization, with all mixed contacts retained.",
        "whole_canonical_pullback_rule": "For any actual canonical finite-jet chart w with invariant map z(w), N_,A=sum_i N_,i z^i_,A and N_,AB=sum_ij N_,ij z^i_,A z^j_,B+sum_i N_,i z^i_,AB. The second term includes the original energy/gradient/vector and density contacts. The invariant variables are not assumed linear or Gaussian. A nonlinear chart's third contacts must likewise be kept.",
        "checks": {
            "whole_first_constraint_chain": (d * first + gradient).applyfunc(s.factor),
            "whole_all_mixed_second_constraint_chain": (
                d * second
                + hessian
                + mixed * first.T
                + first * mixed.T
                + e * first * first.T
            ).applyfunc(s.factor),
            "whole_mixed_second_symmetry": second - second.T,
            "whole_first_directional_chain": s.factor(expanded.coeff_monomial(t)),
            "whole_second_directional_chain": s.factor(expanded.coeff_monomial(t**2)),
            "whole_third_directional_chain": s.factor(expanded.coeff_monomial(t**3)),
        },
        "gates": {
            "all_twelve_density_jet_directions_retained": first.shape == (12, 1),
            "all_144_mixed_root_contacts_retained": second.shape == (12, 12),
            "nonlinear_canonical_map_second_contacts_mandatory": True,
            "third_direction_includes_nonzero_canonical_Czzz_if_present": n3.has(r),
        },
    }


@cache
def physical():
    U = q.R ** -s.Rational(3, 4)
    U1 = s.diff(U, q.N)
    U2 = s.diff(U, q.N, 2)
    v0, v1, v2, n1, n2 = s.symbols(
        "hat_v0 hat_v_first hat_v_second full_N_first full_N_second", real=True
    )
    z = s.Symbol("variation", real=True)
    u0, u1, u2 = s.symbols("whole_U whole_U_N whole_U_NN", real=True)
    first = s.exp(3 * v0) * (3 * u0 * v1 + u1 * n1)
    second = s.exp(3 * v0) * (
        9 * u0 * v1**2 + 6 * u1 * v1 * n1 + u2 * n1**2 + 3 * u0 * v2 + u1 * n2
    )
    polynomial = u0 + u1 * (n1 * z + n2 * z * z / 2) + u2 * n1 * n1 * z * z / 2
    scalar = s.exp(3 * (v0 + v1 * z + v2 * z * z / 2)) * polynomial
    # General two-canonical-direction chain includes nonlinear invariant contacts.
    x, y = s.symbols("canonical_x canonical_y", real=True)
    invariant = s.Matrix([x + y * y, x * y + 2 * y, x * x - y + x * y * y])
    outer = (
        invariant[0] ** 2
        + invariant[0] * invariant[1]
        + invariant[1] * invariant[2]
        + invariant[2] ** 3
    )
    zi = s.symbols("invariant0:3", real=True)
    independent = zi[0] ** 2 + zi[0] * zi[1] + zi[1] * zi[2] + zi[2] ** 3
    bind = dict(zip(zi, invariant, strict=True))
    ng = s.Matrix([s.diff(independent, z).subs(bind, simultaneous=True) for z in zi])
    nh = s.hessian(independent, zi).subs(bind, simultaneous=True)
    mapped1, mapped2 = pullback(
        ng, nh, invariant.jacobian((x, y)), [s.hessian(z, (x, y)) for z in invariant]
    )
    return {
        "whole_original_physical_volume_factor_and_full_lapse_derivatives": [U, U1, U2],
        "whole_complete_physical_volume_first_variation": first,
        "whole_complete_physical_volume_second_variation": second,
        "whole_additional_nonlinear_lapse_second_contact": s.exp(3 * v0) * u1 * n2,
        "whole_physical_second_contact_boundary": "These are complete classical variations of exp(3 vhat) R_full(Nstar)^(-3/4). The full lapse second variation, and the log-scale second variation under nonlinear charts, remain. They are absent from a substitution using only the linear reference lapse. A future canonical ordered-state contraction must include the entire pulled-back Hessian, not only U_NN Var(n_linear). No actual interacting quantum onepoint, subtraction, ordering or Gaussian law for density invariants is defined here.",
        "checks": {
            "whole_original_physical_U_binding": U - q.R ** -s.Rational(3, 4),
            "whole_complete_physical_first_chain": s.factor(
                s.diff(scalar, z).subs(z, 0) - first
            ),
            "whole_complete_physical_second_chain": s.factor(
                s.diff(scalar, z, 2).subs(z, 0) - second
            ),
            "whole_independent_nonlinear_invariant_first_pullback": (
                mapped1 - s.Matrix([s.diff(outer, x), s.diff(outer, y)])
            ).applyfunc(s.factor),
            "whole_independent_nonlinear_invariant_second_pullback": (
                mapped2 - s.hessian(outer, (x, y))
            ).applyfunc(s.factor),
        },
        "gates": {
            "full_nonlinear_lapse_second_contact_not_dropped": second.has(n2),
            "full_metric_chart_second_contact_not_dropped": second.has(v2),
            "S261_shifted_R_binding_not_used": True,
            "classical_chain_not_interacting_mean": True,
        },
    }
