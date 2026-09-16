"""Fixed three-ball dimensional continuation with a retained trace term."""

from functools import cache

import sympy as s

from . import source

RADIAL = s.Symbol("radial_squared", nonnegative=True)
EP = source.EP


def require_direction(direction):
    n = s.Matrix(direction)
    if n.shape != (3, 1):
        raise ValueError("Require a three-dimensional unit direction")
    for value in n:
        source.recoil.exact_real(value)
    if s.simplify(n.dot(n) - 1) != 0:
        raise ValueError("Require an exact unit direction")
    return n


def require_radial(value):
    r = source.recoil.exact_real(value)
    if r < 0 or r > 1:
        raise ValueError("Require radial coordinate in[0,1]")
    return r


def transverse_current(energy, quanta, outgoing, direction, radial=1):
    points, rays, _ = source.recoil.momenta(energy, quanta, outgoing)
    n, r = require_direction(direction), require_radial(radial)
    v = s.Matrix([*(r * n), s.sqrt(1 - r * r)])
    projector = s.eye(4) - v * v.T
    current = s.zeros(4)
    for p in (*points, *rays):
        spatial = s.Matrix([*p[1:4, 0], 0])
        den = p[0] - r * p[1:4, 0].dot(n)
        if den == 0:
            # At a null collinear boundary point select a bounded representative.
            # The angular coefficient is insensitive to this zero-area set.
            if source.recoil.dot(p, p) != 0:
                raise ValueError("Massive Doppler gap lost")
            continue
        projected = projector * spatial
        current += projected * projected.T / den
    return current.applyfunc(s.simplify)


def contractions(current, epsilon=0):
    e = s.sympify(epsilon)
    trace_square = s.trace(current) ** 2
    return s.trace(current * current) - trace_square / (2 + 2 * e), trace_square


def beta_multiplier(epsilon=EP):
    e = s.sympify(epsilon)
    return s.gamma(s.Rational(3, 2) + e) / (s.gamma(s.Rational(3, 2)) * s.gamma(1 + e))


def radial_normalization(epsilon=EP):
    e = s.sympify(epsilon)
    return e * beta_multiplier(e)


def beta_moment(order, epsilon=EP):
    source.require_order(order)
    e = s.sympify(epsilon)
    return s.prod(
        (s.Rational(3, 2) + j) / (s.Rational(3, 2) + e + j) for j in range(order)
    )


def dimensional_first(angular_zero, boundary_zero, trace_boundary):
    return s.sympify(trace_boundary) / 2 + s.Integral(
        s.sqrt(RADIAL) * (angular_zero - boundary_zero) / (1 - RADIAL), (RADIAL, 0, 1)
    )


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = (
            value.applyfunc(s.factor)
            if isinstance(value, s.MatrixBase)
            else s.factor(s.expand_func(value))
        )

    e = s.Symbol("epsilon", positive=True)
    a = s.Rational(3, 2)
    for m in range(1, 9):
        mean = beta_moment(m, e)
        put("beta_moment_boundary_" + str(m), mean.subs(e, 0) - 1)
        put(
            "beta_moment_first_" + str(m),
            s.diff(mean, e).subs(e, 0) + sum(1 / (a + j) for j in range(m)),
        )
    put("beta_multiplier_at_zero", beta_multiplier(0) - 1)
    put(
        "beta_multiplier_log_derivative",
        s.diff(beta_multiplier(e), e).subs(e, 0) - (2 - 2 * s.log(2)),
    )
    put(
        "projector_first_trace_term",
        s.diff(-1 / (2 + 2 * e), e).subs(e, 0) - s.Rational(1, 2),
    )
    put(
        "projector_exact_difference",
        -1 / (2 + 2 * e) + s.Rational(1, 2) - e / (2 * (1 + e)),
    )
    u = s.Symbol("radial_chart", positive=True)
    t = (1 - u * u) ** 2
    put(
        "smooth_radial_substitution",
        -s.diff(t, u) * (1 - u * u) / (1 - t)
        - 4 * (1 - u * u) ** 2 / (u * (2 - u * u)),
    )
    for dim in (5, 6):
        size = dim - 1
        n = s.Matrix([s.Rational(3, 5), s.Rational(4, 5), *([0] * (size - 2))])
        P = s.eye(size) - n * n.T
        put("D" + str(dim) + "_projector_idempotent", P * P - P)
        put("D" + str(dim) + "_projector_trace", s.trace(P) - (dim - 2))
        entries = s.symbols("S0:" + str(size * (size + 1) // 2))
        spatial = s.zeros(size)
        at = 0
        for i in range(size):
            for j in range(i, size):
                spatial[i, j] = spatial[j, i] = entries[at]
                at += 1
        J = s.zeros(dim)
        J[1:dim, 1:dim] = spatial
        J[0, 0] = (n.T * spatial * n)[0]
        J[1:dim, 0] = spatial * n
        J[0, 1:dim] = (spatial * n).T
        eta = s.diag(1, *([-1] * size))
        light = s.Matrix([1, *n])
        put("D" + str(dim) + "_generic_conserved_current", J * eta * light)
        transverse = P * spatial * P
        put(
            "D" + str(dim) + "_generic_covariant_equals_projected",
            s.trace(J * eta * J * eta)
            - s.trace(J * eta) ** 2 / (2 + 2 * e)
            - s.trace(transverse * transverse)
            + s.trace(transverse) ** 2 / (2 + 2 * e),
        )
        put(
            "D" + str(dim) + "_generic_trace_relation",
            s.trace(J * eta) + s.trace(transverse),
        )
    # An actual fully conserved state, not only an abstract transverse tensor.
    w = s.Rational(1, 16)
    quanta = [s.Matrix([w, w, 0, 0]), s.Matrix([w, -w, 0, 0])]
    direction = s.Matrix([0, 1, 0])
    outgoing = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    current = transverse_current(
        s.Rational(5, 4), quanta, outgoing, direction, s.Rational(3, 5)
    )
    points, rays, _ = source.recoil.momenta(s.Rational(5, 4), quanta, outgoing)
    eta = s.diag(1, -1, -1, -1, -1)
    light = s.Matrix([1, 0, s.Rational(3, 5), 0, s.Rational(4, 5)])
    ps = [s.Matrix([*p, 0]) for p in (*points, *rays)]
    J = sum((p * p.T / (p.T * eta * light)[0] for p in ps), s.zeros(5)).applyfunc(
        s.simplify
    )
    put("actual_two_radiation_full_current_Ward", J * eta * light)
    put(
        "actual_two_radiation_covariant_projection",
        s.trace(J * eta * J * eta)
        - s.trace(J * eta) ** 2 / (2 + 2 * e)
        - contractions(current, e)[0],
    )
    return {
        "whole_fixed_ball_measure": "For physical external spatial momenta inR3, t=|n_parallel|^2 has normalized Beta(3/2,e) law inD=4+2e. The remaining angular variable is the ordinary unit two-sphere. Work on this fixed ball, not a varying-dimensional domain.",
        "whole_current": "Embed external vectors intoR4, n_t=(sqrt(t)*n,sqrt(1-t)),P_t=I-n_t*n_t^T,T_t=sum_A P_t*pvec_A*pvec_A^T*P_t/(p_A^0-sqrt(t)*pvec_A.n). F_e=tr(T_t^2)-tr(T_t)^2/(2+2e). Whole momentum conservation identifies this with the covariant full-D current.",
        "whole_dimensional_first": "K1=one_half*mean_S2[tr(T_1)^2]+integral_0^1 sqrt(t)/(1-t)*[mean_S2 F0(t)-K0]dt. The radial difference must be retained, as must the projector derivative.",
        "whole_dimensional_scheme": "Four-dimensional external momenta, full D-dimensional additional-soft projector and angular phase convention ofS296. The finite universal factor is specified in this scheme. No evanescent hard tree/loop terms are replaced by this factor. The interior noninteger-D contraction need not be positive.",
        "checks": checks,
        "gates": {
            "positive_normalized_fixed_ball_measure_not_fictitious_helicities": True,
            "complete_covariant_current_equals_projected_current": True,
            "generic_conserved_tensors_and_actual_recoiled_state_checked": True,
            "boundary_trace_and_radial_finite_terms_both_retained": True,
            "bounded_collinear_representatives_only_on_zero_area_sets": True,
            "full_D_trace_coefficient_not_set_to_D4_early": True,
            "hard_evanescent_scheme_terms_not_inferred": True,
        },
    }
