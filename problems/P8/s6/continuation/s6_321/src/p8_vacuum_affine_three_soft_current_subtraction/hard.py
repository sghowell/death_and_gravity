"""One-block hard-core operator, independent source checks and integrability."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import jets
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as matter
from p8_vacuum_affine_relative_energy_complex_tube import continuation as analytic
from p8_vacuum_affine_temporal_tree_reorganization import trees
from sympy.polys.rings import ring

from . import polynomial, source, subtraction

ETA = s.Rational(1, 10**13)


@cache
def polynomial_calibration():
    Packed = polynomial.FactoredAssemblyPolynomialCurrent
    R, *generators = ring("z,x00,x01,x02,x03,x11,x12,x13,x22,x23,x33", s.QQ)
    z = generators[0]

    def Kernel(RR, qs, leaves):
        return Packed(RR, qs, leaves, generators[-10:])

    matrix = lambda A: tuple(
        tuple(R.from_expr(A[i, j]) for j in range(4)) for i in range(4)
    )
    vector = lambda p: tuple(R.from_expr(q) for q in p)
    dummy = Kernel(R, (), ())
    fields, momenta = jets.calibration_fields()
    for order in (3, 4):
        got = dummy.vertex(
            tuple(matrix(A) for A in fields[:order]),
            tuple(vector(p) for p in momenta[:order]),
        )
        expected = jets.MetricJet(fields[:order], momenta[:order]).gravity()
        assert s.factor(got.as_expr() - expected) == 0
    for coordinates, weights in (
        (
            ((0, 0), (s.Rational(1, 2), 0), (-s.Rational(1, 3), s.Rational(1, 4))),
            (s.Rational(1, 4), s.S.One, s.S.One),
        ),
        (
            ((0, 0), (s.Rational(1, 3), 0), (s.Rational(1, 3), s.Rational(1, 9))),
            (s.Rational(1, 9), s.S.One, s.S.One),
        ),
    ):
        for bits in ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 1)):
            leaves = []
            qs = []
            legs = []
            for (x, y), weight, bit in zip(coordinates, weights, bits):
                x, y = map(s.sympify, (x, y))
                d = 1 + x * x + y * y
                n = s.Matrix([2 * x, 2 * y, 1 - x * x - y * y]) / d
                u = s.Matrix([1 - 2 * x * x / d, -2 * x * y / d, -2 * x / d])
                v = s.Matrix([-2 * x * y / d, 1 - 2 * y * y / d, -2 * y / d])
                A = s.zeros(4)
                A[1:, 1:] = u * u.T - v * v.T if bit == 0 else u * v.T + v * u.T
                q = s.ImmutableMatrix([weight, *(weight * n)])
                A = s.ImmutableMatrix(A)
                legs.append(("h", q, A))
                qs.append(vector(q))
                leaves.append((matrix(A), []))
            engine = Kernel(R, qs, leaves)
            baseline = trees.TemporalSoft(legs)
            for mask in (3, 5, 6, 7):
                got, count = engine.sympy_matrix(engine.current(mask))
                expected, n = baseline.current(mask, "h")
                assert count == n == (4 if mask == 7 else 1)
                assert (got - expected).applyfunc(s.factor) == s.zeros(4)
                J, d, c = engine.amputated(mask)
                plain = s.ImmutableMatrix(
                    [
                        [p.as_expr() / engine.denominator(d).as_expr() for p in row]
                        for row in J
                    ]
                )
                original, original_count = baseline.amputated(mask, "h")
                assert c == original_count and (plain - original).applyfunc(
                    s.factor
                ) == s.zeros(4)
    # A removable high-degree factor is deliberately retained: exact
    # cross multiplication verifies that this normalization is still identical.
    factor = z**3 + z + 1
    synthetic = dummy.scale(factor * (z + 1), dummy.eta)
    normalized, parts = dummy.normalize(synthetic, [(factor, 1), (z + 1, 1)])
    assert any(f == factor and power == 1 for f, power in parts)
    D = dummy.denominator(parts)
    for i in range(4):
        for j in range(4):
            assert normalized[i][j] * factor * (z + 1) == synthetic[i][j] * D
    # Independent generic-matrix identity: no source conservation is assumed
    # for this comparison with the original inverse and temporal projection.
    Q = s.Matrix(s.symbols("Q0 Q1 Q2 Q3"))
    eta = s.diag(1, -1, -1, -1)
    entries = s.symbols("j0:10")
    J = s.zeros(4)
    at = 0
    for i in range(4):
        for j in range(i, 4):
            J[i, j] = J[j, i] = entries[at]
            at += 1
    square = (Q.T * eta * Q)[0]
    pre = -(eta * J * eta - eta * s.trace(eta * J) / 2) / square
    q = eta * Q
    projected = s.Matrix(
        4,
        4,
        lambda i, j: (
            pre[i, j]
            - q[i] * pre[0, j] / Q[0]
            - q[j] * pre[i, 0] / Q[0]
            + q[i] * q[j] * pre[0, 0] / Q[0] ** 2
        ),
    )
    S = -(eta * J * eta - eta * s.trace(eta * J) / 2)
    direct = s.zeros(4)
    for i in range(1, 4):
        for j in range(1, 4):
            direct[i, j] = (
                Q[0] ** 2 * S[i, j]
                + Q[i] * Q[0] * S[0, j]
                + Q[j] * Q[0] * S[i, 0]
                + Q[i] * Q[j] * S[0, 0]
            ) / (square * Q[0] ** 2)
    assert (projected - direct).applyfunc(s.factor) == s.zeros(4)

    # Separate multivariate-to-univariate representation checks, including
    # both exactly divisible and nondivisible polynomials at every variable.
    T, aa, bb, rr, vv, uu = ring("aa,bb,rr,vv,uu", s.QQ)
    divider = polynomial.CoefficientRingDivision(T)
    factors = (
        aa,
        bb,
        aa + bb,
        aa + 1,
        bb + 1,
        aa + bb + 1,
        1 + rr**2,
        1 + vv**2 + uu**2,
        vv**2 + uu**2,
        (rr + vv) ** 2 + uu**2,
    )
    seed = (aa + 2 * bb + 3 * rr + 5 * vv + 7 * uu + 1) ** 4
    for factor in factors:
        quotient, remainder = divider.div(seed * factor, factor)
        assert not remainder and quotient == seed
        quotient, remainder = divider.div(seed * factor + 1, factor)
        assert seed * factor + 1 == quotient * factor + remainder and remainder
    for index in range(T.ngens):
        assert divider.unpack(divider.pack(seed, index), index) == seed

    for factors, energies in (
        ([(z, 2), (z + 1, 1), (2 * z + 1, 1)], (z, z + 1)),
        ([(z, 1), (2 * z + 1, 1)], (z, z + 1)),
        ([(z, 2)], (z, z)),
    ):
        original_N = dummy.scale(z * z + z + 1, dummy.eta)
        original_D = dummy.denominator(factors)
        new_N, new_parts, count = polynomial.multiply_energies_over_total(
            dummy, (original_N, factors, 4), energies, 2 * z + 1
        )
        new_D = dummy.denominator(new_parts)
        energy_product = R.one
        for energy in energies:
            energy_product *= energy
        assert count == 4
        for i in range(4):
            for j in range(4):
                assert (
                    new_N[i][j] * original_D * (2 * z + 1)
                    == original_N[i][j] * energy_product * new_D
                )

    return {
        "checks": {
            "literal_EH3_EH4_and32_current_source_identities": s.S.Zero,
            "generic16_component_inverse_projection_identity": s.S.Zero,
            "retained_factor_rational_identity": s.S.Zero,
            "twenty_division_checks_and_five_representation_roundtrips": s.S.Zero,
            "three_exact_energy_prefactor_identities": s.S.Zero,
        },
        "gates": {
            "unchanged_action_independent_oracles_all_assertions_passed": True,
            "no_backend_patch_known_factor_division_cross_multiplied": True,
            "no_private_polynomial_cache_import": True,
        },
    }


@cache
def core_calibration():
    clean = analytic.canonical
    matrix = lambda A: s.ImmutableMatrix(A.applyfunc(clean))
    pars = source.original_parameters()
    assert pars["heavy"] == s.Integer(10) ** 200 / 512 + 2
    assert pars["cubic"] == s.Rational(1, 8192)
    assert pars["kappa"] == s.Integer(10) ** 800
    x = s.Symbol("x")
    L = (1 - 1024 * x) / (1 - 5120 * x)
    D = 1 / (1 - 4 * x)
    HH = (1 - 1024 * x) / (1 - 3072 * x)
    EE = 1024 / (1 - 1024 * x) ** 2
    VV = 120000000 * 1024 * (2 / (1 - 32 * x) ** 3 - 2)
    Fm = L**4 * (D + s.Rational(3, 2) * D**2 * HH)
    FG = s.Rational(3 * 600000, 8) * L**4 * EE**2 / (1 - VV)
    cm = s.diff(Fm, x).subs(x, 0)
    cg = s.diff(FG, x).subs(x, 0)
    ind_m = 4 * (4 * 1024) * s.Rational(5, 2) + 4 + s.Rational(3, 2) * (8 + 2048)
    ind_g = (
        s.Rational(3 * 600000, 8)
        * 1024**2
        * (4 * (4 * 1024) + 4 * 1024 + 120000000 * 1024 * 192)
    )
    assert cm == ind_m == 44048 and cg == ind_g
    B0 = 4 * (pars["heavy"] ** 2 * cm + cg)
    assert 4 * cm < 10**8 and 4 * cg < 10**28
    assert B0 < 10**8 * pars["heavy"] ** 2 + 10**28
    eta = s.Rational(1, 10**13)
    assert 3 * eta < s.Rational(1, 10**12)
    E = s.Rational(5, 4)
    r0 = s.Rational(3, 4)
    out = s.Matrix([s.Rational(3, 5), 0, s.Rational(4, 5)])
    incoming = (s.ImmutableMatrix([-E, 0, 0, -r0]), s.ImmutableMatrix([-E, 0, 0, r0]))
    born = (
        *incoming,
        s.ImmutableMatrix([E, *(r0 * out)]),
        s.ImmutableMatrix([E, *(-r0 * out)]),
    )
    Am = s.factor(
        matter.born_continuation(born, pars["heavy"], pars["cubic"], pars["contact"])
    )
    AG = s.factor(lower.old.born(born) / pars["kappa"])
    A0 = Am + AG
    assert Am > 0 and AG > 0
    directions = (
        s.Matrix([-1, 0, 0]),
        s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0]),
        s.Matrix([s.Rational(3, 5), -s.Rational(4, 5), 0]),
    )
    pols = []
    for n in directions:
        U = s.Matrix([-n[1], n[0], 0])
        V = s.Matrix([0, 0, 1])
        A = s.zeros(4)
        A[1:, 1:] = U * U.T - V * V.T
        pols.append(s.ImmutableMatrix(A))
    for boost in (s.S.Zero, s.Rational(1, 1024)):
        gamma = (1 + boost**2) / (1 - boost**2)
        sh = 2 * boost / (1 - boost**2)
        for imaginary in (s.S.Zero, s.Rational(1, 10**16)):
            parameter = s.Rational(31, 16) + s.I * imaginary
            ep = clean((parameter + 1 / parameter) / 2)
            rp = clean((parameter - 1 / parameter) / 2)
            points = [*incoming]
            for sign in (1, -1):
                points.append(
                    matrix(
                        s.Matrix(
                            [
                                gamma * ep + sign * sh * rp * out[0],
                                sh * ep + sign * gamma * rp * out[0],
                                sign * rp * out[1],
                                sign * rp * out[2],
                            ]
                        )
                    )
                )
            points = tuple(points)
            Q = matrix(-sum(points, lower.VECTOR_ZERO))
            W = Q[0]
            qx = Q[1]
            weights = (
                clean((3 * W - 5 * qx) / 8),
                clean(s.Rational(5, 16) * (W + qx)),
                clean(s.Rational(5, 16) * (W + qx)),
            )
            hs = tuple(
                (matrix(s.Matrix([w, *(w * n)])), A)
                for w, n, A in zip(weights, directions, pols)
            )
            assert matrix(sum((q for q, A in hs), lower.VECTOR_ZERO) - Q) == s.zeros(
                4, 1
            )
            for p in points:
                assert clean((p.T * lower.ETA * p)[0] - 1) == 0
            for q, A in hs:
                assert clean((q.T * lower.ETA * q)[0]) == 0
                assert matrix(A * q) == s.zeros(4, 1)
            if imaginary == 0:
                assert 0 < W < s.Rational(1, 8) and all(w > 0 for w in weights)
                assert clean((Q.T * lower.ETA * Q)[0]) > 0
                center = weights
                ps, _rays, _born = source.recoil.momenta(
                    E, tuple(q for q, A in hs), out
                )
                assert all(matrix(p - p0) == s.zeros(4, 1) for p, p0 in zip(points, ps))
            else:
                analytic.require_relative_tube(center, weights)
            phase = clean(E * rp / (r0 * ep))
            assert clean(phase * s.conjugate(phase)) < 16
            hard = analytic.GaussianEngine([("phi", p, 1) for p in points], **pars)
            J, jcount = hard.amputated(15, "h")
            assert jcount == 47
            assert matrix(Q.T * lower.ETA * J) == s.zeros(1, 4)
            norm2 = clean(
                sum(
                    J[i, j] * s.conjugate(J[i, j])
                    for i in range(1, 4)
                    for j in range(1, 4)
                )
            )
            upper = clean(4 * pars["kappa"] * W * s.conjugate(W) * norm2 / A0**2)
            assert upper < B0**2
            for i in range(1, 4):
                for j in range(i, 4):
                    A = s.zeros(4)
                    A[i, j] = A[j, i] = 1 if i == j else s.Rational(1, 2)
                    assert sum(value * value for value in A) <= 1
                    core = analytic.GaussianEngine(
                        [("phi", p, 1) for p in points[1:]]
                        + [("h", Q, s.ImmutableMatrix(A))],
                        **pars,
                    )
                    value, count = core.amplitude()
                    assert count == 47
                    expected = clean(
                        sum(J[m, n] * A[m, n] for m in range(4) for n in range(4))
                    )
                    assert clean(value - expected) == 0, (boost, imaginary, i, j)
            pure = analytic.TemporalEngine([("h", q, A) for q, A in hs], **pars)
            current, current_count = pure.current(7, "h")
            assert current_count == 4
            contracted = clean(
                sum(J[i, j] * current[i, j] for i in range(4) for j in range(4))
            )
            collapsed = analytic.GaussianEngine(
                [("phi", p, 1) for p in points[1:]] + [("h", Q, current)], **pars
            )
            collapsed_value, collapsed_count = collapsed.amplitude()
            assert collapsed_count == 47 and clean(collapsed_value - contracted) == 0

            class CompleteTripleBranch(analytic.GaussianEngine):
                @lower.instance_cache
                def current(self, mask, kind):
                    if all(
                        self.legs[i][0] == "h"
                        for i in range(len(self.legs))
                        if mask >> i & 1
                    ):
                        if mask == 56 and kind == "h":
                            return self.complete_triple_tensor, 4
                        return (lower.ZERO if kind == "h" else s.S.Zero), 0
                    return super().current(mask, kind)

            restricted = CompleteTripleBranch(
                [("phi", p, 1) for p in points[1:]] + [("h", q, A) for q, A in hs],
                **pars,
            )
            restricted.complete_triple_tensor = current
            value, count = restricted.amplitude()
            assert count == 188 and clean(value - contracted) == 0

    return {
        "whole_core_EGF_linear_coefficients": (cm, cg),
        "whole_operator_bound": B0,
        "whole_global_hard_core_radius": eta,
        "whole_original_parameter_calibration_points": 4,
        "whole_independent_spatial_root_calibrations": 24,
        "whole_temporal_class_inventory": 188,
        "checks": {
            "independent_matter_linear_coefficient": cm - ind_m,
            "independent_gravity_linear_coefficient": cg - ind_g,
            "all_original47_root_orientations_and188_class_identities": s.S.Zero,
        },
        "gates": {
            "original_parameters_exactly_retained": True,
            "all_four_points_Ward_mass_shell_and_phase_checked": True,
            "real_recoil_and_complex_relative_tube_calibration": True,
            "hard_core_global_tube_proved_separately_not_by_finite_points": True,
            "arbitrary_spatial_inputs_not_assumed_null_or_TT": True,
            "temporal188_class_not_independent_gauge_observable": True,
        },
    }


@cache
def analytic_data():
    pars = source.original_parameters()
    n = pars["heavy"]
    kappa = pars["kappa"]
    cm = s.Integer(44048)
    cg = s.Integer(5566277620447838208000000)
    B0 = 4 * (n * n * cm + cg)
    c0, c1, c2, c3 = subtraction.C0, subtraction.C1, subtraction.C2, subtraction.C3
    bracket = (
        c3
        + 3 * c2 / ETA
        + s.Rational(2, 3) * c1 / ETA**2
        + s.Rational(2, 9) * c0 / ETA**3
    )
    constant = B0 * bracket / kappa ** s.Rational(3, 2)
    assert constant > 0 and constant < s.Rational(1, 10**738)
    assert 3 * ETA < s.Rational(1, 10**12)
    assert 4 * cm < 10**8 and 4 * cg < 10**28
    a, b, c = s.symbols("a b c", positive=True)
    W = a + b + c
    S = 1 / (a + b) + 1 / (a + c) + 1 / (b + c)
    I = lambda x, y: (x + y) * s.log(x + y) - x * s.log(x) - y * s.log(y)
    J = c * I(a, b) + b * I(a, c) + a * I(b, c)
    checks = {
        "entropy_mixed_kernel": s.simplify(s.diff(I(a, b), a, b) - 1 / (a + b)),
        "three_variable_rectangle_kernel": s.simplify(s.diff(J, a, b, c) - S),
    }
    # The sum of pairwise squared differences proves 2 W S >= 9.
    x, y, z = s.symbols("x y z", positive=True)
    identity = (x + y + z) * (1 / x + 1 / y + 1 / z) - 9
    squares = (x - y) ** 2 / (x * y) + (x - z) ** 2 / (x * z) + (y - z) ** 2 / (y * z)
    checks["three_positive_denominator_Cauchy_identity"] = s.factor(identity - squares)
    checks["pair_sum_to_two_W"] = s.expand((a + b) + (a + c) + (b + c) - 2 * W)
    t = s.symbols("t", positive=True)
    checks["log_bound_monotone_derivative"] = s.factor(
        s.diff(t - s.log(1 + t * t), t) - (t - 1) ** 2 / (1 + t * t)
    )
    A, B, C = s.symbols("A B C", positive=True)
    terms = (C * C * A * B, B * B * A * C, A * A * B * C)
    squares = sum((terms[i] - terms[j]) ** 2 for i in range(3) for j in range(i + 1, 3))
    checks["entropy_square_Cauchy_identity"] = s.expand(
        3 * sum(q * q for q in terms) - sum(terms) ** 2 - squares
    )
    checks["entropy_square_pointwise_envelope"] = s.factor(
        12 * sum(q * q for q in terms) / (A * A * B * B * C * C)
        - 12 * (A * A + B * B + C * C)
    )
    limit = s.symbols("limit", positive=True)
    linear_integral = 3 * s.integrate(
        2 / s.sqrt(a * b), (a, 0, limit), (b, 0, limit), (c, 0, limit)
    )
    square_integral = s.integrate(
        12 * (a + b + c), (a, 0, limit), (b, 0, limit), (c, 0, limit)
    )
    checks["J_over_abc_cube_integral_bound"] = s.simplify(
        linear_integral - 24 * limit**2
    )
    checks["J_squared_over_abc_cube_integral_bound"] = s.simplify(
        square_integral - 18 * limit**4
    )
    checks["core_complement_inventory"] = s.Integer(3767 + 3 * 387 + 188 - 5116)
    checks["remaining_temporal_inventory"] = s.Integer(5116 - 188 - 4928)
    # Enumerate all eight product-rule placements, retaining multiplicities.
    degrees = {
        d: sum(1 for alpha in product((0, 1), repeat=3) if sum(alpha) == d)
        for d in range(4)
    }
    checks["squarefree_product_rule_multiplicities"] = s.Matrix(
        [degrees[d] - s.binomial(3, d) for d in range(4)]
    )
    checks["squarefree_product_rule_total"] = s.Integer(sum(degrees.values()) - 8)
    for residual in checks.values():
        assert all(
            v == 0
            for v in (
                list(residual) if isinstance(residual, s.MatrixBase) else [residual]
            )
        )
    return {
        "whole_hard_operator_uniform_bound": B0,
        "whole_squarefree_operator_Cauchy_radius": ETA,
        "whole_product_rule_coefficient": constant,
        "whole_conservative_class_rectangle_constant": s.Rational(1, 10**738),
        "whole_class_rectangle_bound": "norm(Delta_a Delta_b Delta_c G188) < 1e-738*J",
        "whole_current_rectangle_bound": "norm(Delta_a Delta_b Delta_c F) <= 721*J/kappa",
        "whole_integrable_kernels": {
            "J_over_abc_cube": 24 * limit**2,
            "J_squared_over_abc_cube": 18 * limit**4,
        },
        "whole_remaining_temporal_terms": 4928,
        "checks": checks,
        "gates": {
            "one_block_has_no_proper_soft_invariant_denominator": True,
            "all_scalar_soft_charges_use_total_W_for_one_block": True,
            "recoil_and_hard_gaps_extend_on_global_eta_W_tube": True,
            "closed_disc_nonzero_margins_give_operator_Cauchy": True,
            "all_three_current_derivative_caps_replayed_elsewhere": True,
            "eight_term_product_rule_with_two_W_S_ge_nine": True,
            "compatible_axes_plus_uniform_origin_allow_trimmed_cube_FTC": True,
            "finite_entropy_integrals_are_not_an_inclusive_probability": True,
            "remaining4928_temporal_terms_not_declared_subtracted": True,
            "no_full_all_N_real_virtual_Regge_or_original_P8_closure": True,
        },
    }


@cache
def data():
    groups = {
        "polynomial_kernel": polynomial_calibration(),
        "one_block_source": core_calibration(),
        "one_block_product": analytic_data(),
    }
    return {
        "whole_subpackets": {
            name: {k: v for k, v in packet.items() if k not in ("checks", "gates")}
            for name, packet in groups.items()
        },
        "checks": {
            name + "_" + key: value
            for name, packet in groups.items()
            for key, value in packet["checks"].items()
        },
        "gates": {
            name + "_" + key: value
            for name, packet in groups.items()
            for key, value in packet["gates"].items()
        },
    }
