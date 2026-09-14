"""Complete local tensor/projective BRST identities before quantum promotion."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_complement_measure import split as complement
from p8_vacuum_affine_spatial_gauge import gauge as previous_gauge

from .exterior import G, Jets


@cache
def tensor_context():
    jet = Jets()
    fields = {}
    fields["scalar"] = jet.add_tensor("f", 0, ())
    fields["density"] = jet.add_tensor("rho", 0, (), 1)
    fields["covector"] = jet.add_tensor("W", 1, (1,))
    fields["vector_density"] = jet.add_tensor("pW", 1, (-1,), 1)
    fields["metric"] = jet.add_tensor("g", 2, (1, 1), symmetric=True)
    fields["metric_momentum"] = jet.add_tensor("p", 2, (-1, -1), 1, symmetric=True)
    fields["Q"] = jet.add_tensor("Q", 2, (-1, -1), s.Rational(2, 3), symmetric=True)
    fields["shift"] = jet.add_tensor("shift", 1, (-1,))
    for (i,), value in fields["shift"].items():
        jet.rules["shift" + str(i)] += jet.ghost(i, (0,))
    fields["connection"] = jet.add_affine()
    for family, components in fields.items():
        for idx, value in components.items():
            residual = jet.brst(jet.brst(value))
            assert not residual, (family, idx, residual.data())
    for i in range(1, 4):
        assert not jet.brst(jet.sg(i))
    return jet, fields


@cache
def projective():
    jet = Jets()
    full = jet.add_affine()
    trace = {mu: sum((full[a, a, mu] for a in range(4)), G()) for mu in range(4)}
    for mu in range(4):
        wanted = (
            jet.transport(trace[mu])
            + sum((jet.ghost(j, (mu,)) * trace[j] for j in range(1, 4)), G())
            + jet.d(jet.divergence_ghost(), mu)
        )
        assert not jet.brst(trace[mu]) - wanted
    jet.add_affine(compensated=True)
    for idx, value in full.items():
        assert not jet.brst(jet.brst(value)), idx
    for mu in range(4):
        wanted = jet.transport(trace[mu]) + sum(
            (jet.ghost(j, (mu,)) * trace[j] for j in range(1, 4)), G()
        )
        assert not jet.brst(trace[mu]) - wanted
        eta = -jet.d(jet.divergence_ghost(), mu) / 4
        wanted_eta = jet.transport(eta) - sum(
            (
                jet.ghost(j, (mu,)) * jet.d(jet.divergence_ghost(), j) / 4
                for j in range(1, 4)
            ),
            G(),
        )
        assert not jet.brst(eta) - wanted_eta
    return jet, full


@cache
def point_maps():
    jet = Jets()
    gamma = jet.add_tensor("gamma", 2, (1, 1), symmetric=True)
    shift = jet.add_tensor("n", 1, (-1,))
    for (i,) in shift:
        jet.rules["n" + str(i)] += jet.ghost(i, (0,))
    lapse = jet.add_tensor("N", 0, ())[()]
    normal = jet.add_tensor("T", 0, ())[()]
    conformal = jet.add_tensor("C", 0, ())[()]
    W = jet.add_tensor("W", 1, (1,))
    metric = {}
    metric[0, 0] = lapse * lapse - sum(
        (
            conformal * gamma[i, j] * shift[i,] * shift[j,]
            for i in range(1, 4)
            for j in range(1, 4)
        ),
        G(),
    )
    for i in range(1, 4):
        metric[0, i] = metric[i, 0] = -sum(
            (conformal * gamma[i, j] * shift[j,] for j in range(1, 4)), G()
        )
        for j in range(1, 4):
            metric[i, j] = -conformal * gamma[i, j]
    vector = {0: lapse * normal + sum((shift[i,] * W[i,] for i in range(1, 4)), G())}
    vector.update({i: W[i,] for i in range(1, 4)})
    for (a, b), value in metric.items():
        wanted = jet.transport(value) + sum(
            (
                jet.ghost(i, (a,)) * metric[i, b] + jet.ghost(i, (b,)) * metric[a, i]
                for i in range(1, 4)
            ),
            G(),
        )
        assert not jet.brst(value) - wanted, (a, b)
    for a, value in vector.items():
        wanted = jet.transport(value) + sum(
            (jet.ghost(i, (a,)) * vector[i] for i in range(1, 4)), G()
        )
        assert not jet.brst(value) - wanted, a
    eta = s.diag(1, -1, -1, -1)
    linear_metric = {
        (a, b): sum(
            (
                eta[i, b] * jet.ghost(i, (a,)) + eta[a, i] * jet.ghost(i, (b,))
                for i in range(1, 4)
            ),
            G(),
        )
        for a, b in product(range(4), repeat=2)
    }
    for a, b, c in product(range(4), repeat=3):
        variation = sum(
            (
                eta[a, d]
                * (
                    jet.d(linear_metric[d, c], b)
                    + jet.d(linear_metric[d, b], c)
                    - jet.d(linear_metric[b, c], d)
                )
                / 2
                for d in range(4)
            ),
            G(),
        )
        assert not variation - jet.ghost(a, (b, c))
    return jet, metric, vector


@cache
def canonical_action():
    jet = Jets()
    gamma = jet.add_tensor("g", 2, (1, 1), symmetric=True)
    pi = jet.add_tensor("p", 2, (-1, -1), 1, symmetric=True)
    vector = jet.add_tensor("W", 1, (1,))
    pvector = jet.add_tensor("pW", 1, (-1,), 1)
    scalar = [jet.add_tensor("f" + str(a), 0, ())[()] for a in range(4)]
    pscalar = [jet.add_tensor("pf" + str(a), 0, (), 1)[()] for a in range(4)]
    shift = jet.add_tensor("n", 1, (-1,))
    for (i,) in shift:
        jet.rules["n" + str(i)] += jet.ghost(i, (0,))
    h0 = jet.add_tensor("H0", 0, (), 1)[()]
    spatial = {}
    for i in range(1, 4):
        metric = sum(
            (
                pi[j, k] * jet.d(gamma[j, k], i) - 2 * jet.d(gamma[i, j] * pi[j, k], k)
                for j, k in product(range(1, 4), repeat=2)
            ),
            G(),
        )
        matter = sum((pscalar[a] * jet.d(scalar[a], i) for a in range(4)), G())
        vectorpart = sum(
            (
                pvector[j,] * (jet.d(vector[j,], i) - jet.d(vector[i,], j))
                - vector[i,] * jet.d(pvector[j,], j)
                for j in range(1, 4)
            ),
            G(),
        )
        spatial[i] = metric + matter + vectorpart
    theta = sum(
        (pi[i, j] * jet.d(gamma[i, j], 0) for i, j in product(range(1, 4), repeat=2)),
        G(),
    )
    theta += sum((pvector[i,] * jet.d(vector[i,], 0) for i in range(1, 4)), G())
    theta += sum((pscalar[a] * jet.d(scalar[a], 0) for a in range(4)), G())
    action = theta - h0 - sum((shift[i,] * spatial[i] for i in range(1, 4)), G())
    flux = {}
    for j in range(1, 4):
        flux[j] = sum(
            (
                2 * pi[j, k] * gamma[k, i] * jet.ghost(i, (0,))
                for k, i in product(range(1, 4), repeat=2)
            ),
            G(),
        )
        flux[j] += sum(
            (pvector[j,] * vector[i,] * jet.ghost(i, (0,)) for i in range(1, 4)), G()
        )
    residual = jet.brst(action) - sum(
        (jet.d(jet.ghost(i) * action + flux[i], i) for i in range(1, 4)), G()
    )
    assert not residual, residual.data()
    for i in range(1, 4):
        wanted = (
            jet.transport(spatial[i])
            + jet.divergence_ghost() * spatial[i]
            + sum((jet.ghost(j, (i,)) * spatial[j] for j in range(1, 4)), G())
        )
        assert not jet.brst(spatial[i]) - wanted
    Q = jet.add_tensor("Q", 2, (-1, -1), s.Rational(2, 3), symmetric=True)
    chi = {i: sum((jet.d(Q[i, j], j) for j in range(1, 4)), G()) for i in range(1, 4)}
    for i in range(1, 4):
        full = (
            jet.transport(chi[i])
            - sum((jet.ghost(i, (j,)) * chi[j] for j in range(1, 4)), G())
            + s.Rational(2, 3) * jet.divergence_ghost() * chi[i]
        )
        full -= sum(
            (Q[j, k] * jet.ghost(i, (j, k)) for j, k in product(range(1, 4), repeat=2)),
            G(),
        )
        full -= sum(
            (Q[i, j] * jet.d(jet.divergence_ghost(), j) / 3 for j in range(1, 4)), G()
        )
        assert not jet.brst(chi[i]) - full
        assert not jet.brst(full)
    return jet, action, spatial, flux, chi


@cache
def full_projective_doublet():
    jet = Jets()
    full = jet.add_affine()
    for mu in range(4):
        eta = jet.ghost(30 + mu)
        jet.odd_rules[30 + mu] = jet.transport(eta) + sum(
            (jet.ghost(j, (mu,)) * jet.ghost(30 + j) for j in range(1, 4)), G()
        )
    for a, b, mu in full:
        if a == b:
            jet.rules["Gamma" + str(a) + str(b) + str(mu)] += jet.ghost(30 + mu)
    for value in full.values():
        assert not jet.brst(jet.brst(value))
    for mu in range(4):
        g = sum((full[a, a, mu] for a in range(4)), G()) / 4
        shifted_ghost = jet.brst(g)
        assert not jet.brst(shifted_ghost)
        assert not jet.brst(jet.sg(30 + mu))
    return jet, full


@cache
def gauge_fermion():
    jet, _, _, _, chi = canonical_action()
    alpha = s.Symbol("gauge_width", real=True)
    psi = G()
    wanted = G()
    for i in range(1, 4):
        b = jet.field("B" + str(i))
        jet.rules["B" + str(i)] = G()
        bar = jet.ghost(10 + i)
        jet.odd_rules[10 + i] = b
        psi += bar * (chi[i] + alpha * b / 2)
        wanted += b * chi[i] + alpha * b * b / 2 - bar * jet.brst(chi[i])
    assert not jet.brst(psi) - wanted
    assert not jet.brst(wanted)
    return jet, psi, wanted


def scalar(value):
    """Every exterior coefficient has its own independent commuting label."""
    return (
        s.expand(
            sum(
                coef
                * s.Symbol(
                    "wedge_"
                    + "_".join(str(i) + "d" + "".join(map(str, idx)) for i, idx in mon)
                )
                for mon, coef in value.terms.items()
            )
        )
        if value
        else s.S.Zero
    )


def column(values):
    return s.Matrix([scalar(value) for value in values])


@cache
def tensor_identities():
    jet, families = tensor_context()
    checks = {}
    for family, fields in families.items():
        unique = list(
            dict.fromkeys(next(iter(value.terms.values())) for value in fields.values())
        )
        checks["whole_" + family + "_nilpotency"] = column(
            jet.brst(jet.brst(value)) for value in unique
        )
    checks["whole_three_ghost_nilpotency"] = column(
        jet.brst(jet.sg(i)) for i in range(1, 4)
    )
    return {
        "convention": "Left odd derivation, s tensor=Lie_c tensor, sc^i=c^j partial_j c^i. c^0=0; time derivatives of c^i are not zero.",
        "whole_bosonic_transformations": {
            name: value.data() for name, value in jet.rules.items()
        },
        "whole_three_ghost_transformations": {
            str(i): jet.sg(i).data() for i in range(1, 4)
        },
        "checks": checks,
        "gates": {
            "genuine_anticommuting_ghost_jets": True,
            "all_density_weights_and_time_derivatives_retained": True,
            "geometric_rules_not_off_constraint_Dirac_BFV_extension": True,
        },
    }


@cache
def projective_identities():
    jet, full = full_projective_doublet()
    trace = {mu: sum((full[a, a, mu] for a in range(4)), G()) for mu in range(4)}
    compensated, connection = projective()
    checks = {
        "whole_64_semidirect_connection_nilpotency": column(
            jet.brst(jet.brst(value)) for value in full.values()
        ),
        "whole_four_projective_ghost_nilpotency": column(
            jet.brst(jet.sg(30 + mu)) for mu in range(4)
        ),
        "whole_four_trace_and_shifted_ghost_doublets": column(
            jet.brst(jet.brst(trace[mu] / 4)) for mu in range(4)
        ),
        "whole_64_compensated_connection_nilpotency": column(
            compensated.brst(compensated.brst(value)) for value in connection.values()
        ),
    }
    for mu in range(4):
        expected = (
            jet.transport(trace[mu])
            + sum((jet.ghost(j, (mu,)) * trace[j] for j in range(1, 4)), G())
            + jet.d(jet.divergence_ghost(), mu)
            + 4 * jet.ghost(30 + mu)
        )
        checks["whole_trace_inhomogeneous_variation_" + str(mu)] = scalar(
            jet.brst(trace[mu]) - expected
        )
    return {
        "whole_semidirect_connection_rules": {
            name: value.data() for name, value in jet.rules.items()
        },
        "whole_compensated_connection_rules": {
            name: value.data() for name, value in compensated.rules.items()
        },
        "whole_four_shifted_projective_ghosts": {
            str(mu): jet.brst(trace[mu] / 4).data() for mu in range(4)
        },
        "trace_and_compensation": "sF_mu=Lie_c F_mu+partial_mu(div c)+4 eta_mu. F=0 requires eta_mu=-partial_mu(div c)/4, including mu=0. For g=F/4, zeta=sg=eta+Lie_c g+d(div c)/4, and szeta=0.",
        "checks": checks,
        "gates": {
            "uncompensated_spatial_transformation_does_not_preserve_trace_gauge": True,
            "all_four_projective_and_three_spatial_ghosts_retained": True,
            "finite_triangular_ghost_relabeling_not_continuum_regulator": True,
        },
    }


@cache
def physical_point_identities():
    jet, metric, vector = point_maps()
    mr = []
    vr = []
    for (a, b), value in metric.items():
        wanted = jet.transport(value) + sum(
            (
                jet.ghost(i, (a,)) * metric[i, b] + jet.ghost(i, (b,)) * metric[a, i]
                for i in range(1, 4)
            ),
            G(),
        )
        mr.append(jet.brst(value) - wanted)
    for a, value in vector.items():
        vr.append(
            jet.brst(value)
            - jet.transport(value)
            - sum((jet.ghost(i, (a,)) * vector[i] for i in range(1, 4)), G())
        )
    return {
        "whole_physical_metric": {
            str(key): value.data() for key, value in metric.items()
        },
        "whole_physical_vector": {
            str(key): value.data() for key, value in vector.items()
        },
        "dictionary": "h=C gamma; W_0=N T+N^i W_i. C and every original coefficient remain their fixed scalar functions of u,N. This dictionary does not reset the clock, source or state.",
        "checks": {
            "whole_sixteen_ADM_metric_entries": column(mr),
            "whole_four_physical_vector_entries": column(vr),
        },
        "gates": {
            "entire_time_dependent_shift_contact_retained": True,
            "flat_64_connection_inhomogeneous_variation_independently_checked": True,
            "no_coefficient_or_source_change": True,
        },
    }


@cache
def canonical_identities():
    jet, action, spatial, flux, chi = canonical_action()
    _, psi, gauge_action = gauge_fermion()
    residual = jet.brst(action) - sum(
        (jet.d(jet.ghost(i) * action + flux[i], i) for i in range(1, 4)), G()
    )
    return {
        "whole_canonical_action_density": action.data(),
        "whole_spatial_generators": {
            str(i): value.data() for i, value in spatial.items()
        },
        "whole_time_descriptor_spatial_boundary_flux": {
            str(i): value.data() for i, value in flux.items()
        },
        "whole_gauge_fermion": psi.data(),
        "whole_off_gauge_ghost_and_multiplier_action": gauge_action.data(),
        "physical_phase_boundary": "The first two scalar pairs are M1 and H; the last two are auxiliary N,T primaries before restriction. Restricting their vanishing primary momenta leaves the S257 eleven retained positions and twenty-two phases. The shift is a multiplier here, not an asserted unreduced BFV primary sector. H0 denotes the entire S257 scalar-density Hamiltonian, not a replacement action.",
        "checks": {
            "whole_canonical_action_variation_is_retained_spatial_flux": scalar(
                residual
            ),
            "whole_three_off_gauge_conditions_are_nilpotent": column(
                jet.brst(jet.brst(value)) for value in chi.values()
            ),
            "whole_gauge_fermion_action": scalar(jet.brst(psi) - gauge_action),
            "whole_nonlinear_gauge_action_is_BRST_closed": scalar(
                jet.brst(gauge_action)
            ),
        },
        "gates": {
            "every_vector_Gauss_and_matter_generator_term_retained": True,
            "spatial_time_descriptor_flux_not_deleted_pointwise": any(
                bool(value) for value in flux.values()
            ),
            "full_off_gauge_operator_precedes_surface_restriction": True,
        },
    }


@cache
def projective_ghost_block():
    k = s.Matrix(s.symbols("spatial_wavevector1:4", real=True))
    omega = s.Symbol("time_frequency", real=True)
    kk = k.dot(k)
    mixed = -s.Matrix([omega, *k]) * k.T
    M = kk * s.eye(3) + k * k.T / 3
    Minv = (s.eye(3) - k * k.T / (4 * kk)) / kk
    block = (4 * s.eye(4)).row_join(mixed).col_join(s.zeros(3, 4).row_join(M))
    inverse = (
        (s.eye(4) / 4)
        .row_join(-mixed * Minv / 4)
        .col_join(s.zeros(3, 4).row_join(Minv))
    )
    data = complement.matrices()
    Gp, F = map(s.Matrix, (data["projective"], data["gauge"]))
    J = s.Matrix(s.symbols("original_connection_source0:64", real=True))
    projector = s.eye(64) - Gp * F / 4
    Jp = projector.T * J
    compensation = s.Matrix(s.symbols("four_projective_compensators0:4", real=True))
    return {
        "whole_three_nonzero_spatial_wavevector": k,
        "time_frequency": omega,
        "whole_seven_ghost_gauge_matrix": block,
        "whole_seven_ghost_inverse": inverse,
        "whole_original_64_projectively_invariant_source": Jp,
        "determinant": s.Rational(1024, 3) * kk**3,
        "general_triangular_rule": "Before restriction the seven-generator matrix is [[4I,Lie_xi F+d div xi],[0,M(Q)]]. Its finite-block determinant is 256 det M. The normalized trace F/4 has unit projective factor. Spatial translations remain a residual kernel; the displayed inverse requires k^2>0. No separately regularized continuum determinant product is asserted.",
        "source_rule": "G^T J=0 removes the compensating projective shift from the original source pairing. It does not remove the composite connection source or S259 quadratic inverse contact. All differentiation is performed on that full source functional.",
        "checks": {
            "whole_seven_ghost_left_inverse": (block * inverse - s.eye(7)).applyfunc(
                s.factor
            ),
            "whole_seven_ghost_right_inverse": (inverse * block - s.eye(7)).applyfunc(
                s.factor
            ),
            "whole_seven_ghost_determinant": s.factor(
                block.det(method="domain-ge") - s.Rational(1024, 3) * kk**3
            ),
            "whole_original_four_trace_gauge": F * Gp - 4 * s.eye(4),
            "whole_64_projective_source_invariance": Gp.T * Jp,
            "whole_compensator_original_source_pairing": s.expand(
                (Jp.T * Gp * compensation)[0]
            ),
        },
        "gates": {
            "time_frequency_mixed_gauge_block_not_dropped": mixed[0, :]
            != s.zeros(1, 3),
            "projective_normalization_not_confused_with_spatial_dimension": True,
            "nonzero_momentum_and_residual_translations_required": True,
            "original_nonzero_conditional_source_contact_retained": True,
        },
    }


@cache
def zero_modes():
    x = s.Symbol("periodic_coordinate", real=True)
    e1 = G({((1, ()),): s.S.One})
    e2 = G({((2, ()),): s.S.One})
    c = e1 * s.sin(x) + e2 * s.cos(x)
    dc = e1 * s.cos(x) - e2 * s.sin(x)
    cc = c * dc
    source = s.Function("retained_scalar", real=True)(x)
    defect = -cc * s.diff(source, x)
    K = 2

    def bracket(m, n):
        return (m + n, s.I * (n - m)) if abs(m + n) <= K else (m + n, s.S.Zero)

    def jacobi(m, n, p):
        out = 0
        for a, b, c0 in ((m, n, p), (n, p, m), (p, m, n)):
            mode, coef = bracket(a, b)
            _, last = bracket(mode, c0)
            out += coef * last
        return s.expand(out)

    indices = next(
        (m, n, p)
        for m, n, p in product(range(-K, K + 1), repeat=3)
        if jacobi(m, n, p) != 0
    )
    jac = jacobi(*indices)
    parent = previous_gauge.zero_modes()
    return {
        "whole_mean_zero_ghost_example": c.data(),
        "whole_generated_translation_ghost": cc.data(),
        "whole_projected_scalar_nilpotency_defect": defect.data(),
        "finite_Witt_cutoff": K,
        "finite_Witt_Jacobi_counterexample_modes": indices,
        "finite_Witt_Jacobi_defect": jac,
        "retained_parent_spatial_Leibniz_defect": parent[
            "full_finite_SBP_Leibniz_defect"
        ],
        "full_zero_mode_rule": "Write c=c0+cperp, c0=P0c. sc0=P0(c dot grad c); scperp=(1-P0)(c dot grad c). The first expression is generally nonzero even when c0=0. The translation complement in the local slice is not a closed gauge algebra.",
        "checks": {
            "whole_generated_constant_translation": s.trigsimp(scalar(cc + e1 * e2)),
            "whole_projected_nilpotency_is_translation_action": s.trigsimp(
                scalar(defect - e1 * e2 * s.diff(source, x))
            ),
            "both_gauge_basis_functions_mean_zero": s.Matrix(
                [
                    s.integrate(s.sin(x), (x, 0, 2 * s.pi)),
                    s.integrate(s.cos(x), (x, 0, 2 * s.pi)),
                ]
            ),
        },
        "gates": {
            "projected_scalar_nilpotency_defect_nonzero": bool(defect),
            "finite_Fourier_projection_fails_Jacobi": jac != 0,
            "no_finite_mode_diffeomorphism_regulator_inferred": True,
            "residual_group_volume_and_noninvariant_sources_not_omitted": True,
        },
    }
