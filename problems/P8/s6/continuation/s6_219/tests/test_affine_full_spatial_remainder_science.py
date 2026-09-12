"""Independent full spatial current tests, beyond exact report residuals."""

import mpmath as mp
import pytest
import sympy as s
from _full_spatial_reference import (
    expected_nondecaying,
    literal_pair_hamiltonian,
    original_band_difference,
    pair_covariance,
    symplectic,
)
from p8_vacuum_affine_full_spatial_remainder import (
    assembly,
    audit,
    dimension,
    remainder,
    shapes,
    verify,
    vertices,
)
from p8_vacuum_affine_ordered_scalar_symbol import density, geometry

PACKETS = (
    vertices.feature_data,
    vertices.pair_data,
    shapes.data,
    remainder.data,
    dimension.data,
    assembly.data,
    assembly.ward_data,
)


@pytest.mark.parametrize("packet", PACKETS)
def test_complete_exact_scientific_packet(packet):
    data = packet()
    for key, value in data["checks"].items():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(v) == 0 for v in entries), key
    assert all(bool(v) for v in data["gates"].values())


DIRECTIONS = (
    (s.eye(3), s.eye(3)),
    (s.eye(3), s.diag(-1, -1, 2)),
    (s.diag(-1, -1, 2), s.eye(3)),
    (
        s.Matrix([[1, 2, -1], [2, 3, 1], [-1, 1, 2]]) / 11,
        s.Matrix([[2, -1, 1], [-1, -3, 2], [1, 2, 4]]) / 13,
    ),
)


def numeric(M):
    return mp.matrix(
        [
            [
                mp.mpc(str(s.re(v).evalf(mp.mp.dps)), str(s.im(v).evalf(mp.mp.dps)))
                for v in row
            ]
            for row in M.tolist()
        ]
    )


@pytest.mark.parametrize("D,G", DIRECTIONS)
def test_literal_metric_determinant_fieldstrength_and_second_contact(D, G):
    with mp.workdps(65):
        a, mass = mp.mpf("1.4"), mp.mpf(1000)
        k, ell = mp.matrix([3, -2, 7]), mp.matrix([-5, 4, 11])
        left, right = mp.matrix([1, 3, -2, 5, -7, 4]), mp.matrix([-2, 1, 4, -3, 2, 6])
        Dn, Gn = numeric(D), numeric(G)
        call = lambda e, h: literal_pair_hamiltonian(
            e * Dn + h * Gn, k, ell, left, right, a, mass
        )
        first = mp.diff(lambda e: call(e, 0), 0)
        contact = mp.diff(call, (0, 0), (1, 1))
        ks, ls = s.Matrix([3, -2, 7]), s.Matrix([-5, 4, 11])
        Fk = numeric(vertices.physical_map(ks, s.Rational(7, 5), 1000))
        Fl = numeric(vertices.physical_map(ls, s.Rational(7, 5), 1000))
        current = (left.T * Fk.T * numeric(vertices.current_feature(D)) * Fl * right)[0]
        wanted_contact = (
            left.T * Fk.T * numeric(vertices.contact_feature(D, G)) * Fl * right
        )[0]
        assert abs(current + first) < mp.mpf("1e-48")
        assert abs(contact - wanted_contact) < mp.mpf("1e-48")
        if D == s.eye(3) and G == s.eye(3):
            dropped = numeric(vertices.contact_feature(D, G))
            dropped[6, 6] = 0
            wrong = (left.T * Fk.T * dropped * Fl * right)[0]
            assert abs(wrong - contact) > mp.mpf("1e-5")


@pytest.mark.parametrize("D,G", DIRECTIONS)
def test_full_six_phase_space_pair_against_real_covariance(D, G):
    with mp.workdps(60):
        k, ell = s.Matrix([3, -2, 7]), s.Matrix([-5, 4, 11])
        Fk = vertices.physical_map(k, s.Rational(7, 5), 3)
        Fl = vertices.physical_map(ell, s.Rational(7, 5), 3)
        DM = numeric(Fk.T * vertices.current_feature(D) * Fl)
        GM = numeric(Fk.T * vertices.current_feature(G) * Fl)
        J = symplectic(3)
        # Nonzero real symplectic preparation shears make the source pair
        # genuinely complex. An unsqueezed block-diagonal source cannot
        # distinguish these two phase branches, so it is not a negative control.
        S, T = mp.eye(6), mp.eye(6)
        for i, value in enumerate((".2", "-.3", ".1")):
            S[3 + i, i] = mp.mpf(value)
        for i, value in enumerate(("-.4", ".15", ".25")):
            T[3 + i, i] = mp.mpf(value)
        assert mp.norm(S * J * S.T - J) < mp.mpf("1e-50")
        assert mp.norm(T * J * T.T - J) < mp.mpf("1e-50")
        DM, GM = S.T * DM * T, S.T * GM * T
        A = mp.diag([1, 2, 3, 4, 5, 6])
        B = mp.diag([7, 5, 3, 2, 4, 6])
        U, V = mp.expm(J * A * mp.mpf(".13")), mp.expm(J * B * mp.mpf(".21"))
        covariance, kubo, wrong = pair_covariance(U, V, DM, GM)
        assert abs(covariance - kubo) < mp.mpf("1e-45")
        assert abs(covariance - wrong) > mp.mpf("1e-4")


@pytest.mark.parametrize("channel", density.CHANNELS)
@pytest.mark.parametrize("source_jet", (0, 1, 2))
def test_literal_original_mask_full_log_shell_converges(channel, source_jet):
    with mp.workdps(65):
        time, transfer, mass = mp.mpf(".19"), mp.mpf(31), mp.mpf(1000)
        source = [int(r == source_jet) for r in range(5)]
        errors = []
        for K in (mp.mpf("1e8"), mp.mpf("2e8"), mp.mpf("4e8")):
            actual = original_band_difference(channel, time, transfer, mass, K, source)
            expected = expected_nondecaying(channel, time, transfer, mass, K, source)
            errors.append(actual - expected)
        assert abs(errors[0]) > mp.mpf("1e-12")
        for i in (0, 1):
            assert abs(errors[i + 1] / errors[i] - mp.mpf(".5")) < mp.mpf("1e-4")


@pytest.mark.parametrize("channel", density.CHANNELS)
def test_independent_grazing_strip_is_not_optional(channel):
    with mp.workdps(55):
        args = (
            channel,
            mp.mpf(".17"),
            mp.mpf(37),
            mp.mpf(1000),
            mp.mpf("1e6"),
            [1, 0, 0, 0, 0],
        )
        complete = original_band_difference(*args)
        no_strip = original_band_difference(*args, omit_grazing=True)
        assert abs(complete - no_strip) > 1


def test_ordered_linear_shapes_integrated_formal_transpose_not_equal_kernels():
    t = shapes.t
    D, G = (1 - 4 * t * t) ** 5 * (1 + t), (1 - 4 * t * t) ** 6 * (1 - 2 * t)
    left = shapes.predicted_shapes("trace_scalar")[1]
    right = shapes.predicted_shapes("scalar_trace")[1]
    action = lambda row, f: sum(row[r] * s.diff(f, t, r) for r in range(3))
    subs = {shapes.p: 7, shapes.m: 1000}
    correct = s.cancel((D * action(left, G) - G * action(right, D)).subs(subs))
    wrong = s.cancel((D * action(left, G) - G * action(left, D)).subs(subs))
    interval = (t, -s.Rational(1, 2), s.Rational(1, 2))
    assert s.integrate(correct, interval) == 0
    assert s.integrate(wrong, interval) != 0


@pytest.mark.parametrize("channel", density.CHANNELS)
def test_complex_dimension_geometry_majorants_at_nonreal_parameters(channel):
    bounds, _checks = dimension.geometry_bounds()
    points = (
        (s.Rational(3, 20) + s.I / 5, s.I / 100, s.Rational(4, 5)),
        (-s.I / 4, s.Rational(3, 500) + s.I * s.Rational(4, 500), -s.Rational(2, 3)),
    )
    for offset, y, u in points:
        for key, value in geometry.contractions(channel).items():
            actual = value.subs({geometry.d: 3 + offset, geometry.y: y, geometry.u: u})
            square = s.simplify(actual * s.conjugate(actual))
            assert square <= bounds[channel + "_far_" + key] ** 2


@pytest.mark.parametrize("channel", ("unknown", "", 0, True, None))
def test_scalar_channel_rejection(channel):
    with pytest.raises(ValueError):
        shapes.predicted_centered(channel)
    with pytest.raises(ValueError):
        shapes.predicted_shapes(channel)
    with pytest.raises(ValueError):
        shapes.unit_factor(channel)


def test_complex_dimension_denominator_and_float_rejection():
    with pytest.raises(ValueError, match="floor"):
        dimension.rational_bound(
            1 / (1 - 100 * geometry.y), (geometry.y,), (s.Rational(1, 100),)
        )
    with pytest.raises(ValueError, match="exact rational"):
        dimension.rational_bound(s.Float("1.5"), (), ())


@pytest.mark.parametrize("key", tuple(audit.residuals()))
def test_each_complete_reference_exact_residual(key):
    value = audit.residuals()[key]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(v == 0 for v in entries)


@pytest.mark.parametrize("key", tuple(audit.gates()))
def test_each_complete_reference_proof_gate(key):
    assert audit.gates()[key] is True


@pytest.mark.parametrize("_label,call,args", audit.bad_cases())
def test_each_unsupported_input_or_scope_is_rejected(_label, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_all_packets_exact_serializer_and_nonvacuous_counts():
    import json

    for packet in audit.packets().values():
        body = json.dumps(verify.serialize(verify.payload(packet)))
        assert 100 < len(body) < 100000
    assert len(audit.residuals()) == 73
    assert audit.scalar_entry_count() == 457
    assert len(audit.gates()) == 68
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 171
    assert len(audit.matching()) == 75
    assert len({row["id"] for row in audit.matching()}) == 75


def test_exact_zero_germ_slab_and_auxiliary_band_edges():
    assert audit.require_scope(-s.Rational(1, 2))[0] == -s.Rational(1, 2)
    assert audit.require_scope(s.Rational(1, 2))[0] == s.Rational(1, 2)
    assert audit.require_cutoff(2000) == (2000, 1000)
    assert audit.validate_scope(audit.frontier(), audit.matching())


def test_full_weak_reference_response_scope_is_not_reduced_inverse():
    observable = audit.observable()
    assert "constrained/reduced" in observable["remaining"]
    assert "sourced-parent" in observable["remaining"]
    assert "weak" in observable["canonical_spatial"].lower()
    assert "reference" in observable["full_metric"]
