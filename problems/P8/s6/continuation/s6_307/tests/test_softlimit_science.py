"""Independent diagram regrouping, joint limits and uniform integration checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_uniform_radiation_soft_limit import (
    audit,
    currents,
    estimates,
    softlimit,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_all_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize("row", (0, 1, 2))
@pytest.mark.parametrize("pol", (0, 1))
@pytest.mark.parametrize("partition", (0, 1, 2))
def test_original_external_diagrams_and_mandatory_offshell_terms(row, pol, partition):
    v = source.vertices
    ps, k, _ = v.sample(row)
    eps = currents.polarizations(row)[pol]
    left, right = v.PARTS[partition]
    base, shift = currents.grouped_channel(ps, k, eps, left, right)
    # Independent source is the four original scalar-propagator emission diagrams.
    original = v.channel(ps, k, eps, left, right)[0]
    assert s.factor(base + shift - original) == 0
    assert shift != 0 and s.factor(base - original) != 0
    for p in ps:
        assert s.factor(currents.current(-p, k, eps) + currents.current(p, k, eps)) == 0


@pytest.mark.parametrize("row", (0, 1, 2))
@pytest.mark.parametrize("pol", (0, 1))
def test_mixed_pair_current_bounds_on_complete_recoil_state(row, pol):
    v = source.vertices
    ps, k, p0 = v.sample(row)
    eps = currents.polarizations(row)[pol]
    w = k[0]
    for a, b in ((0, 2), (1, 3), (0, 3), (1, 2)):
        tau = -v.dot(p0[a] + p0[b], p0[a] + p0[b])
        paired = currents.current(ps[a], k, eps) + currents.current(ps[b], k, eps)
        assert abs(paired) < 600 * s.sqrt(tau + w * w) / w
        relative = s.sqrt(sum((ps[a][j] + ps[b][j]) ** 2 for j in range(1, 4)))
        assert abs(paired) <= 144 * relative / w
        assert abs(currents.current(ps[a], k, eps)) <= 16 / w
        assert abs(currents.current(ps[b], k, eps)) <= 16 / w


@pytest.mark.parametrize("mode", ("low", "high"))
@pytest.mark.parametrize("den", (100, 10000, 1000000))
@pytest.mark.parametrize("pol", (0, 1))
def test_joint_forward_and_soft_limit_of_original_complete_gravity_tree(mode, den, pol):
    v = source.vertices
    E = s.Rational(5, 4)
    h = s.Rational(1, den)
    c = s.Rational(1, 1000) if mode == "low" else s.S.One
    aux = 2 - c * h
    pair_energy = (aux + 1 / aux) / 2
    w = E - pair_energy**2 / E
    direction = (2 * h / (1 + h * h), 0, (1 - h * h) / (1 + h * h))
    ps, k, p0 = v.recoil(E, w, (s.Rational(3, 5), s.Rational(4, 5), 0), direction)
    delta = min(
        s.S.One,
        -v.dot(p0[0] + p0[2], p0[0] + p0[2]),
        -v.dot(p0[0] + p0[3], p0[0] + p0[3]),
    )
    eps = currents.polarizations(0)[pol]
    whole = v.amplitude(ps, k, eps)
    born = v.born(p0)
    common_soft = sum(currents.current(p, k, eps) for p in p0)
    remainder = s.factor(whole / born - common_soft)
    assert abs(whole / born) * w * s.sqrt(delta + w * w) / delta < estimates.DIRECT
    if mode == "low":
        assert w < s.sqrt(delta) / 192
        assert abs(remainder) < estimates.LOW
    else:
        assert w > s.sqrt(delta) / 192
    assert born > 8 / delta


@pytest.mark.parametrize("x", ("0.125", "0.001", "1e-20"))
@pytest.mark.parametrize("ratio", ("0.001", "0.1", "0.5", "0.9"))
def test_entropy_maximum_and_high_band_integral_in_log_energy(x, ratio):
    with mp.workdps(60):
        x, q = mp.mpf(x), mp.mpf(ratio)
        y = 192 * x * q
        delta = y * y
        a = y / 192
        # This independent inequality even holds outside delta<=1.
        exact = delta * mp.log(x / a)
        assert 0 < exact < (192 * x) ** 2 / 2
        integral = mp.quad(lambda z: delta, [mp.log(a), mp.log(x)])
        assert abs(integral - exact) < mp.mpf("1e-50") * max(mp.mpf(1), abs(exact))
        assert abs(
            (192 * x / mp.sqrt(mp.e)) ** 2 / 2 - (192 * x) ** 2 / (2 * mp.e)
        ) < mp.mpf("1e-50")
        D = mp.mpf(10) ** 17
        actual = mp.quad(
            lambda z: 2 * D * D * delta * delta / (delta + mp.exp(2 * z)),
            [mp.log(a), mp.log(x)],
        )
        assert 0 < actual <= 2 * D * D * exact


@pytest.mark.parametrize(
    "resolution", (s.Rational(1, 8), s.Rational(1, 1000), s.Rational(1, 10**30))
)
@pytest.mark.parametrize(
    "energy,transfer",
    (
        (s.Rational(25, 4), -s.Rational(9, 8)),
        (9, -1),
        (16, -6),
        (9, -s.Rational(1, 10**400)),
        (9, -5 + s.Rational(1, 10**400)),
    ),
)
def test_full_angle_domain_and_both_uniform_rate_bounds(energy, transfer, resolution):
    bound = softlimit.finite_remainder_bound(energy, transfer, resolution, source.KAPPA)
    assert 0 < bound <= s.Integer(10) ** 32 / source.KAPPA
    assert (
        bound <= (2 * 10**14 * resolution + 2 * 10**37 * resolution**2) / source.KAPPA
    )
    assert softlimit.original_relative_bound(
        energy, transfer, resolution
    ) == s.Rational(4, 10**768)
    vanishing = softlimit.original_vanishing_relative_bound(
        energy, transfer, resolution
    )
    assert vanishing > 0
    assert (
        estimates.global_direct_gravity_bound(
            energy, transfer, resolution, source.KAPPA
        )
        > 0
    )


@pytest.mark.parametrize("den", (1000, 10**40, 10**400))
def test_arbitrarily_small_adaptive_low_region(den):
    tau = s.Rational(1, den * den)
    w = s.Rational(1, 192 * den)
    assert estimates.low_remainder_bound(
        9, -tau, w, source.KAPPA
    ) == estimates.LOW / s.sqrt(source.KAPPA)
    with pytest.raises(ValueError, match="low-energy"):
        estimates.low_remainder_bound(9, -tau, 2 * w, source.KAPPA)


@pytest.mark.parametrize(
    "bad", (True, False, 1.0, s.Float(1), "1", None, s.oo, s.nan, s.I, s.Symbol("x"))
)
def test_nonexact_soft_reference_coordinates_rejected(bad):
    with pytest.raises((ValueError, TypeError)):
        softlimit.reference_factor(bad, 0, s.Rational(1, 8))
    with pytest.raises((ValueError, TypeError)):
        estimates.require_domain(9, -1, bad)


@pytest.mark.parametrize("resolution", (0, -1, s.Rational(1, 7), 1))
def test_outside_detector_domain_rejected(resolution):
    with pytest.raises(ValueError):
        softlimit.finite_remainder_bound(9, -1, resolution, source.KAPPA)


@pytest.mark.parametrize(
    "energy,transfer", ((4, -1), (17, -1), (9, 0), (9, -5), (9, 1))
)
def test_forward_or_nonphysical_hard_coordinates_rejected(energy, transfer):
    with pytest.raises(ValueError):
        softlimit.original_relative_bound(energy, transfer, s.Rational(1, 8))


@pytest.mark.parametrize(
    "index,conversion",
    ((-1, 0), (7 / source.KAPPA, 0), (0, 113 / source.KAPPA), (0, -113 / source.KAPPA)),
)
def test_unsupported_original_reference_parameters_rejected(index, conversion):
    with pytest.raises(ValueError):
        softlimit.reference_factor(index, conversion, s.Rational(1, 8))


def test_zero_index_exact_reference_and_unexpanded_resolution():
    x = s.Rational(1, 10**1000)
    assert softlimit.reference_factor(0, 0, x) == 1
    a = 3 / source.KAPPA
    actual = softlimit.reference_factor(a, 0, x)
    assert actual == source.gamma.gamma_factor(a) * x**a
    assert actual != source.gamma.gamma_factor(a) * (1 + a * s.log(x))


@pytest.mark.parametrize("chi", (1, 10, 10000))
def test_exponentially_small_detector_comparison_without_power_expansion(chi):
    with mp.workdps(60):
        kap = mp.mpf(10) ** 800
        alpha = mp.mpf(3) / kap
        logx = -kap * chi
        # Log-domain evaluation avoids losing the meaningful alpha*logx product.
        first = mp.log(mp.mpf(2) * 10**14) + logx
        second = mp.log(mp.mpf(2) * 10**37) + 2 * logx
        lognumerator = first + mp.log1p(mp.exp(second - first))
        log_ratio_upper = mp.log(2) + lognumerator - mp.log(kap) - alpha * logx
        assert abs(alpha * logx + 3 * chi) < mp.mpf("1e-50")
        assert log_ratio_upper < mp.log(4) + 32 * mp.log(10) - mp.log(kap)
        assert 1 + alpha * logx < 0
        assert mp.exp(alpha * logx) > 0


def test_original_and_scoped_frontiers_and_no_allN_promotion():
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 163
    assert len(audit.qualifications()) == 6
    assert audit.matching()[:-1] == audit.previous.matching()
    assert all(
        row == 0
        for value in ROWS.values()
        for row in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert "All-N nonleading radiation" in audit.observable()["not_established"]
    assert (
        "single_real_not_promoted_to_allN_nonleading_control"
        in softlimit.data()["gates"]
    )
