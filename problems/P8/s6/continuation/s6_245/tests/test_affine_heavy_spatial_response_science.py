"""Independent all-transfer scalar remainder, signs, complete cutoffs and scope checks."""

import math
from itertools import pairwise

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_spatial_response import (
    audit,
    domain,
    endpoints,
    time_remainder,
)
from scipy.integrate import quad


@pytest.mark.parametrize("name,value", tuple(audit.residuals().items()))
def test_every_full_exact_identity(name, value):
    assert (
        all(entry == 0 for entry in value)
        if isinstance(value, s.MatrixBase)
        else value == 0
    ), name


@pytest.mark.parametrize("name,value", tuple(audit.gates().items()))
def test_every_written_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_or_claim_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_complete_counts_and_unchanged_frontier():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
    ) == (97, 126, 66, 9, 292)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 101
    assert audit.matching()[:-1] == audit.previous.matching()


@pytest.mark.parametrize("leg", ("k", "l"))
def test_full_six_iterates_not_first_iterate_or_fitted_UV(leg):
    data = domain.direct_frequency(leg)
    rows = data["all_seven_four_jets"]
    assert len(rows) == 7 and rows[1] != rows[2] and rows[2] == rows[6]
    assert data["full_frozen_coefficient_residual"] == s.zeros(5, 1)


@pytest.mark.parametrize("order", range(7))
def test_all_literal_source_time_rows(order):
    assert time_remainder.literal_rows()[order] == time_remainder.recurrence()[order]
    assert len(time_remainder.recurrence()[order]) == order + 1


def test_entire_error_product_retains_all_four_error_factors():
    expression = time_remainder.data()["complete_four_factor_error_identity"]
    errors = sorted(
        (x for x in expression.free_symbols if str(x).startswith("error_")), key=str
    )
    assert len(errors) == 4 and s.diff(expression, *errors) == 1


def test_full_sixth_bulk_and_fifth_endpoint_are_retained_with_their_own_signs():
    data = time_remainder.data()
    assert data["complete_six_endpoint_coefficients"][-1] == 1
    assert data["whole_sixth_bulk_coefficient"] == -1
    assert data["constants"]["whole_fifth_endpoint_numerator"] > 0
    assert data["constants"]["whole_sixth_bulk_numerator"] > 0


def test_all_four_geometries_and_dimension_floors():
    data = domain.data()
    assert len(data["all_eight_full_dimensional_geometry_majorants"]) == 8
    assert len(data["all_48_exact_coefficient_denominator_floors"]) == 48
    assert all(
        value > 0
        for value in data["all_48_exact_coefficient_denominator_floors"].values()
    )
    assert domain.normalized_geometries()["01"] != domain.normalized_geometries()["10"]
    assert domain.far_geometries()["01"] != domain.far_geometries()["10"]


def test_all_three_complete_momentum_regions_and_every_logarithm():
    data = endpoints.constants()
    rows = data["whole_raw_and_all_fifteen_near_terms"]
    assert len(rows) == 16 and sum(im for q, value, im, label in rows) == 5
    assert len(data["all_five_far_terms"]) == 5
    assert data["complete_unexpanded_low_numerator_over_mass_four"] > 0
    assert all(
        q - int(im) <= 4 and q + 1 - int(im) <= 5 and q + 1 <= 6
        for q, value, im, label in rows
    )
    assert data["entire_endpoint_finite_numerator_over_mass_four"] == s.Rational(
        111268759571258000000000000000, 9
    )
    assert data["entire_endpoint_tail_numerator_over_mass_five"] == s.Rational(
        67387527200968092800000000000000, 27
    )


@pytest.mark.parametrize("q", range(7))
@pytest.mark.parametrize(
    "m,p", ((1, 0), (1, 1), (2, 3), (10**98, 1), (10**98, 10**100))
)
def test_full_mass_and_spatial_weight(q, m, p):
    assert (m + p) ** q <= 8 * m**q * (1 + p * p) ** 3


def test_whole_actual_sum_and_original_tail_are_small_without_a_new_inverse():
    data = endpoints.constants()
    assert data["complete_actual_all_transfer_spatial_response_bound"] < s.Rational(
        1, 10**350
    )
    assert data[
        "whole_original_frequency_projected_response_tail_numerator"
    ] < s.Rational(1, 10**250)
    assert data["full_actual_homogeneous_anchor_normalized"] > 0
    assert data["complete_unchanged_scalar_local_UV_difference_bound"] > 0
    with pytest.raises(ValueError):
        audit.require_stage("same_space_inverse_from_graph_smallness")


@pytest.mark.parametrize(
    "label,degree",
    (
        ("detector_time", 0),
        ("detector_spatial", 0),
        ("source_time", 13),
        ("source_spatial", 6),
    ),
)
def test_complete_graph_without_detector_spatial_loss(label, degree):
    assert audit.require_graph(label, degree) == (label, degree)


def test_all_transfer_magnitudes_not_only_a_bounded_Fourier_ball():
    for p in (0, s.Rational(1, 10**100), 1, 10**99, 10**500):
        assert audit.require_domain(0, p) == (0, p)
    assert audit.require_derivative(1) == 1
    assert audit.require_cutoff_squared(4 * state.MASS2) == 4 * state.MASS2
    assert audit.require_projection(
        "complete_unaveraged_subtracted_two_leg_with_full_local_terms"
    )
    with pytest.raises(ValueError):
        audit.require_derivative(2)


def test_whole_scalar_feature_operator_majorant_for_arbitrary_symmetric_D():
    a, b, c, x, y, z = s.symbols("a b c x y z", real=True)
    D = s.Matrix([[a, x, y], [x, b, z], [y, z, c]])
    trace = s.trace(D)
    norm2 = s.trace(D * D)
    assert (
        s.expand(
            3 * norm2
            - trace**2
            - ((a - b) ** 2 + (a - c) ** 2 + (b - c) ** 2 + 6 * (x * x + y * y + z * z))
        )
        == 0
    )
    assert s.sqrt(3) / 2 + 1 < 3


def numerical_identity():
    x = s.Symbol("time", real=True)
    rate = 24 + x + 3 * x * x
    theta = 24 * x + x * x / 2 + x**3
    amplitude = (x + s.Rational(1, 2)) ** 7 * (
        1 + x / 5 + s.I * (s.Rational(3, 10) - 2 * x / 5 + x * x / 10)
    )
    values = [amplitude]
    for _ in range(6):
        values.append(s.factor(s.diff(values[-1] / rate, x)))
    functions = [s.lambdify(x, row, "mpmath") for row in values]
    theta_fn = s.lambdify(x, theta, "mpmath")
    rate_fn = s.lambdify(x, rate, "mpmath")
    with mp.workdps(70):
        end = mp.mpf("0.2")
        start = mp.mpf("-0.5")
        detector = 1 + end / 2 - mp.j * mp.mpf(".7") * end
        direct = (
            mp.conj(detector)
            * mp.exp(mp.j * theta_fn(end))
            * mp.quad(
                lambda t: mp.exp(-mp.j * theta_fn(t)) * functions[0](t), [start, end]
            )
        )
        endpoints = [
            mp.conj(detector) * mp.j * (-mp.j) ** j * functions[j](end) / rate_fn(end)
            for j in range(6)
        ]
        bulk = (
            -mp.conj(detector)
            * mp.exp(mp.j * theta_fn(end))
            * mp.quad(
                lambda t: mp.exp(-mp.j * theta_fn(t)) * functions[6](t), [start, end]
            )
        )
        assert abs(direct - sum(endpoints) - bulk) < mp.mpf("1e-60")
        assert abs(mp.im(endpoints[5])) > mp.mpf("1e-10")
        assert abs(mp.im(bulk)) > mp.mpf("1e-10")
        assert abs(mp.im(direct - (sum(endpoints) - 2 * endpoints[5] + bulk))) > mp.mpf(
            "1e-9"
        )
        assert abs(mp.im(direct - (sum(endpoints) - bulk))) > mp.mpf("1e-9")
        print(
            "FULL_ORDERED_RETARDED_CURRENT",
            mp.nstr(mp.im(direct), 25),
            "j5",
            mp.nstr(mp.im(endpoints[5]), 25),
            "bulk6",
            mp.nstr(mp.im(bulk), 25),
            flush=True,
        )
        print("ENTIRE_SIX_STEP_NUMERICAL_IDENTITY_C0", flush=True)
        return float(abs(direct - sum(endpoints) - bulk))


def test_independent_whole_retarded_identity_with_nonzero_boundary_and_bulk():
    assert numerical_identity() < 1e-60


def removed_fraction(r, p, cutoff):
    if r > cutoff:
        return 1.0
    if p == 0:
        return 0.0
    if r == 0:
        return float(p > cutoff)
    return max(0.0, min(1.0, (1 + (r * r + p * p - cutoff * cutoff) / (2 * r * p)) / 2))


def integral(power, left, right, p, cutoff):
    points = sorted(
        {
            left,
            right,
            *[x for x in (abs(cutoff - p), cutoff, cutoff + p) if left < x < right],
        }
    )
    result = 0.0
    for a, b in pairwise(points):
        term, _error = quad(
            lambda r: r**power * removed_fraction(r, p, cutoff),
            a,
            b,
            epsabs=1e-10,
            epsrel=2e-12,
            limit=200,
        )
        result += term
    return result


def full_power(power, left, right):
    if power == -1:
        return math.log(right / left)
    return (right ** (power + 1) - left ** (power + 1)) / (power + 1)


def mask_diagnostics():
    count = 0
    largest = 0.0
    # Radius is scaled by the actual mass: the tested estimates are homogeneous.
    for p in (0.0, 0.1, 1.0, 4.0, 20.0, 300.0):
        U = 1 + p
        L = 200 * U
        for K in (2.0, 3.0, 10.0, 100.0, 1000.0, 100000.0):
            cutoff = 4 * math.sqrt(K * K - 1)
            assert cutoff > K
            for q in range(5):
                actual = integral(3 - q, 1, L, p, cutoff)
                whole = full_power(3 - q, 1, L)
                assert actual <= whole * (1 + 1e-10)
                bound = whole * 201 * U / cutoff
                assert actual <= bound * (1 + 1e-10), (p, K, q, actual, bound)
                if bound:
                    largest = max(largest, actual / bound)
                count += 1
            # The far integrand r^-2 is evaluated as a finite integral in x=1/r.
            points = sorted(
                {
                    0.0,
                    1 / L,
                    *[1 / r for r in (abs(cutoff - p), cutoff, cutoff + p) if r > L],
                }
            )
            far = 0.0
            for a, b in pairwise(points):
                far += quad(
                    lambda x, p=p, cutoff=cutoff: (
                        1.0 if x == 0 else removed_fraction(1 / x, p, cutoff)
                    ),
                    a,
                    b,
                    epsabs=1e-13,
                    epsrel=2e-12,
                    limit=200,
                )[0]
            assert far <= min(1 / L, 2 / cutoff) * (1 + 1e-9)
            count += 1
            pts = sorted(
                {
                    0.0,
                    1.0,
                    *[x for x in (abs(cutoff - p), cutoff, cutoff + p) if 0 < x < 1],
                }
            )
            low = 0.0
            for a, b in pairwise(pts):
                low += quad(
                    lambda r, p=p, cutoff=cutoff: (
                        r
                        * r
                        * math.sqrt(1 + r * r / 16)
                        * removed_fraction(r, p, cutoff)
                    ),
                    a,
                    b,
                    epsabs=1e-13,
                    epsrel=2e-12,
                )[0]
            full = quad(
                lambda r: r * r * math.sqrt(1 + r * r / 16),
                0,
                1,
                epsabs=1e-13,
                epsrel=2e-12,
            )[0]
            assert low <= full * U / cutoff * (1 + 1e-9)
            count += 1
    # Here r=6,p=8,R=10 gives the included hemisphere u>=0.
    assert 6**2 + 8**2 == 10**2 and math.sqrt(1 + (10 / 4) ** 2) > 2
    unmasked = quad(lambda u: u / 2, -1, 1)[0]
    masked = quad(lambda u: u / 2, 0, 1)[0]
    assert abs(unmasked) < 1e-14 and abs(masked - 0.25) < 1e-14
    assert masked != unmasked
    print(
        "FULL_ORIGINAL_FREQUENCY_TWO_LEG_MASK_QUADRATURE_C0",
        count,
        "max_near_ratio",
        largest,
        flush=True,
    )
    print("ANGULAR_ZERO_UV_ROW_IS_NONZERO_BEFORE_MASK_REMOVAL", masked, flush=True)
    return count


def test_full_two_leg_removed_union_in_all_regions_and_angular_zero_control():
    assert mask_diagnostics() == 252


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def scale(a, c):
    return [c * x for x in a]


def mul(a, b):
    length = min(len(a), len(b))
    return [sum(a[k] * b[j - k] for k in range(j + 1)) for j in range(length)]


def inverse(a):
    out = [1 / a[0]]
    for j in range(1, len(a)):
        out.append(-sum(a[k] * out[j - k] for k in range(1, j + 1)) / a[0])
    return out


def root(a):
    out = [mp.sqrt(a[0])]
    for j in range(1, len(a)):
        out.append(
            (a[j] - sum(out[k] * out[j - k] for k in range(1, j))) / (2 * out[0])
        )
    return out


def derivative(a):
    return [(j + 1) * a[j + 1] for j in range(len(a) - 1)]


def constant(c, length):
    return [mp.mpc(c)] + [mp.mpf(0)] * (length - 1)


def normalized_frequency(time, x, mass, transfer, angle, dimension, leg):
    length = 21
    coord = [time, mp.mpf(1)] + [mp.mpf(0)] * (length - 2)
    vv = add(constant(1, length), mul(coord, coord))
    aa = mul(vv, vv)
    aa2 = mul(aa, aa)
    H = scale(mul(coord, inverse(vv)), 4)
    Ud = add(scale(derivative(H), dimension / 2), scale(mul(H, H), dimension**2 / 4))
    square = add(constant(1, length), scale(aa2, mass**2 * x**2))
    if leg == "l":
        square = add(
            square, constant(transfer**2 * x * x - 2 * transfer * x * angle, length)
        )
    W = root(square)
    for _ in range(6):
        rate = add(mul(derivative(W), inverse(W)), scale(H, -1))
        remainder = add(
            scale(Ud, -1),
            add(
                scale(derivative(rate), -mp.mpf(".5")),
                scale(mul(rate, rate), mp.mpf(".25")),
            ),
        )
        W = root(add(square, scale(mul(aa2, remainder), x * x)))
    full_rate = ((dimension - 1) * H[0] + W[1] / W[0]) / 2
    momentum = -mp.j * W[0] - aa[0] * x * full_rate
    return W[0], momentum, aa[0]


def complex_domain_diagnostics():
    fixtures = []
    for index, center in enumerate(("-0.5", "0", "0.5")):
        for ratio in (0, 1, 4, 100):
            for angle in ("-0.7", "0", "0.8"):
                fixtures.append(
                    (90, mp.mpf(10) ** 6, center, ratio, angle, index + ratio)
                )
    for index, (center, ratio, angle) in enumerate(
        (("-0.5", 0, "-0.7"), ("0", 1, "0"), ("0.5", 100, "0.8"))
    ):
        fixtures.append((1600, None, center, ratio, angle, index + 1))
    largest = mp.mpf(0)
    for precision, mass0, center, ratio, angle, phase in fixtures:
        with mp.workdps(precision):
            mass = mp.sqrt(mp.mpf(10) ** 200 / 512 + 2) if mass0 is None else mass0
            pp = mass * ratio
            uu = mp.mpf(angle)
            phi = mp.mpf(phase) / 7
            x = mp.mpf(".01") / (mass + pp) * mp.exp(mp.j * (phi + mp.mpf(".4")))
            time = mp.mpf(center) + mp.mpf(".99") / 64 * mp.exp(
                mp.j * (phi + mp.mpf(".7"))
            )
            dd = 3 + mp.mpf(".249") * mp.exp(mp.j * (phi + mp.mpf(".2")))
            Wk, Pk, aa = normalized_frequency(time, x, mass, pp, uu, dd, "k")
            Wl, Pl, _ = normalized_frequency(time, x, mass, pp, uu, dd, "l")
            assert mp.re(Wk / aa) > mp.mpf(".5") and mp.re(Wl / aa) > mp.mpf(".5")
            assert abs(Wk) < 2 and abs(Wl) < 2
            assert abs(Pk) < 3 and abs(Pl) < 3
            assert abs(aa / (Wk + Wl)) < 1
            nn = mp.matrix([mp.sqrt(1 - uu * uu), 0, uu])
            ll = -nn + mp.matrix([0, 0, pp * x])
            D = mp.matrix(
                [
                    [mp.mpf(".7"), mp.mpf(".2"), mp.mpf("-.1")],
                    [mp.mpf(".2"), mp.mpf("-.3"), mp.mpf(".15")],
                    [mp.mpf("-.1"), mp.mpf(".15"), mp.mpf(".4")],
                ]
            )
            trace = sum(D[j, j] for j in range(3))
            M = mp.zeros(5)
            M[0, 0] = -trace / 2
            M[4, 4] = trace / 2
            for i in range(3):
                for j in range(3):
                    M[i + 1, j + 1] = (trace / 2 if i == j else 0) - D[i, j]
            Fk = mp.matrix([Pk, *[mp.j * nn[j] for j in range(3)], mass * aa * x])
            Fl = mp.matrix([Pl, *[mp.j * ll[j] for j in range(3)], mass * aa * x])
            pair = (Fk.T * M * Fl)[0] / (2 * aa * mp.sqrt(Wk * Wl))
            norm = mp.sqrt(sum(D[i, j] ** 2 for i in range(3) for j in range(3)))
            assert abs(pair) < 1000 * norm
            largest = max(largest, abs(pair) / norm)
    print(
        "FULL_COMPLEX_DOMAIN_UNEXPANDED_SCALAR_VERTEX_C0",
        len(fixtures),
        "max_pair_norm",
        mp.nstr(largest, 20),
        flush=True,
    )
    return len(fixtures)


def test_full_six_iterate_complex_domain_scalar_vertex():
    assert complex_domain_diagnostics() == 39
