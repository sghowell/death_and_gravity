"""Independent full scalar UV, geometry, dimensional finite part and scope checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_curved_state import state
from p8_vacuum_affine_heavy_spatial_uv import audit, jets, matching, structure


@pytest.mark.parametrize("name,value", tuple(audit.residuals().items()))
def test_every_exact_scalar_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize("name,value", tuple(audit.gates().items()))
def test_every_written_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_or_claim_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_complete_slots_not_only_selected_logarithmic_coefficients():
    rows = jets.endpoint_products()
    assert set(rows) == {
        (geo, j, r)
        for geo in ("00", "01", "10", "11")
        for j in range(5)
        for r in range(j + 1)
    }
    assert sum(map(len, rows.values())) == 140
    assert set(jets.radial_rows()) == {
        (j, r, degree)
        for j in range(5)
        for r in range(j + 1)
        for degree in range(5 - j)
    }
    assert all(len(value) == 5 - key[1] for key, value in rows.items())


def test_full_ring_operations_against_literal_series():
    z = s.Symbol("series_parameter")
    v = (
        s.S.One,
        s.Rational(2, 3),
        s.Rational(7, 5),
        s.Rational(-3, 7),
        s.Rational(1, 11),
    )
    poly = sum(value * z**j for j, value in enumerate(v))
    for actual, expression in (
        (jets.reciprocal(v), 1 / poly),
        (jets.sqrtone(v), s.sqrt(poly)),
        (jets.power(v, 3), poly**3),
    ):
        expected = s.series(expression, z, 0, 5).removeO().expand()
        assert (
            s.expand(sum(value * z**j for j, value in enumerate(actual)) - expected)
            == 0
        )


def test_time_jet_representation_does_not_silently_differentiate_missing_H6():
    with pytest.raises(ValueError, match="H6"):
        jets.derivative(jets.h[-1])
    assert (
        jets.derivative(jets.a**2 * jets.h[2])
        == 2 * jets.a**2 * jets.h[0] * jets.h[2] + jets.a**2 * jets.h[3]
    )


@pytest.mark.parametrize("leg", ("x", None, True, 1, s.Rational(1, 2)))
def test_only_two_internal_legs(leg):
    with pytest.raises((ValueError, TypeError)):
        jets.frequency(leg)


def test_literal_all_spatial_metric_directions():
    data = structure.data()
    M = data["literal_complete_spatial_Hamiltonian_feature_matrix"]
    assert M.shape == (5, 5)
    assert s.factor(M[0, 0] + M[4, 4]) == 0
    assert M[1:4, 1:4] != s.zeros(3)


def test_whole_spatial_contact_cancels_only_in_transfer_difference():
    D = s.diag(1, 0, 0)
    G = s.diag(0, 1, 0)
    BD = D.trace() * s.eye(3) / 2 - D
    BG = G.trace() * s.eye(3) / 2 - G
    second = s.zeros(5)
    second[0, 0] = second[4, 4] = D.trace() * G.trace() / 4
    second[1:4, 1:4] = (BD * BG + BG * BD) / 2
    transfer = s.Symbol("transfer", real=True)
    assert second.subs(transfer, 0) == second
    assert second[0, 0] == s.Rational(1, 4)
    assert second != s.zeros(5)


def test_odd_ordered_endpoints_are_not_zero_in_the_complete_scalar_symbol():
    odd = [
        value
        for (j, _r, degree), value in jets.spatial_rows().items()
        if j % 2 and j + degree == 4
    ]
    assert any(value != 0 for value in odd)
    assert any(jets.endpoint_products()["01", 1, 0])
    assert jets.geometries()[1]["01"] != jets.geometries()[1]["10"]


def test_freezing_dimension_or_omitting_full_counteraction_changes_finite_part():
    data = matching.finite()
    assert any(value != 0 for value in data["complete_invariant_dimensional_slopes"])
    assert data["complete_covariant_evanescent_counteraction"] != 0
    assert s.diff(data["complete_covariant_evanescent_counteraction"], jets.inv[3]) != 0


def test_scalar_pole_is_not_old_vector_pole():
    _, curv = matching.geometry.hessians(*jets.inv)
    scalar = matching.finite()["own_scalar_pole"]
    old_vector = (
        jets.m**2 * curv["R_old"]
        - curv["R_squared"] / 4
        + s.Rational(29, 30) * curv["Ricci_squared"]
        - s.Rational(2, 15) * curv["Riemann_squared"]
    )
    assert (
        s.factor(
            s.diff(scalar - old_vector, jets.m, 2) + s.Rational(8, 3) * curv["R_old"]
        )
        == 0
    )
    assert s.factor(scalar - old_vector) != 0


def test_all_finite_coefficients_and_actual_mass_graph_bound():
    data = matching.finite()
    assert len(data["all_eighteen_fixed_invariant_coefficients"]) == 18
    c1 = data["full_coefficient_sum_mass_squared"]
    c0 = data["full_coefficient_sum_constant"]
    assert c1 == s.Rational(76796107, 23224320)
    assert c0 == s.Rational(9596179133, 464486400)
    assert (c1 * state.MASS2 + c0) / state.KAPPA < s.Rational(1, 10**600)
    assert (
        s.factor(
            data["complete_normalized_local_UV_difference_bound"]
            - (c1 * state.MASS2 + c0) / state.KAPPA
        )
        == 0
    )


def test_actual_log_interval_uses_elementary_exponential_bounds():
    assert state.MASS2 > 10**197 and state.MASS2 < 10**198
    # e<3 follows by j!>=2^(j-1) for j>=2, strictly after j3.
    assert (
        s.Rational(1) + 1 + sum(s.Rational(1, 2) ** (j - 1) for j in range(2, 20)) < 3
    )
    assert sum(s.Rational(7, 3) ** j / s.factorial(j) for j in range(7)) > 10
    assert 197 * 2 == 394 and 198 * s.Rational(7, 3) == 462


def test_original_scope_and_matching_records_unchanged():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 99
    assert "full UV-subtracted" in audit.observable()["remaining"]
    assert "not added again" in audit.observable()["remaining"]
    assert audit.require_parameters(state.MASS2, state.KAPPA) == (
        state.MASS2,
        state.KAPPA,
    )
    assert audit.require_dimension(3) == 3
    assert audit.require_graph("source_time", 2) == ("source_time", 2)


def test_independent_full_W6_endpoint_time_jets():
    assert len(run_full_wkb_diagnostics()) == 3


def test_independent_cartesian_geometry_in_four_dimensions():
    assert run_geometry_diagnostics() == 160


def add(left, right):
    return [x + y for x, y in zip(left, right)]


def scale(row, value):
    return [value * x for x in row]


def mul(left, right):
    length = min(len(left), len(right))
    return [sum(left[k] * right[j - k] for k in range(j + 1)) for j in range(length)]


def inverse(row):
    out = [1 / row[0]]
    for j in range(1, len(row)):
        out.append(-sum(row[k] * out[j - k] for k in range(1, j + 1)) / row[0])
    return out


def root(row):
    out = [mp.sqrt(row[0])]
    for j in range(1, len(row)):
        out.append(
            (row[j] - sum(out[k] * out[j - k] for k in range(1, j))) / (2 * out[0])
        )
    return out


def derivative(row):
    return [(j + 1) * row[j + 1] for j in range(len(row) - 1)]


def constant(value, length):
    return [mp.mpf(value)] + [mp.mpf(0)] * (length - 1)


def background(time, length=27):
    coord = [mp.mpf(time), mp.mpf(1)] + [mp.mpf(0)] * (length - 2)
    v = add(constant(1, length), mul(coord, coord))
    aa = mul(v, v)
    H = scale(mul(coord, inverse(v)), 4)
    return aa, H


def whole_frequency(momentum_squared, mass_squared, aa, H, dimension):
    omega2 = add(
        constant(mass_squared, len(aa)),
        scale(mul(inverse(aa), inverse(aa)), momentum_squared),
    )
    U = add(
        scale(derivative(H), dimension / 2), scale(mul(H, H), dimension * dimension / 4)
    )
    W = root(omega2)
    for _ in range(6):
        rate = mul(derivative(W), inverse(W))
        W = root(
            add(
                omega2,
                add(
                    scale(U, -1),
                    add(
                        scale(derivative(rate), -mp.mpf("0.5")),
                        scale(mul(rate, rate), mp.mpf("0.25")),
                    ),
                ),
            )
        )
    velocity = add(
        scale(W, -mp.j),
        scale(add(scale(H, dimension), mul(derivative(W), inverse(W))), -mp.mpf("0.5")),
    )
    return W, velocity


def whole_endpoints(time, dimension, mass, transfer, cosine, x):
    aa, H = background(time)
    radial = 1 / x
    Wk, Pk = whole_frequency(radial * radial, mass * mass, aa, H, dimension)
    Wl, Pl = whole_frequency(
        radial * radial - 2 * radial * transfer * cosine + transfer * transfer,
        mass * mass,
        aa,
        H,
        dimension,
    )
    norm = scale(inverse(root(mul(Wk, Wl))), mp.mpf("0.5"))
    dot = -radial * radial + radial * transfer * cosine
    tr_num = add(
        scale(mul(Pk, Pl), -1),
        add(scale(mul(inverse(aa), inverse(aa)), -dot), constant(mass * mass, len(aa))),
    )
    U = {
        "trace": scale(mul(norm, tr_num), mp.mpf("0.5")),
        "gradient": scale(mul(norm, mul(inverse(aa), inverse(aa))), radial * radial),
    }
    g = inverse(add(Wk, Wl))
    out = {}
    for name, left, right in (
        ("00", "trace", "trace"),
        ("01", "trace", "gradient"),
        ("10", "gradient", "trace"),
        ("11", "gradient", "gradient"),
    ):
        source = {0: U[right]}
        for j in range(5):
            for r, row in source.items():
                out[name, j, r] = mp.im(
                    mp.j * (-mp.j) ** j * mp.conj(U[left][0]) * g[0] * row[0]
                ) * x ** (1 - j)
            following = {}
            for r, row in source.items():
                product = mul(g, row)
                following[r] = (
                    add(following[r], derivative(product))
                    if r in following
                    else derivative(product)
                )
                following[r + 1] = (
                    add(following[r + 1], product) if r + 1 in following else product
                )
            source = following
    return out


def run_full_wkb_diagnostics():
    import sympy as s
    from p8_vacuum_affine_heavy_spatial_uv import jets

    tt = s.Symbol("fixture_time", real=True)
    HH = 4 * tt / (1 + tt * tt)
    args = (jets.a, *jets.h, jets.d, jets.m, jets.p, jets.u)
    compiled = {
        key: s.lambdify(args, row, "mpmath")
        for key, row in jets.endpoint_products().items()
    }
    with mp.workdps(100):
        maxima = []
        fixtures = (
            ("0.2", "3", "2.7", "0.7", "0.3"),
            ("-0.37", "4", "7", "1.3", "-0.45"),
            ("0.31", "3.4", "0.9", "2.1", "0.61"),
        )
        for fixture in fixtures:
            time, dim, mass, p, u = map(mp.mpf, fixture)
            hh = [s.lambdify(tt, s.diff(HH, tt, j), "mpmath")(time) for j in range(6)]
            values = ((1 + time * time) ** 2, *hh, dim, mass, p, u)
            rows = {key: fun(*values) for key, fun in compiled.items()}
            ratios = []
            for x in (mp.mpf("1e-4"), mp.mpf("5e-5")):
                full = whole_endpoints(time, dim, mass, p, u, x)
                normalized = []
                for key, coeffs in rows.items():
                    expected = sum(c * x**degree for degree, c in enumerate(coeffs))
                    error = abs(full[key] - expected) / x ** len(coeffs)
                    assert error < mp.mpf("1e10"), (fixture, key, mp.nstr(error, 15))
                    normalized.append(error)
                ratios.append(max(normalized))
            assert abs(ratios[1] - ratios[0]) < mp.mpf("0.03") * max(1, ratios[0]), (
                fixture,
                ratios,
            )
            maxima.append(ratios[1])
            print(
                "FULL_W6_ALL_60_ORDERED_ENDPOINT_FIXTURE",
                fixture,
                mp.nstr(ratios[1], 18),
                flush=True,
            )
            # Independently recover individual radial Taylor coefficients from
            # eight FULL positive-radius values, without the compiled ring.
            spacing = mp.mpf("1e-10")
            nodes = list(range(1, 9))
            vandermonde = mp.matrix(
                [[mp.mpf(node) ** degree for degree in range(8)] for node in nodes]
            )
            keys = list(rows)
            samples = [
                whole_endpoints(time, dim, mass, p, u, spacing * node) for node in nodes
            ]
            matrix = mp.matrix([[sample[key] for key in keys] for sample in samples])
            recovered = (vandermonde**-1) * matrix
            coefficient_checks = 0
            for column, key in enumerate(keys):
                for degree, expected in enumerate(rows[key]):
                    actual = recovered[degree, column] / spacing**degree
                    assert abs(actual - expected) < mp.mpf("1e-25") * max(
                        1, abs(expected)
                    ), (fixture, key, degree, mp.nstr(actual - expected, 15))
                    coefficient_checks += 1
            assert coefficient_checks == 140
            print(
                "FULL_W6_INDIVIDUAL_UV_COEFFICIENTS_C0",
                fixture,
                coefficient_checks,
                flush=True,
            )
        wrong = jets.endpoint_products()["01", 1, 0][0:4]
        assert any(value != 0 for value in wrong)
        print(
            "FULL_NUMERICAL_ENDPOINTS_C0", len(fixtures) * 2 * len(compiled), flush=True
        )
        return maxima


def cartesian_average(value, coordinates):
    dimension = len(coordinates)
    total = s.S.Zero
    for powers, coefficient in s.Poly(s.expand(value), *coordinates).terms():
        if any(power % 2 for power in powers):
            continue
        half = sum(powers) // 2
        numerator = s.prod(s.factorial2(power - 1) for power in powers)
        denominator = s.prod(dimension + 2 * j for j in range(half))
        total += coefficient * numerator / denominator
    return s.factor(total)


def run_geometry_diagnostics():
    from p8_vacuum_affine_heavy_spatial_uv import jets

    y, rows = jets.geometries()
    checks = 0
    for dimension in (3, 4, 5, 6):
        n = s.Matrix(s.symbols("n0:" + str(dimension), real=True))
        e = s.eye(dimension)[:, 0]
        ll = y * e - n
        for seed in (1, 2):
            D = s.Matrix(
                dimension,
                dimension,
                lambda i, j, seed=seed: (
                    s.Rational(((i + j + seed) % 5) - 2, 3) if i != j else i + seed
                ),
            )
            G = s.Matrix(
                dimension,
                dimension,
                lambda i, j, seed=seed: (
                    s.Rational(((i + j + 2 * seed) % 7) - 3, 5)
                    if i != j
                    else 2 - i - seed
                ),
            )
            invariants = (
                (D * G).trace(),
                (e.T * D * G * e)[0],
                (e.T * D * e)[0] * (e.T * G * e)[0],
                D.trace() * G.trace(),
                D.trace() * (e.T * G * e)[0],
                G.trace() * (e.T * D * e)[0],
            )
            factors = {
                "00": D.trace() * G.trace(),
                "01": D.trace() * (n.T * G * ll)[0],
                "10": G.trace() * (n.T * D * ll)[0],
                "11": (n.T * D * ll)[0] * (n.T * G * ll)[0],
            }
            substitutions = {jets.d: dimension, **dict(zip(jets.inv, invariants))}
            for label, value in factors.items():
                for degree in range(5):
                    direct = cartesian_average(n[0] ** degree * value, tuple(n))
                    compiled = jets.sphere(jets.u**degree * rows[label]).subs(
                        substitutions, simultaneous=True
                    )
                    residual = s.factor(direct - compiled)
                    assert residual == 0, (dimension, seed, label, degree, residual)
                    checks += 1
    print("LITERAL_CARTESIAN_ORDERED_GEOMETRY_C0", checks, flush=True)
    return checks
