"""Independent covariant Hessians, spectral integrals and retarded extraction."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_flat_spatial_conversion import bounds, conversion, flat
from p8_vacuum_affine_flat_spatial_conversion import tensor as ten

ETA = s.diag(1, -1, -1, -1)


def q_matrix(spin, z, p):
    A = z * mp.eye(3) + p * p.T
    result = mp.zeros(9)
    for u, v in product(range(9), repeat=2):
        a, b = divmod(u, 3)
        c, d = divmod(v, 3)
        trace = A[a, b] * A[c, d] / 3
        result[u, v] = (
            trace if spin == 0 else (A[a, c] * A[b, d] + A[a, d] * A[b, c]) / 2 - trace
        )
    return result


def q_coefficients(spin, p, center=0):
    zero = q_matrix(spin, center, p)
    plus = q_matrix(spin, center + 1, p)
    minus = q_matrix(spin, center - 1, p)
    return zero, (plus - minus) / 2, (plus + minus) / 2 - zero


def density_polynomial(spin, v):
    x = (1 - v * v) / 4
    return (
        (13 + 56 * x + 48 * x * x) / 3840
        if spin == 2
        else (1 - 4 * x + 12 * x * x) / 384
    )


def radial_J(spin, q, n=0, lower=0):
    def integrand(v):
        t = 1 - v * v
        return (
            2
            * v
            * v
            * density_polynomial(spin, v)
            * t**n
            / (mp.pi**2 * (4 + q * t) ** (n + 1))
        )

    return mp.quad(integrand, [lower, 1])


def covariant_full(p, w, lower=0):
    q = (p.T * p)[0]
    z = w - q
    return sum(
        (z * q_matrix(spin, z, p) * radial_J(spin, -z, lower=lower) for spin in (0, 2)),
        mp.zeros(9),
    )


def conversion_coefficients(p, lower=0):
    q = (p.T * p)[0]
    result = [mp.zeros(9) for _ in range(3)]
    for spin in (0, 2):
        H0, H1, H2 = q_coefficients(spin, p, center=-q)
        n = (-q * H0, H0 - q * H1, H1 - q * H2, H2)
        Js = [radial_J(spin, q, j, lower) for j in range(3)]
        for r in range(3):
            result[r] += sum((n[k] * Js[r - k] for k in range(r + 1)), mp.zeros(9))
    return result


def time_subtracted_full(p, w):
    q = (p.T * p)[0]
    result = mp.zeros(9)
    for spin in (0, 2):
        H = q_coefficients(spin, p)
        for r in range(3):

            def integrand(v, spin=spin, r=r):
                t = 1 - v * v
                return (
                    128
                    * v
                    * v
                    * density_polynomial(spin, v)
                    * (t / 4) ** (2 - r)
                    / (mp.pi**2 * (4 + q * t) ** 3 * (4 + (q - w) * t))
                )

            result += w**3 * H[r] * mp.quad(integrand, [0, 1])
    return result


def matrix_norm(A):
    return max(abs(value) for value in mp.eigsy(A, eigvals_only=True))


@pytest.mark.parametrize("direction", ((0, 0, 0), (1, 2, -1), (2, -3, 4), (-1, 1, 3)))
@pytest.mark.parametrize("spin", (0, 2))
def test_full_covariant_projector_at_nonzero_spatial_transfer(direction, spin):
    p = s.Matrix(direction)
    sigma = s.Integer(7)
    P = s.Matrix([s.sqrt(sigma + p.dot(p)), *p])
    theta = ETA - P * P.T / sigma
    direct = s.zeros(9)
    for u, v in product(range(9), repeat=2):
        a, b = divmod(u, 3)
        c, d = divmod(v, 3)
        a += 1
        b += 1
        c += 1
        d += 1
        trace = theta[a, b] * theta[c, d] / 3
        direct[u, v] = (
            trace
            if spin == 0
            else (theta[a, c] * theta[b, d] + theta[a, d] * theta[b, c]) / 2 - trace
        )
    actual = ten.numerator(spin, sigma).subs(dict(zip(ten.P, p))) / sigma**2
    assert actual == direct
    assert sigma * direct * sigma == ten.numerator(spin, sigma).subs(
        dict(zip(ten.P, p))
    )


@pytest.mark.parametrize("direction", ((0, 0, 0), (1, 2, -1), (2, -3, 4), (-1, 1, 3)))
def test_independent_full_four_curvature_Hessians_at_nonzero_transfer(direction):
    frequency = s.Symbol("frequency", real=True)
    momentum = s.Matrix([frequency, *direction])
    D = s.Matrix([[1, 2, -1], [2, -3, 1], [-1, 1, 4]])
    G = s.Matrix([[2, -1, 3], [-1, 1, 2], [3, 2, -2]]) / 7

    def linear_curvature(H):
        h = s.zeros(4)
        h[1:4, 1:4] = H
        R = {}
        for a, b, c, d in product(range(4), repeat=4):
            R[a, b, c, d] = (
                -(
                    momentum[c] * momentum[b] * h[a, d]
                    + momentum[d] * momentum[a] * h[b, c]
                    - momentum[d] * momentum[b] * h[a, c]
                    - momentum[c] * momentum[a] * h[b, d]
                )
                / 2
            )
        Ric = s.Matrix(
            4, 4, lambda a, b: sum(ETA[r, r] * R[r, a, r, b] for r in range(4))
        )
        return R, Ric, s.trace(ETA * Ric)

    A, RA, sa = linear_curvature(D)
    B, RB, sb = linear_curvature(G)
    riem = sum(
        ETA[a, a] * ETA[b, b] * ETA[c, c] * ETA[d, d] * A[a, b, c, d] * B[a, b, c, d]
        for a, b, c, d in product(range(4), repeat=4)
    )
    ric = sum(
        ETA[a, a] * ETA[b, b] * RA[a, b] * RB[a, b]
        for a, b in product(range(4), repeat=2)
    )
    substitutions = {
        ten.Z: frequency**2 - sum(v * v for v in direction),
        **dict(zip(ten.P, direction)),
    }
    vD = s.Matrix(list(D))
    vG = s.Matrix(list(G))
    tensor = (vD.T * ten.numerator(2) * vG)[0].subs(substitutions)
    scalar = (vD.T * ten.numerator(0) * vG)[0].subs(substitutions)
    assert s.expand(2 * (riem - 2 * ric + sa * sb / 3) - tensor) == 0
    assert s.expand(2 * sa * sb - 6 * scalar) == 0


@pytest.mark.parametrize(
    "direction", ((0, 0, 0), ("0.1", "0.2", "-0.1"), (1, 2, -1), (10, -2, 3))
)
@pytest.mark.parametrize("point", ("-1", "1", "3.8", "complex"))
def test_independent_full_matrix_improper_spectral_conversion(direction, point):
    with mp.workdps(65):
        p = mp.matrix([mp.mpf(v) for v in direction])
        w = mp.mpc(1, 1) if point == "complex" else mp.mpf(point)
        F = covariant_full(p, w)
        A = conversion_coefficients(p)
        B = time_subtracted_full(p, w)
        error = F - sum((A[r] * w**r for r in range(3)), mp.zeros(9)) - B
        assert mp.norm(error) < mp.mpf("1e-55") * max(1, mp.norm(F), mp.norm(B))


@pytest.mark.parametrize(
    "direction", ((0, 0, 0), ("0.1", "0.2", "-0.1"), (1, 2, -1), (10, -2, 3))
)
@pytest.mark.parametrize("order", (0, 1, 2))
def test_independent_all_spatial_coefficient_operator_norm(direction, order):
    with mp.workdps(60):
        p = mp.matrix([mp.mpf(v) for v in direction])
        q = (p.T * p)[0]
        A = conversion_coefficients(p)[order]
        assert matrix_norm(A) <= (q + 1) ** 3 / (80 * mp.pi**2)
        if q == 0:
            assert mp.norm(A) == 0


@pytest.mark.parametrize("direction", (("0.1", "0.2", "-0.1"), (1, 2, -1), (10, -2, 3)))
@pytest.mark.parametrize("limit", (4, 8, 100))
@pytest.mark.parametrize("order", (0, 1, 2))
def test_independent_full_conversion_infinite_spectral_tail(direction, limit, order):
    with mp.workdps(60):
        p = mp.matrix([mp.mpf(v) for v in direction])
        q = (p.T * p)[0]
        start = mp.sqrt(1 - mp.mpf(4) / limit)
        tail = conversion_coefficients(p, start)[order]
        assert matrix_norm(tail) <= (q + 1) ** 3 / (9 * mp.pi**2 * limit)


@pytest.mark.parametrize("energy", (2, 3, 10))
@pytest.mark.parametrize("time", ("0.0625", "0.125", "0.1875"))
def test_literal_prepared_retarded_sine_kernel_extraction(energy, time):
    with mp.workdps(65):
        t = s.Symbol("t")
        polynomial = t**7 * (s.Rational(1, 4) - t) ** 7 * (1 + t + t * t)
        derivatives = [
            s.lambdify(t, s.diff(polynomial, t, j), "mpmath") for j in range(7)
        ]
        E = mp.mpf(energy)
        now = mp.mpf(time)
        response = mp.quad(
            lambda u: mp.sin(E * (now - u)) / E * derivatives[0](u), [0, now]
        )
        remainder = (
            -mp.quad(lambda u: mp.sin(E * (now - u)) / E * derivatives[6](u), [0, now])
            / E**6
        )
        boundary = (
            derivatives[0](now) / E**2
            - derivatives[2](now) / E**4
            + derivatives[4](now) / E**6
        )
        assert abs(response - boundary - remainder) < mp.mpf("1e-60")
        assert all(s.diff(polynomial, t, j).subs(t, 0) == 0 for j in range(7))


@pytest.mark.parametrize("order", range(6))
def test_nonlocal_TT_series_has_independently_positive_coefficients(order):
    with mp.workdps(60):

        def integral(sigma):
            x = 1 / sigma
            rho = (
                mp.sqrt(1 - 4 * x)
                * sigma**2
                * (13 + 56 * x + 48 * x * x)
                / (3840 * mp.pi**2)
            )
            return rho / sigma ** (4 + order)

        value = mp.quad(integral, [4, 8, 32, mp.inf])
        from p8_vacuum_affine_flat_tensor_cut import dispersion, projectors

        exact = dispersion.moment(2, 3 + order).subs(projectors.MASS, 1)
        assert value > 0
        assert abs(value - mp.mpf(str(s.N(exact, 65)))) < mp.mpf("1e-55")


@pytest.mark.parametrize("bad", (True, False, 2.0, s.Float(2), 1, 3, "2", None))
def test_nonexact_spin_selector_rejected_before_cache(bad):
    ten.numerator(2)
    ten.centered_coefficients(2)
    with pytest.raises(ValueError):
        ten.numerator(bad)
    with pytest.raises(ValueError):
        ten.centered_coefficients(bad)


@pytest.mark.parametrize("module", (ten, conversion, bounds, flat))
def test_all_scientific_packets(module):
    d = module.data()
    for value in d["checks"].values():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(x) == 0 for x in entries)
    assert all(d["gates"].values())


@pytest.mark.parametrize("spin", (0, 2))
@pytest.mark.parametrize("direction", ((0, 0, 0), (1, 2, -1), (2, -3, 4), (-1, 1, 3)))
@pytest.mark.parametrize("order", (0, 1, 2, 3))
def test_independent_centered_interpolation_matches_literal_frequency_derivatives(
    spin, direction, order
):
    with mp.workdps(65):
        p = mp.matrix(direction)
        q = (p.T * p)[0]
        H0, H1, H2 = q_coefficients(spin, p, center=-q)
        n = (-q * H0, H0 - q * H1, H1 - q * H2, H2)
        polynomial = ten.cubic(spin)
        exact = polynomial.diff(ten.Z, order).subs(
            {ten.Z: -sum(v * v for v in direction), **dict(zip(ten.P, direction))}
        ) / s.factorial(order)
        direct = mp.matrix(
            [[mp.mpf(str(s.N(exact[i, j], 68))) for j in range(9)] for i in range(9)]
        )
        assert mp.norm(n[order] - direct) < mp.mpf("1e-55") * max(1, mp.norm(direct))


def absolute_sine_spectral_majorant(p, lower=4):
    q = (p.T * p)[0]
    theta0 = mp.asin(mp.sqrt(1 - mp.mpf(4) / lower))
    result = mp.zeros(9)
    for spin in (0, 2):
        H = q_coefficients(spin, p)
        for r in range(3):

            def integrand(theta, spin=spin, r=r):
                v = mp.sin(theta)
                t = mp.cos(theta) ** 2
                return (
                    128
                    * v
                    * v
                    * density_polynomial(spin, v)
                    * (t / 4) ** (2 - r)
                    / (mp.pi**2 * (4 + q * t) ** mp.mpf("3.5"))
                )

            result += H[r] * mp.quad(integrand, [theta0, mp.pi / 2])
    return result


@pytest.mark.parametrize(
    "direction", ((0, 0, 0), ("0.1", "0.2", "-0.1"), (1, 2, -1), (10, -2, 3))
)
@pytest.mark.parametrize("lower", (4, 8, 100))
def test_independent_complete_flat_absolute_sine_bulk_and_tail(direction, lower):
    with mp.workdps(60):
        p = mp.matrix([mp.mpf(v) for v in direction])
        q = (p.T * p)[0]
        actual = absolute_sine_spectral_majorant(p, lower)
        assert matrix_norm(actual) <= 7 / (192 * mp.pi**2 * mp.sqrt(lower + q))
        assert matrix_norm(actual) <= mp.mpf(1) / (100 * mp.sqrt(lower))


from p8_vacuum_affine_flat_spatial_conversion import audit


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_exact_residuals(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.cancel(x) == 0 for x in entries)


@pytest.mark.parametrize("case", audit.bad_cases(), ids=lambda case: case[0])
def test_all_scope_mutations_rejected(case):
    _, call, args = case
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_counts_and_guards():
    assert len(audit.residuals()) == 34
    assert audit.scalar_entry_count() == 994
    assert len(audit.gates()) == 35 and all(v is True for v in audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 122


def test_full_flat_norm_keeps_its_local_curved_and_inverse_boundary():
    assert "actual CD state is not changed" in audit.observable()["domain"]
    assert "physical local polynomial" in bounds.data()["full_nonlocal_response"]
    assert "spatially nonlocal" in audit.observable()["boundary"]
    assert "not a fully reduced mixed norm" in bounds.data()["canonical"]
    assert len(audit.frontier()) == 9
