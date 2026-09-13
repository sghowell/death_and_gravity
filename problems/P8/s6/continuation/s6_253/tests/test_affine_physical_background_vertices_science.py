"""Independent parent-path, full constraint and Gaussian operator diagnostics."""

from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_affine_vacuum_domain import family
from p8_vacuum_affine_physical_background_vertices import audit, parent, physical
from p8_vacuum_affine_physical_background_vertices import coupled as c
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import expm_multiply


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_every_complete_exact_identity(name, value):
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    ), name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_every_scoped_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_and_scope_promotion(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


u, X = family.u, family.X
N = s.Symbol("independent_lapse", positive=True)
e = s.Symbol("independent_probe", real=True)
h = (1 + u * u) ** 3
R = 1 + (X - 1) / h
Ru, RX = s.diff(R, u), s.diff(R, X)
energy = (2 + u + u * u) / 100001
pressure = -(1 - 2 * u + u**3) / 100003
F = (
    family.original.data()["original_retuned_tree_scalar"]
    - pressure
    - (energy + pressure) * (X - 1) / 2
)
Rnu = R.subs(X, N**-2)
primitive_N = (-3 * X * Ru * RX * R ** -s.Rational(7, 4) / (2 * N * N)).subs(X, N**-2)
primitive = sum(
    s.diff(primitive_N, N, j).subs(N, 1) * (N - 1) ** (j + 1) / s.factorial(j + 1)
    for j in range(4)
)
B = (-(R ** -s.Rational(3, 4)) * Ru / (2 * N)).subs(X, N**-2) - primitive
Fhat = (R ** -s.Rational(3, 4) * (F + 9 * X * Ru**2 / (16 * R))).subs(
    X, N**-2
) - s.diff(primitive, u) / N
D = Rnu ** s.Rational(1, 4) / N
Z = Rnu ** -s.Rational(3, 4) / N
C3 = Rnu ** s.Rational(3, 4) / 2
Cchi = Rnu ** -s.Rational(1, 4)
full_rows = physical.coefficient_map()["parameters"]
bound_rows = {
    key: tuple(
        x.subs(
            {parent.RHO: energy, parent.PRESSURE: pressure}, simultaneous=True
        ).doit()
        for x in row
    )
    for key, row in full_rows.items()
}


@cache
def direct_physical_fixture(case):
    point = (s.Rational(-1, 4), s.S.Zero, s.Rational(1, 3))[case]
    lapse_probe = [s.Rational(case + j + 1, 7 + j) for j in range(3)]
    scale_probe = [s.Rational(2 * case - j + 1, 11 + j) for j in range(3)]
    matter_probe = [s.Rational(case + 2 * j - 1, 17 + j) for j in range(3)]
    npath = 1 + e * lapse_probe[0]
    ndot, nddot = e * lapse_probe[1], e * lapse_probe[2]
    evaluate = lambda x: x.subs(u, point).subs(N, npath)
    rv = evaluate(Rnu)
    rdot = evaluate(s.diff(Rnu, u)) + evaluate(s.diff(Rnu, N)) * ndot
    rddot = (
        evaluate(s.diff(Rnu, u, 2))
        + 2 * evaluate(s.diff(Rnu, u, N)) * ndot
        + evaluate(s.diff(Rnu, N, 2)) * ndot**2
        + evaluate(s.diff(Rnu, N)) * nddot
    )
    aa = (1 + point * point) ** 2 * s.exp(e * scale_probe[0]) * rv ** s.Rational(1, 4)
    hh = parent.H.subs(u, point) + e * scale_probe[1] + rdot / (4 * rv)
    hhdot = (
        s.diff(parent.H, u).subs(u, point)
        + e * scale_probe[2]
        + rddot / (4 * rv)
        - rdot**2 / (4 * rv**2)
    )
    matter = parent.ell.subs(u, point) + e * matter_probe[1]
    matterdot = s.diff(parent.ell, u).subs(u, point) + e * matter_probe[2]
    dd, zz = evaluate(D), evaluate(Z)
    dn, dnn = evaluate(s.diff(D, N)), evaluate(s.diff(D, N, 2))
    zn, znn = evaluate(s.diff(Z, N)), evaluate(s.diff(Z, N, 2))
    bb, bn, bnn = evaluate(B), evaluate(s.diff(B, N)), evaluate(s.diff(B, N, 2))
    ff, fn, fnn = (
        evaluate(Fhat),
        evaluate(s.diff(Fhat, N)),
        evaluate(s.diff(Fhat, N, 2)),
    )
    theta, ww = -hh * dn + bn / 2, matter * zn
    L0 = -3 * dd * hh**2 + 3 * bb * hh + npath * ff + zz * matter**2 / 2
    Cnn = (
        -3 * hh**2 * dnn + 3 * hh * bnn + 2 * fn + npath * fnn + znn * matter**2 / 2
    ) / 2
    Qn = -3 * hh**2 * dn + 3 * hh * bn + ff + npath * fn + zn * matter**2 / 2
    ddot = evaluate(s.diff(D, u)) + dn * ndot
    zdot = evaluate(s.diff(Z, u)) + zn * ndot
    bdot = evaluate(s.diff(B, u)) + bn * ndot
    cv = -18 * dd * hh + 9 * bb
    cvdot = -18 * ddot * hh - 18 * dd * hhdot + 9 * bdot
    momentum = s.Integer(case + 2)
    zeta = s.Rational(1, 10**6)
    qq = momentum**2 / aa**2
    direct = {
        c.D: dd,
        c.Z: zz,
        c.J: Cnn + 3 * theta**2 / dd - ww**2 / (2 * zz),
        c.K: zeta * evaluate(Cchi) / (npath * aa**2),
        c.th: theta,
        c.w: ww,
        c.c: zz * matter,
        c.Lnv: 3 * Qn + 4 * evaluate(C3 + N * s.diff(C3, N)) * qq,
        c.Vvv: s.Rational(9, 2) * L0 - (cvdot + 3 * hh * cv) / 2,
        c.Vvs: -3 * (zdot * matter + zz * matterdot + 3 * hh * zz * matter),
        c.C: npath * evaluate(C3),
        c.Y: npath * evaluate(Cchi),
        c.Yv: npath * evaluate(Cchi) / aa**2,
        c.r: rv - 1,
        c.rN: evaluate(s.diff(Rnu, N)),
        c.dH: hh - parent.H.subs(u, point),
        c.q: qq,
        physical.a: aa,
    }
    fixture = {
        u: point,
        physical.a: (1 + point * point) ** 2,
        c.P: momentum,
        physical.zeta: zeta,
    }
    fixture.update(zip(physical.eta, lapse_probe))
    fixture.update(zip(physical.vp, scale_probe))
    fixture.update(zip(physical.psi, matter_probe))
    return direct, fixture


@pytest.mark.parametrize("case", range(3))
@pytest.mark.parametrize("key", list(full_rows))
@pytest.mark.parametrize("order", range(3))
def test_independent_complete_parent_physical_path_derivative(case, key, order):
    direct, fixture = direct_physical_fixture(case)
    actual = s.diff(direct[key], e, order).subs(e, 0)
    expected = bound_rows[key][order].subs(fixture, simultaneous=True)
    assert s.factor(actual - expected) == 0


def test_independent_complete_raw_scalar_ADM_and_time_boundary():
    eps = s.Symbol("eps")
    N, a = s.symbols("positive_background_lapse positive_hat_scale", positive=True)
    H, m, q = s.symbols(
        "background_H background_matter_rate momentum_squared", real=True
    )
    n, v, vd, sigma, sd, b = s.symbols("n v vdot sigma sigmadot b", real=True)
    beta, vx, sx, vxx = s.symbols("shift v_x sigma_x v_xx", real=True)
    coefficient_jets = {
        name: tuple(s.Symbol(name + "_" + str(j), real=True) for j in range(3))
        for name in ("M", "B", "F", "U", "C3", "Cchi")
    }
    coefficient = {
        name: sum(jets[j] * (eps * n) ** j / s.factorial(j) for j in range(3))
        for name, jets in coefficient_jets.items()
    }
    lapse = N + eps * n
    volume = 1 + 3 * eps * v + s.Rational(9, 2) * eps**2 * v**2
    curvature_volume = 1 - 2 * eps * v + 2 * eps**2 * v**2
    rate = H + eps * vd - eps**2 * beta * vx
    matter = m + eps * sd - eps**2 * beta * sx
    literal = volume * (
        -3 * coefficient["M"] * rate**2 / lapse
        + 2 * coefficient["M"] * rate * eps * b / lapse
        + coefficient["B"] * (3 * rate - eps * b)
        + lapse * coefficient["F"]
        + lapse
        * coefficient["C3"]
        * curvature_volume
        * (-4 * eps * vxx - 2 * eps**2 * vx**2)
        / a**2
        + coefficient["U"] * matter**2 / (2 * lapse)
        - lapse * coefficient["Cchi"] * curvature_volume * eps**2 * sx**2 / (2 * a**2)
    )
    raw = s.expand(s.diff(literal, eps, 2).subs(eps, 0) / 2)
    spatial = s.expand(
        raw.subs(
            {
                beta * vx: -v * b,
                beta * sx: -sigma * b,
                vx**2: a * a * q * v * v,
                sx**2: a * a * q * sigma * sigma,
                vxx: -a * a * q * v,
            },
            simultaneous=True,
        )
    )
    M0, M1, M2 = coefficient_jets["M"]
    B0, B1, B2 = coefficient_jets["B"]
    F0, F1, F2 = coefficient_jets["F"]
    U0, U1, U2 = coefficient_jets["U"]
    C0, C1, _ = coefficient_jets["C3"]
    MatterC0 = coefficient_jets["Cchi"][0]
    D, DN, DNN = M0 / N, M1 / N - M0 / N**2, M2 / N - 2 * M1 / N**2 + 2 * M0 / N**3
    Z, ZN, ZNN = U0 / N, U1 / N - U0 / N**2, U2 / N - 2 * U1 / N**2 + 2 * U0 / N**3
    L0 = -3 * D * H**2 + 3 * B0 * H + N * F0 + Z * m * m / 2
    theta = -H * DN + B1 / 2
    w = m * ZN
    Cnn = (-3 * H**2 * DNN + 3 * H * B2 + 2 * F1 + N * F2 + ZNN * m * m / 2) / 2
    Qn = -3 * H**2 * DN + 3 * H * B1 + F0 + N * F1 + ZN * m * m / 2
    cv = -18 * D * H + 9 * B0
    expected = (
        -3 * D * vd**2
        + Z * sd**2 / 2
        + 6 * theta * n * vd
        + w * n * sd
        + Cnn * n * n
        + cv * v * vd
        + 3 * Z * m * v * sd
        + 3 * Qn * n * v
        + s.Rational(9, 2) * L0 * v * v
        + b * (2 * D * vd - 2 * theta * n + Z * m * sigma)
        + 2 * N * C0 * q * v * v
        + 4 * (C0 + N * C1) * q * n * v
        - N * MatterC0 * q * sigma * sigma / 2
    )
    assert s.expand(spatial - expected) == 0

    # A separately parameterized entire normalized action avoids identifying the
    # off-reference matter Ward contact or kinetic factors with their clock values.
    D, Z, J = s.symbols("positive_D positive_Z positive_J", positive=True)
    theta, w, c, Lnv, Vvv, Vvs, C, Y = s.symbols(
        "Theta mixed_momentum matter_charge lapse_v potential_vv potential_vsigma curvature matter_gradient",
        real=True,
    )
    pv, ps = s.symbols("pv ps", real=True)
    L = (
        -3 * D * vd**2
        + Z * sd**2 / 2
        + 6 * theta * n * vd
        + w * n * sd
        + (J - 3 * theta**2 / D + w**2 / (2 * Z)) * n * n
        - 3 * c * vd * sigma
        + Lnv * n * v
        + Vvv * v * v
        + Vvs * v * sigma
        + b * (2 * D * vd - 2 * theta * n + c * sigma)
        + 2 * C * q * v * v
        - Y * q * sigma * sigma / 2
    )
    rates = {
        vd: -(pv - 6 * theta * n + 3 * c * sigma - 2 * D * b) / (6 * D),
        sd: (ps - w * n) / Z,
    }
    ham = s.cancel((pv * vd + ps * sd - L).subs(rates, simultaneous=True))
    aux = s.hessian(ham, [n, b]).applyfunc(s.cancel)
    assert aux == s.diag(-2 * J, -2 * D / 3)
    numerator = theta * (pv + 3 * c * sigma) / D - w * ps / Z - Lnv * v
    solution = {n: numerator / (2 * J), b: pv / (2 * D)}
    assert s.cancel(s.diff(ham, n).subs(solution, simultaneous=True)) == 0
    assert s.cancel(s.diff(ham, b).subs(solution, simultaneous=True)) == 0
    reduced = s.cancel(ham.subs(solution, simultaneous=True))
    target = (
        ps * ps / (2 * Z)
        - c * pv * sigma / (2 * D)
        - 3 * c * c * sigma * sigma / (4 * D)
        + Y * q * sigma * sigma / 2
        - (2 * C * q + Vvv) * v * v
        - Vvs * v * sigma
        + numerator**2 / (4 * J)
    )
    assert s.cancel(reduced - target) == 0

    # The source-pinned current reference is exactly the S251 expression after
    # D=Z=1,c=ell,Y=1,C=1/2,Lnv=2E*q+3T,Vvs=0,Vvv=9A/2.
    ell, E, A, T = s.symbols("ell E A Tcorr", real=True)
    at_reference = {
        D: 1,
        Z: 1,
        c: ell,
        Y: 1,
        C: s.Rational(1, 2),
        Lnv: 2 * E * q + 3 * T,
        Vvs: 0,
        Vvv: s.Rational(9, 2) * A,
        w: -ell * E,
    }
    reference = (
        ps**2 / 2
        - ell * pv * sigma / 2
        + (q / 2 - 3 * ell**2 / 4) * sigma**2
        - q * v * v
        - s.Rational(9, 2) * A * v * v
        + (
            theta * pv
            + ell * E * ps
            + 3 * theta * ell * sigma
            - (2 * E * q + 3 * T) * v
        )
        ** 2
        / (4 * J)
    )
    assert s.cancel(target.subs(at_reference, simultaneous=True) - reference) == 0


def test_independent_complete_held_vector_spatial_Fourier_expansion():
    e, cs, sn = s.symbols("epsilon cosine sine", real=True)
    Z, ZN, ZNN, r, rN, rNN, Hdiff, Wbar, P = s.symbols(
        "Z ZN ZNN r rN rNN Hdiff Wbar P", real=True
    )
    n, v, vd, b, A, A0 = s.symbols("n v vd b A A0", real=True)

    # N^x=e*b*sin(Px)/P, v=e*v*cos(Px), W_x=e*A*sin(Px),
    # W_0=Wbar-e*A0*cos(Px). The scalar source is a one-form along dt.
    # N S_normal=r*[3(Hhat-Hclock+vdot-N^x partial_x v)-partial_i N^i].
    local_r = r + e * rN * n * cs + e**2 * rNN * n**2 * cs**2 / 2
    source_numerator = local_r * (
        3 * Hdiff + e * (3 * vd - b) * cs + 3 * e**2 * b * v * sn**2
    )
    normal_numerator = Wbar - e * A0 * cs - e**2 * b * A * sn**2 / P - source_numerator
    local_density = (Z + e * ZN * n * cs + e**2 * ZNN * n**2 * cs**2 / 2) * (
        1 + 3 * e * v * cs + s.Rational(9, 2) * e**2 * v**2 * cs**2
    )
    literal = s.expand(
        s.diff(local_density * normal_numerator**2 / 2, e, 2).subs(e, 0) / 2
    )
    # Twice the spatial average gives equal normalized real-mode cosine and sine norms.
    averaged = s.expand(literal.subs({cs**2: 1, sn**2: 1}, simultaneous=True))
    num0 = Wbar - 3 * r * Hdiff
    num1 = -A0 - r * (3 * vd - b) - 3 * rN * Hdiff * n
    num2 = (
        -b * A / P - 3 * r * v * b - rN * n * (3 * vd - b) - 3 * rNN * Hdiff * n**2 / 2
    )
    target = (
        Z * (num1**2 + 2 * num0 * num2) / 2
        + (ZN * n + 3 * Z * v) * num0 * num1
        + (ZNN * n**2 / 4 + 3 * ZN * n * v / 2 + 9 * Z * v**2 / 4) * num0**2
    )
    assert s.expand(averaged - target) == 0
    assert s.expand(target.subs(Wbar, 3 * r * Hdiff) - Z * num1**2 / 2) == 0
    # Differentiate at fixed light background before restricting to the clock.
    normal_vertex = s.diff(target, Wbar).subs({Wbar: 0, r: 0, Hdiff: 0, Z: 1})
    expected = -b * A / P - rN * n * (3 * vd - b) - (ZN * n + 3 * v) * A0
    assert s.expand(normal_vertex - expected) == 0
    assert s.diff(normal_vertex, n, vd) == -3 * rN

    binding = {
        c.Z: Z,
        c.r: r,
        c.rN: rN,
        c.dH: Hdiff,
        c.P: P,
        c.n: n,
        c.v: v,
        c.vd: vd,
        c.b: b,
        c.longitudinal: A,
        c.A0: A0,
        s.Symbol("held_temporal_vector", real=True): Wbar,
        s.Symbol("Z_N", real=True): ZN,
        s.Symbol("Z_NN", real=True): ZNN,
        s.Symbol("R_NN", real=True): rNN,
    }
    assert (
        s.expand(
            c.held_vector_data()["whole_held_temporal_vector_normal_square"].subs(
                binding, simultaneous=True
            )
            - target
        )
        == 0
    )


@pytest.mark.parametrize("case", range(6))
def test_independent_full_physical_density_three_mode_Dirac_reduction(case):
    ratio = 1 + s.Rational(case + 1, 100)
    R = ratio**4
    lapse = 1 + s.Rational(case - 2, 17)
    D, Z = ratio / lapse, ratio**-3 / lapse
    J = s.Rational(case + 3, 7)
    theta = s.Rational(case - 2, 13)
    w, charge = -s.Rational(case + 1, 31), s.Rational(case + 2, 37)
    lv, vv, vs = (
        s.Rational(case + 2, 41),
        -s.Rational(case + 1, 43),
        s.Rational(case - 2, 47),
    )
    C, Y, Yv = (
        s.Rational(case + 2, 53),
        s.Rational(case + 3, 59),
        s.Rational(case + 4, 61),
    )
    r, rn, hd = R - 1, -s.Rational(case + 1, 67), s.Rational(case - 1, 71)
    K, P, q = s.Rational(case + 2, 10**6), s.Integer(case + 2), s.Rational(case + 3, 73)
    density = s.Integer(10) ** 800 * s.Rational(case + 5, 4) ** 3
    v, sig, A, n, b, A0 = c.v, c.sigma, c.longitudinal, c.n, c.b, c.A0
    vd, sd, Ad = c.vd, c.sd, c.Ad
    lag = density * (
        -3 * D * vd**2
        + Z * sd**2 / 2
        + 6 * theta * n * vd
        + w * n * sd
        + (J - 3 * theta**2 / D + w**2 / (2 * Z)) * n**2
        - 3 * charge * vd * sig
        + lv * n * v
        + vv * v**2
        + vs * v * sig
        + b * (2 * D * vd - 2 * theta * n + charge * sig)
        + 2 * C * q * v**2
        - Y * q * sig**2 / 2
        + K * (Ad - P * A0) ** 2 / 2
        + Z * (A0 + r * (3 * vd - b) + 3 * rn * hd * n) ** 2 / 2
        - Yv * A**2 / 2
    )
    fixture = {
        c.D: D,
        c.Z: Z,
        c.J: J,
        c.K: K,
        c.P: P,
        c.th: theta,
        c.w: w,
        c.c: charge,
        c.Lnv: lv,
        c.Vvv: vv,
        c.Vvs: vs,
        c.C: C,
        c.Y: Y,
        c.Yv: Yv,
        c.r: r,
        c.rN: rn,
        c.dH: hd,
        c.q: q,
    }
    assert (
        s.expand(
            lag - density * c.complete_lagrangian().subs(fixture, simultaneous=True)
        )
        == 0
    )
    rates = s.Matrix([vd, sd, Ad])
    momenta = s.Matrix([c.pv, c.ps, c.pA])
    kinetic = s.hessian(lag, rates)
    linear = s.Matrix([s.diff(lag, x).subs(dict.fromkeys(rates, 0)) for x in rates])
    solution = kinetic.inv() * (momenta - linear)
    ham = s.expand(
        (momenta.dot(rates) - lag).subs(dict(zip(rates, solution)), simultaneous=True)
    )
    expected = density * c.unreduced_hamiltonian().subs(
        fixture, simultaneous=True
    ).subs(
        {c.pv: c.pv / density, c.ps: c.ps / density, c.pA: c.pA / density},
        simultaneous=True,
    )
    assert s.cancel(ham - expected) == 0
    pn, pb, p0 = s.symbols("independent_pn independent_pb independent_p0", real=True)
    positions = [v, sig, A, n, b, A0]
    momentums = [c.pv, c.ps, c.pA, pn, pb, p0]
    bracket = lambda f, g: s.expand(
        sum(
            s.diff(f, x) * s.diff(g, p) - s.diff(f, p) * s.diff(g, x)
            for x, p in zip(positions, momentums)
        )
    )
    constraints = [pn, pb, p0, s.diff(ham, n), s.diff(ham, b), s.diff(ham, A0)]
    dirac = s.Matrix(6, 6, lambda i, j: bracket(constraints[i], constraints[j]))
    aux = s.hessian(ham, [n, b, A0])
    gamma = 1 - 3 * Z * r**2 / (2 * D)
    assert gamma > s.Rational(1, 4)
    assert aux.det() == -4 * density**3 * D * Z * J / (3 * gamma)
    assert dirac.det() == aux.det() ** 2
    assert kinetic.det() * aux.det() == 8 * density**6 * D**2 * Z**2 * J * K
    physical_vars = [v, sig, A, c.pv, c.ps, c.pA]
    cross = s.Matrix(6, 6, lambda i, j: bracket(physical_vars[i], constraints[j]))
    assert (cross * dirac.inv() * cross.T).applyfunc(s.cancel) == s.zeros(6)
    solved = s.solve([s.diff(ham, x) for x in (n, b, A0)], [n, b, A0])
    direct_reduced = ham.subs(solved, simultaneous=True)
    expected_reduced = density * c.reduced_hamiltonian().subs(
        fixture, simultaneous=True
    ).subs(
        {c.pv: c.pv / density, c.ps: c.ps / density, c.pA: c.pA / density},
        simultaneous=True,
    )
    assert s.cancel(direct_reduced - expected_reduced) == 0


@pytest.mark.parametrize(
    "power", (s.Integer(-3), -1, s.Rational(1, 4), s.Rational(3, 4), 0, 1, 2, 3)
)
def test_entire_jet_algebra_against_direct_power_derivatives(power):
    x0 = s.Symbol("independent_positive_x0", positive=True)
    x1, x2 = s.symbols("independent_x1 independent_x2", real=True)
    t = s.Symbol("independent_epsilon", real=True)
    direct = (x0 + t * x1 + t * t * x2 / 2) ** power
    actual = physical.Jet((x0, x1, x2)) ** power
    assert all(
        s.factor(actual.values[j] - s.diff(direct, t, j).subs(t, 0)) == 0
        for j in range(3)
    )
    assert (physical.Jet.constant(0) ** 0).values == (1, 0, 0)


@pytest.mark.parametrize("case", range(4))
def test_complete_eight_mode_product_order_and_CCR(case):
    blocks = []
    for j, dim in enumerate((4, 2, 2, 2, 2, 2, 2)):
        transform = physical.symplectic_fixture(case + j, dim)
        blocks.append(transform * transform.T / 2)
    actual = physical.assemble_product(blocks)
    indices = ([0, 1, 8, 9], [2, 10], [3, 11], [4, 12], [5, 13], [6, 14], [7, 15])
    for i, index in enumerate(indices):
        assert actual.extract(index, index) == blocks[i]
        for other in indices[i + 1 :]:
            assert actual.extract(index, other) == s.zeros(len(index), len(other))
    assert (actual * physical.omega(16) * actual - physical.omega(16) / 4).applyfunc(
        s.cancel
    ) == s.zeros(16)


@pytest.mark.parametrize("case", range(3))
def test_mixed_light_vector_first_vertex_has_nonzero_noise(case):
    H = physical.vertices()["whole_canonical_Hessian"]
    A = H.diff(c.r).subs(
        {
            c.D: 1,
            c.Z: 1,
            c.J: s.Rational(case + 2, 7),
            c.K: s.Rational(1, 10**6),
            c.P: case + 2,
            c.th: s.Rational(case + 1, 11),
            c.w: -s.Rational(1, 13),
            c.c: s.Rational(1, 17),
            c.Lnv: s.Rational(case + 3, 19),
            c.Vvv: 0,
            c.Vvs: 0,
            c.C: s.Rational(1, 2),
            c.Y: 1,
            c.Yv: 1,
            c.r: 0,
            c.rN: -1,
            c.dH: 0,
            c.q: (case + 2) ** 2,
            physical.a: 1,
            physical.zeta: s.Rational(1, 10**6),
            physical.kappa: s.Integer(10) ** 800,
        },
        simultaneous=True,
    )
    result = physical.kernels(A, A, s.eye(6), s.eye(6), s.eye(6) / 2)
    assert result["mean_A"] == 0
    assert s.factor(result["symmetric_noise"]).is_positive is True
    assert A.extract([0, 1, 3, 4], [0, 1, 3, 4]) == s.zeros(4)
    assert A.extract([2, 5], [2, 5]) == s.zeros(2)


@cache
def independent_fock_phase(cutoff):
    annihilate = diags(
        np.sqrt(np.arange(1, cutoff)), 1, shape=(cutoff, cutoff), format="csr"
    )
    one = eye(cutoff, format="csr")
    modes = []
    for j in range(3):
        factors = [annihilate if i == j else one for i in range(3)]
        modes.append(
            kron(kron(factors[0], factors[1], format="csr"), factors[2], format="csr")
        )
    return [(a + a.T) / np.sqrt(2) for a in modes] + [
        (a - a.T) / (1j * np.sqrt(2)) for a in modes
    ]


def independent_weyl_operator(matrix, operators):
    coefficients = np.asarray(matrix, dtype=complex)
    return sum(
        coefficients[i, j]
        * (operators[i] @ operators[j] + operators[j] @ operators[i])
        / 4
        for i in range(6)
        for j in range(6)
    )


def independent_vertex(case, dim=6):
    return s.Matrix(
        dim,
        dim,
        lambda i, j: (
            s.Rational((-1) ** (i + j) * (case + i + j + 1), 31 + i * j + i + j)
            + (s.Rational(i + case + 1, 11) if i == j else 0)
        ),
    )


@pytest.mark.parametrize("case", range(5))
def test_whole_three_mode_Fock_means_ordered_correlations_noise_and_commutator(case):
    T = physical.symplectic_fixture(case, 6)
    L, R = (
        physical.symplectic_fixture(case + 1, 6),
        physical.symplectic_fixture(case + 2, 6),
    )
    V = T * T.T / 2
    A, B, C = [independent_vertex(case + j) for j in range(3)]
    exact = physical.kernels(A, B, L, R, V, C)
    operators = independent_fock_phase(5)
    a = independent_weyl_operator(T.T * L.T * A * L * T, operators)
    b = independent_weyl_operator(T.T * R.T * B * R * T, operators)
    c = independent_weyl_operator(T.T * L.T * C * L * T, operators)
    vacuum = np.zeros(5**3, dtype=complex)
    vacuum[0] = 1
    mean_a, mean_b = np.vdot(vacuum, a @ vacuum), np.vdot(vacuum, b @ vacuum)
    connected = np.vdot(vacuum, a @ (b @ vacuum)) - mean_a * mean_b
    commutator = -1j * np.vdot(vacuum, a @ (b @ vacuum) - b @ (a @ vacuum))
    assert abs(mean_a - complex(exact["mean_A"])) < 2e-11
    assert abs(mean_b - complex(exact["mean_B"])) < 2e-11
    assert abs(connected - complex(exact["ordered_connected_Wick"])) < 2e-10
    assert abs(connected.real - float(exact["symmetric_noise"])) < 2e-10
    assert abs(commutator - complex(exact["retarded_observable_before_step"])) < 2e-10
    assert (
        abs(-np.vdot(vacuum, c @ vacuum) - complex(exact["effective_second_contact"]))
        < 2e-11
    )


@pytest.mark.parametrize("case", range(3))
def test_whole_eight_mode_explicit_four_index_contraction(case):
    dim = 16
    T = physical.symplectic_fixture(case, dim)
    L, R = (
        physical.symplectic_fixture(case + 1, dim),
        physical.symplectic_fixture(case + 2, dim),
    )
    A, B = s.eye(dim), 2 * s.eye(dim)
    A[0, 7] = A[7, 0] = s.Rational(case + 1, 13)
    B[1, 15] = B[15, 1] = s.Rational(case + 2, 17)
    result = physical.kernels(A, B, L, R, T * T.T / 2)
    W = result["whole_unequal_time_W"]
    aa = [(i, j, A[i, j]) for i in range(dim) for j in range(dim) if A[i, j] != 0]
    bb = [(k, l, B[k, l]) for k in range(dim) for l in range(dim) if B[k, l] != 0]
    literal = sum(
        v * w * (W[i, k] * W[j, l] + W[i, l] * W[j, k]) / 4
        for i, j, v in aa
        for k, l, w in bb
    )
    assert s.expand(literal - result["ordered_connected_Wick"]) == 0


@pytest.mark.parametrize("case", range(3))
def test_whole_three_mode_CTP_overlap_at_two_occupation_cutoffs(case):
    O = np.asarray(physical.omega(6), dtype=float)
    frequencies = np.array([1.0, 1.3, 1.7])
    H0 = np.diag(np.r_[frequencies, frequencies])
    A = np.asarray(independent_vertex(case), dtype=float) * 0.035
    B = np.asarray(independent_vertex(case + 1), dtype=float) * 0.02
    Hplus, Hminus = H0 + A + B / 2, H0 - A + B / 2
    duration = 0.08 + case * 0.01
    relative = expm(-duration * O @ Hminus) @ expm(duration * O @ Hplus)
    alpha = (
        relative[:3, :3] + relative[3:, 3:] + 1j * (relative[3:, :3] - relative[:3, 3:])
    ) / 2
    # Short continuously connected path, with no logarithm-cut crossing.
    expected = np.exp(-np.log(np.linalg.det(alpha.conjugate())) / 2)
    overlaps = []
    for cutoff in (5, 7):
        ops = independent_fock_phase(cutoff)
        plus, minus = (
            independent_weyl_operator(Hplus, ops),
            independent_weyl_operator(Hminus, ops),
        )
        vacuum = np.zeros(cutoff**3, dtype=complex)
        vacuum[0] = 1
        left = expm_multiply(-1j * duration * minus, vacuum)
        right = expm_multiply(-1j * duration * plus, vacuum)
        overlaps.append(np.vdot(left, right))
    assert abs(overlaps[1] - expected) < 2e-11
    assert abs(overlaps[0] - overlaps[1]) < 2e-9


@pytest.mark.parametrize("case", range(3))
def test_whole_three_mode_time_Duhamel_and_second_contact(case):
    O = np.asarray(physical.omega(6), dtype=float)
    A = np.asarray(independent_vertex(case), dtype=float) / 7
    C = np.asarray(independent_vertex(case + 1), dtype=float) / 11
    readout = np.asarray(independent_vertex(case + 2), dtype=float)
    V = np.eye(6) / 2
    duration = 0.27 + case / 100
    H0 = lambda t: np.diag(
        [1 + t / 11, 1.2 + t / 13, 1.4 + t / 17, 1 + t / 19, 1.1 + t / 23, 1.3 + t / 29]
    )
    f = lambda t: np.sin(2 * t) + 0.3
    g = lambda t: 0.5 + t * t

    def rhs(t, y):
        base, first, second = y.reshape(3, 6, 6)
        G, GA, GC = O @ H0(t), f(t) * O @ A, g(t) * O @ C
        return np.array(
            [G @ base, G @ first + GA @ base, G @ second + 2 * GA @ first + GC @ base]
        ).ravel()

    initial = np.array([np.eye(6), np.zeros((6, 6)), np.zeros((6, 6))]).ravel()
    result = solve_ivp(
        rhs, [0, duration], initial, rtol=2e-12, atol=2e-14, dense_output=True
    )
    assert result.success
    S0, S1, S2 = result.y[:, -1].reshape(3, 6, 6)

    def solve(epsilon):
        fun = lambda t, y: (
            O
            @ (H0(t) + epsilon * f(t) * A + epsilon**2 * g(t) * C / 2)
            @ y.reshape(6, 6)
        ).ravel()
        value = solve_ivp(fun, [0, duration], np.eye(6).ravel(), rtol=1e-12, atol=1e-14)
        assert value.success
        return value.y[:, -1].reshape(6, 6)

    epsilon = 0.001
    plus, minus = solve(epsilon), solve(-epsilon)
    assert np.max(np.abs((plus - minus) / (2 * epsilon) - S1)) < 2e-8
    assert np.max(np.abs((plus - 2 * S0 + minus) / epsilon**2 - S2)) < 2e-7
    target = np.trace(readout @ (S1 @ V @ S0.T + S0 @ V @ S1.T)) / 2
    pulled_readout = S0.T @ readout @ S0

    def commutator_integrand(t):
        S = result.sol(t).reshape(3, 6, 6)[0]
        pulled = S.T @ (f(t) * A) @ S
        return (
            np.trace((pulled_readout @ O @ pulled - pulled @ O @ pulled_readout) @ V)
            / 2
        )

    integral = quad(commutator_integrand, 0, duration, epsabs=1e-12, epsrel=1e-12)[0]
    assert abs(integral - target) < 2e-10
    contact = g(duration) * np.trace(C @ S0 @ V @ S0.T) / 2
    response = f(duration) * np.trace(A @ (S1 @ V @ S0.T + S0 @ V @ S1.T)) / 2
    current = lambda eps, S: (
        np.trace((f(duration) * A + eps * g(duration) * C) @ S @ V @ S.T) / 2
    )
    finite = (current(epsilon, plus) - current(-epsilon, minus)) / (2 * epsilon)
    assert abs(finite - contact - response) < 2e-8
    assert abs(contact) > 1e-4
