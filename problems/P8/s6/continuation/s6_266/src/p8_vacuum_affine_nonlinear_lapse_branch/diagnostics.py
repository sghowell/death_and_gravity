"""Independent modest-parameter complete-source root fixtures, not a retuning."""

from functools import cache

import numpy as np
import sympy as s
from scipy.integrate import quad
from scipy.optimize import brentq

from . import source as q


@cache
def kernels():
    N = q.N
    expressions = {
        "C": q.CONSTRAINT,
        "C_N": s.diff(q.CONSTRAINT, N),
        "C_NN": s.diff(q.CONSTRAINT, N, 2),
        "T": q.TEMPORAL,
        "T_N": s.diff(q.TEMPORAL, N),
        "normal": q.NORMAL,
        "primitive_N": q.PRIMITIVE_N,
    }
    atoms = sorted(
        set().union(
            *(
                expr.atoms(s.Derivative) | expr.atoms(s.Function)
                for expr in expressions.values()
            )
        ),
        key=str,
    )
    symbols = s.symbols("full_numeric_source_jet0:" + str(len(atoms)), real=True)
    mapping = dict(zip(atoms, symbols, strict=True))
    inputs = (N, *q.COORDS, q.mu, *symbols)
    lifted = [expr.xreplace(mapping) for expr in expressions.values()]
    assert all(expr.free_symbols <= set(inputs) for expr in lifted)
    return atoms, tuple(expressions), s.lambdify(inputs, lifted, "numpy", cse=True)


@cache
def fixture(case):
    if not isinstance(case, int) or isinstance(case, bool) or not 0 <= case < 5:
        raise ValueError("Use one of five explicitly finite diagnostics")
    N = q.N
    X = N**-2
    order = s.Integer(8 + 2 * case)
    localizer = s.Integer(17 + case)
    T = X**order / (X**order + (1 - X) ** order)
    bump = order * X**2 * s.exp(-order * X**2)
    B = T + (1 - T) * bump
    heavy = X**2 * (1 - X) ** 8 * s.exp(-localizer * X**2)
    R = 1 + B * (X - 1) + order * heavy
    Ruu = -6 * B * (X - 1)
    # Nonzero profile and vacuum/heavy terms are diagnostic coefficients,
    # not assertions that a changed fixture is the actual P8 parent.
    F = (
        q.F_TREE
        + (case + 1) * s.Rational(1, 100000) * T * (1 + 2 * (X - 1))
        + (1 - T) * (1 + X * X)
        + heavy / 7
    )
    j = (
        s.Integer(10) ** 100
        * s.Rational(case + 1, 100)
        * ((X - 1) ** 3 * s.exp(-localizer * X * X) + T / s.Integer(100))
    )
    bind = {}
    for f, value in [(q.R, R), (q.F, F), (q.Ruu, Ruu), (q.j, j)]:
        for k in range(4):
            bind[s.diff(f, N, k)] = s.diff(value, N, k)
    mass = s.Rational(2 + case, 1000)
    rule = {**bind, q.mu: mass}
    atoms, names, kernel = kernels()
    coefficients = [bind[atom] for atom in atoms]
    assert all(expr.free_symbols <= {N} for expr in coefficients)
    source_kernel = s.lambdify(N, coefficients + [R, F, Ruu, j], "numpy", cse=True)

    def evaluated(n, values):
        jets = source_kernel(n)
        return kernel(n, *values, float(mass), *jets[: len(atoms)])

    def select(index):
        return lambda n, *values: evaluated(n, values)[index]

    functions = {name: select(index) for index, name in enumerate(names)}
    functions["source"] = lambda n: source_kernel(n)[len(atoms) :]
    functions["rule"] = rule
    return functions


def vector(case, scale=1e-5):
    rng = np.random.default_rng(820266 + case)
    values = rng.uniform(-scale, scale, size=12)
    # Physical invariant subset for every nominal root diagnostic.
    values[5:11] = np.abs(values[5:11])
    return values


def root(case, values):
    f = fixture(case)["C"]
    return brentq(lambda n: float(f(n, *values)), 0.999, 1.001, xtol=5e-15, rtol=1e-14)


def primitive(case, n):
    f = fixture(case)["primitive_N"]
    zero = np.zeros(12)
    return quad(lambda x: float(f(x, *zero)), 1, float(n), epsabs=2e-14, epsrel=2e-13)[
        0
    ]


def independent_hamiltonian(case, n, values):
    data = fixture(case)
    p, G, dp, ph, eta, sh, el, ma, wm, gm, gh, curv = values
    R, F, _, j = data["source"](n)
    mass = (2 + case) / 1000
    # Direct independent numeric transcription keeps every original channel.
    normal = (
        -3 * (p - (R - 1) * G) ** 2 / (4 * R**0.25)
        - R ** (-0.75) * F
        + G**2 * R**0.75 / 2
        + 2 * sh * R ** (-0.25)
        + ((0.1 + dp) ** 2 + ph**2) * R**0.75 / 2
        + R ** (-0.25) * (gm + gh) / 2
        + R ** (-0.75) * (mass * eta**2 / 2 - j * eta / 10.0**100)
        - R**0.75 * curv / 2
        + el * R**0.25 / 2
        + ma * R**0.25 / 4
        + R ** (-0.25) * wm / 2
    )
    return n * normal + primitive(case, n)
