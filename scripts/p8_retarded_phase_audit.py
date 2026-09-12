"""Independent finite Fock and real-covariance orientation diagnostics.

No frozen endpoint helper is imported. Quadratic products acting on the vacuum
need at most two quanta, so the three-level fixture gives the exact CCR result.
"""

import sympy as s

ann = s.zeros(3)
ann[0, 1] = 1
ann[1, 2] = s.sqrt(2)
cre = ann.T
Q = (ann + cre) / s.sqrt(2)
P = -s.I * (ann - cre) / s.sqrt(2)
J = s.Matrix([[0, 1], [-1, 0]])
C = s.eye(2) / 2
U = s.Matrix(
    [[s.Rational(3, 5), s.Rational(4, 5)], [-s.Rational(4, 5), s.Rational(3, 5)]]
)
z = [U[0, 0] * Q + U[0, 1] * P, U[1, 0] * Q + U[1, 1] * P]


def observable(M, fields):
    return sum(
        (M[i, j] * fields[i] * fields[j] / 2 for i in range(2) for j in range(2)),
        s.zeros(3),
    )


bases = (s.diag(1, 0), s.diag(0, 1), s.Matrix([[0, 1], [1, 0]]))
for i, D in enumerate(bases):
    for j, G in enumerate(bases):
        hd, hg = observable(D, z), observable(G, [Q, P])
        literal = s.expand(s.I * (hd * hg - hg * hd)[0, 0])
        covariance = -s.trace(D * U * (J * G * C + C * G * J.T) * U.T) / 2
        u = s.Matrix([1, -s.I]) / s.sqrt(2)
        det_ann = (u.T * U.T * D * U * u)[0]
        src_ann = (u.T * G * u)[0]
        annih_pair = s.im(s.conjugate(det_ann) * src_ann).expand(complex=True)
        assert s.simplify(literal - covariance) == 0
        assert s.simplify(literal - annih_pair) == 0
        assert s.simplify(literal + s.im(det_ann * s.conjugate(src_ann))) == 0
print(
    "Exact finite Fock / covariance / annihilation orientation: 9 cases PASS",
    flush=True,
)

dr, di, gr, gi, theta = s.symbols("dr di gr gi theta", real=True)
ad, ag = dr + s.I * di, gr + s.I * gi
physical = s.expand(s.im(s.conjugate(ad) * ag * (s.cos(theta) + s.I * s.sin(theta))))
frozen = s.expand(-s.im(s.conjugate(ad) * ag * (s.cos(theta) - s.I * s.sin(theta))))
assert s.expand(physical - frozen - 2 * (dr * gi - di * gr) * s.cos(theta)) == 0
assert (physical - frozen).subs({dr: 1, di: 0, gr: 0, gi: 1, theta: 0}) == 2
print("Wrong-phase discrepancy =", s.factor(physical - frozen), flush=True)

GJ = s.symbols("g0:8", real=True)
AJ = s.symbols("a0:8", real=True)
HJ = s.symbols("h0:8", real=True)


def derivative(expr):
    return s.expand(
        sum(s.diff(expr, seq[j]) * seq[j + 1] for seq in (GJ, AJ, HJ) for j in range(7))
    )


rows = [AJ[0] * HJ[0]]
for _ in range(6):
    rows.append(derivative(GJ[0] * rows[-1]))
for n in range(1, 7):
    plus = sum(s.I * (-s.I) ** j * GJ[0] * rows[j] for j in range(n))
    assert (
        s.expand(
            derivative(plus) - s.I * plus / GJ[0] - rows[0] + (-s.I) ** n * rows[n]
        )
        == 0
    )
print("Positive detector-sharp phase: all 6 retarded ODE identities PASS", flush=True)
zr, zi = s.symbols("zr zi", real=True)
for j in range(6):
    old = -s.im(-s.I * s.I**j * (zr + s.I * zi))
    new = s.im(s.I * (-s.I) ** j * (zr + s.I * zi))
    assert s.expand(new - (-1) ** j * old) == 0
print("Correct current multiplies frozen endpoint j by (-1)^j: all 6 PASS", flush=True)


def dimensional_regression():
    """Optional exact check of the frozen tracefree evanescent error."""
    import sys
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    sys.setrecursionlimit(4000)
    sys.set_int_max_str_digits(0)
    for source in sorted((repo / "problems/P8").rglob("src")):
        sys.path.insert(0, str(source))
    sys.path.insert(0, str(repo / "problems/P8/s6/continuation/s6_160/src"))
    from p8_vacuum_affine_dimensional_spatial_symbol import density, geometry, jets
    from sympy.core.random import seed

    seed(0)
    rows = density.invariant_spatial_coefficients()
    odd = {
        key: row for key, row in rows.items() if key[0] % 2 and any(v != 0 for v in row)
    }
    assert set(odd) == {(1, 0, 3)}
    wanted = (-s.Rational(31, 420), s.Rational(1, 14), s.Integer(0))
    for value, coefficient in zip(odd[1, 0, 3], wanted):
        assert s.factor(value.subs(geometry.d, 3)) == 0
        assert (
            s.factor(
                s.diff(value, geometry.d).subs(geometry.d, 3)
                - coefficient * jets.p**2 * (3 * jets.t**2 + 1)
            )
            == 0
        )
    print("Exact tracefree evanescent phase correction: PASS", flush=True)
    print(
        "Corrected minus frozen Ffinite = p^2(3t^2+1)(-31T+30V)Gamma0/(420pi^2)",
        flush=True,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dimensional", action="store_true")
    if parser.parse_args().dimensional:
        dimensional_regression()
