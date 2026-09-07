"""Independent Fraction/Taylor polarization and polynomial replay.

No SymPy or primary-module imports. Only the frozen S6.20 arithmetic
primitive/family and older sparse polynomial primitive are reused; no old
constraint, scalar-health, vector-speed or certificate verdict is imported.
Quadratic Hessians are recovered by polarization, not symbolic derivatives.
"""

from fractions import Fraction as Q
from itertools import product
from math import comb

from p8_composite_modes.independent import NAMES, Poly
from p8_variable_beta.independent import Jet, derivative, exact, family


def transpose(matrix):
    return [list(row) for row in zip(*matrix, strict=True)]


def matmul(a, b):
    return [[sum((x*y for x, y in zip(row, column, strict=True)), Q())
             for column in zip(*b, strict=True)] for row in a]


def inverse(matrix):
    n = len(matrix)
    if not n or any(len(row) != n for row in matrix):
        raise ValueError("An exact nonempty square matrix is required")
    rows = [[exact(x) for x in row]+[Q(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            raise ValueError("A singular matrix has no inverse")
        rows[j], rows[pivot] = rows[pivot], rows[j]
        factor = rows[j][j]
        rows[j] = [x/factor for x in rows[j]]
        for i in range(n):
            if i != j:
                coefficient = rows[i][j]
                rows[i] = [x-coefficient*y for x, y in zip(rows[i], rows[j], strict=True)]
    return [row[n:] for row in rows]


def hessian_by_polarization(function, dimension):
    basis = [[Q(i == j) for j in range(dimension)] for i in range(dimension)]
    diagonals = [function(row) for row in basis]
    result = [[None]*dimension for _ in range(dimension)]
    for i in range(dimension):
        result[i][i] = 2*diagonals[i]
        for j in range(i):
            value = function([x+y for x, y in zip(basis[i], basis[j], strict=True)])-diagonals[i]-diagonals[j]
            result[i][j] = result[j][i] = value
    return result


def scalar_jet_fixture(time, lapse, momentum_squared):
    """Independent complete regular Hamiltonian value/first-jet fixture."""
    time, c, K = map(exact, (time, lapse, momentum_squared))
    if not (abs(time) <= Q(1, 10) and (c == 1 or 2 < c <= 4) and K > 0):
        raise ValueError("Outside the nonzero-mode frozen local family")
    f = family(Jet(time, 1), c)
    a, b, y, h, hp, w, n, kin, P = (f[name] for name in
                                    ("a", "b", "y", "h", "h_u", "chi_speed", "nbar", "kbar", "P"))
    C = a**5*y*y*P/(2*(c+y))
    r, alpha = 1+2*b*K/(3*a**3*P), h/(3*c*a**3*P)
    # Only value and first derivative of these differentiated jets are used.
    rp = Jet(r.first, 2*r.second_coefficient)
    ap = Jet(alpha.first, 2*alpha.second_coefficient)

    def transformed(values):
        R, e, z, pr, pe, px, p0 = values
        qg, qf, pg, pf, old_pe = r*R+e/3+alpha*pr, R+alpha*p0, p0, pr-r*p0, pe-p0/3
        sg, sf = 3*qg, 3*qf+e
        peg = a**3*w*z-old_pe
        clock = -h*pg+a**3*n*sg+a**3*y*P*(sg-sf)+2*a*K*qg
        H = (3*peg*peg-2*pg*peg)/(4*a**3)+c*(3*old_pe**2-2*pf*old_pe)/(4*b**3)
        H += K*old_pe**2/(4*C)+(clock-w*px)**2/(2*a**3*kin)+px**2/(2*a**3)
        H += -sg*clock+a**3*n*sg**2/4+a*K*z*z/2
        H += -a**3*(hp+Q(3, 2)*h*h)*sg**2-b**3/c*(-hp+Q(3, 2)*h*h)*sf**2
        H += -a*K*qg*qg-c*b*K*qf*qf
        H += a**3*(f["b0"]+c*f["b1"])*sg**2
        H += a**3*y*f["b1"]*(2*(2*qg+qf)**2+(2*qg+qf+e)**2)+c*b**3*f["b4"]*sf**2
        H += -rp*p0*R-ap*p0*pr+(r*ap-alpha*rp)*p0*p0/2
        return H

    matrix = hessian_by_polarization(transformed, 7)
    D = matrix[6][6]
    p0 = [-matrix[6][i]/D for i in range(6)]
    reduced = [[matrix[i][j]+matrix[i][6]*p0[j] for j in range(6)] for i in range(6)]
    return {"u": time, "c": c, "K": K,
            "M0": [[x.value for x in row] for row in matrix],
            "M1": [[x.first for x in row] for row in matrix],
            "H0": [[x.value for x in row] for row in reduced],
            "H1": [[x.first for x in row] for row in reduced],
            "p0": [[x.value for x in p0], [x.first for x in p0]],
            "D": [D.value, D.first]}


def scalar_fixtures():
    cases = [(0, 4, 1), (0, 3, Q(7, 2)), (0, 1, 5),
             (Q(1, 100), 4, 7), (Q(-1, 100), 4, 7),
             (Q(1, 10), Q(201, 100), Q(3, 2))]
    return [scalar_jet_fixture(*case) for case in cases]


def vector_fixtures():
    result = []
    for u, c, K in ((0, 4, 1), (Q(1, 10), 3, Q(7, 2)),
                    (Q(-1, 10), Q(201, 100), 5), (0, 1, 200)):
        c, K = exact(c), exact(K)
        f = family(Jet(exact(u), 1), c)
        a, b, y, P = (f[name] for name in ("a", "b", "y", "P"))
        A, B, C, D = a**3/4, b**3/(4*c), a**5*y*y*P/(2*(c+y)), a**3*y*P/4
        # Solve the literal two-shift matrix, without the advertised Schur formula.
        shift_hessian = [[(A.value*K+C.value), -C.value], [-C.value, (B.value*K+C.value)]]
        source = [[A.value*K], [Q()]]  # E_g'=1,E_f'=0.
        shifts = matmul(inverse(shift_hessian), source)
        sg, sf = shifts[0][0], shifts[1][0]
        literal = A.value*K*(1-sg)**2+B.value*K*sf**2+C.value*(sf-sg)**2
        kinetic = A*B*C*K/(A*B*K+(A+B)*C)
        if literal != kinetic.value:
            raise ValueError("Independent vector stationarity disagreed with harmonic formula")
        theta = kinetic.first/(2*kinetic.value)
        pump = kinetic.second_coefficient/kinetic.value-theta**2
        result.append({"u": exact(u), "c": c, "K": K, "kinetic": literal,
                       "bare": D.value*K/literal, "pump": pump,
                       "speed_squared": (c+y.value)/(2*y.value),
                       "mass_squared": y.value*P.value*(1+c/y.value**3)})
    return result


def margin_polynomial():
    """Derive the cleared D margin by polynomial arithmetic in w=u²,c."""
    w, c = Poly.variable("x"), Poly.variable("c")
    d = 1+w
    F = 6400*(1-w)-800*c*(1-w)*d**12-c*d**2
    t = d**3*(c*d**4-2)
    # alpha=u*t/[48(1-w)]; differentiate before clearing denominators.
    alpha_numerator = (1-w)*t+2*w*((1-w)*derivative(t, "x")+t)
    return (c*(1-w)**2*F*(8+c*d**12)*Q(1, 2)
            -38400*w*(1-w)**2*d**4*(c*d**4+2)**2
            +8*w*(1-w)*d**4*F*(c*d**4-2)
            -c*d**6*F*alpha_numerator*Q(1, 2)-3*c*(1-w)**2*d**6*F)


def margin_coefficients():
    """Independent affine substitution and binomial Bernstein conversion."""
    raw = margin_polynomial()
    iw, ic = NAMES.index("x"), NAMES.index("c")
    result = {}
    for name, origin, width in (("c1", Q(1), Q()), ("positive", Q(2), Q(2))):
        terms = {}
        for power, coefficient in raw.terms.items():
            degree_w, degree_c = power[iw], power[ic]
            for j in range(degree_c+1):
                value = coefficient*Q(1, 100)**degree_w*comb(degree_c, j)*origin**(degree_c-j)*width**j
                terms[degree_w, j] = terms.get((degree_w, j), Q())+value
        terms = {power: value for power, value in terms.items() if value}
        nw, nc = max(p[0] for p in terms), max(p[1] for p in terms)
        values = []
        for i, j in product(range(nw+1), range(nc+1)):
            values.append(sum((value*Q(comb(i, p), comb(nw, p))*Q(comb(j, q), comb(nc, q))
                               for (p, q), value in terms.items() if p <= i and q <= j), Q()))
        if any(value <= 0 for value in values):
            raise ValueError("Independent continuous constraint margin lost positivity")
        result[name] = tuple(values)
    return result


def checks():
    scalar, vector, margin = scalar_fixtures(), vector_fixtures(), margin_coefficients()
    if any(item["D"][0] >= -Q(1, 8) for item in scalar):
        raise ValueError("An actual scalar fixture violated the proved denominator bound")
    return {"scalar_phase_jet_fixtures": [{key: f[key] for key in ("u", "c", "K", "D")} for f in scalar],
            "vector_literal_shift_jet_fixtures": vector,
            "independent_Bernstein_counts": {name: len(values) for name, values in margin.items()},
            "independent_Bernstein_minima": {name: min(values) for name, values in margin.items()}}
