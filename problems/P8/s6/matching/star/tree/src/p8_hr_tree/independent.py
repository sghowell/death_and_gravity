"""Independent Fraction/first-jet engine: no SymPy or new primary imports.

Only immutable Poly and Dual arithmetic primitives are reused. Literal
edge stress derivatives and incidence Gaussian elimination are independent
of the primary symbolic and leaf-stripping implementations.
"""

from fractions import Fraction as Q

from p8_composite_modes.independent import NAMES, Poly
from p8_trimetric.independent import Dual


def edge_lagrangian(beta, ni, ai, nj, aj):
    b0, b1, b2, b3, b4 = beta
    return -2*(b0*ni*ai*ai*ai+b1*(nj*ai*ai*ai+3*ni*ai*ai*aj)
               +3*b2*(nj*ai*ai*aj+ni*ai*aj*aj)
               +b3*(3*nj*ai*aj*aj+ni*aj*aj*aj)+b4*nj*aj*aj*aj)


def polynomial_derivative(polynomial, name):
    """Literal coefficient differentiation in the independent sparse ring."""
    index = NAMES.index(name)
    result = {}
    for monomial, coefficient in polynomial.terms.items():
        power = monomial[index]
        if power:
            reduced = list(monomial)
            reduced[index] -= 1
            result[tuple(reduced)] = coefficient*power
    return Poly(result)


def stress_fixtures():
    cases = [((-3, 1, 0, 0, -1), Q(1), Q(1)),
             ((1, -2, Q(3, 2), -4, 5), Q(2), Q(3, 2)),
             ((0, 1, -1, 1, 0), Q(1), Q(2)),
             ((0, 1, -Q(1, 2), 0, Q(1, 2)), Q(1), Q(3)),
             ((2, 0, 0, 0, -3), Q(3, 2), Q(2, 3))]
    result = []
    for beta, y, c in cases:
        ni, ai, nj, aj = Q(2), Q(3), 2*c, 3*y
        vals = (ni, ai, nj, aj)
        gradients = []
        for index in range(4):
            jets = tuple(Dual(value, Q(i == index)) for i, value in enumerate(vals))
            gradients.append(edge_lagrangian(beta, *jets).tangent)
        ri, pi = -gradients[0]/ai**3, gradients[1]/(3*ni*ai**2)
        rj, pj = -gradients[2]/aj**3, gradients[3]/(3*nj*aj**2)
        poly = 2*(beta[1]+2*beta[2]*y+beta[3]*y*y)
        if ri+pi != (y-c)*poly or rj+pj != (c-y)*poly/(c*y**3):
            raise ValueError("Independent literal edge first-jet stresses failed")
        result.append({"beta": tuple(map(Q, beta)), "y": y, "c": c, "P": poly,
                       "rho_i": ri, "pressure_i": pi, "rho_j": rj, "pressure_j": pj,
                       "null_i": ri+pi, "null_j": rj+pj})
    return result


def identities():
    y, c, p, hi, hj = (Poly.variable(name) for name in ("y", "c", "p", "x", "z"))
    # Explicit cleared derivative of the cubic lapse density; rho_i'=3P.
    ci = 3*p*y*(c*hj-hi)+3*hi*(y-c)*p
    cj_cleared = -3*p*c*(c*hj-hi)+3*hj*c*(c-y)*p
    # cj_cleared=c² y³ C_j, from derivative of the inverse-ratio density.
    g0, g1, r, rp, h, hp = (Poly.variable(name) for name in ("A", "B", "y", "u", "x", "v"))
    a, ap = g0+g1*r*r, 2*g1*r*rp
    full = -2*g0*hp-2*g1*(r*r*hp-r*h*rp)
    # Derive both cleared null stresses from independent four-coframe
    # elementary-action differentiation, not from their claimed factors.
    ni, ai, nj, aj = (Poly.variable(name) for name in ("A", "B", "C", "z"))
    beta = tuple(Poly.variable(f"beta{i}") for i in range(5))
    lagrangian = edge_lagrangian(beta, ni, ai, nj, aj)
    ni_null = -3*ni*polynomial_derivative(lagrangian, "A")+ai*polynomial_derivative(lagrangian, "B")
    nj_null = -3*nj*polynomial_derivative(lagrangian, "C")+aj*polynomial_derivative(lagrangian, "z")
    expected_null = 6*(aj*ni-nj*ai)*(beta[1]*ai*ai+2*beta[2]*ai*aj+beta[3]*aj*aj)
    # ni_null = 3 Ni ai³ I_i and nj_null = 3 Nj aj³ I_j.
    return {"literal_i_Bianchi_cleared": ci-3*p*c*(y*hj-hi),
            "reciprocal_Bianchi_two_lapses": ci+cj_cleared,
            "literal_i_null_coframe_polynomial": ni_null-expected_null,
            "literal_j_null_coframe_polynomial": nj_null+expected_null,
            "reciprocal_null_one_lapse": ni_null+nj_null,
            "conditional_two_vertex_weighted_null": full-(-2*a*hp+h*ap)}


def incidence_fixtures():
    """Independent rational row reduction, not the primary pruning route."""
    cases = [(2, ((0, 1),)), (4, ((0, 1), (1, 2), (2, 3))),
             (4, ((1, 0), (0, 2), (3, 0))),
             (6, ((0, 1), (2, 1), (2, 3), (4, 2), (4, 5)))]
    result = []
    for n, edges in cases:
        columns = n-1
        original = [[Q((i == left)-(i == right)) for left, right in edges] for i in range(n)]
        matrix = [row[:] for row in original]
        rank = 0
        for column in range(columns):
            pivot = next(i for i in range(rank, n) if matrix[i][column])
            matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
            divisor = matrix[rank][column]
            matrix[rank] = [entry/divisor for entry in matrix[rank]]
            for i in range(n):
                if i != rank:
                    coefficient = matrix[i][column]
                    matrix[i] = [x-coefficient*y for x, y in zip(matrix[i], matrix[rank], strict=True)]
            rank += 1
        if rank != columns or any(matrix[-1]):
            raise ValueError("Independent finite-tree full column rank failed")
        fluxes = tuple(Q((-1)**i*(i+1), i+2) for i in range(columns))
        rhs = tuple(sum(row[i]*fluxes[i] for i in range(columns)) for row in original)
        result.append({"vertices": n, "edges": edges, "rank": rank, "fluxes": fluxes, "divergences": rhs})
    return result


def component_fixtures():
    cases = [((1,), (1,), (1,), (2,), Q(0), (0,), 0),
             ((1, 0, 2), (1, 2, 3), (1, 3, 2), (1, 0, 2), Q(-2, 3), (0, 1, 2), 0),
             ((1, 1, 1), (1, 1, 1), (1, Q(1, 3), 1), (2, 0, 2), Q(1, 3), (0, 2), 0),
             ((0, 3, 2), (2, 1, Q(3, 2)), (3, 1, 2), (0, 1, 2), Q(5, 7), (0, 1, 2), 1),
             ((1, 0, 7), (1, 2, 9), (1, 1, 5), (3, 0, 4), Q(2), (0,), 0)]
    result = []
    for gs, ys, cs, ns, h, component, root in cases:
        gs, ys, cs, ns = (tuple(map(Q, values)) for values in (gs, ys, cs, ns))
        ajet = Dual(Q(0))
        yp = {}
        for i in component:
            # The kinematic scale ratio is differentiated first in the Dual ring.
            yp[i] = (cs[i]-ys[i])*h
            ratio = Dual(ys[i], yp[i])
            ajet += gs[i]*ratio*ratio
        source = sum(cs[i]*ys[i]**3*ns[i] for i in component)
        lam, eta = ajet.tangent/(2*ajet.value), source/(2*ajet.value)
        hp = lam*h-eta
        full = -2*sum(gs[i]*cs[i]*ys[i]**3*((Dual(h, hp)/Dual(ys[i], yp[i])).tangent/cs[i]) for i in component)
        rate_bound = max(abs(yp[i]/ys[i]) for i in component)
        if full != source or ajet.value < gs[root] or gs[root] <= 0 or abs(lam) > rate_bound:
            raise ValueError("Independent component clock/source/denominator replay failed")
        result.append({"G": gs, "y": ys, "c": cs, "nulls": ns, "H": h, "component": component, "root": root,
                       "A": ajet.value, "A_fixed_prime": ajet.tangent, "weighted_source": source,
                       "lambda": lam, "eta": eta, "Hprime": hp, "component_log_rate_bound": rate_bound,
                       "yprimes": yp})
    return result


def checks():
    polys = identities()
    if any(not value.is_zero() for value in polys.values()):
        raise ValueError("Independent coefficientwise HR-tree identity failed")
    return {"coefficientwise_identities": dict.fromkeys(polys, "0"),
            "literal_edge_stress_fixtures": stress_fixtures(),
            "Gaussian_incidence_fixtures": incidence_fixtures(),
            "component_first_jet_fixtures": component_fixtures()}
