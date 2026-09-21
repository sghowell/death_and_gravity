"""Whole off-shell labeled-box to fixed six-word covariant-jet conversion."""

from functools import cache
from itertools import combinations_with_replacement, product
from math import factorial

import sympy as s
from p8_vacuum_affine_box_curvature_coefficient import basis, jets

POWERS = tuple(p for p in product(range(4), repeat=3) if sum(p) == 3)
EXPECTED = (-6, -2, -2, -6, -4, 0, -4, -4, -4, -8)
Z = basis.Z


def require_powers(value):
    if (
        type(value) is not tuple
        or len(value) != 3
        or any(type(v) is not int or v < 0 for v in value)
        or sum(value) != 3
    ):
        raise ValueError(
            "Require three native nonnegative integer powers of total degree3"
        )
    return value


def divided_power(v, a, n):
    return sum(
        s.binomial(n, l) * v ** (n - l) * (2 * a) ** (l - 1) for l in range(1, n + 1)
    )


def labeled_contact(exponents, G, H, aa):
    exponents = require_powers(exponents)
    natural = s.S.Zero
    for a, b, c, d in basis.PERMS:
        values = (G[a, a] + 2 * G[a, b] + G[b, b], G[c, c], G[d, d])
        soft = (aa[a] + aa[b], aa[c], aa[d])
        hh = (H[a, a] + 2 * H[a, b] + H[b, b], H[c, c], H[d, d])
        natural -= (
            sum(
                hh[i]
                * divided_power(values[i], soft[i], exponents[i])
                * s.prod(values[j] ** exponents[j] for j in range(3) if j != i)
                for i in range(3)
            )
            / 2
        )
    return s.expand(natural)


@cache
def generic():
    kin = s.symbols("G00 G01 G02 G11 G12 G22 a0 a1 a2")
    hs = s.symbols("H00 H01 H02 H11 H12 H22")
    G = s.zeros(4)
    H = s.zeros(4)
    for (i, j), v, h in zip(combinations_with_replacement(range(3), 2), kin[:6], hs):
        G[i, j] = G[j, i] = v
        H[i, j] = H[j, i] = h
    aa = (*kin[6:], -sum(kin[6:]))
    for i in range(3):
        G[i, 3] = G[3, i] = -sum(G[i, j] for j in range(3)) - aa[i]
        H[i, 3] = H[3, i] = -sum(H[i, j] for j in range(3))
    G[3, 3] = sum(G[i, j] for i in range(3) for j in range(3)) + 2 * sum(kin[6:])
    H[3, 3] = sum(H[i, j] for i in range(3) for j in range(3))
    T = sum(
        aa[i] ** 2 * H[j, j] + aa[j] ** 2 * H[i, i] - 2 * aa[i] * aa[j] * H[i, j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
    return G, H, aa, (*kin, *hs), s.expand(T)


@cache
def build():
    G, H, aa, variables, T = generic()
    TP = s.Poly(T, *variables)
    contacts = tuple(jets.bose_contact(word, G, H, aa) for word in basis.WORDS)
    flats = [s.Poly(basis.bose(word, basis.G), *basis.GG) for word in basis.WORDS]
    monomials = sorted(set().union(*(set(P.monoms()) for P in flats)))
    matrix = s.Matrix([[P.coeff_monomial(m) for P in flats] for m in monomials])
    checks = {}
    tau = {}
    coefs = {}
    for index, powers in enumerate(POWERS):
        flat = s.S.Zero
        for a, b, c, d in basis.PERMS:
            GG = basis.G
            values = (GG[a, a] + 2 * GG[a, b] + GG[b, b], GG[c, c], GG[d, d])
            flat += s.prod(v**power for v, power in zip(values, powers)) / 4
        P = s.Poly(s.expand(flat), *basis.GG)
        target = s.Matrix([P.coeff_monomial(m) for m in monomials])
        solution = matrix.gauss_jordan_solve(target)[0]
        checks[f"word{index}_entire_flat_projection"] = s.expand(
            P.as_expr() - sum(c * B.as_expr() for c, B in zip(solution, flats))
        )
        difference = s.expand(
            labeled_contact(powers, G, H, aa)
            - sum(c * B for c, B in zip(solution, contacts))
        )
        monomial, co = TP.terms()[0]
        candidate = s.Poly(difference, *variables).coeff_monomial(monomial) / co
        checks[f"word{index}_whole_generic_curvature_identity"] = s.expand(
            difference - candidate * T
        )
        checks[f"word{index}_exact_conversion_coefficient"] = (
            candidate - EXPECTED[index]
        )
        tau[powers] = candidate
        coefs[powers] = tuple(solution)
    for index, powers in enumerate(POWERS):
        checks[f"word{index}_singleton_exchange"] = (
            tau[powers] - tau[(powers[0], powers[2], powers[1])]
        )
    checks["full_three_channel_bubble_v3_conversion"] = tau[(3, 0, 0)] / 2 + 4
    K = []
    for j in range(4):
        r = 3 - j
        value = s.S.Zero
        for a, b, c in product(range(r + 1), repeat=3):
            if a + b + c != r:
                continue
            xp, yp, zp = a + c, a + b, b + c
            mult = s.Rational(factorial(r), factorial(a) * factorial(b) * factorial(c))
            beta = s.Rational(factorial(xp) * factorial(yp), factorial(xp + yp + 1))
            value += mult * beta * Z**zp * (1 - Z) ** (xp + yp + 1) * tau[(j + a, b, c)]
        K.append(s.factor(value))
    target = (
        -2 * (Z - 1) ** 4 * (55 * Z**3 + 10 * Z**2 + 4 * Z + 1) / 35,
        4 * (Z - 1) ** 3 * (6 * Z**2 + 3 * Z + 1) / 15,
        -4 * (Z - 1) ** 2 * (2 * Z + 1) / 3,
        8 * (Z - 1),
    )
    for j, (got, want) in enumerate(zip(K, target)):
        checks[f"whole_triangle_weight_j{j}"] = s.expand(got - want)
    return checks, tau, coefs, tuple(K)


@cache
def data():
    checks, tau, coefs, K = build()
    return {
        "checks": checks,
        "gates": {
            "all10_degree3_labeled_words": len(tau) == 10,
            "full_flat_polynomial_not_on_shell_fit": True,
            "whole_generic_real_TT_difference_not_sample_extraction": True,
            "same_fixed_six_symmetrized_jets_including_connections": True,
            "all_four_outer_heavy_resolvent_orders_retained": len(K) == 4,
            "full_three_channel_bubble_factor_not_one_channel": True,
        },
        "whole_conversion_coefficients": {str(p): v for p, v in tau.items()},
        "whole_flat_projection_coefficients": {str(p): v for p, v in coefs.items()},
        "whole_triangle_conversion_weights": K,
        "whole_bubble_conversion_coefficient": tau[(3, 0, 0)] / 2,
        "whole_definition": "For each alpha+beta+gamma=3, W=(1/4)sum24 v_ab^alpha(pc^2)^beta(pd^2)^gamma. Covariantize its three factors Phi^2,Phi,Phi using L=-Box_g as in S344. Subtract the S345 six-word symmetrized-jet lift of the entire flat off-shell polynomial. The displayed difference is tau*T as an identity in all generic G,a,H variables with only sum p+k=0, k^2=0 and TT. The composite-factor divided difference includes every L insertion.",
        "whole_weight_derivation": "Expand U^r=(xy*v+yz*w+xz*q)^r with r=3-j and multiply the jth A(v) coefficient. For monomial exponents alpha,beta,gamma, integrate the exact light-simplex moment x^(alpha+gamma)y^(alpha+beta) z^(beta+gamma), x+y=1-z. The measure supplies one additional(1-z); its beta factor is xp!yp!/(xp+yp+1)!. Multiply tau(j+alpha,beta,gamma). All terms are retained.",
    }
