"""Literal curvature contraction, all-64 inverse and curved commutator.

Source signature (-+++), derivative index last. Only the quadratic
rolling connection Hessian is specialized to the original p=1/2.
"""
from functools import cache
from itertools import combinations

import sympy as sp
from p8_affine import connection

PAIRS = tuple(combinations(range(4), 2))
K = sp.symbols("k0:4", real=True)
SIGNS = connection.SIGN
PAIRING = sp.ImmutableMatrix(sp.diag(*(SIGNS[a]*SIGNS[b] for a, b in PAIRS)))


def _matrix(value):
    return sp.ImmutableMatrix(value.applyfunc(sp.factor))


def two_form(values):
    result = sp.zeros(4)
    for value, (a, b) in zip(values, PAIRS, strict=True):
        result[a, b], result[b, a] = value, -value
    return result


@cache
def flat_map():
    """C=(F13-F14)_[mu,nu], with the conventional factor 1/2."""
    rows = {}
    for mu in range(4):
        for nu in range(4):
            row = sp.zeros(1, 64)
            for a in range(4):
                row[connection.index(a, nu, mu)] += K[a]
                row[connection.index(a, nu, a)] -= K[mu]
                row[connection.index(nu, a, mu)] -= SIGNS[nu]*SIGNS[a]*K[a]
                row[connection.index(nu, a, a)] += SIGNS[nu]*SIGNS[a]*K[mu]
            rows[mu, nu] = row
    full = sp.Matrix.vstack(*((rows[a, b]-rows[b, a])/2 for a, b in PAIRS))
    return _matrix(full)


@cache
def flat_inverse():
    data = connection.quotient()
    embedding = data["embedding"]
    n = flat_map()*embedding
    w = sp.zeros(60, 6)
    for indices, block in zip(data["blocks"], data["block_matrices"], strict=True):
        piece = block.subs(connection.P, sp.Rational(1, 2)).inv()*n[:, list(indices)].T
        for i, index in enumerate(indices):
            w[index, :] = piece[i, :]
    w = _matrix(w)
    lower = two_form(sp.symbols("Z0:6", real=True))
    raised = sp.diag(*SIGNS)*lower
    expected = sp.Matrix([K[c]*raised[a, b] for a, b, c in connection.INDICES])
    expected_map = expected.jacobian(tuple(lower[a, b] for a, b in PAIRS))
    m = data["hessian"].subs(connection.P, sp.Rational(1, 2))
    full_m = connection.quadratic()["hessian"].subs(connection.P, sp.Rational(1, 2))
    return {"M": m, "N": _matrix(n), "W": w, "embedding": embedding,
            "D": _matrix(n*w), "covariant_gradient_map": _matrix(expected_map),
            "gradient_lift_residual": _matrix(embedding*w*PAIRING-expected_map),
            "quotient_Euler_residual": _matrix(m*w-n.T),
            "full_64_Euler_residual": _matrix(full_m*embedding*w-flat_map().T),
            "projective_residual": _matrix(flat_map()*connection.quadratic()["gauge"])}


@cache
def commutator():
    """C[nabla Z]=[nabla_a,nabla_mu]Z^a_nu-(mu<->nu)."""
    acceleration, spatial = sp.symbols("H_dot_plus_H_squared H_squared", real=True)
    riemann = {}
    for a, b, c, d in sp.utilities.iterables.cartes(range(4), repeat=4):
        value = sp.S.Zero
        if a == 0 and b > 0:
            value = acceleration*((c == 0 and d == b)-(d == 0 and c == b))
        elif a > 0 and b == 0:
            value = acceleration*((c == 0 and d == a)-(d == 0 and c == a))
        elif a > 0 and b > 0:
            value = spatial*(int(a == c and b == d)-int(a == d and b == c))
        riemann[a, b, c, d] = value
    ricci = sp.Matrix(4, 4, lambda b, d: sum(riemann[a, b, a, d] for a in range(4)))
    z = sp.symbols("Z0:6", real=True)
    lower = two_form(z)
    raised = sp.diag(*SIGNS)*lower

    def contracted(mu, nu):
        return sum(ricci[r, mu]*raised[r, nu] for r in range(4))-sum(
            riemann[r, nu, a, mu]*raised[a, r] for a in range(4) for r in range(4))

    result = sp.Matrix([contracted(mu, nu)-contracted(nu, mu) for mu, nu in PAIRS])
    endomorphism = _matrix(result.jacobian(z))
    return {"acceleration": acceleration, "spatial": spatial,
            "Riemann": riemann, "Ricci": ricci, "endomorphism": endomorphism,
            "expected": 2*(acceleration+spatial)*sp.eye(6),
            "scalar_curvature": sum(SIGNS[i]*ricci[i, i] for i in range(4))}


@cache
def off_clock_symbol():
    """Conditional constant-coefficient response, not a physical mode verdict.

Rotations put a spatial covector along z; its frequency stays arbitrary.
This checks why the null on-clock identity is not an open-tube identity.
"""
    omega, spatial_k = sp.symbols("omega spatial_k", real=True)
    data = connection.quotient()
    n = flat_map().subs(dict(zip(K, (omega, 0, 0, spatial_k), strict=True)))*data["embedding"]
    response = sp.zeros(6)
    for indices, block in zip(data["blocks"], data["block_matrices"], strict=True):
        rows = n[:, list(indices)]
        response += rows*block.inv(method="DM")*rows.T
    response = _matrix(response)
    p = connection.P
    electric = (2*p-1)**2*(3*p+1)/(2*p*(8*p**2-1))
    magnetic = (2*p-1)**2/(8*p**2)
    expected = spatial_k**2*sp.diag(-electric, -electric, 0, 0, magnetic, magnetic)
    return {"omega": omega, "spatial_k": spatial_k, "response": response,
            "expected": _matrix(expected), "electric_coefficient": electric,
            "magnetic_coefficient": magnetic,
            "additional_propagating_mode_or_ghost_claim": False}


@cache
def checks():
    data, curved, other = flat_inverse(), commutator(), off_clock_symbol()
    return {"projective_kernel": data["projective_residual"],
            "all_block_Euler": data["quotient_Euler_residual"],
            "all_64_Euler": data["full_64_Euler_residual"],
            "literal_covariant_gradient_lift": data["gradient_lift_residual"],
            "all_covector_components_null_Schur": data["D"],
            "FLRW_commutator_endomorphism": _matrix(curved["endomorphism"]-curved["expected"]),
            "FLRW_R_over_three": sp.factor(curved["scalar_curvature"]/3
                                             -2*(curved["acceleration"]+curved["spatial"])),
            "off_clock_conditional_Schur_formula": _matrix(other["response"]-other["expected"]),
            "off_clock_no_frequency_in_response": other["response"].diff(other["omega"])}
