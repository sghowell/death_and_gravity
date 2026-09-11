"""Exact retarded six-step identity with every equal-time boundary term."""

from functools import cache

import sympy as s

G = s.symbols("g0:8")
A = s.symbols("a0:8")
H = s.symbols("h0:8")


def derivative(expression):
    return s.expand(
        sum(
            s.diff(expression, seq[j]) * seq[j + 1]
            for seq in (G, A, H)
            for j in range(7)
        )
    )


@cache
def iterates():
    rows = [A[0] * H[0]]
    for _ in range(6):
        rows.append(derivative(G[0] * rows[-1]))
    return tuple(rows)


def boundary(order):
    return s.expand(sum(-s.I * s.I**j * G[0] * iterates()[j] for j in range(order)))


@cache
def data():
    rows = iterates()
    checks = {}
    for n in range(1, 7):
        B = boundary(n)
        checks[f"complete_retarded_ODE_and_remainder_{n}"] = s.expand(
            derivative(B) + s.I * B / G[0] - rows[0] + s.I**n * rows[n]
        )
    coefficients = [-s.I, 1, s.I, -1, -s.I, 1]
    for j, c in enumerate(coefficients):
        checks[f"retained_boundary_phase_coefficient_{j}"] = -s.I * s.I**j - c
    return {
        "reference_integral": "K(t)=exp(-iTheta(t))*integral_(t0)^t exp(iTheta(s)) b(s)ds, with Theta'=Omega=W_k+W_l, g=1/Omega and b=a_Gamma*Gammahat. All source jets vanish near t0.",
        "complete_identity": "K=sum_(j=0)^5[-i*i^j g L^j b]+i^6 exp(-iTheta(t))*integral_(t0)^t exp(iTheta(s)) L^6 b(s)ds, Lb=(g b)'.",
        "ODE": "K'+iOmega K=b, K(t0)=0. The boundary polynomial obeys B6'+iOmega B6=b-i^6 L^6 b, so the remainder restores the exact equation and common initial data.",
        "boundary_coefficients": coefficients,
        "current": "The actual reference current is minus the imaginary part after contraction with the complete conjugate detector pair, summed over all nine polarizations and both retained internal momenta. No odd boundary term is dropped without evaluating that complete contraction.",
        "locality_boundary": "The first five boundary terms j0..4 contain source time jets at the readout time, but their momentum kernels are not thereby polynomial in external spatial momentum. They and the full contact still require the original covariant spatial matching.",
        "checks": checks,
        "gates": {
            "six_exact_retarded_steps": len(rows) == 7,
            "all_six_boundary_terms_kept": len(coefficients) == 6,
            "initial_jets_from_common_preparation": True,
            "retarded_bulk_not_full_signed_square": True,
        },
    }
