"""Exact coupled local complex structure and canonical propagation identities.

This chooses initial free modes on one local window.  It is not a global
vacuum prescription or the instantaneous ground state of Hnorm=E-P^T Omega Y.
"""

from functools import cache

import sympy as sp


@cache
def checks():
    angle = sp.Symbol("rotation_parameter", real=True)
    w1, w2 = sp.symbols("frequency1 frequency2", positive=True)
    cosine, sine = (1-angle**2)/(1+angle**2), 2*angle/(1+angle**2)
    O = sp.Matrix([[cosine, -sine], [sine, cosine]])
    R = O*sp.diag(w1, w2)*O.T
    W = O*sp.diag(w1**2, w2**2)*O.T
    U = O*sp.diag(1/sp.sqrt(2*w1), 1/sp.sqrt(2*w2))*O.T
    P = -sp.I*O*sp.diag(sp.sqrt(w1/2), sp.sqrt(w2/2))*O.T
    residuals = {"orthogonal_initial_basis": O.T*O-sp.eye(2),
                 "positive_frequency_square": R*R-W,
                 "coordinate_CCR": U*U.H-U.conjugate()*U.T,
                 "momentum_CCR": P*P.H-P.conjugate()*P.T,
                 "mixed_CCR": U*P.H-U.conjugate()*P.T-sp.I*sp.eye(2),
                 "initial_energy_matrix": (P.H*P+U.H*W*U)/2-R/2}
    # Generic real coupled Hamiltonian, with no diagonalization of W(t).
    w0, w01, w11, omega, dw0, dw01, dw11 = sp.symbols(
        "w0 w01 w11 omega dw0 dw01 dw11", real=True)
    W = sp.Matrix([[w0, w01], [w01, w11]])
    Wdot = sp.Matrix([[dw0, dw01], [dw01, dw11]])
    Omega = sp.Matrix([[0, omega], [-omega, 0]])
    identity, null = sp.eye(2), sp.zeros(2)
    generator = (-Omega).row_join(identity).col_join((-W).row_join(-Omega))
    symplectic = null.row_join(identity).col_join((-identity).row_join(null))
    energy = W.row_join(null).col_join(null.row_join(identity))
    energy_dot = Wdot.row_join(null).col_join(null.row_join(null))
    target = (Wdot+Omega*W-W*Omega).row_join(null).col_join(null.row_join(null))
    residuals["Hamiltonian_preserves_CCR"] = generator*symplectic+symplectic*generator.T
    residuals["full_covariant_energy_rate"] = energy_dot+generator.T*energy+energy*generator-target
    return {f"{name}_{i}{j}": sp.simplify(matrix[i, j])
            for name, matrix in residuals.items()
            for i in range(matrix.rows) for j in range(matrix.cols)}


def negative_controls():
    # With a nonzero connection, replacing the prescribed covariant P by the
    # ordinary Ydot changes the initial datum.  A symmetric added connection
    # would also destroy the canonical skew-generator and energy identities.
    o, y0, y1, p0, p1 = sp.symbols("omega y0 y1 p0 p1", real=True)
    return {"P_is_not_Ydot": [-o*y1, o*y0],
            "dropping_the_mixed_free_Hamiltonian": -o*(p0*y1-p1*y0),
            "one_scalar_species_omits_a_CCR_channel": sp.Matrix([[0, 0], [0, sp.I]])}


def exponential_enclosure():
    """Independent rational Taylor enclosure for the full-window exp(4/11)."""
    exponent = sp.Rational(4, 11)
    partial = sum(exponent**n/sp.factorial(n) for n in range(10))
    first_tail = exponent**10/sp.factorial(10)
    upper = partial+first_tail/(1-exponent/11)
    if not (1 < partial < upper < sp.Rational(11, 7) < 2):
        raise ValueError("Independent free-energy exponential enclosure failed")
    return {"exponent": str(exponent), "lower": str(partial), "upper": str(upper),
            "rational_energy_upper": "11/7"}
