"""Full physical pair readout, including the additional trace-constraint channel."""

import sympy as s

from . import geometry


def amplitudes(O1, O2, W1, W2, P1, P2, mass_a_over_r, sector):
    norm = 1 / (2 * s.sqrt(W1 * W2))
    pp = P1 * P2
    if sector == "TT":
        return ((mass_a_over_r**2 - pp) * norm, norm)
    if sector == "TL":
        return (mass_a_over_r * (O2 - pp / O2) * norm,)
    if sector == "LT":
        return (mass_a_over_r * (O1 - pp / O1) * norm,)
    if sector == "LL":
        return (
            (O1 * O2 - mass_a_over_r**2 * pp / (O1 * O2)) * norm,
            -pp * norm / (2 * O1 * O2),
        )
    raise ValueError("Require an actual ordered polarization sector")


def readout(khat, ellhat, ek, el, D, O1, O2, W1, W2, P1, P2, mass_a_over_r, sector):
    tau = s.trace(D)
    BD = D - tau * s.eye(D.rows) / 2
    amp = amplitudes(O1, O2, W1, W2, P1, P2, mass_a_over_r, sector)
    projected = (ek.T * BD * el)[0]
    value = amp[0] * projected
    if sector == "TT":
        value += amp[1] * (ek.T * geometry.magnetic(khat, ellhat, D) * el)[0]
    if sector == "LL":
        value += amp[1] * tau * (ek.T * khat)[0] * (ellhat.T * el)[0]
    return value
