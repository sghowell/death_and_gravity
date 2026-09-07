"""Exact Poeschl--Teller connection with the two physical-u phase factors.

This is a classical ODE transfer, not a particle-production statement. The
left/right power coordinates have opposite radial orientations, so their
2x2 transfer has determinant -1, not the canonical Cauchy determinant.
"""

from functools import cache

import sympy as sp

from .exact import MU, parameters

PHASE_SCALE = 4*sp.sqrt(2)


def coefficients():
    mu = MU
    aa = (sp.gamma(1-sp.I*mu)*sp.gamma(-sp.I*mu)
          /(sp.gamma(sp.Rational(3, 2)-sp.I*mu)*sp.gamma(-sp.Rational(1, 2)-sp.I*mu)))
    bb = sp.I/sp.sinh(sp.pi*mu)
    return {"mu": mu, "A": aa, "B": bb,
            "A_abs_squared": sp.coth(sp.pi*mu)**2,
            "B_abs_squared": sp.csch(sp.pi*mu)**2,
            "right_Jost": "exp(i*mu*z)*2F1(-1/2,3/2;1-i*mu;(1-tanh(z))/2)",
            "left_Jost_asymptotic": "A*exp(i*mu*z)+B*exp(-i*mu*z)"}


def power_transfer(delta, kbar=1):
    """(L+,L-) -> (R+,R-), powers |u|^(1/2 +/- i*mu)."""
    delta, _ = parameters(delta, kbar)
    d = coefficients()
    phase = sp.exp(2*sp.I*MU*sp.log(PHASE_SCALE/sp.sqrt(delta)))
    return sp.Matrix([[-sp.conjugate(d["B"]), sp.conjugate(d["A"])*phase],
                      [d["A"]/phase, -d["B"]]])


def unitary_power_coordinates():
    """Real cosine/sine coefficients -> normalized complex power coefficients.

    Q=sqrt(r)*(C1 cos(mu log r)+C2 sin(mu log r))
     =sqrt(r)/sqrt(2)*(Z+ r^(i mu)+Z- r^(-i mu)).
    Hence Z-=conjugate(Z+) for real data and the norm is unchanged.
    """
    return sp.Matrix([[1, -sp.I], [1, sp.I]])/sp.sqrt(2)


def real_transfer(delta, kbar=1):
    unitary = unitary_power_coordinates()
    return sp.diag(sp.eye(2), unitary.conjugate().T*power_transfer(delta, kbar)*unitary)


@cache
def checks():
    w, mu = sp.symbols("w mu", positive=True)
    f, fp, fpp = sp.symbols("F Fp Fpp")
    wp = -2*w*(1-w)
    wpp = sp.diff(wp, w)*wp
    transformed = wp**2*fpp+(wpp+2*sp.I*mu*wp)*fp+3*w*(1-w)*f
    hypergeometric = w*(1-w)*fpp+(1-sp.I*mu-2*w)*fp+sp.Rational(3, 4)*f
    a, ac, b, bc, phase = sp.symbols("A Abar B Bbar phase", nonzero=True)
    left_from_right = sp.Matrix([[b, ac*phase], [a/phase, bc]])
    right_from_left = sp.Matrix([[-bc, ac*phase], [a/phase, -b]])
    product = left_from_right*right_from_left-sp.eye(2)*(a*ac-b*bc)
    # Phase exponents follow by changing u=epsilon*x before identifying the
    # LEFT positive-frequency power. That power is the opposite right one.
    eps, scale, radius = sp.symbols("epsilon scale radius", positive=True)
    pplus, pminus = sp.Rational(1, 2)+sp.I*mu, sp.Rational(1, 2)-sp.I*mu
    expected = sp.exp(2*sp.I*mu*sp.log(eps/scale))
    scaling = sp.expand_log(sp.log(eps/scale), force=True)
    u = unitary_power_coordinates()
    return {
        "hypergeometric_ODE_pullback": sp.expand(transformed-4*w*(1-w)*hypergeometric),
        "gamma_B_reflection": sp.simplify((sp.pi/sp.sin(sp.pi*sp.I*mu))/(-sp.pi)-sp.I/sp.sinh(sp.pi*mu)),
        "Wronskian_modulus": sp.simplify(sp.coth(sp.pi*mu)**2-sp.csch(sp.pi*mu)**2-1),
        "power_connection_inverse": sum(sp.expand(v)**2 for v in product),
        "orientation_determinant": sp.expand(right_from_left.det()+a*ac-b*bc),
        "physical_scaling_phase": sp.simplify(sp.exp((pplus-pminus)*(sp.log(eps)-sp.log(scale)))
                                              -expected.subs(sp.log(eps/scale), scaling)),
        "unitary_power_basis": sum(sp.simplify(v)**2 for v in u.conjugate().T*u-sp.eye(2)),
        "positive_frequency_reverses_at_left": sp.diff(mu*sp.log(-radius), radius)-mu/radius,
    }


def controls():
    """Exact nonzero controls; none is promoted to a low-energy verdict."""
    d = coefficients()
    phase = sp.Symbol("phase", nonzero=True)
    aa = sp.Symbol("A", nonzero=True)
    # Even B=0 has a phase-dependent off-diagonal transmission amplitude.
    return {"nonzero_reflection_squared": d["B_abs_squared"],
            "transmission_not_small": d["A_abs_squared"]-1-d["B_abs_squared"],
            "B_zero_phase_derivative": sp.diff(aa*phase, phase),
            "wrong_left_frequency_label_changes_phase": 2*MU,
            "canonical_Cauchy_symplectic_claim_from_power_determinant": False}
