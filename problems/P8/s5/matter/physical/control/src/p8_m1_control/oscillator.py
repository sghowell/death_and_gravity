"""Exact volume/kinetic whitening, retained gyroscopic connection and energy.

The normalized Lagrangian is |Ydot+Omega Y|^2/2-Y^T W Y/2,
with W=q*I/ell^2+M/ell^2. Omega is antisymmetric, not discarded.
Off-diagonal quantities are stored as sqrt(ratio) times rational coefficients.
"""

from functools import cache

import sympy as sp

from . import model as m


@cache
def derive(chart):
    base, f = m.coefficients(chart), m.factors(chart)
    alpha, beta = base["alpha"], base["beta"]
    symmetric, antisymmetric = (beta+beta.T)/2, (beta-beta.T)/2
    delta = f["delta"]
    dimensions = (delta, 0)
    R = sp.Matrix(2, 2, lambda i, j: m.canonical(base["gamma"][i, j]
                  +m.derivative(symmetric[i, j], 1+dimensions[i]+dimensions[j])
                  +3*m.H*symmetric[i, j]-alpha[i, j]/m.z))
    for value in R:
        m.no_high_frequency_pole(value)
    d1, d2, chi, ratio = (f[key] for key in ("d1", "d2", "chi", "ratio"))
    f1, f2, h, log_ratio = (f[key] for key in ("f1", "f2", "h", "log_root_ratio_derivative"))
    k = m.canonical(antisymmetric[0, 1]/d2)
    omega = m.canonical(k+h/2)
    # Actual Omega12=sqrt(ratio)*omega; Omega21=-Omega12.
    M11 = m.canonical((R[0, 0]-2*chi*R[0, 1]+chi**2*R[1, 1])/d1
                      -f1**2-m.derivative(f1, 1)+ratio*(-h**2-2*h*k+omega**2))
    M22 = m.canonical(R[1, 1]/d2-f2**2-m.derivative(f2, 1)+ratio*omega**2)
    M12 = m.canonical((R[0, 1]-chi*R[1, 1])/d2-h*f2+k*(f1-f2)
                      -(m.derivative(h, 1)+log_ratio*h)/2)
    # The symmetric momentum boundary maps the ORIGINAL phase momenta to the
    # covariant normalized momentum P=Ydot+Omega Y. This is not merely a free
    # velocity replacement in interacting vertices.
    S11 = m.canonical(f1-(symmetric[0, 0]-2*chi*symmetric[0, 1]+chi**2*symmetric[1, 1])/d1)
    S22 = m.canonical(f2-symmetric[1, 1]/d2)
    S12 = m.canonical(h/2-(symmetric[0, 1]-chi*symmetric[1, 1])/d2)
    # Covariant derivative of the mass. [Omega,M] is symmetric, and the
    # leading q*I commutes exactly with the connection.
    cov11 = m.canonical(m.derivative(M11, 2)+2*ratio*omega*M12)
    cov22 = m.canonical(m.derivative(M22, 2)-2*ratio*omega*M12)
    cov12 = m.canonical(m.derivative(M12, 2)+log_ratio*M12+omega*(M22-M11))
    regular = {"mass11": M11, "mass22": M22, "mass12_factor": M12,
               "connection_factor": omega, "covariant_mass11": cov11,
               "covariant_mass22": cov22, "covariant_mass12_factor": cov12}
    for value in regular.values():
        m.no_high_frequency_pole(value)
    return {"factors": f, "unwhitened_mass_residual": R, **regular,
            "momentum_boundary11": S11, "momentum_boundary22": S22,
            "momentum_boundary12_factor": S12}


def physical_matrices(chart, x_value, q_value):
    """Exact point evaluation in local ell=1 units, with explicit chart checks."""
    x_value, q_value = sp.Rational(x_value), sp.Rational(q_value)
    if abs(x_value) > 1 or q_value < 1000:
        raise ValueError("Require |x|<=1 and q>=1000")
    if chart == "unitary" and abs(x_value) < sp.Rational(1, 9):
        raise ValueError("Outside the certified unitary exterior")
    if chart == "gamma" and abs(x_value) > sp.Rational(1, 4):
        raise ValueError("Outside the certified gamma core")
    mapping = {m.x: x_value, m.z: 1/q_value, m.l: (1-x_value*x_value)**sp.Rational(11, 2)/10}
    data = derive(chart)
    evaluate = lambda value: sp.simplify(value.subs(mapping, simultaneous=True))
    root = sp.sqrt(evaluate(data["factors"]["ratio"]))
    off, connection = root*evaluate(data["mass12_factor"]), root*evaluate(data["connection_factor"])
    mass = sp.Matrix([[evaluate(data["mass11"]), off], [off, evaluate(data["mass22"])]])
    return {"mass": mass, "potential": q_value*sp.eye(2)+mass,
            "connection": sp.Matrix([[0, connection], [-connection, 0]])}


def tensor():
    mass = -6-24*m.x*m.x
    return {"mass": mass, "covariant_mass": m.derivative(mass, 2), "connection": sp.Integer(0)}
