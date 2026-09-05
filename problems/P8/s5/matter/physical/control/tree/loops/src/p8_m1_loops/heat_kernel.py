"""Real-scalar local heat coefficient and frozen-variable nonclosure checks.

The source theorem is Vassilevich hep-th/0306138, equations (1.16),
(1.18), (2.32), (3.6), (4.28), and section 2.3.  A2 here is the local
coefficient of s**2 after stripping (4*pi*s)**(-n/2), not his integrated a4.
Only chi is integrated out; the metric and clock are classical backgrounds.
"""

from functools import cache

import sympy as sp

R, RIC2, RIEM2, BOX_R = sp.symbols("R Ric2 Riem2 BoxR")
E, BOX_E, OMEGA2 = sp.symbols("E BoxE Omega2")
XI = sp.Symbol("xi", real=True)
HBAR = sp.Symbol("hbar", nonnegative=True)
LOOP = 1/(4*sp.pi)**2
EULER = RIEM2-4*RIC2+R**2
WEYL2 = RIEM2-2*RIC2+R**2/3


def universal_a2():
    """Equation (4.28) local integrand, including its BoxR representative."""
    return (60*BOX_E+60*R*E+180*E**2+12*BOX_R+5*R**2
            - 2*RIC2+2*RIEM2+30*OMEGA2)/360


def scalar_a2(xi=0):
    """Constant xi in DE=-nabla_E**2+xi*R; xi=0 is the frozen M1 action."""
    xi = sp.sympify(xi)
    return sp.expand(universal_a2().subs({E: -xi*R, BOX_E: -xi*BOX_R, OMEGA2: 0}))


def bulk_coefficients(xi=0):
    """Four-dimensional C2/E4/R2/BoxR coefficients before any bulk quotient."""
    xi = sp.sympify(xi)
    return {"C2": sp.Rational(1, 120), "E4": -sp.Rational(1, 360),
            "R2": (xi-sp.Rational(1, 6))**2/2,
            "BoxR": sp.Rational(1, 30)-xi/6}


def decomposition_checks():
    coefficients = bulk_coefficients(XI)
    reconstructed = (coefficients["C2"]*WEYL2+coefficients["E4"]*EULER
                     + coefficients["R2"]*R**2+coefficients["BoxR"]*BOX_R)
    minimal_ricci_basis = EULER/180+RIC2/60+R**2/120+BOX_R/30
    completed_square = ((RIEM2-RIC2)/180+(XI-sp.Rational(1, 6))**2*R**2/2
                        + (sp.Rational(1, 30)-XI/6)*BOX_R)
    return {"weyl_euler_basis": sp.expand(scalar_a2(XI)-reconstructed),
            "minimal_ricci_basis": sp.expand(scalar_a2()-minimal_ricci_basis),
            "nonminimal_completed_square": sp.expand(scalar_a2(XI)-completed_square)}


@cache
def convention_checks():
    """Residue and constant-scale bookkeeping, with all signs defined here.

    A positive auxiliary proper-time IR regulator extracts only the UV residue;
    it is not a mass added to M1 or a claim that its full determinant exists.
    Under GammaE=-i*GammaL|Wick for P8 +---, a curvature-squared local
    coefficient obeys cE=-cL.  Counterterms have the opposite pole sign.
    """
    eps, regulator, mu = sp.symbols("eps regulator mu", positive=True)
    proper_time_a2 = -HBAR*mu**(2*eps)*sp.gamma(eps)*regulator**(-eps)/(2*(4*sp.pi)**(2-eps))
    residue = sp.limit(eps*proper_time_a2, eps, 0, dir="+")
    n = sp.Symbol("n")
    epsilon_pole = -HBAR*LOOP/(2*eps)
    n_pole = HBAR*LOOP/(n-4)
    zeta0, zeta_prime = sp.symbols("zeta0 zeta_prime")
    renormalized = -HBAR*zeta_prime/2-HBAR*sp.log(mu**2)*zeta0/2
    rescaling = mu*sp.diff(renormalized, mu)
    euclidean_local = sp.Symbol("cE")
    # dt=-i dtE, and the two metric-sign flips in R2/C2 cancel.
    lorentzian_local = -euclidean_local
    continued_action_coefficient = -sp.I*lorentzian_local
    return {"real_scalar_proper_time_residue": sp.simplify(residue+HBAR*LOOP/2),
            "n_minus_four_vs_epsilon_pole": sp.simplify(n_pole.subs(n, 4-2*eps)-epsilon_pole),
            "constant_scale_rescaling": sp.simplify(rescaling+HBAR*zeta0),
            "wick_local_coefficient_sign": sp.simplify(-sp.I*continued_action_coefficient-euclidean_local)}


def convention_table():
    unit = HBAR*LOOP
    return {"dimension": "n=4-2*epsilon",
            "Euclidean_effective_action": "GammaE=(hbar/2)*Tr log(DE/mu^2)",
            "Euclidean_pole_coefficient_per_A2_over_epsilon": str(-unit/2),
            "Euclidean_counterterm_coefficient_per_A2_over_epsilon": str(unit/2),
            "Euclidean_pole_coefficient_per_A2_over_n_minus_4": str(unit),
            "P8_Lorentzian_pole_coefficient_per_A2_over_epsilon": str(unit/2),
            "P8_Lorentzian_counterterm_coefficient_per_A2_over_epsilon": str(-unit/2),
            "constant_scale_GammaE_derivative_per_A2": str(-unit),
            "constant_scale_GammaL_derivative_per_A2": str(unit),
            "Euclidean_compensating_beta_per_bulk_coefficient": str(unit),
            "P8_Lorentzian_compensating_beta_per_bulk_coefficient": str(-unit),
            "continuation": "GammaE=-i*GammaL|continued; dt=-i*dtE; gL|continued=-gE; cE=-cL",
            "locality_only": "no global Euclidean CD determinant or global state is assumed"}


@cache
def nonclosure_witness():
    """Off-shell acceleration Hessian inside phi=t, N=1, X=1, a>0.

    No metric redefinitions or equations of motion enter this test.  A total
    derivative cannot change the nonzero fourth-order Euler-Lagrange symbol.
    """
    a = sp.Symbol("a", positive=True)
    adot, addot, a3, a4 = sp.symbols("adot addot a3 a4", real=True)
    jets = (a, adot, addot, a3, a4)

    def dt(expression):
        return sum(sp.diff(expression, jets[j])*jets[j+1] for j in range(4))

    ricci = -6*(addot/a+adot**2/a**2)
    lagrangian = sp.expand(a**3*ricci**2/72)
    el = sp.factor(sp.diff(lagrangian, a)-dt(sp.diff(lagrangian, adot))
                   + dt(dt(sp.diff(lagrangian, addot))))
    f, g = sp.symbols("f g")
    strict_restriction = f*a**2*addot+g*a*adot**2
    return {"a": a, "adot": adot, "addot": addot, "a3": a3, "a4": a4,
            "lagrangian": lagrangian, "acceleration_hessian": sp.diff(lagrangian, addot, 2),
            "euler_lagrange": el, "fourth_derivative_coefficient": sp.diff(el, a4),
            "strict_class_acceleration_hessian": sp.diff(strict_restriction, addot, 2)}


def field_redefinition_checks():
    """Illustrative Einstein+free-chi trade; emphatically not the CD EOM.

    delta g^ab=(alpha R^ab+beta R g^ab)/M^2 removes the Ricci-basis
    curvature term from -M^2 R/2, but induces nonminimal matter operators.
    """
    alpha, beta, mass, y, ruu = sp.symbols("alpha beta M Y Ruu")
    original = RIC2/60+R**2/120
    delta_einstein = -alpha*RIC2/2+(alpha/4+beta/2)*R**2
    delta_matter = (alpha*ruu-(alpha/2+beta)*R*y)/(2*mass**2)
    choice = {alpha: sp.Rational(1, 30), beta: -sp.Rational(1, 30)}
    induced = sp.expand(delta_matter.subs(choice))
    einstein_on_shell = induced.subs({ruu: y**2/mass**2, R: y/mass**2})
    return {"gravity_term_traded": sp.expand(original+delta_einstein.subs(choice)),
            "induced_matter_terms": sp.expand(induced-ruu/(60*mass**2)-R*y/(120*mass**2)),
            "illustrative_einstein_on_shell_Y_squared": sp.expand(einstein_on_shell-y**2/(40*mass**4))}


def controls():
    """Positive controls and nonzero omissions; these do not redefine M1."""
    c2 = sp.Symbol("C2")
    a = nonclosure_witness()["a"]
    # Setting R=0 need not set Weyl to zero.  Work modulo Euler and BoxR.
    return {"conformal_R2_coefficient": bulk_coefficients(sp.Rational(1, 6))["R2"],
            "conformal_C2_coefficient": bulk_coefficients(sp.Rational(1, 6))["C2"],
            "R_flat_bulk_Weyl_term": c2/120,
            "Ricci_flat_local_is_Euler": sp.expand(scalar_a2().subs({R: 0, RIC2: 0, BOX_R: 0})-RIEM2/180),
            "omit_real_scalar_half_residue_error": -HBAR*LOOP/2,
            "omit_R2_acceleration_hessian_error": -a,
            "flip_E_sign_at_conformal_coupling_error": sp.Rational(1, 18),
            "using_Euclidean_as_Lorentzian_beta_error_per_R2": HBAR*LOOP/36}
