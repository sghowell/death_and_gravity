"""Literal quadratic ADM action in a non-H-dividing clock/spatial chart.

Units M=tau=1; time is u=T/tau and K=(tau*k_com)^2. Spatial metrics are
a² diag(exp(2qg),exp(2qg),exp(2qg+2eg)) and its f analogue. The common
spatial gauge is eg=0 and the time gauge is delta_phi=0, permitted by the
strictly positive actual clock speed. e denotes the f logarithmic scalar
shear, e=-k_com² E_f, not the transverse-vector shear.

The EH time boundary -[a³h Sg²-b³h Sf²/c]' has been removed before the
Legendre map; all its background derivatives remain in the bulk action.
The sourced clock and truly free chi are both retained. No background H
is set to zero before differentiating any coefficient.
"""

from functools import cache

import sympy as sp
from p8_variable_beta import background

K = sp.Symbol("K", positive=True)
QG, QF, E, Z = sp.symbols("q_g q_f e z", real=True)
PG, PF, PE, PX = sp.symbols("p_g p_f p_e p_chi", real=True)
VG, VF, VE, VZ = sp.symbols("q_g_prime q_f_prime e_prime z_prime", real=True)
NG, NF, BG, BF = sp.symbols("n_g n_f B_g B_f", real=True)


@cache
def derive():
    bg = background.derive()
    c, a, b, y, h, hp, w, n, kin, b0, b1, b4, P = (
        bg[key] for key in ("c", "a", "b", "y", "h", "h_u", "chi_speed",
                            "nbar", "kbar", "b0", "b1", "b4", "P"))
    sg, sf = 3*QG, 3*QF+E
    difference = sg-sf
    shift = a**5*y**2*P/(2*(c+y))
    potential = (a**3*(b0+c*b1)*sg**2
                 +a**3*y*b1*(2*(2*QG+QF)**2+(2*QG+QF+E)**2)
                 +c*b**3*b4*sf**2)
    p_eg = a**3*w*Z-PE
    clock_numerator = -h*PG+a**3*n*sg+a**3*y*P*difference+2*a*K*QG
    clock_momentum = (clock_numerator-w*PX)/sp.sqrt(kin)
    hg = (3*p_eg**2-2*PG*p_eg)/(4*a**3)
    hf = c*(3*PE**2-2*PF*PE)/(4*b**3)
    hs = K*PE**2/(4*shift)
    hm = ((clock_numerator-w*PX)**2/(2*a**3*kin)+PX**2/(2*a**3)
          -sg*clock_numerator+a**3*n*sg**2/4+a*K*Z**2/2)
    he = (-a**3*(hp+sp.Rational(3, 2)*h**2)*sg**2
          -b**3/c*(-hp+sp.Rational(3, 2)*h**2)*sf**2
          -a*K*QG**2-c*b*K*QF**2)
    constraint = -h*PF-2*c*b*K*QF+c*a**3*P*difference
    hamiltonian = hg+hf+hs+hm+he+potential
    return {**bg, "K": K, "coordinates": (QG, QF, E, Z),
            "momenta": (PG, PF, PE, PX), "state": (QG, QF, E, Z, PG, PF, PE, PX),
            "S_g": sg, "S_f": sf, "volume_difference": difference,
            "C_shift": shift, "potential_H2": potential,
            "p_Eg": p_eg, "clock_constraint_numerator": clock_numerator,
            "p_phi": clock_momentum, "H_g": hg, "H_f": hf,
            "H_shift": hs, "H_matter": hm, "H_EH_pump_gradient": he,
            "H": hamiltonian, "C_f": constraint,
            "B_g": (3*p_eg-PG)/(2*a**3*K),
            "EH_time_boundary": -a**3*h*sg**2+b**3*h*sf**2/c}


@cache
def lagrangian():
    """Uneliminated four lapse/shift equations in the same clock chart."""
    d = derive()
    a, b, c, h, hp, w, n, sg, sf, y = (
        d[key] for key in ("a", "b", "c", "h", "h_u", "chi_speed", "nbar",
                           "S_g", "S_f", "y"))
    lg = a**3*(-3*VG**2+6*h*NG*VG-3*h**2*NG**2+3*h**2*NG*sg
               +(hp+sp.Rational(3, 2)*h**2)*sg**2-2*K*(VG-h*NG)*BG)
    lg += a*K*(QG**2+2*NG*QG)
    lf = b**3/c*(-3*VF**2-2*VF*VE-2*h*NF*(3*VF+VE)-3*h**2*NF**2
                 +3*h**2*NF*sf+(-hp+sp.Rational(3, 2)*h**2)*sf**2
                 -2*K*(VF+h*NF)*BF)
    lf += c*b*K*(QF**2+2*NF*QF)
    lm = a**3*(VZ**2/2+w*VZ*(sg-NG)+n*(sg**2/2-NG*sg+NG**2)/2
               -w*K*BG*Z)-a*K*Z**2/2
    lp = (-d["potential_H2"]-2*a**3*d["b0"]*NG*sg
          -2*a**3*y*d["b1"]*NG*(6*QG+3*QF+E)
          -2*c*a**3*d["b1"]*NF*sg-2*c*b**3*d["b4"]*NF*sf
          +d["C_shift"]*K*(BF-BG)**2)
    total = lg+lf+lm+lp
    return {"L": total, "L_g": lg, "L_f": lf, "L_matter": lm, "L_potential": lp,
            "velocities": (VG, VF, VE, VZ), "auxiliaries": (NG, NF, BG, BF),
            "auxiliary_equations": tuple(sp.diff(total, x) for x in (NG, NF, BG, BF))}
