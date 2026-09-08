"""Positive-mass Liouville reduction without changing the physical metric."""
from functools import cache

import sympy as sp

a, N, am, bm, k, m = sp.symbols("scale lapse mass_a mass_b comoving_momentum mass", positive=True)
U, U1, U2, b, Acurv = sp.symbols("mass_potential mass_potential_first mass_potential_second pump_rate pump_curvature", real=True)


def sector(value):
    if type(value) is not str or value not in ("T", "L"):
        raise ValueError("Require native physical sector T or L")
    return value


def data(value):
    return _data(sector(value))


@cache
def _data(value):
    q = k**2/a**2
    if value == "T":
        rate = N/a
        g2 = a/N
        frequency2 = N**2*(q+m**2*bm)
        pump2 = sp.Integer(1)
        mass_potential = a**2*m**2*bm
    else:
        rate = N*sp.sqrt(bm/am)/a
        g2 = a**3*m**2*am*q/(N*(q+m**2*am))
        frequency2 = N**2*bm*(q/am+m**2)
        mass_potential = a**2*m**2*am
        pump2 = m**2*a**2*sp.sqrt(am*bm)*k**2/(k**2+mass_potential)
    return {"acoustic_rate": rate, "physical_g_squared": g2, "physical_frequency_squared": frequency2,
            "acoustic_pump_squared": pump2, "acoustic_mass_potential": mass_potential,
            "checks": {"positive_canonical_pump_"+value: sp.simplify(g2*rate-pump2),
                       "unit_acoustic_momentum_"+value: sp.simplify(frequency2/rate**2-k**2-mass_potential)}}


def longitudinal_potential():
    denominator = k**2+U
    remainder = (b*U1+U2/2)/denominator-3*U1**2/(4*denominator**2)
    return {"scalar_principal": U-Acurv, "momentum_remainder": remainder,
            "exact": U-Acurv+remainder}


@cache
def checks():
    out = {}
    for value in ("T", "L"):
        out.update(data(value)["checks"])
    denominator = k**2+U
    log_pump = b-U1/(2*denominator)
    log_pump_derivative = Acurv-b**2-U2/(2*denominator)+U1**2/(2*denominator**2)
    out["exact_longitudinal_pump_potential_split"] = sp.factor(
        U-log_pump_derivative-log_pump**2-longitudinal_potential()["exact"])
    # A separate direct derivative of B=A*k/sqrt(k^2+U).
    t = sp.Symbol("parameter_time", real=True)
    A0, A1, A2 = sp.symbols("pump0 pump1 pump2", positive=True)
    B = (A0+A1*t+A2*t**2/2)*k/sp.sqrt(k**2+U+U1*t+U2*t**2/2)
    curvature = (sp.diff(B, t, 2)/B).subs(t, 0)
    out["independent_longitudinal_pump_second_derivative"] = sp.factor(
        curvature-(U-longitudinal_potential()["exact"]).subs({b: A1/A0, Acurv: A2/A0}))
    # Time change d sigma=f dt and physical v -> sqrt(f)*v.
    f, f1, f2, g, g1, g2 = sp.symbols("rate rate1 rate2 g g1 g2", positive=True)
    v0, p0, v1, p1 = sp.symbols("v0 p0 v1 p1", complex=True)
    transformed = sp.ImmutableMatrix([[sp.sqrt(f), 0], [0, 1/sp.sqrt(f)]])
    symplectic = sp.ImmutableMatrix([[0, 1], [-1, 0]])
    out["acoustic_canonical_Wronskian_preserved"] = sp.ImmutableMatrix(
        transformed.T*symplectic*transformed-symplectic)
    out["acoustic_pair_Wronskian_explicit"] = sp.simplify(
        (sp.sqrt(f)*v0)*(p1/sp.sqrt(f))-(p0/sp.sqrt(f))*(sp.sqrt(f)*v1)-(v0*p1-p0*v1))
    C = (g+g1*t+g2*t**2/2)*sp.sqrt(f+f1*t+f2*t**2/2)
    Ccurv = ((sp.diff(C, t, 2)-f1*sp.diff(C, t)/f)/(f**2*C)).subs(t, 0)
    target = g2/(f**2*g)+f2/(2*f**3)-3*f1**2/(4*f**4)
    out["acoustic_pump_retains_time_change_Schwarzian"] = sp.simplify(Ccurv-target)
    vp = sp.Symbol("physical_mode_time_derivative", complex=True)
    fv = sp.sqrt(f)*v0
    fvd = (sp.sqrt(f)*vp+f1*v0/(2*sp.sqrt(f)))/f
    rate_C = (g1/g+f1/(2*f))/f
    out["moving_acoustic_momentum"] = sp.simplify(fvd-rate_C*fv-(vp-g1*v0/g)/sp.sqrt(f))
    weight_A, weight_B, physical_frequency2 = sp.symbols("readout_A readout_B physical_frequency_squared", real=True)
    old_readout = (weight_A*p0**2+weight_B*physical_frequency2*v0**2)/2
    new_readout = f*(weight_A*(p0/sp.sqrt(f))**2
                    +weight_B*(physical_frequency2/f**2)*(sp.sqrt(f)*v0)**2)/2
    out["physical_readout_retains_acoustic_rate"] = sp.simplify(old_readout-new_readout)
    return out
