"""Exact scalar-clock spatial-chart transformation of the frozen CD action."""
from functools import cache

import sympy as sp
from p8 import adm as old_adm
from p8_affine import dictionary

u = dictionary.u
X, s, N = sp.symbols("X s N", positive=True)
K, KK, V = sp.symbols("K_hat Khat_ij_squared V", real=True)


@cache
def coefficients():
    target = dictionary.target()
    f = -target["f"].subs(dictionary.x, -X)
    a3 = -target["alpha3"].subs(dictionary.x, -X)
    original = old_adm.coefficients()
    mapping = {original["X"]: X, original["F2"]: f,
               original["F2X"]: sp.diff(f, X), original["A1"]: 0, original["A3"]: a3}
    B, C, D, E = [sp.factor(original[name].subs(mapping, simultaneous=True)) for name in ("B", "C", "D", "E")]
    omega = -sp.log((dictionary.lift()["h"]-1+X)/dictionary.lift()["h"])/4
    F = target["scalar_F"].subs(dictionary.x, -X)
    braid = -target["braiding"].subs(dictionary.x, -X)
    if braid != 0:
        raise ValueError("This exact CD branch must have the frozen zero metric braiding function")
    return {"f": f, "f_phi": sp.diff(f, u), "B": B, "C": C, "D": D, "E": E,
            "omega": omega, "omega_X": sp.diff(omega, X), "omega_phi": sp.diff(omega, u), "F": F}


@cache
def transform():
    data = coefficients()
    sub = {X: s**2}
    B, C, D, fphi, ox, op = [data[name].subs(sub) for name in
                             ("B", "C", "D", "f_phi", "omega_X", "omega_phi")]
    Q = s*op+2*s*ox*V
    Kphysical, KKphysical = K+3*Q, KK+2*Q*K+3*Q**2
    literal = B*(Kphysical**2-KKphysical)+C*Kphysical*V+D*V**2+2*s*fphi*Kphysical+data["F"].subs(sub)
    Lv = sp.factor(sp.diff(literal, V).subs(V, 0))
    expected_Lv = 12*s**2*fphi*ox
    Blinear = sp.factor(4*B*s*op+2*s*fphi)
    Fnew = data["F"].subs(sub)+6*B*s**2*op**2+6*s**2*fphi*op
    expected = B*(K**2-KK)+Blinear*K+Fnew+expected_Lv*V
    B4 = -data["f"]
    delta = -2*X*data["omega_X"]
    # R3 conformal transformation and its spatial IBP are both kept.
    gradient = data["E"]+2*B4*delta**2+4*(B4-2*X*sp.diff(B4, X))*delta
    volume = ((dictionary.lift()["h"]-1+s**2)/dictionary.lift()["h"])**(-sp.Rational(3, 4))
    return {"literal": literal, "B": B, "Blinear": Blinear, "Fnew": Fnew,
            "volume": volume, "Lv": Lv, "primitive_s": volume*Lv,
            "trace_velocity_residual": sp.factor(sp.diff(literal, K, V)),
            "lapse_velocity_squared_residual": sp.factor(sp.diff(literal, V, 2)),
            "remaining_linear_V": sp.factor(Lv-expected_Lv),
            "full_time_action_residual": sp.factor(literal-expected),
            "full_spatial_lapse_gradient_residual": sp.factor(gradient)}


@cache
def clock_lapse():
    from p8_affine_aligned import dynamics
    bg = dynamics.old.background()
    data = transform()
    # I_s=volume*Lv with I(u,1)=0 removes the linear normal-s derivative.
    I0 = sp.Integer(0)
    I_N = sp.factor(-data["primitive_s"].subs(s, 1))
    I_NN = sp.factor((sp.diff(data["primitive_s"], s)+2*data["primitive_s"]).subs(s, 1))
    I = I0+I_N*(N-1)+I_NN*(N-1)**2/2
    Iphi = sp.diff(I, u)
    volume, B, Blinear, Fnew = [data[name].subs(s, 1/N) for name in ("volume", "B", "Blinear", "Fnew")]
    coefficient = volume*B*sp.Rational(2, 3)
    linear = volume*Blinear-I
    potential = volume*Fnew-Iphi/N
    p = -2*bg["H"]
    H = N*((p-linear)**2/(4*coefficient)-potential+bg["ell"]**2/(2*volume))
    first = sp.factor(sp.diff(H, N).subs(N, 1))
    second = sp.factor(sp.diff(H, N, 2).subs(N, 1))
    return {"I_N": I_N, "I_NN": I_NN, "first": first, "second": second,
            "actual_background_lapse_equation": first,
            "actual_secondary_lapse_Jacobian": sp.factor(second+2*bg["J"])}


@cache
def checks():
    data = transform()
    names = ("trace_velocity_residual", "lapse_velocity_squared_residual", "remaining_linear_V",
             "full_time_action_residual", "full_spatial_lapse_gradient_residual")
    out = {name: data[name] for name in names}
    out.update({name: clock_lapse()[name] for name in
                ("actual_background_lapse_equation", "actual_secondary_lapse_Jacobian")})
    return out
