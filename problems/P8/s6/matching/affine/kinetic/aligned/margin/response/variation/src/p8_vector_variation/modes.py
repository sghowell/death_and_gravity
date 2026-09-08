"""Full homogeneous mass-tensor modes and derivative-free canonical variation."""
from functools import cache

import sympy as sp

N, scale, mass2, q = sp.symbols("N scale mass_squared q", positive=True)
am, bm = sp.symbols("a_mass b_mass", positive=True)
n, zeta, alpha, beta = sp.symbols("lapse_variation physical_log_scale_variation alpha beta", real=True)
parameter = sp.Symbol("variation_parameter", real=True)


def clean(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.factor))


@cache
def canonical():
    sigma, velocity, temporal = sp.symbols("sigma sigma_dot temporal", real=True)
    L = scale**3*(q*(velocity-temporal)**2/(2*N)+mass2*am*temporal**2/(2*N)
                 -N*mass2*bm*q*sigma**2/2)
    solution = q*velocity/(q+mass2*am)
    longitudinal_g2 = scale**3*mass2*am*q/(N*(q+mass2*am))
    longitudinal_bare = N**2*(bm*q/am+mass2*bm)
    transverse_g2 = scale/N
    transverse_bare = N**2*(q+mass2*bm)
    return {"longitudinal_L": L, "temporal": solution,
            "transverse": {"g_squared": transverse_g2, "bare_frequency_squared": transverse_bare},
            "longitudinal": {"g_squared": longitudinal_g2, "bare_frequency_squared": longitudinal_bare},
            "temporal_Euler": sp.factor(sp.diff(L, temporal).subs(temporal, solution)),
            "longitudinal_reduced_action": sp.factor(
                L.subs(temporal, solution)-longitudinal_g2*(velocity**2-longitudinal_bare*sigma**2)/2)}


@cache
def first_variations():
    out = {}
    substitutions = {N: 1+parameter*n, scale: scale*sp.exp(parameter*zeta),
                     q: q*sp.exp(-2*parameter*zeta),
                     am: 1+parameter*alpha*n, bm: 1+parameter*beta*n}
    z = q/(q+mass2)
    targets = {"transverse": ((zeta-n)/2, -2*q*zeta+(2*(q+mass2)+mass2*beta)*n),
               "longitudinal": (((1+2*z)*zeta+(alpha*z-1)*n)/2,
                                -2*q*zeta+((2+beta)*(q+mass2)-alpha*q)*n)}
    for kind in ("transverse", "longitudinal"):
        data = canonical()[kind]
        r = sp.factor(sp.diff(sp.log(data["g_squared"].subs(substitutions, simultaneous=True)), parameter).subs(parameter, 0)/2)
        frequency = sp.factor(sp.diff(data["bare_frequency_squared"].subs(substitutions, simultaneous=True), parameter).subs(parameter, 0))
        out[kind] = {"delta_log_g": r, "delta_bare_frequency_squared": frequency,
                     "normalization_identity": sp.factor(r-targets[kind][0]),
                     "frequency_identity": sp.factor(frequency-targets[kind][1])}
    return out


@cache
def symplectic_variation():
    g = sp.Symbol("canonical_g", positive=True)
    d, dp, bare, r, rp, rpp, variation = sp.symbols(
        "d d_prime bare_frequency_squared r r_prime r_second delta_bare_frequency_squared", real=True)
    original = sp.Matrix([[0, g**-2], [-g**2*bare, 0]])
    delta_original = sp.Matrix([[0, -2*r/g**2], [-g**2*(variation+2*r*bare), 0]])
    C = sp.Matrix([[g, 0], [d*g, 1/g]])
    Cp = sp.Matrix([[d*g, 0], [g*(dp+d**2), -d/g]])
    oscillator = sp.Matrix([[0, 1], [-bare+dp+d**2, 0]])
    deltaC = sp.Matrix([[r*g, 0], [g*(rp+d*r), -r/g]])
    S = sp.Matrix([[r, 0], [rp+2*d*r, -r]])
    Sp = sp.Matrix([[rp, 0], [rpp+2*dp*r+2*d*rp, -rp]])
    changed = Sp+S*oscillator-oscillator*S+C*delta_original*C.inv()
    target = sp.Matrix([[0, 0], [-variation+rpp+2*d*rp, 0]])
    symplectic = sp.Matrix([[0, 1], [-1, 0]])
    return {"physical_canonical_matrix": original, "physical_canonical_variation": delta_original,
            "oscillator_map": C, "map_variation": deltaC, "oscillator_frequency_variation": variation-rpp-2*d*rp,
            "map_is_symplectic": clean(C.T*symplectic*C-symplectic),
            "time_dependent_map_gives_oscillator": clean(Cp*C.inv()+C*original*C.inv()-oscillator),
            "map_variation_retained": clean(deltaC*C.inv()-S),
            "two_causal_variation_forms_agree": clean(changed-target),
            "original_variation_has_no_metric_derivatives": not delta_original.has(rp, rpp)}


@cache
def cauchy_response():
    f, fc, fd, fcd = sp.symbols("f f_conjugate f_prime f_conjugate_prime")
    # G(u,s)=[f(s)f*(u)-f*(s)f(u)]/i, with W(f,f*)=i.
    coincident = sp.expand((f*fc-fc*f)/sp.I)
    slope = (f*fcd-fc*fd)/sp.I
    positive = sp.Symbol("positive_frequency", positive=True)
    time = sp.Symbol("time_separation", real=True)
    flat = (sp.exp(sp.I*positive*time)-sp.exp(-sp.I*positive*time))/(2*positive*sp.I)
    return {"retarded_G_diagonal_zero": coincident,
            "retarded_G_diagonal_slope_is_Wronskian_over_i": sp.factor(slope-(f*fcd-fd*fc)/sp.I),
            "retarded_G_diagonal_normalization": sp.factor(slope.subs(fcd, (sp.I+fd*fc)/f)-1),
            "flat_retarded_G_sign": sp.simplify(sp.expand_complex(flat)-sp.sin(positive*time)/positive)}


@cache
def checks():
    data = canonical()
    out = {name: data[name] for name in ("temporal_Euler", "longitudinal_reduced_action")}
    for kind, values in first_variations().items():
        out.update({kind+"_"+name: values[name] for name in ("normalization_identity", "frequency_identity")})
    s = symplectic_variation()
    out.update({name: s[name] for name in ("map_is_symplectic", "time_dependent_map_gives_oscillator",
               "map_variation_retained", "two_causal_variation_forms_agree")})
    out.update(cauchy_response())
    return out
