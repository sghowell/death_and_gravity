"""Three constrained orthonormal vector modes, with longitudinal phases."""
from functools import cache

import sympy as sp

mass = sp.Symbol("mass", positive=True)
cosine, sine = sp.symbols("angle_cosine angle_sine", real=True)


def label(value):
    if type(value) is not int or value not in (0, 1):
        raise ValueError("Require native time or momentum index 0 or 1")
    return value


def angular(value):
    value = sp.expand(value)
    return sp.Poly(value, sine).rem(sp.Poly(sine**2+cosine**2-1, sine)).as_expr().expand()


def data(time, momentum):
    return _data(label(time), label(momentum))


@cache
def _data(time, momentum):
    suffix = str(time)+"_"+str(momentum)
    return {"k": sp.Symbol("physical_k_"+suffix, positive=True),
            "omega": sp.Symbol("omega_"+suffix, positive=True),
            "vT": sp.Symbol("vT_"+suffix),
            "vL": sp.Symbol("vL_"+suffix),
            "pL": sp.Symbol("pL_"+suffix)}


def direction(momentum):
    momentum = label(momentum)
    return sp.Matrix([0, 0, 1]) if momentum == 0 else sp.Matrix([sine, 0, cosine])


def polarization(time, momentum):
    time, momentum = label(time), label(momentum)
    return _polarization(time, momentum)


@cache
def _polarization(time, momentum):
    d = data(time, momentum)
    e1 = sp.Matrix([1, 0, 0]) if momentum == 0 else sp.Matrix([cosine, 0, -sine])
    e2 = sp.Matrix([0, 1, 0])
    result = [sp.ImmutableMatrix([0, *(e*d["vT"])]) for e in (e1, e2)]
    result.append(sp.ImmutableMatrix([d["k"]*d["pL"]/(mass*d["omega"]),
                                      *(sp.I*d["omega"]*direction(momentum)*d["vL"]/mass)]))
    # These are a^(3/2) times the physical orthonormal components.
    return tuple(result)


def wightman(momentum):
    momentum = label(momentum)
    return sp.ImmutableMatrix(sum((left*sp.conjugate(right).T for left, right in
                                  zip(polarization(0, momentum), polarization(1, momentum))),
                                 sp.zeros(4)))


@cache
def checks():
    out = {}
    for momentum in (0, 1):
        unit = direction(momentum)
        out["unit_direction_"+str(momentum)] = angular((unit.T*unit)[0]-1)
        x, y = data(0, momentum), data(1, momentum)
        W = wightman(momentum)
        P = sp.eye(3)-unit*unit.T
        target = sp.zeros(4)
        target[0, 0] = x["k"]*y["k"]*x["pL"]*sp.conjugate(y["pL"])/(mass**2*x["omega"]*y["omega"])
        for i in range(3):
            target[0, i+1] = -sp.I*x["k"]*y["omega"]*unit[i]*x["pL"]*sp.conjugate(y["vL"])/(mass**2*x["omega"])
            target[i+1, 0] = sp.I*x["omega"]*y["k"]*unit[i]*x["vL"]*sp.conjugate(y["pL"])/(mass**2*y["omega"])
            for j in range(3):
                target[i+1, j+1] = P[i, j]*x["vT"]*sp.conjugate(y["vT"])+x["omega"]*y["omega"]*unit[i]*unit[j]*x["vL"]*sp.conjugate(y["vL"])/mass**2
        out["constrained_Wightman_"+str(momentum)] = (W-target).applyfunc(angular)
        out["transverse_projector_rank_"+str(momentum)] = angular(sp.trace(P)-2)
        out["transverse_projector_idempotent_"+str(momentum)] = (P*P-P).applyfunc(angular)
        zero_momentum = {x["k"]: 0, y["k"]: 0, x["omega"]: mass, y["omega"]: mass,
                         x["vL"]: x["vT"], y["vL"]: y["vT"]}
        isotropic = sp.diag(0, *([x["vT"]*sp.conjugate(y["vT"])]*3))
        out["regular_isotropic_zero_momentum_"+str(momentum)] = (W.subs(zero_momentum)-isotropic).applyfunc(angular)
    H, wave, omega, vv, pp = sp.symbols("H wave omega v p", real=True)
    aa, mm = sp.symbols("a m", positive=True)
    d = H*(sp.Rational(1, 2)+wave**2/omega**2)
    U0 = wave*pp/(aa**sp.Rational(3, 2)*mm*omega)
    replacements = {wave: -H*wave, omega: -H*wave**2/omega,
                    aa: H*aa, pp: -omega**2*vv-d*pp}
    derivative = sum(sp.diff(U0, key)*value for key, value in replacements.items())
    longitudinal_space_divergence = -wave*omega*vv/(aa**sp.Rational(3, 2)*mm)
    out["on_shell_Proca_divergence"] = sp.factor(-derivative-3*H*U0+longitudinal_space_divergence)
    out["flat_longitudinal_polarization_norm"] = sp.factor((-wave**2+omega**2)/mm**2-1).subs(omega**2, mm**2+wave**2).expand()
    return out
