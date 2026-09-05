"""Full CD/M1 phase-space vertices, all nonexceptional rational momenta.

The scalar pairs use the mixed-boundary convention of pinned p8.gamma.
Velocity substitution uses the full two-by-two momentum Hessian. Canonical
configuration normalization, amplitudes and interaction control remain open.
"""

from functools import cache

import sympy as sp
from p8_physical import geometry
from p8_physical import jets as j
from p8_physical.momentum import ZeroMomentumConstraint
from p8_physical.vertices import Leg, nonexceptional, polarization, tensor_basis
from p8_s5.lapse_series import stationary_series

from . import background, momentum, quadratic

__all__ = ["Leg", "hamiltonian_kernel", "tensor_basis"]


def fields(context, legs, bg, chart, time_point, representation, *, mixed_boundary=True):
    if representation not in ("phase", "velocity"):
        raise ValueError("Use phase or velocity representation")
    scalar, scalar_p, matter, matter_p = (context.jet() for _ in range(4))
    scalar_dot, matter_dot = context.jet(), context.jet()
    tensor, tensor_p, tensor_dot = j.zeros(context), j.zeros(context), j.zeros(context)
    for index, leg in enumerate(legs):
        amplitude = context.leg(index)
        if representation == "velocity" and leg.kind in ("s", "m", "s_dot", "m_dot"):
            q = sum(sp.Rational(k)**2 for k in leg.wave)
            response = quadratic.response(time_point, chart, q)
            column = 0 if leg.kind in ("s", "s_dot") else 1
            is_velocity = leg.kind.endswith("_dot")
            mat = response["alpha"] if is_velocity else response["beta"]
            scalar_p += mat[0, column]*amplitude
            matter_p += mat[1, column]*amplitude
            if leg.kind == "s":
                scalar += amplitude
            elif leg.kind == "m":
                matter += amplitude
            elif leg.kind == "s_dot":
                scalar_dot += amplitude
            else:
                matter_dot += amplitude
        elif leg.kind == "s":
            scalar += amplitude
        elif leg.kind == "p":
            scalar_p += amplitude
        elif leg.kind == "m":
            matter += amplitude
        elif leg.kind == "P":
            matter_p += amplitude
        elif leg.kind in ("t", "pi", "t_dot"):
            matrix, norm = polarization(leg)
            update = [[matrix[i, k]*amplitude/(norm if leg.kind == "pi" else 1) for k in range(3)] for i in range(3)]
            if leg.kind == "t":
                tensor = j.madd(tensor, update)
            elif leg.kind == "pi":
                tensor_p = j.madd(tensor_p, update)
            elif representation == "velocity":
                tensor_dot = j.madd(tensor_dot, update)
                tensor_p = j.madd(tensor_p, j.mscale(update, sp.Rational(1, 4)))
            else:
                raise ValueError("Velocity leg requested in phase representation")
        else:
            raise ValueError("Use phase s,p,m,P,t,pi or velocity s,s_dot,m,m_dot,t,t_dot legs")
    generator = context.jet()
    if chart == "unitary":
        zeta, pv = scalar, scalar_p
    elif chart == "gamma":
        def rescale_modes(jet, inverse=False):
            result = {}
            for mask, coefficient in jet.data.items():
                q = sum(k*k for k in context.wave[mask])
                if q == 0:
                    raise ZeroMomentumConstraint("Gamma canonical swap requires nonzero momentum")
                factor = 1/(2*q) if inverse else -2*q
                result[mask] = coefficient*context.convert(factor)
            return j.Jet(context, result)
        zeta, pv = rescale_modes(scalar_p, True), rescale_modes(scalar)
        generator = -bg["H"]*scalar*scalar_p
    else:
        raise ValueError("Use unitary or gamma phase chart")
    # Inverse of the canonical boundary 3*a^3*l*zeta*m; (a^3*l)'=0.
    old_p = pv+3*bg["l"]*matter if mixed_boundary else pv
    matter_density = bg["l"]+matter_p+(3*bg["l"]*zeta if mixed_boundary else 0)
    return {"zeta": zeta, "tensor": tensor, "old_scalar_p": old_p, "tensor_p": tensor_p,
            "matter": matter, "matter_density": matter_density, "generator": generator,
            "symplectic": scalar_p*scalar_dot+matter_p*matter_dot+j.contract(tensor_p, tensor_dot)}


def stationary_hamiltonian(A, L, Q):
    """The pinned quartic formula, without expanding the unused third lapse jet.

    L and Q may already contain nonlinear physical geometry. Their minimum
    physical degrees are one and two, so the invariant truncation is sufficient.
    Expanding n3 as a physical Fourier jet is unnecessary for H through degree4.
    """
    n1 = -L[1]/A[2]
    force2 = A[3]*n1**2/2+L[2]*n1+Q[1]
    orders = (A[0], L[0], Q[0]-L[1]**2/(2*A[2]),
              A[3]*n1**3/6+L[2]*n1**2/2+Q[1]*n1,
              A[4]*n1**4/24+L[3]*n1**3/6+Q[2]*n1**2/2-force2**2/(2*A[2]))
    return {"lapse": (n1, -force2/A[2]), "hamiltonian": orders}


@cache
def stationary_bridge_checks():
    A, L, Q = sp.symbols("a:5"), sp.symbols("l:4"), sp.symbols("q:3")
    expected, actual = stationary_series(A, L, Q), stationary_hamiltonian(A, L, Q)
    return {**{f"Hamiltonian_{k}": sp.cancel(actual["hamiltonian"][k]-expected["hamiltonian"][k])
               for k in range(5)},
            **{f"lapse_{k+1}": sp.cancel(actual["lapse"][k]-expected["lapse"][k]) for k in range(2)}}


def invariant_reduction(context, invariants, jets):
    sigma, rho, eta, shear2, z = invariants
    L = tuple(jets["B"][k]*sigma+jets["C"][k]*rho+jets["L"][k]*eta for k in range(4))
    Q = tuple(jets["D"][k]*sigma**2+jets["E"][k]*shear2+jets["M"][k]*eta**2+jets["Z"][k]*z
              for k in range(3))
    if jets["A"][1] != 0 or jets["A"][2] == 0:
        raise ValueError("Stationary background/nonzero lapse Hessian required")
    result = stationary_hamiltonian(jets["A"], L, Q)
    return sum(result["hamiltonian"], context.jet()), result


def construction(legs, time_point=None, chart="unitary", representation="phase", *, constraint_order=None,
                 mixed_boundary=True, subtract_clock_velocity=True):
    if not 1 <= len(legs) <= 4:
        raise ValueError("Only first through fourth physical order is supported")
    nonexceptional(legs)
    background.point(time_point)
    context = j.Context([leg.wave for leg in legs], parameter=background.rtime if time_point is None else None)
    bg = background.functions(time_point)
    f = fields(context, legs, bg, chart, time_point, representation, mixed_boundary=mixed_boundary)
    zeta, tensor = f["zeta"], f["tensor"]
    geo = geometry.derive(context, zeta, tensor)
    mom = momentum.derive(context, geo, zeta, tensor, f["old_scalar_p"], f["tensor_p"], bg["H"],
                          f["matter_density"], f["matter"], order=constraint_order)
    volume_inverse = geo["volume"].power(-1)
    mixed = j.mscale(j.mmul(mom["momentum"], geo["metric"]), volume_inverse)
    trace = j.trace(mixed)
    sigma = trace+3*bg["H"]
    shear2 = j.trace(j.mmul(mixed, mixed))-trace**2/3
    eta = f["matter_density"]*volume_inverse-bg["l"]
    z = sum(geo["inverse"][i][k]*f["matter"].derivative(i)*f["matter"].derivative(k)
            for i in range(3) for k in range(3))
    invariants = (sigma, geo["curvature"], eta, shear2, z)
    for value in invariants:
        if not value.homogeneous(0).is_zero():
            raise ValueError("A perturbative invariant has nonzero background")
    if not shear2.homogeneous(1).is_zero() or not z.homogeneous(1).is_zero():
        raise ValueError("A weight-two invariant has a linear term")
    density, lapse_result = invariant_reduction(context, invariants, background.lapse_jets(time_point))
    h = geo["volume"]*density-2*bg["H"]*j.contract(mom["momentum"], geo["metric"])
    h -= bg["Adot"]*(6*zeta+3*zeta**2-j.contract(tensor, tensor)/2)
    if subtract_clock_velocity:
        h -= bg["l"]*f["matter_density"]
    h += f["generator"]
    return {"context": context, "hamiltonian": h, "momentum": mom, "geometry": geo, "fields": f,
            "invariants": invariants, "lapse_result": lapse_result, "background": bg}


@cache
def hamiltonian_kernel(legs, time_point=None, chart="unitary", representation="phase"):
    """Labelled integrated coefficient; repeated fields use distinct legs, no 1/n!."""
    result = construction(legs, time_point, chart, representation)
    ctx = result["context"]
    return {"kernel": sp.factor(result["hamiltonian"].coefficient(ctx.full)),
            "symplectic": sp.factor(result["fields"]["symplectic"].coefficient(ctx.full)),
            "constraint_checks": result["momentum"]["checks"], "chart": chart,
            "representation": representation, "normalized": False}
