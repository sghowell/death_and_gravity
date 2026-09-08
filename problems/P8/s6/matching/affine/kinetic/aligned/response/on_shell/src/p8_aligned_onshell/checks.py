"""Physical readout and moving canonical-boundary consistency checks."""
from functools import cache

import sympy as sp

from . import phase


@cache
def canonical_swap():
    data = phase.system()
    old, bg = phase.old, phase.old.background()
    pv = sp.Symbol("p_v", real=True)
    original_states = sp.Matrix([old.v, old.matter, pv, phase.PS])
    c = 2*data["a"]**3*data["q"]
    # Recover the unswapped density before adding the moving boundary.
    original_H = (data["H"]+bg["H"]*old.shift*phase.PB).subs(
        {old.shift: -pv/c, phase.PB: c*old.v}, simultaneous=True)
    original_A = phase.clean(data["symplectic"]*sp.hessian(original_H, tuple(original_states)))
    mapping = sp.Matrix([-pv/c, old.matter, c*old.v, phase.PS])
    transform = mapping.jacobian(original_states)
    actual = phase.clean((transform.diff(phase.u)+transform*original_A)*transform.inv())
    wrong = phase.clean(transform*original_A*transform.inv())
    return {"residual": phase.clean(actual-data["A"]),
            "omitted_time_boundary_residual": phase.clean(wrong-data["A"])}


@cache
def residuals():
    data = phase.system()
    old, bg = phase.old, phase.old.background()
    state = phase.STATES
    flow = data["A"]*state

    def derivative(value):
        return sp.diff(value, phase.u)+sum(sp.diff(value, item)*rate for item, rate in zip(state, flow, strict=True))

    N, eps = sp.symbols("N eps", positive=True)
    omega = -sp.log(1+(N**-2-1)/bg["h"])/4
    literal_K = (3*bg["H"]+eps*(3*old.vd+data["q"]*old.shift))/(1+eps*old.n)
    linear_K = sp.diff(literal_K, eps).subs(eps, 0)
    expected = linear_K.subs({old.vd: bg["theta"]*old.n-bg["ell"]*old.matter/2,
                             old.n: data["lapse"]}, simultaneous=True).subs(old.n, data["lapse"])
    jets = phase.initial_jets()
    return {"actual_shift_constraint_after_moving_boundary": phase.clean(
                derivative(data["physical_v"])-bg["theta"]*data["lapse"]+bg["ell"]*old.matter/2),
            "actual_matter_momentum_readout": phase.clean(flow[1]-phase.PS/data["a"]**3+bg["w"]*data["lapse"]),
            "spatial_chart_metric_map": phase.clean(sp.diff(omega, N).subs(N, 1)-bg["delta"]),
            "ADM_trace_readout": phase.clean(expected-data["trace"]),
            "prepared_initial_lapse_jets": jets["preparation_residual"],
            "independent_original_canonical_swap": canonical_swap()["residual"],
            "canonical_phase_symplectic_identity": phase.clean(
                data["A"]*data["symplectic"]+data["symplectic"]*data["A"].T)}
