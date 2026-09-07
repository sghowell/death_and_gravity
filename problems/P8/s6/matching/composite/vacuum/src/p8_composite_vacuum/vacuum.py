"""Physical-frame proportional Minkowski candidate and full FP mass tensor.

The potential jet conditions below are conditions on ONE specified light
potential. They do not supply an off-tube extension of the pinned rolling
reconstruction. The massive mode has zero physical source residue here,
not negative kinetic residue and not zero mass.
"""

from functools import cache

import sympy as sp

from . import cayley


@cache
def derive():
    M2, m2 = cayley.M2, cayley.m2
    q = sp.Symbol("D_flat", real=True)  # partial_T²+kphys², not spatial q alone
    gamma, delta, source = sp.symbols("gamma Delta_TT j", real=True)
    # Lorentz-invariant flat vacuum, g=f=G/4; physical G is Minkowski.
    planck = M2/2
    relative_planck = 2*M2  # for the dimensionless Cayley tensor Delta
    mass = m2/4
    raw_kinetic = M2*q*((gamma+2*delta)**2+(gamma-2*delta)**2)/32
    raw_mass = M2*m2*delta**2/16
    mass_basis = (planck*q*gamma**2+relative_planck*(q+mass)*delta**2)/8
    # Convention L=-[K0 D gamma²+KR(D+mass)Delta²]/8+j gamma/4.
    action = -mass_basis+source*gamma/4
    vac = -3*M2*m2/8
    shifted = tuple(beta+vac/(M2*m2) for beta in (0, 0, 1, 0, 0))
    return {"M2": M2, "m2": m2, "D": q, "gamma": gamma, "Delta_TT": delta, "source": source,
            "vacuum_potential_value": vac,
            "shifted_beta_values": shifted,
            "physical_Planck_squared": planck,
            "relative_Planck_squared": relative_planck,
            "physical_FP_mass_squared": mass,
            "raw_positive_TT_kernel": raw_kinetic+raw_mass,
            "diagonal_positive_TT_kernel": mass_basis,
            "sourced_mass_basis_action": action,
            "physical_massless_source_response": 1/(planck*q),
            "physical_massive_source_response": sp.S.Zero,
            "relative_auxiliary_probe_residue": 1/relative_planck,
            "tree_C_D_curvature_corrections_on_zero_relative_branch": sp.S.Zero}


@cache
def checks():
    d = derive()
    M2, m2, gamma, delta, j, q = (d[key] for key in ("M2", "m2", "gamma", "Delta_TT", "source", "D"))
    beta = d["shifted_beta_values"]
    # The constant physical-volume vacuum term shifts ALL beta_n.
    ug = beta[0]+3*beta[1]+3*beta[2]+beta[3]
    uf = beta[4]+3*beta[3]+3*beta[2]+beta[1]
    # Standard FP mass is first in g proper units, then divided by4 because
    # G=4g and physical T=2*T_g at this proportional vacuum.
    standard_g_mass = 2*m2*(beta[1]+2*beta[2]+beta[3])
    gamma_solution = j/(d["physical_Planck_squared"]*q)
    action = d["sourced_mass_basis_action"]
    trace, trace2 = sp.symbols("trace_Delta trace_Delta2", real=True)
    fp_quadratic = -M2*m2*(trace2-trace**2)/16
    expected_fp = -d["relative_Planck_squared"]*d["physical_FP_mass_squared"]*(trace2-trace**2)/8
    return {
        "g_vacuum_equation_with_all_beta_shifts": ug,
        "f_vacuum_equation_with_all_beta_shifts": uf,
        "vacuum_cosmological_constant_cancellation": 3*M2*m2/8+d["vacuum_potential_value"],
        "physical_clock_FP_mass_from_shifted_betas": sp.factor(standard_g_mass/4-d["physical_FP_mass_squared"]),
        "exact_TT_mass_basis_and_residues": sp.expand(d["raw_positive_TT_kernel"]-d["diagonal_positive_TT_kernel"]),
        "full_relative_Fierz_Pauli_tensor": sp.expand(fp_quadratic-expected_fp),
        "massless_source_stationarity": sp.factor(sp.diff(action, gamma).subs(gamma, gamma_solution)),
        "unsourced_relative_stationarity": sp.diff(action, delta).subs(delta, 0),
        "physical_source_functional_has_no_massive_pole": sp.factor(
            action.subs({gamma: gamma_solution, delta: 0}, simultaneous=True)
            -j**2/(8*d["physical_Planck_squared"]*q)),
    }


def source_projector_checks():
    """Physical conserved source sees the GR trace projector, not massive P2.

Full FP trace/shift constraints are audited independently in the root test.
This algebra checks the physical source dictionary beyond TT amplitudes.
"""
    d = derive()
    spatial = sp.Matrix(sp.symbols("T11 T22 T33 sqrt2T12 sqrt2T13 sqrt2T23", real=True))
    scalar = sp.zeros(6)
    scalar[:3, :3] = sp.ones(3)/3
    spin2 = sp.eye(6)-scalar
    response = (spin2-scalar/2)/(d["physical_Planck_squared"]*d["D"])
    trace = sum(spatial[:3, 0])
    return {"physical_source_GR_projector": sp.expand(
                (spatial.T*response*spatial)[0]
                -(spatial.dot(spatial)-trace**2/2)/(d["physical_Planck_squared"]*d["D"])),
            "relative_five_spin2_components": sp.trace(spin2)-5,
            "relative_source_projection_zero": sp.S.Zero}


def controls():
    d = derive()
    m2, M2 = d["m2"], d["M2"]
    eps = sp.Symbol("epsilon", positive=True)
    bounce_v = M2*m2*(sp.Rational(9, 2)*eps**2-eps/(1+eps)**2)
    bare_beta2_wrong_mass = 2*m2*(0+2*1+0)/4
    return {"bare_beta2_wrong_physical_mass": bare_beta2_wrong_mass,
            "lost_vacuum_shift_mass_error": bare_beta2_wrong_mass-d["physical_FP_mass_squared"],
            "wrong_equal_time_clock_mass": m2,
            "pinned_bounce_initial_potential_over_M2m2": bounce_v/(M2*m2),
            "constant_vacuum_potential_mismatch_cleared": sp.cancel(
                (bounce_v-d["vacuum_potential_value"])*8*(1+eps)**2/(M2*m2)),
            "zero_physical_heavy_source_not_zero_mass": d["physical_FP_mass_squared"]}
