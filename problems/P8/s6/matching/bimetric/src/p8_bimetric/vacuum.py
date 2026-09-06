"""Minkowski potential expansion, mass eigenstates and positive residues."""

import sympy as sp

from .background import MF2, MG2, NU

D = sp.Symbol("D", real=True)  # Flat Lorentzian D=partial_t^2+k^2, not energy.
HG, HF, H0, HM, SOURCE = sp.symbols("h_g h_f h_0 h_m j", real=True)
MPHI2 = sp.Symbol("m_phi_squared", positive=True)


def spectrum():
    total = MG2+MF2
    return {"M_squared": total, "alpha_squared": MF2/MG2,
            "relative_mode_kinetic_coefficient": MG2*MF2/total,
            "m_FP_squared": NU*(1/MG2+1/MF2),
            "massless_source_residue": 1/total,
            "massive_source_residue": MF2/(MG2*total),
            "massless_coordinate": (MG2*HG+MF2*HF)/total,
            "relative_coordinate": HF-HG,
            "physical_metric_from_modes": H0-MF2*HM/total,
            "other_metric_from_modes": H0+MG2*HM/total}


def potential_checks():
    """General matrix trace expansion, not just a TT mass fit.

    G=eta^-1*h_g,F=eta^-1*h_f. Cyclic trace handles noncommuting matrices.
    sqrt(det(I+eG)) and sqrt((I+eG)^-1(I+eF)) are expanded to e^2.
    """
    eps = sp.Symbol("epsilon")
    tg, tf, tgg, tgf, tff = sp.symbols("trG trF trG2 trGF trF2", real=True)
    volg = 1+eps*tg/2+eps**2*(tg**2/8-tgg/4)
    volf = 1+eps*tf/2+eps**2*(tf**2/8-tff/4)
    trace_root = 4+eps*(tf-tg)/2+eps**2*((tgg-tgf)/2-(tff-2*tgf+tgg)/8)
    # sqrt(g)*e4(sqrt(g^-1*f))=sqrt(f) on the chosen positive branch.
    potential = sp.expand(-3*volg+volg*trace_root-volf)
    fp = (tff-2*tgf+tgg-(tf-tg)**2)/8
    return {"Minkowski_zero_potential": potential.coeff(eps, 0),
            "Minkowski_no_tadpole": potential.coeff(eps, 1),
            "full_Fierz_Pauli_relative_mass": sp.expand(potential.coeff(eps, 2)-fp),
            "unit_TT_mass_normalization": fp.subs({tg: 0, tf: 0, tgg: HG**2, tgf: HG*HF, tff: HF**2})-(HG-HF)**2/8}


def mass_basis_checks():
    data = spectrum()
    total, reduced = data["M_squared"], data["relative_mode_kinetic_coefficient"]
    mapping = {HG: data["physical_metric_from_modes"], HF: data["other_metric_from_modes"]}
    quadratic = MG2*D*HG**2+MF2*D*HF**2+NU*(HG-HF)**2
    diagonal = total*D*H0**2+reduced*(D+data["m_FP_squared"])*HM**2
    return {"positive_mass_basis_diagonalization": sp.factor(quadratic.subs(mapping, simultaneous=True)-diagonal),
            "physical_source_couples_to_both_modes": sp.factor((SOURCE*HG).subs(mapping)-SOURCE*H0+SOURCE*MF2*HM/total),
            "kinetic_weighted_mass_eigenvalues": sp.factor(sp.det(sp.Matrix([[MG2*D+NU, -NU], [-NU, MF2*D+NU]]))-MG2*MF2*D*(D+data["m_FP_squared"])),
            "source_residues_sum": sp.factor(data["massless_source_residue"]+data["massive_source_residue"]-1/MG2)}


def canonical_matter():
    """Explicit parent matter: massive canonical clock plus free M1 chi.

    Both fields are at constant vacuum values, V=V_phi=0. Their first-order
    stress then vanishes, so there is no vacuum quadratic metric-matter mix.
    This says nothing about the existence of the required rolling CD state.
    """
    return {"clock_kinetic": sp.Integer(1), "clock_mass_squared": MPHI2,
            "M1_kinetic": sp.Integer(1), "M1_mass_squared": sp.Integer(0),
            "potential": "V(phi)=m_phi_squared*(phi-phi_v)^2/2",
            "quadratic_metric_matter_mix_at_constant_zero_energy_vacuum": sp.Integer(0)}


def controls():
    data = spectrum()
    return {"negative_spring_is_tachyonic": data["m_FP_squared"].subs(NU, -1),
            "zero_spring_has_no_heavy_gap": data["m_FP_squared"].subs(NU, 0),
            "negative_second_Einstein_coefficient_is_not_healthy": sp.Integer(-1),
            "physical_metric_not_massless_eigenfield": -MF2*HM/(MG2+MF2)}
