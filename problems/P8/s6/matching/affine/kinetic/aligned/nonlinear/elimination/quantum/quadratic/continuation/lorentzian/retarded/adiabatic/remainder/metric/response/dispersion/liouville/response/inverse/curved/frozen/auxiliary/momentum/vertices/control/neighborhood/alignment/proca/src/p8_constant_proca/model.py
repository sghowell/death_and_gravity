"""Actual new constant-mass Hamiltonian and full quadratic coefficient bridge."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model as old
from p8_zero_source import coefficients as previous_coefficients
from p8_zero_source import model as previous


@cache
def data():
    g=old.generic()
    before=previous.data()["new_Hamiltonian"]
    after=before.subs({g["gt"]:1,g["gs"]:1},simultaneous=True)
    delta=old.N*((1-g["gt"])*old.j**2/(2*g["U"])
                 +g["eo"]*(1-1/g["gs"])*old.vector/2)
    return {"new_Hamiltonian":after,"previous_source_free_Hamiltonian":before,
            "new_minus_previous_Hamiltonian":delta,
            "new_temporal_constraint_solution":-old.j/g["U"],
            "new_temporal_Hessian":-old.N*g["U"],
            "new_source_free_mass_matrix":sp.diag(1,-1,-1,-1),
            "canonical_Proca_mass_squared":sp.Integer(10**6)}


def powers(index,degree):
    out=[0]*len(old.VARIABLES)
    out[index]=degree
    return tuple(out)


@cache
def jets():
    result=dict(previous_coefficients.jets())
    c=old.coefficients()
    changed={powers(old.VARIABLES.index(old.j),2):old.N/(2*c["U"]),
             powers(old.VARIABLES.index(old.vector),1):old.N*c["e_omega"]/2}
    for key,value in changed.items():
        result[key]=tuple(sp.factor(sp.diff(value,old.N,n).subs(old.N,1)) for n in range(5))
    return result


@cache
def checks():
    d=data()
    c=old.coefficients()
    before,after=previous_coefficients.jets(),jets()
    out={"literal_constant_mass_Hamiltonian_difference":sp.factor(
        d["new_Hamiltonian"]-d["previous_source_free_Hamiltonian"]-d["new_minus_previous_Hamiltonian"]),
        "unchanged_all_homogeneous_light_Hamiltonians":sp.factor(d["new_minus_previous_Hamiltonian"].subs({old.j:0,old.vector:0})),
        "new_complete_Hamiltonian_vector_parity":sp.factor(d["new_Hamiltonian"].subs(old.j,-old.j)-d["new_Hamiltonian"]),
        "actual_old_temporal_mass_at_every_clock_point":sp.factor(c["gamma_t"].subs(old.N,1)-1),
        "actual_old_spatial_mass_at_every_clock_point":sp.factor(c["gamma_s"].subs(old.N,1)-1)}
    for key,row in before.items():
        name="_".join(map(str,key))
        out["full_clock_invariant_coefficient_bridge_"+name]=sp.factor(after[key][0]-row[0])
        # Only the first-order invariants enter the linear lapse force.
        # The squared vector invariant starts at physical degree two;
        # its changed N derivative is intentionally NOT asserted equal.
        if sum(key)==1 and key[old.VARIABLES.index(old.vector)]==0:
            out["clock_linear_lapse_force_bridge_"+name]=sp.factor(after[key][1]-row[1])
    zero=(0,)*len(old.VARIABLES)
    for n in range(5):
        out[f"unchanged_full_background_lapse_jet_{n}"]=sp.factor(after[zero][n]-before[zero][n])
    polynomial=sp.Poly(sp.expand(d["new_Hamiltonian"]),*old.VARIABLES)
    actual_Gauss=polynomial.coeff_monomial(old.j**2).subs(old.specialization(),simultaneous=True)
    actual_vector=polynomial.coeff_monomial(old.vector).subs(old.specialization(),simultaneous=True)
    for n in range(5):
        out[f"literal_new_Gauss_coefficient_derivative_{n}"]=sp.factor(
            after[powers(2,2)][n]-sp.diff(actual_Gauss,old.N,n).subs(old.N,1))
        out[f"literal_new_spatial_vector_coefficient_derivative_{n}"]=sp.factor(
            after[powers(6,1)][n]-sp.diff(actual_vector,old.N,n).subs(old.N,1))
    return out
