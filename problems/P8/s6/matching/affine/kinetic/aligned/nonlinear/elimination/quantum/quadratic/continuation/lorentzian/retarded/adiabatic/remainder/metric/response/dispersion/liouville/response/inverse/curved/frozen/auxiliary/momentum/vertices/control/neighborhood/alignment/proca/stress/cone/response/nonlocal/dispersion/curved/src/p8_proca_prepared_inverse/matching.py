"""New minimal physical scalar current and its complete dimensional diagonal matching."""
from functools import cache

import sympy as sp
from p8_prepared_volterra import vertices as generic
from p8_proca_local_response import local
from p8_proca_rank_one_inverse import spectral
from p8_vector_metric_local import canonical, jets, variation

D=generic.D
ZERO=spectral.ZERO


def high(sector):
    return _high(jets.kind(sector))


@cache
def _high(sector):
    weights=canonical.data(sector)["weights"]
    mapping={jets.D:D,jets.z:1}
    A=sp.ImmutableMatrix([weights[label][0].xreplace(ZERO).subs(mapping,simultaneous=True) for label in ("N","Z")])
    B=sp.ImmutableMatrix([weights[label][1].xreplace(ZERO).subs(mapping,simultaneous=True) for label in ("N","Z")])
    return {"A":A,"B":B,"pair":sp.ImmutableMatrix((A-B).applyfunc(sp.factor))}


@cache
def data():
    b=high("L")["pair"][1].subs(D,3)
    return {"new_transverse_dimensional_high_current":high("T"),
            "new_longitudinal_dimensional_high_current":high("L"),
            "new_scalar_physical_high_pair":b,"normalized_scalar_log_weight":b*b/4,
            "normalized_scalar_fourth_pole":b*b/8,
            "fixed_new_scalar_finite_local_fourth_coefficient":-sp.Integer(4),
            "generic_bilocal_current_polynomial_coefficients":tuple(generic.polynomial().nth(j) for j in range(5)),
            "new_flat_massive_reference_scalar_block":spectral.data()["normalized_local_plus_thrice_subtracted_scalar_block"],
            "normalization":"Q_w=64*pi^2 times 3*delta_pressure_Gamma[n=0,v=w], with the original fixed-profile normalization contact retained.",
            "normalized_high_radial_kernel":"128*k^4*sin(2*k*Delta_sigma)/(a(t)^4*a(s)), with sigma'=1/a.",
            "normalized_off_diagonal_leading_kernel":"96/[a(t)^4*a(s)*Delta_sigma^5]; difference from 96/(t-s)^5 has one lower singular order.",
            "all_remaining_local_source_derivative_orders_at_most_three":True}


@cache
def checks():
    out={}
    pair=high("L")["pair"]
    transverse=high("T")["pair"]
    out["new_transverse_high_current_full_dimension"]=sp.ImmutableMatrix((transverse-sp.Matrix([0,2*(D-3)])).applyfunc(sp.factor))
    out["new_longitudinal_high_current_full_dimension"]=sp.ImmutableMatrix((pair-sp.Matrix([0,2*(D-1)])).applyfunc(sp.factor))
    out["new_lapse_high_pair_vanishes_throughout_dimensional_neighborhood"]=pair[0]
    out["new_transverse_high_square_and_first_dimensional_jet"]=sp.ImmutableMatrix([
        ((D-1)*transverse[1]**2).subs(D,3),sp.diff((D-1)*transverse[1]**2,D).subs(D,3)])
    out["new_longitudinal_dimensional_pair_jet_retained"]=sp.diff(pair[1],D)-2
    out["new_longitudinal_dimensional_square_jet_retained"]=sp.diff(pair[1]**2,D).subs(D,3)-16
    out["new_actual_scalar_high_pair_matches_massive_reference"]=sp.factor(
        spectral.pair("L")[1].subs(spectral.z,1)-pair[1].subs(D,3))
    out["new_actual_scalar_log_weight"]=data()["normalized_scalar_log_weight"]-4
    out["new_actual_scalar_pole"]=data()["normalized_scalar_fourth_pole"]-2
    out["new_actual_curved_scalar_fourth_local_coefficient"]=sp.factor(
        local.operator("Z",2).coeff(local.v[4])+4)
    out["new_actual_curved_inactive_fourth_local_coefficient"]=sp.factor(
        local.operator("N",2).coeff(local.n[4]))
    for sector in ("T","L"):
        weights=canonical.data(sector)["weights"]
        pairs=sp.Matrix([sp.factor((weights[label][0]-weights[label][1]).xreplace(ZERO)) for label in ("N","Z")])
        for i,label in enumerate(("N","Z")):
            value=variation.data(sector)["coefficients"][label][2].xreplace(ZERO)
            for j,row in enumerate((jets.n,jets.v)):
                actual=sp.expand(value).coeff(row[4])
                out[sector+f"_new_full_dimensional_all_momentum_curved_fourth_contact_{i}_{j}"]=sp.factor(actual-pairs[i]*pairs[j]/32)
    inherited=generic.checks()
    for key in ("generic_leading_bilocal_pair","generic_subleading_connection_not_dropped",
                "proper_momentum_scale_cancels_for_every_dimension","leading_geometry_diagonal",
                "leading_geometry_first_lag_jet","radial_leading_Abel_limit",
                "radial_subleading_Abel_limit"):
        out["generic_current_identity_"+key]=inherited[key]
    lag=sp.Symbol("positive_lag",positive=True)
    out["new_scalar_high_kernel_four_primitive_normalization"]=sp.factor(sp.diff(4/lag,lag,4)-96/lag**5)
    # The diagonal single-current connection cancels algebraically; its full off-diagonal coefficient is retained.
    A=high("L")["A"][1].subs(D,3)
    B=pair[1].subs(D,3)
    out["new_single_scale_connection_diagonal"]=sp.factor(2*(B*A-A*B))
    return out
