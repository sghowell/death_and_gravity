"""Exact principal-embedding, both-orientation and actual Proca constraint tests."""
import json

import pytest
import sympy as sp
from p8_affine import verify as original
from p8_proca_spectator_repair import audit, pencil, proca, repair


def zero(value):
    return all(v==0 for v in value) if isinstance(value,sp.MatrixBase) else value==0


def test_exact_named_and_scalar_identity_counts():
    rows=audit.residuals()
    assert len(rows)==31
    assert sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in rows.values())==40
    assert all(zero(value) for value in rows.values())


@pytest.mark.parametrize("name",tuple(audit.residuals()))
def test_each_exact_identity(name):
    assert zero(audit.residuals()[name])


def test_all_actual_gap_and_explicit_scope_gates():
    assert len(audit.gates())==16
    assert all(value is True for value in audit.gates().values())


def test_canonical_comparison_keeps_actual_matter_clock_mixing():
    d=pencil.light()
    K,H,T=(d[key] for key in ("literal_light_kinetic","physical_matter_frame_light_gradient",
        "old_kinetic_canonical_congruence"))
    assert K.has(pencil.modes.ell)
    assert all(H.has(v) for v in (pencil.modes.ell,pencil.c))
    assert T[1,0]!=0
    assert sp.factor(T.det())==1/sp.sqrt(pencil.modes.kc*pencil.modes.km)


def test_all_heavy_and_cross_terms_are_present_before_clock_restriction():
    d=pencil.extension()
    P=d["full_characteristic_pencil_anchor"]
    restricted=d["clock_direction_pencil"]
    heavy={s for s in P.free_symbols if str(s).startswith("heavy_")}
    assert len(heavy)==21
    assert not restricted.has(*heavy)
    assert all(restricted.has(d[key][0,0]) for key in
        ("light_kinetic_change","light_mixed_change","light_gradient_change"))


def test_arbitrary_mixed_term_changes_opposite_orientations_with_opposite_signs():
    d=pencil.extension()
    value=d["clock_direction_pencil"]
    mixed=d["light_mixed_change"][0,0]
    assert sp.diff(value.subs(pencil.s,1),mixed)==2
    assert sp.diff(value.subs(pencil.s,-1),mixed)==-2


@pytest.mark.parametrize("speed",(sp.Rational(1000001,1000000),2,3))
def test_mixed_only_control_repairs_one_direction_and_not_the_other(speed):
    d=repair.data()
    f=d["positive_direction_only_mixed_repair_pencil"].subs(pencil.c,speed)
    assert f.subs(pencil.s,1)==0
    assert f.subs(pencil.s,-1)<0
    assert f.subs(pencil.s,-speed**2)==0
    assert -speed**2<-1


def test_real_light_gradient_change_is_not_excluded_by_screening():
    d=repair.data()
    assert d["algebraic_gradient_change_that_can_repair_principal_pair"][0,0]==1-pencil.c**2
    assert d["repaired_canonical_gradient"]==sp.eye(2)
    assert d["algebraic_repair_not_a_covariant_action_or_UV_construction"] is True


BAD=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
    sp.Symbol("free"),sp.sqrt(2),[],-1)


@pytest.mark.parametrize("bad",BAD)
def test_relative_norm_domains_reject_bad_input_in_both_arguments(bad):
    assert repair.budget()["strict_nonrepair_margin"]>0
    for pair in ((bad,0),(0,bad)):
        with pytest.raises((TypeError,ValueError)):
            repair.budget(*pair)


@pytest.mark.parametrize("pair",((repair.GAP,0),(0,repair.GAP),
    (repair.GAP/2,repair.GAP/2),(repair.GAP,repair.GAP)))
def test_nonstrict_or_negative_nonrepair_margins_rejected(pair):
    with pytest.raises(ValueError):
        repair.budget(*pair)


def test_valid_rational_budget_retains_positive_margin_without_bounding_mixed_terms():
    d=repair.budget(repair.GAP/5,repair.GAP/3)
    assert d["strict_nonrepair_margin"]==7*repair.GAP/15
    assert d["no_smallness_assumption_on_finite_mixed_time_space_or_heavy_blocks"] is True
    assert d["bounds_are_conditional_inputs_not_claimed_actual_UV_matching_errors"] is True


@pytest.mark.parametrize("bad",BAD[:-3]+(sp.Rational(1,2),[],-1))
def test_finite_species_domain_rejected_by_both_entry_points(bad):
    assert proca.species(1)==1
    for call in (proca.species,proca.characteristic):
        with pytest.raises((TypeError,ValueError)):
            call(bad)


@pytest.mark.parametrize("n",(0,1,3,1000))
def test_finite_species_product_keeps_one_fast_clock_factor(n):
    d=proca.characteristic(n)
    assert d["number_of_physical_Proca_polarizations"]==3*n
    assert d["number_of_luminal_physical_modes"]==3+3*n
    assert d["total_number_of_physical_modes"]==4+3*n
    assert d["normalized_principal_characteristic"]==(pencil.c**2-pencil.s**2)*(1-pencil.s**2)**(3+3*n)


def test_literal_vector_action_varies_every_metric_first_jet():
    d=proca.data()
    H=d["inverse_metric_first_jet"]
    assert H==H.T
    assert len(H.free_symbols)==10
    assert all(d["full_relevant_metric_vector_order_action"].has(v) for v in H.free_symbols)
    assert d["metric_volume_first_jet"]!=1


def test_actual_constant_canonical_mass_normalization_not_assumed_unit():
    d=proca.data()
    assert d["actual_literal_unscaled_Maxwell_coefficient"]==sp.Rational(1,10**6)
    assert d["actual_existing_canonical_Proca_mass_squared"]==10**6


def test_longitudinal_constraint_and_finite_mass_are_retained():
    d=proca.data()
    assert all(d["temporal_constraint_solution"].has(v) for v in (proca.k,proca.m))
    assert d["longitudinal_kinetic_coefficient"]==proca.m**2/(proca.k**2+proca.m**2)
    assert d["longitudinal_kinetic_coefficient"].is_positive is True
    assert d["tangent_frame_frozen_massive_dispersion_squared"]==proca.k**2+proca.m**2


def test_all_invalid_exact_domains_are_rejected():
    d=audit.controls()
    assert d["groups"]["relative_norm_budget"]["rejected_inputs"]==34
    assert d["groups"]["finite_species"]["rejected_inputs"]==30
    assert d["rejected_inputs"]==64


def test_all_exact_report_data_serialize_without_float_rounding():
    blocks=(pencil.light(),pencil.extension(),repair.data(),proca.data())
    serialized=json.loads(json.dumps(original.serialize(blocks)))
    assert len(serialized)==4
    assert serialized[0]["actual_clock_squared_speed_excess_lower"]==str(repair.GAP)
