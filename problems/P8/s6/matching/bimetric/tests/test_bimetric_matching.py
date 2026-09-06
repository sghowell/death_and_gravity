import pytest
import sympy as sp
from p8_bimetric import matching, vacuum
from p8_bimetric.background import MF2, MG2, NU


@pytest.mark.parametrize("check", [vacuum.potential_checks, vacuum.mass_basis_checks,
                                  matching.source_matching_checks, matching.projector_checks,
                                  matching.physical_residue_checks, matching.fierz_pauli_source_checks,
                                  matching.flat_curvature_checks, matching.direct_curvature_checks,
                                  matching.remainder_checks])
def test_exact_vacuum_and_full_source_identities(check):
    assert all(sp.simplify(value) == 0 for value in check().values())


def test_positive_kinetic_mass_and_physical_poles():
    data = vacuum.spectrum()
    for key in ("M_squared", "relative_mode_kinetic_coefficient", "m_FP_squared",
                "massless_source_residue", "massive_source_residue"):
        assert data[key].is_positive
    matter = vacuum.canonical_matter()
    assert matter["clock_kinetic"] == matter["M1_kinetic"] == 1
    assert matter["clock_mass_squared"].is_positive
    assert matter["M1_mass_squared"] == 0
    assert matter["quadratic_metric_matter_mix_at_constant_zero_energy_vacuum"] == 0


def test_exact_equal_planck_example():
    point = {MG2: 1, MF2: 1, NU: 1}
    v, m = vacuum.spectrum(), matching.tt_schur()
    assert v["M_squared"].subs(point) == 2
    assert v["m_FP_squared"].subs(point) == 2
    assert v["massless_source_residue"].subs(point) == sp.Rational(1, 2)
    assert v["massive_source_residue"].subs(point) == sp.Rational(1, 2)
    assert m["c_C"].subs(point) == sp.Rational(1, 4)
    assert m["beta_C"].subs(point) == sp.Rational(1, 8)


@pytest.mark.parametrize("ratio", [-sp.Rational(1, 2), sp.Rational(1, 2)])
def test_kernel_remainder_exact_spectral_disk_examples(ratio):
    data = matching.tt_schur()
    point = {MG2: 2, MF2: 3, NU: 5, matching.D: sp.Rational(5, 3)*ratio}
    error = abs(data["kernel_remainder"].subs(point))
    eta = sp.Rational(1, 2)
    dabs = abs(point[matching.D])
    bound = MF2**3*dabs**3/(NU**2*(1-eta))
    assert error <= bound.subs(point)
    relative = error/abs(((MG2+MF2)*matching.D).subs(point))
    assert relative <= eta**2/(1-eta)


def test_trace_source_detects_non_TT_matching():
    p = matching.projectors()
    trace_source = sp.Matrix([1, 1, 1, 0, 0, 0])
    assert p["P2"]*trace_source == sp.zeros(6, 1)
    assert (trace_source.T*(p["P2"]-p["P0"]/2)*trace_source)[0] == -sp.Rational(3, 2)
    assert matching.tt_schur()["c_R_flat_quadratic_tree"] == 0
    # The negative off-shell constraint entry is not the on-shell residue:
    # the independent null-source calculation gives two positive squares.
    assert all(value == 0 for value in matching.physical_residue_checks().values())


def test_source_and_wrong_field_omission_controls():
    controls = matching.controls()
    assert controls["wrong_field_equation_changes_four_derivative_kernel"] == (MG2+MF2)**2/NU
    assert controls["massless_and_massive_trace_projectors_differ"] == sp.Rational(1, 6)
    assert controls["auxiliary_inverse_pole_not_physical_massive_pole"] == NU/MG2
    assert controls["discarding_massive_source_contact_loses_exchange"] != 0
    assert controls["physical_g_not_the_massless_eigenmode"] != 0
    assert controls["zero_curvature_R2_tree_coefficient_not_radiative_closure"] == 0


def test_outside_class_spectrum_controls():
    controls = vacuum.controls()
    assert controls["negative_spring_is_tachyonic"].is_negative
    assert controls["zero_spring_has_no_heavy_gap"] == 0
    assert controls["negative_second_Einstein_coefficient_is_not_healthy"] < 0
