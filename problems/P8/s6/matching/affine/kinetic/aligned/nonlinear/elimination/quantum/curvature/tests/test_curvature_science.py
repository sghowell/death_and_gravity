"""Direct tensor, physical geometry and independent spectral controls."""
import pytest
import sympy as sp
from p8_vector_curvature import geometry, heat, sphere


def test_exact_curvature_and_physical_metric_identity_chains():
    for module in (heat, geometry, sphere):
        for value in module.checks().values():
            assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_general_four_dimensional_curvature_not_just_constant_curvature():
    assert heat.explicit_tensor_traces()["curvature_parameter_count"] == 20


def test_full_clock_reduction_does_not_erase_off_clock_mass_variations():
    from p8_aligned_quantum import potential
    assert all(value == 0 for value in heat.clock_reduction().values())
    assert potential.clock_jets()["pole_weight"]["N_first"] != 0


def test_finite_contact_determinant_cannot_be_dropped_without_a_prescription():
    d = heat.cochain_control()
    missing_contact = sp.factor(d["K"].det()*(d["Delta0"]+heat.m2*sp.eye(3)).det()
                               -(d["Delta1"]+heat.m2*sp.eye(4)).det())
    assert missing_contact.subs(heat.m2, 2) != 0


def test_massive_scalar_subtraction_is_not_massless_Maxwell_or_three_scalars():
    d = heat.coefficients()
    assert d["proca"][0] == 3
    assert sp.expand(d["proca"][2]-3*d["scalar"][2]) != 0
    assert sp.expand(d["proca"][4]-(d["vector"][4]-2*d["scalar"][4])) != 0


def test_actual_rolling_curvature_is_not_de_sitter():
    d, u = geometry.rolling(), geometry.u
    assert d["R"].subs(u, 0) == 24
    assert d["R"].subs(u, sp.sqrt(sp.Rational(5, 7))) == 49
    assert sp.diff(d["R"], u) != 0
    assert sp.limit(d["R"], u, sp.oo) == 0
    assert d["boxR"] != 0


@pytest.mark.parametrize("point", (0, sp.Rational(1, 2), 1, 2, 10))
def test_continuous_envelopes_have_exact_point_controls(point):
    d = geometry.rolling()
    for name, bound in geometry.continuous_bounds().items():
        assert abs(d[name].subs(geometry.u, point)) <= bound["absolute_upper"]


def test_sphere_constant_mode_must_not_be_silently_dropped():
    d = sphere.spectral()
    assert d["leading"].coeff(d["t"], 0) == -sp.Rational(11, 30)
    coexact_only = d["leading"]+1
    assert coexact_only.coeff(d["t"], 0) == sp.Rational(19, 30)
    assert d["trace_remainder_over_t_upper_for_0_lt_t_le_1"] > 0
    assert all(value is True for value in sphere.proof_checks().values())


def test_scale_bound_is_local_pole_not_a_finite_loop_claim():
    d = geometry.scale_bounds(1000)
    assert d["scalar_curvature_over_mass_squared"] == sp.Rational(49, 10**6)
    assert d["finite_curved_or_in_in_remainder"] is False


@pytest.mark.parametrize("bad", (True, 1.0, "1", sp.oo, sp.I, 0, -1))
def test_invalid_mass_scale_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        geometry.scale_bounds(bad)
