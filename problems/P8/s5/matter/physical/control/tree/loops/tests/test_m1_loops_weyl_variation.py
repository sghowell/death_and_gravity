"""Exact tensor contractions showing why background Weyl-flatness is insufficient."""

from fractions import Fraction as F
from itertools import product


def linearized_invariants(wave_covector):
    """Use h_11=f, h_22=-f and factor out the common f'' amplitude."""
    metric = (1, -1, -1, -1)
    perturbation = (0, 1, -1, 0)

    def derivative(i, j, a, b):
        return wave_covector[i]*wave_covector[j]*(perturbation[a] if a == b else 0)

    riemann = {}
    for a, b, c, d in product(range(4), repeat=4):
        riemann[a, b, c, d] = F(
            derivative(c, b, a, d)+derivative(d, a, c, b)
            - derivative(c, a, d, b)-derivative(d, b, a, c), 2)
    ricci = {(b, d): sum(metric[a]*riemann[a, b, a, d] for a in range(4))
             for b, d in product(range(4), repeat=2)}
    scalar = sum(metric[a]*ricci[a, a] for a in range(4))
    ricci_squared = sum(metric[a]*metric[b]*value**2 for (a, b), value in ricci.items())
    riemann_squared = sum(metric[a]*metric[b]*metric[c]*metric[d]*value**2
                          for (a, b, c, d), value in riemann.items())
    return {
        "scalar": scalar,
        "ricci_squared": ricci_squared,
        "riemann_squared": riemann_squared,
        "weyl_squared": riemann_squared-2*ricci_squared+scalar**2/3,
        "euler": riemann_squared-4*ricci_squared+scalar**2,
        "nonzero_curvature_components": sum(value != 0 for value in riemann.values()),
    }


def test_static_tt_wave_has_nonzero_weyl_second_variation():
    result = linearized_invariants((0, 0, 0, 1))
    assert result["scalar"] == 0
    assert result["ricci_squared"] == F(1, 2)
    assert result["riemann_squared"] == 2
    assert result["weyl_squared"] == 1
    assert result["euler"] == 0
    full = (result["riemann_squared"]-result["ricci_squared"])/180+result["scalar"]**2/72
    wrongly_drop_weyl = -result["euler"]/360+result["scalar"]**2/72
    assert full-wrongly_drop_weyl == F(1, 120)


def test_null_wave_zero_scalar_invariants_do_not_mean_zero_tensor():
    result = linearized_invariants((1, 0, 0, -1))
    assert result["scalar"] == 0
    assert result["ricci_squared"] == 0
    assert result["riemann_squared"] == 0
    assert result["weyl_squared"] == 0
    assert result["euler"] == 0
    assert result["nonzero_curvature_components"] > 0
