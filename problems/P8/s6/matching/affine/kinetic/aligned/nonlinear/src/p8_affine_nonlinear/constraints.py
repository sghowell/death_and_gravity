"""Auxiliary block rank and spatial first-class count on a local regular branch."""
from functools import cache

import sympy as sp
from p8_affine_aligned import alignment

from . import adm, trace


@cache
def boundaries():
    density, I, Is, Iphi, ds, N, K = sp.symbols(
        "sqrt_hat I I_s I_phi advected_s N K_hat", real=True)
    # The spacetime divergence divided by sqrt_hat is
    # I_s*(s_dot-N^i D_i s)+I_phi+N*K_hat*I.
    divergence = density*(Is*ds+Iphi+N*K*I)
    replacement = -density*Iphi-N*density*K*I
    return {"full_linear_lapse_velocity_boundary": sp.expand(density*Is*ds-divergence-replacement)}


@cache
def block_inverse():
    D, E = sp.MatrixSymbol("D", 2, 2), sp.MatrixSymbol("E", 2, 2)
    zero = sp.ZeroMatrix(2, 2)
    matrix = sp.BlockMatrix([[zero, -D], [D.T, E]])
    inverse = sp.BlockMatrix([[D.T.inv()*E*D.inv(), D.T.inv()], [-D.inv(), zero]])
    # Materialize the collapsed product before subtracting the identity.
    # Re-collapsing an unevaluated block-diagonal identity minus Identity(4)
    # gives a wrong 2I-I result in the installed SymPy; keep an exact control.
    residual = sp.block_collapse(matrix*inverse).as_explicit()-sp.eye(4)
    return {"residual": residual,
            "does_not_require_an_inverse_of_secondary_secondary_operator": True}


@cache
def proof_checks():
    bg = alignment.parent.old.background()
    numerator, denominator = sp.fraction(bg["J"])
    poly = sp.Poly(numerator, adm.u)
    lower = sp.Rational(1199, 800)/(sp.Rational(5, 4)**18)
    gamma_correction = 3*sp.Rational(1, 100)/(8*sp.Rational(9, 40))
    gamma_lower = sp.Rational(18, 19)-gamma_correction
    return {"temporal_mass_correction_upper": gamma_correction == sp.Rational(1, 60),
            "effective_temporal_gamma_lower": gamma_lower == sp.Rational(1061, 1140),
            "effective_temporal_gamma_above_nine_tenths": bool(gamma_lower > sp.Rational(9, 10)),
            "positive_even_J_numerator": bool(poly.TC() == 1199 and all(
                powers[0] % 2 == 0 and coefficient > 0 for powers, coefficient in poly.terms())),
            "actual_J_denominator": sp.factor(denominator-800*(1+adm.u**2)**18) == 0,
            "compact_J_lower_above_one_over_40": bool(lower > sp.Rational(1, 40)),
            "joint_clock_auxiliary_Jacobian_exact": trace.clock_auxiliary()["actual_joint_auxiliary_Jacobian"] == sp.zeros(2),
            "uniform_clock_secondary_pivot_above_one_over_20": bool(2*lower > sp.Rational(1, 20)),
            "auxiliary_constraint_block_has_local_inverse": block_inverse()["residual"] == sp.zeros(4),
            "complete_timelike_phase_count_with_free_matter": sp.Rational(2*(6+1+3+4+1)-2*6-4, 2) == 7,
            "count_not_nonlinear_energy_or_global_wellposedness": True}


@cache
def checks():
    out = dict(boundaries())
    # The abstract inverse is an operator-composition identity, not a
    # commuting determinant argument for differential secondary brackets.
    out["four_auxiliary_constraint_block_inverse"] = sp.Matrix(block_inverse()["residual"])
    return out
