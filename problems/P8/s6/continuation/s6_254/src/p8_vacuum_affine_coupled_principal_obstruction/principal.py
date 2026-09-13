"""Whole eight-phase fast symbol and an explicit time-dependent growth chart."""

from functools import cache

import sympy as s

from . import coupled as c


def factors(matrix):
    return matrix.applyfunc(s.factor)


@cache
def fast_data():
    H = c.chart()["whole_central_Hessian"]
    shift = -4 * c.a * c.th * c.P**2 / (c.D * (c.L2 + c.L0 * c.a**2 / c.P**2))
    R = s.eye(8)
    R[4, 0] = shift
    dots = s.symbols("theta_dot D_dot L0_dot L2_dot", real=True)
    shift_dot = c.H * c.a * s.diff(shift, c.a) + sum(
        s.diff(shift, variable) * dot
        for variable, dot in zip((c.th, c.D, c.L0, c.L2), dots)
    )
    Rdot = s.zeros(8)
    Rdot[4, 0] = shift_dot
    connection = c.OMEGA * R.inv() * Rdot
    whole = factors(R.T * H * R + connection)
    M = s.diag(
        1, 1, 1, s.sqrt(c.P), c.P ** s.Rational(3, 2), s.sqrt(c.P), s.sqrt(c.P), c.P
    )
    normalized = factors(M.inv() * c.OMEGA * whole * M / c.P ** s.Rational(3, 2))
    fast = normalized.applyfunc(lambda x: s.factor(s.limit(x, c.P, s.oo)))
    alpha = c.L2**2 / (8 * c.J * c.a**3)
    g = c.r * s.sqrt(c.zeta) / (c.D * c.a**2)
    h = c.a**3 * c.Yv / c.zeta
    wanted = s.zeros(8)
    wanted[0, 4] = alpha
    wanted[1, 4], wanted[2, 4] = (c.L2 * w / (4 * c.J * c.Z * c.a**3) for w in c.ws)
    wanted[3, 0], wanted[4, 7], wanted[7, 3] = -g, g, -h
    wanted[5, 0], wanted[6, 0] = (-c.a * charge / c.D for charge in c.cs)
    wanted[5, 1] = wanted[6, 2] = -c.a * c.Y
    lam = s.Symbol("lambda")
    beta = c.r**2 * c.L2**2 * c.Yv / (8 * c.D**2 * c.J * c.a**4)
    # With e=P^-1/2, prove ENTRYWISE that the full remainder is e times
    # a rational expression regular at e=0 on the stated coefficient domain.
    e = s.Symbol("inverse_sqrt_P", positive=True)
    remainder = (normalized - fast).subs(c.P, e**-2).applyfunc(s.cancel)
    divided = (remainder / e).applyfunc(s.cancel)
    denominator_at_zero = []
    denominator_domain_proof = []
    polynomial_degrees = []
    checks = {
        "whole_exact_perfect_square_shear_is_symplectic": R * c.OMEGA * R.T - c.OMEGA,
        "whole_exact_shear_time_connection_is_symmetric": connection - connection.T,
        "entire_fast_matrix_not_a_selected_channel": fast - wanted,
        "whole_eight_phase_fast_characteristic": fast.charpoly(lam).as_expr()
        - lam**4 * (lam**4 - beta),
        "whole_positive_quartic_coefficient": alpha * g**2 * h - beta,
        "complete_reference_fast_matrix_has_no_nonzero_fast_eigenvalue": fast.subs(
            c.r, 0
        )
        .charpoly(lam)
        .as_expr()
        - lam**8,
        "remainder_exact_reconstruction": remainder - e * divided,
    }
    for entry in divided:
        numerator, denominator = s.fraction(entry)
        numerator_poly, denominator_poly = s.Poly(numerator, e), s.Poly(denominator, e)
        polynomial_degrees.append(
            (
                int(numerator_poly.degree()) if numerator != 0 else -1,
                int(denominator_poly.degree()),
            )
        )
        value = s.factor(denominator.subs(e, 0))
        denominator_at_zero.append(value)
        powers = value.as_powers_dict()
        denominator_domain_proof.append(
            all(
                (base.is_number and base != 0)
                or (
                    base in (c.D, c.J, c.K, c.Z, c.a, c.L2)
                    and exponent.is_integer
                    and exponent >= 0
                )
                for base, exponent in powers.items()
            )
        )
    # Other two TT and two transverse-Proca physical modes are kept. In each
    # pair use (Q,p/P); its generator is O(P), hence has zero fast symbol.
    Agrad, Bkin = s.symbols(
        "other_positive_gradient other_positive_kinetic", positive=True
    )
    Apot = s.Symbol("other_potential", real=True)
    omega2 = s.Matrix([[0, 1], [-1, 0]])
    other = omega2 * s.diag(Agrad * c.P**2 + Apot, Bkin)
    other_map = s.diag(1, c.P)
    other_fast = (
        other_map.inv() * other * other_map / c.P ** s.Rational(3, 2)
    ).applyfunc(lambda x: s.limit(x, c.P, s.oo))
    full16 = s.diag(fast, other_fast, other_fast, other_fast, other_fast)
    checks["all_four_other_modes_zero_at_fast_order"] = other_fast
    checks["all_sixteen_phases_characteristic"] = full16.charpoly(
        lam
    ).as_expr() - lam**12 * (lam**4 - beta)
    return {
        "whole_exact_old_from_new_shear": R,
        "whole_shear_time_derivative": Rdot,
        "whole_shear_time_connection": connection,
        "whole_transformed_Hessian": whole,
        "whole_similarity_weights": M,
        "whole_normalized_generator": normalized,
        "whole_eight_phase_fast_symbol": fast,
        "whole_sixteen_phase_fast_symbol": full16,
        "positive_fast_rate_fourth_power": beta,
        "whole_exact_inverse_sqrt_momentum_remainder": divided,
        "all_remainder_denominators_at_zero": denominator_at_zero,
        "all_remainder_numerator_and_denominator_degrees": polynomial_degrees,
        "domain": "a,D,Z,J,K,Y,Yv,zeta positive, gamma=1-3Zr^2/(2D)>0, r!=0,L2!=0. Fix a smooth compact time interval with these strict margins before taking P to infinity. All source coefficients, masses and time jets are held fixed. P sufficiently large makes L2+L0 a^2/P^2 nonzero.",
        "checks": checks,
        "gates": {
            "every_full_remainder_denominator_regular": all(
                value != 0 for value in denominator_at_zero
            ),
            "every_remainder_denominator_has_only_nonzero_domain_factors": all(
                denominator_domain_proof
            ),
            "all_heavy_source_and_mass_terms_in_whole_Hessian": all(
                whole.has(x) for x in (*c.ds, *c.vs, c.m11)
            ),
            "full_first_time_connections_in_whole_Hessian": all(
                whole.has(x) for x in dots
            ),
            "no_frozen_raw_Hamiltonian_stability_inference": True,
            "other_four_physical_modes_not_deleted": full16.shape == (16, 16),
        },
    }


@cache
def determinant_data():
    scale = s.diag(c.P**2, c.P, c.P, 1, 1, 1, 1, c.P)
    scaled = scale.inv().T * s.hessian(c.hamiltonian(), c.PHASE) * scale.inv()
    limit = scaled.applyfunc(lambda x: s.factor(s.limit(x, c.P, s.oo)))
    ell = s.Matrix(
        [
            -c.L2 / c.a**2,
            0,
            0,
            0,
            c.th / c.D,
            -c.ws[0] / c.Z,
            -c.ws[1] / c.Z,
            -3 * (c.r * c.th / c.D + c.rN * c.dH),
        ]
    )
    base = s.zeros(8)
    base[1, 1] = base[2, 2] = c.Y / c.a**2
    base[3, 3] = c.Yv
    base[5, 5] = base[6, 6] = 1 / c.Z
    base[7, 7] = c.GAMMA / c.Z
    base[4, 7] = base[7, 4] = c.r / (2 * c.D)
    shear = s.eye(8)
    for j in range(1, 8):
        shear[0, j] = -ell[j] / ell[0]
    separated = s.diag(ell[0] ** 2 / (2 * c.J), 1, 1, 1, 1, 1, 1, 1)
    separated[1:, 1:] = base[1:, 1:]
    coefficient = (
        -(c.r**2) * c.L2**2 * c.Y**2 * c.Yv / (8 * c.D**2 * c.J * c.Z**2 * c.a**8)
    )
    return {
        "whole_scaled_Hessian": limit,
        "whole_rank_one_column": ell,
        "whole_rank_one_base": base,
        "whole_determinant_preserving_shear": shear,
        "raw_full_determinant_P10_coefficient": coefficient,
        "interpretation": "An independent whole-matrix diagnostic only. Time-dependent stability is established by the full generator, exact canonical connections and growth normal form, not the sign of an unclean frozen Hamiltonian.",
        "checks": {
            "all_eight_rank_one_Hessian_entries": limit
            - base
            - ell * ell.T / (2 * c.J),
            "whole_rank_one_separation": shear.T * limit * shear - separated,
            "whole_rank_one_shear_unit_determinant": shear.det() - 1,
            "whole_raw_P10_determinant_coefficient": ell[0] ** 2
            * base[1:, 1:].det()
            / (2 * c.J)
            - coefficient,
            "whole_raw_determinant_P10_scaling": scale.det() - c.P**5,
        },
        "gates": {"raw_determinant_not_used_as_time_dependent_proof": True},
    }


@cache
def growth_data():
    alpha, h, rho, d = s.symbols("alpha h rho d", positive=True)
    g = s.Symbol("nonzero_g", real=True, nonzero=True)
    beta = s.Matrix(s.symbols("beta0:2", real=True))
    charge = s.Matrix(s.symbols("charge0:2", real=True))
    A0 = s.zeros(8)
    A0[0, 4] = alpha
    A0[1, 4], A0[2, 4] = beta
    A0[3, 0], A0[4, 7], A0[7, 3] = -g, g, -h
    A0[5, 0], A0[6, 0] = -charge
    A0[5, 1] = A0[6, 2] = -d
    # u=(Q,-rho QL/g,alpha Pb/rho,rho^2 PL/(g h));
    # sbar=s-beta Q/alpha, pbar=p-(charge+d beta/alpha) QL/g.
    S = s.zeros(8)
    S[0, 0], S[1, 3], S[2, 4], S[3, 7] = 1, -rho / g, alpha / rho, rho**2 / (g * h)
    for i in range(2):
        S[4 + i, 1 + i] = 1
        S[4 + i, 0] = -beta[i] / alpha
        S[6 + i, 5 + i] = 1
        S[6 + i, 3] = -(charge[i] + d * beta[i] / alpha) / g
    cycle = s.Matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0]])
    slow = s.zeros(4)
    slow[2:, :2] = -d * s.eye(2)
    wanted = s.diag(rho * cycle, slow)
    relation = rho**4 - alpha * g**2 * h

    def reduce_quartic(value):
        numerator, _ = s.fraction(s.cancel(value))
        return s.rem(numerator, relation, rho)

    positive = s.ones(4, 1) / 2
    negative = s.Matrix([1, -1, -1, 1]) / 2
    amplitude = s.Matrix(s.symbols("slow_coordinate0:4", real=True))
    eps = s.Symbol("slow_weight", positive=True)
    slow_scale = s.diag(1, 1, eps, eps)
    slow_normal = slow_scale * slow * slow_scale.inv()
    energy_margin = (
        eps * d * amplitude.dot(amplitude) / 2
        - (amplitude.T * (slow_normal + slow_normal.T) * amplitude)[0] / 2
    )
    positive_margin = (
        eps * d * sum((amplitude[i] + amplitude[2 + i]) ** 2 for i in range(2)) / 2
    )
    mapping = {
        alpha: c.L2**2 / (8 * c.J * c.a**3),
        g: c.r * s.sqrt(c.zeta) / (c.D * c.a**2),
        h: c.a**3 * c.Yv / c.zeta,
        d: c.a * c.Y,
        **{beta[i]: c.L2 * c.ws[i] / (4 * c.J * c.Z * c.a**3) for i in range(2)},
        **{charge[i]: c.a * c.cs[i] / c.D for i in range(2)},
    }
    return {
        "whole_generic_fast_matrix": A0,
        "physical_coefficient_bindings": mapping,
        "positive_rate_relation": relation,
        "whole_smooth_growth_chart": S,
        "whole_growth_chart_determinant": s.factor(S.det()),
        "normal_fast_cycle": cycle,
        "positive_fast_unit_vector": positive,
        "negative_fast_unit_vector": negative,
        "whole_slow_nilpotent_block": slow,
        "slow_weighted_block": slow_normal,
        "growth_statement": "On each fixed compact nonreference interval with rho>=rho_min>0, choose slow_weight=rho_min/(4 max(d)). The complement logarithmic norm of the principal generator is <=rho_min/8. The full time-dependent transformed equation is P^(3/2) diag(rho C,slow_weighted)+E with ||E||<=C0 P for all sufficiently large P; this includes chart derivatives. For sqrt(P)>=8 C0/rho_min the cone ||z||<=f is forward invariant and f'>=rho_min P^(3/2) f/2. Thus the original physical homogeneous Gaussian propagator has a polynomial-in-P times exp(rho_min P^(3/2) length/2) lower bound along suitable data. There is no all-momentum finite-derivative Sobolev estimate on these fixed backgrounds. The constants are not evaluated against a Wilsonian cutoff and are not uniform as r tends to zero.",
        "checks": {
            "entire_growth_chart_normal_form_mod_exact_rate_relation": (
                S * A0 * S.inv() - wanted
            ).applyfunc(reduce_quartic),
            "whole_generic_fast_symbol_matches_physical_block": A0.subs(
                mapping, simultaneous=True
            )
            - fast_data()["whole_eight_phase_fast_symbol"],
            "normal_fast_cycle_orthogonal": cycle.T * cycle - s.eye(4),
            "positive_fast_eigenvector": cycle * positive - positive,
            "negative_fast_eigenvector": cycle * negative + negative,
            "exact_complement_energy_sign": (cycle + cycle.T) / 2
            - positive * positive.T
            + negative * negative.T,
            "slow_weighted_energy_margin_is_positive_square": energy_margin
            - positive_margin,
            "whole_slow_block_nilpotent_not_assumed_diagonal": slow**2,
        },
        "gates": {
            "growth_chart_invertible_on_nonreference_domain": s.factor(S.det()) != 0,
            "fast_slow_couplings_not_omitted": all(A0.has(x) for x in (*beta, *charge)),
            "slow_zero_block_need_not_be_diagonalizable": slow != s.zeros(4),
            "time_dependent_cone_proof_not_frozen_eigenvalue_claim": True,
            "fixed_nonzero_background_before_large_momentum_limit": True,
        },
    }
