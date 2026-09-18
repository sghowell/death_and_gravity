"""Original scalar vertex, prescriptions and nonzero phase audit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree
from p8_vacuum_affine_minimal_gravity_radiation import vertices


@cache
def data():
    eta = vertices.ETA
    etas = -eta
    p = s.Matrix(s.symbols("p0:4"))
    b = s.Matrix(s.symbols("b0:4"))
    q = s.Matrix(s.symbols("q0:4"))
    hv = s.symbols("h0:10")
    h = s.zeros(4)
    index = 0
    for i in range(4):
        for j in range(i, 4):
            h[i, j] = h[j, i] = hv[index]
            index += 1
    ps, bs, qs, hs = -p, -b, -q, -h
    mu, L, delta = s.symbols("mu L delta", positive=True)
    dot = lambda a, b: (a.T * eta * b)[0]
    dots = lambda a, b: (a.T * etas * b)[0]
    pair = lambda A, B: sum(A[i, j] * B[i, j] for i in range(4) for j in range(4))
    check = {}

    def put(name, value):
        if isinstance(value, s.MatrixBase):
            value = value.applyfunc(s.factor)
            assert value == s.zeros(*value.shape), (name, value)
        else:
            value = s.factor(value)
            assert value == 0, (name, value)
        check[name] = value

    stress_ss = ps * bs.T + bs * ps.T - etas * (dots(ps, bs) - mu)
    put(
        "original_general_stress_matches_opposite_signature",
        stress_ss - vertices.stress(p, b, mu),
    )
    put(
        "original_scalar_vertex_i_stress_matches_Sen_minus_i_vertex",
        -s.I * pair(hs, stress_ss) - s.I * pair(h, vertices.stress(p, b, mu)),
    )
    put(
        "massive_Feynman_prescription_preserved",
        -s.I / (-L + mu - s.I * delta) - s.I / (L - mu + s.I * delta),
    )
    put(
        "massless_Feynman_prescription_preserved",
        -s.I / (-L - s.I * delta) - s.I / (L + s.I * delta),
    )
    P = lambda a, b: (a.T * h * b)[0]
    Ps = lambda a, b: (a.T * hs * b)[0]
    D, Db = dot(p, q), dot(b, q)
    put("one_leg_leading_soft_sign", s.cancel(Ps(ps, ps) / dots(ps, qs) - P(p, p) / D))
    B = P(p, p) * Db / D - P(p, b)
    Bs = Ps(ps, ps) * dots(bs, qs) / dots(ps, qs) - Ps(ps, bs)
    put("pair_derivative_contraction_changes_sign", s.cancel(Bs + B))
    put("pair_invariant_changes_sign", dots(ps, bs) + dot(p, b))
    put("Doppler_squared_log_argument_preserved", dots(ps, q) ** 2 - D**2)
    # A nonforward exact Born configuration compares directly to the frozen current.
    E = s.Rational(5, 4)
    r = s.Rational(3, 4)
    mom = [
        s.Matrix([-E, 0, 0, -r]),
        s.Matrix([-E, 0, 0, r]),
        s.Matrix([E, 3 * r / 5, 0, 4 * r / 5]),
        s.Matrix([E, -3 * r / 5, 0, -4 * r / 5]),
    ]
    Q = s.Matrix([1, 1, 0, 0])
    eps = s.diag(0, 0, 1, -1) / s.sqrt(2)
    put("original_Born_momentum_conservation", sum(mom, s.zeros(4, 1)))
    for j, k in enumerate(mom):
        put("original_mass_shell_" + str(j), dot(k, k) - 1)
    direct = sum((k.T * eps * k)[0] / dot(k, Q) for k in mom)
    put(
        "original_complete_leading_current_correspondence",
        pair(eps, tree.soft_current(mom, Q)) - direct,
    )
    sen = sum(((-k).T * (-eps) * (-k))[0] / dots(-k, -Q) for k in mom)
    put("original_Born_current_Sen_convention", sen - direct)
    # Forward hard-angle LIMIT of the coefficient, not a value of the hard pole.
    forward = [
        s.Matrix([-E, 0, 0, -r]),
        s.Matrix([-E, 0, 0, r]),
        s.Matrix([E, 0, 0, r]),
        s.Matrix([E, 0, 0, -r]),
    ]
    den = [dot(k, Q) for k in forward]
    ep = lambda a, b: (a.T * eps * b)[0]
    soft = sum(ep(k, k) / d for k, d in zip(forward, den))
    put("forward_leading_current_zero", soft)
    cl = s.S.Zero
    real_by_invariant = {}
    for i, a in enumerate(forward):
        for j, b in enumerate(forward):
            if i == j:
                continue
            z = dot(a, b)
            B = ep(a, a) * den[j] / den[i] - ep(a, b)
            real_by_invariant[abs(z)] = real_by_invariant.get(abs(z), s.S.Zero) + B
            if a[0] * b[0] > 0:
                cl += z * (2 * z * z - 3) / (z * z - 1) ** s.Rational(3, 2) * B
    for z, value in real_by_invariant.items():
        put("forward_quantum_group_zero_" + str(z), value)
    put(
        "forward_imaginary_coefficient_retained",
        cl / (8 * s.pi) + s.Rational(9843, 9000) / (s.sqrt(2) * s.pi),
    )
    assert cl != 0
    return {
        "checks": check,
        "gates": {
            "original_scalar_vertex_trace_mass_and_i_preserved": True,
            "signature_change_not_loop_complex_conjugation": True,
            "pair_invariant_and_derivative_signs_both_mapped": True,
            "imaginary_log_phase_not_discarded": cl != 0,
            "forward_coefficient_limit_not_hard_amplitude_point_value": True,
        },
        "whole_convention_map": "All-outgoing P8 k_i and omega*q map to Sen p_i=-k_i,ksoft=-omega*q,etaSS=-eta,A_SS=-A,sqrt(8*pi*G)=1/sqrt(kappa).",
        "whole_forward_limit_phase": s.factor(s.I * cl / (8 * s.pi)),
        "whole_primary_theorem_boundary": "Sahoo-Sen1808.03288v3 e1.17int/eqgrsoft and Krishna-Sahoo2308.16807v2 FIRST line S1_loop_gr and section4.3 eq:oneloop-soft_gr-thm. Fixed other-leg kinematics; arbitrary mass/spin with tree gauge-fixed EFT for massless legs. Neither uniform O(1) remainder nor conjectural omega-log/two-loop terms are inferred.",
    }
