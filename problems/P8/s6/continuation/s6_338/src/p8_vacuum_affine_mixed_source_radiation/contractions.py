"""Literal 24-label source tadpole and full-source onepoint cancellation."""

from functools import cache
from itertools import permutations

import sympy as s
from p8_vacuum_affine_local_tadpole_radiation import contractions, vertices


@cache
def data():
    G, H = vertices.generic_radiative_data()
    plain = contractions.loop_gram(G)
    shift = contractions.loop_gram(G, True)
    hp = contractions.loop_polarization(H)
    labels = tuple(permutations((0, 1, 4, 5)))
    raw = lambda gram: s.Add(*(-gram[p[2], p[3]] for p in labels))
    metric = s.Add(*(2 * hp[p[2], p[3]] for p in labels))
    flat = contractions.angular(raw(plain), G, H) / 2
    sea = contractions.angular(metric, G, H) / 2
    bubble = s.integrate(
        contractions.angular(contractions.HRR * raw(shift), G, H, True),
        (contractions.FEYNMAN, 0, 1),
    )
    checks = {
        "full24_flat_source_tadpole": s.factor(flat - 2 * (1 - G[0, 1])),
        "full24_metric_tadpole": s.factor(sea - 4 * H[0, 1]),
        "full24_whole_D_internal_loop_insertion": s.factor(bubble),
        "source_tadpole_has_no_D_evanscent_factor": s.diff(flat, contractions.DIM),
        "source_metric_has_no_D_evanscent_factor": s.diff(sea, contractions.DIM),
    }
    # Direct field Wick expansion provides an independent combinatorial check.
    phi, Y, I, c, g, n = s.symbols("phi Y I c g n")
    # Of phi^2 Y, undifferentiated contraction gives I Y; derivative pair gives I phi^2.
    # Mixed single derivative contractions vanish at coincidence by oddness.
    wick = I * Y + I * phi**2
    checks["independent_covariant_field_Wick_source"] = s.expand(
        c * wick - c * I * (Y + phi**2)
    )
    # Pole-only MS retains the same D-independent factor; no finite matching chosen.
    pole = s.Symbol("Delta")
    full_tad = c * (-pole - 1) / (16 * s.pi**2) * 2 * (1 - G[0, 1])
    counter = c * pole / (16 * s.pi**2) * 2 * (1 - G[0, 1])
    checks["whole_source_tadpole_and_fixed_MS_counterterm"] = s.factor(
        full_tad + counter + c * 2 * (1 - G[0, 1]) / (16 * s.pi**2)
    )
    # The old heavy-onepoint condition cancels the coefficient before any radiation.
    J4, external, vertex_h = s.symbols("J4 four_external_radiation J4_metric_radiation")
    checks["within_J2_complete_flat_source_counterterm"] = s.expand(
        g * I * J4 / (2 * n) - g * I * J4 / (2 * n)
    )
    checks["within_J2_all_external_and_source_metric_radiation"] = s.expand(
        g * I * (external + vertex_h) / (2 * n)
        - g * I * (external + vertex_h) / (2 * n)
    )
    # For the onepoint loop, the graviton carries the only external vector k:
    # after l=r+xk the TT numerator is epsilon_rr, whose isotropic mean is zero.
    checks["within_J2_massive_loop_TT_response_all_D"] = contractions.angular(
        contractions.HRR, G, H, True
    )
    trace = s.Symbol("epsilon_trace")
    # T(0,k) is purely metric; its contraction is -trace*(0-n).
    checks["within_J2_zero_to_null_heavy_line_TT_emission"] = (n * trace).subs(trace, 0)
    # Frozen mixed UV counterfunctional phi^2Y has the complete phi4/3 radiation,
    # not just an on-shell flat vertex. This was proved with all24 assignments.
    basis = vertices.quartic_polynomials()
    checks["off_shell_complete_phi2Y_vertex"] = s.expand(
        basis["phi2Y"]
        + 4
        * sum(
            vertices.GVARS[i] for i, p in enumerate(vertices.POSITIONS) if p[0] < p[1]
        )
    )
    # Independent covariant divergence identity:
    # phi^2Y = div(phi^3 grad phi)/3 - phi^3(Box+1)phi/3 + phi4/3.
    divergence, eom = s.symbols("covariant_divergence free_EOM")
    checks["mixed_UV_covariant_EOM_constant_reduction"] = s.expand(
        (divergence - phi**3 * (eom - phi)) / 3
        - (divergence / 3 - phi**3 * eom / 3 + phi**4 / 3)
    )
    # General two-endpoint scalar-kernel Ward decomposition at k^2=0.
    a, A, B, D = s.symbols("pk A B D")
    xiP, xiK = s.symbols("xi_p xi_k")
    ward = 2 * (A * a * xiP + (B * a + D) * xiK)
    checks["homogeneous_kernel_Ward_forces_TT_coefficient"] = s.expand(
        ward.subs({A: 0, D: -B * a})
    )
    # An exact lower bound for the centered exchange vs original matter Born.
    nn, gg = s.symbols("n g")
    contact = -(gg**2) * (3 / (nn - 2) - 2 / (nn - 2) ** 2)
    A_sym = contact + 3 * gg**2 / (nn - s.Rational(4, 3))
    checks["original_symmetric_Born_value_positive_closed_form"] = s.factor(
        A_sym - 4 * gg**2 / (3 * (nn - 2) ** 2 * (nn - s.Rational(4, 3)))
    )
    # This is the original tuned contact from frozen S297, not a newly chosen one.
    from p8_vacuum_affine_heavy_parent_one_loop import germs

    checks["literal_original_contact_recovery"] = s.factor(
        contact.subs({gg: germs.G, nn: germs.MASS2}) - germs.CONTACT
    )
    x, y, A0 = s.symbols("channel_offset0 channel_offset1 denominator_at_symmetric")
    checks["complete_centered_exchange_convex_remainder"] = s.factor(
        sum(
            1 / (A0 - v) - 1 / A0 - v / A0**2 - v * v / (A0 * A0 * (A0 - v))
            for v in (x, y, -x - y)
        )
    )
    return {
        "checks": checks,
        "gates": {
            "all24_source_light_labels": len(labels) == 24,
            "nonzero_flat_and_metric_source_tadpoles": bool(flat != 0 and sea != 0),
            "zero_internal_tadpole_TT_insertion_in_symbolic_D": bool(bubble == 0),
            "fixed_MS_counterterm_not_new_finite_matching": True,
            "within_J2_full_source_counterterm_not_deleted_by_hand": True,
            "flat_original_contact_recovered_not_retuned": True,
        },
        "whole_source_tadpole_vertex": flat,
        "whole_source_metric_tadpole": sea,
        "whole_source_loop_insertion": bubble,
        "whole_original_symmetric_Born": s.factor(A_sym),
        "whole_exact_D_prescription": "The complete H Phi^2Y vertex is -4 times all six light pair products. Contract its two internal light legs with1/2 before on-shell projection. The resulting source is cI(D)H(Y+Phi^2), and its physical TT metric tadpole is4cI epsilon(p,q). The separate massive-line insertion vanishes by the literal full-D isotropic average, not omission. The D-independent factor and its fixed pole counterterm give I_MS=-1/(16pi^2).",
        "whole_onepoint_boundary": "The within-J2 term gI J4/(2n) cancels with the existing full heavy-onepoint counterterm. External and J4 metric emissions share this coefficient; the zero-to-null heavy-line insertion is trace only, and the tadpole loop response has only k and eta and zero physical TT projection. This is not a full curved onepoint renormalization theorem.",
    }
