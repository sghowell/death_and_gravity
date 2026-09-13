"""First heavy absorptive coefficient and outgoing leading-rate sign, not a full width."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import model

S, M, g, beta, x = s.symbols("s heavy_mass trilinear beta Feynman_x", positive=True)
WIDTH_LOWER = 7 * 999 * model.G2 / (1000 * 32 * 22 * model.MASS2)
WIDTH_UPPER = model.G2 / (96 * model.MASS2)


@cache
def data():
    y = s.Symbol("cut_root", real=True)
    phase = s.sqrt(1 - 4 / S)
    xp, xm = (1 + phase) / 2, (1 - phase) / 2
    mixeddisc = (S - M * M - 1) ** 2 - 4 * M * M
    amix, bmix = (
        (S + 1 - M * M - s.sqrt(mixeddisc)) / (2 * S),
        (S + 1 - M * M + s.sqrt(mixeddisc)) / (2 * S),
    )
    c = s.Symbol("positive_absorptive_coefficient", positive=True)
    displacement = -s.I * c
    checks = {
        "equal_mass_cut_interval_length": s.simplify(xp - xm - phase),
        "equal_mass_upper_log_cut_endpoint": s.cancel(
            (1 - S * y * (1 - y)).subs(y, xp)
        ),
        "equal_mass_lower_log_cut_endpoint": s.cancel(
            (1 - S * y * (1 - y)).subs(y, xm)
        ),
        "mixed_Kallen_threshold_factorization": s.expand(
            mixeddisc - (S - (M + 1) ** 2) * (S - (M - 1) ** 2)
        ),
        "mixed_log_cut_interval_length": s.cancel(bmix - amix - s.sqrt(mixeddisc) / S),
        "mixed_log_cut_endpoint": s.simplify(
            (y * M * M + 1 - y - y * (1 - y) * S).subs(y, amix)
        ),
        "heavy_bubble_identical_loop_weight": g
        * g
        / (16 * s.pi**2)
        * s.Rational(1, 2)
        * s.pi
        * beta
        - g * g * beta / (32 * s.pi),
        "light_mixed_distinct_loop_weight": g * g / (16 * s.pi**2) * s.pi * beta
        - g * g * beta / (16 * s.pi),
        "full_two_body_decay_phase_space": s.integrate(
            beta / (32 * s.pi**2), (x, -1, 1)
        )
        * 2
        * s.pi
        - beta / (8 * s.pi),
        "identical_decay_rate_matches_heavy_cut": s.Rational(1, 2)
        / M
        * s.Rational(1, 2)
        * g
        * g
        * beta
        / (8 * s.pi)
        - g * g * beta / (32 * s.pi * M),
        "outgoing_formal_pole_displacement_sign": displacement + s.I * c,
        "heavy_width_over_mass": (g * g * beta / (32 * s.pi * M)) / M
        - g * g * beta / (32 * s.pi * M * M),
    }
    return {
        "heavy_first_self_energy": "Pi_H,first(s)=g^2 B0(s;1,1)/(32pi^2), up to real local mass/kinetic counterterms. The one-half loop factor comes from the two identical internal light lines.",
        "heavy_first_absorptive_coefficient": "For real s>4, Im Pi_H,first(s+i0)=g^2 sqrt(1-4/s)/(32pi)>0. The first-sheet logarithm has negative imaginary part on its negative argument, and the minus logarithm gives this positive sign.",
        "mixed_light_first_absorptive_coefficient": "For real s>(M_H+1)^2, Im Pi_phi,first(s+i0)=g^2 sqrt[(s-(M_H+1)^2)(s-(M_H-1)^2)]/(16pi s). There is no identical factor for a heavy-light pair. At this order the mixed cut lies outside the named low disk; this is not the full quantum threshold statement.",
        "leading_heavy_decay_rate": "Gamma_H,first=g^2 sqrt(1-4/M_H^2)/(32pi M_H). The phase-space formula includes1/(2M_H) and exactly one identical final1/2!. It agrees with Im Pi_H,first(M_H^2)/M_H.",
        "outgoing_pole_scope": "For the inverse s-M_H^2+Pi_H, the formal first outgoing resonance displacement on the continued second sheet has imaginary part -Im Pi_H,first(M_H^2). This sign is a perturbative decay coefficient, not an exact pole/width, a first-sheet complex pole, a controlled resonance resummation or a theorem that the full amplitude has no other states.",
        "actual_width_ratio_lower": WIDTH_LOWER,
        "actual_width_ratio_upper": WIDTH_UPPER,
        "actual_first_width_range": "10^-208<Gamma_H,first/M_H<10^-207. The actual beta exceeds999/1000 and3<pi<22/7. The narrow formal leading coefficient does not license substituting a stable delta pole into an exact quantum dispersion relation.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "heavy_light_decay_open": model.MASS2 > 4,
            "actual_heavy_beta_lower": 1 - 4 / model.MASS2 > s.Rational(999, 1000) ** 2,
            "actual_first_width_positive_lower": WIDTH_LOWER > s.Rational(1, 10**208),
            "actual_first_width_small_upper": WIDTH_UPPER < s.Rational(1, 10**207),
            "actual_width_brackets_ordered": 0 < WIDTH_LOWER < WIDTH_UPPER,
            "unstable_heavy_not_exact_stable_atom": True,
        },
    }
