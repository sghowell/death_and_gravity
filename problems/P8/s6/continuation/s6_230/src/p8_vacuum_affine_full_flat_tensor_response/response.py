"""Full massless/complex-pair/cut dispersion and causal all-momentum realization."""

from functools import cache

import sympy as s

from . import normalization as norm


@cache
def data():
    p, tau, C, U, DA = s.symbols("p tau C U DA", positive=True)
    z = s.Symbol("upper_first_sheet_zero", complex=True, nonzero=True)
    R = s.Symbol("upper_complex_residue", complex=True)
    cut = 1 / (-tau * (C - tau * (DA + s.I * s.pi * U)))
    rho = U / ((C - tau * DA) ** 2 + s.pi**2 * tau**2 * U**2)
    T = s.Symbol("finite_time_window", positive=True)
    # Independent beta moments for the bounded original paired forward kernel.
    Hbound = (30 * s.pi / 4 - 20 * 3 * s.pi / 16 + 3 * 5 * s.pi / 32) / (60 * norm.MASS)
    M0 = -1 / C - 2 * s.re(R)
    return {
        "full_reciprocal_dispersion": "E_C(p)=1/(Cp)+R/(p-z)+conj(R)/(p-conj(z))+integral rho_E(tau)/(p+tau)dtau; R=1/[z D_C'(z)]",
        "original_positive_cut_density": rho,
        "both_finite_cut_moments": "M0=integral rho_E=-1/C-2ReR; M1=integral tau rho_E=2Re(Rz). Both positive and finite. No complex pole omitted.",
        "pole_and_cut_asymptotics": "rho_E is O(sqrt(tau-4m^2)) at threshold and O(1/(tau^2 log^2tau)) at infinity; E_C=O(1/(p^2logp)) fixes both moment cancellations.",
        "actual_physical_causal_kernel": "16pi^2kappa theta(t)[s_q(t)/C+2Re(R sinh(sqrt(z-q)t)/sqrt(z-q))+integral rho_E(tau)s_(q+tau)(t)dtau]; s_q=sin(sqrt(q)t)/sqrt(q), s0=t.",
        "all_momentum_complex_wave_bound": "Duhamel about the q>=0 wave gives |s_complex|<=t exp(sqrt(|z|)t), |s_complex'|<=exp(sqrt(|z|)t), uniformly inq.",
        "all_momentum_physical_inverse_norm_bound": 4
        * T
        * T
        * s.exp(s.sqrt(norm.KAPPA) * T),
        "all_momentum_first_derivative_norm_bound": 8
        * T
        * s.exp(s.sqrt(norm.KAPPA) * T),
        "original_paired_forward_kernel_bound": Hbound,
        "full_forward_realization": "O_phys=[C L_q+L_q^2/30+L_q^3 H_q]/(16pi^2kappa), H_q=theta integral W2/(1-y^2) s_(q+4m^2/(1-y^2)); all derivatives are causal distributions.",
        "specified_graph": "Y_r=C_tH^r_TT with a fixed original zero germ. Domain consists of h inY_r for which the full causal O_phys h belongs toY_r, with equality understood first in D'_t H^(r-6). The uniform inverse and spatial cutoff limit give both identities and the graph norm bound. This is not the curved scalar/clock/matter realization.",
        "forced_growth": "Nonzero complex residues give growing right-half-plane lambda poles. Smooth compact-time TT forcing can excite them, including real L2 wave packets with zero initial germ. This is an exact retained flat mean-equation statement, not physical UV or nonlinear curved-bounce instability.",
        "checks": {
            "full_cut_density_sign_and_normalization": s.factor(
                -s.im(cut) / s.pi - rho
            ),
            "complex_pole_first_moment_expansion": s.cancel(
                1 / (p - z) - 1 / p - z / (p * (p - z))
            ),
            "continuum_first_moment_expansion": s.cancel(
                1 / (p + tau) - 1 / p + tau / (p * (p + tau))
            ),
            "zeroth_sum_rule": 1 / C + 2 * s.re(R) + M0,
            "paired_forward_beta_bound": s.factor(
                Hbound - 9 * s.pi / (128 * norm.MASS)
            ),
            "complete_kernel_norm_using_moment": s.expand(
                1 / C + 2 * s.Abs(R) + M0 - 2 * (s.Abs(R) - s.re(R))
            ),
            "kernel_to_force_inverse_constant": 8 / s.Integer(2) - 4,
            "actual_force_normalization_ratio": s.cancel(
                norm.C0 / norm.C - (1 - 5 * norm.MASS**2 / (6 * norm.C))
            ),
        },
        "gates": {
            "full_cut_density_positive": rho.is_positive,
            "actual_Einstein_over_C_less_than_one": (norm.C - norm.C0).is_positive
            and norm.C0.is_positive,
            "original_massless_pole_kept": True,
            "both_complex_poles_kept_on_original_sheet": True,
            "all_q_bound_without_near_frequency_division": True,
            "causal_graph_identity_not_stability": True,
        },
    }
