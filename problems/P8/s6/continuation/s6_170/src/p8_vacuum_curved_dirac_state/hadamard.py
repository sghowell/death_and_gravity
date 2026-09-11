"""All-order projector comparison with a local external Hadamard reference."""

from functools import cache

import sympy as s


@cache
def data():
    a, b, c, d, q = s.symbols("a b c d q")
    P = s.diag(1, 0)
    D = s.Matrix([[a, b], [c, d]])
    H = s.diag(q, -q)
    proj = P * D + D * P - D
    comm = H * D - D * H
    checks = {
        "leading_projection_diagonal_11": proj[0, 0] - a,
        "leading_projection_diagonal_22": proj[1, 1] + d,
        "leading_projection_offdiagonal_12": proj[0, 1],
        "leading_projection_offdiagonal_21": proj[1, 0],
        "Sylvester_positive_gap_12": comm[0, 1] - 2 * q * b,
        "Sylvester_negative_gap_21": comm[1, 0] + 2 * q * c,
        "zero_diagonal_and_offdiagonal_imply_zero": s.trace(
            D.subs({a: 0, b: 0, c: 0, d: 0})
        ),
    }
    return {
        "external_primary_reference": {
            "url": "https://arxiv.org/pdf/2108.11630",
            "version": "v3, 11 January 2022",
            "locations": "Definition2.3; hypotheses(H1-H3),(M) in4.5; Proposition6.8 and its construction; Theorem6.9; Proposition3.8",
            "role": "A local pure Hadamard reference with order-zero projectors, the correct principal energy branch and a smoothing evolution defect. This external theorem is applied, not re-proved or formalized.",
        },
        "hypothesis_witness": "Use U=(-1,1)xR^3 and a smooth comparison scale a_tilde=1+zeta(t)(a(t)-1), 1<=a_tilde<=25, zeta=theta(t+2)theta(2-t), theta(s)=E(s)/(E(s)+E(1-s)), E(s)=0 for s<=0 and exp(-1/s) otherwise. The comparison is exactly physical on U, globally hyperbolic and bounded relative to dt^2+dx^2. Set u_tilde=0. The original smooth scalar mass has bounded derivatives for each fixed tau. Cosmic-time monotonicity gives causal compatibility. This only witnesses theorem hypotheses, not a change of physical metric or a different physical state.",
        "conventions": "Gamma0=-i gamma_plus^0, Gammaj=i gamma_plus^j, beta_source=gamma_plus^0 and m_source=-M. The source curved slash-plus-mass equals the physical i gamma^mu nabla_mu-M. Rescaling a^(3/2) gives H_source=-H_physical, so the source spectral labels exchange and the physical occupied negative-energy projector remains fixed.",
        "finite_order_symbols": "On each compact time interval P_N is an exact selfadjoint order-zero classical projector with the negative principal branch. Its exact evolution defect is the conjugate of [g_N sigma_axis,P_minus], order -N-1. Smooth high-momentum cutoffs handle the zero-momentum helicity chart; bounded low-momentum multipliers are Sobolev smoothing.",
        "formal_uniqueness_proof": "Compare P_N with the external reference pi(t). If their difference starts in order -k, k>=1, projection identities force its leading diagonal energy blocks to vanish. The evolution equations in order1-k force its off-diagonal blocks to commute with the principal Hamiltonian; the nonzero gap2p/a forces them to vanish. Induction through k=N+1 gives P_N-pi in Psi^(-N-2). This works for both rank-two energy blocks of the four-component spinor.",
        "same_in_out_state_identified": "The exact in/out projectors differ from P_N by at most the half-line integral|g_N|=O(p^(-N/4)). For each desired Sobolev order choose N sufficiently large, then its finite high-p threshold. Combining with the fixed reference pi makes each difference W^-infinity. Local evolution is bounded on every Sobolev space and time derivatives cost finite powers of momentum. Thus the spacetime two-point differences are smooth, and these SAME actual curved in/out states are Hadamard globally by propagation.",
        "excluded_shortcuts": "Neither finite N20 nor in-minus-out finiteness alone proves Hadamard. The scale factor is not asymptotically static, so the S6.166 flat in/out theorem is not imported unchanged. No infinite-order constants are inferred from one fixed m*tau gap. No momentum derivative bound of the exact state is needed for the Sobolev multiplier argument.",
        "checks": checks,
    }
