# Complete covariance bridge and reflected ordered response

## Canonical current tangent

At a finite reference-mode projector write the real symmetrized covariance as C and the canonical symplectic matrix as J. With A=J M0, the full source tangent is

C1'=A C1+C1 A^t+J M_G C+C(J M_G)^t, C1(initial)=0.

The full current tangent is
-(1/2)tr(M_D C1+M_DG C).
It includes both reverse source blocks and the entire instantaneous second-Hamiltonian contact.

An independent Wick calculation checks the sign and factor. For propagation U(t,s), the unsymmetrized two-point matrix is W=U(C+iJ/2); the reversed ordering uses U(C-iJ/2), not an independently conjugated source. Expanding the connected quadratic commutator gives exactly the covariance source above after multiplication by i and the current factor. This matrix identity does not depend on a chosen phase of the mode or on stationarity.

## Full Fourier ordering, including odd channels

For a real Cartesian ADM component basis A define

V_A(t;p,q)=F_p(t)^t M_A(t)F_q(t), q=P-p.

The connected two-quadratic Wick contraction has factor1/2 from the two identical internal contractions. Before integrating momentum and multiplying by a(t)³a(s)³, the causal current memory is

(i/2)theta(t-s) [
V_A(t;p,q)conjugate(V_B(s;p,q))
-conjugate(V_A(t;-p,-q))V_B(s;-p,-q) ].

This is the FULL ordered form. Complex external coefficients are inserted bilinearly only after component kernels are formed. Conjugating a complex supplied source coefficient is wrong.

The physical annihilation feature has components(S',ip_i S/a,mS). Simultaneous internal momentum reflection acts by R=diag(1,-I3,1), WITHOUT conjugating the time mode. Lapse and spatial-metric features obey R M R=+M; shift features obey R M R=-M. Thus even-even and odd-odd components reduce to minus the imaginary part of the product, whereas mixed even-odd components reduce to plus i times its real part. The latter Fourier kernels need not be real; they obey the correct opposite-momentum reality identity.

The independent seven-mode lattice calculation verifies all ordered N,Q,beta channels, both reflected transfers and nontrivially squeezed states against the exact canonical covariance tangent. It also checks the full single-mode contact. A prepared lapse evolution, with nonzero source at readout, separately verifies the finite-difference current including its nonzero final contact. The deliberately same-imaginary-part odd-channel rule fails strongly.

## Research context, not a borrowed curved-state estimate

[Anderson, Molina-Paris and Mottola](https://arxiv.org/html/gr-qc/0209075), Section III, equations21-25, derive the causal stress response together with its local metric contact and explain why their sum uses the same covariant counterterms as the mean equation. That supports the structure being checked here. Their flat-space application is not a theorem about this curved massive state, these quantitative constants, or the full P8 parent. The feature, ordering and covariance identities above are derived independently.
