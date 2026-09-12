# S6.224: full flat scalar quotient reference inverse

Keep the source-pinned S223 shear factor, its subthreshold pole, the S86 trace factor, the S200 full nonzero-transfer covariant projectors and the original finite Hessian. The target is an isolated flat reference in the two-dimensional spatial scalar quotient, not the full curved quantum-force system.

For each Fourier direction let q=|P|^2 and Q=2zeta I-2c Pi, Pi=Phat Phat^T. Put D=partial_t and

S=(D^2+2q/3)zeta-D^2 c/3,
W=-q zeta-D^2 c,
B=[[D^2+2q/3,-D^2/3],[-q,-D^2]].

The full spatial covariant projectors give p^2 Pi0=12 S_D S_G and p^2 Pi2=(8/3)W_D W_G, p=lambda^2+q. Including the original finite coefficients, the normalized quotient reference is

A_q=B^T diag(Ftrace(D^2+q),(8/3)F2(D^2+q)) B.

The complete three-source reference in(n,zeta,b) has a gauge kernel. Both source and detector time legs are reduced explicitly; it is not called a full three-source inverse.

The causal B inverse is built from s_q(t)=sin(sqrt(q)t)/sqrt(q), with its continuous s_0=t extension:

R_q=[[s_q,-s_q/3],[s_q-t,-s_q/3-2t/3]].
R_q(0)=0, ||R'_q||<5/2, ||R_q(t)||<5t/2.

The shifted positive inverse spectral measures give primitives Jtrace,q and J2,q with bounds1/2 and60 uniformly inq>=0. J2,q includes its undamped pole. For Jq=diag(Jtrace,q,(3/8)J2,q), ||Jq||<=45/2. The ordinary convolution

E_q=R'_q*Jq*R_q^T

has both causal inverse identities for A_q and obeys

||E_q(t)||<=375t^3/16,
||E_q*||<=375T^4/64<6T^4.

For every realr the latter is a uniform C_tH^r-to-C_tH^r and L2_tH^r bound for the two amplitude coordinates and their dual sources, with no spatial derivative loss or transfer cutoff. A first-time-derivative bound is also proved. The graph domain includes the complete initial boundary, not just the equation at strictly positive times.

These are normalized reference units. The amplitude-coordinate physical Hessian inverse restores64pi^2; the S222 force-normalized flat reference also restoreskappa. Its bound is375pi^2 kappa T^4, not kappa smallness. The amplitude-to-Frobenius metric embedding has squared norm below14, which must be retained if translating source/output norms.

A full curved/state/contact/tree/matter normal form and its controlled remainder are not supplied here. Neither are the S222 coupled inverse, a finite-coupling/nonlinear parent remainder, stability, cutoff or original V/G/B/P8 closure. All frozen inputs and previously scoped objectives are unchanged.
