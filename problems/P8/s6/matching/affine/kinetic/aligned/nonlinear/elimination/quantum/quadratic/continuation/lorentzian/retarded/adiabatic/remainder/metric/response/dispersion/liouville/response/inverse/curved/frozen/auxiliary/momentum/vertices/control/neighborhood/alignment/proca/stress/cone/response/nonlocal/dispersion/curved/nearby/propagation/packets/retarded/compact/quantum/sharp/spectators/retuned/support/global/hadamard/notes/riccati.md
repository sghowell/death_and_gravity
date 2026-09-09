# All-orders negative-frequency graph

Use the fresh canonical gamma variables and write their Hamiltonian
matrix as mathcal H(k)=-Omega M_Y(k). Set rho=1/k. The EXACT scaled
Hamiltonian is h(rho)=mathcal H/k=h0+rho*h1+rho^2*h2.
All h_j are real symmetric smooth matrices on a neighborhood of the
compact initial slab. The leading diagonal blocks are K^-1 and the
positive spatial block stated in notes/canonical.md.

A graph p=R x is dynamically invariant precisely when

    rho R' + R h_pp R + R h_pq + h_qp R + h_qq = 0.

Seek R~sum rho^j R_j, with R0=-i B0. At order n>=1 the only terms
containing the new R_n are

    L(R_n)=R0 K^-1 R_n + R_n K^-1 R0.

Every other term S_n depends on already determined R_0,...,R_(n-1),
their first time derivatives, and the complete h0,h1,h2. Thus
L(R_n)=-S_n. For an arbitrary symmetric forcing F set Ftilde=V^T F V.
The UNIQUE inverse is

    X=V^-T Xtilde V^-1,
    Xtilde_ij=i*Ftilde_ij/(omega_i+omega_j).

This identity is checked with an arbitrary symmetric matrix. All
denominators are positive frequency SUMS, separated from zero on
the entire initial slab. Symmetry of the Riccati expression and
of this inverse gives complex-symmetric R_n by induction.
No small clock--matter difference denominator occurs in this graph
recursion. Subsequent diagonal phase propagation is a separate step.

This is an all-integer induction: expanding a finite polynomial
through a chosen order produces S_n; the displayed invertible
linear map produces the next coefficient, with no unsolved condition.
The implementation computes this same rule. Native residual tests
through order four exercise two coupled fields, time-dependent
off-diagonal transport, and both lower Hamiltonian blocks.
The actual global center jets independently check its first
nontrivial derivative contribution. These finite tests support the
algebra; the induction above, not a finite order cutoff, proves
existence of the formal sequence.

To obtain a genuine symbol, choose a real smooth cutoff eta which
vanishes below one and equals one above two. Choose successively
large radii K_j and form

    R_B(u,k)=R0(u)+sum_(j>=1) eta(k/K_j)*k^-j*R_j(u).

For each j choose K_j so every time/frequency seminorm of order at
most j/2, measured as a symbol of order -N for N<=j/2, of that summand
is at most 2^-j. This is possible because j>N and all coefficient
derivatives are bounded on a slightly larger compact slab. Fixing any
finite seminorm leaves only finitely many uncontrolled early summands;
the remaining geometric tail converges. Applying the same argument
after subtracting the first N formal terms shows the required
asymptotic expansion and its differentiated remainders.
A diagonal enumeration of the seminorms gives this construction to
every order. Real cutoffs preserve complex symmetry.

Finite products and time derivatives respect these asymptotics.
The Riccati defect of R_B is therefore in S^-infinity, uniformly with
all derivatives on the initial slab. Since -Im R0=B0 is positive with
a compact positive lower eigenvalue, -Im R_B remains positive for
all k above some finite K>=64. This K is an existence threshold,
not an interacting EFT validity bound. Write R_B=A-i B with A,B
real symmetric and B>0.

Let F solve x'=(M_xx+M_xp R_B)x, F(0)=I, and set
Y_app=(I,R_B)^T F. The leading configuration generator is
-i*k*K^-1 B0, diagonalized by V with two positive frequencies.
A smooth leading symmetrizer and Gronwall give polynomial (in fact
bounded at zeroth derivative) real-k transfer on the compact slab.
All time and frequency derivatives have polynomial bounds by
differentiating the equation. The residual of Y_app under the full
canonical equation is the smoothing Riccati defect times F.
Duhamel with the complete polynomially bounded transfer shows that
the exact solution with the same initial graph differs by S^-infinity,
with all derivatives. No lower Hamiltonian term is discarded.

For context, the scalar pseudodifferential construction in
[Gérard--Wrochna, Appendix A.3](https://arxiv.org/pdf/1209.2604)
uses an asymptotic Riccati strategy. Its Klein--Gordon theorem is not
applied to this system: the coupled matrix inverse, actual canonical
bridge, positivity, and continuation are proved here.
