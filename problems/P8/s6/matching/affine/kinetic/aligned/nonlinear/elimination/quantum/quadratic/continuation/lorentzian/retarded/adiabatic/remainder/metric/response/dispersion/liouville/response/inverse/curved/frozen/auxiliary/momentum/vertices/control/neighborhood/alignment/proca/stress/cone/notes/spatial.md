# Independent spatial reduction and Euler-first principal limit

The four raw scalar phase channels are curvature v, its current
canonical momentum p, matter chi and matter momentum pi_chi.
Generate all ten unordered pairs and their reversed orders from
the full physical metric, scalar curvature and the three spatial
momentum constraints, using two independent labelled Fourier legs.
The matter current is retained. No isolated Proca or minisuperspace
principal matrix replaces this calculation.

Use a formal first time jet with p0=eta*p0' and b=eta*b_u.
Every even coefficient has zero first time derivative at the
central slice. Differentiate the full lapse-dependent Hamiltonian
at fixed p0' and P before their background equations are imposed.
The lapse Hessian is

    h_NN=(16874985N^4-14244998N^2-8625003)/(2000000N^(7/2)).

The full linear lapse force contributes -force^2/(2h_NN) to
the reduced quadratic Hamiltonian. For two labelled legs this
coefficient is -force_1*force_2/h_NN. Project its exact first
time jet before rational division; higher time powers cannot
contribute. This is an algebraic simplification, not a dropped
physical term or a change to the polynomial backend.

Keep the actual moving metric term -2H_hat*pi:g, the canonical
trace boundary -(H_tilde'+3H_hat H_tilde)(6v+3v^2) with
H_tilde=-p0/2, and the matter background shift. At the central
slice the trace boundary contributes 3p0' to the vv Hessian.
The comparison with the previous central spatial matrix includes
this term explicitly.

Use q=k_com^2/R^2 and the regular gamma transformation

    v=P_b/(2q), p=-2q*b,

retaining its moving generator -H_hat*b*P_b. All ten pair
reversal identities use simultaneous Fourier momentum reversal.
Write the resulting canonical Hamiltonian as
(P^T A P+2P^T B Q+Q^T C Q)/2. Its velocity coefficients are

    alpha=A^-1, beta=-A^-1 B, gamma=C-B^T A^-1 B.

The actual central time jets are formed BEFORE taking q to
infinity. Since H_hat=0 at the slice, the principal gradient is

    G=lim_(q->infinity) (gamma+beta')/q,

not gamma/q. The leading beta itself vanishes at the slice,
but its derivative does not. Its only q-leading entry is
beta'_bb/q=4N^(3/2)Fp'/3.

Independently, the leading lapse force in gamma variables is
c*P_b+D*pi_chi+e*q*b, with c=1/(2N^(3/2)), D=-P*c and
e=-2Fp/3. The bare scalar p^2 cancels against the solved shear
term. The resulting high-q kinetic inverse acting on (c,D)
gives beta_bb/q=4N^(3/2)Fp/3. Differentiating at the central
data reproduces the full generated result.

The other principal potential entries are N^(3/2)P and
N^(3/2). Actual fixed q matrices, first time jets and limiting
matrices are serialized. Their rational q limits have nonzero
leading denominators uniformly on the stated compact N domain.
Consequently this is a regular local characteristic calculation,
not a finite-q phase-velocity interpretation.
