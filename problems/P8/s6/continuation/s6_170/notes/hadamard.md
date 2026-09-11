# Identifying the same actual in/out states as Hadamard

The external input is [Gerard and Stoskopf](https://arxiv.org/pdf/2108.11630),
v3 of11 January2022: Definition2.3, hypotheses in4.5,
Proposition6.8, Theorem6.9 and Proposition3.8.
It supplies a local pure Hadamard reference with order-zero energy
projectors whose evolution defect is smoothing. Its Dirac mass may
be a smooth Hermitian endomorphism. The theorem is applied, not
proved or formalized here.

## Hypotheses and conventions

The physical geometry is not asymptotically static. The S6.166
flat in/out theorem is therefore not invoked unchanged.

Here is an explicit witness for the local theorem hypotheses.
Define E(s)=exp(-1/s) for s>0 and zero otherwise,
theta(s)=E(s)/(E(s)+E(1-s)), and
zeta(t)=theta(t+2)theta(2-t). Then zeta is smooth, equals one
on[-1,1], vanishes outside[-2,2], and lies in[0,1].
The comparison metric with scale
a_tilde=1+zeta(t)[a(t)-1] has1<=a_tilde<=25,
bounded metric/inverse derivatives relative to dt^2+dx^2,
and the trivial spin structure. Every constant-time slice is Cauchy.
It equals the actual metric on U=(-1,1)xR^3. Since cosmic time is
strictly monotone on causal curves, a curve between points of U
cannot leave this time slab and return; causal compatibility holds.
Set u_tilde=0. The original SAT8 mass is bounded with bounded
derivatives of every finite order for the fixed tau>0; its spatial
derivatives vanish. The surface t=0 has Euclidean bounded geometry.
These directly verify(H1-H3),(M). The comparison metric is only
a witness for theorem applicability, not a changed physical
matter frame or a separately prepared physical state.

Use the checked S6.166 Clifford dictionary in an orthonormal
tetrad: Gamma0=-i gamma_plus^0, Gammaj=i gamma_plus^j,
beta_source=gamma_plus^0 and m_source=-M. Covariant derivatives
give the same physical curved Dirac operator. The a^(3/2)
Cauchy normalization makes H_source=-H_physical.
Exchange the source spectral labels, retaining physical negative
energy as the occupied subspace.

## Our all-order comparison

On a compact time interval, let P_N be the exact projector obtained
by conjugating the negative diagonal projector with the mass
rotation and N further exact rotations. At high p all coefficients
have classical inverse-momentum expansions with bounded derivatives.
Specifically g0 is order-1, e0 is order1, each successive connection
g_j is order-j-1, and the rotations are order0.
The evolution defect of P_N is a conjugate of the commutator with
g_N sigma_axis, hence order-N-1. P_N is an EXACT orthogonal
projector, not a truncated polynomial masquerading as a state.

Let pi(t) denote the fixed reference projector from the external
construction. Both principal projectors agree. If
D=P_N-pi first appears at order-k, k>=1, their projection
identities give at leading order

    P0 D + D P0 - D =0.

Its two diagonal energy blocks vanish. Subtracting the evolution
equations at order1-k gives [h1,D]=0, provided k<=N+1.
The off-diagonal blocks vanish because the principal energies
+p/a and-p/a have a nonzero gap. Therefore D has one lower order.
Induction gives P_N-pi in Psi^(-N-2).
The argument works for the rank-two energy blocks of the full
four-component spinor; it does not assume nondegenerate helicity
eigenvectors or a global smooth helicity frame.

The exact in/out Cauchy projectors c_in/out differ from P_N
by the corresponding half-line transition estimate. At each t,

    ||c_in/out(t,p)-P_N(t,p)||<=integral_halfline|g_N|
                                  =O(p^(-N/4)).

For any desired Sobolev smoothing order, choose N large enough
and then its finite high-momentum threshold. A bounded
Fourier multiplier of sufficiently rapid decay maps H^-s to H^s;
bounded low-momentum multipliers do so as well. This argument
does not require momentum derivatives of the exact state.
Together with the pseudodifferential comparison above, it gives

    c_in/out(0)-pi(0) in W^-infinity.

The same conclusion holds for the complementary covariance.
Local exact evolution preserves every Sobolev space; each time
derivative costs a finite number of momentum powers.
Thus the spacetime covariance difference from the Hadamard
reference has a smooth kernel. The same actual in/out states
are Hadamard near t=0 and globally by propagation of singularities.

N20 supplies numerical bounds only. Hadamard follows from the
arbitrary-order argument, not finite energy, not a finite
adiabatic order, and not the difference of two unidentified
states. No quantitative absolute local stress follows here.
