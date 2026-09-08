# S6.43 working scope: prepared leading nonlinear vector response

Status: calculation in progress; no certificate yet. Original P8 OPEN.
This is not full nonlinear heavy elimination.
The literal S6.42 action is unchanged. This checkpoint concerns its
leading, second-order forced vector response about the actual bounce,
not a nonlinear solution or an alternative parent.

1. Derive the source term in the actual constrained canonical equation
   from the nonzero second-order S6.42 source. Keep the full temporal
   vector equation before eliminating it. For a scalar source S2,
   test J=(g_L²*S2)'/g_L, rather than replacing the force by S2 itself.
2. Bound the actual time-dependent frequency and its first two
   derivatives at fixed comoving momentum, including normalization
   derivatives. Use a declared compact interval inside [-1/2,1/2],
   length <=1, 0<k_com²<=1 and 0<zeta<=1/20000.
3. Impose zero heavy initial data and prepared forcing J=J'=0 at the
   left endpoint. For a C² forcing with all three jet norms bounded
   by A, obtain an explicit retarded error against the local response
   zeta*J. This is an absolute norm estimate, not a pointwise relative
   error at zeros, a Fourier band or an arbitrary-state claim.
4. Supply a sufficient source preparation S2=S2'=S2''=0 initially
   and a C³ source-jet bound E. Restore the physical spatial-vector
   and temporal readouts, not just the canonical variable. State what
   derivative control is additional to S6.42's spatial source bound.

The exact coefficient ODE fixes the previously symbolic quadratic term:

    S2=-(2/h)*n*delta_K-[6H/h+27H/(8h²)]*n².

On [-1/2,1/2], C³ light-jet envelopes epsilon_n, epsilon_K give the
sufficient source envelope E=514*epsilon_n*epsilon_K+(132027/8)*epsilon_n².
For Fourier modes these are summable modewise time-jet envelopes, so
the products use convolution. The momentum restriction is on each
output source mode, not an incoming light mode. See notes/light-jets.md.

With the specified preparation and domain, the canonical error obeys
||v-zeta*J||<=zeta*A/80. The physical spatial-vector error against
zeta*sqrt(q)*(S2'+2rho*S2) is <=6*zeta*E; the temporal error against
S2 is <=500*zeta*E. These are normalized amplitudes, with physical
factor 1/tau. Sufficient preparation is n=n'=n''=0 at the left endpoint.

The source class is a declared perturbative forcing class; no claim
that every such source comes from an on-shell light solution is made.
If the actual light evolution satisfies the stated preparation and jet
bounds, the leading vector response obeys the derived estimate.
The homogeneous spatial mode and the transverse block have no scalar
source at this order and are separate zero-data controls.

No full nonlinear remainder, nonlinear secondary-constraint rank,
induced quantum action, loop bound, stationary gap, interacting cutoff,
vacuum/finite-gravity V/G/B gate or original P8 closure is asserted.
