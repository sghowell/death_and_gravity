# Complete MS vertex as a paired difference plus its anchor

Let S(l)=(m+i slash l)^-1 and
B_b(k)=1/(k^2+b), with b=1 for the scalar
chord and b=0 for the gauge chord. The short
arc kernel is

 K_b(q,l,p)=B_b(q-l) S(l) S(l+p).

For gauge exchange cap the ordered fermion
product by gamma_mu on both sides and sum
mu. The scalar relative vertex has internal
coupling +Y, the gauge relative vertex has
-aCf. The sign follows by differentiating
the proper inverse self-energy: scalar
Sigma=-Y integral B S, gauge
Sigma=+aCf integral B gamma S gamma, and
the exact identity partial_m S=-S^2.

At zero external momenta the kernel is
K_b(0,l,0)=B_b(l) S(l)^2. Use the identical
dimensional regulator and subtract these
kernels before integration. Their leading
large-l degree-zero vertex divergences cancel.
Add back the full finite proper MS vertex
at zero momentum, not a freely chosen value:

 Gamma_MS(q,p)/y = ups_MS
   + the signed scalar/gauge integral of
     [K_b(q,l,p)-K_b(0,l,0)].

The S6.139 fixed-scale mass derivative gives

 ups_MS=[Y(J0+2R)-6aCf]/Q,
 |ups_MS| <= (2Y+6aCf)/Q.

This is a proper vertex anchor, not the
canonical bare-Yukawa reference map. In
particular the gauge finite value is -6,
and not the result of resetting mu=m before
the derivative. The proper UV pole is
(-Y+4aCf)/Q times the MS pole; the total
proper counterterm is its negative.

The subtracted kernel is UV convergent for
each external vertex momentum. Its difference
can therefore be evaluated in four dimensions;
finite epsilon-times-pole terms of the proper
subgraph are all retained in ups_MS. The
whole two-loop graph has only a degree-zero
overall divergence after this pairing.
Any remaining dimensional finite product
with that overall pole is a local Phi^4
contact and has zero b2. For positive soft
degrees the bounds below prove convergence
before the regulator is removed.

The subtraction and local anchor are used
together. No unregulated divergent kernel
is assigned a finite value. Other finite
field and parameter conversions, including
S6.138's outer MS conversion, remain separate.
