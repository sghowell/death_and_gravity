# Literal all-D TT tree and complete physical polarization sum

Write n=D-2 for the transverse Euclidean dimension. In the pair center
of mass, let r be the scalar momentum transverse to the graviton axis
and delta=e^2-p_z^2. For arbitrary symmetric traceless transverse H,K,

 A(H,K)=[delta tr(HK)-4r^T HK r
                     +4(r^T H r)(r^T K r)/delta]/kappa
       =tr(H Q K Q)/(kappa delta), Q=2rr^T-delta I_n.

The proof starts with the exact Einstein GammaGamma density at
g=eta+2h/sqrt(kappa). Every first/second connection and density term
multilinear in the three distinct graviton legs is retained. Both
scalar exchanges and the complete Phi^2 h^2 contact are added.
The internal scalar stress through the harmonic projector is
W=a b^T+b a^T+eta*2mu/(D-2).
The discarded Einstein divergence has zero total external momentum,
as in S280; no cosmological canonical boundary is being discarded.

Rotational covariance, parity and the two TT trace conditions leave
exactly the three bilinear invariants displayed above. Three controls
span them for r!=0: H=K=offdiag(2,3), offdiag(1,2), diag(2,-1,-1),
taking r along the first transverse direction.
The r=0 identity follows by rational continuation.

Here is the dimension degree bound required to use finite component
evaluations. Each cubic term contains one W, one H and one K.
Every nonzero momentum and TT component is supported on the displayed
external subspace. A spectator index can therefore occur in at most
one metric trace: two independent spectator traces would require a
second W or a spectator component of H,K or momentum, all absent.
The dependence on the number of spectators is affine. W's only
dimension denominator is D-2, retained symbolically throughout.

The code evaluates five and six COMPONENT dimensions with a symbolic
trace dimension D, then forms value5+(D-5)(value6-value5).
The resulting whole tree equals the target for symbolic D.
The degree argument makes this a full rational identity, not an
unqualified numerical scan over dimensions. Independent tests also
evaluate the full literal graphs in D7 and D8.

This agrees with the D-dimensional covariant organization in
[Bjerrum-Bohr et al.,section4](https://arxiv.org/html/1908.09755).
The original action, not its printed prefactor, fixes normalization.

For the symmetric-traceless projector
P_ij,kl=(delta_ik delta_jl+delta_il delta_jk)/2-delta_ij delta_kl/n,
the rank is n(n+1)/2-1 and P^2=P. The complete two-projector sew is

 [(tr QL QR)^2+tr QL QR QL QR]/2
 -2 tr QL^2 QR^2/n +tr QL^2 tr QR^2/n^2,

divided by kappa^2 deltaL deltaR.
This follows by expanding P=I_sym-|I><I|/n and cyclic traces.
The code independently expands every dyad trace word and checks full
matrix projector contractions. No dilaton or antisymmetric tensor
square is used. Its n2 limit is the entire S282 helicity-phase numerator.
