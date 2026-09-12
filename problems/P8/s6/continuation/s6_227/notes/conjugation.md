# The actual local remainder becomes a bounded middle operator

Use the exact factorization and BOTH causal inverse identities. In their actual order,

Omega=Z Rloc Y
=A Y+Z A^T+(ZD*)H(DY)+(ZD*)JY+Z J^T(DY)+Z V0 Y.

Only adjacent matching inverse pairs are canceled. In particular, A,H,J,V0 are time-dependent and are not commuted throughY orZ. The two derivative operators are controlled on their correct sides by notes/reciprocity.md.

Put g=sigma-20>0 and C0=5/2. Using the actual coefficient bounds gives

||Omega||
<=2(110m^2) C0/g^2
 +(185m^2)C0^2/g^2
 +2(286m^2)C0^2/g^3
 +(160m^4)C0^2/g^4
=6825m^2/(4g^2)+3575m^2/g^3+1000m^4/g^4.

Every term is a bounded causal operator on the SAME weightedL2_tH^r space, uniformly inq and for every realr. The original q terms were absorbed into the exact curvature factors before taking norms. No spatial derivative loss is hidden in a bound on a multiplication byq.

For sigma=sigma_M=2m exp(32+8M), M>=0, the elementary boundexp32>=1+32+32^2/2=545 givesg>=100m at the unchangedm1000. Hence

||Omega||<=6825543/40000000<1/4.

This is a deliberately loose rational upper bound. The exact bound decays as sigma grows, but no unweighted smallness is inferred.

## Domain extension

Initially perform the identity on smooth causal functions with compact spatial momentum, as an equality of finite-order distributions including the initial boundary. Each of the six right-hand terms extends boundedly to weightedL2_tH^r by the displayed estimates. This definesOmega on that completed source space.

The factor identities then hold by approximation in the finite-order Sobolev-valued distribution topologies. For every x in the completed source space and s=Yx, B*Omega x=Rloc s as a causal distribution. In particular B*Omega B=Rloc on this pullback domain. The complete inverse uses the additional middle/output graph conditions in notes/inverse.md; no claim of a maximal unrestricted distributional graph is made. No source-time derivative atom is lost in the extension. This is a new proved mapping of the ACTUAL finite local remainder. It is not a statement that every derivative-losing weak response admits the same conjugation estimate.
