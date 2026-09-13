# Complete current source germs and lapse derivatives

The entire S174 affine reduction is CD+M1 plus

    kappa sqrt(-g)[-zeta F(W)^2/4+(W-S)^2/2],

with the S238 H sector and both fixed reference-profile corrections retained.
The56 algebraic affine-complement variables and four projective directions
are not thereby assigned a nonlinear quantum determinant prescription.
The unchanged constants are kappa=10^800, zeta=10^-6 and the original H mass.

Write h=(1+u^2)^3, T=X^1024/[X^1024+(1-X)^1024] and let b(X) denote the
original bump. The exact clock expressions are

    Rclock=1+(X-1)/h,
    Fclock=Ftree+A+Bprofile(X-1),
    A=-Ptotal, Bprofile=-(rhototal+Ptotal)/2.

rhototal and Ptotal are names for the original fixed normalized Proca
profile plus the full S240 H profile divided by kappa. They are fixed
functions of u, not live quantum means. Their exact bindings, the complete
R/F differences, both heavy corrections and all three vacuum constants
are stored in parent.source_germs. Explicitly,

    Rfull-Rclock=(1-T)[-(X-1)(1-b)/h]+delta_R,
    Ffull-Fclock=(1-T)[exp(-u^4)(1-b)(Fvac-Ftree)
                     -(A+Bprofile(X-1)+PVtotal)]+delta_F.

Both heavy corrections have the complete (X-1)^8 exponential switch factor.
The ordinary complement switch has order1024 and denominator1 at X=1.
Consequently the entire differences have an eighth-order clock zero and
their fixed time derivatives preserve it. This licenses the finite clock
jets used here; it does not replace either full function away from the tube.
The exact test replaces only the common inverse switch denominator and
checks an exact roundtrip before expanding. Replacing an arbitrarily
flattened occurrence of the whole switch would not be a valid identity.

The physical metric and affine-clock ADM metric obey

    g=C ghat+Ddis du^2, C=R^-1/2, Ddis=(1-C)/X.

The lapse N and contravariant shift are unchanged, X=N^-2 and
a_phys=R^-1/4 a_hat. The complete scalar ADM coefficient families are

    M=R^1/4, U=R^-3/4, C3=R^3/4/2, Cchi=R^-1/4,
    B=-U sqrt(X) R_u/2-I,
    Fhat=U[F+9X R_u^2/(16R)]-I_u/N,
    I_s=3 U X R_u R_X/(2R), s=1/N, I(u,1)=0.

Cchi is not U: the matter gradient also contains the inverse spatial
metric. The original S174 chart calculation cancels lapse velocities and
spatial lapse-gradient terms before restriction. One must retain the
primitive in later physical variation. In particular

    I_N=0, I_NN=-3 h'/h^3

on the clock, while I_NNN and I_NNNN are also nonzero. parent.lapse_jets
lists every lapse derivative0..4 for all six families. The primitive is
formed by integrating its derivative jet with the fixed zero constant.
R_X in that derivative loses one clock order, while integration restores
the required primitive order. The eighth-order full difference therefore
suffices for every coefficient and mixed time/lapse derivative in this
degree-four background/fluctuation calculation.

The full chart source is a one-form along du, with normal component

    Snormal=(R-1)(Khat-3Hclock(u)/N), Hclock=4u/(1+u^2).

Hclock is a fixed coefficient, NOT the varied background Hhat. On the
reference S0=S1=0, but S2=-4 delta*n*(3vdot-b), delta=1/(2h).
Thus -W.S begins cubically and S^2/2 quartically. Those terms leave the
reference Gaussian unchanged but contribute to its first and second
background vertices. The complete source square is mandatory.
