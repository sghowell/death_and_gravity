# Exact sourced Proca equation and causal inverse

Use the physical +--- metric and one-forms, raising indices with that
metric where convenient. Write K W = nabla_nu F(W)^{nu mu} and
D W = grad(div W). Variation of the actual action gives

    O W = S,   O=I+zeta K.

Antisymmetry and torsion-free metric compatibility imply div K=0:
the commutator contracts the antisymmetric tensor with symmetric Ricci.
Also K grad=0 because the curl of a scalar gradient vanishes. Therefore
KD=DK=0, including on a curved metric. Define

    N=I+zeta D,   L=I+zeta(K+D).
    ON=NO=L.

The principal symbol of (K+D) is that of the metric wave operator times
the identity. Any curvature endomorphism is lower order; its convention
does not affect this statement. L/zeta is normally hyperbolic for
positive constant zeta. This does not say that O itself is normally
hyperbolic.

On a smooth globally hyperbolic spacetime let G_L be a retarded or advanced
Green operator for L. The existence, uniqueness, two-sided compact-source
identities and causal support follow from Corollary3.4.3 and Definition3.4.1
of [Baer, Ginoux and Pfaeffle, 0806.1036v1](https://arxiv.org/pdf/0806.1036).
The source is a compact smooth section of the one-form bundle.

N commutes with L; uniqueness with causal support gives NG_L=G_L N
on compact sources. G_O=N G_L is a two-sided Green operator for O.
For example O G_O=L G_L=I; G_O O=G_L N O=G_L L=I.
A homogeneous O solution with retarded or advanced support also solves L
and is zero by uniqueness, so these are the unique causal inverses.
One may similarly commute K through G_L by the same support argument.
Consequently

    W=G_O S,    W-S=-zeta G_L K S,
    div W=div S.

The last equation is the sourced constraint. Replacing it by div W=0
would remove a genuine longitudinal source term. Causal support and
smoothness for compact sources are qualitative conclusions; no unweighted
spacetime L2 norm, chosen quantum state or stress bound is supplied.

The flat Fourier matrices are checked in all components, using
p^2=p^T eta p and the mixed-index outer product p p^T eta.
Independent coordinate variation and dense inversion tests supplement
the symbolic identities. Off the mass pole,

    O=(1-zeta p^2)I+zeta p p^T eta,
    O^-1=(I-zeta p p^T eta)/(1-zeta p^2).

This rational expression must be assigned its actual causal, Feynman or
other boundary prescription when used as a distribution. The distinctions
are not settled by algebra alone.
