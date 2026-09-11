# A fixed mass reference, with a named high-field extension

Write chi=f_R(Phi), R=10^300, and Y=y^2. The active
opposite masses are mF +/- y chi, with three colors
per sign and N=6 total active color/flavor copies.
The twelve inert flavors have constant mass; all their
vacuum constants are included in the same fixed
zero-field vacuum reference, not discarded.

The paired MS potential is

    -N/(2Q) sum_{sigma=+,-} (mF+sigma ychi)^4
        [log((mF+sigma ychi)^2/mF^2)-3/2].

Its zero-field second derivative is fF(0)=4NYmF^2/Q.
The earlier on-shell calculation fixes fF(1), including
the full tadpole and momentum kernel. Its complete
positive moment series gives

    fF(1)-fF(0)=(2NY/Q)[2/3-sum_{n>=1} a_n/mF^(2n)],
    a_n=3(n!)^2/[n(2n+3)(2n+1)!].

For mF>=36 this lies strictly between zero and
4NY/(3Q). No cancellation of huge rounded decimals
is used. The external light pole invariant is one;
restoring dimensions places m_Phi^2 in that upper
bound and in the corresponding energy allowance.

Define GY14-SAT8-MR by replacing, at formal loop one,
the finite potential counterterm -fF(1)Phi^2/2 by
-fF(1)chi^2/2. The coefficient is the inherited pole
anchor, not fitted to a background energy. The MS
field and quartic prescription gets no extra finite
subtraction. The minimal fermion pole functional
already uses chi; the fixed vacuum reference remains.

This is an explicit additional high-field EFT choice,
not a rewrite of a frozen functional. In particular,
the unsaturated -fF(1)Phi^2/2 diverges on the unbounded
clock even though the mass and determinant stay
bounded. A bounded mass alone was not enough.

The first finite difference follows from

    chi^2=Phi^2-Phi^10/(4R^8)+5Phi^18/(32R^16)+... .

Its loop-one degree-ten vertex has weighted excess
(10-2)/2+1=5. Connected graphs obey

    L_eff=1-E/2+sum_v[(n_v-2)/2+j_v]+map_degree.

All other weights and the physical-map extra degree
are nonnegative; forest contraction preserves the
weight. For E=0,1,2,4 the loop lower bounds are
6,6,5,4, respectively. Thus all named through-two-loop
pole/residue/four-point, vacuum/H and low-cut data
remain identical. Degree-ten data can change at one
loop; the complete quantum functional is not equated.

Use chi_D and R_D=mu^-epsilon R and the inherited
regulated mass-reference coefficient before Laurent
finite products. The degree8k+2 coefficient has
dimension2-8k+8k epsilon as required. No untracked
positive-epsilon pole product is dropped. The
physical unit-residue source dictionary remains
locally the inherited one. This energy is specified
in the MS field coordinate, not the whole canonically
re-expressed scalar/heavy stress.

The extension has one formal loop power. The complete
classical zero-fermion action match is unchanged, but
quantum target matching and a clock solution do not
follow.
