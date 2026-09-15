# Positive core probability and fully coupled regulator comparisons

Write SH,EH,BH,DeltaH for the complete order-zero Hamiltonian bounds
of notes/operators.md; S1,E1,B1,Delta1 for each first homogeneous
derivative; and B2 for each second derivative. Use SF,EF,BF,DeltaF
and BF1 for the centered physical-volume symbol. All are exact
positive rationals, including each remainder.

Each of the four S6.273 solutions has its own actual scalar phase
factored from the interaction state phi. The original g0 force stays
in the classical equations. Factoring does not replace a live
background by the reference or erase an energy derivative.

## Positive POVM probability

The SAME original seed has strictly positive outside-core coherent
probability bounded by tail=1e-4000. In the unchanged whitened frame
its exact Gamma48 tail is exp(-t) sum(k=0..47)t^k/k!, t=R^2/2.
S6.270's complete Chernoff estimate gives tail<exp(-1e39)<1e-4000,
using e^3>10 and1e39>12000. No new tail distribution is chosen.
For any allowed live Y path,
bounded unitarity gives ||phi(u)-psi0||<=BH|u|. Its outside-core coherent
effect satisfies0<=E<=I, hence

<phi,E phi> <= tail+2 BH T <1e-6.

This is a positive POVM probability, not sharp simultaneous phase
support, a zero tail or a projection. Both effect and state are
transported by the SAME fixed reference unitary.

## Direct cutoff and ordering differences

At fixed Y, the difference of original cutoff symbols has a P6
coherent part of norm<=2SH whose action on the seed is
<=2SH sqrt(tail). The FULL Weyl remainder has norm<=2EH even though
it is not core-supported. In a unitary Duhamel comparison its action
on the evolving comparison state can be bounded using displacement
BH|u|. Integration gives

ecut = 2SH T sqrt(tail)+2EH T+BH^2 T^2.

For any homogeneous derivative, the corresponding expectation
difference on paths within BH T of the seed is at most
2S1 tail+2E1+4B1 BH T. The full canonical force conversion gives

Fcut = 1e104/kappa * (2S1 tail+2E1+4B1 BH T).

For ordering, the direct terms are
eorder=T DeltaH and Forder=1e104 Delta1/kappa.
The same conservative bounds cover calibrated and Weyl cutoff
comparisons without changing either definition.

## Solve BOTH feedback rows

Let y be the supremum of the five-coordinate Y difference and z the
supremum of phase-factored state difference. With the original complete
homogeneous Lipschitz bound Lhom=1e112, put

Bq=1e104 B1/kappa, Lq=1e106(B1+B2)/kappa,
a=T(Lhom+Lq), b=2T Bq, c=5T B1.

The exact coupled integral equations give

y <= a y+b z+T Fdirect,
z <= c y+edirect.

Both cross rows are necessary; the two backgrounds are not assumed
identical. Their inverse denominator1-a-bc is>99/100, giving

y <= (b edirect+T Fdirect)/(1-a-bc),
z <= c y+edirect.

Using the above direct terms yields the stated cutoff and ordering
ceilings. M1 is still reconstructed from its full cyclic quadrature;
no new sup-norm claim for that separate coordinate is made.
The uniqueness proved in S6.273 identifies these with its existing
solutions, not trajectories of a replacement model.

## Compare the complete physical volume

The complete parameter derivative of
exp(3alpha)[F0(Y)I+f(Y)] has row sum bounded by

12+4*2*12*NZ+6BF+10BF1 <1e8.

Here |a^3|<2, |U_N|<12, |N_z|<1e4, all four N-dependent homogeneous
directions remain, and centered-symbol derivatives use BF1. At fixed
Y, a state change cancels the scalar identity F0 I exactly; its bound
is4BF z, not a multiple of the uncentered identity norm.

At fixed ordering, the readout operator difference on the evolving
state is bounded by
4SF tail+4EF+8BF BH T. The full heat remainder and actual displacement
remain. At fixed cutoff, the ordering readout difference is2DeltaF.
Consequently

Vcut <= 1e8 ycut+4BF zcut+4SF tail+4EF+8BF BH T <1e-514,
Vorder <= 1e8 yorder+4BF zorder+2DeltaF <1e-566.

The five-coordinate/state ceilings are respectively
1e-580/1e-12 and1e-630/1e-65.

Distinct live Y paths may give large differences in their scalar g0
phases. Restore each solution's OWN exp[-i integral g0(Y)] and the
SAME fixed reference unitary. The vector-norm comparison is only for
phi, not those un-factored vectors. All physical volume and POVM
expectations remain phase invariant. Every scalar-center classical
force, heavy mode and original covariance contact is retained.
