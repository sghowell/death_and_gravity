# S6.233: separate classical model and complete tree remainder

## Model and original target

The model name is V2S-T1. Its only fields are a canonical light scalar phi and a canonical heavy scalar H on Minkowski spacetime with signature+---. Its classical Lagrangian is

L=(partial phi)²/2-phi²/2+(partial H)²/2-M_H² H²/2
  +(g/2)H phi²+(C/24)phi4.

Use the original target lambda=10^-600, gamma=1024*10^-800 and scalar mass1. Set D=2lambda/gamma, M_H²=D+2, g²=gamma D4=1/2^26 and C=-g²(3/D-2/D²). The contact is part of THIS new model; it is not a new renormalization choice for a frozen original theory.

The target is the complete S177 four-scalar on-shell tree, unchanged by the nonconstant low vacuum jets of S182. Matching this target does not identify the models' full independent coefficient functions, higher-field interactions, covariant dictionaries or states. In particular V2S-T1 does not contain the original Proca field or DHOST curvature sector.

## Classical result

The full kinetic matrix is identity. The vacuum mass matrix is diag(1,D+2). An exact positive-square identity proves the potential coercive and nonnegative for D>1, with unique minimum phi=H=0. A second identity at the actual D>2 gives V>=(phi²+H²)/6 globally. These are classical flat-potential and vacuum statements. No quantum vacuum or cosmological stability theorem is inferred.

## Complete massive tree result

The literal cubic and quartic vertex factors give

A_V2S=C+g²[1/(M_H²-s)+1/(M_H²-t)+1/(M_H²-u)].

For s+t+u=4, set xi_i=channel_i-2. The exact shifted identities sum xi=-2 and sum xi³=3stu-8, together with the finite geometric identity, give

A_V2S=2lambda sum xi²+3gamma stu-8gamma
        +gamma sum xi4/(D-xi).

Every original mass and potential term remains. The higher remainder and the heavy pole are explicit.

On the physical angular interval, s>=4, t,u<=0 and abs(xi_i)<=s-2. Below the heavy s-channel pole all denominators in the remainder are positive. The actual original tree obeys A_original>=lambda(s-2)²>0. Hence, with r=(s-2)/D<1,

0<(A_V2S-A_original)/A_original<=6r²/(1-r).

At s<=10^196, r<32/625, and the rational bound6144/370625 is below1/60. This is uniform in every physical scattering angle. It is a comparison of two complete tree amplitudes, not the full error premise in S231.

## Forward and first elastic observables

At t0 and v=s-2, the exact remainder is

16gamma/(D+2)+2gamma D v4/(D²-v²).

Thus b20_tree=4lambda and b21_tree=-3gamma match exactly, while b40_tree=gamma²/lambda and the forward constant shift16gamma/(D+2) are retained. The model supplies the S232 tree positive atom but is not thereby exact unitary.

The full angular amplitude is rational and has infinitely many nonzero even partial waves for s>4. Closed integrals of A and A², with removable massive-threshold limits, are derived without an l0/l2 truncation. With the same labelled-identical convention as S231,

t0_tree=beta integral A/(64pi),
rho_first=beta integral A²/(64pi).

The named window is below the heavy single pole and heavy-pair threshold. Phi parity forbids the mixed phi-H two-body channel. The first elastic coefficient therefore has the stated comparison to the original coefficient. The real amplitude, physical renormalized masses/residues and higher quantum terms remain uncomputed.

## Certification boundary

All17 source inputs are immutable once frozen. Exact algebra, rational margins and independent high-precision diagnostics are distinguished from the written continuum inequalities. The latter are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only full P8 regression uses the audited exact-GCD adapter.

No original action, state, pole, regulator, counterterm or primitive verdict is overwritten. Full quantum V, finite-gravity G and common-parent bounce B remain open.
