# Complete first light two-point coefficient and finite local conditions

All calculations are for the separate V2S-T1 model. Define the one-particle-irreducible insertion as iPi. The geometric Dyson series has first correction

[i/(s-1)] iPi [i/(s-1)]=-iPi/(s-1)²,

so the inverse denominator is s-1+Pi. This sign is also checked against the Euclidean Hessian and the absorptive rate.

Labelled Wick counting gives the mixed heavy-light bubble factor1, the two-identical-light heavy bubble factor1/2, and the quartic light tadpole factor1/2. With two gHphi²/2 vertices, there are eight routes placing one external light leg on each vertex and four placing both on the same vertex; the expansion and vertex denominators give respectively1 and1/2. The latter is the one-point-reducible heavy-source tadpole, explicitly cancelled by the one-point counterterm in notes/gaussian.md.

Use dimensional regularization and MSbar at mu1. In our conventions

A0_MS(m²)=m²[1-log(m²/mu²)],
B0_MS(s;a,b)=-integral_0^1 Log[xb+(1-x)a-x(1-x)s-i0]dx.

The Wick-rotated Euclidean tadpole is -A0/(16pi²), while the two-propagator integral continues to B0/(16pi²). Consequently the remaining complete first light coefficient is

Pi_MS(s)=-C/(32pi²)+g² B0_MS(s;1,M_H²)/(16pi²).

No heavy-light diagram is omitted because the pure-heavy determinant is constant. The light parity forbids odd-light vertices, and the displayed tadpole and bubble exhaust the two-point topologies at this loop order. Higher-loop diagrams are not excluded.

The heavy-one-point-zero condition and the light on-shell conditions are new, explicit finite conditions for this separate model. Add the local inverse counterterm

-Pi_MS(1)-(s-1)Pi_MS'(1).

It is representable by the ordinary light mass and kinetic operators. It is not an unreported change of any original affine prescription. It makes Pi_OS(1)=Pi_OS'(1)=0 by definition at this order. The renormalized tree parameters retain their S233 values; no finite four-point condition is chosen.

Set F=(1-x)²+xM_H² and alpha=x(1-x)/F. Subtracting the constant and linear terms of the FULL Feynman-parameter function gives

Pi_OS(s)=g²/(16pi²) integral[-Log(1-alpha(s-1))-alpha(s-1)]dx.

The log branch is continued from s1. This representation retains the full nonpolynomial remainder. Omitting either local subtraction fails its respective anchor condition.

Since1<=F<=M_H²<10^198 and log10<7/3,

abs(Pi_MS(1))<[-C+924g²]/(288)<10^-7.

The exact rational exponential partial sum through order8 at7/3 exceeds10, independently checking the logarithm upper bound. The bound pi>3 is used in the denominator. The finite mass coefficient can be nonzero despite the very small low-energy four-point tree.

Moreover0<alpha<=(1-x)/M_H² gives

0<Pi_MS'(1)<g²/(32pi²M_H²)<10^-207.

The un-subtracted formal first pole shift is -Pi_MS(1); its higher-order remainder is not bounded here. The on-shell pole and unit residue are imposed finite-order conditions, not a construction of the exact physical asymptotic field. Full real four-point matching, its finite vertex conditions and all-loop control remain necessary.
