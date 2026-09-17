# Probability subtraction and the necessary state-transport connector

For a given physical polarization pair, write the two state-correct
soft faces U,V and their common overlap O as in the formulation.
S313 subtracts the amplitude B2=U+V-O. The exact complex identity is
 |B2|^2-|U|^2-|V|^2+|O|^2
 =2Re[(U-O)*conjugate(V-O)].
Consequently the probability inclusion-exclusion density
 r2=|F2|^2-|U|^2-|V|^2+|O|^2
equals the S313 signed density plus this extra interference.

The S313 bound |G-G00|<B*W, with actual B<10^-350, has continuous
single-soft axes. Passing to those axes gives
|G(0,b)-G00|<10^-350*b and its symmetric counterpart. Thus
|U-O|<10^-350/a, |V-O|<10^-350/b.
The extra density is absolutely bounded by2*10^-700/(ab).
Multiplying by ab da db and integrating the triangle a+b<=x yields
10^-700*x^2. The same full polarization/angle/Bose factor as S313 is
1/(8pi^4)<1. Therefore
 TV(R2;x)<2*10^-725*x+2*10^-652*x^2.
This is finite uniformly in a common lower seed cutoff.

## Exact density reconstruction

After polarization sums let K_a(sigma)=sum|J_a(sigma)|^2/kappa.
Let f1a,f1b,f2 be the complete phase-weighted physical tree densities,
and K_a0,K_b0 the elastic single-soft densities. Then
 r1a=f1a-K_a0, r1b=f1b-K_b0,
 r2=f2-K_a(sigma_b)f1b-K_b(sigma_a)f1a+K_a0*K_b0.
The elastic leading two-real density plus the S309 one-marked class is
 K_a0*K_b0+K_a(sigma_b)r1b+K_b(sigma_a)r1a.
Adding r2 is NOT f2. The exact missing density is
 C2=[K_a(sigma_b)-K_a0]K_b0+[K_b(sigma_a)-K_b0]K_a0.
Adding C2 makes the entire expression identically f2.
The mismatch without C2 is a nonzero general polynomial, and actual
original recoil currents have nonzero state transport.

C2 by itself is not another finite signed two-real measure: a state
change of order b in its first current still leaves a single-soft
logarithm in a. Its same-state virtual contribution must be retained.
The connector
 E1=int dB1[P_sigma_b(x-b)-P0(x-b)]
provides the prescribed completion. At a common positive seed cutoff,
its one-additional-real term is C2 after symmetric relabeling; its
no-additional-real virtual terms are the associated difference of
the same-state virtual factors. The entire-sum difference has a
finite cutoff/removal limit by notes/series.md.

## Meaning and limit of tree-sector matching

The statement of matching is about differential physical D4 TREE
densities BEFORE integration. It can be audited with formal radiation
tags: one-real seeds have degree1, two-real seeds degree2, and each
additional real emission degree1. For a tree coefficient take the
zero-virtual part of every factor. Hold the original hard Born and
all numeric coefficients fixed; this is not a variation of kappa
in the hard exchange. Degree0 is the normalized Born, degree1 is
K0+r1=f1, and degree2 is the explicit polynomial identity above.

This is not a fixed-N dimensional-regulator limit inside a divergent
integral. The completed soft factors are defined by the separately
proved entire sums and their ordered limits. Nor does the real
identity determine the unknown finite hard five-point loop or
four-point higher-loop amplitude. Their evanescent dimensional
pieces remain separate matching obligations.

Six original full434/47-face states independently replay the identity
between the probability and amplitude subtractions, including both
phase ratios and the factor1/2 from the two unnormalized plus tensors.
No sampled amplitude replaces the general density identity.
