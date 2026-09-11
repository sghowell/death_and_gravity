# Uniform removed-momentum bound

The memory cutoff removes max(|k|,|l|)>K. Split only the
EXTERNAL momentum P=k+l at K/2; it is not cut off.

For |P|<=K/2, the removed union forces both internal legs
>K/2. The complete memory-difference tail is bounded by

    36L(Cr Ce J4_tail(K/2)+Ce^2 J10_tail(K/2)),
    J4_tail(K/2)<=A^4/(pi^2 K)<A^4/(9K),
    J10_tail(K/2)<100/K^7.

Use L<2(1+|P|^2) and K>=1000. The low-external1/K
numerator is
72[Cr Ce A^4/9+Ce^2*100/1000^6]<1e26.

For |P|>K/2, it is wrong to force both legs large.
Keep the full internal bound and use
L/(1+|P|^2)<1/(100K), established for all K>=1000.
Its1/K numerator is
(36/100)[Cr Ce J4+Ce^2 J10]<1e21.

The exact sum of these two memory numerators is still<1e26.
The contact has one internal momentum and contributes
<1e14/K^2. Adding its exact numerator divided by1000 gives
a complete1/K coefficient<1e27.

The same unit-time and Fourier Cauchy inequalities as in
notes/memory.md prove

    |R-R_K|<1e27 M[D]M[Gamma]/K.

The estimates apply with the retarded factor present. They
justify an absolutely convergent continuum difference and
its quantitative error, not a renormalized limit of the
reference term or a physical momentum cutoff.
