# Full two-leg memory and distinct contact regulator limit

The removed memory region is max(|k|,|l|)>K with P=k+l.
For |P|<=K/2, this forces BOTH |k| and |l| above K/2.
The previous full pair-weight calculation therefore gives

    36 Cref^2 B^2 L J10_tail(K/2).

Since nu>=|k|/Amax,

    J10_tail(K/2)<=64 Amax^10/(7pi^2 K^7)<100/K^7.

After L<2(1+|P|^2), the recorded low-transfer numerator is
72 Cref^2 B^2*100 for K^-7. Using K>=m converts it to a
1/K numerator below 3e16. The independent numerical tests integrate
the complete removed union, including its moving angular boundary;
they do not test only a radial subset.

For |P|>K/2 use the complete internal integral and
L/(1+|P|^2)<1/(100K). The resulting 1/K numerator is
(36/100) Cref^2 B^2 J10<1e10. This keeps high-low pairs and
all external momenta.

Add the full one-mode contact tail from contact.md, whose
1/K numerator is below 1. The exact recorded sum is below 1e17.
Thus

    |R_alpha-R_alpha,K|<1e17 M[D]M[Gamma]/K.

The complete memory plus contact bound is below 2e12 in those norms.
The quantities are exactly zero if either test is zero. K is removed
as a computational regulator, not declared a physical cutoff.
No equality between the one-leg and two-leg regions is imposed.
