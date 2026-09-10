# Exact local jet and substitution majorants

Expand R and L4_trunc as polynomials in independent field
jets, and take the sum of absolute numerical coefficients
separately in each actual positive coupling. The exact groups are

    C_R=c/6+13lambda+(434/3)gamma,
    C_L=(13/6)c+81lambda+729gamma.

This grouping is a safe triangle majorant even when different
coupling contributions cancel inside a monomial. Native checks
verify cubic and quartic field homogeneity. The selected exact
parameters give C_R<10^-403.

Each coordinate derivative of a cubic field-jet polynomial
differentiates three factors. Its coefficient sum therefore
grows by at most three. Repeating this argument bounds every
derivative of order j<=6 by 3^j C_R, on jets of phi through
order four+j. The independent first derivatives are also
checked directly. Hence the ten-jet cube controls the
substitution into the sixth-derivative quartic action.

If every original jet through order ten has absolute value
at most B, then every mapped jet needed by L4_trunc changes
by at most C6 B^3, where C6=3^6 C_R<10^-400.
The exact free-action remainder has bound

    [(4*(3C_R)^2+C_R^2)/2] B^6
      =(37/2) C_R^2 B^6.

For a quartic monomial the substitution changes the product
by at most B^4[(1+C6 B^2)^4-1]. Summing coefficients yields

    E_field(B)=(37/2)C_R^2 B^6
               +C_L B^4[(1+C6 B^2)^4-1].

At B=1 the free part is below 10^-805 and the total is
below 10^-800. All these are exact rational calibrations.
The derivative bound is a finite-polynomial induction,
not an assumption about arbitrary inverse differential maps.

This is initially pointwise. It does not by itself bound a
spacetime integral on an infinite domain. The common-class
proof supplies the separate L2 and Fourier-support
hypotheses needed to integrate every error term.
