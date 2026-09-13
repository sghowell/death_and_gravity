# Complete first-loop coefficient intervals

At v=t=0, define b20=[v²t0]A, b21=[v²t1]A and b40=[v4t0]A. These are partial_v²/2, partial_t partial_v²/2 and partial_v4/24 respectively. There is no extra factorial in the t derivative.

The explicit complete n^-2 through n^-5 amplitude gives

delta_b20,approx=g4/(16pi²)[(3ln n-157/45)/n4+(24ln n-7969/315)/n5],
delta_b21,approx=g4/(16pi²)(-3ln n+13/21)/n5,
delta_b40,approx=0.

The first two amplitude orders are constant. The next two have degrees2 and3, so neither contributes a v4 coefficient.

Let R=10^12 g4/(16pi² n6) be the FULL bidisk remainder norm. Cauchy coefficient estimates on radii1/4 give

|e20|<16R, |e21|<64R, |e40|<256R.

Thus delta_b20 and delta_b21 equal the displayed approximations plus their stated errors, while delta_b40=e40. The fixed OS4 contact has zero derivative in each observable; no new finite condition is imposed.

## Actual signs and low-coefficient tolerances

The exact mass n obeys10^197<n<10^198. The elementary bounds2<ln10<7/3 imply394<ln n<462. For the lower bound, the factorial-series estimate e<3 gives e²<9<10. For the upper bound, the positive exponential partial sum through degree8 at7/3 exceeds10. No floating-point logarithm certificate is needed.

The actual rational comparisons, including the COMPLETE error terms, give

3*394-157/45-16*10^12/n²>1000,
24*394-7969/315>0,
1400+11100/n+16*10^12/n²<1500,

and

3*394-13/21-64*10^12/n>1000,
3*462-13/21+64*10^12/n<1500.

Consequently

1000g4/(16pi²n4)<delta_b20<1500g4/(16pi²n4),
-1500g4/(16pi²n5)<delta_b21<-1000g4/(16pi²n5).

The first b20 shift is positive and the first b21 shift negative. Their signs do not rely on a truncated leading term alone.

The unchanged separate tree matches b20=4lambda=2g²/D³ and b21=-3gamma=-3g²/D4, where D=n-2. Using pi²>9,

0<delta_b20/(4lambda)<1500g²/(288n)<10^-203,
0<delta_b21/(-3gamma)<1500g²/(432n)<10^-203.

## Required higher coefficient in the specified truncation

The separate tree has b40=gamma²/lambda=2g²/D5. Its complete first-loop shift obeys

abs(delta_b40)/(gamma²/lambda)
 <8*10^12 g² D5/(pi² n6)
 <8*10^12 g²/(9n)
 <10^-192.

No sign for delta_b40 follows. Nevertheless the complete tree-plus-first-loop b40 lies strictly between(1-10^-192) and(1+10^-192) times gamma²/lambda, so that explicitly truncated coefficient stays positive.

The heavy exchange contribution is retained. There is no subtraction of unknown heavy positive weight, exact Gram-saturation claim or inference that finite matching proves full dispersion/UV properties. All results here concern the complete first-loop corrections with the existing first-order light mass/residue conditions. Omitted loops and exact physical coefficients remain separate.
