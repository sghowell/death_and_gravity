# Finite original-parameter bound

Do not substitute an asymptotic limit for an inequality.
The exact original expression can be written

    c_core=-(69n^13+R(n,log n))/(315n^4(n-2)^4(n-1)^9).

Every term in R has n-degree at most12 and log power0 or1.
For n>=10^6 and0<=log n<=462, the exact rational tail budget is
the sum over ALL those terms of

    abs(coefficient)*462^logpower/(10^6)^(13-degree).

Thus |R|/n^13<=tail. With d0=(1-2/10^6)^4(1-1/10^6)^9,
the denominator divided by315n^17 lies between d0 and1.
Exact rational comparisons prove

    (69-tail)/315>1/5, (69+tail)/(315*d0)<1/4.

Therefore -1/(4n^4)<c_core<-1/(5n^4) throughout this finite domain.

Original10^6<n<10^198 is checked exactly. The degree8 positive Taylor
polynomial for exp(7/3) exceeds10 by rational arithmetic, proving
log10<7/3, hence log n<462. No floating logarithm comparison is needed.

Restore16pi^2>144 and the prior positive Born bound A0>4g^2/n^3:

    chi_core<0,
    |chi_core|<g^4/(576n^4),
    |chi_core|/A0<g^2/(2304n)<10^-207.

This is a coefficient bound only, not an above-threshold Taylor-error
bound, a total physical radiation bound, or a constraint on extra chi.
