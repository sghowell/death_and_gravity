# Complete scalar moments and derivative bounds

## Absolute moments and scalar master bounds

Let c0=1/32 and J[p,r]=integral_0^1 t^r/(c0+n0 t)^p dt.
The exact antiderivatives give

    J[1,0] < 500/n0;
    J[2,0] < 32/n0;
    J[2,1] < 500/n0^2;
    J[3,1] < 16/n0^2;
    J[3,2] < 500/n0^3.

Here log(1+32 n0)<200 log 10<500. No sampled floating-point inequality
is used in this logarithmic estimate.
Using Jacobian product <4 and the triangle measure modulus <=2 yields

    |C_off| < 8*64*500/n0 = 256000/n0.

Using box measure modulus <=4t yields

    |D_ordered_off| < 16*64^2*500/n0^2 = 32768000/n0^2.

Weighted moments can be bounded with explicit powers of
|xi|,|1-xi|<=2, |z|<=2t, and |1-z|<=2(1-t).
Their numerator weights must be displayed; an unweighted scalar integral
bound alone does not bound an arbitrary weighted or tensor integral.


Set U=c+n t. Exact primitives for J10,J20,J21,J31,J32 respectively are

 log(U/c)/n;
 -1/(n U);
 [log(U/c)+c/U]/n^2;
 [-1/U+c/(2U^2)]/n^2;
 [log(U/c)+2c/U-c^2/(2U^2)]/n^3.

Differentiate each before taking endpoints. With y=c/(c+n) in(0,1),
the nonlogarithmic J32 residual is-(1-y)(3-y)/2<=0, and the J31
residual is-y(2-y)/2<=0. This proves the stated upper bounds.
The positive exponential partial sum through degree6 proves
exp(7/3)>10, hence log10<7/3 and log(1+32n)<500.

For every external multiindex alpha, the fixed polydisk and Cauchy
formula give alpha!*4096^|alpha| times the master bound. Additional mass
derivative r multiplies by r!*(128/n0)^r. Joint analyticity here is for
the continued physical boundary germ; there is no assumption that the
principal first-sheet function agrees on the other side of the cut.
