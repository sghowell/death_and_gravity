# The complete analytic source on a real clock neighborhood

No leading-germ replacement is used for the estimates.
On 9/10<=X<=11/10 set q=(1-X)/X; |q|<=1/9. For the fixed even n=1024,

    T=1/(1+q^n), w=nX^2 exp(-nX^2), B=T+(1-T)w,
    R=1+B(X)(X-1)/(1+u^2)^3.

Thus 0<=T,w,B<=1 and 1+q^n>=1. The exact first and second derivatives
of T, with q'=-1/X^2 and q''=2/X^3, are bounded in coefficients.py
by explicit rational numbers below1. Their rational identities
extend from positive q to negative q for the frozen integer exponent;
both sides of X=1 are independently tested against the full function.

For the full exponential, y=nX^2>0 and exp(y)>y^3/6, exp(y)>y^4/24.
Applying these to w'=2nX(1-nX^2)exp(-nX^2) and
w''=2n(1-5nX^2+2n^2X^4)exp(-nX^2) gives rational bounds below1.
Consequently |B'|<=2 and |B''|<=4. The complete coefficients satisfy

    |R-1|<=|X-1|, |R_X|<=6/5, |R_XX|<=22/5,
    |R_u|<=3|R-1|, |R_Xu|<=18/5, |R_uu|<=15|R-1|.

The last three follow from the EXACT u dependence (1+u^2)^(-3),
not a bounded-u approximation. In particular
R_u=-6u(R-1)/(1+u^2) and
R_uu=(-6+42u^2)(R-1)/(1+u^2)^2.
The existing global R>1/2 bound remains valid.

For E=Boxu-3H(u)X and Z=u^mu Hess(u)_mu_nu u^nu, the complete
S6.174 source is exactly

    S_mu=partial_mu u (R-1) Q,
    Q=E/X+3R_u/(4R)+C Z,
    C=-1/X^2+3R_X/(2RX),

where the last denominator is2*R*X, not a new coefficient symbol.
The full bounds give |C|<=424/81<6, |C_u|<15 and |C_X|<32.

All three defects R-1,E,Z vanish on the unchanged clock, and R_u
has the R-1 factor. This proves the double zero of S without
discarding a Taylor remainder. It also retains spatially nonclosed
sources. For psi=epsilon f with f=(t+t^2/2)c(x) near the bounce,
direct differentiation of the complete source gives

    S2_0=3c(x)^2, S2_i=0, partial_i S2_0=6c partial_i c.

A compact smooth cutoff can realize these local jets. Conversely,
a homogeneous source is closed and has the compatible retarded
solution W=S exactly. That algebraic homogeneous comparison is not
used to insert a finite support-volume norm on R^3.
