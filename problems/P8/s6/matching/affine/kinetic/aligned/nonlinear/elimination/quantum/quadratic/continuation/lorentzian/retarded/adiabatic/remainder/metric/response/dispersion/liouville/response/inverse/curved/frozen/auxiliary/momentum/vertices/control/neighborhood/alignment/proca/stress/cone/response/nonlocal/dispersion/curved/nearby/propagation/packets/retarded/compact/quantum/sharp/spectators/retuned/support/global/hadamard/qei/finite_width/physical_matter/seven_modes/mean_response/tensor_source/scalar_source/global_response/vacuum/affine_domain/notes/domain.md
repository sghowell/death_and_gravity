# A new connected analytic domain with the full rank factor

For even n>=1024, h=(1+u^2)^3, define

    T=X^n/[X^n+(1-X)^n],
    w=nX^2 exp(-nX^2),
    B=T+(1-T)w, S=(1-T)(1-w),
    R=1+B(X-1)/h.

The denominator of T is positive for all real X. Both T
and w lie in [0,1], and w<=1/e<3/8. Thus 0<=B<=1.
Use all real u and -1/(4n)<X<6/5.

For negative X, T<1/(16n), w<=1/(16n), giving
R>1-5/(32n). For 0<=X<=1/4, T<=3^-n<1/40,
so B<2/5 and R>3/5. On 1/4<=X<=1/2, put y=nX^2>=64.
Since y exp(-y) decreases there, w<=64 exp(-64)<1/100.
Together with T<=1/2 this gives B<101/200 and R>497/800.

On 1/2<=X<=1/2+1/(8n), the logarithm of the step odds is
bounded above by

    n log[(1+1/(4n))/(1-1/(4n))]
      <=1/[2(1-1/(4n))]<1.

Hence the odds are below e<3, T<3/4 and B<301/400,
so R>499/800. For the remaining positive X, use
R>=min(1,X). Consequently

    R>=1/2+1/(8n),  2R-1>=1/(4n)>0,  R<=6/5.

These are per-order bounds. The tensor floor stays above
1/2 throughout the family, but the full quotient-factor
lower bound is not uniform as n tends to infinity.

The earlier S6.108 polynomial switch gives, at u=0,

    R(0)=R(1)=1,
    R(1/4)=1/4+(3/4)(15/16)^n<1/2,
    R(1/2)=1/2+(1/2)(3/4)^n>1/2.

Continuity therefore forces two distinct R=1/2 loci,
one in (0,1/4), one in (1/4,1/2). A positive reduced
tensor factor does not imply a regular full affine quotient.
The old sourced solution remains compatible at rank loss;
that does not restore uniqueness or prove a propagating ghost.

For the new switch S has a zero of order n at X=1, whereas
near X=0 it is 1-nX^2+O(X^4)+O(X^n). R, its exceptional
Ia completion and the lower scalar therefore preserve the
same finite clock jets and massive vacuum quartics as
S6.108, while avoiding that particular rank defect.
