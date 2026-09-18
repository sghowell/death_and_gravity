# Massive kernel and original convention audit

The original scalar stress, including its trace, mass term and vertex i,
is compared directly with frozen S304. Massive and massless Feynman
prescriptions are preserved by the map in source.md; changing metric
signature is not complex conjugating a loop. Frozen S295 independently
checks the leading Born current. The pair invariant and its Lorentz
derivative contraction both change sign in the source convention map.

For t>=1, use r(t)=2int_0^1dv/[t+1-(t-1)v^2], with denominator>=2.
It equals acosh(t)/sqrt(t^2-1) on t>1 and supplies its removable endpoint.
Differentiation on this fixed compact interval gives r(1)=1,r'(1)=-1/3,
r''(1)=4/15. The signed exterior kernel f_signed(z)=sign(z)f(|z|),
f(t)=2(2t^2-1)r(t), has derivative f'(|z|), even across pair orientations.
At t1, f'=22/3 and f''=16/5.

The equivalent hyperbolic integrand for f' is
(4t^2+8t*cosh(u)+2)/(t+cosh(u))^2>0. Since r'<=0 and r<=1,
0<f'<=8t<=56 on[1,7]. Fixed-interval derivative bounds give
|r'|<=1/3 and |r''|<=4/15. Hence
|f''|<=8+16*7/3+2*97*4/15<161. No finite sample supplies this proof.

For same-time pairs, c=t(2t^2-3)/(t^2-1)^(3/2),
c'=3/(t^2-1)^(5/2). The outgoing pair has t>=29/16, incoming>=17/8.
The identity4(t^2-1)^3-t^2(2t^2-3)^2=3t^2-4 shows0<c<2 there;
squaring the positive c' bound gives c'<1.

In the forward hard-angle LIMIT at E5/4,spatial momenta+/-3z/4,
q=(1,1,0,0),A=diag(0,0,1,-1)/sqrt2, the leading current and each
real-pair group vanish but C=-3281sqrt2*i/(6000pi) is nonzero.
This is a kinematic coefficient limit, not a value assigned to the
divergent forward hard amplitude. The imaginary phase cannot be discarded.
