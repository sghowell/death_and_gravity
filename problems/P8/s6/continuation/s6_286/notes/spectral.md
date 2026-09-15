# Entire positive spectral measure and first-sheet domain

Let A=Delta_ab(0), w=xy. Both triangles have A>0. For n>=mu>0:
A_light-4mu*w=nz+mu(1-z)^2 v^2>=0;
A_heavy-4n*w=nz(1-z)+mu z^2+n(1-z)^2 v^2>=0.
Hence their first-sheet supports begin no earlier than4mu and4n.

This is a proof of the full Feynman integral domain, not a scan of a
finite number of external momenta. Setting y=z/(1-z) gives
A_light/(1-z)^2=mu+n y+n y^2,
A_heavy/(1-z)^2=n+n y+mu y^2.
They are monotone on y>=0 and their lowest values occur at y0.
The first support thresholds are therefore precisely4mu and4n.
Any algebraic pseudothreshold on another sheet is not silently added.

The imaginary part at sigma>4a is obtained with
Im[1/(A-w(sigma+i0))]=pi delta(A-w sigma).
The v root and its full Jacobian give the positive density
2pi c/sigma int_0^ymax dy y^2/(1+y)^3 /
 sqrt[1-4(a+n y+b y^2)/sigma],
where(a,b)=(mu,n) or(n,mu), c=g^2/(16pi^2), and
4(a+n ymax+b ymax^2)=sigma. The endpoint square root is integrable.

Equivalently the Feynman measure pushes forward under sigma=A/w.
The identity
(1/w)t/[sigma(sigma-t)]=1/(A-wt)-1/A
proves the entire once-subtracted positive Stieltjes representation.
The inverse moments are finite; the first is f1prime(0)>0.
This is the generated three-point form-factor measure, not a claim
that it equals the full four-point unitarity density.

Put r=abs(t)/(4mu)<1. Positivity and w/A<=1/(4mu) imply
abs(f1(t))<=abs(t)f1prime(0)/(1-r),
abs(f1(t)-t f1prime(0))<=abs(t)^2 f1prime(0)/[4mu(1-r)].
All higher Taylor coefficients are positive. At t=-tau<=0,
f1(-tau)/(-tau) is between f1prime(0)/(1+tau/(4mu)) and
f1prime(0). The zero-transfer limit is controlled for these massive
graphs. Nothing here provides that limit for internal graviton loops.
