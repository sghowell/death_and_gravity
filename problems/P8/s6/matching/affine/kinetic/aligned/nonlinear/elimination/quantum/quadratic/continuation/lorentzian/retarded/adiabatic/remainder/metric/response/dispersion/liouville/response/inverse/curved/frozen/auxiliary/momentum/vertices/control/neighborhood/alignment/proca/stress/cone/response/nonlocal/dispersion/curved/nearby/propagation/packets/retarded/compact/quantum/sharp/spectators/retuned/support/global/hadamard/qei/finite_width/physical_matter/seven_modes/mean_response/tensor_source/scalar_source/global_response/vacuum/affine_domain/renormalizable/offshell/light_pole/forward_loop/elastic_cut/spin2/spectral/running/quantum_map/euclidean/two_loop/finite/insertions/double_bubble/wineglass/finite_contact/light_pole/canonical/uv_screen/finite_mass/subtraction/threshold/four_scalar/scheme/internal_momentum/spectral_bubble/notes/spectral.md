# A direct spectral identity, without an assumed high-energy contour

Write T=4mF^2>1, A=x(1-x), tau=mF^2/A>=T and C=2NY/Q.
For interior x the inherited parameter kernel is

    g_x(s)=-(T-s) log(1-s/tau).

Its complete on-shell subtraction satisfies

    -g_x,R(s)/(s-1)^2
       = integral_tau^infinity (u-T)/[(u-1)^2(u-s)] du.

To prove this directly, multiply the right side by -(s-1)^2 and
differentiate twice. The second derivative of
(s-1)^2/[(u-1)^2(u-s)] is 2/(u-s)^3. Its elementary u integral gives

    -1/(tau-s)-(tau-T)/(tau-s)^2 = g_x''(s).

Both mass and residue anchors vanish at one. These properties identify
the analytic functions on Re(s)<T. The integrand and its differentiated
forms have integrable large-u bounds locally uniformly on this domain;
no contour-at-infinity hypothesis is assumed. The x endpoints A=0
contribute zero by the tau-to-infinity limit.

For fixed u, the allowed x interval has endpoints
(1 minus/plus beta)/2, beta=sqrt(1-T/u). Its length is beta.
Thus the complete density is (u-T) beta=u beta^3 and

    -f_R(s)/(s-1)^2 = integral_T^infinity w(u)/(u-s) du,
    w(u)=C u beta^3/(u-1)^2 >=0.

The real below-threshold integral permits nonnegative interchange of the
parameters; integrable complex majorants give the same analytic identity.

In a formal expansion of the canonical propagator about 1/(1-s), the
first correction is -f_R(s)/(1-s)^2, exactly this positive continuum.
It is not the exact inverse of a resummed one-loop functional. In fact,
for u>=2T, beta^3>=1/(2sqrt(2))>1/4 and w(u)>=C/(4u). The total continuum
weight diverges logarithmically when C>0 even though its propagator
integral converges. Do not infer a finite normalized exact spectral
measure, all-order reflection positivity or a complete Lorentzian spectrum.
