# Entire elastic angular convolution

Let mu=m^2, Q=s-4mu, V0=s^2-4mu*s+2mu^2, with4mu<s<n.
For t=-Q(1-x)/2,u=-Q(1+x)/2 the exact tree is

 A(x)=P/(1-x^2)+B/(v^2-x^2)+a0+a2 P2(x),
 P=4V0/(kappa Q), v=1+2n/Q, B=4g^2 v/Q,
 a0=C+g^2/(n-s)+[-2s+6mu-2mu^2/s+Q^2/(6s)]/kappa,
 a2=-Q^2/(6kappa s).

The independently checked invariant form is
Agrav=-sum_channels[2mu^2-2mu*z-other1*other2]/(kappa*z),
Amatter=C+g^2 sum_channels1/(n-z). No interference is removed.

Temporarily replace the first denominator by w^2-x^2, w>1. This is an
auxiliary angular regulator, NOT a graviton mass or detector resolution.
For unit external axes a,b and z=a.b, define the normalized spherical
average and T_w(x)=1/(w^2-x^2). Feynman parametrization of two linear
resolvents and direct spherical integration of a squared resolvent yield

 I(w,v,z)=integral_0^1 dt /
 [(tw+(1-t)v)^2-|ta+(1-t)b|^2]
 =log[(c+Delta)/sqrt((w^2-1)(v^2-1))]/Delta,
 c=wv-z, Delta^2=c^2-(w^2-1)(v^2-1).
 J(w,v,z)=<T_w(a.n)T_v(b.n)>=[I(w,v,z)+I(w,v,-z)]/(2wv).

For w,v>1 the Feynman denominator is strictly positive because the
convex energy exceeds1 while the convex unit-vector norm is at most1.
The physical real sheet has c>0 and Delta>=0. The coincident limit
Delta0 is1/c, including w=v,z1. The certificate checks the full
quadratic primitive, its discriminant and the endpoint cross ratio;
isolated primitive singular coordinates are filled by this positive
integral's continuous limit. No chosen endpoint value substitutes for
the whole integral.

Azimuthal averaging gives <P2(b.n)>_az=P2(a.n)P2(z).
Consequently <T_w>=Q0(w)/w and
<T_w(a.n)P2(b.n)>=Q2(w)P2(z)/w, where
Q0=log((w+1)/(w-1))/2 and Q2=((3w^2-1)Q0-3w)/2.
Legendre orthogonality supplies <P2(a.n)P2(b.n)>=P2(z)/5.

Thus the ENTIRE regulated cut integrand average is
P^2 J(w,w,z)+2PB J(w,v,z)+B^2 J(v,v,z)
+2P[a0 Q0(w)/w+a2 Q2(w)P2(z)/w]
+2B[a0 Q0(v)/v+a2 Q2(v)P2(z)/v]+a0^2+a2^2 P2(z)/5.
Its physical two-Phi phase-space factor is beta/(32pi), with
beta=sqrt(1-4mu/s). Identical-pair1/2! and optical1/2 are both included.
Numerical quadratures and finite-parameter diagnostics are independent
checks, not the proof of these continuum identities.
