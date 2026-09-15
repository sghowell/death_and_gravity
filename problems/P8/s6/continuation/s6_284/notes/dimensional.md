# Whole D-dimensional angular integrals and common scalar masters

Let h=s/q,q=s-4mu,z=1+2t/q and n=D-2. Put
rL^2=1-x^2,rR^2=1-y^2,rL.rR=z-xy,
deltaL=h-x^2,deltaR=h-y^2. Half the TT trace numerator in the
polarization note is N_D. The entire graviton cut is the normalized
D-dimensional two-body phase times
q^2/(256pi kappa^2)*E[N_D/((h-x^2)(h-y^2))].

The unit sphere has spatial dimension d=n+1. Rotational invariance gives
E[x^2]=1/d,E[xy]=z/d,E[x^2y^2]=(1+2z^2)/(d(d+2)).
Conditioning on x gives E[y^2|x]=z^2x^2+(1-z^2)(1-x^2)/n.
For L0=E[1/(h-x^2)], this yields
E[xy/(h-x^2)]=z(hL0-1) and
E[y^2/(h-x^2)]=[(n+1)z^2-1](hL0-1)/n+(1-z^2)L0/n.

Exact polynomial division in BOTH axes expresses the whole integral as
C0 J0+C1 J1+CL L0+CP. The implementation constructs all four rational
coefficients, checks the unaveraged division and the complete n2 limit.
No uncomputed moment of higher degree remains.

For I(z)=E[1/((sqrt(h)-x)(sqrt(h)-y))], partial fractions and
axis-sign invariance give
J0=(I(z)+I(-z))/(2h), J1=(I(z)-I(-z))/2 in every dimension.
The normalized scalar cuts carry the same full D-dimensional phase:
Im I4(s,t)=4pi*phase/(s*q)*I(z),
Im C00mu(s)=-2pi*phase/q*L0, Im B00=pi*phase.
Thus inside overall1/(16pi^2 kappa^2),

 box(s,t)=q^4 C0/128+s*q^3 C1/128,
 C00mu=-q^3 CL/32, B00=q^2 CP/16.

The box is exactly [V(t)+2mu^2 EP/(1+EP)]^2 at n=2+2EP.
The same graph appears in the massive cut; it is not counted again.

For the complete massive elastic polynomial set
VD=s^2-4mu*s+2mu^2+2mu^2 EP/(1+EP),P=4VD/q,
b=-q^2/(4s),
a=-2s+6mu-2mu^2/s+q^2/(4s)-2mu^2 EP/[s(1+EP)].
The tree is P/(1-x^2)+a+b*x^2. Conditioning and the same sphere
moments give C0mumu=-2VD(a+bz^2) and

 Bmm=1/2{2Pb[-z^2+(1-z^2)/n]+a^2+2ab/d
                              +b^2(1+2z^2)/(d(d+2))}.

This retains the whole polynomial and single-resolvent interference,
not just its soft pole. Both species agree with S283 at EP0.
The exact common dimensional phase cancels when extracting coefficients;
it is not discarded from the amplitude. Its raw Gamma/4pi conversion
remains the one fixed in S283.
