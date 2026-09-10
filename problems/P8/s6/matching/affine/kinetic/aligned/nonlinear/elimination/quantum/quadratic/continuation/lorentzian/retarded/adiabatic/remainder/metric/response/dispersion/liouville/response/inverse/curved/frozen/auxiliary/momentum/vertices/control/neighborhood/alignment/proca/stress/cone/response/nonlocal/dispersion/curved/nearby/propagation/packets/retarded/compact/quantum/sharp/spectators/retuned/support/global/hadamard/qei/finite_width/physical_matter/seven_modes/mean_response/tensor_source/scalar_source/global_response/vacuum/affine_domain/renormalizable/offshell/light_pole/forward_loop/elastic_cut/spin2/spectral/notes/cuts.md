# Both parameter cuts and independent spin-two normalization

Let mA denote the mass squared of the two stress-adjacent lines and
mB that of the other line. With chi the other-line parameter and z
the unit pair parameter, the denominator is

D=(1-chi)mA+chi mB-chi(1-chi)
  -(1-chi)^2 z(1-z)(T+i0).

The positive prescription gives Im(1/D)=pi delta(Re D).
The two pair roots each have Jacobian
(1-chi)^2 T beta_chi. Integrating them gives

rho_A(T)=g/(8pi sqrt(T)) integral_0^chi_minus
 chi^2 dchi/sqrt[(T-4)chi^2-2(T-2mA+2mB-2)chi+T-4mA].

Here rho means Im F, not Im F divided by pi. For T>4mA the
quadratic is positive at chi=0 and equals -4mB at chi=1, so its
first root lies inside the integration interval. Use only the
positive-square-root branch below that first root.

Put A=T-4, C=T-4mA,
Z=(T-2mA+2mB-2)/sqrt(A C), chi=sqrt(C/A)y. For both actual
mass assignments, Z>1 above the respective threshold. The universal
integral is

integral_0^(Z-sqrt(Z^2-1)) y^2 dy/sqrt(y^2-2Zy+1)
 =Q2(Z)
 =(3Z^2-1)log[(Z+1)/(Z-1)]/4-3Z/2.

To verify it, set w=Z-y and use the anchored primitive

(3Z^2-1)log[(w+sqrt(w^2-Z^2+1))/sqrt(Z^2-1)]/2
 +(w/2-2Z)sqrt(w^2-Z^2+1).

Its derivative and both endpoints are checked. At the upper
endpoint w=Z the radical is one. All logarithm arguments are
positive; squaring (Z+1)/sqrt(Z^2-1) gives (Z+1)/(Z-1), which
establishes the required real-log identity. The angular derivation
also checks that identity by derivative and an anchor at Z=2.

The two actual densities are therefore

rho_L(T)=g Q2(1+2M/(T-4))/(8pi sqrt(T(T-4))),       T>4,

rho_H(T)=g(T-4M) Q2((T-2M)/sqrt((T-4)(T-4M)))
         /[8pi sqrt(T)(T-4)^(3/2)],               T>4M.

Each is zero below its own threshold and continuously zero at it.
For the heavy argument, the squared numerator minus squared
denominator is 4(T+M^2-4M)>0 above T=4M.

## Independent tree-vertex and phase-space derivation

Let pe,pa be external and intermediate pair momenta in the
transfer center-of-mass frame. The canonical spin-two stress has
the relative momentum factor pa^2/pe^2. With S=1+iT, the
two-particle phase measure per cosine is beta/(16pi).
The identical-pair 1/2! and the factor 2 Im are separate:

Im F=beta/(64pi)*(pa^2/pe^2) integral_-1^1 A_tree(z) P2(z) dz,
P2(z)=(3z^2-1)/2.

For the light pair the full actual tree amplitude is
C0+g/(B-kz)+g/(B+kz), with B=M+(T-4)/2 and k=(T-4)/2.
C0 includes every angle-independent contact and direct-channel
exchange; its spin-two projection is zero.

For the heavy pair the actual conversion amplitude is
g/(B_H-K_H z)+g/(B_H+K_H z), with
B_H=T/2-M and K_H=sqrt((T-4)(T-4M))/2.
There is no independent Phi^2 H^2 contact in this polynomial model.

The anchored Legendre integral gives 4g Q2(B/k)/k for either
even exchange sum. Substitution in the phase-space expression
reproduces the triangle prefactors exactly. The heavy channel here
uses free intermediate heavy lines at the retained perturbative
order, not exact stable asymptotic particles.

## Stable positive evaluation

Expanding 1/(Z-z) uniformly for Z>1 and integrating its even
Legendre moments gives

Q2(Z)=sum_(n>=1) [2n/((2n+1)(2n+3))] Z^(-2n-1).

Every coefficient is positive and at most 2/15. After N terms,
the remainder is strictly positive and at most

2/[15 Z^(2N+3)(1-Z^-2)].

The exact rational enclosure implements this expression. The
software accepts N from one through 32; that is an implementation
limit, not a physical cutoff. Direct evaluation of the difference
of large logarithmic terms is not used for the actual hierarchy.

The leading term 2/(15Z^3) and its explicit tail show the
five-halves threshold power of both densities at the actual M>1.
No massless cut or exact interacting threshold is inferred.
