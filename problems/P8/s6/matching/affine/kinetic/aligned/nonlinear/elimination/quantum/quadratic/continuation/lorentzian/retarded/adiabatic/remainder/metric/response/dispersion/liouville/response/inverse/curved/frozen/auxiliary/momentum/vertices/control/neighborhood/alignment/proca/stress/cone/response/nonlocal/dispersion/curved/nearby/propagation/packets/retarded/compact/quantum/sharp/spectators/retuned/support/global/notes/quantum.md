# New global reduced scalar state and commutator

The original density evolution is symplectic and, on every finite time
strip, has real-momentum polynomial degree at most nine. Its original
generator is polynomial in spatial momentum of degree four.
Differentiated Duhamel induction gives a sufficient degree
9+13*|alpha| for every spatial Fourier derivative: the new derivative
term adds at most nine degrees from the outer transfer and four from
a generator derivative to the lower-order bound. Compact momenta are
smooth by the original entire generator. Reversing time gives the same
finite bounds with absolute duration. Thus global evolution preserves
Schwartz tests at all finite times, with continuous compact-time bounds.

At u=0 choose mu(k)=sqrt(1+|k|^2) and the ordered covariance

    C0=hbar/(2*kappa)*[mu^-1 I,i I;-i I,mu I]=B*B^dagger,
    B=sqrt(hbar/(2*kappa))*[mu^-1/2 I;-i*mu^1/2 I].

For positive hbar,kappa this is positive and tempered, with
C0-C0^T=i*hbar*Omega/kappa. Its Wick extension defines a reduced
quadratic Gaussian state. Transport using the GLOBAL retuned U_rho(t,0)
at both endpoints; neither nearby state nor old fast-clock evolution is
transferred. Gram factorization, Schwartz preservation and exact
symplecticity prove the resulting global state properties.

Its state-independent relational commutator multiplier is
-i*hbar*U_rho(t,s)[1,3]/kappa. The normalized physical Kubo factor
i*kappa*a(s)^3/hbar gives exactly G_hat=a(s)^3*U_rho(t,s)[1,3].
The unscaled source j=kappa*J retains an extra 1/kappa in its response.
The complete global support theorem therefore implies
[O(f),O(a^3 J)]=0 for compact matter-spacelike supported tests at any
finite source/output times. Smooth positive volume multiplication does
not change test support.

This construction supplies neither Hadamard wavefront behavior nor
renormalized stress, a unitary implementation on one global Fock space,
interacting microcausality or a semiclassical Einstein solution.
