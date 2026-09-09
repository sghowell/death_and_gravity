# Positive infrared-completed state of the exact global equation

At u=0 and k>=K use the Borel-summed graph R=A-i B, B>0.
Its canonical ordered Cauchy covariance is

    C_R=hbar/(2*kappa) *
        [B^-1, B^-1*A+i*I;
         A*B^-1-i*I, A*B^-1*A+B].

Choose the positive-diagonal Cholesky factor B=L L^T and T=L^-T.
Then C_R is the Gram product of

    sqrt(hbar/(2*kappa)) * [T; (A-i*B)*T].

This proves positivity without requiring A and B to commute.
Direct matrix identities give C_R-C_R^T=i*hbar*Omega/kappa
and Hermiticity. Its range is the negative-frequency graph.
The matrix square roots and inverses have smooth classical
order-zero symbols because B has a positive leading matrix and
a uniform positive high-frequency lower bound.

Map to the original density phase with sqrt(k)*E^-1 at u=0.
This map is symplectic, so the density covariance retains exactly
the same CCR. Its high-frequency symbol order is at most three.
No unnormalized use of E is allowed.

For low k use the fresh INITIAL positive prescription
C_mu=hbar/(2*kappa)*[mu^-1 I,iI;-iI,mu I], mu=sqrt(1+k^2).
Choose a smooth radial cutoff 0<=chi<=1 which is zero for k<=K
and one for k>=2K, and set

    C_initial=chi*C_high+(1-chi)*C_mu.

Only the already well-defined region k>=K is used for C_high.
Positivity and the CCR are preserved by convexity. A cutoff built
from squares of a smooth partition also supplies smooth Gram
factors if desired. The difference from the high-frequency graph
construction has compact frequency support, hence a smooth
spacetime kernel after exact evolution. There is no singular
choice of |xi| at xi=0 in this completed covariance.

Define the new two-point function using the SAME complete global
density equation on both legs:

    W(t,s,k)=U(t,0,k) C_initial(k) U(s,0,k)^T.

It is a bisolution, is of positive type, and has the exact commutator
i*hbar U(t,s,k) Omega/kappa. On every compact time strip the
transfer has polynomial degree nine, so the two-point multiplier
has a degree-21 majorant. All differentiated multipliers are also
polynomially bounded. It defines a positive spacetime distribution
on compact-time spatial-Schwartz tests and in particular on compact
spacetime tests. The CCR uncertainty condition is the same Gram
positivity just proved. Wick pairing of this covariance constructs
a quasifree state on the corresponding bosonic field algebra.

This is NOT the previously chosen global mu state promoted by a
name change: only its low-frequency initial prescription is used,
and the high-frequency initial covariance is replaced by the
all-orders graph construction. Different admissible Borel choices
are allowed; the theorem is existence of such a state, not a
numerically unique vacuum or a preferred Fock representation.

The complete response and therefore quadratic relational
microcausality are unchanged. Positivity, CCR and exact evolution
alone do not prove the short-distance condition; the separate
wavefront argument supplies it.
