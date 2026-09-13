# Complete joint scalar/longitudinal-Proca reduction and measure

On a homogeneous light background the source one-form is closed. The
classically aligned lift Wbar=Sbar has F(Wbar)=0 and Wbar-Sbar=0. This is a
background slice, not a restriction of arbitrary vector perturbations.
Choose a real Fourier phase Wi=A_L sin(Px), W0 fluctuation=-A0 cos(Px).
Then the electric combination is Adot_L-P A0. Set

    r=R-1, rN=R_N, dH=Hhat-Hclock,
    K=zeta*Cchi/(N*a_hat^2), Yv=N*Cchi/a_hat^2.

To the entire base Lagrangian in adm.md add

    K(Adot_L-P A0)^2/2
      +Z[A0+r(3vdot-b)+3rN*dH*n]^2/2-Yv*A_L^2/2.

This is the complete aligned quadratic square, not its first-order cross
term alone. Its three-velocity Hessian is diag(-6D gamma,Z,K), where

    gamma=1-3Zr^2/(2D)=1-3(R-1)^2/(2R).

For1/2<R<2, gamma>1/4, as follows from the displayed exact factorization
in coupled.data. Work where D,Z,J,K are positive and gamma>0. This
ensures the stated velocity/auxiliary pivots, not a uniform all-momentum
perturbed physical stability theorem.

Perform the full three-velocity Legendre transformation. The auxiliary
Hamiltonian Hessian on(n,b,A0) has

    Dnb=DbA0=0, Dbb=-2D/3,
    DA0A0=-Z/gamma,
    DnA0=-3Z(r Theta/D+rN dH)/gamma.

Eliminating A0 in this matrix leaves lapse pivot-2J. Thus

    det Daux=-4DZJ/(3gamma),
    det Avelocity*det Daux=8D^2Z^2JK.

The cancellation of gamma belongs to the whole Gaussian measure. It does
not justify independently regularized factorization into unrelated sectors.

For normalized momenta(pv,ps,pA), define

    L=Theta*(pv-3r P pA+3c sigma)/D
        -w ps/Z-Lnv v-3rN dH P pA.

The constraints include n=L/(2J), b=pv/(2D); the full A0 solution is also
stored before those substitutions. The entire reduced Hamiltonian is

    h=ps^2/(2Z)-c(pv-3r P pA)sigma/(2D)-3c^2 sigma^2/(4D)
      +Y q sigma^2/2-(2Cq+Vvv)v^2-Vvs v sigma+L^2/(4J)
      +r P pA pv/(2D)-3r^2 P^2 pA^2/(4D)
      +pA^2/(2K)+P^2 pA^2/(2Z)+Yv A_L^2/2.

Every entry is checked against the direct Legendre transform and successive
auxiliary solves. At r=dH=0,D=Z=1 this equals the WHOLE S251 scalar
Hamiltonian plus ordinary longitudinal Proca. First mixed r vertices and
second source contacts are nonzero despite that reference block split.

For W=kappa*a_hat^3, restore physical momenta Pi=W*p. If Caux denotes
the three secondary constraints, their physical Poisson bracket is the
nonzero lower block Laux. The full six-constraint Dirac matrix and inverse
are

    M=[[0,-Dphysical],[Dphysical^T,Laux]],
    M^-1=[[Dphysical^-T Laux Dphysical^-1,Dphysical^-T],
          [-Dphysical^-1,0]], Dphysical=W Daux.

The code constructs Laux from the whole Hamiltonian rather than setting it
to zero. det M=(det Dphysical)^2. The density
4W^3DZJ/(3gamma) cancels the auxiliary delta-function Jacobian. Physical
variables have zero brackets with the primary auxiliary momenta; the zero
lower-right inverse block therefore leaves their complete canonical bracket
unchanged. Independent twelve-dimensional Poisson calculations check this
at the actual kappa normalization, including a zero-Theta fixture.

For x=(v,sigma,A_L,pv,ps,pA), the canonical map is

    Cphase=sqrt(kappa)*diag(1,1,sqrt(zeta),
                           a_hat^3,a_hat^3,a_hat^3/sqrt(zeta)),
    Cphase Omega Cphase^T=W Omega,
    Hcanonical=W Cphase^-T hessian(h) Cphase^-1.

Its full time derivative cancels the weighted3Hhat momentum damping. The
canonical symplectic generator and the inherited time boundary are retained.
Before field rescaling, the joint configuration density squared is
8W^6D^2Z^2JK. Rescaling the four scalar/auxiliary configuration variables
by sqrt(kappa), and the two vector variables by sqrt(kappa*zeta), gives
8a_hat^18D^2Z^2JK/zeta^2. At the reference this is8a^16J/zeta. The
reference outer reduced configuration determinant is

    2a^7J/[Theta^2(1+zeta P^2/a^2)].

That last chart requires Theta nonzero. The canonical phase and joint
Dirac system remain regular at its crossing. All Gaussian phases are
continued jointly from the chosen finite-regulator canonical evolution.
This is an already gauge-fixed physical quadratic measure, not a nonlinear
covariant, ghost, affine-complement or field-map determinant prescription.
