# Full spatial tensor numerator

The full covariant transverse tensor at timelike P is
theta_ab=eta_ab-P_a P_b/s. With eta=diag(+1,-1,-1,-1) and fixed
spatial p, its spatial block is -I-p p^T/s. Therefore the complete
spatial cut from S6.199 is

    K(s,p)=rho2(s)Q2(s,p)/s^2+rho0(s)Q0(s,p)/s^2,

where Q0 and Q2 are the polynomial tensors in FORMULATION.md.
COM projectors alone are insufficient at nonzero spatial transfer.

The flat compact C^2 Hessian restricted to spatial directions is
Q2(z,p); the R^2 Hessian is 6Q0(z,p). These statements follow by
inserting an arbitrary spatial h into the four-dimensional linearized
Riemann tensor and contracting all Lorentz indices. Independent tests
do that calculation for noncollinear nonzero p and symbolic frequency.

Define N_i=z Q_i. It is cubic in z and a finite polynomial in p.
Thus F=sum N_i integral rho_i/(s^3(s-z)) has the required full
absorptive cut and no artificial z=0 projector singularity.
For each finite z off the massive cut the integral is absolutely
convergent, since rho=O(s^2). Multiplying out N before discussing
regularity is essential.

The finite operator acts on all symmetric spatial tensors. The
9-by-9 representation also includes the annihilated antisymmetric
subspace; it does not add physical polarizations. At nonzero p the
restricted spatial projectors are not Euclidean orthogonal COM
projectors. The norm proof uses their complete congruence and trace
terms, not a unit-projector assumption.

The physical local tensor/contact polynomial remains separate. The
specified F is a covariant nonlocal representative, not an assertion
that the full physical finite action has zero local coefficients.
