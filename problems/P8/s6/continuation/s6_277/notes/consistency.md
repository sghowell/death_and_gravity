# First and second preservation exclude smooth passage through this family

Set D_u=partial_u+L_H with the full normalized density flow of canonical.md.
The algebraic lapse constraint is C=H_N=0. Along any regular solution,

    0=dC/du=C_N N_u+K, K=D_u C.

At the fold C_N=0. Complete clock differentiation and the literal canonical
bracket give

    K=p[-3/25+3(P-rho)]+3P1/2-rho1.

If K is nonzero, no finite N_u preserves C. If K is zero, the denominator
3/25-3(P-rho) is strictly positive and fixes

    p_c=(3P1/2-rho1)/(3/25-3(P-rho)), |p_c|<21epsilon.

No new free gauge multiplier is inferred at this exceptional point.

At C=C_N=K=0, second preservation requires

    a lambda^2+b lambda+c=0, lambda=N_u,
    a=C_NN, b=2 D_u(C_N), c=D_u^2 C.

The normalized-weight detail matters. Off the constraint,

    K_N=D_u(C_N)+L_C C=D_u(C_N)-C C_p.

Only on C=0 does this become D_u(C_N). A generic independent polynomial
test checks the nonzero off-surface term. In forming D_u^2 C, the exact
C_uu, both mixed flow terms L_H(C_u), L_(H_u)(C), and
L_H(C_u+L_H C) are all retained. The heavy pair is evaluated at zero
only after these differentiations.

For clarity the ENTIRE coefficients at this homogeneous fold are

    a=(11391+1400P-1600rho)/400,
    b=-(600Pp-700P1-1200rho*p+800rho1-12627p)/200,
    c=-(20000P^2-36000P*rho-36000Pp^2-9100P-4500P1p-3000P2
         +16000rho^2+45000rho*p^2+4940rho+2000rho2
         +55485p^2-120528)/2000.

Every actual profile jet has magnitude at most epsilon=2/10^400.
At the sole compatible p use |p|<21epsilon. For each of the three
polynomials, bound every nonconstant monomial by its absolute coefficient
times the product of these radii. The exact finite sums give

    28<a<29, |b|<1, c>60,
    b^2-4ac<1-4*28*60=-6719.

The zero-profile diagnostic is
3(18985lambda^2+40176)/2000; it is not substituted for the actual bound.
Thus the actual second equation has no real lambda.

Only a C1 lapse is needed for this contradiction. If x_u=F(u,N,x)
with smooth F and C1 N, then x is C2. Taylor-expand C(u,N(u),x(u))
to second order at the fold. The vanishing coefficient C_N removes any
need for N_uu; the remaining quadratic coefficient is the one above.
The argument does not assume an unproved C2 extension of N.

The full raw N,T block gives the same obstruction. At the fold it is

    [[-9p^2,3p],[3p,-1]], null vector(1,3p).

Complete the temporal square around the exact Tstar. Projecting the full
secondary time-preservation equations onto this null vector cancels both
Tstar time and phase derivatives and gives precisely K. The negative
second condition therefore is not created by eliminating a regular
temporal variable.

C_s=3 and the nonzero Pbeta imply a regular constraint gradient. A rank
change of the auxiliary Poisson matrix is distinct from a redundant or
ineffectively squared constraint. A nonvanishing rescaling of C multiplies
the first condition, and at C=C_N=K=0 its second quadratic, by the same
nonzero factor. A regular lapse-coordinate change maps finite lambda
affinely with nonzero slope. Neither creates a real solution.

The role of repeated consistency conditions is described by
[J. David Brown, section VII](https://arxiv.org/html/2201.06558).
His illustrative model is different; no theorem about quantum resolution
or our endpoint is imported from it. Here all displayed conditions are
derived from the original parent itself.

This excludes C1 clock-time passage through THIS homogeneous zero-heavy,
zero-curvature fold family for every real trace density p. It does not
classify all inhomogeneous folds, choose a quantum constraint prescription,
or exclude the original corrected finite hybrid.
