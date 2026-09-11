# Differentiating the complete infinite weighted tail

Use W(epsilon,t,k)=omega(epsilon,t,k)Delta(epsilon,t,k).
The weight also varies with epsilon. The full product
rule, with |omega1|<=3nu_minus and |omega2|<=6nu_minus,
gives

    ||W||<=B0 nu_minus^-9,
    ||W1||<=B1 nu_minus^-8,
    ||W2||<=B2 nu_minus^-7,

where
B0=3C0,
B1=3(C1+C0/K),
B2=3(C2+2C1/K+2C0/K^2),
and C_a are the covariance displays in covariance.md.

The proof band nu_minus>=K=1e16 is independent of
epsilon and retains the same lower frequency
nu_minus^2=m^2+b|k|^2, b=(99/100)/(25/16)^2.
It is not an epsilon-dependent physical cutoff.
For d^3k/(2pi)^3, the exact radial Jacobian and
b^-3/2<4, pi^2>9 give

    integral ||W_a|| d^3k/(2pi)^3
      <[4 B_a/(18(6-a))] K^(-6+a),  a=0,1,2.

Every upper radial endpoint is infinity. The resulting
bounds are respectively below1e-65,1e-47,1e-28.
No fixed numerical momentum grid replaces the tail.

For each finite k the exact Hamiltonian initial-value
problem, positive-root frame and finite reference
are smooth in epsilon. The same holds jointly in
(t,epsilon) on compact parameter and time intervals.
The displayed common integrable majorants justify
differentiating the high-band Bochner integral twice.
They also give continuity of those derivatives,
uniformly in t: finite-k joint continuity is uniform
on compact sets, and dominated convergence applies
to the uniform-in-time norms.

Thus the matrix-valued high-band tail integral T
is C2 in epsilon as a C0(I)-valued function. Taylor's
formula with its integral remainder proves

    ||T(epsilon)-T(0)-epsilon T'(0)||C0(I)
      <=epsilon^2*1e-28/2,  |epsilon|<=1/100.

The inequality is strict for nonzero epsilon; at zero
the remainder and its bound are both exactly zero.

The entire statement concerns omega times the
actual-minus-finite-reference balanced covariance.
The reference itself has UV-divergent local terms.
Correctly subtracting those terms in the original
covariant prescription, including all contacts,
is still required. The low band and full physical
stress vertices also remain separate obligations.
This finite-amplitude tail bound is not a full
gravitational feedback or quantum-background result.
