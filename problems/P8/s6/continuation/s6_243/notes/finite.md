# Actual scalar finite UV difference in the unchanged prescription

Let f2 be the complete spatial quadratic power coefficient and f4r(d) the logarithmic coefficient at independent source order r=0,1,2. The other lower radial grades and the full d-dependent r3,r4 logarithmic rows vanish exactly. Define

F2=f2(3)/(2pi²),
F4r=f4r(3)/(2pi²),
Sr=partial_d f4r(3)/(2pi²).

All are six-invariant expressions. Extract the tensor invariant coefficients BEFORE taking the dimension derivative, and keep detector/source ordering fixed. The generic time jets are specialized to the actual a,H only after extraction.

Use the SAME fixed comoving split r=m and radial MSbar normalization. At d=3-2epsilon the full sphere measure is
C(d)=2^(1-d)pi^(-d/2)/Gamma(d/2), C(3)=1/(2pi²).
For the displayed scale factor
C(3-2epsilon)/C(3) times exp(gamma_E epsilon)[mu²/(4pi)]^epsilon m^(-2epsilon),
its logarithmic derivative at epsilon0 is2-2log2-log(m²/mu²). This follows from the actual Gamma function at3/2 and is checked independently by differentiating the full expression. The physical prescription here fixes mu=1.

The two analytic radial tails are

integral_m^infinity r^(1-2epsilon)dr=-m^(2-2epsilon)/(2-2epsilon),
integral_m^infinity r^(-1-2epsilon)dr=m^(-2epsilon)/(2epsilon).

They are first defined in their respective convergent half-planes and analytically continued. Since f4r(3-2epsilon)=f4r(3)-2epsilon partial_d f4r(3)+O(epsilon²), the finite part of the entire continued tail is

-m² F2 Gamma0/2
+sum_r[(1-log2-ell/2)F4r-Sr]Gamma_r,
ell=log n.

Subtract the FULL covariant scalar pole counteraction
H_pole(d)/(64pi² epsilon)
before taking the physical limit. H_pole(3)=16 sum f4r(3)Gamma_r cancels the entire pole. Its d-dependence leaves the additional term

partial_d H_pole(3)/(32pi²).

The final finite expression is the sum of these terms, without any fitted finite target. An independent calculation multiplies the full analytic expression by epsilon, checks its value0 vanishes and takes its derivative at0; it gives the same formula. All eighteen invariant/source coefficients are serialized, not only a pole or tracefree restriction.

The geometric Hessian includes the full a^(d-2) and a^(d-4) volume factors and the dimension-dependent curvature identities. In particular a physical Euler boundary must not be set to zero before dimension differentiation. Both the mode slopes and evanescent counteraction are nonzero. Borrowing the old vector pole changes even the mass-curvature coefficient and fails the scalar match.

For the full coefficients f0,f1,f2, transposition swaps UD and UG. With the actual dt dP pairing and compact prepared variations, exact identities are

f2=transpose(f2),
f1-2dt transpose(f2)+transpose(f1)=0,
f0-transpose(f0)+dt transpose(f1)-dt² transpose(f2)=0.

These are local formal Green identities; they do not make a retarded nonlocal response symmetric. Every complete finite coefficient vanishes at P=0 because the result is a spatial difference, not because the full homogeneous response vanishes.

This is the local finite term in a direct endpoint/Laurent representation of the SAME scalar prescription. It is not an additional freely chosen action to add to S241's whole finite heat Hessian. A later assembly must relate the full subtracted kernels in the two representations and count the prescription exactly once.
