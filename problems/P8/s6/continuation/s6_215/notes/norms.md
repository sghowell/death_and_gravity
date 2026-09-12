# Explicit known-piece and conditional full-response estimates

V04[D] is the L2 spacetime Hilbert norm of (1-Delta)^2 acting on the lapse, shift and symmetric spatial tuple. U138[G]^2 sums through13 time derivatives the squared L2 norm of (1-Delta)^4 acting on that tuple. Matrix norms are Frobenius. These are unreduced coordinate norms, not canonical scalar norms.

Both time primitives have L2 norm<=1 on the unit slab, by Cauchy-Schwarz and integration. For j>=1, (If)^(j)=f^(j-1). On |z-t|<=1/4 with real |t|<=1/2, |z|<=3/4 and |1+z^2|>=7/16. Thus |H(z)|<7, |a(z)^-2|<28, and Cauchy gives H_j<=7j!4^j, c_j<=28j!4^j. Put B_r(b)=sum_j0..r binom(r,j)b j!4^j.

Every eta derivative through13 is individually bounded by U. The spatial H7 norm of every chi derivative through13 is at most[1+B12(28)]U, including j0 because B0<=B12. Each Qsyn time derivative in H6 is at most[3+4B13(7)+2B12(28)]U; use 2sqrt3<4 for its isotropic term. The Hilbert sum of14 derivatives costs sqrt14<4. The exact resulting coefficient62408287126207901212 is below1e20.

The detector j0 calculation in H1 gives1+4*7+2(1+28)=87<100; V04 safely controls the extra spatial derivatives. Hence

    M[Qsyn_D]<=100 V04[D],
    Z136[Qsyn_G]<=1e20 U138[G].

Orthogonal trace/tracefree/scalar projections are contractions. The initial-only domain is justified in boundary.md.

For the one-point Ward terms use a<2, |H|<=2, |H'|<=4, c<=1 and the original rho/P derivative bounds1e30 through order1. Then ||E||F<4e30, ||E'||F<2e31. In the appropriate input norm eta<=1 and chi<=2, so xi<=3 and its full spacetime derivative matrix<=4. The ADM gauge tuple has norm<13 because n,beta equal the original components and ||Qxi||<=12. The same estimates hold with the spatial derivatives used below.

The first metric chart has norm<=6 times its ADM tuple. The second has pointwise bilinear norm<=30 times the two tuple norms: separate component coefficients sum to2+8+6+6+4=26<30. The source Lie-density term is bounded by(20+32+32)e30=84e30. Its complete source Ward coefficient is[6*84+30*4*13]e30=2064e30.

For Gsyn, Qsyn and its spatial gradient have bounds13U and Qsyn' has bound29U. Thus hG and its spatial gradient are<=52U and hG'<=324U; use376U for the spacetime derivative. Spacetime Cauchy gives detector Lie coefficient4(3*376+2*4*52)e30 and chart coefficient30*4*13^2 e30. Summing both ordered Ward terms gives28520e30 VU<1e35 VU.

The matched tracefree current has coefficient2e95 M Z136. After projection its contribution plus the full one-point correction is below1e118 VU. This is only a known-piece bound.

If, in the same prescription and initial-only domain, each missing scalar form obeys |Rab(d,g)|<=Lab M[d]Z136[g], the conditional full bound is

    |Rfull(D,G)| <= [1e22(2e95+Ltt+Lt0+L0t)+1e35] V04[D] U138[G].

No finite Lab is asserted. Derivative loss remains; this is not a reduced inverse, contraction, finite-amplitude or quantum stability theorem.
