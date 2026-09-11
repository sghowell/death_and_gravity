# S6.196 formulation

## Claim

Retain the actual canonical Gaussian Proca sector, mass1000,
kappa1e800, the fixed profile, all-order prepared pure Hadamard
state and unchanged finite stress prescription. The CD clock is
a=(1+t^2)^2 on I=[-1/2,1/2]. Real smooth compact spatial tracefree
source Gamma and readout D have supports inside the open slab,
with sup time(supp Gamma)<inf time(supp D). They need not be
homogeneous and have no spatial Fourier cutoff.

Let T[f]=integral dt dx a^3 f_ij T_ij in the physical frame,
psi_f=(T[f]-<T[f]>)Omega, and P_K^(2) retain both created
momenta |k|,|l|<=K. With the unchanged S186 norm

    N[f]^2 = sum(j=0..3) ||partial_t^j f||_L2(dt dx;F)^2
             + ||grad_x f||_L2(dt dx;F)^2,

the discarded vector obeys ||(1-P_K^(2))psi_f||^2
<1e52 N[f]^2/K for K>=1000 and nonzero f.

The complete weak first-order current response on these supports is

    B(D,Gamma) = i<[T[D],T[Gamma]]>/4.

All metric and fixed local finite-curvature contacts vanish by
support, without changing the prescription. Its magnitude is
<5e49 N[D]N[Gamma]. A common orthogonal projection gives
|B-B_K|<5e51 N[D]N[Gamma]/K, so all internal and external
momentum regulators are removed. If either test is zero, each
corresponding quantity is exactly zero, not a strict inequality.

For canonical tensor h=sqrt(kappa)gamma/2 the two chain factors
give |B_can|<2e-750 N[D]N[Gamma] and regulator error
<2e-748 N[D]N[Gamma]/K. At computational K1e16 the latter is
<2e-764 N[D]N[Gamma]. Mutually spacelike supports have zero
full response by complete Proca locality.

## Boundary

This is a weak bilinear first-order limit on separated supports.
It is not operator-norm differentiability of the full metric
evolution, a finite-amplitude quantum remainder, or a response
on overlapping times or the diagonal. In particular it does not
supply the same-space estimate needed for a feedback inverse.
K is not a physical cutoff. Full mixed constraints, interacting
background, stability, parent matching, V/G/B and original P8
remain OPEN.
