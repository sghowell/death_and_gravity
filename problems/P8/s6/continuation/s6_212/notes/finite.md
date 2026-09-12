# Original radial normalization, finite subtraction and explicit norm

Let d=3-2epsilon. For the spatial difference, the high-band radial UV expression is

    M(d) Z(epsilon)[f2(d)m^(2-2epsilon)/(2epsilon-2)
                       +f4(d)m^(-2epsilon)/(2epsilon)],

where the source-time slots are implicit,
M(d)=2pi^(d/2)/[Gamma(d/2)(2pi)^d], and

    Z(epsilon)=[exp(EulerGamma)mu^2/(4pi)]^epsilon.

This is the unchanged MSbar normalization. To check it independently, integrate a radial monomial omega^(1-2n)z^k/4 using the beta integral. After extracting m^(4-2n)/(64pi^2), the dimension-dependent factor is

    2sqrt(pi)(d/2)_k/Gamma(n+k-1/2),

times Gamma(epsilon+n-2)exp[epsilon(EulerGamma-ell)]. It is exactly the previously fixed radial implementation for every monomial, not a new scale choice.

At epsilon0 the power term is-m^2 F2/2. The normalized-sphere derivative contributes-partial_d f4(3)/(2pi^2). Combining the angular factor, Z and lower-tail mass factor contributes(1-log2-ell/2)F4. In particular the split is the original fixed comovingm, not a*m.

Subtract the original pole Hpole(d)/(64pi^2 epsilon) with fixed scalar weights and continued metric contractions/volume. Its finite contribution to the subtracted expression is+partial_d Hpole(3)/(32pi^2). This proves the finite coefficient displayed in FORMULATION.md. The original pair-band conversion remains separate; no finite conversion artifact is added after S209's complete cancellation.

The continued volume contributes log(a) times the physical pole. Together with the other evanescent terms it makes the first source-time coefficient equal the proper-time derivative of the second one. A missing volume logarithm produces a nonzero Green defect. The evanescent odd endpoint and dimension-dependent scalar probe also give explicit nonzero omission controls.

## Norm

Use Z24^2=sum_(r=0)^2 integral dt d^3P/(2pi)^3 (1+p^2)^4 |partial_t^r Gamma_hat|_F^2.

At ell0,m1000 multiply each invariant coefficient by pi^2*a and replace L=log(1+t^2), B=log2. Every result is an exact rational polynomial in t,p,L,B, with p degree at most4. On the original slab use|t|<=1/2,0<=L<=1/4,0<B<1,a>=1 and pi^2>9. Sum absolute monomial coefficients with these weights.

The exact nine-row bound is

    2905633954513/66355200 < 1e5.

Since|T|,|V|,|W|<=||D||F||Gamma||F, Plancherel and spacetime Cauchy-Schwarz give the stated finite coefficient norm. The symbol is polynomial in P at the origin; both canonical metric insertions multiply the bound by4/kappa, yielding4e-795.

This bound concerns only the local UV finite difference. Combining it with the homogeneous anchor and finite actual remainders, and proving the complete dimension-limit interchange, is still required. It is not a completed full response, mixed inverse or nonlinear estimate.
