# Leading vacuum germ: complete inverse below its pole

At u=0 and vanishing gradient, r=(R-1)/X^2=-n. Substituting the canonical
u=Phi/sqrt(kappa) into the COMPLETE regular source and taking its lowest
field degree gives

    S4_mu=(n/kappa^2) partial_mu Phi [Z_Phi-X_Phi Box Phi].

This has four fields and five derivatives. The source is even under
Phi -> -Phi. Its higher terms have even field degrees six and above.
The exact switch, exponential and rational functions are not a finite
polynomial, so the complete source is not assigned a finite Fourier
support radius by this argument.

For S4 alone, convolution gives support radius at most4 on the stated
Fourier-unit field class. Products of the jet supremum bounds with the
outer gradient L2 norm give each component at most32n kappa^-2 U and
the positive four-component norm at most64n kappa^-2 U.

For |p_E|<=4, |p_M^2|<=16 and the Euclidean matrix norm of
p p^T eta is |p_E|^2<=16. The EXACT, unexpanded inverse from operator.md
then obeys, for 0<zeta<=1/2000,

    ||O^-1|| <=(1+16zeta)/(1-16zeta)<=63/62,
    ||I-O^-1|| <=32zeta/(1-16zeta)<=1/62.

The Lorentz pairing has the same absolute Cauchy-Schwarz bound as the
positive component norm because eta has norm one. Therefore the leading
degree-eight source functional satisfies

    |Delta S_L[S4]| <= kappa/124 (64n kappa^-2 U)^2
       =1024 n^2/(31 kappa^3) U^2 <10^-2392 U^2.

There is no mass pole on this finite support, so the causal/in-out
boundary prescriptions agree there. This is not a derivative expansion
of the propagator. It remains ONLY the degree-eight germ functional,
not an estimate of the remainder from all higher source field degrees.

## Why a generic global real-time L2 inverse bound is unavailable

Set m^2=1/zeta. In a small Fourier neighborhood of the positive-energy
mass surface choose a smooth transverse vector packet, for example
v=(p1/p0,1,0,0) times a compact scalar bump. Then p^T eta v=0.
Restrict the bump to p_M^2-m^2 in [epsilon,2epsilon] and p0 bounded
away from zero. Reflect it at negative momentum with conjugate values
to obtain a real Schwartz spacetime source. On these positive-measure
bands the inverse acts as multiplication by 1/(1-zeta p_M^2);
its norm is at least1/(2zeta epsilon). Normalize the packet in L2.
Letting epsilon decrease proves that the inverse has no finite global
unweighted L2 operator norm over arbitrary sources. No delta function
on the pole is needed for this control.

This is NOT a proof that the actual nonlinear image S(Phi) contains
these source packets. It is neither a no-go theorem for this parent nor
an exclusion of any V/G/B row. It explains why the Euclidean contraction
bound cannot simply be imported to the full real-time functional.
A bound adapted to the actual image, a finite/weighted time domain,
a state prescription and/or contour control remains to be developed.
