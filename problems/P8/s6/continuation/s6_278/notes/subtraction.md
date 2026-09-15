# Exact low-cut subtraction and the transfer-channel coefficient

Use the logarithms on their continuous principal sheets. The identity

    log(-z/nu^2)
      =[log(-z/nu^2)-log((L-z)/nu^2)]+log((L-z)/nu^2)

is an identity of the displayed logarithms; it does not use an unchecked
logarithm-of-a-product simplification across a cut. For0<z<L, only the
first logarithm is discontinuous. For z>L, both have the same imaginary
jump and their difference has none. Thus the written low piece carries
exactly the original discontinuity below L.

Multiply by EACH full channel polynomial and sum s,t,u. Keep the real
local polynomial unchanged. This specifies the finite low-cut convention
and its exact add-back, rather than erasing a cut and pretending its
analytic remainder was the original amplitude. Polynomial ambiguities of
a dispersive subtraction must be matched to this convention explicitly.

For L>4m^2, the resulting M1-sector logarithms have all their branch
points above the massive crossing center. With D=L-2m^2, the s and u
polynomials at t0 are14m^4 plus or minus6m^2 v plus v^2.
Direct differentiation gives

    b20=-c[2log(D/nu^2)+log(L/nu^2)-12m^2/D-14m^4/D^2].

Here b20 refers ONLY to the displayed nonlocal sector after subtraction;
the unchanged local b20 must still be added.
Differentiating with respect to L gives the full formula in FORMULATION.
Equivalently,

    d_L b20=-2 Im A_s(L,0)/[pi(L-2m^2)^3]-c/L.

The last term is from the crossed transfer logarithm. The derivative is
strictly negative for positive m^2 and L>4m^2; this is a fact about this
cut convention, not an RG flow or a positivity verdict.

For the named m^2=nu^2=1,L=10^196 use0<log(L-2)<log L<588,
since log10<3, and12/(L-2)+14/(L-2)^2<1. The entire bracket is below2000
in absolute value, and pi^2>1 proves the recorded rational bound.
It is far below4lambda, but it bounds neither independent local terms
nor the other sectors or all higher orders. The cap is a scale in a
specified subtraction, not evidence for a physical cutoff.
