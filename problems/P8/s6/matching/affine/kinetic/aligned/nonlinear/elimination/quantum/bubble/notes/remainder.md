# Finite subtracted loop: analytic domain and quantitative bound

This is a quadratic vector mass-insertion loop on a stationary,
constant reference metric, at the isotropic on-clock mass. The external
coefficient h is frozen at any point of |u|<=1/2. Neither a global
Minkowski solution nor a curved/in-in transfer is assumed. Notation,
normalization and the full numerator are in [the pole note](pole.md).

## Finite logarithmic coefficient from the actual integral

Divide the d-dimensional radial integral by 1/(16*pi^2), use
epsilon=(4-d)/2 and the modified minimal-subtraction scale, and set
ell=log(Delta/mu_bar^2). For n=0,1,2 it is

    J_n=Delta^n exp[epsilon*(EulerGamma-ell)]
         Gamma(n+2-epsilon)/Gamma(2-epsilon)*Gamma(epsilon-n).

The code computes the residue and finite Laurent coefficient of this
expression, not just the Gamma residue. They are, respectively,

    J_0: 1/epsilon - ell,
    J_1: -2Delta/epsilon + Delta*(2ell-1),
    J_2: 3Delta^2/epsilon + Delta^2*(-3ell+2).

Thus each logarithmic coefficient is minus its pole. Corrections
from d-dependent angular factors and the dimensional continuation
of the mass trace are analytic at epsilon=0. They multiply the simple
poles into finite local polynomials of momentum degree at most four.
They cannot change the nonlocal logarithmic coefficient. Likewise
the remaining finite pieces displayed above are such polynomials.

The bubble's finite nonlocal part is consequently

    +1/(64*pi^2) integral_0^1 dx F(x,k) log(Delta/mu_bar^2),

up to those local polynomials, with the Gaussian -1/4 and minus-log
sign both retained. No assertion that the full unsubtracted finite
polynomial vanishes is needed. Taylor subtraction through external
momentum degree four removes that polynomial and the scale-dependent
polynomial F(x,k)*log(m2/mu_bar^2).

## Exact subtraction

With t=x(1-x), the homogeneous parameter pieces are

    f0=m2^2*(T^2+2U)/8,
    f2=m2*[C/2+t*(K*(T^2-2U)/4+T*B)],
    f4=t^2*[K^2*(T^2+2U)/8+2K*C+K*T*B+B^2]-t*K*C/2.

Their sum is F(x,k), checked before integration. Under k -> sqrt(z)k,
let v=t*K/m2. Subtracting the z^0,z^1,z^2 Taylor polynomial from
(f0+z*f2+z^2*f4)*log(1+z*v) leaves exactly

    f0*[log(1+z*v)-z*v+(z*v)^2/2]
      +z*f2*[log(1+z*v)-z*v]+z^2*f4*log(1+z*v).

At z=1, integrating this expression and dividing by 64*pi^2 defines
R_loop(k). The formula is independent of finite local choices of
momentum degree at most four. Higher local UV matching operators,
which are not calculated here, are not removed by this subtraction.

## Independent all-order isotropic finite control

The denominator reduction for A=cI leaves a scalar bubble with
coefficient c^2*(3m2^2+m2*K+K^2/4), plus local tadpoles/scaleless
pieces. Its finite logarithmic contribution therefore uses that
coefficient times integral log(1+t*K/m2). The angular derivation
instead gives c^2*[3m2^2+m2*K*(1/2+6t)+K^2*(10t^2-t/2)]
inside the logarithmic integral. Their subtracted functions agree.

To check all orders, use B_j=integral_0^1 t^j dx=j!^2/(2j+1)!.
For every integer n>=3, divide the Taylor coefficient difference
by the nonzero common factor c^2*(-1)^(n+1)*m2^(2-n)*B_(n-2).
Then B_(n-1)/B_(n-2)=(n-1)/[2(2n-1)] and
B_n/B_(n-2)=n(n-1)/[4(2n-1)(2n+1)] reduce the difference to
an exactly vanishing rational function of n. This is checked
symbolically, with independent direct x integrations at n=3,4,5,8.
Analyticity inside the massive logarithm's convergence domain then
identifies the functions, not just these sampled coefficients.

## Continuous complex-domain envelope

For complex Euclidean momentum, distinguish the bilinear K=sum k_mu^2
from kappa=sum |k_mu|^2. Use the principal logarithm analytically
continued from k=0. For the actual positive diagonal A, alpha>=beta,

    |K|<=kappa, |B|<=alpha*kappa, |C|<=alpha^2*kappa.

Let rho=kappa/m2<=1. Since 0<=t<=1/4, |v|<=t*rho<=1/4,
so this ball is strictly inside the first massive branch singularity.
The convergent log series and absolute geometric tail imply

    |log(1+v)| <= |v|/(1-|v|),
    |log(1+v)-v| <= |v|^2/[2(1-|v|)],
    |log(1+v)-v+v^2/2| <= |v|^3/[3(1-|v|)].

For the actual alpha=4/(9h), beta=28/(81h), direct triangle bounds give

    |f0| <= m2^2/h^2 * 904/2187,
    |f2| <= m2*kappa/h^2 * [8/81+(2032/2187)*t],
    |f4| <= kappa^2/h^2 * [(3640/2187)*t^2+(8/81)*t].

Here (T^2-2U) is positive in the actual direction. For example its
normalized quarter is 592/2187, and h^2*T*alpha=1440/2187;
their sum yields the f2 envelope. The f4 t^2 envelope is the sum
(904+864+1440+432)/2187, and its remaining absolute term is 8t/81.
These inequalities do not require real k or a sign for K*C.

Taking the three log tails above, replacing their denominators by
1-rho/4, and integrating yields the coefficient

    C_abs=integral_0^1 [ (904/2187)*t^3/3
            +(8/81+(2032/2187)*t)*t^2/2
            +((3640/2187)*t^2+(8/81)*t)*t ] dx
         =4852/229635 < 1/40.

In particular integral t^2=1/30 and integral t^3=1/140 are exact.
Uniform absolute convergence justifies termwise bounds and the x
integral throughout the closed ball. Consequently

    |R_loop(k)| <= m0^4/(64*pi^2*h^2)
                      * C_abs*rho^3/(1-rho/4)
                 <= kappa^3/(1920*pi^2*h^2*m0^2).

The second inequality uses 1-rho/4>=3/4 and C_abs<1/40;
it is strict for kappa>0. Both sides vanish at kappa=0.

## Physical scales and what is not transferred

With Lp=M*tau, R=m0*tau and Q=kappa*tau^2, require positive Lp,R
and 0<=Q<=R^2. Since h>=1 and pi>3 (the inscribed hexagon bound),

    |R_loop|/(M^2/tau^2) <= Q^3/(17280*Lp^2*R^2).

For Lp=10^12,R=1000,Q=1 the right side is
1/(17280*10^30). It is the bound for this subtracted quadratic
Fourier kernel. For a real Euclidean source supported within such
a momentum band, |Gamma_sub^(2)| is additionally bounded by that
kernel supremum times integral d^4k/(2*pi)^4 |r(k)|^2. There is
no pointwise energy or nonzero bounce-energy lower bound implied.

The complex ball includes low-momentum Wick-continued k_0=i*omega
with omega^2+|k_space|^2<=m0^2, providing a stationary analytic
continuation of this component. It is not a causal evolving-background
estimate. In particular this calculation does not bound curvature,
other mass insertions, metric/complement/light/gravity loops, finite
higher matching operators, a full in-in stress tensor, an interacting
cutoff, or the quantum-corrected constrained matter cone. Small norm
cannot by itself establish the sign of a saturated luminal cone
inequality. Original V/G/B and original P8 remain open.
