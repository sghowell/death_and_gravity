# Literal source and regular homogeneous response

Set all classical vector components to zero and restrict the
literal local quadratic scalar density to zero spatial
derivatives. Do this before dividing any scalar-shift equation
by q. The resulting density divided by a³ is

    L=-3v'²+(J+w²/2-3Theta²+delta_J)n²
        +6Theta*n*v'+s'²/2+w*n*s'-3ell*v'*s
        +F_v*v+F_n*n.

Here s is the free-matter perturbation, not sqrt(-x), and
delta_J=4epsilon/h². The metric point chart gives physical
curvature v+delta*n at first order, delta=1/(2h).
Direct contraction of 1/2 T^{mu nu} delta g_{mu nu} gives

    F_v=3p_quantum, F_n=3delta*p_quantum-rho_quantum.

Dropping the lapse-dependent physical spatial variation would
miss 3delta*p. The actual clock Euler source is already fixed
by the same action and S6.54's
J_clock=-[rho'+3H(rho+p)]. Diffeomorphism compatibility is
therefore retained; no new clock source is chosen independently.

## Constraints before the crossing

Let p_v, p_m be the momenta divided by a³. Direct Legendre
transformation, with J_e=J+delta_J, gives

    H=-(p_v+3ell*s)²/12+p_m²/2
        +[Theta(p_v+3ell*s)-w*p_m-F_n]²/(4J_e)-F_v*v,
    n=[Theta(p_v+3ell*s)-w*p_m-F_n]/(2J_e).

There is no inverse Theta or spatial momentum. The negative
unreduced gravitational trace coefficient is not promoted to
a ghost statement about inhomogeneous physical modes.

Canonical-density evolution gives the friction -3H on
normalized momenta. Using ell'=-3H ell, the equations imply

    [a³(p_m+3ell*v)]'=0.

Thus zero initial independent canonical perturbations imply
p_m=-3ell*v. They do not imply p_m stays zero. Nor is n an
independent initial datum: n(u0)=-F_n(u0)/(2J_e(u0)).

Set S=p_v+3ell*s, Y=4a³ S. The remaining system is

    [v,Y]' = M [v,Y] + f,
    M=[ B       A/(4a³) ]
      [ -4a³ C     -B  ],
    f=[ -Theta*F_n/(2J_e),
         4a³(F_v+3ell*w*F_n/(2J_e)) ],

where A=Theta²/(2J_e)-1/6,
B=3Theta*w*ell/(2J_e),
C=9ell²[1+w²/(2J_e)].
The code derives every entry from the full action rather than
postulating this favorable normalization. The remaining matter
displacement obeys s'=-3ell*v-w*n.

## Continuous bounds and uniqueness

On I, the inherited all-time J*h²>=1199/800 gives J_e>3/8.
Also 1<=a³<4, |H|,|Theta|<=2, |ell|<=1/10,
|w|<=1/20 and 0<delta<=1/2. In particular
h<3 makes |Lambda|<=1/2 and w=-ell*Lambda.

In the maximum norm the two row sums of M are bounded by

    283/200, 928/625,

both below 3/2. With |rho|,|p|<=eta0,
|F_n|<=5eta0/2, |F_v|<=3eta0; both force components are
below 49eta0. Continuous linear-ODE existence and uniqueness
apply on the entire interval, with zero initial v,Y. The
fundamental matrix norm is at most exp(3|u-s|/2).
An exact exponential-series tail bound proves exp(3/2)<5.
Since the interval length is one,

    ||(v,Y)|| <=49eta0*(exp(3/2)-1)/(3/2)<131eta0,
    ||(v,Y)'|| <246eta0, |v'|<193eta0.

Reconstruction yields

    n=c_v*v+c_Y*Y-F_n/(2J_e),
    c_v=3w*ell/(2J_e), c_Y=Theta/(8a³J_e).

The coefficient sum is at most 103/150, giving |n|<94eta0.
Integration of the matter equation with s(u0)=0 gives
|s|<=44eta0. Signed sources are allowed; no positive-energy
or equation-of-state assumption is used.

## Time derivative control

Absolute polynomial coefficient bounds, retaining the positive
(1+u²) denominators, give

    |J_e'|<185, |Theta'|<7, |w'|<4/3,
    |ell'|<=3/5, |delta'|<=3/2.

The J_e bound includes |delta_J'|<=24epsilon for the whole
allowed epsilon interval. Direct differentiation then gives

    |c_v'|+|c_Y'| <=77792/225<346,
    |F_n'| <=9eta0/2+5eta1/2.

Combining these with the phase and source bounds gives
|n'|<48000eta0+4eta1. The constants are deliberately
conservative but explicit, not fitted to sampled trajectories.

The theorem applies to C1 source functions in the reduced
equations. The actual fixed-background vector source has the
needed regularity by S6.55. Its all-order preparation additionally
allows smooth fixed-background readouts: for each desired time
derivative, compare to a sufficiently high WKB order. Products
have the common full asymptotic series and rapidly decreasing
mixing. Subtracting orders zero, two and four leaves an integrable
order-omega^-5 tail after that fixed derivative; finite local
matching terms are smooth. This supplies qualitative regularity,
not numerical higher-derivative or functional-response bounds.
