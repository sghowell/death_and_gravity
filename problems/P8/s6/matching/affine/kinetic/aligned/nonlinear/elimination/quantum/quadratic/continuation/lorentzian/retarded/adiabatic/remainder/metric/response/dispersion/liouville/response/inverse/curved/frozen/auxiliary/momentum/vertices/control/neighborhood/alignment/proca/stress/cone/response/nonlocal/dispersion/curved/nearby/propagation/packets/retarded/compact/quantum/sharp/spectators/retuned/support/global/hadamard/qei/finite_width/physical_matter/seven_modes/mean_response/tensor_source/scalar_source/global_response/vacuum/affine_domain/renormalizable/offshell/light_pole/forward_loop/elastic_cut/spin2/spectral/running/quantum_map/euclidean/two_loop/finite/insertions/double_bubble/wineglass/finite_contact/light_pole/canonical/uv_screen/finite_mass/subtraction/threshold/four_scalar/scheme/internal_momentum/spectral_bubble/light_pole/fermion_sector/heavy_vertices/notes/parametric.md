# Full unequal-mass parameter domain and master bounds

On |s-2|<=1 let u=4-s and t=0. S6.113 derives the full
box parameter polynomial from the actual routes
0,P,p1,p2. Assign the spectral mass squared v>=T>=16
to the first light line, with parameter alpha0. Write
H=alpha2+alpha3 and Lp=alpha0+alpha1, Lp+H=1.
The polynomial becomes

    Delta=Lp^2+M H-s alpha0 alpha1-u alpha2 alpha3
          +(v-1)alpha0.

Because Re(s),Re(u)<=3 and the pair products are at
most Lp^2/4 and H^2/4,

    Re Delta >= 1/4+(M-1)H+(v-1)alpha0+H(1-H)/2
              >=1/4+(M-1)H+(v-1)alpha0 > 0.

For the other insertion position replace alpha0 by
alpha1. Triangle and bubble faces inherit this bound.
Repeated heavy routes remove a potentially negative
pair term. In the t channel the light-pair transfer
is zero and the heavy-route difference carries s or
u; the same lower bound is conservative. External
heavy factors are analytic for M>=24, with
|C(z)|<=V=L+g/(M-3).

This is continuation of complete Euclidean diagrams
through their Feynman parameters. It does not assume
a shift of a real loop contour after substituting
complex external momenta.

Strip off the common four-dimensional loop factor
1/Q, Q=16 pi^2. For a three-propagator triangle the
Gamma factor is Gamma(1)=1; for a four-propagator box
it is Gamma(2)=1. Fix the marked parameter x. The
remaining simplex volumes are 1-x and (1-x)^2/2.
Dropping only the positive heavy-mass part of the
gap, put a=1/4 and k=v-1. Then

    |J3(v)| <= integral_0^1 (1-x)/(a+kx) dx
              <= log(4v-3)/(v-1) <=2 log(4v)/v,

    |J4(v)| <= (1/2) integral_0^1 dx/(a+kx)^2
              <=2/(v-1) <=4/v.

These bounds intentionally sacrifice heavy-mass
suppression; the spectral-mass suppression already
suffices for the actual candidate. They never expand
the internal momentum relative to either mass.

For the paired bubble define J2(v,z)=B(v,z)-B(v,0).
The frozen S6.135 half-plane derivative estimate gives
|B'(v,z)|<=1/[2(v-4)] along the segment from zero to z.
Since |z|<=3 in the relevant channels,

    |J2(v,z)|<=3/[2(v-4)].

The closed disc has a strict uniform parameter gap,
and these finite pieces have a holomorphic extension
to a neighborhood of it.
