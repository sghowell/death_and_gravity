# Actual kernel at every Euclidean momentum

Keep the on-shell subtraction of S6.112, with light pole mass and
residue one. Let y>=0 be Euclidean momentum squared, M the actual
heavy mass squared, g the cubic coupling squared and

    d(x)=xM+(1-x)^2,
    b(x)=x(1-x)/d(x),  0<=x<=1.

The complete heavy denominator is retained. The one-loop result is

    Pi_R(-y)=g/(16pi^2) integral_0^1
             [b(x)(y+1)-ln(1+b(x)(y+1))] dx.

Define Q(y)=Pi_R(-y)/(y+1) and

    f(u)=1-ln(1+u)/u = integral_0^1 u t/(1+u t) dt.

The integral supplies the continuous f(0)=0. Native code checks
its primitive, anchor, derivative and endpoints. For u>0,

    u/[2(1+u)] < f(u) < min(u/2,1).

The lower and linear upper gaps are explicit positive integrals.
Moreover f'(u)=[ln(1+u)-u/(1+u)]/u^2>0: the numerator vanishes
at zero and has derivative u/(1+u)^2. Also f(u) tends to one.

Therefore Q(y)>0 and is strictly increasing, with

    sup_(y>=0) Q(y) = alpha
      = g/(16pi^2) integral_0^1 b(x) dx.

The limit follows by dominated convergence: b is positive in the
open interval, integrable and independent of y, and 0<=f<=1.
The value alpha is approached, not attained at finite y.

Differentiating the original unsubtracted self-energy at the
mass-one shell gives precisely alpha. Independently, summing
both actual stress triangles at transfer zero and integrating
their pair parameter yields the same integrand. Thus alpha is
the magnitude of the already fixed kinetic counterterm, not a
new subtraction constant. The renormalized gravitational charge
remains one, not alpha.

Since b(x)<=(1-x)/M, alpha<g/(32pi^2 M)<g/(288M).
The high-Euclidean-momentum estimate is weaker than the local
quadratic remainder near the mass shell, but is uniform without
a momentum cutoff or a heavy derivative expansion. It says
nothing by itself about complex timelike contours or the
interacting heavy resonance.
