# Whole minimal matter endpoint, counterterms and curved anchors

Use Gamma=2PP F1+Q F2, Q=qq-eta*t, P.q=0,
P^2=mu-t/4. All endpoint contributions are first order in
the formal matter-loop marker. F1=1,F2=-1/2 at tree level.

## Both full massive stress triangles

With z in[0,1] and v in[-1,1], the active/spectator masses are
(a,b)=(mu,n) or(n,mu). Let alpha=(1-z)v and

Delta=(1-z)a+zb-z(1-z)mu
 -(1-z)^2(1-v^2)t/4.

After the literal momentum shift, the active momenta are
zP+(alpha-1)q/2-ell and zP+(alpha+1)q/2-ell.
The entire isotropically averaged stress numerator is

2z^2 PP+z alpha(Pq+qP)+(alpha^2-1)qq/2
 +eta[(2/D-1)<ell^2>+a-z^2 P^2-(alpha^2-1)t/4].

The middle Pq term is present before integration. The denominator
is even in v, so the explicit symmetric v reflection removes
only that term. The remaining integral may then use0<=v<=1.
Independent off-diagonal PP and diagonal-difference qq projections
recover z^2 and(alpha^2-1)/2. The complete four-dimensional raw
and even tensors are checked, not only selected scalar coefficients.

The projected triangles are UV finite; their D4 projections are
therefore appropriate after the complete regulated Ward terms
are retained. With c=g^2/(16pi^2), their F1 is the frozen S286
function, including its entire on-shell subtraction:

f1(t)=c integral dzdv (1-z)z^2
 sum_(a,b)[1/Delta_t-1/Delta_0].

The complete triangle F2 plus OS kinetic term is

f2_triangle,OS(t)=c integral dzdv(1-z)[(1-z)^2v^2-1]/(2Delta_t)
 summed over both(a,b), plus Pi_gprime(mu)/2.

At t0, reverse z for the second triangle. Including the full
kinetic counterterm, the combined one-dimensional integrand is
exactly-1/[3F], F=mu(1-z)^2+nz. Thus

f2_triangle,OS(0)=-g^2 integral_0^1 dz/F/(48pi^2).

This constant is retained in the vertex, although its contribution
to the fully crossed endpoint b20 is zero.

## Entire dimensional quartic and H-mixing bubble

For the bubble shift k=ell+xq put y=x(1-x),
Delta=mu-y t. Before the D limit the full isotropic eta
coefficient is

(2/D-1)[A0(Delta)+Delta B0(Delta)]+(mu+y t)B0(Delta).

The exact dimension-dependent identity
A0=2Delta B0/(D-2) makes this2y t B0.
The qq coefficient is-2y B0, so the complete bubble tensor is
-2y(qq-eta*t)B0. No evanescent finite term is discarded.

The quartic graph has its identical internal factor1/2.
The literal vertex/propagator/loop phases give the stripped
coefficient+C/(32pi^2) multiplying that tensor.
The H-metric bubble has the same factor with C replaced by g.
Its heavy propagator and external H Phi Phi vertex multiply it
by g/(n-t). Hence the complete bubble coefficient in F2 is

-(C+g^2/(n-t)) integral_0^1 y B0(mu-y t)dx/(16pi^2).

It includes the direct quartic and H-mixing contribution, not
three direct heavy exchanges. The two crossed heavy exchanges
belong to the active-light triangle, as the cut proof checks.

## Whole Ward and covariant counterterm completion

The full two-triangle identity, both metric-cubic contacts and
their signs are the frozen S286 off-shell Ward calculation:
q.Gamma_loop=p Pi(r^2)-r Pi(p^2).
On shell this is-Pi(mu)q, not zero before mass matching.

For deltaZ=-Pi_prime(mu), delta_m^2=Pi(mu)-mu Pi_prime(mu),
the literal covariant countervertex is

Gamma_ct=deltaZ[pr+rp-eta(p.r)]+delta_m^2 eta
        =-Pi_prime(mu)[2PP-Q/2]+Pi(mu)eta.

Every component and its Ward contraction+Pi(mu)q are checked.
Its PP and qq pieces are the full F1 subtraction and
+Pi_prime(mu)/2 in F2. The remaining eta term cancels the
unsubtracted on-shell Ward term, including finite values.
For the quartic, the entire metric tadpole-Pi_C eta cancels
the OS mass term+Pi_C eta. The original heavy-onepoint
counterfunctional cancels its flat constant density and all
onepoint-reducible attachments after metric variation, as
explained in notes/source.md. The transverse H response remains.

Lorentz covariance of the two-scalar metric endpoint leaves
PP,qq,eta and the parameter-odd Pq structure. The latter has
been removed only by its proved integration symmetry. With
the two independent projections and the whole Ward identity,
the full remaining eta component is fixed. Internal gravitational
gauge-fixing graphs are not needed for this matter-only loop
sector; they remain in separate gravitational-loop sectors.

## Renormalized representation without choosing curved data

Set

Bbar(t)=-integral_0^1 dx y log(1-t y/mu), Bbar(0)=0.

At the old flat MS scale1 write Bconst=B0_MS(mu)/6=-log(mu)/6.
The complete parameterized minimal endpoint is

F1=1+f1(t),
F2=-1/2+f2_triangle,OS(t)
 -(C+g^2/(n-t))Bbar(t)/(16pi^2)+ell+h/(n-t).

For the covariant coefficients c_RPhi2 R_old Phi^2+c_RH R_old H,
the literal linear curvature is R_old^(1)=2Q:hmetric/sqrt(kappa).
Differentiating the two identical Phi fields, or one H, and
stripping-i/sqrt(kappa) gives the factors-4 and-2. Thus

ell=-4c_RPhi2-C Bconst/(16pi^2),
h=-2g c_RH-g^2 Bconst/(16pi^2).

At actual mu1 the displayed loop constant Bconst vanishes,
but neither curved coefficient is fixed. The expression is
a decomposition into known functions and free matching
coordinates, not a prescription setting them to zero.

The analytic functions are continued from below threshold.
Physical values mean explicit t+i delta,delta->0+ limits
of the whole function, including the formal H denominator.
Near t=n the fixed-loop expression requires separate
width/resummation treatment; it is not an exact stable
heavy atom. Additional higher-derivative EFT matching and
other loop sectors remain outside this minimal endpoint.
