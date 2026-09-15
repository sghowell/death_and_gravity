# Closed finite forward coefficient and explicit uniform remainder

## Entire D-dependent center derivatives

At s=2mu+v,t=0,u=2mu-v, a single-channel function contributes its
second derivative at2mu to the coefficient of v^2: the two Taylor
halves add. The t-channel value and all constant leg/local terms
have no v dependence. Set y=x(1-x), A0=mu(1-2y)>=mu/2.

Differentiate the complete parameter integrands before expanding in e.
For the proper term the exact second channel derivative is
T2(e)=integral_0^1[
 (1-2e*y)A0^(e-1)
 -mu^2(e-1)(e-2)*y^2*A0^(e-3)/(1+e)]dx.
For the endpoint it is
E2(e)=integral_0^1[
 -2e*y^2*A0^(e-1)
 +2mu(2+e)*e(e-1)*y^3*A0^(e-2)/(1+e)]dx.
The entire raw coefficient is therefore
C Gamma(-e)(4pi nu^2)^(-e)[-2T2(e)+E2(e)]/(16pi^2 kappa).
The complete positive parameter gap justifies differentiation and
Taylor expansion on the compact domain. No small mu/n expansion
or physical-root singular quadrature enters this crossing-center proof.

## Exact finite integrals

At zero e, the proper integral is
T2(0)=integral[1/A0-2mu^2*y^2/A0^3]dx=3pi/(8mu)=B2.
Its first derivative is
T2'(0)=integral[
 ln(A0)/A0-2y/A0+mu^2*y^2(5-2ln A0)/A0^3]dx.
The endpoint gives E2(0)=0 and
E2'(0)=-integral[2y^2/A0+4mu*y^3/A0^2]dx
      =-(1-pi/4)/mu.

For the full proper integral use reflection and x=(1+z)/2.
The coefficient of ln mu is
fB(z)=(1+6z^2+z^4)/(1+z^2)^3.
A primitive is
3atan(z)/2-z(1-z^2)/[2(1+z^2)^2],
so its integral is3pi/8. The remaining rational term integrates to
1-3pi/16. The logarithmic integral is
integral_0^1 fB(z)ln[(1+z^2)/2]dz
=3pi ln2/8-3Catalan/2+1/4-pi/16.

To derive it put z=tan theta,0<=theta<=pi/4.
Then fB dz=(3/2-cos4theta/2)dtheta.
The defining Catalan integral
Catalan=-integral_0^(pi/4)ln(tan theta)dtheta,
together with the sine doubling identity and
integral_0^(pi/2)ln(sin theta)dtheta=-pi ln2/2,
gives integral ln cos theta=-pi ln2/4+Catalan/2.
Integration by parts gives
integral cos4theta ln cos theta=1/4-pi/16.
Every integral is convergent; no finite endpoint constant is introduced.
Thus
T2'(0)=[3pi ln(2mu)/8-3Catalan/2+5/4-pi/4]/mu.
The endpoint rational integral independently equals S6.290's full
quartic contribution; it has not been counted twice.

## Retained formal soft division

Use only the already inherited S6.278 convention
logW_E=Bsoft*(E/nu)^(2e)/(8pi^2 kappa e).
At the coefficient linear in C, its v^2 term is
C*B2*(E^2/nu^2)^e/(8pi^2 kappa e).
Subtract it from the complete raw graph coefficient. The finite result is
C/(8pi^2 kappa mu)*
 [3pi/8*(EulerGamma+ln(mu/(2pi E^2))-1)
  -3Catalan/2+7/4].
The regulator scale cancels, while the specified E dependence remains.
Its derivative in ln E equals the inherited resolution-flow coefficient.
The sign of this subset depends on E; no E is selected to force a sign.

This is a computation within a named conditional analytic convention.
It is not a construction or proof of physical soft-factor
factorization/unitarity, nor a new detector/finite-subtraction choice.

## Evaluated uniform regulator bound

Now restrict to original mu=nu^2=1,1/4<=E^2<=1 and0<=e<=1/8.
Let
R_e=Gamma(1-e)(4pi)^(-e)[2T2(e)-E2(e)]
    -2B2*(E^2)^e.
R0=0. The complete finite-e soft-divided result is C*R_e/(16pi^2 kappa e).
Taylor's theorem bounds its difference from the finite limit by
e*abs(C)*sup|R''|/(32pi^2 kappa).

Here is an explicit bound for that supremum. On0<=x<=1,
0<=y<=1/4,1/2<=A0<=1, |ln A0|<1, A0^(e-1)<=2,
y^2/A0^2<=1/4 and y^3/A0<=1/32.
Write g(e)=(e-1)(e-2)/(1+e)=e-4+6/(1+e).
Then |g|<=2, |g'|<=5, |g''|<=12 on the stated interval.
The proper integrand is A0^(e-1)*(1-2e*y-g*y^2/A0^2).
The bracket and its first/second derivatives are bounded by(2,2,3).
After the exponential product rule, its integrated derivative bounds
are T2:(4,8,18).

For the endpoint write h(e)=(e+2)(e-1)/(1+e)=e-2/(1+e),
with |h|<=2, |h'|<=3, |h''|<=4.
Its integrand is A0^(e-1)*e[-2y^2+2h*y^3/A0].
The bracket including e has bounds(1/32,1/2,1/2).
The product rule then gives E2:(1,2,4), strictly larger than the
corresponding elementary estimates.

The Gamma integral split at1 yields, for j=0,1,2,
|Gamma^(j)(1-e)|<
 integral_0^1 t^(-1/8)|ln t|^j dt
 +integral_1^infinity t^j exp(-t)dt
 <j!/(7/8)^(j+1)+j!.
These bounds are respectively below(3,3,5).
Since0<ln(4pi)<3 and(4pi)^(-e)<=1, the Gamma/scale product and
its first two derivatives are bounded by(3,12,50).
Hence the product with2T2-E2 has second derivative below
50*9+2*12*18+3*40=1002.
The subtracted soft factor adds at most12 since
2B2<3, |ln E^2|<2 and(E^2)^e<=1.
Therefore sup|R''|<1014<1024. The complete regulator error is
below32e*abs(C)/(pi^2 kappa)<8e*abs(C)/kappa.
These are explicit uniform constants, not an unevaluated supremum.
They bound this regulator calculation, not a physical detector error.

## Actual finite magnitude

For the same original resolution range,
-3<ln[1/(2pi E^2)]<0,0<EulerGamma<1,0<Catalan<1 and2<pi<4.
The absolute finite bracket is below6+7/4<8.
Thus abs(b20_C)<abs(C)/(4kappa).
The unchanged original contact obeys
4g^2/n-abs(C)=g^2(n-4)^2/[n(n-2)^2]>0,
and its exact parameters give g^2/(kappa*n)<10^-1005.

Both estimates cover the whole stated E range and do not fix a sign.
Finite constant quartic and constant curvature anchors have zero b20
at this insertion; higher matching and the full source do not follow.
