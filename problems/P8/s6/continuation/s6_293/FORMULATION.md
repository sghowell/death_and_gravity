# Formulation and limited claim

## Fixed theory and domain

Keep the original affine parent, flat vacuum, canonical mass and
coupling conventions. C denotes the original nondifferentiated
four-scalar tree amplitude, not a newly selected coupling. Work to
one loop, linear in C and in1/kappa, in the inherited harmonic metric
gauge and D=4+2e. The inverse scalar propagator is p^2-mu+Sigma.
Use the Feynman normal sheet and the raw negative-Euclidean C0
convention of S6.283.

The analytic center has mu>0, s=2mu+v,t=0,u=2mu-v. The physical
two-scalar cut has s>4mu,t,u<0 and0<e<=1/8. Its singular angular
integral is regulated before removal of e. The explicit uniform
error and numerical magnitude bounds below require original
mu=nu^2=1 and1/4<=E^2<=1. These are bounds over a stated range of
the inherited formal resolution parameter; no physical detector
value is selected.

## Whole selected amplitude

For each a in(s,t,u), set y=x(1-x), A_a=mu-a*y-i0,
V_D(a)=(a-2mu)^2-2mu^2/(1+e).
Define M_a(e)=integral A_a^e dx and J_a(e)=M_a(e-1). Then
T_e=sum_a[V_D(a)J_a(e)/2+(a-2mu)M_a(e)]
    +mu^(1+e)/(1+e),
E_e=sum_a[a+2mu/(1+e)]integral y*A_a^e dx.
The complete known graph sum is
deltaA_C=C*Gamma(-e)(4pi nu^2)^(-e)*[-2T_e+E_e]
         /(16pi^2 kappa).
Six pair triangles, four contacts, four external factors, both
endpoint assignments in all three channels, the full on-shell
mass terms and the entire flat onepoint condition are accounted for.

T0 is the complete four-leg massive soft coefficient Bsoft.
The selected proper-plus-external UV terms cancel separately.
The matter endpoint has E0=5mu/3, a local momentum-independent UV
pole in the retained OS/curvature organization. Finite constant
quartic/curvature anchors remain unassigned; they have zero b20
at this single insertion. Higher-derivative matching is not fixed.

## Exact known forward coefficient

Let b20 mean the coefficient of v^2 at the stated crossing center.
Subtract only the linear C coefficient of the already specified
S6.278 analytic soft factor:
C*B2*(E^2/nu^2)^e/(8pi^2 kappa e),
B2=3pi/(8mu).
The finite known coefficient is
b20_C(E)=C/(8pi^2 kappa mu)*
 [3pi/8*(EulerGamma+ln(mu/(2pi E^2))-1)
  -3Catalan/2+7/4].
All first e coefficients are retained. nu cancels. E remains explicit,
and changing it can change the sign of this subset.

On the original domain above, the difference between the complete
finite-e soft-divided coefficient and this limit is strictly below
8e*abs(C)/kappa. Its finite magnitude is below abs(C)/(4kappa),
and the unchanged original parameters imply abs(b20_C)<10^-1005.
The uniform remainder is proved with numerical majorants, not an
unevaluated derivative supremum.

## Nonclaims

S6.278's analytic soft-factor convention remains conditional on
physical factorization; this computation does not construct a
unitary dressed/inclusive observable. The C-linear cut is an
interference term, not a positive measure. All other couplings,
higher-EFT matching, Coulomb/other physical infrared contributions,
fixed-transfer Regge and all-loop bounds, state/domain/bounce
matching and original V/G/B/P8 remain open.
