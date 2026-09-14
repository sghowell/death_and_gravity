# Onepoint contact and physical volume cancellation

The physical spatial metric is h_phys=C(u,N)gamma with C=R-1/2.
The full density is

rho_V=a^3 exp(3v) C(u,N)^(3/2).

Under the complete spatial transformation,
delta v=xi dot grad v+div xi/3 and delta N=xi dot grad N.
Consequently

delta rho_V =
rho_V[3delta v+(3/2)(C_N/C)delta N]
=div(rho_V xi).

This identity uses the actual full conformal function, not its reference
value. On the periodic domain the integral is invariant. A nonperiodic
domain retains the boundary flux. No time boundary or lapse equation is
replaced by a gauge condition.

At a homogeneous reference set c1=(3/2)C_N/C. The relevant second-order
density expansion includes

3v+c1 n+(9/2)v^2+3c1 vn+c2 n^2.

The second n derivative c2 is retained in the full expansion but its
gauge variation at this order vanishes because delta n1=0. For the
unchanged current clock jet c1=-6/(1+u^2)^3, hence c1(0)=-6.
Only a first lapse jet is used, licensed by the order-eight exact source
factorization; the full off-clock C is not redefined.

The gauge flow gives delta v1=3v_perp/4 and the mean second-order jets

<delta v2>=-9 Cvv_perp/4,
<delta n2>=-9 Cvn_perp/4.

The onepoint part of the density variation is therefore
-27 Cvv/4-9c1 Cvn/4.
The two-point part is
9<v delta v1>+3c1<n delta v1>
=27 Cvv/4+9c1 Cvn/4.
Their sum vanishes exactly. The cancellation is independent of the
full coupled covariance values; neither Cvn nor either onepoint contact
is set to zero. It holds for symmetrized finite polynomial/Weyl jets as
specified in reference.md, without selecting a nonlinear quantum
ordering for the original source functional.

The exact periodic orbit supplies a nonperturbative CLASSICAL check:
rho_V,pulled=C(N(phi))^(3/2)exp(3v(phi))phi'. The change of variables
y=phi(x) proves equality of the entire original and pulled integrals.
Both tensor profiles and the full lapse dependence are present. The
numerical tests reconstruct the inverse map and the entire 3x3 metric
independently and check the same volume.

Dropping the nonlinear coordinate mean would create a false nonzero
gauge variation of physical volume. Conversely, this cancellation does
not say that every physical observable is already renormalized, that
all S260 Ward defects vanish, or that the interacting P8 mean exists
and lies on the reference branch.
