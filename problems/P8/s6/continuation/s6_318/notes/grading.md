# All-order Einstein null grading

## 3. All-order Einstein vertex grading

Temporarily replace delta by a formal variable e in these decompositions.
A spatial rotation by pi in the transverse plane sends e to -e:
transverse momenta and one-longitudinal-index tensors change sign,
while the transverse and double-longitudinal tensor blocks do not.

For a transverse-transverse root test tensor the Einstein vertex is
therefore EVEN in e. Its constant term vanishes for every valence:
all constant child fields are arbitrary transverse tensors with
collinear momenta W_A*k, and the root has the compensating collinear
momentum. Such fields lie in the Rosen metric family

ds^2=2 du dv+gamma_AB(u) dx^A dx^B.

For arbitrary invertible gamma_AB(u), only Ricci_uu can be nonzero:

Ricci_uu=-tr(gamma^-1 gamma'')/2
         +tr(gamma^-1 gamma' gamma^-1 gamma')/4.

Ricci_AB=0, the scalar curvature vanishes, and the Gamma-Gamma action
contraction itself vanishes on the entire family. Thus all its
transverse-root multilinear action coefficients vanish, at every order.
This is a geometric identity, not a finite-N extrapolation.

Consequently the transverse-root source vertex starts at e^2.
A transverse-longitudinal root is ODD in e and starts at e.
A longitudinal-longitudinal root needs no vanishing.

The S315 Einstein bound r!*32^r*L^2*product(field norms), with canonical
factor kappa^(1-r/2), extends to the coefficient-l1 Banach algebra.
Every step of its proof uses submultiplication and triangle/Cauchy bounds,
which also hold for polynomial coefficient-l1 matrix and momentum norms.

For a root with k children, r=k+1, the common source-component budget is

B_pi=16*W^2*(k+1)!*32^(k+1)*13^k
      *kappa^(-(k-1)/2)
      *product_(A in pi)[(W/W_A)*N(H_A,Q_A)].

For |e|<=1 and evaluation e=delta, the complete sum B=sum_pi B_pi obeys

|Rxx|,|Rxy|,|Ryy| <=delta^2*B,
|Rxz|,|Ryz|       <=delta*B,
|Rzz|             <=B.

Component bounds follow by testing the vertex on unit-Frobenius root
tensors; offdiagonal normalization only improves these conservative
bounds. Every complete root partition is retained.


## Why the all-order premise is legitimate

The transverse metric block in the Rosen family is arbitrary and need
not be traceless or solve the uu Einstein equation. Only its transverse
Euler components and its Gamma-Gamma action contraction vanish. The
independent Ricci calculation deliberately retains the generally nonzero
uu component. No nonlinear vacuum condition has been silently imposed
on each child tensor.

The reflection argument is a proper spatial rotation, not a heuristic
parity choice for one helicity. It applies to arbitrary real or complex
symmetric tensor coefficients and to every covariant action coefficient.
Momentum conservation holds throughout the auxiliary polynomial grading.

The coefficients of the auxiliary polynomial are held fixed at the
configuration being estimated. Lower currents need not solve their
field equations at other values of that auxiliary variable: the
grading is an identity of the covariant vertex for arbitrary tensors.
Only its evaluation at the actual angular variance uses the complete
current's Noether conservation. The finite metric jets are multilinear
polynomials in their field inputs, so no new inverse in the auxiliary
variable is introduced by this coefficient-norm argument.

Only the COMPLETE source invokes the Noether identity. Individual
vertices need not be conserved, and an off-shell collapsed subtree is
never treated as a physical free external wave. These distinctions
preserve the S316 obstruction rather than assuming it away.
