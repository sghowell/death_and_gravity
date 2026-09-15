# Entire dimension-dependent covariant local projection

Write X=partial_m Phi partial^m Phi and Y=Phi^2. The eight
coordinates(a,b,c,d,e,f,r,w) multiply

(a X^2+b mu Y X+c mu^2 Y^2)/kappa^2
 +(d R_old X+e Ricci_old_mn partial^m Phi partial^n Phi
    +f mu R_old Y)/kappa+r R_old^2+w Weyl^2.

All are arbitrary renormalized or raw matching coordinates as
specified by their use; none is set to zero. Evaluating the first
counterterm variation of the on-shell action requires its value
on the zeroth-order classical solution: the variation of the
zeroth-order action itself is zero with the fixed on-shell
external data. This is a one-insertion calculation. It is not
an actual change of integration variables or of the original
physical matter frame.

## Conserved sourced metric and Euler kernel

For the ordinary metric convention the scalar stress is
T_mn=partial_m Phi partial_n Phi-eta_mn(X-muY)/2.
Its trace is-(D-2)X/2+DmuY/2. The conserved sourced metric is
gamma_mn=2[T_mn-eta_mn T/(D-2)]/(kappa q^2). Its literal linear
curvature in the original R_old convention gives

Ricci_old_mn=-(partial_m Phi partial_n Phi-muY eta_mn/(D-2))/kappa,
R_old=(-X+DmuY/(D-2))/kappa.

The code independently varies every symmetric metric component
in D4,D5,D6 at q=(omega,0,...), constructs the Ricci and Riemann
tensors, and verifies the complete conserved response. Arbitrary
symmetric spatial T spans the conserved source at timelike q.
Lorentz covariance and polynomial continuation extend the tensor
identity away from this frame; the written D-dependent contraction
above gives the general D formula, not an inference from three
integer dimensions.

For a general metric perturbation let H be its spatial block.
The complete quadratic invariants at the same timelike q are

Riemann^2=omega^4 tr(H^2),
Ricci^2=omega^4[(tr H)^2+tr(H^2)]/4,
R^2=omega^4(tr H)^2.

Hence the whole quadratic Euler kernel E4=Riemann^2-4Ricci^2+R^2
vanishes for every D, including all h00 and h0i components.
These are polynomial momentum identities, so they extend also
to the null locus without division by q^2. This does NOT claim
that Euler is topological in general dimension. The independent
one-counterterm tree counting in notes/source.md explains why
only this quadratic kernel is needed here.

Use Weyl^2=E4+4(D-3)Ricci^2/(D-2)
-D(D-3)R^2/[(D-1)(D-2)]. The sourced invariants reduce to

kappa^2 Ricci^2=X^2-2muXY/(D-2)+Dmu^2Y^2/(D-2)^2,
kappa^2 R^2=(-X+DmuY/(D-2))^2,
kappa Ricci_mn partial^mPhi partial^nPhi=-X^2+muXY/(D-2).

## Literal scalar vertices and full map

Let Phi=sum_i phi_i exp(ik_i x), with four incoming on-shell
k_i^2=mu and sum_i k_i=0. The full Gram matrix has
k1.k2=(s-2mu)/2 and the t/u crossed entries. Extracting the
phi1 phi2 phi3 phi4 coefficient gives, including all permutations,

X^2 -> 2(s^2+t^2+u^2-4mu^2),
mu Phi^2 X -> 8mu^2,
mu^2 Phi^4 ->24mu^2.

Thus Phi^2 X has the same four-point insertion as mu Phi^4/3.
The literal extraction, not an assumed off-shell EOM equivalence,
checks that statement.

The output amplitude is[A2 sum(channel^2)+A0 mu^2]/kappa^2,
where(A2,A0)^T=M_D(a,b,c,d,e,f,r,w)^T and the full rows are

A2:[2,0,0,-2,-2,0,2,2(D-3)(3D-4)/((D-2)(D-1))],
A0:[-8,8,24,16(D-1)/(D-2),8(D-1)/(D-2),
16(D+1)/(D-2),32(2D-1)/(D-2)^2,
-32D(D-3)/((D-2)(D-1))].

At D4 these become[2,0,0,-2,-2,0,2,8/3] and
[-8,8,24,24,12,40,56,-64/3]. The(a,c)minor is48.
The six explicit columns in local.nullspace_columns obey
M_D N_D=0; their rows(b,d,e,f,r,w) have determinant1.
Consequently the rank is exactly2 and these six columns span
the whole kernel for the stated D domain D>2. This is an
on-shell degeneracy, not a determination of off-shell coefficients.

At fixed t, s=2mu-t/2+v and u=2mu-t/2-v, so the v^2
coefficient is2A2/kappa^2. At D4 it is

b20_local=4(a-d-e+r+4w/3)/kappa^2.

## Evanescent and physical-frame boundary

D=4+2EP must be retained before multiplying raw1/EP terms.
The EP-linear columns of(r,w) in the two output rows are

[[0,26/9],[-80,-160/9]].

These generate finite constants and are not discarded. The full
matrix derivative is supplied for all columns. Using its D4
value is appropriate only for coefficients already renormalized
in the specified convention; it is not a raw-dimensional shortcut.

The linear EOM/insertion simplification does not authorize a finite
nonlinear field redefinition or an all-order quantum equivalence.
The order-by-order distinction and the accompanying source/matching
changes are discussed by Criado and Perez-Victoria,
[Field redefinitions in effective theories at higher orders](https://arxiv.org/html/1811.09413).
Here no such change is made: higher-field terms, the measure,
curved off-shell coefficients and the original bounce frame remain.
