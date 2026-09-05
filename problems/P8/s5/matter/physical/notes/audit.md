# Independent quadratic target, regressions and claim boundary

## Independent full quadratic density

Expanding the original unreduced metric/matter action, including background
pressure `F0=-2*Hdot-3*H^2-l^2/2`, gives after the paired mixed boundary

`sigma1=(p+3*l*s)/2`, `eta1=P`, `rho1=4*q*v`, `z2=q*s^2`.

Solving the linear matter-sourced shift gives integrated scalar
`shear2_2=p^2/24`. In particular it is p, not p+3*l*s, that occurs in this
shear invariant: the matter source cancels the extra longitudinal contribution.
The combination `-sigma1^2/3+2*shear2_2` therefore gives the terms
`-l*p*s/2-3*l^2*s^2/4` which a matter-free reduction would miss.
Background volume and canonical boundary terms cancel using the independent
pressure equation. The full scalar quadratic Hamiltonian is

`H2=-q*v^2+P^2/2+(q/2-3*l^2/4)*s^2-l*p*s/2`

` +[-Theta*p/2+Lambda*q*v+w*P/2-3*l*Theta*s/2]^2/J`.

Applying the gamma swap and its generator gives

`H2_gamma=-P_b^2/(4*q)+P^2/2+(q/2-3*l^2/4)*s^2`

` +q*l*b*s-H*b*P_b`

` +[q*Theta*b+Lambda*P_b/2+w*P/2-3*l*Theta*s/2]^2/J`.

The tensor block is `2*PiTT:PiTT+q*gammaTT:gammaTT/8`, with no scalar/tensor
quadratic cross. Every independent scalar Hessian entry in both charts,
tensor entry, mixed entry and homogeneous unitary tadpole is compared with
the nonlinear construction at symbolic compact time: 35 exact residuals.
Separate literal rational-time checks cover the bounce, an interior point
and both compact tail limits. An independent review also compared all forty
scalar entries at two times in both charts and all scalar velocity entries.

## Finite-q inverse Hessian and the principal limit

For the scalar momentum Hessian A, exact determinants are

`det A_unitary=Theta^2/(2*J)`,

`det A_gamma=(q*Lambda^2-J0)/(2*J*q)`, `J0=J+w^2/2`.

With `D=q*Lambda^2-J0`, the gamma kinetic matrix `K=A^-1/2` is

`K=[[q*J0/D,-q*Lambda*w/(2*D)],`

`   [-q*Lambda*w/(2*D),(q*Lambda^2-J)/(2*D)]]`.

For J>0 its positive chart requires D>0. Its q->infinity limit recovers the
earlier principal matrix, but the finite-q lower-right entry is not generally
1/2. At the bounce the velocity pole is q=6. An exact test checks a finite
phase coefficient there while the velocity API rejects that pole and the
non-positive q=5 chart. No low-q stability verdict follows from this chart test.

## Omission and symmetry controls

The report contains explicit nonzero errors for deleting the mixed boundary,
deleting the matter background-velocity subtraction, using the vacuum instead
of matter-sourced momentum constraint, and retaining only W1 in a quartic
interaction. The last gives error `-3/100` for the specified bounce fixture.
Previous invariant geometry/lapse errors remain pinned by S5.5 and S5.2.D.

The principal fixture momenta are

`TRI=((2,3,0),(0,2,4),(-2,-5,-4))`,

`QUAD=((2,3,0),(0,2,4),(3,-2,1),(-5,-3,-5))`.

Every proper-subset squared transfer is greater than six, so the bounce gamma
velocity chart is valid for all external and quartic internal momenta. Tests
check exact parity, label permutation and reality, invalid polarization and
zero-transfer rejection. A mixed quartic fixture has nonzero contributions
from both internal tensor polarizations; an orthogonal rational change of TT
basis preserves each partition's tensor total and the full vertex.

At the bounce the all-matter-velocity quartic has

`L4=15739581640021/33058771200000`,

`-H4(P0)=-3269548239/12940480000`.

Each of its three scalar partitions has a nonzero off-diagonal contact. The
independent arbitrary symmetric two-scalar completion checks its sign; a
separate review checked the unordered-partition factor with a two-species toy.
These are exact regression fixtures, not phenomenological observables.

## Boundaries of evidence

Generic Fourier/York and weighted stationary identities plus the written
reduction lemmas justify the formal construction for arbitrary admitted
inputs. Selected fixtures are not an exhaustive enumeration of momentum space.
The all-time fixed-momentum phase regularity result is the finite-algebra lemma,
not a sampled extrapolation; no expanded symbolic quartic majorant is promoted.

All prior sources are hash-pinned and unchanged. This checkpoint derives no
coupled canonical mode normalization, positive-frequency evolution estimate,
scattering amplitude, interacting cutoff, radiative correction, UV positivity
condition, global nonlinear stability or full P8 completion. In particular,
none of the prior D-only hard-tree constants is silently transferred to M1.
