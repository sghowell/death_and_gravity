# Regular constraints, physical maps, and the unresolved scalar symbol

## 1. Action, clocks, and quadratic variables

Use exactly the S6.20 action and solution in the formulation. In this note
`M=tau=1`; a prime is `d/du`. Restoring dimensionless spacetime coordinates
and scalar fields gives the positive overall action factor `(M*tau)^2`.
No sign/rank conclusion depends on that factor. The P8(b) curvature
convention is `R_B=-6(H'+2H^2)` with `G_B00=3H^2` and Einstein term
`-M^2 R_B/2`; it is not the opposite A-track curvature convention.

Write the two scalar spatial metrics as

```
h_g=a² diag(exp(2qg),exp(2qg),exp(2qg+2eg)),
h_f=b² diag(exp(2qf),exp(2qf),exp(2qf+2e)).
```

The common time gauge is `delta_phi=0`; the common spatial gauge is
`eg=0`. These require `phi'=sqrt(k)>0` and a nonzero Fourier mode, not
`H!=0`. Here `e=-K E_f` is the longitudinal logarithmic scalar shear.
The lapse perturbations are multiplicative, `N_g=1+ng`,
`N_f=c(1+nf)`, while `N_i^z=partial_z B_i`. Set

```
w=chi'=1/(10a³), n=phi'²+w², k=n-w²,
Sg=3qg, Sf=3qf+e, Jv=Sg-Sf, P=2b1,
Cs=a^5*y²*P/[2(c+y)], y=b/a.
```

All background functions here are the actual sourced-clock profiles, not
free clock replacements. `beta2=beta3=0`, but all `beta0,beta1,beta4`
terms and their time dependence remain. No perturbation equation is
obtained by setting the background `h` to zero before differentiation.

The exact intrinsic spatial scalar curvature is

```
sqrt(h) R3 = a exp(q-e)[-4q_zz-2q_z²+4e_z q_z].
```

Including the lapse and Einstein factor `1/2`, its quadratic Fourier
coefficient is `aK(q²+2nq)`: the shear dependence cancels after spatial
integration by parts. Expanding the lapse-retaining extrinsic curvature
density gives the kinetic/constraint terms in `action.lagrangian()`.
The time boundary removed from the Einstein action is

```
F_EH=-a³h Sg²+b³h Sf²/c.
```

Its derivative contributes the retained bulk coefficients `h'+3h²/2`
and `-h'+3h²/2`. The literal potential follows by expanding the diagonal
positive square roots; the relative-shift coefficient follows the
independent two-dimensional root. The root-authored audit directly
expands the Einstein, both scalar and potential densities before the
Legendre transform. It checks the time boundary, intrinsic curvature and
the complete lapse/shift elimination at rational actual backgrounds.

Gauge fixing the clock does not turn its coupling into an external
source. The clock equation follows from the common diffeomorphism
identity once the retained metric and free-chi equations hold, because
the actual clock speed is nonzero. No separate clock-stress conservation
is imposed; the beta(phi) exchange in the frozen action is preserved.

## 2. Exact Hamiltonian and remaining lapse constraint

Let `(pg,pf,pe,px)` be momenta conjugate to `(qg,qf,e,z)`, with
`z=delta_chi`. After the common spatial momentum constraint,

```
pEg=a³ w z-pe,
A0=-h pg+a³ n Sg+a³ yP Jv+2aK qg,
p_phi=(A0-w px)/sqrt(k).
```

The clock and relative-shift constraints give the following Hamiltonian,
still with the f lapse as a multiplier:

```
Htotal=H+nf Cf,
Cf=-h pf-2c bK qf+c a³P Jv,
H=(3pEg²-2pg pEg)/(4a³)+c(3pe²-2pf pe)/(4b³)+K pe²/(4Cs)
  +(A0-w px)²/(2a³k)+px²/(2a³)-Sg A0+a³n Sg²/4+aK z²/2
  -a³(h'+3h²/2)Sg²-(b³/c)(-h'+3h²/2)Sf²-aK qg²-cbK qf²
  +a³(b0+cb1)Sg²+a³yb1[2(2qg+qf)²+(2qg+qf+e)²]+cb³b4 Sf².
```

This is derived by a full quadratic Legendre transform and an auxiliary
Schur complement, not by guessing a constraint at the bounce. In
particular no `nf²` survives. The independent audit retains all four
lapse/shift equations and reproduces `H+nf Cf`, including nonzero positive
and negative times and the `c=1` control.

## 3. A non-H-dividing canonical reduction and its strict domain

Put

```
r=1+2bK/(3a³P), alpha=h/(3c a³P),
qg=Q+rR+e/3+alpha PR,  qf=R+alpha p0,
pg=p0, pf=PR-rp0, pe=PE-p0/3.
```

The matter pair is unchanged. This is an invertible canonical map for
every finite `K` in the stated background domain and gives
`Cf=3c a³P Q`. The old one-form is the new one-form plus `dF` and an
explicit time term, where

```
F=alpha p0 PR-r alpha p0²/2,
G=-r' p0 R-alpha' p0 PR+(r alpha'-alpha r')p0²/2.
```

Thus the new Hamiltonian is the transformed old Hamiltonian **plus G**.
Keeping only the frozen-time symplectic matrix would miss this term.
The first constraint is `Q=0`; its preservation is `dHnew/dp0=0`.
The coefficient of `p0` in that equation is exactly

```
D=-(1+c/y³)/(6a³)+h²(1+y/c)²/(a³k)
  -2h²y/(3c²a³P)+alpha'.
```

It has no `K` dependence. Cancellation of the apparent `K` terms uses
the **actual** f null equation, not the old constant-beta Bianchi lock.
`reduction.checks()` compares this expression with the literal transformed
Hamiltonian Hessian. At the center `D=-5/24` for every admitted `c`.

The stronger bound is `D<-1/8` on the whole local box. To prove it define
`d=1+u²` and

```
Fclock=6400(1-u²)-800c(1-u²)d^12-cd²=100c d^14 k>0,
Den=24c(1-u²)²d^6 Fclock>0,
Num=Den*(-D-1/8).
```

`Num` is even in `u`. Substituting `u²=x/100` gives degree27 in `x`.
For `c=1` all28 Bernstein coefficients on `0<=x<=1` are strictly
positive. For `c=2+2z` it has bidegree(27,3), and all112 Bernstein
coefficients on the unit square are strictly positive. Their exact
rational minima are in the report. The convex-hull property proves the
continuous-domain bound; no finite sampling inference is used. A separate
Fraction/polynomial implementation derives the cleared numerator using
`alpha=u d³(cd⁴-2)/[48(1-u²)]`, performs the affine substitutions, and
reproduces all140 coefficients. The `c=2` polynomial edge is only a
cleared-boundary device, not a regular finite-action solution.

The pair `(Q,dHnew/dp0)` is therefore second class. Solving it leaves
exactly the three canonical pairs `(R,e,z;PR,PE,px)`, with a finite
analytic Hamiltonian for every finite `K>0`. Preservation of the second
constraint determines `nf`; it does not introduce another mode. This is
a quadratic constraint-rank result, not a proof of nonlinear constraint
rank, scalar positivity or well-posedness uniform in momentum.

## 4. Physical observables, lapse recovery and finite-band continuation

Let `Z=(R,e,z,PR,PE,px)`, `Z'=A Z`, `A=Jcan H2`, and let
`p0=L Z` be the solved secondary momentum. The spatial gauge equation
recovers the actual g shift:

```
Bg=(3pEg-pg)/(2a³K),
xi=a² Bg=3(a³w z-PE)/(2aK).
```

All `p0` terms cancel. For a common time shift `xi_time`, the ungauged
combination `a²(Bg-Eg')` shifts by `+xi_time`, while the matter
perturbations shift by minus their background velocities. Consequently
the invariant observable coordinates used here are

```
q=(xi, chi_B, Psi_f_B),
delta_phi_B=sqrt(k)*xi, chi_B=z+w xi,
Psi_f_B=qf-h xi=R+alpha LZ-h xi.
```

The first coordinate is the physical clock perturbation divided by its
nonzero background speed. `R` is not relabeled a physical curvature.
The full recovery also gives

```
ng=Sg-(A0-wpx)/(a³k), Bf=Bg-pe/(2Cs),
nf=-[Hnew,Q|constraints+L'Z+L A Z]/(3c a³P).
```

The last formula is the secondary consistency equation. It divides
`P`, not `h`. The new module reconstructs all old phase entries and their
derivatives and substitutes these rows into every original lapse/shift
equation, every momentum definition and all eight unreduced Hamilton
equations, including secondary-lapse preservation. The residuals vanish
exactly. Merely checking the lapse/shift equations at `h=0` would not
detect a wrongly omitted `nf`; the full phase equations do.
For constant `c`, the lapse Bardeen variables are `ni_B=ni+xi'`, and
`Psi_g_B=qg+h xi`; the remaining f shift in g Newtonian gauge is
`Bf_B=Bf-c²xi/b²`. This is an actual metric/source dictionary, not a
change of physical matter frame.

For `q=O Z` the exact Cauchy and acceleration maps are

```
C = stack(O, O'+O A),
q'' = [O''+2O'A+O(A'+A²)] C^-1 (q,q').
```

Every derivative is taken before time evaluation. At `u=0` the term
`alpha L''` vanishes, but `alpha' L` and `2alpha' L'` do not. At nonzero
time the implementation retains the full second jet. The canonical
conserved antisymmetric bilinear is explicitly `Z1^T Jcan Z2`; in the
physical Cauchy chart it is `X1^T C^-T Jcan C^-1 X2`. This sign convention
is fixed independently of how one names the canonical two-form.

At the center define

```
F_K=c²(c-2)K²+4c(c-2)(c-8)K+1536(8-c).
det C=5F_K/[12K²(c-2)(6400-801c)].
```

For `2<c<=4`, `F_K` is positive for all `K>0`: its minimum in real K
is `4(8-c)[384-(c-2)(8-c)]>0`. Thus this is a genuine finite-momentum
Cauchy chart at the bounce. Its equal-time `q,q` bracket vanishes there.
Its `q,q'` bracket has positive final diagonal `c/240`, positive
2-by-2 determinant `100F_K/[K²c(c-2)(6400-801c)]`, and positive second
diagonal numerator `200c(c-2)K²-c(c-2)K+384`. This proves a finite-time
algebraic property of the conserved bilinear, not a frequency/residue
or scalar-health theorem.

For each fixed compact momentum band `0<Kmin<=K<=Kmax<infinity` and
each fixed admitted positive-branch `c`, compactness and continuity give
a neighborhood of the bounce on which `C` remains invertible. The
regular canonical first-order equation exists across the entire local
background interval because its coefficients are continuous there.
No explicit band-independent neighborhood, propagator bound or scalar
energy estimate is supplied.

## 5. Why the apparent scalar frozen symbols are not health verdicts

At `c=4` the physical Cauchy determinant has exact time jets

```
det C|0=5(K²-4K+192)/(2397K²),  (det C)'|0=0,
(det C)''|0=(19975K³-48440K²+1644744K-313661952)/(5745609K²).
```

In particular `(det C)''/[2K det C] -> 5/6`. The map changes on a
momentum-dependent time scale near the bounce; its second derivative
does not have a uniform order-zero bound in `K`. Away from the bounce
the physical observable coordinates need not commute. An instantaneous
Hamiltonian momentum block is therefore not this observable kinetic
matrix, and an instantaneous frozen observable polynomial is not yet a
uniform physical characteristic theorem.

For reproducibility the code retains the apparent center polynomial
`(v²-1)²[v²-(c-1)(9c²-22c+144)/120]`, up to a nonzero overall sign.
It is recorded **only as a nonuniform diagnostic**, not as three
certified scalar propagation speeds. Earlier ordinary `H`-dividing or
different phase-coordinate frozen computations need not reproduce it:
they omit or redistribute derivative-map terms whose momentum orders
are not uniform near `H=0`. No health inference is made from any of
those polynomials, nor from their apparent positive or complex roots.

The limitation is mathematical, not merely cautious terminology. For
the exact free oscillator `x''+Kx=0`, set `q=(1+Ku²)x`. Its exact
transformed equation is

```
q''-2(f'/f)q'+[K+2(f'/f)²-f''/f]q=0, f=1+Ku².
```

At the center its frozen squared frequency is `-K` although the original
equation has `+K`. The two equations have the same exact solutions under
an invertible map. A derivative map with a nonuniform large-K time jet
invalidates the proposed inference. The replay verifies this control and
the actual Cauchy-map jet, and retains the full physical maps for future
analysis. A uniformly weighted first-order or pseudodifferential scalar
reduction remains a substantive research task; it is not a user blocker.

## 6. Complete vector action, sign, clock and inner limit

Independently use the transverse unit-Jacobian spatial pullback

```
h/a² = [[1,0,0],[0,1,E_z],[0,E_z,1+E_z²]].
```

Direct matrix contraction gives
`tr(Kij²)-(tr Kij)²=-6H²+(E_z'-shift_z)²/(2N²)`.
The positive relative spatial root has trace
`y sqrt(4+(Delta E_z)²)` and determinant `y²`; literal elementary
polynomials give potential `-N a³ mu(Delta E_z)²/4`, where
`mu=2y[b1+b2(c+y)+b3cy]=yP` in this action. These calculations fix the
transverse-component normalization without importing a TT or paper
vector prefactor.

For each of the two transverse components the action is

```
L=A K(Eg'-sg)²+B K(Ef'-sf)²+Cs(sf-sg)²-Dv K(Ef-Eg)²,
A=a³/4, B=b³/(4c), Dv=a³mu/4.
```

Literal two-shift elimination yields

```
Kv=AB Cs K/[AB K+(A+B)Cs],
omega_bare²=Dv K/Kv=m_alg²+cV² K/a²,
m_alg²=mu(1+c/y³), cV²=(c+y)/(2y).
```

On `2<c<=4`, `Cs>0` and hence `Kv>0` for every `K>0` on the local
interval. Since `c>y`, the full-parent vector principal cone is wider
than the actual matter cone. For every fixed c, the high-K limit of
`Kv` is the smooth nonzero `Cs`, so its normalization pump and time
derivatives are bounded on compact time intervals and do not change
this vector principal coefficient. This uniformity argument is not
silently transferred to the scalar map.

For the `c=1` control `Cs<0`; at the center
`Kv=-128K/[3(K-192)]`, negative for `K>192`. This is a formal
high-momentum wrong-sign kinetic result. The finite-K shift singularity
at `K=192` is not ignored or claimed regular. Neither this result nor
the wider positive-branch vector cone is by itself a light-only,
below-cutoff matching exclusion.

On the positive-kinetic branch let `V=sqrt(2Kv)(Ef-Eg)` and
`theta=Kv'/(2Kv)`. Up to the explicit boundary `-(theta V²)'/2`, the
canonical action is `(V'²-[omega_bare²-theta'-theta²]V²)/2`.
At `c=4,u=0` it has

```
Kv=8K/[3(K+16)], omega_bare²=3K/2+24,
theta'+theta²=-(47K-96)/[3(K+16)],
omega_can²=(9K²+382K+2112)/[6(K+16)]>0.
```

This last positive instantaneous coefficient is not a rolling spectral
gap theorem. The normalization correction is explicitly retained.

For `c=2+epsilon²`, `u=epsilon x`, fixed `K>0` and compact x intervals,
rewrite `Kv=AB K/[A+B+AB K/Cs]`. Its denominator has a smooth positive
extension at `epsilon=0`, where `Kv=K/5`. The canonical normalization
pump is finite and its scaled contribution vanishes; the exact scaled
equation tends to `V_xx+80V/(1+8x²)=0`. Uniform coefficient convergence
on each such compact x set gives ordinary first-order Cauchy-propagator
convergence for convergent initial data. This is only the shrinking
physical window, not a fixed-CD-window, large-momentum, adiabatic,
retarded-source or rolling-gap theorem.

## 7. Verification and remaining obligations

The source-hashed replay pins and replays S6.20. It checks symbolic
identities, the exact continuous Bernstein margin, independently
polarized Fraction/Taylor phase matrices and literal vector-shift
fixtures, and the separately authored action/constraint audit. It does
not equate matching formulas to independent evidence; actual output
bridges compare all retained phase coefficients at the advertised
fixtures. The real-analysis continuity arguments are written here, not
proof-assistant formalized.

Nothing changes the physical matter frame or the free-chi action, and
no scalar degree of freedom is dropped. What remains unresolved is the
uniform scalar kinetic/gradient/characteristic problem and its physical
finite-band control. Nonlinear health, source matching, cutoff, loops,
UV positivity and an original C/D operator realization are not inferred.
The separate variable-lapse G1 profiles have different second jets and
are not covered. The completed scoped photon objective and frozen
linear classification are unchanged; original P8 remains open.
