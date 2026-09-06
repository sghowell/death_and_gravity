# S6.3.beta1: a specified two-metric parent screen

This gate tests one new parent family, not a modification of any frozen P8
action. It establishes a positive quadratic Minkowski spectrum and a
source-preserving flat quadratic tree-level Weyl coefficient. The same
parent cannot support the CD bounce on its stipulated regular FLRW branch.
Thus this is an exact failed matching screen, not completion of S6 or P8.

## Frozen physical frame and exact parent

Use signature `+---`, `R=-6(Hdot+2H^2)` on spatially flat physical FLRW, and

\[
S=-\frac{M_g^2}{2}\int\sqrt{|g|}R_g
  -\frac{M_f^2}{2}\int\sqrt{|f|}R_f
  -\mathfrak m^4\int\sqrt{|g|}\sum_{n=0}^4\beta_n e_n(\sqrt{g^{-1}f})
  +S_{\rm can}[g].
\]

`M_g^2,M_f^2,nu=mathfrak m^4 beta1` are strictly positive and
`beta0=-3 beta1`, `beta4=-beta1`, `beta2=beta3=0`. The positive square-root
branch is used. An explicit matter choice is a canonical clock with
`V(phi)=m_phi^2(phi-phi_v)^2/2`, `m_phi^2>0`, and free canonical M1 matter
`chi`, both coupled only to `g`. The vacuum has constant fields, `phi=phi_v`,
zero vacuum energy, and `f=g=eta`. The background obstruction also holds
for any classical positive field-metric canonical `g` matter with arbitrary
smooth potential: it uses only `rho+p>=0`.

The physical matter metric remains `g`. The old DHOST action is not appended
to its Einstein action, and no derivative metric redefinition is used to
identify a different geometry as physical. A flat vacuum spectrum does not
construct a rolling clock, interpolate to CD, or bound interactions.

## Exact claims

For a common flat spatial foliation, require positive finite scale factors
`a,b`, positive finite lapses `N_g,N_f`, and smoothness through the putative
bounce. All four scale/lapse equations are derived before choosing physical
cosmic time `N_g=1`. The undivided Bianchi equation is

\[
N_g\dot b-N_f\dot a=0.
\]

At a stationary physical slice, it implies `H_f=0`; the full `f` lapse
equation uniquely fixes `y=b/a=1`. The two acceleration equations give

\[
-2(M_g^2+M_f^2)\dot H_g=\rho+p\geq0.
\]

Consequently a nondegenerate CD bounce, requiring `Hdot_g(0)=4/tau^2`, is
impossible in this parent family. This is independent of a mass-gap or
derivative expansion. A quantitative finite-window mismatch is proved in
[the proof](notes/parent-screen.md): moving the stationary slice does not
allow arbitrarily small physical `H,Hdot` errors on `[-tau/2,tau/2]`.

At the constant Minkowski vacuum the two spin-2 kinetic coefficients and
all physical on-pole residues are positive, with

\[
M^2=M_g^2+M_f^2,\qquad
m_{\rm FP}^2=\nu\left(M_g^{-2}+M_f^{-2}\right)>0.
\]

Keeping the prescribed physical source, both true massive-mode integration
and the metric Schur calculation give the same response. Including the
conserved-source scalar/trace projector, not just TT waves, fixes the
flat quadratic four-derivative pure-metric representative, modulo boundaries,
to `c_C C^2+c_R R^2` with

\[
c_C=\frac{M_f^4}{4\nu}>0,\qquad c_R=0
\quad\text{(flat quadratic tree order only).}
\]

The coefficient `c_R=0` is not nonlinear, cosmological, or quantum matching.
In particular it cannot remove the isolated-matter R² counterterm of S5.9.
The flat rational-kernel remainder is controlled in its stated spectral
disk; that is not a rolling-background, spatial-momentum, or strong-coupling
cutoff estimate.

## Evidence and exclusions

The new read-only report pins and replays S6.2 and S5.11, including their
replayable lineage. Exact algebra and independent covariant/lapse/source
tests support the [written proof](notes/parent-screen.md) and
[source audit](notes/sources.md). It reports a failed parent match even
though the separate flat quadratic vacuum test passes.

This gate does not exclude other parameter families, non-flat or non-common
spatial foliations, singular/nonpositive square-root branches, noncanonical
or nonminimal matter, quantum NEC violation, other physical-frame proposals,
or general UV completions. It neither resums a curvature truncation to infer
ghosts nor claims parent nonlinear/cosmological health from TT residues.
There is no common parent CD solution on which to transfer the S5.10/S5.11
finite-band claims. Those certificates remain unchanged and conditional on
their own candidate, data, frame, and matching-coefficient assumptions.
