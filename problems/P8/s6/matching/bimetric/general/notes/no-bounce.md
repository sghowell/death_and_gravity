# Proof: regular common-flat HR cannot realize a nondegenerate CD bounce

The [formulation](../FORMULATION.md) fixes the action, matter frame, regular
positive-root branch and matter hypotheses. All parameters `beta0,...,beta4`
are free real constants. This is a background theorem: it needs no
assumption about the stability of a corresponding vacuum.

## 1. Full two-lapse equations

Write `G=M_g^2>0`, `F=M_f^2>0`, `m4=mathfrak m^4>0`,
`y=b/a>0`, `c=N_f/N_g>0` and `H_i=adot_i/(N_i a_i)`,
`D_i=N_i^{-1}d/dt`. The square-root generating polynomial is
`(1+zc)(1+zy)^3`. After the exact Einstein time boundaries, the full
coordinate-volume Lagrangian is

\[
\begin{split}
L={}&-\frac{3Ga\dot a^2}{N_g}-\frac{3Fb\dot b^2}{N_f}
 +\frac{a^3K_g}{2N_g}-N_ga^3V_g
 +\frac{b^3K_f}{2N_f}-N_fb^3V_f\\
&-m4\{\beta_0N_ga^3+\beta_1(N_fa^3+3N_ga^2b)
 +3\beta_2(N_fa^2b+N_gab^2)\\
&\hspace{35mm}+\beta_3(3N_fab^2+N_gb^3)+\beta_4N_fb^3\}.
\end{split}\tag{1}
\]

For canonical matter `K_i` denotes a coordinate-time positive kinetic
quadratic form in that sector's own fields. It is independent of the
lapse in the variation. The primary `g`-only case sets `K_f=V_f=0`.
For two independent sectors,

\[
\rho_i=\frac{K_i}{2N_i^2}+V_i,\quad
p_i=\frac{K_i}{2N_i^2}-V_i,\quad
\rho_i+p_i=\frac{K_i}{N_i^2}\ge0. \tag{2}
\]

The theorem only uses the corresponding classical NEC and separate
conservation, so (2) is a constructive sufficient matter class rather
than a required choice of potential. No physical field is shared between
the two actions.

The interaction density and pressure read

\[
\begin{split}
\rho_{gv}&=m4(\beta_0+3\beta_1y+3\beta_2y^2+\beta_3y^3),\\
p_{gv}&=-m4[\beta_0+\beta_1(2y+c)+\beta_2(y^2+2cy)+\beta_3cy^2],\\
\rho_{fv}&=m4(\beta_4+3\beta_3/y+3\beta_2/y^2+\beta_1/y^3),\\
p_{fv}&=-m4[\beta_4+\beta_3(2/y+1/c)
                    +\beta_2(y^{-2}+2/(cy))+\beta_1/(cy^2)].
\end{split}\tag{3}
\]

Literal Euler variation of `N_g,N_f,a,b`, before choosing cosmic time,
gives all four equations

\[
\begin{split}
3GH_g^2&=\rho_g+\rho_{gv},&
G(2D_gH_g+3H_g^2)&=-p_g-p_{gv},\\
3FH_f^2&=\rho_f+\rho_{fv},&
F(2D_fH_f+3H_f^2)&=-p_f-p_{fv}.
\end{split}\tag{4}
\]

Their null combinations are

\[
-2G D_gH_g=\rho_g+p_g+\rho_{gv}+p_{gv},\qquad
-2F D_fH_f=\rho_f+p_f+\rho_{fv}+p_{fv}. \tag{5}
\]

Fixing either lapse before varying or retaining only its Friedmann
equation can lose one of (4). The code and independent covariant tests
retain both lapse equations and both accelerations.

## 2. Null-stress cancellation and undivided Bianchi identities

Define the real polynomial

\[
P(y)=m4(\beta_1+2\beta_2y+\beta_3y^2).
\]

Direct factorization of (3) gives

\[
\rho_{gv}+p_{gv}=(y-c)P(y),\qquad
\rho_{fv}+p_{fv}=\frac{c-y}{cy^3}P(y). \tag{6}
\]

Their positively weighted sum vanishes. This identity is independent of
`beta0,beta4` and does not require `P` to be nonzero.

The matter conservation equations in each metric and the metric Noether
identities imply separate interaction balances. With
`B=N_g bdot-N_f adot` they are

\[
\begin{split}
D_g\rho_{gv}+3H_g(\rho_{gv}+p_{gv})
 &=\frac{3P(y)}{N_g^2a}B=0,\\
D_f\rho_{fv}+3H_f(\rho_{fv}+p_{fv})
 &=-\frac{3P(y)}{N_gN_fa\,c\,y^3}B=0.
\end{split}\tag{7}
\]

All denominators in (7) are nonzero on the stipulated regular branch.
The resulting constraint is `P(y) B=0`, not a quotient by `P` or by a
Hubble parameter. Independent Fraction coefficient algebra checks both
forms of (7) and all of (3)--(6).

## 3. Exhaustive stationary-point argument

Let `t0` be a physical stationary slice, `H_g(t0)=0`. Set `N_g=1` only
after deriving (4); overdots below are physical `g` cosmic-time
derivatives. Smoothness implies continuity of `y(t)` and `P(y(t))`.

**Case I: `P(y(t0))!=0`.** There is an open neighborhood of `t0` in
which `P` stays nonzero. Equation (7) gives `bdot=c adot` there, hence
`H_f=H_g/y` without division by `H_g`. Differentiating only on this
neighborhood gives

\[
D_fH_f=\frac{\dot H_g}{cy}
        -\frac{H_g\dot y}{cy^2}
        =\frac{\dot H_g}{cy}\quad\hbox{at }t_0. \tag{8}
\]

This is the origin of the lapse/ratio weights; replacing it by
`D_fH_f=Hdot_g` would be incorrect. Substitute (6) and (8) into (5),
multiply the `f` equation by `cy^3`, and add:

\[
\boxed{-2(G+Fy^2)\dot H_g
   =(\rho_g+p_g)+cy^3(\rho_f+p_f)\ge0.} \tag{9}
\]

The inertia coefficient and matter weight are positive. Thus
`Hdot_g(t0)<=0`. No `f` lapse constraint, condition `y=1`, choice of
vacuum or assumption on the sign of `P` is needed.

**Case II: `P(y(t0))=0`.** At this one point the `g` interaction null
stress in (6) vanishes. The `g` equation in (5) immediately gives

\[
\boxed{-2G\dot H_g=\rho_g+p_g\ge0.} \tag{10}
\]

This argument leaves the `f` Hubble parameter, acceleration, lapse and
matter arbitrary. In particular it does not infer `H_f=0` or use (8).
It therefore applies to a simple root, a double root, an isolated crossing,
a branch-switching point, an algebraic-branch interval, or the identically
zero polynomial `beta1=beta2=beta3=0` with arbitrary cosmological terms
`beta0,beta4`. No analyticity, fixed branch on a whole interval, or bound
on the number of switching events is assumed.

The cases are exhaustive at every stationary slice and imply the same
pointwise sign. They exclude all nondegenerate physical bounces in the
stipulated common-flat sector. The proof does not purport to classify every
possible degenerate extremality with `Hdot_g=0`.

For the primary `g`-only matter case, (9) reduces to
`-2(G+Fy^2)Hdot_g=rho_g+p_g>=0`. The separate-sector extension is stronger
than merely adding an arbitrary hidden fluid: its own NEC and conservation
are required. A shared doubly coupled field generally has a different
source/Bianchi structure and is outside this proof.

## 4. CD and robust physical-time mismatch

For `a_CD=(1+(t/tau)^2)^2`,

\[
H_{CD}=\frac{4t}{\tau^2+t^2},\qquad
\dot H_{CD}(0)=\frac4{\tau^2}>0.
\]

Either (9) or (10) contradicts an exact CD branch. The obstruction does
not depend on a derivative expansion, a mass gap, or the value of
`M*tau`.

The inherited mismatch is also stable against shifting a proposed
stationary slice. On `I=[-tau/2,tau/2]` in physical `g` cosmic time,

\[
H_{CD}(\pm\tau/2)=\pm\frac8{5\tau},\qquad
\dot H_{CD}\ge\frac{48}{25\tau^2}.
\]

For `w=(t/tau)^2∈[0,1/4]`, the exact lower-bound remainder is

\[
\tau^2\dot H_{CD}-\frac{48}{25}
 =\frac{(1-4w)(52+12w)}{25(1+w)^2}\ge0.
\]

Suppose a parent obeyed endpoint Hubble errors strictly below
`8/(5tau)` and uniform derivative error strictly below
`48/(25tau^2)` throughout `I`, while remaining within the stated regular
domain. Its Hubble parameter would have opposite endpoint signs and
positive derivative throughout `I`. The intermediate-value theorem would
give a zero with positive derivative, forbidden by the exhaustive cases.
Hence those two errors cannot both hold. The comparison uses the prescribed
physical metric and its cosmic time, not a redefined metric or clock.

## 5. Independent arithmetic and adversarial controls

The independent `polynomial.py` has no SymPy dependency. It multiplies
the four elementary factors of `(1+zc)(1+zy)^3`, retains independent
couplings and metric jets, and varies the resulting Laurent polynomial
with exact `Fraction` coefficients. Monomials in `a,b,N_g,N_f` alone
may have negative powers. The four normalized Euler expressions have
respectively 7, 7, 11 and 11 terms. Thirteen coefficient identities verify
the full equations, interaction stress, weighted cancellation and both
Bianchi relations. A separate bridge compares the full polynomial action
and all four equations to the symbolic implementation. This is not
evaluation on a finite sample or an interpolation argument.

The checks also include:

- A positive-ratio stationary proportional configuration `y=c=2`, with
  `beta0=-6,beta1=1,beta2=beta3=0,beta4=-1/8` and zero matter. It obeys
  all four equations, demonstrating that the old beta1-specific conclusion
  `y=1` is not a general premise. No spectrum is asserted for this control.
- A jet with `G=F=m4=a=b=N_g=1`, `c=2`, `adot=bdot=0`,
  `addot=1/10,bddot=1/5`, the old beta1 coefficients, and
  `K_g=4/5,V_g=-2/5,K_f=V_f=0`. Both lapse equations and the `g`
  acceleration hold with nonnegative matter null stress. The omitted `f`
  acceleration has exact error `3/5`.
- A forbidden `g` NEC jet and, separately, a forbidden `f` NEC jet,
  each satisfying all four equations algebraically. Their null stresses
  are respectively `-2/5` and `-1/2`. They are not claimed to extend
  to healthy parent solutions; they test the hypothesis boundary.
- On a constant-ratio algebraic chart, `H_f=H_g/c` at `y=1`, so
  `D_fH_f=Hdot_g/c^2` differs from the inapplicable dynamical expression
  `Hdot_g/c`. At `c=2,Hdot_g=1` the difference is `-1/4`.
  The root proof never makes this inference.
- Setting a lapse to zero makes the metric determinant zero. Formal
  continuation to `c=-1` makes `cy^3` negative and can turn a weighted
  pair of nonnegative matter null stresses into a negative sum. This is
  solely a check that the positivity step needs its hypotheses: the
  original `sqrt|f|` positive-root action is not obtained by analytically
  continuing its positive-lapse formulas to a negative lapse.

These controls do not supply physical counterexamples outside the theorem.
The coefficient engine rejects rounded, nonfinite and boolean coefficient
inputs, disallows inverse non-monomial polynomials, and does not silently
evaluate a negative scale/lapse power at zero.

## 6. Precise conclusion

This closes the common-flat, regular, arbitrary-parameter HR background
screen under the specified separate classical NEC matter assumptions.
The immutable beta1 vacuum/source calculation remains a separate fact.
Nothing here establishes healthy arbitrary-beta spectra, interaction or
quantum bounds, nonlinear curvature-squared matching, or an alternative
parent construction.

Non-flat spatial curvature, non-common/non-bidiagonal geometries, singular
or different square-root branches, quantum or NEC-violating stress,
nonminimal/derivative/shared matter couplings, and changed physical-frame
prescriptions are outside this theorem. The known interaction NEC
anti-correlation must not be confused with a statement that physical
matter violates the NEC. The source audit describes the supporting
literature without claiming global novelty. General S6 matching and
both tracks of P8 remain separate research obligations.
