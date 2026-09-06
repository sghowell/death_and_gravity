# P8-S6.10.COMPOSITE — large-asymmetry coefficient limit and homogeneous response

This child pins S6.9's `../certificates/composite-light.json`, SHA256
`36c899d2b82957f2fc83df8e7bc66582cbf77829fbed8b64724d2d39c2086e48`,
and replays its lineage without modification.

## Scope and model

Keep the beta2 composite-matter action with constant
\[
G=F=M^2,\quad m_4=M^2m^2,\quad\alpha=\beta=1,\quad
(\beta_0,\ldots,\beta_4)=(0,0,1,0,0).
\]
The physical metric is the composite metric in the pinned \(+---\)
convention. It is not the original CD/M1 matter frame. The positive-root
common-flat pressure branch reconstructs a canonical scalar potential
for each finite member of the family; that potential need not be
positive and is not a single prescribed free M1 action.

Use physical \(u=mT\), \(A_e(0)=1\), and
\[
\varepsilon=Y^{-1}>0,\quad y_0=Y,\quad
\bar\rho_0=9/Y^2,\quad
h(u)=\frac{11\varepsilon u}{1+11\varepsilon u^2/4},\quad
m\tau=\frac2{\sqrt{11\varepsilon}}.
\]
Bars denote division of densities by \(M^2m^2\), and \(h=H_e/m\).
The exact scaled variables are
\[
e=1/y=\varepsilon z,\quad \eta=\bar\rho/e^2,\quad
j=h/e,\quad\mathfrak a=h'/e .
\]
The symbol \(\mathfrak a\) is a rate variable, not a metric scale factor.

## Claims

1. The full scaled CD reconstruction ODE and regularized TT canonical
   coefficient functions extend analytically to \(\varepsilon=0\)
   in a neighbourhood of the explicit limiting trajectory on every
   fixed compact forward interval \([0,L]\).
2. The limiting flow is exactly
   \[
   z=12\cosh u+2\sinh u-11,\quad
   v=z'/z,\quad\mathfrak a=11/z,\quad j=11u/z,\quad
   R=(z'-11u)/z>0,\quad\eta=3(R^2-1).
   \]
   For \(u\geq0\), \(z\geq1\) and
   \(1\leq v\leq\sqrt{140/19}<3\).
3. The complete time-dependent limiting relative normalization gives
   \[
   \bar N_R=\frac{3(\mathfrak a+1)^2+3v^4-2v^2}{4v^2},\qquad
   \bar m_{\rm alg}^2=v^2-v,
   \]
   with bars on these coefficients denoting division by \(m^2\).
   More generally at a centre with fixed \(\eta_0>0,\mathfrak a_0\geq0\),
   the same leading formulas imply
   \(\bar N_R>\tfrac34\bar m_{\rm alg}^2\). This is a fixed-scaling
   mass-led hierarchy screen, not a general matching exclusion.
4. In the exact limiting fixed-charge \(k=0\) sector, take
   \(J=0,l(0)=0,H(0)=1,H'(0)=0\). Then \(l=0\) and
   \[
   H''=(\bar N_R-\bar m_{\rm alg}^2)H,\qquad
   \bar N_R-\bar m_{\rm alg}^2\geq\tfrac14,\qquad
   H(u)\geq\cosh(u/2).
   \]
   This uses the actual full time-dependent limiting equation.
5. The physical composite tensor projection has the finite relative
   coefficient \(R_{\rm src}\to-\sqrt v\). It retains the limiting
   response and has a nonzero linearized physical tidal derivative.
   With matched canonical/charge data the finite-\(Y\) coefficients,
   solutions and composite projection converge on each fixed \([0,L]\).
   The convergence constants exist by the written compact-domain proof;
   no numerical \(Y_*(L)\) or growing-window error estimate is supplied.

## Essential boundaries

Only forward intervals are covered by the explicit uniform bounds.
The limiting density \(\eta\) may become negative; the physical canonical
null source remains positive for sufficiently large finite \(Y\).
The individual \(g\) lapse and scale vanish at \(\varepsilon=0\), and the
scalar clock also degenerates in that limit. The limit is a regularized
coefficient problem, not a new regular two-metric parent solution.

The heavy response is a specified linear homogeneous \(k=0\) solution.
It is not inevitable for zero heavy data, not a finite-\(k\) dispersion
or scalar/vector result, and not a finite-\(Y\) all-time instability
theorem. Neither a CD-duration window \(u=O(\varepsilon^{-1/2})\) nor a
logarithmically growing window is covered. The limits \(Y\to\infty\)
and \(L\to\infty\) are not interchanged.

Unit canonical mode columns may be rescaled to keep the original metric
perturbations small at finite \(Y\); no uniform finite-amplitude
nonlinear perturbation claim follows from their normalization. No
ghost/cutoff/UV verdict, exclusion of all nonadiabatic light-only
reductions, original-frame matching, or completion of S6/P8 is claimed.
