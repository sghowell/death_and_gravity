# Exact large-asymmetry limit and fixed-interval response

## 1. Scope, scales and state

Use the fixed action, physical composite metric and parameter family in
[FORMULATION.md](../FORMULATION.md). Work in physical \(u=mT\); a prime
below always denotes \(d/du\). The rate variables and squared frequencies
are divided by \(m\) and \(m^2\) respectively. The Planck coefficient
\(M\) is retained only when reconstructing an actual metric perturbation.

This gate concerns a singular large-asymmetry limit of regular
finite-parameter backgrounds. It does not regard the limiting individual
\(g\) metric as a regular parent solution. It proves exact limiting
equations and convergence on every fixed compact forward interval.
Nothing in the proof is uniform on an interval whose length grows with
\(\varepsilon^{-1}\).

For clarity first let the fixed initial scaled density and rate be
\(\eta_0>0,a_0>0\), with \(y_0=\varepsilon^{-1}\),
\(\bar\rho_0=\eta_0\varepsilon^2\), and
\[
D_\varepsilon=1+a_0\varepsilon u^2/4,\quad
A_e=D_\varepsilon^2,\quad
h=\frac{a_0\varepsilon u}{D_\varepsilon},\quad
(m\tau)^2=\frac4{a_0\varepsilon}.
\]
The full evolution/growth theorem then fixes \(\eta_0=9,a_0=11\).
The centre-only coefficient theorem includes \(a_0=0\) as an algebraic
endpoint; it is not a finite-duration CD history.

## 2. Full scaled ODE, with no division by a vanishing background Hubble rate

Starting from the pinned pressure-branch reconstruction, set
\[
e=1/y,\quad \eta=\bar\rho/e^2,\quad j=h/e,\quad
\mathfrak a=h'/e,\quad
\bar X=\sqrt{1+\eta e(1+e)^3/3},\quad
\bar R=\sqrt{1+\eta(1+e)^3/3}.
\]
On the regular positive branch \(X=\bar X/e\), \(Y_f=-e\bar R\).
Here \(Y_f\) denotes the reconstruction's negative \(f\)-Hubble root,
not the initial ratio \(Y=\varepsilon^{-1}\).
Define
\[
\lambda=\frac{y'}y
=\frac{(1+e)[-\bar X(\bar R+j)+je^2\bar R]}{\bar X+e\bar R},
\qquad
\kappa=ec=\frac{\bar X-je(1+e)}{\bar R+j(1+e)}.
\]
Then the full source-aware equations are exactly
\[
e'=-e\lambda,\qquad
\eta'=-3j\left[\eta e+\frac2{(1+e)^2}\right]+2\eta\lambda,\qquad
j'=\mathfrak a+j\lambda.
\]
The matter identity follows by differentiating \(\eta=\bar\rho/e^2\);
holding the scaled density \(\eta\) constant during time differentiation
would be incorrect.

Use \(e=\varepsilon z\) to remove the parameter-dependent initial scale:
\[
j=\frac{a_0u}{D_\varepsilon z},\qquad
\mathfrak a=
\frac{a_0(1-a_0\varepsilon u^2/4)}{D_\varepsilon^2z},\qquad
z'=-z\lambda,\quad z(0)=1,\quad\eta(0)=\eta_0.
\]
This is an exact two-dimensional nonautonomous ODE for \(z,\eta\).
Its right-hand side is analytic wherever its explicit radicands and
denominators remain positive. The individual lapses in physical time are
\[
N_g=\frac e{e+\kappa},\qquad
N_f=\frac{\kappa}{e+\kappa},\qquad N_g+N_f=1.
\]
The positive-root identities are verified against the immutable
reconstruction, not against a quotient of Hubble rates at the bounce.

The canonical null source is
\[
\bar\rho+\bar p=e\left[\eta e+\frac2{(1+e)^2}\right].
\]
Its bracket tends to \(2\). Thus positivity of the kinetic scalar does
not require the scaled energy density \(\eta\) itself to stay positive.
When the bracket is positive and the lapse/root conditions hold, a local
canonical scalar and reconstructed potential follow from
\(\chi'=M\sqrt{\bar\rho+\bar p}\) and
\(V=M^2m^2(\bar\rho-\bar p)/2\).
Each finite member has an invertible scalar clock on the interval;
uniform clock/potential-derivative bounds at the endpoint are not claimed.

## 3. The limiting flow is integrable

At \(\varepsilon=0\) define
\[
R=\sqrt{1+\eta/3}>0,\qquad v=R+j=1/\kappa.
\]
The full scaled equations imply
\[
\lambda=-v,\quad R'=1-Rv,\quad
j'=\mathfrak a-jv,\quad
v'=1+\mathfrak a-v^2,\quad z'=vz.
\]
The prescribed CD history also gives
\(\mathfrak a'=h''/e+\mathfrak a\lambda\).
On each fixed interval \(h''/e\to0\), exactly as checked from
\(D_\varepsilon\). Hence
\[
\mathfrak a'=-\mathfrak av,\qquad
\mathfrak az=a_0,\qquad z''=z+a_0.
\]
With \(v_0=\sqrt{1+\eta_0/3}\),
\[
z=(1+a_0)\cosh u+v_0\sinh u-a_0,\qquad
v=z'/z,\quad\mathfrak a=a_0/z.
\]
Moreover \((jz)'=a_0\), so
\[
j=a_0u/z,\qquad R=(z'-a_0u)/z,\qquad\eta=3(R^2-1).
\]
These reconstruct both limiting density and lapse ratio, not just an
isolated equation for a frequency.

For \(\eta_0=9,a_0=11,v_0=2\), set
\[
z=12\cosh u+2\sinh u-11.
\]
On \(u\geq0\), \(z\geq1\), \(z'>0\), and the first integral is
\[
z'^2=z^2+22z-19,\qquad v^2=1+\frac{22}{z}-\frac{19}{z^2}.
\]
Since \(z\geq1\), \(v^2\geq1\). The exact square identity
\[
\frac{140}{19}-v^2=\frac{(11z-19)^2}{19z^2}\geq0
\]
therefore proves \(1\leq v\leq\sqrt{140/19}<3\).
The numerator of \(R\) obeys
\[
(z'-11u)'=z\geq1,\qquad z'(0)=2,
\]
so \(R>0\) on every forward interval. The limiting radicands and lapse
ratio are regular in scaled variables there.

These are not all-real-time bounds. At
\(u=-\operatorname{atanh}(1/6)\), \(z'=0\) and
\(z=\sqrt{140}-11<1\), so the forward proof cannot be reused backwards.
Also
\[
R-1=\frac{11(1-u)-10e^{-u}}z,
\]
which is negative at \(u=1\). Thus the limiting density becomes negative;
this is allowed by the formulation and does not make the scaled null
source negative.

## 4. Complete canonical coefficient limit

The S6.9 weighted/relative map is retained:
\[
l=f_\Sigma(w_1h_g+w_2h_f),\qquad
H=f_R(h_f-h_g),\qquad
w_1=\frac{\kappa e^2}{1+\kappa e^2},\quad w_2=\frac1{1+\kappa e^2}.
\]
The physical-time kinetic normalizations are
\[
f_\Sigma^2=\frac{M^2A_e^3}{4}
\frac{(e+\kappa)(1+\kappa e^2)}{\kappa(1+e)^3},
\qquad
f_R^2=\frac{M^2A_e^3e^2}{4}
\frac{e+\kappa}{(1+e)^3(1+\kappa e^2)}.
\]
Although \(f_R\to0\), its logarithmic derivatives are regular after
separating \(e'/e=-\lambda\). Write the two displayed rational factors
as \(S_\Sigma\) and \(S_R\), so
\[
\theta_\Sigma=\tfrac32h+\tfrac12(\log S_\Sigma)',\qquad
\theta_R=\tfrac32h-\lambda+\tfrac12(\log S_R)',\qquad
N_i=\theta_i'+\theta_i^2.
\]
The exact moving-field coefficient, algebraic stiffness and physical
source projection are
\[
\omega=\frac{e\sqrt\kappa}{1+\kappa e^2}
\left(\lambda-\frac{\kappa'}{2\kappa}\right),
\]
\[
\bar m_{\rm alg}^2=
\frac{(1-e)(1-\kappa)(1+\kappa e^2)}
{(1+e)(e+\kappa)^2},
\qquad
R_{\rm src}=\frac{e\kappa-1}{(1+e)\sqrt\kappa}.
\]
All canonical time boundaries and the phase definition
\(p_l=l'-2\omega H,\ p_H=H'\) remain those of S6.9.

In the limit \(h=0,\kappa=1/v,\lambda=-v\), direct differentiation
using the actual \(v,\mathfrak a\) flow gives
\[
\theta_\Sigma=N_\Sigma=\omega=0,\quad
\theta_R=\frac{3v^2-1-\mathfrak a}{2v},
\]
\[
N_R=\frac{3(\mathfrak a+1)^2+3v^4-2v^2}{4v^2},\qquad
\bar m_{\rm alg}^2=v^2-v .
\]
Here \(N_R\) is already in \(m^2\) units. In particular, differentiating
the normalization only after freezing the initial ratios would lose
the first term in the time evolution.

At a general centre with \(v>1,\mathfrak a\geq0\), the exact leading
hierarchy margin is
\[
N_R-\tfrac34\bar m_{\rm alg}^2
=\frac{3(\mathfrak a+1)^2+v^2(3v-2)}{4v^2}>0.
\]
This rules out making that particular ratio small within the fixed
large-asymmetry scaling. It does not rule out every nonadiabatic
light-only reduction. The often useful tuning
\(\mathfrak a=3v^2-1\) sets \(\theta_R=0\) but leaves
\[
N_R=(15v^2-1)/2,\qquad
(\bar m_{\rm alg}^2)'=2v^2(2v-1).
\]
For the chosen initial state these are \(59/2\) and \(24\); the second
mass-squared derivative is \(-34\), and \(\bar m_{\rm alg}^2=2\).

The smaller centre orders are also checked:
\[
\theta_\Sigma=\frac e2(1+\mathfrak a-3v)+O(e^2),\quad
N_\Sigma=\frac e2(v-3)+O(e^2),\quad
\omega=-\frac e{\sqrt v}\theta_R+O(e^2),\quad
\omega'=-\frac e{\sqrt v}N_R+O(e^2).
\]
They do not cancel the \(O(1)\) relative normalization curvature.

## 5. Actual time-dependent limiting homogeneous response

The exact \(k=0\) charge inherited from the full action is
\[
J=f_\Sigma(p_l-\theta_\Sigma l),\qquad J'=0.
\]
At \(J=0\), the full finite-parameter equations are
\[
H''+( \bar m_{\rm alg}^2-N_R)H=0,\qquad
l'=\theta_\Sigma l+2\omega H.
\]
These are not equations obtained by freezing a diagonal entry. For the
normalized data
\[
H(0)=1,\quad H'(0)=0,\quad l(0)=0,\quad J=0
\]
the phase condition is \(p_l(0)=0\). At finite parameter it need not
mean \(l'(0)=0\). In the limit \(\theta_\Sigma=\omega=0\), so \(l=0\).

Let \(D_H=N_R-\bar m_{\rm alg}^2\) in the exact limiting trajectory:
\[
D_H=\frac{3(\mathfrak a+1)^2}{4v^2}
    +\frac{-v^2+4v-2}{4}.
\]
Because \(1\leq v<3\),
\[
(-v^2+4v-2)-1=(v-1)(3-v)\geq0,
\]
and therefore \(D_H\geq1/4\) at every forward time.

The equation \(H''=D_HH\) with the stated data has \(H>0\): a putative
first zero is impossible because until that time \(H''\geq0\) and
\(H'\geq0\). Set \(W=H-\cosh(u/2)\). It obeys
\[
W''-\tfrac14W=(D_H-\tfrac14)H,\qquad W(0)=W'(0)=0.
\]
The retarded kernel is \(2\sinh[(u-s)/2]\), positive for \(u>s\).
The exact integral formula thus gives
\[
W(u)=\int_0^u2\sinh[(u-s)/2](D_H(s)-\tfrac14)H(s)\,ds\geq0 .
\]
Consequently \(H(u)\geq\cosh(u/2)\) for the actual full limiting
time-dependent equation. With zero heavy initial data instead, the
homogeneous solution is identically zero. The result does not assert
that a source or every initial state excites this column.

## 6. Physical projection and the linear-amplitude guard

Linearizing the composite spatial metric gives
\[
\gamma_{\rm eff}=\frac{h_g+yh_f}{1+y}
=\frac{l+R_{\rm src}H}{f_\Sigma}.
\]
The exact expression for \(R_{\rm src}\) above follows from
\((w_1-1/(1+y))/\sqrt{w_1w_2}\). In the limit
\[
f_\Sigma=M/2,\qquad R_{\rm src}=-\sqrt v,\qquad
\gamma_{\rm eff}=-\frac2M\sqrt v\,H.
\]
Thus the normalized relative response is not projected away. In the
fixed common-field gauge,
\[
|\gamma_{\rm eff}(u)|\geq \frac2M\cosh(u/2).
\]
A common constant homogeneous TT coordinate shift cannot remove the
time dependence. In particular the difference between its value at
\(u\) and at zero satisfies
\[
|\gamma_{\rm eff}(u)-\gamma_{\rm eff}(0)|
\geq\frac2M\max\{0,\cosh(u/2)-\sqrt2\}.
\]
For a local curvature check, \(v(0)=2,v'(0)=8,v''(0)=-54\) and
\(H''(0)=55/2\). Therefore
\[
\gamma_{\rm eff}'(0)=-\frac{4\sqrt2}{M},\qquad
\gamma_{\rm eff}''(0)=-\frac{20\sqrt2}{M}.
\]
The physical limiting composite background is Minkowski at fixed \(u\).
The nonzero second derivative gives a nonzero linearized physical tidal
curvature; it is not a removable constant tensor coordinate shift.
Cosmic-time derivatives restore one or two powers of \(m\).

There is an essential amplitude qualification. The unit canonical
column \(H(0)=1\) has
\[
\delta=h_f-h_g=H/f_R=O(\varepsilon^{-1})
\]
on a fixed interval. It is a normalized solution of the linear equations,
not a uniform small-metric-perturbation family. At each finite
\(\varepsilon\), multiply the entire linear solution by an independent
amplitude, for example \(\varepsilon\epsilon_{\rm lin}\), with
\(\epsilon_{\rm lin}\) small enough for that fixed interval. The
individual metric perturbations then stay small, while the absolute
composite perturbation is \(O(\varepsilon\epsilon_{\rm lin})\) and tends
to zero. The normalized growth ratios and nonzero normalized tidal
derivatives survive linear rescaling. No unsuppressed absolute response
for fixed small individual metrics, finite nonlinear instability, or
uniform nonlinear approximation follows.

## 7. Fixed-interval finite-parameter convergence

Fix any finite \(L\geq0\) before taking \(\varepsilon\to0^+\).
The explicit limiting trajectory on \([0,L]\) has \(z\geq1\),
\(R>0\), and \(v\geq1\), with finite upper bounds on every continuous
coefficient. In particular the limiting radicands are \(1\) and \(R^2\);
the limiting numerator/denominator of \(\kappa\) are \(1\) and \(v\).
A compact tube around this trajectory can therefore be chosen with
strict positive margins for all these quantities and for \(z\).

The exact two-dimensional scaled vector field is analytic on this tube
for sufficiently small \(\varepsilon\). It and its required derivatives
have finite uniform bounds there. Let \(K_L\) bound its state derivative
and \(B_L\) bound its parameter derivative. With equal initial \(z,\eta\),
the standard integral estimate, proved directly by Gronwall, is
\[
\sup_{0\leq u\leq L}|(z_\varepsilon,\eta_\varepsilon)
 -(z_0,\eta_0)|
\leq B_L\varepsilon\frac{e^{K_LL}-1}{K_L},
\]
with the continuous interpretation \(B_L\varepsilon L\) if \(K_L=0\).
Choose \(\varepsilon\) so this bound is smaller than half the tube
radius. A first-exit argument then guarantees existence throughout
\([0,L]\) and closes the estimate. No numerical values of these compact
derivative bounds or of the required parameter threshold are asserted.

All regularized canonical coefficient formulas, including their
required derivatives, are analytic functions of the state and
parameter on the same tube. They converge uniformly at rate
\(O_L(\varepsilon)\). The factors \(N_g\) and \(a_g\) are positive for
every sufficiently small positive parameter but tend to zero; their
inverses are not included in a claimed uniform regular-geometry norm.
Similarly
\[
(\bar\rho+\bar p)/e\longrightarrow2
\]
uniformly, so the finite-parameter canonical scalar retains positive
kinetic energy even after the density changes sign.

Apply variation of constants to the first-order system for
\((H,H',l)\) at \(J=0\). Its coefficients are precisely the regularized
ones just controlled, and its canonical initial data are fixed.
Finite compact-interval propagator bounds give
\[
\sup_{[0,L]}\bigl|(H_\varepsilon,H_\varepsilon',l_\varepsilon)
 -(H_0,H_0',0)\bigr|=O_L(\varepsilon).
\]
The composite source coefficient and \(f_\Sigma^{-1}\) also converge
uniformly, so the normalized physical composite projection and its
required time derivatives converge on the same fixed interval.
This proves more than pointwise centre jets, but it supplies no numerical
finite-\(Y\) error threshold.

The argument does not interchange \(L\to\infty\) with
\(\varepsilon\to0\). The physical CD duration in \(u\) grows like
\(\varepsilon^{-1/2}\); even a logarithmically growing interval requires
new quantitative estimates. The limiting growth theorem is therefore
not relabeled as an all-time finite-\(Y\) instability or a controlled
finite-band EFT exclusion.

## 8. Verification boundary

The source-hashed replay checks the full scaled source/clock algebra,
the exact integrable trajectory, derivative/normalization identities,
positive polynomial decompositions, and physical projection. A separate
Fraction engine checks the coefficient identities and exact initial
fixtures. Independent covariant/ODE and amplitude audits protect the
source, clock, inverse map and time-dependent comparison.

The positivity, comparison and compact-interval convergence arguments
are written analytic proofs supported by those exact identities, not
Lean-formalized analytic theorems. No expanded numerical propagator
majorant, finite-\(k\) spectrum, scalar/vector health, strong-coupling
scale, quantum state, original-M1-frame matching, or general UV verdict
is established. Other initial data, other HR coefficients and possible
nonadiabatic light-only reductions remain research options.
