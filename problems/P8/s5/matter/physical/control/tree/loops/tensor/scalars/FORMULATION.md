# P8-S5.11.CD: constant-Weyl candidate, coupled scalar response

## Claim and immutable input

This is a new gate for the **same** candidate

\[
S_{\rm CD/M1}+c_C\int\sqrt{|g|}\,C_{abcd}C^{abcd},\qquad c_C=\text{constant},
\]

in the prescribed physical matter metric. It pins and replays S5.10.CD,
certificate SHA256
`05a66f42c280e099d73ef4629a0358cf84d1eba4af5eb692dc4ff83ecbd739d8`.
No old action, formulation or certificate is changed. Flat FLRW is an exact
background of the candidate because its Weyl tensor and the full first
variation of the added functional vanish, not because variation is restricted
to FLRW in advance.

The gate establishes the scalar quadratic Weyl functional; regular first-order
stationary elimination of the old lapse/shift; its derivative phase map and
physical lapse/shift/spatial-metric reconstruction; and a coupled, rank-one
quadratic Hamiltonian correction. It gives finite-band bounds for the chosen
first-order two-scalar Hamiltonian and for a chosen first-order phase pullback.
These are not exact higher-derivative branch approximation or quantum causality
theorems.

## Variables and normalization

Signature is `+---`. The physical ADM scalars in clock gauge and scalar spatial
gauge are `n`, `zeta_phys`, and `b=a^2*beta_shift`, with
`g_0i=-a^2 partial_i beta_shift`. The old auxiliary spatial map gives
`zeta_phys=v_old+delta*n`, `delta=1/(2*(1+u^2)^3)`. Neither `v_old` nor the gamma
canonical coordinate is silently identified with the corrected physical metric
or shift.

Write the physical matter perturbation as `delta chi=M*s`. The old normalized
phase density is `p*v_dot+P_s*s_dot-h0` with overall `M^2*a^3`. In cosmic units
`beta=c_C/M^2`. In fixed centre units
`sigma=(t-t0)/ell0`, `ell0=tau*sqrt(1+u0^2)`, the corresponding coefficient is
`epsilon0=c_C/(M*ell0)^2`. Scalar momentum densities include the expansion
term `-3H*p`, not a constant-volume canonical convention.

The scalar profile is a real spatial Fourier mode of squared spatial norm one.
For a literal cosine the spatial average is one half; restoring unit norm
multiplies it by two. Both coupled scalar species remain. The TT sector was
treated separately in S5.10 with its stated unit polarization convention.

## Exact first-order result

Let `q=k_com^2/a^2`, `w=l*(3delta-1)=-l*Lambda` on CD, and use the pinned
`h0`, `J>0`, `theta`, and `Lambda=1-3delta`. The physical scalar lensing variable
is

\[
\mathcal L=(1-\delta)n-v_{\rm old}+\dot b-Hb,
\quad \Delta S_S^{(2)}=\frac{4c_C}{3}\int dt\,a^3q^2\mathcal L^2.
\]

After the old partial Legendre transform the auxiliary Hessian is
`diag(-2J,-2q^2/3)`. Thus the phase construction is regular for every finite
time and `q>0`, including the gamma crossing. It does not divide by `theta`
or `Lambda`. The old solutions are `n0=-R/J`, `b0=-p/(2q)`. Define

\[
\mathcal L_0=2\delta n_0-2(v+Hb_0),\qquad
E_p=\dot p+3Hp-2q(v+\Lambda n_0).
\]

The full off-shell map and reduced Hamiltonian are

\[
v_{\rm old}=v+\frac\beta3(E_p-4q\mathcal L_0)+O(\beta^2),\qquad
h_{\rm red}=h_0-\frac{4\beta}{3}q^2\mathcal L_0^2.
\]

The `E_p` term can be dropped only on the perturbative leading branch.
Physical reconstruction on that branch is

\[
\begin{split}
v_{\rm old}&=v-\tfrac43\beta q\mathcal L_0+O(\beta^2),\\
n_{\rm phys}&=n_0-\tfrac83\beta\delta q^2\mathcal L_0/J+O(\beta^2),\\
b_{\rm phys}&=b_0+4\beta D_t\mathcal L_0+O(\beta^2),\\
\zeta_{\rm phys}&=v_{\rm old}+\delta n_{\rm phys}.
\end{split}
\]

`D_t` includes all CD coefficient derivatives and the full old coupled phase
flow. The pre-mixed auxiliary metric momentum is `p+3l*s`; the normalized
physical matter momentum density is `l+P_s+3l*v_old`. These boundary terms are
retained, and `a^3*l` is constant.

## Quantitative operational contract

On every centre window `abs(t-t0)<=ell0/100`, set `a(t0)=1` and restrict
`10^11<=k_com*ell0<=4*10^12`. Use one of the pinned overlapping scalar charts
throughout each window. The full S5.7 whitening, gyroscopic connection and
symmetric momentum boundary are retained, with positive free energy norm
`N_t(X)^2=P^T P+Y^T W(t)Y` in **fixed dimensionless time units**.

If an independent matching assumption supplies
`epsilon=abs(c_C)/(M*tau)^2<=10^-116`, then the chosen reduced evolution with
equal reduced-canonical centre data differs from the old evolution by at most
`10^115*epsilon` in relative old free-energy norm. Its first Born/Taylor
remainder is at most `10^229*epsilon^2`. The chosen reduced momentum Hessian
is positive, with eigenvalues greater than `1/2`, on this domain only.

The branch phase map in the same norm differs from identity by at most
`10^84*epsilon` and has a Neumann inverse. Treating that first-order map as an
exact **chosen pullback**, the same-pulled-old-phase-data evolution bound is
`10^116*epsilon`. Equal reduced data, equal pulled old `(v,s,p,P_s)` data and
equal physical metric/lapse/shift data are different prescriptions. No
physical-metric same-data error theorem is asserted by relabelling them.

For independently supplied `abs(c_C)<=1` and S5.8's `M*tau=10^324`, the reduced
evolution bound is `10^-533`, its Born remainder bound is `10^-1067`, and the
phase-map bound is `10^-564`. These deliberately loose numbers certify a
conditional small correction, not a predicted matching coefficient or cutoff.

## Remaining scope

The proof [notes/scalar-response.md](notes/scalar-response.md) distinguishes
formal first-order equivalence, exact evolution of a specified finite-order
representative, and the unsolved higher-derivative branch problem. No bound is
placed on unknown order-beta-squared operators or on an exact extra-mode
initial-value problem. There is no unbounded-momentum/front-velocity, global
fixed-comoving, full nonlinear/vector/loop, corrected quantum-state or complete
candidate health verdict. S5.9 running is not a bound on this finite matching
coefficient; omitted nonlocal, anomaly, state, clock and graviton terms remain
open. P8 is not completed by this gate.
