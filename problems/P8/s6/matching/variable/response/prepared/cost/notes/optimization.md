# Exact source-cost characterization and its attainable relaxation

This is a written Hilbert-space and Fourier proof, with exact replay of
the finite algebra and rational enclosures. The optimizer is characterized
by an operator formula; no numerical solution of that formula, numerical
band optimum or proof-assistant formalization is asserted. Companion
notes prove the actual loading and construction constants. This file does
not replace those proofs by arbitrarily supplied Gramian entries.

## 1. The physical loading constraint and the specified norm

Write ell=1/100, a=-1/50, b=-1/100. All time derivatives below are with
respect to the physical g clock u=T/tau. For one unit-normalized real TT
polarization, the physical state is Z=(g,g',f,f'). Its actual equations are

```
Z'=A_delta,K(u) Z+B sigma,   B=(0,2,0,0)^T.
```

The factor two follows from the literal source action and
`K_g=a_g^3/8`: the g equation has right side `a_g^3 sigma/4` before
division by K_g. It is not a source coupled solely to the canonical light
field. The delta=0 coefficients are smooth on the closed punctured J.
Let Phi_0,K be their homogeneous fundamental solution and let
`P0^-(b)` be the full physical endpoint matrix of S6.21. Define

```
h_K(s)=P0^-(b)^-1 Phi_0,K(b,s) B,
C_K sigma=integral_J h_K(s) sigma(s) ds.
```

Equivalently, because the full normalized outer wave matrix propagates
the delta=0 homogeneous equation, the kernel is
`h_K(s)=P0^-(s)^-1 B`. In canonical/outer variables it is

```
h_K(s)=Psi_-(|s|)^-1 (0,j_l(s),0,-sqrt(|s|)j_H(s)/mu)^T.
```

The minus sign is u<0. Both j weights are the physical g-source weights.
This identity does not substitute the physical g observable for the light
canonical coordinate, or delete its heavy projection. The target
`y=Lv`, `L=[I2;0]`, is in this full four-dimensional endpoint frame.
At the formal regular-light delta=0 origin the physical jet is
`(g,g')=(2/sqrt(5))(v_even,v_odd)`; a derivative in proper T has another
factor tau^-1. Thus unit v is not silently called unit physical g.

Set H=L2(J;R). The actual full Gramian is
`G_K=C_K C_K*=integral_J h_K h_K^T`. It is positive definite. For example,
if an adjoint covector is orthogonal to every g-source, its second
physical component p2 vanishes on J. The adjoint equation for p2 makes
p1 vanish; the equation for p1 then makes the fourth component vanish,
because the actual intermetric spring coupling is nonzero; its remaining
equation makes p3 vanish. Hence no nonzero covector annihilates the
range. The companion construction proves surjectivity constructively and
the independent continuous-K calculation gives the quantitative bounds.

The source-free analytic target is not imposed by resetting data at b.
All four data begin at zero at a and are created by the actual source.
Every source optimized below is fixed independently of delta; only the
punctured loading map used to choose it is a limit.

## 2. The exact Fourier operator

Extend sigma by zero to the real line and use

```
sigmahat(omega)=integral_R exp(-i omega u) sigma(u) du,
||sigma||2^2=(1/(2pi)) integral_R |sigmahat(omega)|^2 domega.
```

Plancherel and Fubini on a bounded band give

```
E_< = <sigma,B_Omega sigma>,
B_Omega(s,t)=sin(Omega(s-t))/(pi(s-t)),
E_> = <sigma,R_Omega sigma>,  R_Omega=I-B_Omega.
```

The diagonal kernel is Omega/pi; there is no additional factor of two.
The integral operator is positive and compact, with a continuous kernel, and
`0<=B_Omega<=I`. Its trace is `Omega ell/pi`. For every finite Omega its
operator norm is strictly below one: an eigenvector with eigenvalue one
would have both compact time support and Fourier support in a bounded
interval. The Fourier transform of an L2 function supported on finite J
is entire (differentiate the L1 integral on every compact complex set).
Vanishing on an open real interval would make it identically zero.
Compactness then excludes norm one. Thus R_Omega is strictly coercive.
This argument proves a positive gap, not its numerical value at a large
Omega. For Omega=0, B_Omega=0 and R_Omega=I exactly.

For the calibrated range, positivity and the trace bound instead give the
explicit continuous inequalities

```
||B_Omega||<=Omega ell/pi<=Omega/300<=1/3,
R_Omega >=alpha I,   alpha=1-Omega/300>=2/3.          (1)
```

The rational trace ceiling is strict at positive Omega because pi>3.
These inequalities measure the Fourier energy of the actual zero-extended
physical source. They do not label a compact window times a low-frequency
carrier as bandlimited. In particular, no nonzero admitted source can
have more than one third of its energy in this whole calibrated band.

## 3. Exact optimizer and dual, including active budgets

In this section write C=C_K, y=Lv, R=R_Omega and G=CC*. For any nonzero
target define

```
sigma_min=C* G^-1 y,
s_min^2=y^T G^-1 y.
```

Orthogonal projection onto the affine space `C sigma=y` proves that
sigma_min is its unique minimum-L2-norm member. A budget S<s_min is
infeasible. At S=s_min it is the only L2-feasible source. For S>s_min the
L2 cost problem has a unique optimizer: its feasible set is weakly compact,
the positive quadratic cost is weakly lower semicontinuous, and strict
coercivity gives strict convexity.

For lambda>=0 put

```
R_lambda=R+lambda I,
G_lambda=C R_lambda^-1 C*,
sigma_lambda=R_lambda^-1 C* G_lambda^-1 y.           (2)
```

G_lambda is positive definite. Completing the square on the affine
constraint proves

```
min_(C sigma=y) <sigma,R_lambda sigma>
 =y^T G_lambda^-1 y.
```

Consequently the exact dual for S>s_min is

```
p_Omega(S,y)=sup_(lambda>=0)
             [y^T G_lambda^-1 y-lambda S^2].        (3)
```

There is no unproved minimax exchange: the affine constraint is surjective
and sigma_min is strictly inside the norm ball, so the usual separating
hyperplane/Lagrange multiplier argument applies. Alternatively, (2)
directly constructs the primal optimizer. Choose lambda=0 if
`||sigma_0||<=S`; otherwise choose its unique positive value with
`||sigma_lambda||=S`. Indeed differentiating
`R_lambda sigma_lambda=C* p_lambda` and `C sigma_lambda=y` gives

```
d/dlambda ||sigma_lambda||^2
 =-2 <sigma_lambda',R_lambda sigma_lambda'> <=0.     (4)
```

Its limit as lambda tends to infinity is s_min^2. If the derivative ever
vanishes then sigma_lambda'=0, sigma_lambda belongs to range C*, and is
sigma_min; the stationary condition then makes it the same optimizer for
every lambda. Otherwise (4) is strictly negative. Thus the construction
covers the active-budget case without falsely claiming strict monotonicity
when R=I. At S=s_min the dual value is understood as a possible limit
lambda->infinity; attainment of its multiplier is not automatic.

For a tolerated loading error `||C sigma-y||<=e`, a useful exact dual
under strict feasibility is

```
sup_(lambda>=0,p in R4)
 [2p.y-2e||p||-p^T G_lambda p-lambda S^2].           (5)
```

Without strict feasibility every such pair still gives a valid lower
bound. Equation (5) follows by minimizing the quadratic Lagrangian and
using the support function of the Euclidean error ball. It is not an
instruction to project away the actual heavy state at the endpoint.

## 4. H0^2 infimum versus attained L2 minimum

The H0^2 norm is not bounded in the problem as posed: only the L2 source
budget is bounded. Thus H0^2 is not a closed feasible subset of L2, and it
would be incorrect to assert its minimum is attained. Even the elementary
one-moment problem `integral_J sigma=1`, Omega=0, has the unique L2 optimizer
sigma=1/ell, which does not have zero endpoint traces. At its minimum
budget the H0^2 feasible set is empty.

Let D=C_c^infinity(J). It is dense in H. Since C(D) is a dense linear
subspace of finite-dimensional R4, it equals R4. Hence there is a bounded
finite-rank right inverse R_s:R4->D. The actual H0^2 Hermite right inverse
in the companion proof suffices as well for H0^2 approximation; no
C-infinity claim about that particular Hermite source is required.
For any affine-feasible sigma and approximants d_n in D, the correction

```
tau_n=d_n+R_s(y-C d_n)
```

has exactly the same loading and converges to sigma in L2. It follows
that every S>s_min contains a strictly feasible smooth/H0^2 source.
If the L2 optimizer is at the budget boundary, first mix it with such a
strictly feasible source with weight theta>0, then approximate and correct
as above inside the remaining norm slack. Letting theta and the correction
error go to zero proves

```
inf_(H0^2, C sigma=y, ||sigma||<=S) E_>(sigma)
 =min_(L2, C sigma=y, ||sigma||<=S) E_>(sigma),
                                                     S>s_min.          (6)
```

At S=s_min equality of feasible sets must be checked separately: only
sigma_min is available, and it must itself belong to H0^2. On J the h
kernels are smooth through the endpoints. The optimizer equation
`(1+lambda)sigma_lambda=B_Omega sigma_lambda+h^T p_lambda` makes the L2
optimizer smooth there. It belongs to H0^2 precisely when its values and
first derivatives vanish at both endpoints. No such four trace identities
are presumed for the actual model. Therefore (6) is an infimum statement,
not a hidden claimed H0^2 minimizer.

The approximation can be made constructive. Let a cutoff chi_h lie in
[0,1], equal one away from endpoint strips of total length at most 2h,
and have zero value/first derivative at the endpoints. If the smooth
optimizer has sup norm U, then

```
||chi_h sigma-sigma||2 <= U sqrt(2h),
||chi_h sigma+R_s C[(1-chi_h)sigma]-sigma||2
 <=(1+||R_s|| ||C||) U sqrt(2h).                    (7)
```

Choose h so this is smaller than the available strict norm slack and
desired cost tolerance. Quadratic cost continuity gives
`|E_>(f)-E_>(g)|<=||f-g||2 (||f||2+||g||2)` because ||R||<=1.
For an active budget use the preceding strict-interior mixture first.
The cutoff and correction have finite H0^2 norms, but those derivative
norms need not stay bounded as the optimum is approached. If a fixed H2,
peak or derivative budget is desired, that is an additional optimization
constraint and is not supplied by (6).

No uncontrolled regularity is needed in (7). If `g_- I<=G` and
`||h^(j)||infinity<=H_j`, the optimizer equation and a bounded-band
Cauchy--Schwarz estimate give, for j=0,1,2,

```
||sigma_lambda^(j)||infinity
 <=sqrt[Omega^(2j+1)/(pi(2j+1))] ||sigma_lambda||2/(1+lambda)
   +H_j ||y||/g_-.                                 (8)
```

Here `G_lambda >=G/(1+lambda)` bounds its multiplier. These are available
finite suprema on the actual punctured coefficient interval. A polynomial
cutoff with controlled derivatives and the explicit Hermite right inverse
then makes (7) a finite construction; no band projection is treated as a
causal compactly supported control.

## 5. Cost envelopes and finite Fourier certification

Let `g_- I<=G<=g_+ I`. Equations (1)-(3) imply

```
||y||^2/g_+ <=s_min^2<=||y||^2/g_-,
alpha s_min^2 <=p_Omega(infinity,y)<=s_min^2,
||sigma_0||^2 <=s_min^2/alpha.                       (9)
```

The upper cost uses sigma_min. The optimizer norm estimate follows from
its cost being no larger than that competitor's. For a budget satisfying
`S^2>||y||^2/(alpha g_-)`, lambda=0 is feasible with strict slack, so all
bounds in (9) apply to the actual budget and to the H0^2 infimum (6).

The companion continuous-Gramian bound is `g_-=9/(4*10^12)`,
`g_+=16/25`. Thus the generic full-target norm floor is `5||y||/4`.
The independently derived physical symplectic source kernels for the
regular-light components, `-a_g^3 g_odd/2` and `a_g^3 g_even/2`, improve
the norm floor to `19||v||`, and to `900|v_even|` for the even component.
These are actual source moment bounds, not chosen rescalings of a control.
With unit v, Omega<=100 and S=10^6, (9) gives

```
||sigma_0||^2<=2*10^12/3,
S^2-2*10^12/3=10^12/3>0,
722/3 <p_Omega(S,Lv)<=4*10^12/9,
p_Omega(S,L(1,0))>540000.                           (10)
```

Strict lower signs in (10) come from the companion strict norm floors.
The broad upper bound is not advertised as sharp. No value of the
minimizer or of p_Omega is computed by this arithmetic.

For an explicitly constructed Hermite source with
`||sigma_H||2<R2||y||` and `||sigma_H''||2<D2||y||`, zero-extension H0^2
and Plancherel give

```
E_>(sigma_H)<=min(R2^2,D2^2/Omega^4)||y||^2,
                                                     Omega>0,         (11)
```

and `R2^2||y||^2` at Omega=0. This is an actual feasible upper bound when
its norm budget is allowed, not a claim that it fits S=10^6. The companion
values are R2=2*10^8 and D2=4*10^14. Their coarse derivative tail estimate
does not establish predominantly low-band preparation. A finite peak
bound is also available: `||sigma_H||infinity<=ell^(3/2)||sigma_H''||2/pi`,
by the endpoint and Poincare inequalities. Sharper construction bounds
may be used without changing the source norm definition.

For later improvement of the Fourier part use

```
B_N(s,t)=(Omega/pi) sum_(n=0)^N
             (-1)^n Omega^(2n)(s-t)^(2n)/(2n+1)!.
```

Real sine Taylor's theorem, expanding through the zero even coefficient
of degree 2N+2, bounds the kernel remainder by
`(Omega/pi)(Omega ell)^(2N+2)/(2N+3)!`. The Schur test therefore gives

```
||B_Omega-B_N|| <=(Omega ell/pi)(Omega ell)^(2N+2)/(2N+3)!.
```

This uses no assumption that alternating terms decrease. For N=4 and
Omega<=100 it is at most `1/(3*11!)`. Consequently

```
||C B_Omega C*-C B_N C*||
 <=||C||^2/(3*11!) <=(16/25)/(3*11!).                (12)
```

Using centered moments `m_j=integral_J(s-(a+b)/2)^j h(s) ds`, the finite
matrix `C B_N C*` is the sum of

```
(Omega/pi) (-1)^(n+j) Omega^(2n) binomial(2n,j)/(2n+1)!
 *m_(2n-j) m_j^T,  0<=n<=N, 0<=j<=2n.
```

Thus N=4 needs nine actual moment vectors, degrees zero through eight.
The exact weights and continuous remainder are replayed. Actual band
moments are not computed in this module; they cannot be replaced by the
unweighted Gramian. No numerical band optimum is inferred from (12).

A useful future finite-matrix lower witness needs only `G_< =C B_Omega C*`.
For any covector z set `c=|z.y|`, `a_z^2=z^T G_< z`,
`b_z^2=z^T(G-G_<)z`, `n_z^2=a_z^2+b_z^2`. Orthogonal Fourier decomposition
and Cauchy--Schwarz give

```
c<=a_z sqrt(S^2-E_>)+b_z sqrt(E_>).
```

If c>S n_z the declared preparation is impossible. If c>a_z S and
c<=S n_z, every feasible source obeys the stronger bound

```
E_> >=[ (b_z c-a_z sqrt(n_z^2 S^2-c^2))/n_z^2 ]^2.  (13)
```

The numerator is positive precisely at that threshold. The complementary
root `(a_z c+b_z sqrt(n_z^2 S^2-c^2))/n_z^2` verifies the circle and
constraint identities directly. This is a rigorous source-aware witness,
not a claim that its Cauchy equalities are realizable by a compact source.

## 6. Fixed-source finite-delta propagation and physical units

Let sigma be a single source chosen through C_K, independently of delta,
and let `||C_K sigma-Lv||<=e`, `||sigma||2<=S`. The pinned S6.21/S6.23
actual maps obey

```
||C_delta,K sigma-C_K sigma||<=3*10^6 delta ||sigma||1,
||T_delta||<=42,
||H_delta^side-L||<=200 delta,
T_delta H_delta^-=H_delta^+.
```

Here C_delta,K uses the actual P_delta endpoint map; it is not C_K with
its physical weights frozen. The smooth-source estimates extend to all
L2(J) sources, including the H0^2 controls used here, as follows. Choose
`sigma_n in C_c^infinity(J)` converging to sigma in L2, hence also in L1
because J has finite length. On the fixed punctured preparation interval
the actual first-order coefficients, physical endpoint maps and retarded
kernels are continuous and bounded, uniformly over the stated delta and
K box, including its limiting delta=0 preparation operator. The Duhamel
integral therefore defines a bounded L1-to-endpoint map and
`C_delta,K sigma_n -> C_delta,K sigma` (likewise for C_K). Pass to the
limit in the smooth-source loading estimate with its same constant. The
subsequent source-free transfer is the finite-dimensional bounded map
T_delta, so its estimate and the prepared-sector comparison also pass to
the limit. For L2 sources this defines the usual absolutely continuous
forced solution, smooth after the source has ended; the H0^2 subclass
also has its stated classical endpoint traces. No source derivative enters
any of these kernel estimates, and no source/state is reset in this
extension. Since `||sigma||1<=sqrt(ell) S=S/10`,

```
||H_out,delta-Lv|| <=42e+delta(8600||v||+12600000 S),
||H_out,delta-H_delta^+v||
 <=42e+delta(8400||v||+12600000 S).                  (14)
```

The first combines the fixed-light error and the source-loading error.
The second instead compares the incoming source to H_delta^-v and uses
exact prepared-sector transport; it is not the same target. To report a
physical endpoint vector, multiply both target and error by the same
actual `P_delta^+`, with its separately justified norm bound. A physical
time derivative then includes tau^-1. No canonical-coordinate error is
silently declared a raw g-derivative error.

At S=10^6, ||v||=1, e=0 and delta=10^-17, both right sides are below
1/1000. At delta=10^-9 the coarse bound is not small. For a desired error
eta>42e, (14) instead requires the explicit further restriction
`delta<=(eta-42e)/(8600||v||+12600000 S)` for the Lv target, together with
the frozen delta range. This is a bound, not a manufactured cutoff.

The dimensional source is `Pi(T)=M^2 sigma(T/tau)/tau^2`. Therefore

```
Pihat(omega_T)=(M^2/tau) sigmahat(tau omega_T),
integral |Pi|^2 dT=(M^4/tau^3)||sigma||2^2,
physical out-of-band energy=(M^4/tau^3)E_>,
physical frequency edge=Omega/tau,
integral |Pi|dT=(M^2/tau)||sigma||1.
```

The same multiplier acts on the total and band energies, so their ratio
is invariant. This is a source norm per declared Fourier/polarization
normalization, not the energy density of additional physical matter.
The source has zero temporal components and transverse, traceless spatial
components; conservation is the literal TT probe identity on the same
background. Its amplitude is differentiated at zero before any delta
limit. Finite-norm controls here prove no finite-amplitude backreaction,
quantum-state selection or positivity theorem.

## 7. Boundary of the matching conclusion

This theorem makes source preparation and its temporal spectrum explicit.
It does not make the full four-mode response equal to a two-mode EFT for
arbitrary sources. The cost problem permits signed sources and coherent
cancellations, unlike the near-positive pulse class of S6.27. The
optimizer selects a correlated endpoint state; it does not show that a
generic low-frequency source or a vacuum selects that state.

The unwindowed low-pass projection of a compact source has a noncompact
past and future. Evolving it from newly zero data inside J would change
the problem. Equations (1)-(14) use actual compact controls, not that
replacement. Temporal Fourier concentration is a specified source
diagnostic in proper time, not a stationary spectral theorem for the
rolling dynamics. An RMS derivative scale, a heavy gap, a Wilsonian
cutoff and a nonlinear quantum error remain separate obligations.
No long-history theorem outside the recorded local branch, finite-gravity
positivity result, original C/D matching or original P8 closure is claimed.
