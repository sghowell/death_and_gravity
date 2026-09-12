# Mixed-order prepared inverse and variable-pivot regularity

For a fixed finite Pmax and real r, let H^r_ball be the Fourier-supported subspace of H^r with support in |P|<=Pmax. All real spatial Sobolev norms are equivalent on this fixed ball, with constants that may depend on Pmax and the orders. Use the maximum of the four component C(I;H^r_ball) norms and the same norm with exponential weight exp(-lambda(t+1/2)). Each source has a common nonempty initial zero neighborhood. Equivalently work for each fixed positive germ length, then take their union. No claim of spatial Schwartz-output preservation is needed.

These are the original scalar AMPLITUDE norms, with b the shift divergence. They do not assert that the shift vector beta proportional to P b/P^2 has the same unrestricted infrared H^r norm. First prove the identities on smooth Fourier data supported away from P0, then extend in the amplitude norms by the S219/S222 original weak bound and the new uniform ball estimates. Such data are dense in the ball amplitude space. The inverse-gradient shift is still a well-defined tempered distribution: on a bounded three-dimensional ball, multiplication by1/P sends L2 locally to L1 by Cauchy-Schwarz. No stronger shift-vector norm is claimed.

Let Ktrace and K2 be the original causal inverse kernels. The first is half-line L1; the full pole-plus-cut shear kernel is only finite-window L1, which is all that is used. Define

Kmax=max(||Ktrace||L1(0,1),(3/8)||K2||L1(0,1)),
b0=max(15625/6144,16Kmax/(9gamma),1).

Both L0 and its inverse, as well as their transposes, have infinity norm4/3. Thus

B0=diag(A^-1,
 gamma^-1 L0^-1 diag(Ktrace,(3/8)K2)(L0^T)^-1,
 -1)

has weighted continuous-time operator norm at most b0. It is a two-sided inverse of Fref on prepared causal distributions. The eta block is time-dependent multiplication, not a convolution. All gamma^-1=64pi^2 kappa factors remain.

Let C=C0(Pmax) be a valid COMPLETE row-sum remainder majorant. The integrable weak-log estimate and

integral_0^1 exp(-lambda tau)(1-log tau)dtau
 <=(2+log lambda)/lambda<=2/sqrt(lambda), lambda>=1,

give ||B0 V||_lambda<=2b0 C/sqrt(lambda). Choose

lambda=(4b0 C+1)^2.

Then ||B0V||<1/2, including C0. No physical small-coupling assumption is used. The convergent causal inverse is

E_adapt=(I+B0V)^-1 B0 mathcal I.

The order is essential. On normalized adapted forcing f,

||E_adapt f||_C <=2exp(lambda)b0 ||mathcal I f||_C.

This is a finite bound with unevaluated C and Kmax, possibly enormous, not a small physical inverse or stability estimate. Removing the weight costs exp(lambda); gamma^-1 and the original source derivatives/density normalization are not suppressed. No bound uniform as Pmax tends to infinity is obtained.

## Smoothness with the nonstationary pivot retained

The first construction gives a continuous solution of (Fref+V)x=f0. Causality preserves the exact original zero germ. To prove smoothness, shift both time arguments on a shortened output interval. For sufficiently small epsilon, the changed lower source interval lies in the unchanged zero germ. Convolution by the original Ktrace,K2 commutes with these prepared shifts.

The only nonstationary reference coefficient is A(t) in the eta row. Set A_epsilon(t)=A(t+epsilon), V_epsilon(t,s)=V(t+epsilon,s+epsilon), and let P_eta denote the eta component projection. The difference quotient satisfies

(Fref+V)Delta_epsilon x
 =Delta_epsilon f0
 -(A_epsilon-A)/epsilon P_eta S_epsilon x
 -(V_epsilon-V)/epsilon S_epsilon x.

Apply the bounded ordered inverse of Fref+V. The A difference converges uniformly to A'. Away from the diagonal the V difference converges uniformly in the essential Fourier-multiplier norm; near the diagonal the uniform weak-log derivative majorant gives convergence in the integral-kernel operator norm. The underlying smooth high-frequency coefficients and the uniformly absolutely convergent comparison remainder justify this uniform convergence. This proves a continuous first derivative, including one-sided endpoint limits.

Inductively,

(Fref+V)x^(n)
 =f0^(n)-sum_(j=1)^n binom(n,j)
   [A^(j)P_eta x^(n-j)+V_j x^(n-j)],

V_j=(Dt+Ds)^j V.

Every fixed derivative has a finite bound in terms of finitely many Cj, derivatives of A and forcing norms. No analytic-in-order bound is claimed. Omitting A' would be wrong even though the two middle reference channels are stationary convolutions.

The original forward factors preserve smooth prepared histories: their paired representations are the original local delta terms minus D^2 times a bounded causal spectral primitive. Prepared derivatives can be moved onto the input, retaining all initial distributions. The remainder identity then gives the full forward operator on the same smooth Fourier-supported class. Apply rowwise D4,D4,D4,D2 to the normal-form equation to recover T_adapt x=f. Conversely every smooth prepared solution gives the bounded normal-form equation and is unique by the original inverse identities. There is no order reduction, final-time condition, deleted pole or free choice of new higher-derivative branch data.

Independent finite noncommuting matrix diagnostics check both inverse products, row-primitives, variable reference blocks and RIGHT density multiplication. They detect reversed products and omitted coefficient commutators. These are algebraic controls, not numerical certificates for the continuum kernel.
