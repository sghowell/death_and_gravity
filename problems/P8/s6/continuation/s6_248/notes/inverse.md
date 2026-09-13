# Ordered mixed-row inverse and smooth prepared regularity

Fix Pmax<infinity and r real. Work with scalar amplitudes in C(I;H^r_ball)^4, maximum component norm, and its weight exp(-lambda(t+1/2)). The smooth domain is the union over common nonempty initial zero-germ lengths of C^infinity prepared histories in this space. No derivatives of the Fourier multiplier in P or spatial Schwartz-output preservation are required.

S247 proves that the q0 total trace and shear factor inverses K_trace,total and K_2,total are both ordinary half-line L1 kernels. This uses the reciprocal of the COMPLETE two-mass sum, both physical thresholds and the two-sided interior trace cusp. The old isolated Proca pole is not deleted from its own theorem, but is not a pole of this summed reciprocal.

Set Kmax=max(||K_trace,total||L1,(3/8)||K_2,total||L1), and
b0=max(15625/6144,16Kmax/(9gamma),1).
The L1 norms can be taken on the half line; the finite interval needs no larger constants. Both L0 and its inverse, and their transposes, have infinity norm4/3. Thus the two-sided causal reference inverse

B0=diag((-6delta²)^-1,
 gamma^-1 L0^-1 diag(K_trace,total,(3/8)K_2,total)(L0^T)^-1,
 -1)

has weighted continuous-time operator norm at most b0. The first block is variable multiplication, the middle block is convolution, and the last block is multiplication. All physical gamma^-1=64pi²kappa factors remain.

Let C=C0(Pmax) be a majorant for the ENTIRE remainder in remainder.md. For lambda>=1,
integral_0^1 exp(-lambda tau)(1-log(tau))dtau
 <=(2+log(lambda))/lambda<=2/sqrt(lambda).
Consequently ||B0 V||_lambda<=2b0 C/sqrt(lambda).
Choose lambda=(4b0 C+1)². Then
||B0V||_lambda<=2b0 C/(4b0 C+1)<1/2,
also when C0.

The complete ordered inverse is

E_adapt=(I+B0 V)^-1 B0 diag(I4,I4,I4,I2).

For normalized adapted forcing f it satisfies
||E_adapt f||_C<=2exp(lambda)b0 ||diag(I4,I4,I4,I2)f||_C.
The constants are finite but unevaluated and may be enormous. Exponential weighting does not prove a small physical perturbation, a convergent physical Born expansion or stability. No bound uniform in Pmax is asserted.

## Regularity without commuting the variable pivot

The continuous solution preserves the original zero germ by causality. Translate both time arguments by a small epsilon on a shortened interval; the moved lower segment lies in that zero germ. Stationary convolution by both complete reference kernels commutes with these prepared translations. The eta pivot A(t)=-6delta(t)² does not.

Writing V_epsilon(t,s)=V(t+epsilon,s+epsilon), the difference quotient equation is

(Fref+V)Delta_epsilon x
 =Delta_epsilon f0
 -(A(t+epsilon)-A(t))/epsilon P_eta x(t+epsilon)
 -(V_epsilon-V)/epsilon S_epsilon x.

The complete weak-log derivative bounds imply convergence of V's difference quotient in integral-kernel operator norm: use uniform smooth convergence off the diagonal and the integrable uniform bound close to it. The high-radius expansion and integrable actual-state error in remainder.md justify those uniform statements for each fixed derivative. Apply the bounded ordered inverse to obtain the first continuous derivative, including one-sided endpoint limits.

Inductively the n-th derivative obeys
(Fref+V)x^(n)=f0^(n)
 -sum_(j=1)^n binom(n,j)[A^(j)P_eta x^(n-j)+V_j x^(n-j)],
V_j=(Dt+Ds)^j V.
Every fixed smooth seminorm is finite in terms of finitely many C_j, pivot derivatives and source seminorms. There is no assertion of analyticity in n.

The original forward factors preserve this smooth prepared class using S247's paired representation
F_i,total=-c_i delta-D²(G_Proca,i+G_H,i).
The spectral primitives G are bounded on compact time windows, and prepared derivatives move onto the full source without discarding initial distributions. Together with the normal form this gives the actual full forward operator on the same class.

Applying D4,D4,D4,D2 to the normal-form equation yields the original adapted equation. Conversely every smooth prepared solution yields that normal-form equation, so uniqueness follows from the complete reference inverse and contraction. No order reduction, new higher-derivative branch data, final-time condition or reset state is introduced.
