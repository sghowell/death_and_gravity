# Arbitrary-order short-distance control for the specified matrix state

This is a written homogeneous-mode proof on the fixed compact reference slab, not a formalization by a finite Riccati test. It establishes the oriented two-cone short-distance condition for the reduced physical free reference. It does not identify a covariant gauge-fixed Green-hyperbolic operator or its renormalized stress.

## Actual symbol hypotheses and the two chart transitions

The current S241 coefficient map is retained in full. Its exact bounds give F>1/100, Jc<100 and Jc-F>1/1000. Therefore c_s^2=F/Jc lies strictly between0 and1, c_s>1/100, and1-c_s^2>1e-5. Since a<=25/16, the positive roots after factoring k=|P| are omega_s=c_s/a and omega_t=1/a. They satisfy omega_s>4/625 and omega_t-omega_s>8/(25*10^5). The sum of any two positive roots is bounded away from zero. These statements include the tiny FULL existing profile; it is not set to zero.

All exact reference functions are smooth. The Proca and H profiles are the fixed renormalized means of their already constructed Hadamard states; their smoothness, not merely five numerically bounded derivatives, is the needed input. On any fixed compact subinterval every specified finite derivative is bounded. No analytic-in-order or uniform mass bound is required.

After the exact symmetric-boundary shear in state.md, the outer coefficients K,A_B,V are smooth and independent of q except the explicit qG. Centrally, the complete rational coefficients have K=K0+O(q^-1), A_B,V=O(1), with nonvanishing Delta at q>=4096. Expanding a rational function with denominator bounded away from zero at q=infinity gives an asymptotic symbol to ANY finite depth. Taylor's remainder formula bounds each fixed t and k derivative, including q=k^2/a^2. Retain all its coefficients and all full A,Tcorr terms.

In either principal pair choose R as recorded in uv.py. Direct full algebra gives R^T K0 R=I and R^T G R=diag(c_s^2,1). For example R_outer=[[Theta/sqrt(2Jc),0],[-w/sqrt(2Jc),1]]. Its central analogue has E in the first entry and -ell E in the second. Nonzero pivots are needed only in their certified charts.

The symplectic change Q=R a^-3/2 x, p=R^-T a^3/2 pi has a smooth time connection. The whole canonical equations in this frame take

x'=A(t,k)x+B(t,k)pi,
pi'=-C(t,k)x-A(t,k)^T pi,

where B0=I, C2=diag(omega_s^2,omega_t^2), B and C are symmetric, and A,B,C have full classical inverse-k expansions. The change's time derivative is part of A and the lower C, not an omitted term.

The cleaned central/outer map, scaled by diag(sqrt(k)I,k^-1/2 I), is a classical orderzero symbol with bounded inverse and principal diag(-E/Theta,1,-Theta/E,1). This is checked on the ENTIRE finite-q map, not assumed from the principal pencils. Both overlaps are compact and keep E,Theta nonzero. Hence it preserves the bounded high-frequency phase norms and the sign of the leading symplectic polarization. The two exact switches introduce neither a new state nor a source.

The original canonical matrix is smooth at every finite k and gamma crossing. The known S241 tame estimate is polynomial in k times a finite constant on this slab. Conjugation by the displayed polynomial-order maps preserves such a bound. Differentiating its ODE any fixed number of times in k or t adds at most a fixed polynomial loss, by Duhamel and the polynomial/rational symbol coefficients. Thus arbitrarily deep asymptotic error can absorb every fixed derivative loss. A uniform-in-order estimate or convergent formal series is not used.

## Positive Lagrangian recursion and exact comparison modes

Write the positive-frequency Lagrangian graph pi=Z x, with Z^T=Z and negative imaginary part. It obeys

Z'+C+A^T Z+ZA+ZBZ=0.

Seek Z~k Z1+Z0+k^-1 Z_-1+..., with Z1=-i diag(omega_s,omega_t). At each successive order the unknown symmetric matrix X appears ONLY through -i(Omega_f X+X Omega_f), Omega_f=diag(omega_s,omega_t). Its inverse divides each entry by the POSITIVE SUM omega_i+omega_j. All remaining terms involve already determined coefficients, their derivatives and the full lower A,B,C coefficients. This gives a unique smooth formal negative-frequency graph to every finite depth. No same-branch eigenvalue difference is needed for this graph recursion.

After any finite truncation, choose k large enough that -Im Z_N>0. At a fixed reference time in the preparation chart, set x_N=[2(-Im Z_N)]^-1/2 and pi_N=Z_N x_N. The resulting full complex frame F_N has i F_N^dagger Omega F_N=I and F_N^T Omega F_N=0 EXACTLY, not just asymptotically. Evolve these Cauchy data with the full original canonical equation. The exact frame remains normalized and defines a legitimate positive pure comparison state.

The formal transport within the negative graph is a two-by-two first-order symbol with principal -ik diag(omega_s,omega_t). Its off-diagonal amplitude transport can be diagonalized recursively using the strictly nonzero DIFFERENCE omega_s-omega_t, retaining every connection and lower coefficient. Diagonal terms are integrated into smooth symbol amplitudes multiplying phases exp[-ik integral omega_i dt]. This is a finite-order construction with a residual at any prescribed inverse power; derivatives of the amplitudes and residual remain symbols. It is a local computation in one time variable: at each order the off-diagonal algebraic equation has a uniform nonzero gap, and the diagonal transport is an ordinary smooth scalar ODE.

Normalize the finite construction to the exact initial frame above. Variation of constants and the already established polynomial tame bound make the difference between the exact and finite construction smaller than any prescribed k^-R, with any chosen finite collection of derivatives, by FIRST choosing its depth sufficiently large. This is valid across the finite chart cover. On an overlap the exact symplectic transition sends the negative graph to a negative graph and satisfies the same full transformed equation. The leading sign and uniqueness of the recursive coefficients identify the two formal graphs; the remaining mismatch is as small as the chosen depth. This supplies matching through both switches without an instantaneous reinitialization.

The two tensor oscillators are the one-frequency version: positive principal root1/a, smooth full coefficients, and the same positive-sum recursion. The degeneracy between their roots is harmless because the rotation sectors are EXACTLY block diagonal. This argument does not infer a general coupled equal-root theorem.

## Smooth sampling controls the anomalous matrix, not just diagonal modes

Express the ACTUAL selection form in a sufficiently deep exact normalized comparison frame F_N. Its number-conserving Hermitian block is the integral of F_N^dagger E_t F_N, and its anomalous symmetric block is the integral of F_N^T E_t F_N. All cross entries are retained. The leading normalized amplitudes and strict positive K,G imply that the number block has eigenvalues between c k and C k on the preparation interval, for positive fixed c,C and sufficiently large k. The nonzero smooth bump preserves this bound after integration. These constants are finite and need not be small.

Every finite symbolic term of the anomalous integral has the form integral f(t)^2 b_ij(t,k) exp[-ik integral(omega_i+omega_j) dt] dt. Its phase derivative has absolute value at least2k omega_min. Repeated integration by parts with the full smooth f has NO boundary term. Every fixed derivative of its symbol amplitude and inverse phase derivative is bounded at the required order. Each integration gains one inverse power. Momentum derivatives can add powers of bounded time intervals and derivatives of symbol amplitudes, and time derivatives of later exact propagation cost only fixed powers. Deepen the finite construction before fixing the desired decay. The exact-comparison error is controlled by the preceding tame bound. Thus the WHOLE anomalous block is O(k^-R) for every R and every chosen finite derivative family.

This argument does not require the number-conserving matrix to be diagonal. A realified positive Hermitian number block commutes with the standard Omega, and its Gaussian ground covariance is exactly I/2 regardless of its mixing. In the same normalized comparison basis write the real preparation matrix as M_0+Delta, where M_0 commutes with Omega, c k I<=M_0<=C k I, and Delta is the anomalous part. The explicit ground formula is smooth on SPD matrices and homogeneous of degree0. Divide both matrices by k: the segment between them stays in a fixed compact positive cone for large k, so the derivative of the formula is bounded. Repeated derivatives use the inverse, positive determinant square root and positive scalar frequency sum; none divides by a frequency difference. Consequently V(M_0+Delta)-I/2 is rapidly decreasing with the same arbitrary-order control.

For one tensor mode the formula sqrt(det M) M^-1/2 has the same properties. Hence the ACTUAL selected scalar and tensor covariances differ from arbitrarily deep comparison covariances by an arbitrarily prescribed inverse power. No finite-N comparison state is identified with the actual state, and no finite cutoff is part of its definition.

## Infrared existence, distributions and oriented frequency condition

On every bounded momentum set the original complete canonical ODE, preparation Gramian and its positive ground formula are smooth in P^2, with finite bounds and strict positivity. In particular the preparation norm's positive mass-scale term makes both tensor and scalar covariances bounded at P0. The three-dimensional measure d^3P is therefore locally integrable. The bounded TT projector may depend on direction there, but is still integrable; no singular global polarization frame is introduced. Physical shift-vector reconstruction has at worst the existing1/|P| amplitude, whose squared size remains locally integrable in three dimensions. This does not authorize an infrared-unrestricted scalar shift-potential field, whose stronger inverse Laplacian would need a separate test-function domain.

At high momentum, the exact covariances and propagators have polynomial bounds and the preceding asymptotics. Their full Fourier kernels therefore define operator-valued distributions. Positivity and CCR follow from the finite-mode positive covariance integrated on arbitrary complex compact test vectors; this is a statement about the algebraic quasifree state and does not assume a single globally unitarily equivalent Fock representation at all times.

Fix the transform convention fhat(tau,xi)=integral exp[-i(tau t+xi.x)] f(t,x) dt dx. An annihilation phase exp[-ik integral omega dt+iP.x] has first-leg time covector tau=-omega |xi|. Let C_ann be the union of the two nonzero shells tau=-|xi|/a and tau=-c_s|xi|/a. The selected field two-point functions obey WF(W) contained in C_ann times(-C_ann), with the associated component polarizations. This explicit sign convention avoids importing a label for future frequency from a different Fourier convention.

To see the inclusion directly, localize away from these shells and integrate the finite positive-phase Fourier parametrices by parts. Their amplitudes are symbols and their phase gradients are separated from the excluded directions. The exact remainder can be made smoother than any prescribed order by increasing the construction depth first. Thus it contributes no wavefront direction outside the displayed set. This establishes the claim without assuming that a finite-order WKB truncation is Hadamard. The exact CCR identifies the singular part with the corresponding frequency half of the actual reduced commutator. Tensor modes contribute only their speed1 shell.

Any TWO states on this SAME reduced free algebra satisfying this oriented condition differ smoothly: their two-point difference equals its transpose by CCR, while the wavefront sets of the difference and its transpose lie in opposite, disjoint products of cones. Their intersection is empty because every root has nonzero time sign. Local pseudodifferential chart reconstructions with the displayed symbols preserve this inclusion; they do not convert it into a proof of a covariant gauge-fixed field algebra or a local renormalized stress prescription.

## Precise literature boundary

[Fewster, Definition5.2 and Theorem5.4](https://arxiv.org/html/2503.12537) supply a generalized multiple-cone framework for decomposable Green-hyperbolic operators. Its existence clause presupposes an existing Hadamard state; it does not construct the current state. The proof above uses the explicit homogeneous mode system, not an unverified assignment of Green-hyperbolic differential-operator hypotheses to a reduced constraint chart.

[Fewster and Strohmaier, Theorems2/4](https://arxiv.org/html/2510.11492v2) emphasize exact positivity and CCR in constructing states from parametrices. Their normally hyperbolic principal symbol is not the present two-speed matrix. Neither that result nor the scalar Olbermann theorem is imported as a coupled-state theorem here.

Normal ordering relative to this state on its reduced free algebra is now meaningful. A full curved covariant counterfunctional, gauge/constraint measure, renormalized physical metric/light means, quantitative loop bounds and the interacting same-state bounce remain separate calculations.
