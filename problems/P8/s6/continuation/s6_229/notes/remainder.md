# Actual strong remainder uniformly on each bounded external-momentum ball

Fix Pmax<infinity. Internal loop momenta still range to infinity. For the pair ellvec=-kvec+Pvec, split the radial integral into a bounded region and k>R, where R exceeds2Pmax and a fixed background/mass threshold. The split is an estimate of the original integral, not a new regulator, state, or physical EFT cutoff.

On k>R, both k and |ellvec| are uniformly large. The original canonical frequencies, pumps, constrained two-form vertices and normalized mode expansions are smooth rational functions of these momenta and background jets. For each specified finite number of time derivatives, the original all-order prepared state admits a sufficiently deep finite WKB/Riccati comparison, as in S228. The variation-of-constants error is arbitrarily inverse-power small at physical dimension3. Uniformity follows from |ellvec|>=k/2 and |P|<=Pmax. No derivative in P of the selected cutoff-summed state is used.

An essential finite-transfer correction to the phase is

omega_k+omega_ell =2k/a-(n dot P)/a+O(k^-1).

Consequently the leading oscillation has residual factor exp(-i(n dot P)Delta_sigma). This is orderzero in k, NOT an inverse-k-small error. Retain it as a smooth time-triangle amplitude. It equals1 on the diagonal and all its fixed time and simultaneous-time derivatives are bounded uniformly on the fixed ball. Its difference from1 has a lag factor tau=t-s. Angular integration preserves these uniform bounds.

The leading two-form matrix is the all-dimensional matrix in notes/matching.md. For physicald3 its geometric factor relative to the flat proper-time leading singularity is

tau^5/[a(t)^4 a(s)Delta_sigma^5],

which equals1 on the diagonal. The residual finite-transfer phase, this ratio, and all subleading polarization, mass and gradient terms therefore leave a difference of causal singular order at most4 after the leading order-five distribution is matched. This statement keeps the finite-transfer massive reference difference; the chosen reference is q0 with the ORIGINAL mass, and its finite-P corrections are not omitted.

For every fixed simultaneous derivative Dt+Ds, differentiating the phase inserts k[sigma'(t)-sigma'(s)], which contains a lag factor compensating the added radial power. The orderzero finite-P phase has bounded simultaneous derivatives on the ball. Deepen the original finite comparison as needed for the chosen derivative count. The absolutely integrable error remains uniformly integrable; no bound uniform in derivative order is asserted.

In the original analytic-dimensional neighborhood, the unchanged W4 plus dimension-independent Borel correction has mixing O(k^-6). The leading k^(4+Re(d)-3) factor leaves an undifferentiated error O(k^(-2+|Re(d)-3|)), locally uniformly integrable. The original dimensional finite part may therefore be taken before the physicald3 higher-regularity argument. The all-dimensional invariant match and fixed physical contractions exclude an extra unfixed fourth-order contact. Evanescent jets, original finite curvature action, normalization and second-vertex contacts are retained. In particular the finite constants -4 and -1/30 are unchanged.

On the bounded internal region, the original massive full covariance and constrained feature estimates are bounded, including k0 and ell0. Individual polarization directions need not have a continuous coordinate chart there; the complete physical sums, with bounded projectors and positive mass, give an integrable bound. All fixed time derivatives are uniformly bounded on this region. Thus it gives an ordinary smooth-time contribution, without needing angular or external-momentum derivatives of the state.

After subtracting the matched fourth-order reference, the remaining fixed causal finite parts have off-diagonal powers tau^-n with n<=4, logarithms, ordinary remainders, and local derivatives of order at most3. Four OUTPUT primitives of tau^-n have the form constant times tau^(4-n)log(tau) plus ordinary polynomials. For a variable smooth amplitude, expand about the diagonal to the necessary finite depth; the residual vanishing order makes the remaining integral ordinary. The same reasoning applies to each fixed simultaneous derivative. Every original finite extension remains in this operation. Preparation removes only the original lower zero germ, not a final source interval or an upper retarded contact.

The local integration formula, for j<n, is

I^n[c_j(t)D^j h](t)
 =integral ds sum_(r=0)^j (-1)^r binom(j,r)
   (t-s)^(n-1-j+r)c_j^(r)(s)/(n-1-j+r)! h(s).

For j=n there is additionally c_n(t)h(t), and the ordinary kernel has the same sum with r=1,...,n. It follows by n-fold prepared integration by parts, keeping coefficient order and the upper delta contribution. Apply it with n4 to the first three rows, and n2 to the matter row. The variable eta pivot therefore contributes its multiplication A(t)=-6delta(t)^2 AND its full lower kernel; it is not commuted through I4.

Combining the complete nonlocal scalar response, both ordered Ward legs, the distinct clock contact, all classical and retuning terms, and the matter invariant yields the actual normal form

mathcal I T_adapt = Fref+V,
mathcal I=diag(I4,I4,I4,I2),
Fref=diag(A(t), gamma L0^T diag(Ftrace,(8/3)F2)L0, -1).

Here Ftrace,F2 act by their original q0 proper-time convolution. For every fixed j, the COMPLETE row-sum matrix multiplier kernel satisfies

ess sup_(|P|<=Pmax) max_i sum_l
 |(Dt+Ds)^j V_il(t,s;P)| <= Cj(Pmax)(1+|log(t-s)|)

with finite actual constants Cj. These constants include the full local coefficient and kappa factors and are not evaluated numerically here. The estimates are uniform down to small nonzero external momentum because the complete Ward expression has no inverse P, the direction projectors are bounded and the high-k split was uniform on a ball. The literal P0 homogeneous constraint sector remains a separate mode problem; it is a measure-zero set in the admitted Fourier L2 space.

This strong normal form is proved directly from the actual mode expansion and fixed singular extension. It is not inferred from the old derivative-losing weak response estimate.
