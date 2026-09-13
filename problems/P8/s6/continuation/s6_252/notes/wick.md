# Continuum reference Wick noise and a time-retarded extension family

Work on the open slab I=(-1/2,1/2), with compact smooth spacetime tests. The state is exactly FREE-REF-S251, not a cutoff state or an instantaneous vacuum. This note treats specified finite-jet bilinears of its reduced canonical scalar fields and the two TT fields. A physical covariant stress still requires its own complete vertex and subtraction construction.

## Products for the full two-cone state

Every scalar matrix entry of W has first-leg wavefront time component strictly negative and second-leg component strictly positive, in the convention of S251. The first leg lies on one of the two shells tau=-|xi|/a or tau=-c_s|xi|/a. Both positive root magnitudes are bounded below by4/625 after extracting |P|. Tensor entries have the speed1 shell. Differentiation and the stated classical-symbol reconstructions do not introduce a reverse orientation. The low-momentum contributions are smooth: their compact momentum integrals and all coordinate/time derivatives are integrable by the existing infrared bounds. The scalar shift potential remains outside an unrestricted infrared domain.

Consider an entrywise product W_ij(x,y) W_kl(x,y). At any possible wavefront sum the first time component is the sum of two strictly negative numbers. It cannot vanish. The same remains true for arbitrary finite field derivatives and smooth coefficient multipliers. The distribution multiplication criterion therefore defines the whole product on I times R^3 squared, including coincident base points. Conjugation is performed AFTER that ordered product when taking its real or imaginary part. The criterion would not automatically justify an opposite-orientation product W times conjugate(W).

The same reasoning constructs reference-normal-ordered finite Wick polynomials. In a contraction graph for an ordered product of such polynomials, no edge begins and ends at one vertex. Choose the smallest vertex label incident to a nonempty set of contractions. All its incident oriented contractions leave for larger labels, so their time components have the same strict sign and their sum is nonzero. This proves the no-zero-wavefront-sum condition for each contraction graph, retaining arbitrary component mixing and both cones. It is not an assumption that our operator is scalar Klein-Gordon.

One can use smooth approximate-identity operators on both field legs, preserving the common Gaussian state, to construct these products on its finite-polynomial GNS domain. Mixed smoothing parameters have the same microlocal limits. The norm squared of the difference of two regularized Wick vectors is a sum of their contraction distributions and tends to zero. Thus the regularized vectors are Cauchy; finite products have the analogous limit. Positivity passes from these same-state smooth approximations.

For O_A(x)=:z_i A_ij z_j:/2, with finite field derivatives applied to their respective legs, the connected ordered kernel is the full Wick expression of vertices.md. Its real part N is a positive symmetric distribution: for each finite real test combination its pairing is the limiting squared norm of the centered Hermitian operator applied to the reference state. No momentum cutoff or diagonal-value evaluation of an unrenormalized W is needed for this connected kernel. Replacing a fixed bilinear by itself plus a c-number changes its mean but not this connected noise. Changing the operator's actual coupling matrices is a different operation and is not covered by that invariance.

The general multiplication criterion is [Brunetti, Fredenhagen and Koehler, Theorem2.6](https://arxiv.org/pdf/gr-qc/9510056). Their scalar Wick-field theorem is not imported as a theorem for our matrix system: the contraction-graph, sign and smoothing argument above checks the needed product hypotheses directly.

## Complete phase symbol order, including uncleaned momenta

A finite-order negative-frequency comparison frame in either clean chart has position amplitudes of order k^-1/2 and clean momentum amplitudes of order k^1/2, k=|P|. Its entire principal normalization is positive, and the full all-order transport and arbitrarily regular exact-comparison error are supplied by S251. All statements below are uniform on each fixed compact time subinterval and for every specified finite set of derivatives.

The exact full K,G are symbols of q-order0, while EVERY entry of the complete B is of q-order at most1 in both charts. This is checked on the whole rational matrices, not on a selected principal entry. Since q=k^2/a^2, restoring the symmetric-boundary shear gives uncleaned momenta of order at most k^3/2. Outer positions remain order k^-1/2. Centrally, the original position is Pi_b/(2a^3 q) and the original momentum is-2a^3 q Q_b, so the same original-position/order-minus-one-half and original-momentum/order-three-halves bounds hold. No division by Theta is used at crossing. The two TT amplitudes satisfy the same conservative bound; their coordinate-free projector is order0.

The full equal-time covariance is therefore a classical symbol with matrix entries of order at most3. To obtain its local unequal-time form, start the unique formal negative graph at the second time s and use the complete separated same-sign transport to time t. The all-order state condition identifies its Cauchy polarization with the same formal graph. This gives, to any desired finite accuracy,

W(t,x;s,y)=sum_i integral exp(iP.(x-y)-i|P| integral_s^t omega_i(r)dr)
a_i(t,s,P) d^3P,

with each full matrix amplitude a_i of order at most3, plus a remainder as regular as required. All finite amplitude derivatives satisfy their symbol estimates. This is a local parametrix of the fixed state, not a new state chosen at s. Starting the transport at s eliminates any need to postulate unrelated preparation-memory phases. The exact chart transformations and unique formal negative graph give the same statement across the finite cover.

For every fixed collection of distribution seminorms, choose a sufficiently deep finite parametrix first. Its remainder and all needed derivatives then contribute smoother terms than the bounds below. Neither convergence of the formal inverse-k series nor a uniform-in-order constant is used.

## A conservative diagonal scaling bound

Let d_A,d_B be the TOTAL number of field derivatives in the two fixed bilinears. Their ordered product has a double-mode representation. On a dyadic internal shell |P|+|R| of size rho, its integration volume has order rho^6 and its complete amplitude has order at most rho^(6+d_A+d_B). Consequently the conservative power is

omega_bound=12+d_A+d_B.

This is deliberately not a minimal counterterm or physical power-counting claim. Its uniform distributional scaling estimate can be seen directly. Use center coordinates X and four normal coordinates Y=(t-s,x-y) for the diagonal. Localize a normal test function to an annulus away from Y0 and scale Y to lambda Y. For every required finite tangential derivative, the complete dyadic term obeys a seminorm bound

|T_rho,lambda(phi)| <= C_N rho^omega_bound (1+lambda rho)^(-N).

Here is why the decay is uniform. Where the normal time is separated from zero, the double phase has time frequency equal to the SUM of the two positive frequencies, bounded below by a fixed positive multiple of rho. Repeated time integration by parts gains inverse lambda*rho. Near equal times on the annulus, the spatial separation is nonzero; restrict to a sufficiently narrow equal-time region. The spatial derivative of the phase with respect to at least one high internal momentum is then separated from zero, since the speeds are uniformly bounded. Integration by parts in that high momentum supplies the same gain. Split the shell into regions where P or R is the high variable, so no angular derivative at a zero momentum is used. Amplitude derivatives are bounded symbols. Tangential derivatives of the oscillatory phase cost lambda*rho, because its time integral runs from s to t; choose N larger to absorb them.

For the time-retarded expression, multiply the commutator by theta(t-s) off the diagonal. On the first region this factor is smooth. On the equal-time spatially separated region the two-point kernels are smooth and the high-momentum integration by parts does not differentiate theta. Thus the same annular estimate holds. Gauge-reduced TT projection can produce a spatially nonlocal equal-time commutator; this causes no difficulty for multiplying a smooth off-diagonal kernel by theta, and it is one reason no spacelike-microcausality claim is made here.

Sum dyadic shells below and above rho=lambda^-1, choosing N>omega_bound. This gives a bound C lambda^-omega_bound. Extra harmless logarithms do not raise scaling degree. The constants depend on the compact background patch and finitely many test seminorms, not on lambda. The same construction with sufficiently many tangential derivatives gives the uniform microlocal scaling bound.

Each displayed phase at the diagonal has opposite covectors on its two legs. In center/normal coordinates the tangential phase covector tends to zero relative to the normal covector as Y tends to zero. The closure of the off-diagonal wavefront set at the diagonal is therefore conormal. Theta contributes only the corresponding time-normal covector. This verifies the required conormal closure and tangential regularity, rather than inferring them from an unspecified finite-order distribution.

For the complete seven-parameter scalar Hamiltonian family, the entire first and second vertices have spatial momentum degree at most4. They can be written as finite differential bilinears of total order at most4, with their full ordering retained. Thus the conservative bound is20. Tensor vertices have lower degree and obey this common bound as well. Nothing asserts that all allowed degree20 terms actually occur.

## Extension is not a physical counterfunctional

The finite scaling degree and conormal/tangential bounds meet the hypotheses of [Brunetti and Fredenhagen, Theorems5.3 and6.9](https://arxiv.org/pdf/math-ph/9903028). Apply their normal-coordinate Taylor subtraction, with a smooth compact cutoff equal to1 near Y0 and a smooth partition of unity along the diagonal. Subtract normal test jets through floor(omega_bound-4). The dyadic estimate makes the subtracted shell series convergent; the finitely many removed jets have arbitrary smooth coefficient functionals along the diagonal. This yields time-retarded extensions with the same bound.

Two such extensions differ by a finite sum c_alpha(X) partial_Y^alpha delta^4(Y), with |alpha|<=floor(omega_bound-4). For the displayed coefficient family the coarse upper limit is16. The off-diagonal expression is zero at t<s, and adding diagonal terms preserves that time-retarded support. No unique set of coefficients has been chosen here.

This is ONLY the extension of the specified reduced free bilinear response distributions. These coefficients have not been shown to arise from a single local covariant variational counterfunctional, to satisfy the parent Ward identities, to match the existing physical vacuum finite conditions, or to be numerically small. The finite-loop-order source ancestry, normal-ordering convention and existence of a time-retarded distribution do not supply those extra facts.

In particular, setting a reference Wick mean to zero is not a calculation of zero physical metric stress. The gauge/constraint and physical-vertex correspondence still matters. This note provides a continuum free-noise input and a controlled extension family, not the missing interacting P8 theory.
