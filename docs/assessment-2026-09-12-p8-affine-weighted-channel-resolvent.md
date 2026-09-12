# P8 continuation: weighted original-channel and curved reference inverses

S6.226 proves a quantitative weighted inverse for the original trace/shear factors with a stated bounded causal correction, then composes it with the actual S225 curved coordinate inverses. It does not assert that the complete actual curved quantum remainder belongs to that bounded class. Original P8 remains OPEN.

## Uniform first-sheet logarithmic lower bound

LetA0=-Ftrace andA2=-F2 be the original factors, including the fixed finite coefficients and the original shear pole. The existing trace bound isReA0>=16/15. The full first-sheet shear bound proved here isReA2>=-172/225, from the actual cut, finite threshold, uniform large-circle asymptotics and the harmonic minimum principle.

For both channels, on the first-sheet exterior|p|>=4m^2,

ReAi(p)>=(13/80)log(|p|/(4m^2))-4.

The proof treats both cut banks, the inner circle and a large closing circle. On the cut, u=log(tau/(4m^2))>=16 gives a direct lower logarithmic bound from the actual positive polynomials; for0<=u<=16 the global channel bounds supply the required margin. The first-sheet identity

atanh(1/sqrt(d))=[Log zeta+2Log(1+sqrt(d))]/2,
zeta=p/(4m^2), d=1+1/zeta,

gives uniform complex exterior asymptotics by analytic expansion in1/zeta. A positive-real-axis expansion alone would not establish the required contour bound. The two original leading logarithms are2 and13/60, both larger than13/80.

## One weight controls every temporal and spatial frequency

Forlambda=sigma+i omega, q=P^2>=0 andp=lambda^2+q, the exact identity is

|p|^2-sigma^2(sigma^2+omega^2+q)
 =(q-omega^2)^2+sigma^2(q+omega^2)>=0.

Thus|p|>=sigma^2 uniformly inq andomega. The right-half-plane Laplace contour never intersects the first-sheet cut. For a stated matrix boundM>=0, choose

sigma_M=2m exp(32+8M),
d_M=(32+13M)/5.

The original diagonal inverse, including the factor3/8 in the shear channel, has norm at most

||Kdiag||<=5/(32+13M)

on weightedL2_tH^r, all spatial momenta and every realr. Plancherel is applied to the original causal distribution, using its bounded spectral primitive to justify the Laplace transform and completion. No different Green function or pole prescription is introduced.

## Ordered bounded correction and its graph

LetV(t) be a measurable, spatially independent2x2 matrix with essential-sup norm at mostM. It may mix channels and need not be symmetric. Then

||Kdiag V||,||V Kdiag||<=5M/(32+13M)<5/13<1/2.

The convergent causal Neumann series gives

K_V=(I+Kdiag V)^-1 Kdiag=Kdiag(I+V Kdiag)^-1,
||K_V||<=5/(32+8M).

The equality is an ordered series identity, not a commutation assumption. Both original forward/inverse distribution identities retain the complete initial boundary. The middle graph consists ofu in the weighted source space for which(Fdiag+V)u is an ordinary source there. Existence follows by applying the original forward distribution to the bounded equation; the original inverse gives uniqueness and the reverse identity on that graph.

The proof works for any causal bounded middle operator once its actual norm and domain compatibility are established. The certified comparison class here is explicitly norm-bounded time multiplication. Differential order, formal matching or a derivative-losing weak estimate does not put an unknown physical remainder in that class.

## Composition with the actual curved coordinates

Using the unchanged S225 coordinatesBc and their causal inversesY,Z, the reference is

A_V=Bc*(Fdiag+V)Bc,
E_V=Y K_V Z.

The inherited coordinate bounds give

||E_V||<=125/[4(32+8M)(sigma_M-20)^2(sigma_M-108)^2].

The compatible composed graph requiresx=Bc s andy=(Fdiag+V)x in the weighted source space andBc*y an ordinary source there. Equivalently,s=Yx in that graph. This makes products with merely measurableV well-defined; it is not an unrestricted maximal distributional realization or the full S222 graph.

The physical conformal force inverse retains RIGHT multiplication by64pi^2 kappa a^4. Since the actuala^4<6, its norm is at most

12000pi^2 kappa/[(32+8M)(sigma_M-20)^2(sigma_M-108)^2].

Removing the weight on a finite interval costs exp(sigma_M T). This cost and kappa are not suppressed.

## An explicit check against confusing existence with stability

The bounded comparison

Vcmp=diag(0,(8/3)A2(4m^2))

has norm below208/315<1, using the original positive radial weight. Yet the corrected shear channel has a simple zero atp=4m^2. Atq0 this is a growing Laplace polelambda=2m. The weighted theorem withM1 still applies, at a contour to its right.

The pole is retained. This is not a change of the physical prescription or a claim that the actual curved system has that pole. It is a concrete demonstration that weighted causal existence does not imply unweighted smallness or stability.

## Independent checks and completed fresh validation

Independent exact checks verify the global/exterior constants, complex-frequency inequality, channel normalizations, noncommuting inverse products and physical density order. High-precision radial comparisons include eight complex closed-form points and132 additional exterior/contour cases. The growing-pole comparison is checked with the original mass and weight.

Two early private probe comparisons needed factoring of the entire exact difference rather than comparison of differently factored expressions; no formula or physical coefficient changed. Corrected private science passed265 tests. Preflight checked17 inputs and20 fields in0.42 seconds. Final private science passed265 in1.71 seconds, and repository science265 in1.81 seconds. Formatting and lint passed before freezing.

Fresh original-SymPy ordinary replay passed290 tests in2311.92 seconds. Standalone CLI replay passed. Complete P8 regression passed46058 tests in4491.24 seconds, exit code0. Its705-file snapshot SHA is

`d790eec4acdcccc9bf7c14ef247c2a1d4ff1c2a67e60d86f497dc0d6ef9d5b39`.

The frozen S219 helper-directory PYTHONPATH allowance was retained in the full launcher. Only full regression used the audited exact-GCD adapter:128 original tuple comparisons passed; final counters were36072 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

The38960-character native report was transferred in four checked chunks. Its20 fields contain35 named identities,35 scalar entries,25 gates,nine controls and183 rejected inputs. Every frozen input hash and the native report SHA were verified:

`98ac59825e67545f5afd2ce5a430d0c381f49811d85ee5bf805a3939b8c9de3d`.

The exact21-file publication manifest contains17 inputs, report, this assessment, CLAIMS and README. Pending finite-local and actual homogeneous successors, private coupled research and unrelated P4/P9 changes are excluded. The complex analysis and functional arguments are written proofs, not FORMALIZED.

## Remaining work

The next finite-local checkpoint proves a compatible bound for the actual original local remainder; a further checkpoint addresses the actual homogeneous nonlocal shear problem. Both are undergoing fresh verification. The complete nonzero-transfer curved mass/state/contact response, coupled scalar/clock/matter graph and nonlinear/finite-coupling common-parent problem remain further work. No original primitive or prior matching status is closed by this reference result, and no user intervention is currently required.
