# Pole-retaining inverse of the actual homogeneous total equation

LetK2 be the ORIGINAL S223 retarded inverse ofF2(Dt^2), including its isolated subthreshold pole and complete cut. The proved spectral representation gives

K=||K2||L1(0,1)<infinity.

The full kernel is not inL1 on the half-line: its undamped pole remains. Only the finite unit window is used below. No numerical K is invented.

LetC be a valid finite majorant for the COMPLETE actualV_total in notes/remainder.md and notes/primitives.md, and put beta=KC. The norm of convolution byK2 in the weighted continuous space is at mostK. Thus

||K2 V_total||<=2beta/sqrt(lambda).

Choose the explicit finite weightlambda=(4beta+1)^2. Then

2beta/(4beta+1)<1/2,

with exact margin1/[2(4beta+1)]. This is a mathematical time weight depending on proved finite kernel constants, not a retuning of a physical coefficient, state or source.

The bounded inverse ofI+K2V_total is its norm-convergent Neumann series. Every term is causal, so its limit is causal. For normalized forcingg, define

h=(I+K2V_total)^-1 K2 I4g.

The order is essential. K2V_total is not commuted toV_total K2 or throughI4. Independent noncommuting finite causal matrices verify the correct products and detect both wrong orders; they are diagnostics, not a continuum discretization certificate.

The bound on the original unit interval is

||h||C0<=2exp(lambda)K||I4g||C0.

For canonical forcingf, g=-16pi^2 kappa f, so

||h||C0<=32pi^2 kappa exp(lambda)K||I4f||C0.

C itself includes the kappa-dependent classical coefficients. The exponent and source normalization are retained. This is constructive in the actual finite but unevaluatedC andK; it is not a numerical smallness or stability bound.

The original causal identitiesF2*K2=delta=K2*F2 include the complete initial distributions. ApplyingF2 to the bounded equation gives
(F2+V_total)h=I4g in causal distributions. Conversely, applyingK2 to this normal-form equation recovers the bounded equation and its unique solution. No homogeneous solution or pole is discarded by hand. Smoothness and equivalence to the unintegrated actual equation are established in notes/regularity.md.
