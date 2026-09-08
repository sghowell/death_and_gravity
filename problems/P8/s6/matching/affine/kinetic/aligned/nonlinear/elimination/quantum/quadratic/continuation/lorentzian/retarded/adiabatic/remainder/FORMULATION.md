# P8-S6.66: finite homogeneous prepared mass response

Date: 2026-09-08. Original P8 remains OPEN.

Keep the physical clock, S6.55 vector state and S6.60 fixed profiles.
Use the same independent mass source as S6.64 and S6.65:

    a_m=1+alpha n, b_m=1+beta n,
    alpha=4/[9(1+u^2)^3], beta=28/[81(1+u^2)^3].

The source is spatially homogeneous, smooth, compactly supported
in time strictly after u0=-1/2, and zero on an initial neighborhood.
Bounds are on I=[-1/2,1/2], with m=m0 tau>=1000 and L=M tau>0.
This is the K=0 mass-only response, not a simultaneous metric lapse
variation or a full gravitational feedback map.

## Acceptance gates

1. Vary the actual eighth-order reference by the exact linearized
   Riccati equation and recover S6.65 through fourth adiabatic order.
2. Derive continuous C10 source bounds for its varied residual,
   normalization, phase and differentiated canonical readout.
3. Retain the nonzero selected-state initial mixing and bound the
   exact mixing variation with zero independent initial response.
4. Prove the varied reference-minus-subtraction tail is integrable.
5. Justify the common-dimensional finite limit and source derivative,
   then combine the finite remainder with the S6.65 local component.
6. Give a continuous C10-to-C0 bound, pin sources and pass regressions.

The subtraction order remains 0,2,4. The eighth-order reference
is an auxiliary comparison, not a new subtraction prescription or
a replacement of the exact all-order selected state.

## Scope

The finite response is the exact Gaussian selected-state linear
response in this prepared homogeneous mass-source sector, with the
same local prescription and no added finite Y^2 matching coefficient.
It is retarded in time. The bound uses the maximum over source
time derivatives zero through ten and is not a no-loss feedback norm.

It does not establish the arbitrary-spatial-momentum finite kernel,
full metric or second mass-source blocks, arbitrary varied Cauchy
states, coupled stability/cones, interactions, cutoff, finite Wilson
matching, common UV parent, V/G/B or original P8 closure.
