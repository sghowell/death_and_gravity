# P8: a controlled prepared leading vector response

Original P8 remains open. The scoped P8(a) photon objective, linear
classification, physical matter frame and adopted V/G/B contract are
unchanged. This follows the [source-aligned quadratic checkpoint](assessment-2026-09-07-p8-source-aligned-modes.md).

## New result

[S6.43](../problems/P8/s6/matching/affine/kinetic/aligned/response/FORMULATION.md)
keeps that action unchanged and studies its actual second-order heavy
response with zero heavy initial data. Differentiating the original
coefficient ODE fixes the nonlinear source explicitly:

    S2=-(2/h)*n*delta_K-[6H/h+27H/(8h²)]*n².

After the complete temporal constraint, the actual canonical forcing
is J=(g_L²*S2)'/g_L, not S2 alone. The fixed-comoving frequency and
its first two derivatives obey continuous bounds W>=19985, |W'|<=78,
|W''|<=712 on the declared domain. The energy argument then gives

    ||v-zeta*J|| <= zeta*A/80,

for a prepared C² forcing with all three jet norms <=A. The physical
spatial-vector error is <=6*zeta*E and the temporal error against S2
is <=500*zeta*E, with physical factor 1/tau. Here E bounds the source
and its first three time derivatives, and
A=(815/2)*sqrt(zeta)*|k_com|*E.

The exact domain is a forward interval inside [-1/2,1/2], output
0<k_com²<=1, 0<zeta<=1/20000, zero heavy data, and initial
S2=S2'=S2''=0. A sufficient light preparation is n=n'=n''=0.
The homogeneous temporal and unforced transverse components are
separate controls. Unprepared forcing, arbitrary heavy states,
amplitude-only forcing bounds and out-of-window momenta are not included.

Continuous coefficient bounds give the sufficient source envelope

    E=514*epsilon_n*epsilon_K+(132027/8)*epsilon_n².

The light envelopes control derivatives through order three. In
Fourier space the products are convolutions with summable modewise
envelopes. The momentum restriction is on the source's output mode;
an incoming cosine produces zero and double momentum. Both the
canonical map and the physical readout normalization are retained.

## What this does and does not resolve

This supplies a quantitative retarded estimate for the actual leading
block under explicit conditions. It is stronger than a source norm or
a positive canonical frequency alone. It does not bound the higher
perturbative orders, the complete induced action, quantum loops, an
interacting cutoff, a stationary gap or the adopted V/G/B UV gates.

The next calculation tests whether a nonzero solution of the original regular linear light
system realizes the preparation and derivative bounds through the
bounce. Missing nonlinear secondary constraints and full nonlinear
control remain separate. No new user choice or external authorization
is needed for these calculations; original P8 is not closed.

## Verification

The frozen report pins 15 local sources and fully rebuilds S6.42 and
its ancestry. It checks 34 exact scalar identities, 30 continuous/interface
proof checks and 34 rejected inputs. The final scientific suite passes
48 tests in 22.93 seconds. The ordinary suite passes **70 tests in
156.16 seconds**, without the broad GCD adapter, and the separate
seeded CLI passes. The full **4399-test P8 regression passes in
771.28 seconds**, with no certified checkpoint excluded. The exact
GCD adapter passes 128 original normalized-tuple comparisons and
records 5602 domain fallbacks and 6509 exact descents, using the
unchanged [per-test seeded recipe](assessment-2026-09-07-p8-constant-two-trace-family.md#verification).

These are exact symbolic identities and written continuous proofs,
with independent manufactured-solution, constraint, convolution and
readout controls. They are not proof-assistant formalization or external
peer review.
