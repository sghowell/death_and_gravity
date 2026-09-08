# P8: source alignment passes the quadratic cone test

Original P8 remains open. The completed scoped P8(a) photon objective,
linear CD/M1 classification and adopted V/G/B contract are unchanged.
This continues the [mass/cone obstruction](assessment-2026-09-07-p8-retuned-trace-cone.md)
without modifying that frozen action or its rejected verdict.

## New result

[S6.42](../problems/P8/s6/matching/affine/kinetic/aligned/FORMULATION.md)
replaces the curled trace T by W=T-B(phi,x)dphi, keeping the complete
source-centered mass, original scalar functions and physical metric.
The coefficient B=-15H(1+x)/(8h) is local, not a prescribed time source.

The exact shifted stationary source is

    Sstar_normal=(4p²-1)*(K_hat-3H*s)+(3/2)*s*Q_lower.

Its full background and first variation vanish, as independently
checked from the all64 source. The triangular temporal-variable map
has unit Jacobian and introduces no lapse velocity. The quadratic
physical action, including all lapse, shift, vector-temporal and
free-matter constraints, is exactly old CD/M1 plus positive Proca.
All three scalar and both transverse-vector principal cones are the
original matter cone, and the old regular first-order crossing remains.
The tensor block is unchanged. This removes the actual S6.41 cone
obstruction at quadratic rolling order, not by weakening the contract.

The actual time-dependent vector normalizations give

    Omega²=q+1/zeta-g''/g >= q+1/zeta-15 >= q+1985,
    0<zeta<=1/2000, q=k_com²/a².

This holds at all real times and momenta, with a separate homogeneous
coordinate-vector chart at zero momentum. Fixed comoving, not physical,
momentum is differentiated; nonunit physical scales are tested.

## What is still missing

The source has a nonzero second variation and a cubic vector-light
coupling. The old coefficient ODE gives explicit normal-source and
spatial-electric bounds on the original tube, but not an inverse
operator or an omitted-action error. The canonical frequency floor
is not a stationary scattering gap or an interacting cutoff.

The next [working calculation](../problems/P8/s6/matching/affine/kinetic/aligned/response/FORMULATION.md)
derives a prepared leading retarded response with explicit forcing
derivatives, initial data, finite time/momentum domain and physical
vector readouts. Its assumptions must not be confused with a general
on-shell light evolution or a full nonlinear elimination theorem.

Full nonlinear secondary constraints, controlled nonlinear and loop
remainders, an interacting cutoff, and the adopted vacuum/finite-gravity
V/G/B UV gates remain research tasks. None currently requires a new
user choice or external authorization. This is not P8 closure.

## Verification

The frozen source-aligned report covers 17 sources, 35 named identities
comprising 60 scalar entries, 24 continuous/interface proof checks and
34 rejected inputs. The scientific suite passes 68 tests in 43.22
seconds. The complete ordinary suite passes **88 tests in 166.35
seconds**, without the broad GCD adapter; the separate seeded CLI passes.
Its source and parent hashes are read-only and all frozen ancestry is
fully rebuilt. The full **4329-test P8 regression passes in 750.49
seconds**, with no certified checkpoint excluded. The exact GCD adapter
passes 128 original normalized-tuple comparisons and records 5599
domain fallbacks and 6509 exact descents, using the unchanged
[per-test seeded recipe](assessment-2026-09-07-p8-constant-two-trace-family.md#verification).

These are exact symbolic checks, continuous written proofs and
independent internal calculations, not proof-assistant formalization
or external peer review.
