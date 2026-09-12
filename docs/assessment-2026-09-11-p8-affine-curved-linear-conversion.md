# P8 continuation: actual curved linear cutoff-conversion coefficient

S6.210 evaluates the complete actual curved linear coefficient in the original sharp-band UV conversion. It retains every physical sector, the longitudinal curvature shift and both source-time derivatives. Original P8(b), V/G/B and P8 remain OPEN.

## Result and boundary

For p=|P|, T=tr(D Gamma), V=(D phat).(Gamma phat), W=(phat.D.phat)(phat.Gamma.phat), define

    L=18T-12V-W, M=42T+4V+3W,
    S=298T-420V+63W, C=2T+4V-W.

On the unchanged CD clock with a=(1+t^2)^2, H=a'/a, Q=H'+2H^2 and m=1000,

    A3=-p L/(512pi^2 a),
    A1=p/pi^2[-m^2 a M/1024+p^2 S/(24576a)
                   -a Q C/1024+a L(Gamma''+H Gamma')/2048].

The source tensor in each form is the indicated time jet. Both multipliers extend continuously as zero at P0. The full actual WKB calculation, not a substituted flat physical state, fixes the curvature and time-derivative terms.

The actual first WKB coefficients obey P1_T(z1)=0 and P1_L(z1)=-Q/2. The mixed amplitude's imaginary second-degree term is retained. The complete j1/d1 slot evaluates to zero, while its higher nonzero slot and finite odd endpoints remain. Source differentiation precedes detector/source coincidence.

Seven exact general-tensor azimuthal identities reconstruct all five tracefree components. The first/second source derivatives have the proper-time divergence form, with an explicit Green-identity omission control. Coefficient bounds are1/100 for A3 in the stated gradient norm and1e4 for A1 in Z23. Their canonical displays are4e-802 and4e-796.

These are coefficient bounds, not cutoff-uniform response bounds: K^3 and K still grow. The complete quadratic and finite conversion cancellations and all-P1e40 X46/K tail from S6.209 are unchanged. Full fixed spatial matching, current assembly and reduced inverse/background work are separate.

## Verification

All initial scientific probes and formulas passed. The initial independent suite passed198 tests in48.96 seconds; expanded science passed264 in94.52 seconds. Integrated science passed447 in93.87 seconds; final preflight passed447 in89.90 seconds; final private science passed447 in96.12 seconds; repository science passed447 in92.96 seconds. Unused-variable/import lint corrections were made before freezing; final lint passed.

Independent original-SymPy ordinary replay passed472 tests in2384.56 seconds. Original-SymPy CLI replay passed. Full regression passed40818 tests in3818.15 seconds with final exit code0, retaining all673 captured files. Snapshot path-list SHA:

    7d15618a37f9b9baf4730ee533b100532bc5672d73206ba53ae4a98c23767fd1

The full-only exact GCD adapter passed128 original tuple comparisons. Final counters were34443 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct science, ordinary and CLI retain original SymPy.

The38107-character native report was transferred losslessly in four chunks and verified against eighteen frozen source hashes. Its twenty fields contain49 named identities,92 scalar entries,40 gates,nine controls and132 rejected inputs. Report SHA:

    7aa73d597f35e68b485e1aed2f0e6ed8c3a72648d992f1f555a8d5fed02d3612

Independent actual full W8 constrained fields, all nine pairs, multiple clocks/directions/noncommuting tensors, separated transfer scales, fixed-detector source derivatives, complete angular quadrature and two Cauchy resolutions check the result. The continuous Fourier-norm and regulated-current arguments are not FORMALIZED. No frozen scientific byte changed.

This exact22-file publication excludes S6.211/S6.212 and unrelated P4/P9 work.

## Continuing frontier

S6.211 is frozen with the complete actual spatial UV difference, full contact cancellation only in that difference, and the original curved pole identity. Its integrated science passed288 tests in338.99 seconds, final preflight288 in259.97 seconds, final private288 in339.74 seconds and repository288 in345.79 seconds. Its38001-character native report and all eighteen sources are verified. Independent ordinary, CLI and full replays are still in progress.

S6.212 is frozen with the original dimensional spatial UV finite part, including genuine general-dimensional magnetic geometry, all actual source/mode jets, fixed-source invariant derivatives, the evanescent odd endpoint, continued Euler and volume terms, and the fixed MSbar convention. Expanded independent science passed129 tests in224.63 seconds; integrated338 in235.55 seconds; final preflight338 in112.13 seconds; final private338 in234.42 seconds. Its native report and repository replay are underway. It is not yet a published completed-replay checkpoint. The finite local coefficient alone does not establish the full assembled response.

Private S6.213 has an explicit proposed anchored assembly using the same fixed homogeneous current plus the known spatial remainder difference and the new finite UV difference. The dimension-limit proof, Sobolev lift and complete original-regulator estimates remain under development; no private assembly result is certified here.

Original remaining obligations include the complete spatial/scalar/metric inverse, finite-coupling background and stability, remaining parent measure/sectors/loops, threshold/heavy/cutoff and omitted-order matching, canonical vacuum cuts/contour/truncation and finite-gravity IR/Regge. Scoped P8(a) and A.20-A.23 are unchanged. No user-intervention blocker is identified.
