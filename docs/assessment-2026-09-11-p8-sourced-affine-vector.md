# P8 continuation: exact sourced affine vector and scoped bounds

S6.175 keeps the same CD-REG-AFFINE-ISO classical parent of S6.174.
Original P8(b), V/G/B and P8 remain OPEN.

## What is established

The complete nonlinear source and its source-squared contact are retained.
Writing K for the Maxwell differential operator and D for grad-div,
the covariant identities KD=DK=0 factor the constrained Proca operator:

    O=I+zeta K, N=I+zeta D, L=ON=NO=I+zeta(K+D).

The normally hyperbolic reduction and checked global hyperbolicity of the
actual CD geometry give unique compact-source causal Green operators.
They preserve div(W)=div(S); this is not a transverse-source assumption.
Qualitative Green existence is not a state or norm estimate.

If dS=0, W=S solves the full vector equation and its extra first metric
and light variations vanish. Every homogeneous light-target solution
therefore has an exact classical parent lift with compatible vector data.
This does not cover arbitrary vector initial data or assert retarded
preparation from a noncompact infinite past.

At one common regulator and one stated boundary prescription, the exact
Gaussian Schur identity keeps both the source contribution and determinant.
The determinant is light-field independent only at fixed metric.
A retarded inverse in a single-branch quadratic action instead varies to
the retarded/advanced symmetric kernel; a causal response needs the state
and closed-time-path construction.

## Distinct quantitative observables

For real Schwartz canonical fields with Fourier support in the Euclidean
unit ball and Fourier L1 norm at most one, the COMPLETE analytic source
has the bound

    ||S||_2 <= 256 n kappa^(-2) ||Phi||_2.

The proof keeps the rational switch, exponential factor and all derivatives;
it is not a finite Taylor truncation. For n=1024 and kappa=1e800 the
separately defined real Euclidean Gaussian observable obeys

    0 <= Delta S_E < 1e-2389 ||Phi||_2^2.

This is not obtained by declaring the Lorentzian field class Wick rotated.
For the leading degree-four vacuum source S4, its finite support lies below
the vector pole for 0<zeta<=1/2000. The complete Lorentzian inverse gives

    |Delta S_L[S4]| < 1e-2392 ||Phi||_2^2.

The nonlinear complete source is not band-limited merely because Phi is.
The generic real-time inverse is unbounded near its physical mass shell,
but this is NOT a proof that the actual nonlinear source image contains
the bad packets, and is not a parent or whole-row exclusion.

## Verification

Final private science passed233 tests in7.97 seconds after two private
normal-form/roundoff controls and exact-serialization checks were corrected.
Fresh repository science passed233 in8.22 seconds.
Ordinary replay passed258 in2097.56 seconds; independent native CLI passed.
The complete captured P8 snapshot passed30002 tests in3277.86 seconds,
with final exit code0 and all603 captured test files unchanged.

Path-list SHA-256:

    5cec5789d51db9436489896a82b0e88ffa6bc14e8f4cdb3baf68eb15a140e07e

The full-only exact GCD adapter passed128 original tuple comparisons.
Final counters:34147 domain fallbacks,7340 exact descents,88 mixed fallbacks.
Native, direct science, ordinary and CLI used original SymPy.

The native41984-character report was transferred losslessly in four chunks
and independently checked against all18 frozen scientific source hashes.
Its20 fields contain47 named identities,213 scalar entries,19 proof gates,
9 controls and114 rejected inputs. Report SHA-256:

    151f0642185a494dd1a49e77360349c5105596fbf4485cdb911586ef0d480b94

Independent controls cover constrained matrix inverses, source signs and
contacts, homogeneous solutions, finite-regulator translations, full
analytic coefficient bounds and the on-shell boundary. These are exact
algebra and written proofs, not FORMALIZED. No frozen byte was changed.
Exact22-file staging excludes S6.176 and unrelated P4/P9 work.

## Next frontier

S6.176 has passed private preflight and207 science tests and is frozen for
repository/native checks. It specifies the conditional sourced Proca
Hadamard family, its canonical measure and exact clock operator/readouts,
and a fixed covariant prescription for the ordinary-vector clock stress.
It does not import the old scalar stress-canceling profiles or a full
interacting parent.

The next construction holds the complete independent canonical scalar
functions fixed along an actual gravitational-decoupling family through
the S6.174 base point. Explicit gravity-dependent terms must be tracked,
not deleted at the level of diagrams.

Full real-time response, interacting state and measure, physical cutoff,
omitted-order/background bounds, vacuum contour/cuts/truncation and
finite-gravity IR/Regge control remain research obligations. Scoped P8(a)
and A.20-A.23 are unchanged. No user-intervention blocker is identified.
