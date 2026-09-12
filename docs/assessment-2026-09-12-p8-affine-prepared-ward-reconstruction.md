# P8 continuation: conditional prepared Ward reconstruction

S6.215 supplies the prepared geometric reduction of a covariant Gaussian metric response to synchronous spatial inputs, together with explicit source/readout norm estimates and three remaining ordered scalar channels. Original P8 remains OPEN.

**Publication limitation:** the [phase erratum](assessment-2026-09-12-p8-retarded-phase-erratum.md) withdraws S213's physical current identification. Therefore S215's frozen claim to an already known tracefree input, and its application of the1e118 known-piece bound, are **not currently certified**. The new prepared Ward identities and the conditional implication from a correctly matched tracefree input remain the result published here. The frozen report is preserved as a historical record with this explicit limitation.

## Prepared source and detector maps

For g=diag(1,-a^2 I), lapse n, shift beta and Q, a vector xi=(eta,chi) produces

    n_xi=eta', beta_xi=chi'-a^-2 grad eta,
    Q_xi=2H eta I+sym grad chi.

Take the source primitive Iminus f=integral_(left)^t f and the detector primitive Iplus f=-integral_t^(right) f. In either case set eta=I n, chi=I(beta+a^-2 grad eta), and Qsyn=Q-2H eta I-sym grad chi. No division by H or spatial momentum occurs, so these maps remain regular at the bounce and at zero transfer.

The source vector vanishes at the unchanged initial germ. The detector vector vanishes at the final boundary. They are deliberately different maps, not a symmetric retarded Hessian construction.

## Ordered Ward identities and contacts

Use the actual contravariant weight-one one-point density E, with E00=-a^3 rho/2 and Eii=-a pressure/2. Its full Lie derivative is

    Lie_xi E=xi.partial E-Jac(xi)E-E Jac(xi)^t+div(xi)E.

The density conservation identity is a^3[rho'+3H(rho+pressure)]=0. The complete first metric variation h_D and second chart variation q_DG are those of S214, including lapse/lapse and shift/spatial nonlinear metric contacts even where the independent ADM Hamiltonian has zero corresponding second vertex.

For a covariantly renormalized prepared response, source covariance and detector conservation give separate ordered identities

    R(D,G_xi)=integral h_D:Lie_xi E+E:q(D,G_xi),
    R(D_xi,G)=integral -E:Lie_xi h_G+E:q(D_xi,G).

Their differentiated time flux contains both delta E and E:h_G. It vanishes at the final boundary from the detector vector and initially from the unchanged source/current germ. It is not permissible to delete either flux term or to assume a sharp finite computational band is separately diffeomorphism invariant.

The response equals R(Dsyn,Gsyn) plus these two Ward corrections. This structural assertion is conditional on existence of the original covariant distributional response; it is not a construction of the missing scalar operator or strong stress-norm differentiability.

## Ordered scalar and domain obligations

At nonzero P use unit trace I/sqrt3 and scalar tracefree(3ee^t-I)/sqrt6. Beyond a corrected full tracefree block there are three kernels: trace/trace, trace/scalar and scalar/trace. The last two are not interchangeable forward retarded kernels. Rotation covariance makes both cross entries zero at P0; the homogeneous trace anchor remains necessary.

The source primitive generally does not preserve a compact final support. The written audit shows that the original retarded estimates use the source's zero initial germ, not a final source cutoff; no uncontrolled cutoff is inserted inside the slab. The detector norm has no time derivatives, allowing the stated initial zero-output-domain adjustment. This domain analysis remains an input for a corrected response assembly, not a restoration of the withdrawn one.

## Quantitative conditional estimates

Let V04[D]^2=||(1-Delta)^2(n,beta,Q)||L2^2 and U138[Gamma]^2=sum_(r=0)^13||(1-Delta)^4 partial_t^r(n,beta,Q)||L2^2.

Exact primitive and Cauchy estimates give the synchronous detector factor below100 and source factor below10^20. The full one-point Ward correction is below10^35 V04 U138. If a correctly identified tracefree current obeys the stated2e95 input bound, the corresponding tracefree-plus-Ward contribution is below10^118 V04 U138. That input's physical identification remains withdrawn pending reassembly; the conditional inequality is not a certified existing known piece.

If all three missing scalar entries also admit constants L_ab in the required norms, the complete conditional bound is

    [1e22(2e95+sum L_ab)+1e35] V04 U138.

No L_ab is supplied by this checkpoint, and no reduced canonical scalar normalization or inverse follows.

## Verification and immutable record

Core kinematics, ordered Ward, scalar projection and norm checks passed. Independent science passed101 tests in15.61 seconds, integrated289 in16.63 seconds, preflight289 in15.62 seconds, final private289 in16.47 seconds and repository289 in16.37 seconds. Tests include direct inverse-metric/determinant volume currents, full source/readout Ward variations, integrated action reconstruction, boundary/contact counterexamples, normalized arbitrary-axis rotations and primitive/derivative majorants. Lint passed before freezing.

Fresh original-SymPy ordinary replay passed314 tests in2346.07 seconds and CLI passed. Complete captured P8 regression passed42377 tests in4248.35 seconds with final exit code0, retaining683 files. Snapshot path-list SHA:

    320c0ec5d6f4b9b5563754584f9649e3ba974f9682e4756a2288b71eab4d318c

The full-only exact-GCD adapter passed128 original tuple comparisons; final counters were35543 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native/direct/ordinary/CLI retain original SymPy.

The38745-character native report was transferred in four checked chunks. Its twenty fields contain38 named identities,100 scalar entries,38 gates,nine controls and147 rejected inputs. All eighteen frozen source hashes and the report SHA were verified:

    bba0d9e34fd6645444a52940587dced5500746a8fb49124434cce1ad819059a5

Exactly22 publication files exclude private S216 and unrelated P4/P9 work. These passing tests are historical reproducibility evidence; they do not fix the independently discovered ancestor phase error. Continuous covariant/distributional/Sobolev arguments are written proofs, not FORMALIZED.

## Continuing work

Private S216 now explicitly corrects the phase, all three scalar poles and the full six-invariant local finite part, with467 integrated tests passing. Its local finite bound does not yet restore S213 or S215's known-input application. Correct finite endpoint/bulk/regulator assembly, scalar anchor and complete scalar remainders remain next, followed by reduced inverse, quantum background and original V/G/B obligations. No user-intervention blocker has been identified.
