# P8 continuation: complete ordered spatial Gaussian response

S6.219 completes the three ordered nonzero-transfer scalar kernels and assembles the full reference spatial Gaussian response. Original P8 remains OPEN. The [homogeneous anchor](assessment-2026-09-12-p8-affine-homogeneous-trace-anchor.md), [corrected tracefree current](assessment-2026-09-12-p8-affine-corrected-spatial-current.md) and [original finite local term](assessment-2026-09-12-p8-affine-ordered-scalar-symbol.md) enter with their stated normalizations. The [phase erratum](assessment-2026-09-12-p8-retarded-phase-erratum.md) remains in force; old S207/S212/S213 physical formulas are not restored.

## Full trace and constraint geometry

The ten-component physical feature retains electric, magnetic, mass and temporal-constraint terms. For B_Q=Q-tr(Q)I/2, the current feature is diag(-B_Q,-B_Q,tr(Q)/2,B_Q), and the contact includes the symmetrized two-matrix products and tr(D)tr(G)/4. Its complex Frobenius bound follows from
||B_Q||²=||Q||²-|tr(Q)|²/4 and |tr(Q)|²<=3||Q||².
Both longitudinal legs and all mixed constraint terms remain. Independent literal ADM determinant/inverse differentiation checks the first and mixed second vertices, including arbitrary noncommuting matrices with nonzero trace.

The complete real covariance-flow derivative agrees with the nine annihilation-pair contractions, including the positive detector phase, negative source phase and +Im convention. The ordered endpoint multiplier is i(-i)^j. The complete sixth-order bulk uses the same branch. An independent fixture with nonzero symplectic preparation shears makes the wrong-phase control nondegenerate.

## Original cutoff, dimensional continuation and all remainders

All three scalar channels are evaluated separately. The unaveraged odd nonlogarithmic geometry is nonzero; only the complete centered total parity is available. Forty-five raw shell rows include q0..4 and every source jet. Thirty nondecaying grades give the original cubic and linear conversion coefficients. The quadratic and finite terms cancel only after the full hemisphere and grazing strip are retained. The q4 logarithmic shell is nonzero and contributes to the controlled 1/K remainder.

In particular, the cubic trace/trace coefficient is -p/(16 pi² a), while each cross coefficient is +p/(64 pi² a), before the stated orthonormal scalar conversions. The linear cross coefficients are distinct. Their local formal transpose is checked; equality of the two retarded cross kernels is not asserted.

The full complex-dimensional calculation checks 96 geometry and inverse-basis majorants, 35 grouped transverse moments and 30 exact normalized-to-original reconstructions. The common |d-3|<=1/4 domain, both phase rates, all four normalized W8 defects and the far/near angular powers are retained. These estimates extend the actual longitudinal/trace geometry, not just the old tracefree labels.

The full feature transfers the existing state, time and contact comparisons with pair constant below1e9 and error constant below1e15. The corrected known-piece bound is5e48 with tail5e53/K. The complete shape remainder is bounded by1e40/K. No state, regulator, contact or finite counterterm is changed.

## Assembled response and scope

With the S218 full homogeneous anchor and S216 original finite local term, the anchored assembly is
J_ren(P)=H0+Known(P)-Known(0)+F_finite(P).
It has the uniform full-spatial weak bound2e95 M[D] Z136[G] and original-regulator error3e54 M[D] Z136[G]/K. The normalized trace/trace bound is2e95 and each ordered cross bound is3e49. Both cross anchors vanish at P=0 by rotation covariance, which does not equate the nonzero-transfer kernels.

The canonical displays are8e-705 and12e-746/K in their declared weak norms. Thirteen source time derivatives and six spatial derivatives remain. These are not bounds for a reduced inverse or a nonlinear solution.

The full prepared Ward reconstruction now has its spatial input, giving the1e118 V04[D] U138[G] bound, with all local density and ADM chart contacts. The reference identification uses S176's vanishing source and first source variation at the reference. It does not make the finite-amplitude affine parent source free. Only the precisely stated S215 matching row is advanced; the remaining matching obligations and all nine primitive statuses are unchanged.

## Independent checks and fresh immutable validation

Independent science passed37 tests in84.96 seconds. The original phase-control fixture had33 passes and four failed controls because its unsqueezed source was real; adding preparation shears before freezing made the control sensitive without changing the physical current. An integrated run later had351 passes and one case-sensitive text assertion failure; that assertion was corrected before freezing. Final private science passed352 tests in83.25 seconds and repository science passed352 in82.83 seconds. Lint passed.

The first fresh ordinary replay failed collection with ModuleNotFoundError for the frozen independent helper, and the first full replay also stopped at collection after660.82 seconds. Neither is counted as a pass. Importlib collection did not add the helper directory to sys.path. Frozen scientific, proof and test bytes were left unchanged. The fresh retry explicitly supplied
`PYTHONPATH=problems/P8/s6/continuation/s6_219/tests`
to the ordinary and full launch commands. This is an import-path allowance only.

The corrected-launch ordinary replay passed377 tests in2361.01 seconds. CLI replay passed. Full regression passed43980 tests in4391.51 seconds with final exit code0. Its691-file snapshot SHA is
`44fed37ff3d6cdae6a860c8a63e8dd78c61699e8bbb68490c28566d09c626d1b`.
Only full regression used the audited exact-GCD adapter:128 original tuple comparisons passed, with final counters36063 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

The60711-character native report was transferred in six checked chunks. Its20 fields contain73 named identities,457 scalar entries,68 gates,nine controls and171 rejected inputs. All20 frozen source/proof/test hashes and the report SHA were verified:
`df188e09486c7666c0dcf9c811702d1be096849acf993c8efae18bfffc9532a6`.
The exact24-file publication manifest includes those20 inputs, report, this audit, CLAIMS and README. S220, private propagator work and unrelated P4/P9 changes are excluded. Infinite-mode, dimensional, angular, Sobolev and Ward arguments are written proofs, not FORMALIZED.

## Remaining work

The full quantum constraint inverse, common functional-space derivative control, finite nonlinear sourced-parent remainder, quantum-corrected background, physical cutoff/heavy-sector control and original V/G/B remain unresolved. The reduced classical scalar and infrared pullback successor is undergoing separate validation. Scoped P8(a) and A.20–A.23 are unchanged. There is no user-intervention blocker; continuation is active.
