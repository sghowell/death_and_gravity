# P8 continuation: complete four-point loop and explicit value matching

S6.235 computes the complete first four-light loop of the separately named V2S-T1 model. It then defines V2S-T1-OS4 by one explicit finite value-matching condition. No original affine/DHOST/Proca action, state or subtraction prescription is overwritten. Original P8 remains OPEN.

## Complete representation and independent normalization checks

With n=M_H²=D+2 and A_s=C+g²/(n-s), the full first-loop amplitude is

A1_base=(16pi²)^-1 sum_channels[A_s² B_MS(s)/2+2A_s g² Cbar(s)+g4(Dbar(s,t)+Dbar(s,u))].

All three bubbles, all triangles and all six ORDERED boxes remain. In Dbar(s,t), the first channel cuts two light lines and the second two heavy lines; exchanging the arguments is not an identity. The full two-field Gaussian elimination and the fourth derivative of the light trace logarithm independently fix the diagram multiplicities.

For equal external virtuality e, the complete parameter denominators are

Q_C=1-z+nz-ez(1-z)-s(1-z)²xi(1-xi),
Q_D=Q_C-tz²eta(1-eta).

Cbar integrates (1-z)/Q_C, and Dbar integrates z(1-z)/Q_D², with all parameters in [0,1]. These are Feynman boundary values across cuts, not absolute-value integrals over real double poles. On shell e=1; the independent zero-momentum off-shell check has e=s=t=u=0.

Only the bubbles have the UV pole Delta sum A_s²/(32pi²). Combined insertions of deltaC=-3C²Delta/(32pi²), deltag=-CgDelta/(32pi²) and deltan=g²Delta/(32pi²) cancel it exactly, including the simple and double heavy denominators. The heavy one-point and light on-shell quadratic conditions remain those of S234 and add no direct first-order four-point contact.

For 4<s<n, the complete forward discontinuity agrees with beta integral A_tree²/(64pi), with the full rational angular dependence and infinitely many even waves. Independent triangle delta-function and heavy-mass differentiation checks reproduce the box cuts. The heavy pole and heavy-pair threshold remain outside this open comparison window.

The [primary scalar-integral conventions](https://arxiv.org/abs/0709.1075) and [functional one-loop matching treatment](https://arxiv.org/abs/1604.01019) supply notation and context; the actual action, diagrams, signs and bounds are derived in the frozen package.

## Two distinct value checks

An independent constant-background full two-field Hessian gives the complete off-shell zero jet. Its triangle and box integrals reduce to

B21=(n log n-n+1)/(n-1)²,
B22=((n+1)log n-2(n-1))/(n-1)³=-partial_n B21.

The resulting zero jet is negative at the actual parameters and below one millionth of the positive classical quartic margin. It is not substituted for the on-shell matching point, and no off-shell zero-jet counterterm is adopted.

At the on-shell subthreshold symmetric point s=t=u=4/3, every parameter denominator exceeds 2/3+D z. This point is not a real scattering configuration. Complete integral bounds imply

A1_base(s0,s0,s0)<-5985g4/(64pi²D²)<0.

The complete tree value is positive but extremely small: 4g²/[D²(3D+2)]. The absolute loop/tree ratio exceeds 10^190. This demonstrates a mismatch in a finely cancelled observable in this particular unadjusted scheme, not strong coupling or a UV no-go.

The new prescription V2S-T1-OS4 sets deltaC_fin=-A1_base(s0,s0,s0)>0. It fixes ONE first-loop value, not derivatives or an entire physical-angle amplitude. An independent complete upper bound gives

0<deltaC_fin/(24q)<14793g²(D+2)/[768(D-1)]<10^-6.

Thus the CONTACT-ONLY classical comparison retains more than (1-10^-6)q in its completed-square quartic and remains coercive with a unique origin. This does not establish the full quantum effective potential or the bare vacuum after all one-point and quadratic counterterms.

A 260-digit diagnostic finds A1_base approximately -1.00470478309049137*10^-409 and deltaC_fin/(24q) approximately 6.43011061*10^-8. Numerical quadrature is a check on the written integral bounds, not their certificate.

## Validation and completed publication gate

The frozen package contains 17 source inputs and 20 report fields, 55 named identities, 55 scalar entries, 25 proof gates, nine controls and 270 rejected unsupported inputs. It adds one matching record for 91 records; all nine original primitive statuses remain unchanged.

Private science first passed 389 tests in 7.52 seconds. Expanded-versus-factored comparison assertions were normalized before the first science run. A private orchestration guard caught a wrong temporary path before command execution; correcting that path did not change the scientific expressions. There was no failed scientific run. Final lint and format checks passed. Final preflight verified the exact 17-source/20-field manifest and all counts in 0.44 seconds, with 389 science tests passing in 6.02 seconds. Repository science passed 389 in 6.13 seconds.

The 50454-character native report was transferred in five checked chunks. All source hashes, the report bytes and all 20 certificate fields were checked independently. Report SHA:

`1277412d19d9b022eb81c3728251d1836d88e4b047c8ad80794107c2aff7c911`.

Fresh ordinary original-SymPy replay passed 414 tests in 2464.58 seconds. Standalone CLI replay passed. Full P8 regression passed 49002 tests in 4572.43 seconds, exit code 0. The 723-file snapshot SHA is

`e0665969a594297d2ec4d31f6775e2869c87985d287cd0b13742d83f035f639d`.

Only full regression used the audited exact-GCD adapter and frozen S219 helper-directory allowance. All 128 original tuple self-checks passed; final counters were 36163 domain fallbacks, 7340 exact descents and 94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy. The collection used 545 static namespace ancestors. These are exact algebra and written proofs, not FORMALIZED proofs.

Publication is restricted to the exact 21-file manifest, with every frozen source and report hash rechecked. Successors and unrelated P4/P9 changes are excluded. Full angular real remainders and low-coefficient matching are successor work; all-loop errors, physical UV completion, finite-gravity Regge control and the common-parent nonlinear bounce remain open. No user intervention is required.
