# Full finite local reference inverse with all local terms retained

Keep the original flat scalar factorsFdiag and, optionally, the S226 bounded channel matrixV of norm at mostM. At the same weight sigma_M, the original causal inverseKdiag has norm at most5/(32+13M), while the newly controlled actual local conjugate has norm below1/4. Therefore

||Kdiag(V+Omega)||
<=5(M+1/4)/(32+13M)<5/13<1/2.

The operatorV+Omega is causal and bounded in the same weighted space. Omega is generally a NONLOCAL coordinate conjugate of a local differential operator; it is not being passed off as time multiplication. The bounded-causal Neumann argument applies directly to this proved operator bound, with coefficient order intact.

Define

Kfinite=(I+Kdiag(V+Omega))^-1 Kdiag.

Its norm is at most

1/[d_M-M-1/4]=20/(123+32M).

As in S226, the bounded equation and the original forward distribution prove both graph inverse identities forFdiag+V+Omega. No instantaneous term, shear pole, initial atom or original finite coefficient is changed.

The complete finite local reference is exactly

Afinite=B*(Fdiag+V)B+Rloc=B*(Fdiag+V+Omega)B.

Its inverse isEfinite=Y Kfinite Z. Both inverse products cancel in the stated order and retain the original completeRloc. The resulting norm is

||Efinite||<=125/[(123+32M)(sigma_M-20)^4],

uniformly in all comoving momenta on weightedL2_tH^r for every realr, on the unchanged conformal slabT<=1.

For the physical force reference, output normalization is1/(64pi^2 kappa a^4). Its inverse isEfinite M_[64pi^2 kappa a^4], on the RIGHT. The actuala^4<6 gives the bound

48000pi^2 kappa/[(123+32M)(sigma_M-20)^4].

Unweighting on[0,T] costs exp(sigma_M T). The amplitude-to-metric norm factors inherited from S224 are also retained if physical Frobenius conventions are used. The displayed number is not an unweighted physical-feedback estimate or stability bound.

## The precise compatible composed graph

Write H=weightedL2_tH^r and Fmiddle=Fdiag+V+Omega. Its middle graph is

Dom(Fmiddle)={x in H : Fmiddle x is an ordinary element of H as a causal distribution}.

For the full reference, use the explicit pullback graph

Dom(Afinite)={s=Yx : x in Dom(Fmiddle), y=Fmiddle x in H, B*y in H},
Afinite s=B*y.

All equalities here are causal distribution equalities including the complete initial boundary. Equivalently, this graph requires Bs=x in H, the middle output y in H, and B*y in H. This requirement matters: for merely measurable V, the expression VBs need not be defined for an arbitrary distribution Bs. We do NOT enlarge the domain to every s in H whose informal complete operator expression appears to be in H.

For any source f in H, put y=Zf and x=Kfinite y. Then y is in H, B*y=f, x is in the middle graph, Fmiddle x=y, and s=Yx is in the stated full graph. This proves the right inverse. Conversely, for s in the graph with output f, causal uniqueness for B* gives y=Zf. The middle graph inverse gives x=Kfinite y, and s=Yx gives the left inverse and uniqueness. Thus the displayed inverse is onto precisely this compatible graph, and its graph norm ||s||_H+||Afinite s||_H is controlled by the displayed inverse bound plus1.

This is not a bounded forward self-map. Maximality among every possible distributional realization, graph density, and invariance of the full S222 quantum-force graph are not established. The explicit graph also makes precise the corresponding curvature-composition scope inherited from S225/S226; their frozen inputs are not altered.
