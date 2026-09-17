# Source ownership and coupling conventions

S315 supplies the exact mixed coefficients of the unchanged covariant
action at every finite metric order and the terminating labeled rooted
tree evaluator. S310 supplies an independent frozen lower-multiplicity
implementation. S311 supplies the exact24 conserved-pair coefficients.
All three are imported read-only, with the immediate S315 and direct
S311 native reports rebuilt and hash checked during native acceptance.

The external Phi root is amputated. An internal h edge carries
-trace_reverse(R)/P^2, where trace_reverse(R)=eta*R*eta
-eta*trace(eta*R)/2. The Einstein cubic with the canonical free fields
carries1/sqrt(kappa). These signs and normalizations are the same as
the inherited action, not selected from the desired cancellation.

The complete source includes all matter and Einstein trees. For
four Phi legs, its formal inventory has only contact, heavy-exchange
and pure-gravity sectors. C and g in the inventory are formal markers;
they are not substitutions for the original physical couplings.

The new arbitrary-N statements are algebraic identities for a finite
tree source, plus a conditional current norm inequality. They are not
an all-N S-matrix construction. No UV matching coefficient is chosen,
no hard loop is dropped as negligible, and no quantum-state condition
is supplied by the combinatorics.
