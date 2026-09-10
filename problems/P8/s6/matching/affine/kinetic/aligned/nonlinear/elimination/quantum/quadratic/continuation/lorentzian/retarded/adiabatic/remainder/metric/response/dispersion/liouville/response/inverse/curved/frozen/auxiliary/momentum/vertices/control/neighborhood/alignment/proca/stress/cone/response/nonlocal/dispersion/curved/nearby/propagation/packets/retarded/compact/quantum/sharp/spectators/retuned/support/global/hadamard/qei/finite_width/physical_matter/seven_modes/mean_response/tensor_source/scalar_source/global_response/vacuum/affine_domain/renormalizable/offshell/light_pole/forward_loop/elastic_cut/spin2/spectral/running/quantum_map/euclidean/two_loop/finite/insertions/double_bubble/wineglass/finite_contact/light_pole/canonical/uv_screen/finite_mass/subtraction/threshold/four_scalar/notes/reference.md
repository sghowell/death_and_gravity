# Reference, diagram ownership and canonical field factors

The selected nonlocal effective action consists of the
scalar tree action plus the complete one-loop fermion
determinant, with S6.131's fixed local references.
It has exact Phi parity and no direct H Yukawa.

At one fermion loop, four external Phi legs can attach
directly to the four-vertex determinant term. The
two-vertex term is a self-energy insertion. A connected
scalar tree with four external Phi has only a contact
or an H exchange, never an internal Phi exchange:
each interaction contains an even number of Phi legs.
Thus the two-vertex fermion term enters this amplitude
only through external pole/mass normalization. Higher
determinant vertices cannot give four external legs
without adding a scalar loop or a forbidden odd branch.
An H insertion joining two points of the fermion cycle
requires scalar lines and creates an additional loop.
Gauge exchange, fermion parameter counterterms inside
the box and direct H quantum corrections likewise
first enter later loop orders. Inert fermions supply
only the already-fixed vacuum constant at this order.

The pole mass is fixed to one, so s+t+u=4 remains the
same external-shell relation. With
Phi_can=sqrt(kappa) Phi_reference, every four-point
term scales by kappa^-2. The reference tree amplitude
and the box therefore combine as

    A_selected=(A0+A_F)/kappa^2.

Here the tree quartic is L/kappa^2. If instead one
displays the full zero-momentum potential quartic
(L+v4)/kappa^2 in the tree, that tree is
(A0-v4)/kappa^2, and the remaining box must be
(A_F+v4)/kappa^2. The two expressions are identical.
Adding the entire box to the potential-quartic tree
would double count -v4/kappa^2.

Introduce a formal loop variable h before expanding:

    (A0+h A_F)/(1-h f'(1))^2
      =A0+h[A_F+2f'(1)A0]+O(h^2).

The center coefficient increment is
b2_F+8lambda f'(1), not another LSZ copy on top of
already normalized total canonical vertices.

The rational S6.131 slope interval and the complete
box error give explicit intervals for both the formal
increment and the exactly normalized selected
functional. Generic interval arithmetic handles a
negative lower numerator with the correct denominator
endpoint; small valid examples are allowed to be
inconclusive. At the actual named boundary the lower
increment is positive and the normalized coefficient
lies above 4lambda by less than a relative 10^-204.

Exact constant normalization of a finite selected
functional does not bound all loop orders. In
particular the old scalar calculation uses different
finite references for its light bubble and canonical
vertices. Neither its numerical couplings nor its
complete two-loop budget are silently identified with
full MSbar matching of the present candidate.
