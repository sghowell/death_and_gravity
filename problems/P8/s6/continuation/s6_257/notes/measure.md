# Full second-class density at a fixed finite regulator

Take the canonical system after the complete boundary transformations.
At a fixed finite regulator, collect all lapse and normal-vector
auxiliaries in a, with primary momenta p_a and secondary equations
C=partial_a H. On a single regular root branch set D=partial_a C.
Keep every cell and field dependence of D. The secondary bracket
E={C,C} is generally nonzero and may contain spatial couplings.

For constraints (p_a,C), their complete matrix is

    M=[[0,-D],[D^T,E]],
    M^-1=[[D^-T E D^-1,D^-T],[-D^-1,0]],
    det M=(det D)^2.

Block multiplication proves the inverse at arbitrary finite dimension.
To compute the determinant, interchange the two block rows; the resulting
block-triangular determinant and row-permutation sign give the displayed
square. No assumption E=0 or independent-cell factorization is required.

The positive second-class density is sqrt(det M)=|det D|. Integrating
delta(p_a) and then delta(C) on one regular local root branch contributes
1/|det D|. These cancel exactly at the same finite regulator, leaving
the retained canonical Liouville measure. Multiple roots require their
separate branch contributions and cannot be identified by this local
statement. A zero pivot is outside the theorem.

Retained canonical functions commute with p_a. Their two bracket vectors
therefore have zero primary slots; the lower-right block of M^-1 is zero.
The second-class Dirac correction vanishes even with nonzero E. This
leaves 11 dynamic coordinates and 22 phases: metric six, vector three,
and matter two. At the classical level three regular spatial gauges and
their three first-class constraints remove six phases, leaving the
known eight physical modes. This count does not compute a quantum gauge
determinant or prove its preservation by a finite regulator.

Keep the spatial constraint delta functions, gauge conditions and their
first-class determinant as a separate unevaluated factor. The finite
summation-by-parts diagnostic is not a discrete Leibniz rule or a proof
of first-class diffeomorphism closure. Affine-complement/projective-gauge
quantization is a further separate task. A remaining nonlinear momentum
integral is not generally Gaussian; no product of simple square roots
of its velocity Hessian is asserted.

Perform the cancellation before a continuum limit or separate functional
determinant regularization. The result does not set every quantum
determinant, ordering term, local counterterm or anomaly to zero.

The finite canonical distinction is consistent with
[Bellorin and Droguett, section 5.1](https://arxiv.org/pdf/1912.06749).
Their theory-specific operator is not imported into P8. The need to
check a regular Hamiltonian in unitary gauge is discussed by
[Langlois and Noui, section VI](https://arxiv.org/pdf/1512.06820).
Here the source-pinned full pivots, not an external model, supply that
check. Relevant source pages were rendered and visually inspected.

