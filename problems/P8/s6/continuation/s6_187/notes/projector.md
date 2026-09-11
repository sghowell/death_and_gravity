# TT projection and a genuinely nonempty compact class

For k nonzero let n=k/|k| and P=I-n n^T. On symmetric
matrices the spatial Fourier TT projector is

    Pi(A)=P A P-(1/2)P tr(PA).

It is self-adjoint and idempotent in the Frobenius inner
product: the transverse two-dimensional block is projected
orthogonally onto its traceless part. Its operator norm
is1. It commutes with all time/spatial derivatives and
with the actual scalar-coefficient wave operator L.
A bounded representative at k=0 changes no stated L2 norm.
The energy proof itself has no division by k.

General TT projection does NOT preserve compact spatial
support. Where it is used, the S6.186 stress smearing
extends by continuity in
N(f)^2=sum_(j=0..3)||partial_t^j f||L2^2+||grad f||L2^2.
Smooth compact tests are dense in the corresponding
time-compact mixed Sobolev space; the bounded Fourier
projector and its commuting derivatives preserve this
norm. No false compact-support assertion is required.

There is also a local construction of actual compact TT
tests. Put Q_ij=-Delta delta_ij+partial_i partial_j, with
polynomial Fourier symbol |k|^2 delta_ij-k_i k_j, and set

    psi_ij=sum_kl[Q_ik Q_jl-(1/2)Q_ij Q_kl] A_kl

for any real compact smooth symmetric tensor potential A.
The derivatives commute. k^T Q=0, Q^2=|k|^2 Q and
tr Q=2|k|^2 show that psi is symmetric, transverse and
traceless. Local differentiation preserves compact support.
The operator is not identically zero: for the A11 symbol
and k along z its output is diag(k_z^4/2,-k_z^4/2,0).
A nonzero differential operator acts nontrivially on some
compact smooth potential, as can also be seen by a bump
equal to a suitable polynomial near an interior point.

Set q=L psi. Then q is compact TT and has both required
weighted moments by integration by parts. A nonzero
compact psi cannot have L psi=0, because its vanishing
data before its support imply psi=0 by wave uniqueness.
Hence this construction gives genuinely nonzero detectors.
It is not a restriction to homogeneous spatial momentum.
