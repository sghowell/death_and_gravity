# Relational matter observable and reconstruction

Use the local background-dependent scalar

    O=chi-chi_background(phi).

It vanishes on the actual nearby homogeneous solution. At linear
order its perturbation is
delta chi-(chi_background'/phi_background')delta phi.
The original clock derivative is1, and in unitary gauge this is
exactly the reconstructed chi coordinate. The linear observable
is invariant under the time-gauge shift; homogeneity also makes
the residual spatial-gauge variation vanish.

This is a linear relational scalar on a fixed classical
background, not a statement about nonperturbatively dressed local
observables in quantum gravity. A formal covariant test source
sqrt(-g) J O gives matter force J and clock force
-J chi_background'/phi_background'. Their background Noether
combination cancels. Because O_background=0, there is no
quadratic source-times-metric term. This source identity does
not itself construct a compact retarded preparation.

## Nonzero clock pole

In S6.89's triangular characteristic basis chi=sigma-ell b,
the principal chi response is exactly

    ell^2/(kappa_c omega^2-q G_c)
      +1/(kappa_m omega^2-q G_m).

Equivalently the clock pole weight in the normalized frequency
denominator is ell^2/kappa_c. The actual real bounds give

    ell^2/kappa_c >
    [9/1000]/(1+10^-6)^2 *4(.24)^2/3
    =691200000/1000002000001 >1/2000.

The distinct luminal pole cannot cancel this nonzero clock pole.
This is a classical principal response, not a fundamental
quantum spectral-positivity or UV argument.

## Field reconstruction error

For each Fourier mode the exact observable is

    O_hat=R^(-3/2) Y_2/k.

For either outgoing branch put c_j=R^(-3/2) S_(2,j).
The complex annuli and explicit basis give the safe real-time
bounds 10^-4<|c_j|<1000. For the clock, |ell|>1/25 and
sqrt(2|omega_c kappa_c|)<40 suffice before the volume factor;
the matter projection is also nonzero.

Take leading modal amplitude a(k) and pure branch j.
If ||U-U_lead||_infinity<=epsilon |a(k)|, then the pointwise
Fourier error relative to the leading O amplitude is bounded by

    [4*1000/(10^-4)] epsilon=4e7 epsilon.

The factor k cancels between the exact and leading observable.
This cancellation must precede comparing the field norms: the
scaled phase Y and the physical relational field have different
k powers.

For k in [K/2,2K], packets.md gives
epsilon<=(1+2e21/K)(4T*1e38/K)+2e21/K.
The resulting relative scalar bound is below 10^44/K.
It is uniform over the whole band. Multiplying by any prescribed
Fourier envelope and using Plancherel therefore gives the same
relative L2 bound. The two conjugate carrier supports are disjoint,
so forming real packets preserves the estimate.

One can take a(k)=k O_initial_hat/c_j(initial) to normalize the
leading initial relational field. The included near-identity
correction changes that field by the same controlled error; it
is not discarded or silently called zero. Finite spectral support
makes the exact initial phase and all reconstructed linearized
fields smooth with finite Sobolev norms of every order. Scaling
their amplitude is allowed in the linear theory, but does not
prove a nearby exact nonlinear solution or quantum resolution.
