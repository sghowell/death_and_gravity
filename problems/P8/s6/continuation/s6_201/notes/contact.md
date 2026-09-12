# Complete one-mode contact correction

The local second-metric vertex contains
H=(D Gamma+Gamma D)/2, with positive mass sign. In Hamiltonian
energy-feature order its full matrix is diag(H,H,H,0), not one
selected polarization or trace. Its norm is bounded by
||D||F||Gamma||F. The final zero block records the direct unimodular
constraint vertex, not omission of the constrained physical readout.

For every mode the alpha-weighted covariance differs from the unit
one by exactly b_k v_k v_k^dagger. The norm bound
||v_k||<=4000 sqrt(nu) therefore gives, after the half-Hamiltonian
factor and sum over all three modes,

    (3/2)4000^2 B^2 nu^-11 ||D||F||Gamma||F.

The complete massive radial integral is

    J11=8 Amax^3/(315pi^2 m^8),

using integral_0^infinity y^2/(1+y^2)^(11/2) dy=16/315.
The resulting contact pairing is below
1e-3||D||L2||Gamma||L2 for nonzero norms. It retains the
noncommuting product, longitudinal contribution and spatial convolution.

The removed contact band is ONE internal momentum. Its tail obeys

    J11_tail(K)<=Amax^11/(16pi^2 K^8)
                 <Amax^11/(144K^8).

With the complete mode coefficient, K>=m makes the contact tail
less than ||D||L2||Gamma||L2/K. This contact remains present even
when external transfer exceeds twice the memory band.
