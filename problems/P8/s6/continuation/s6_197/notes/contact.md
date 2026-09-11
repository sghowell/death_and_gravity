# Complete second-metric contact difference

S195's exact spatial Hamiltonian gives the second vertex
with H=(D Gamma+Gamma D)/2 and POSITIVE mass sign.
Its local current term is -<H_DGamma>. Keep the full
noncommuting product, Fourier convolution and all3 modes.

The ten-row energy feature vector is the physical readout
multiplied by a^(3/2), up to a fixed permutation and signs.
Its contact quadratic matrix is diag(H,H,H,0), with
operator norm at most||D||F||Gamma||F. The last row is the
retained temporal-constraint readout: zero direct shear
vertex does not mean that a polarization has been removed.

For each mode, v=vref+e, and the exact covariance difference is

    vv^dagger-vref vref^dagger
    =vref e^dagger+e vref^dagger+ee^dagger.

The S186 bounds give ||vref||<=2F sqrt(nu),
||e||<=F R nu^(-11/2), F=4000, R<2e6.
The Hamiltonian factor1/2 and all3 polarizations bound the
complete difference by

    (3/2)F^2[4R nu^-5+R^2 nu^-11] ||D||F||Gamma||F.

Because R/nu^6<1, its coefficient is
<(15/2)F^2(2e6)nu^-5=2.4e14 nu^-5.
The exact radial integral is

    J5=integral nu^-5 d^3k/(2pi)^3=A^3/(6pi^2 m^2).

Using pi>3 yields the complete contact pairing bound
<1e8 ||D||L2||Gamma||L2. This is a full one-momentum
integral, not a two-created-mode overlap.

Its discarded tail satisfies
J5_tail(K)<=A^5/(4pi^2 K^2), giving error
<1e14 ||D||L2||Gamma||L2/K^2.
The contact can survive at high external transfer even when
a two-mode propagation overlap is empty. It is neither
dropped nor incorporated by changing a counterterm.
