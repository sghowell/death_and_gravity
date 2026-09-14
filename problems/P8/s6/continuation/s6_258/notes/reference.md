# Reference TT dictionary and mandatory nonlinear contacts

Linearize gamma=a_hat^2 exp(2v)exp(t), trace t=0. The density Q=exp(-t)
is independent of v, and chi=-div t at first order. The existing
TT-plus-volume reference therefore satisfies this spatial gauge
without a new Gaussian state selection.

For a general tracefree symmetric Fourier perturbation h at k!=0
define P=I-kk^T/k^2 and

    xi = i M0(k)^-1 h k,
    delta h = i(k xi^T + xi k^T - (2/3)I k.xi),
    h_TT = P h P - P Tr(P h)/2.

Direct substitution gives h+delta h=h_TT, h_TT k=0 and Tr h_TT=0.
The same coordinate transformation necessarily changes scalar volume:

    delta v = i k.xi/3 = -(k^T h k)/(4k^2).

That contact vanishes for an already TT h, but not for a generic
tracefree h. A change of shape variables does not authorize discarding
the volume component of the physical field map.

## Two distinct waves at second order

Let ell1=e1, ell2=e2 and take A_yz=A_zy=1, B_xz=B_zx=1, all other
entries zero. Both A and B are tracefree and transverse to their own
wavevector. In

    t = s A cos x1 + r B cos x2,
    Q = exp(-t),

the mixed derivative in s,r has Fourier coefficient at ell=e1+e2

    Q12 = (AB+BA)/8 = (E_xy+E_yx)/8.

Its divergence coefficient is i Q12 ell != 0. Thus the naive
two-wave family is not in the nonlinear gauge, although both
individual first-order waves are TT. The full off-gauge ghost
block and the gauge-surface-only block differ at this coefficient.

At the homogeneous reference the required mixed coordinate correction is

    xi12 = -i M0(ell)^-1 Q12 ell,
    delta h12 = i(ell xi12^T+xi12 ell^T-(2/3)I ell.xi12),
    Q12_restored = Q12 - delta h12.

Then Q12_restored ell=0 exactly. The corrected full ghost block equals
the surface formula only after this transformation. Its scalar-volume
contact is

    delta v12 = i ell.xi12/3 = 1/32.

These are mixed derivatives/Fourier coefficients with the above cosine
normalization; no extra factorial is silently inserted. Independent
full matrix-exponential quadrature extracts the mixed Fourier
coefficient before applying the correction.

The full local slice exists by the implicit-function argument; this
second-order construction illustrates a term in its expansion rather
than claiming a globally convergent coordinate series. Keep its
coordinate change in the entire physical source pullback, canonical
momenta and primitive temporal boundary. On an off-shell reference,
second field-map variations multiply one-point terms. The parent
already showed why these contacts cannot be erased by identifying the
fixed quantum mean with an unforced classical reduction.
