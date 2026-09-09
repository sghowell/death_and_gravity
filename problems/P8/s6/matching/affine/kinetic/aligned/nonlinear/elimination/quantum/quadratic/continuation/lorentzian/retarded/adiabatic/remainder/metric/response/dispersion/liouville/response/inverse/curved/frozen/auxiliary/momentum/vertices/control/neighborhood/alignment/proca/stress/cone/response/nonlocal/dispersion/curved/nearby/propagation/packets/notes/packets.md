# Exact normal-form error and three-dimensional packet localization

Write Y=SU and U=(I+Z/k)V. The exact identity
i[Lambda,Z]+B0=0, including B0's zero diagonal, gives

    V'=ik Lambda V
       +(1/k)(I+Z/k)^-1 [B0 Z+W+W Z/k-Z']V.

The entire Z' term is retained. For k>=4*10^31,
||(I+Z/k)^-1||<=2 and the bracket with this inverse has norm
less than10^38. The leading diagonal phase has norm1 on real time.
Across I=[-T/2,T/2], of length T=10^-7, variation of constants gives

    ||V-V_lead||/||V_initial|| <=exp(T*10^38/k)-1
                                 <=2T*10^38/k.

The last inequality uses T*10^38/k<=1/4 and the complete
exponential series bound. It is not an exponential of the carrier
or the large crude transport norm. Returning to U adds the
near-identity correction. In the carrier band [K/2,2K],

    epsilon_U <=(1+2e21/K)(4T*10^38/K)+2e21/K.

The observable reconstruction in observable.md bounds its
relative field error by 10^44/K.

## Explicit finite-band real data

Take K=10^72 and sigma=10^24. With the unitary Fourier convention,

    f_sigma(x)=sqrt(sigma/pi) sin(sigma x)/(sigma x)

is the inverse transform of the normalized constant on
[-sigma,sigma]. Its L2 norm is1. The real three-dimensional
leading packet is

    sqrt(2) cos(K x1) product_i f_sigma(x_i).

Its two Fourier cubes centered at +/-K e1 have coordinate
half-width sigma and lie in K/2<|k_vector|<2K.
Their disjointness gives the stated real normalization exactly.
The finite-frequency ODE and its initial data define exact
linearized classical solutions mode by mode, not a ray-only
approximation. Homogeneity preserves their finite spectral support.

Choose outgoing -clock or -matter mode for the positive carrier,
with its conjugate for the negative carrier. Let
s_j(u)=integral from the initial slice to u of omega_j(v) dv.
The leading positive-carrier phase is exp(-i|k_vector|s_j).

To compare with a rigid translation along x1, use

    |k_vector|-k1=k_perp^2/(|k_vector|+k1)
                    <=2sigma^2/K

on the positive cube. The real ray length is below2T, so
|exp(i phase_error)-1|<=|phase_error| gives additional relative
L2 error at most4T sigma^2/K=4*10^-31.
It is uniform across the full three-dimensional cube, not a
one-dimensional packet used as a finite-energy three-dimensional
solution. Combining this with the field reconstruction gives
total relative error delta<10^-27.

## Nonzero tails and normalization of the exact field

Across this half interval the parent ray-separation lower bound is

    d=99/(4*10^14).

Choose initial/comparison cubes of half-width r=d/8.
The two-sided sinc bound is
integral outside [-r,r] |f_sigma|^2 <=2/(pi sigma r).
Using pi>2, a three-coordinate union bound and the real-carrier
factor2 gives leading outside-cube L2 mass

    epsilon_tail <=6/(sigma r)<10^-9.

The tails are explicitly nonzero. If f is the translated leading
field and g the exact field, with
||g-f||<=delta ||f||, then

    ||g_outside||^2/||g||^2
      <=(sqrt(epsilon_tail)+delta)^2/(1-delta)^2
      <=8(epsilon_tail+delta^2)<10^-8.

The second inequality uses delta<1/2. The final exact field norm,
not only the leading field norm, appears in the denominator.
The same estimate applies on the initial slice with its retained
near-identity correction.

The actual clock-ray center is more than d ahead of the matter
ray center on the final common clock slice. The final clock cube
has x1 minimum greater than s_matter+d-r, whereas the physical
matter causal future of the initial cube has x1 maximum at most
s_matter+r. Since d>2r the two regions are disjoint.

This is controlled finite-frequency classical packet transport
in the relational field, with tails and errors retained. It is
not strict compact support, compact local preparation, a retarded
commutator theorem, or a signal below an interacting EFT cutoff.
Larger declared decimal carrier exponents divisible by three
between72 and120 improve every error and tail bound with
sigma=K^(1/3); no wider physical EFT domain is inferred.
