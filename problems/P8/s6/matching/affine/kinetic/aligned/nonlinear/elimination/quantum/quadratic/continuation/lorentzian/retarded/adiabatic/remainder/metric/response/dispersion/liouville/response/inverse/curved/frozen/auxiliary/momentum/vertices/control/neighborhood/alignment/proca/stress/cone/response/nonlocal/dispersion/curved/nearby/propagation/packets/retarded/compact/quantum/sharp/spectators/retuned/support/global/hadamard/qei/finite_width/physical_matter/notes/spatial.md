# Complete quadratic spatial reconstruction and its domain

Use the same spatial canonical gauge as S6.76-77,

    hat_g=I+2v I+t,

where t and its canonical momentum are symmetric, tracefree and
transverse. This is the LINEAR spatial metric coordinate, not
exp(2v)I. Its inverse, volume, curvature, connection and momentum
constraints are reconstructed before perturbative truncation.

The canonical phase has fourteen channels in this order:

    (v,p,chi,Pchi,T,PT,S,PS,Wx,Wy,Wz,Px,Py,Pz).

At wave direction x the two tensor bases are diag(0,1,-1) and
E_yz=E_zy=1. Their canonical momentum bases are half of these,
since E:E=2. The actual scalar and vector spatial generators,
including the Proca divergence term, enter the York recursion.
For proper-subset momenta all required quadratic constraint
components are solved. The validated inputs have zero background
and only linear phase fields; malformed gauges are rejected.

At zero total output a two-leg context solves first-order York.
The missing homogeneous second-order York correction is flat
tracefree. It cannot enter dp2 through its flat trace, and it
cannot enter the other invariants until higher order. Therefore
it cannot change the quadratic lapse force or matter observable.
This is not an inversion of a nonzero homogeneous source.

The matrices are defined by external momenta (k,0,0),(-k,0,0),
k>0, with rho2=(1/2) Q^T M Q. They are the zero-output coefficients
of the POINTWISE scalar normal density. They are not multiplied
by a physical or hatted volume measure. A local bilinear kernel at
nonzero output, or a perturbed-volume integral, is different.

Three-leg fixtures use (k,0,0),(0,k,0),(-k,-k,0). The coefficients
with masks 3,5,6 have nonzero output; the appropriate order-two
York correction is included. All fields, both tensor directions
and vector sectors are mixed with exact rational amplitudes,
with TT tensors transverse to their own wavevector. Checks run
at all three exact times -1/100,0,1/100. In each of the nine
quadratic output coefficients the n2 contribution is nonzero.

The universal invariant reconstruction follows by the written
Taylor/constraint argument, not by extrapolating these fixtures.
The zero-output matrix covers every center phase polarization
for a fixed nonzero wave direction; isotropic background rotation
transports the basis. No zero-frequency gauge inverse or locality
of the reduced constraint variables is inferred.
