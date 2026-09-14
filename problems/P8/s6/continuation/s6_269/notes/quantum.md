# Evaluated same-seed coherent tail and two-volume comparison

Keep the exact pure original seed psi0 and its full finite
unit-CCR covariance V. Whiten z=S w with V=SS^T/2.
The S268 normalized coherent resolution is

    Q_V(a)=(2pi)^(-d) integral a(z)
                    |D(z)psi0><D(z)psi0| dz.

Its analysis map is an isometry. In particular Q_V is
positive and unital, ||Q_V(a)||<=||a||infinity, and
Q_V(a)*Q_V(a)<=Q_V(|a|^2). No finite occupation cutoff
or finite-dimensional representation of the CCR is used.

The probability density obtained by this POVM in the SAME
seed has covariance2V. In w it is the standard Gaussian
on2d real phase coordinates, not covarianceI/2.
Consequently t=||w||^2/2 has the gamma(d,1) distribution.
This probability law belongs to a coherent POVM. It is
not a joint sharp measurement of noncommuting q and p.

At the actual d=48 and core R=10^20,
t0=R^2/2=5*10^39. The exact positive outside mass is

    T48(t0)=exp(-t0) sum_(k=0)^47 t0^k/k! >0.

The code checks its COMPLETE polynomial derivative
after factoring out the exponential. It does not
equate floating underflow with zero probability.
The independent160-digit arbitrary-exponent mpmath
test also verifies the actual strictly positive tail.

For t0>d the gamma Chernoff estimate gives

    T_d(t0)<=exp[-t0+d+d log(t0/d)].

The first five positive terms in the exponential
series give e>8/3. Exact rational exponentiation gives
(8/3)^100>10^40, while t0/d<10^40.
Thus log(t0/d)<100 and the exponent above is less
than-10^39. We obtain the evaluated strict bound

    0<T48(t0)<exp(-10^39).

Smooth cutoffs equal1 on the core and supported
strictly inside radius2*10^20 exist, for example a
radial bump of ||w||^2 with a smooth flat transition.
Translations act by rotations preserving the exact
seed covariance, hence preserve ||w|| and this
cutoff. This is not a state projection: S268's compact
effect has norm strictly below1 and no nonzero state
has exact finite phase support.

The full nonlinear physical symbol F and actual
background F0 satisfy notes/invariants.md uniformly.
For any such 0<=chi<=1,

    Fext=F0+chi(F-F0),  ||Fext-1||infinity<10^-255.

Positivity and unital contractivity imply

    ||Q_V(Fext)-I||<10^-255,
    |<psi0,Q_V(Fext)psi0>-1|<10^-255.

These operator bounds also hold on the reducing
common zero-translation-charge sector. They refer
to the UNCALIBRATED positive coherent volume.
They neither assert positivity of Q_V((1-D)Fext)
nor evaluate that distinct readout's ordering error.

For two admissible smooth cutoffs chi1,chi2 with
the same core, physical F and actual background F0,

    delta F=(chi1-chi2)(F-F0).

This is zero on the core and has sup norm below
2*10^-255 everywhere. The isometry/Schwarz inequality
and the exact POVM mass give

    ||Q_V(delta F)psi0||
       <=2*10^-255 sqrt(T48(t0))
       <2*10^-255 exp(-10^39/2).

The expectation uses the mass rather than its square
root, so

    |<psi0,Q_V(delta F)psi0>|
       <2*10^-255 exp(-10^39).

The same seed is used on both sides; there is no
conditioning, re-minimization or coherent-window
replacement of the initial state. The bound is
evaluated at the BOUNCE ONLY, for two defined finite
volume regulators. Initial POVM concentration does
not imply later-time support or a Hamiltonian
comparison with an undefined singular operator.

No derivative of chi enters this particular
uncalibrated volume norm argument. In contrast,
S268's first-Weyl-calibrated interaction and volume
DO require every derivative of chi and the full
mixed covariance. None of those distinct ordering
or time-evolution bounds is claimed here.
