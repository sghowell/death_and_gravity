# Actual three-polarization scalar-force covariance

Use the positive-frequency Minkowski vacuum already specified
for the new S6.182 action, not the CD Cauchy state transported
to a different background. The physical metric is held fixed.
Its canonical vector has mass m=1000 and the same three-mode
phase measure and covariant prescription. Turning on the
prepared classical source adds a smooth coherent mean,
but leaves the connected covariance unchanged.

The fluctuating part of the unintegrated scalar force paired
with a real test eta is exactly

    F_eta-<F_eta>=-sqrt(kappa)m integral A_centered^mu DS_mu[eta].

The deterministic source contact and fixed scalar coefficient
have no connected contribution to this LINEAR field smearing.
This is not a coincident quadratic stress-tensor observable.

For momentum p=(omega,q,0,0), omega=sqrt(m^2+q^2), take physical
polarizations (q/m,omega/m,0,0),(0,0,1,0),(0,0,0,1).
Their Euclidean matrix sum P=-eta+p p^T/m^2 is positive
semidefinite, with eigenvalues0,1,1,1+2q^2/m^2.
Thus P<= (1+2|k|^2/m^2)I at arbitrary spatial momentum
by spatial rotation. The temporal component and longitudinal
polarization are included. Minkowski index sign changes
are orthogonal for the Euclidean norm used in this bound.

For any real smooth spatially Schwartz test four-vector f
compact in a time interval of length T, the complete
positive-frequency covariance is its mass-shell Fourier
quadratic form with measure d^3k/[(2pi)^3 2omega].
Cauchy-Schwarz in time and Plancherel in space give

    C_A(f,f)<=T/(2m)||f||L2(dt dx)^2
             +T/m^3||grad_spatial f||L2(dt dx)^2.

Indeed 1/(2omega)<=1/(2m), and
|k|^2/(omega m^2)<=|k|^2/m^3. This integrates the entire
mass shell; no momentum cutoff or finite-mode sum replaces
the continuum inequality. The norms are finite for the
stated smooth test class, so no coincident subtraction or
additional noise counterterm is needed for this smearing.

Set f=DS_eta. The full source estimates give
||f||L2<=C||eta||J3 and ||grad f||L2<=62C||eta||J3.
For T<=1,

    Var(F_eta)<=kappa C^2[m/2+3844/m] ||eta||J3^2.

At kappa0 this coefficient is below4e-2380; the standard
deviation is below2e-1190||eta||J3. Variance scales as
kappa0/kappa and standard deviation as its square root.

This is an exact conditional Gaussian connected force
bound at fixed scalar history and fixed Minkowski metric.
It does not bound metric stress noise, a stochastic
metric solution or interacting light/gravity loops.
