# Complete covariance and the prescribed lower polynomial

For a flat minimal scalar mode use canonical variables with
G0=omega I2 and vacuum covariance I2/2, omega=sqrt(k^2+nu).
For a homogeneous spatial metric perturbation gamma, the balanced
Hamiltonian variations are

X=[omega^2 tr(gamma)/2-k.gamma.k]/omega,
Y=-omega tr(gamma)/2.

The full Hamiltonian also has an instantaneous second metric vertex;
denote its detector contribution by C. It is retained throughout.

For Laplace frequency lambda the covariance equation is

lambda deltaSigma=omega[J,deltaSigma]+(J deltaG+(J deltaG)^T)/2.

Solving all entries gives U=omega(Y-X)/(lambda^2+4omega^2),
V=-U and cross=lambda(Y-X)/[2(lambda^2+4omega^2)]. The complete
detector current is

omega(YD-XD)(Y-X)/[2(lambda^2+4omega^2)]+C.

The derivative marker multiplies lambda. By uniqueness of this exact
linear recurrence, adiabatic orders0,2,4 are precisely the Taylor
coefficients through p^2, where p=lambda^2. The entire second vertex is
in the zero-order term. Pointwise in internal momentum,

1/(sigma+p)-1/sigma+p/sigma^2-p^2/sigma^3
=-p^3/[sigma^3(sigma+p)].

Thus the actual-minus-fourth-order term is exactly the three-subtracted
nonlocal cut, with no residual lower polynomial. The arbitrary symmetric
homogeneous metric spans both conserved timelike sectors. Lorentz
covariance fixes their invariant continuation.

For H the original scalar subtraction and the entire dimensionally varied
pole counteraction fix the remaining finite heat polynomial. The same
Proca statement was proved in S193/S230 with all physical modes. One can
also read these finite Taylor coefficients from the gapped covariant
determinant: volume, R, and quadratic curvature terms have derivative
orders0,2,4. At quadratic metric order every higher heat coefficient
starts at six derivatives or has at least three curvatures, and hence
cannot change these Taylor coefficients. The mass gap gives analytic
Taylor coefficients below threshold; the covariance identity supplies a
convergent remainder, not a use of an asymptotic heat series as a full
determinant.

This argument identifies the H/Proca prescriptions already specified by
the source. It does NOT choose a Phi finite curvature prescription.
Phi's cut and subtraction remainder have the same form, but its lower
R,Weyl^2,R^2 polynomial remains explicit and source-unmatched.

For each positive density define

D_i(z)=z^3 integral_(4nu)^infinity rho_i(sigma)/
                         [sigma^3(sigma-z)] dsigma.

Since rho=O(sigma^2), the integral is locally dominated on the slit plane.
It is analytic and has zero Taylor coefficients through degree2.
For |z|<=2nu, |sigma-z|>=sigma/2 gives

|D_i|<=2|z|^3 I_i,
|D_i-z^3 I_i|<=2|z|^4 J_i.

These bounds use the entire massive density, not a finite-energy cap.
