# Full tensor projection and normalization

At total timelike momentum P with s=P^2, set
theta_ab=eta_ab-P_a P_b/s,

    P0_abcd=theta_ab theta_cd/3,
    P2_abcd=(theta_ac theta_bd+theta_ad theta_bc)/2-P0_abcd.

The full stress-pair cut is conserved and symmetric. Rotation averaging
in the center-of-mass frame therefore gives a P2+b P0. This follows
from the decomposition of symmetric spatial tensors into trace and
traceless parts. Their ranks are 1 and 5, so, summing all nine pairs,

    b = sum |tr Tij|^2/3,
    a = (sum ||Tij||_F^2-sum |tr Tij|^2/3)/5.

The displayed COM amplitudes are real after the two creation phases
are multiplied, so ordinary squares equal absolute squares. Exact
evaluation gives

    a=2(13k^4+40k^2 m^2+30m^4)/15,
    b=(4k^4+4k^2 m^2+3m^4)/3.

Substitute k^2=s/4-m^2. This yields
a=(13s^2+56m^2s+48m^4)/120 and
b=(s^2-4m^2s+12m^4)/12.
Independent full 9-by-9 angular quadrature checks the entire tensor,
not only its trace or one transverse channel.

The two-body measure is
d^3k/((2pi)^3 2E_k) d^3l/((2pi)^3 2E_l)
times (2pi)^4 delta4(P-k-l).
Radial integration of delta(sqrt(s)-2E) gives beta/(8pi);
the remaining angular average is the one just evaluated.
Connected Wick contraction supplies exactly 2. Thus the positive-energy
stress Wightman cut is 2 beta(a P2+b P0)/(8pi).

If sigma_T is the Kallen-Lehmann density, W=2pi sigma_T at positive
energy. The actual current has the additional factor 1/4. Hence
rho=W/(8pi), giving the two densities in FORMULATION.md.
An independent oscillator has
i<[q(t)^2,q(0)^2]>/4=sin(2Et)/(4E^2) and positive Laplace response
1/(2E(4E^2+u^2)). This fixes the sign and spectral Jacobian.

With Fourier convention exp(i P.x), z=(omega+i0)^2-|p|^2,
the retarded spectral denominator is (sigma-z)^-1; at positive energy
Im chi(s+i0)=pi rho(s). Theta support in the two-body measure sets both
densities to zero through 4m^2. For x=m^2/s in [0,1/4],
1-4x+12x^2=12(x-1/6)^2+2/3>0, and the spin2 polynomial is positive.
These are full external-metric Gaussian cut densities, not full parent
interacting scattering cuts.
