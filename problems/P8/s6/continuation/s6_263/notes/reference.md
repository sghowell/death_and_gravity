# One unchanged reference energy and its full time boundary

The reference here is the original fixed coefficient reference from
S221/S253, including current Jc, A and Tcorr, not a newly selected
classical solution, quantum state or vacuum. The literal complete
S254 central Hamiltonian specializes to the same weighted Legendre
Hessian. In the light scalar order (Qb,sigma_M1,Pb,p_M1), the checked
blocks are

    Hqq = a_ref^3 YY,
    Hqp = -source^T M^(-T),
    Hpp = M^(-1) XX M^(-T)/a_ref^3.

The -H_ref Qb Pb term from the original time-dependent canonical
map is present in this identity. In the weighted reduced Lagrangian,
write B=S+A_skew, where S=(B+B^T)/2 and A_skew=(B-B^T)/2.
The fixed reference momentum split is

    p_old = Pi + a_ref^3 S y,
    Pi = a_ref^3 (K ydot + A_skew y).

This is a symplectic momentum shear. S is zero at the bounce, but
its full time derivative is not zero there. No derivative is removed
by setting its instantaneous value to zero.

The eight entries of S/q and Sdot/q are bounded by complete
polynomial majorants after multiplying by E^4 Delta^2. The original
reference domain gives |E|>=1/4 and Delta>=1/8; q=P^2/a_ref^2.
All original reference jets and the nonzero A,Tcorr enter the
polynomials. Their exact term counts, reconstruction residuals and
entry bounds are recorded. Both entry families are below10^40.
Using a_ref<=2 and |H_ref|<=2 gives

    |(a_ref^3 S/P^2)'|
      <=2*(max|Sdot/q|+6 max|S/q|)<10^40,
    |a_ref^3 S/P^2|<10^-20

on the stated symmetric slab. This fixed shear is applied to both
the reference and the nearby trajectory. Its coefficients are never
reselected from the latter.

The original scalar comparison energy is

    E_s = (ydot^T K ydot + q y^T G y)/2.

K and G are not assumed to commute. Their original operator
bounds are10^-8<=K,G<=10^4. In the weighted clean coordinates
xq=P y, xp=Pi,

    ydot = K^(-1)(xp/a_ref^3-A_skew xq/P).

The original gyro bound is |A_skew|<=2*10^18, so its scaled value
is at most2*10^18/P_min. Young inequalities give the displayed
positive lower and upper constants, safely enclosed by
10^-12<=Q_s<=10^12 when E_s=x^T Q_s x/2.
The complete nonstationary S221 energy identity yields the safe
root-energy rate10^28. Its nonsymmetric lower matrix and skew
velocity term have not been dropped.

Add the positive energy of the other SIX oscillator pairs: one
heavy scalar, one longitudinal Proca, two tensors and two transverse
Proca. Their normalized q/p coefficients are, respectively,

    ((a P^2+a^3 mH^2)/(P^2+mH^2), a^-3),
    (a, a^-3+1/(a zeta P^2)),
    (a, a^-3),
    (a^-1+a/(zeta P^2), a^-1).

On1<=a_ref<=2 each is between1/8 and8 plus the shown tiny
positive mass term. The actual heavy mass is kept. For any pair
H=(A q^2+B p^2)/2, its exact derivative along Hamilton's equations
is (Adot q^2+Bdot p^2)/2: the mass frequency cancels, rather than
being bounded as a raw perturbative Lipschitz constant. These
reference rates are much smaller than the safe scalar rate.

The direct sum is one complete16-phase positive reference energy,
with eigenvalues in[10^-12,10^12] and root-energy rate10^28.
This is a mathematical norm for a classical linear equation, not
an interacting quantum mean or an ordering prescription.
