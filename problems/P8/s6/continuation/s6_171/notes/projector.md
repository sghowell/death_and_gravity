# Complete twenty-frame physical projector remainder

Use exactly the actual CD metric and the exact free in/out Hadamard
states of S6.170. At fixed real t set

    a=(1+t^2)^2, q=p/a, L=sqrt(t^2+tau^2), D=1+t^2,
    E=sqrt(q^2+m0^2), m0=.99m, K20=163840,
    B=Delta/tau*w(t/tau)+mL/D, w(x)=(1+x^8)^(-9/8),
    A=16qB/E^3, r=K20/(LE).

The inherited complex disks and actual gap give A<1/4,r<1/2.
Every real exact frequency e_j>=omega>=E, and
|theta_j|=|atan(g_j/e_j)|<=A r^j. There are twenty rotations j=0,...,19,
alternating x and y axes. Their exact product maps the positive z vector
to n20; the physical negative projector is (1-n20.sigma)/2.
It is a comparison projector, not an exact evolving or Hadamard state.

For the internal j>=1 product, its spherical distance from z is at most
the sum of absolute angles. Thus |v_z-1|<=2A^2r^2.
Only even j>=2 x rotations change v_y, and each changes that coordinate
by at most its absolute angle: |v_y|<=Ar^2/(1-r^2)<=4Ar^2/3.
The outer j0 rotation therefore gives

    |n20_z-cos(theta0)|<=4A^2r^2.

For pressure, separate j1, since j0 leaves x unchanged. In the j>=2
product, |v_x|<=Ar^3/(1-r^2), |v_z-1|<=2A^2r^4.
Writing q0=g0/omega and b=dot(q0)/(2omega), the exact next ratio is

    q1=b/(1+q0^2)^(3/2).

Cauchy's inequality on the initial disk gives
|b|<1024 A/(LE)<=Ar: the complex |q0| bound is (10/9)A,
the disk radius is L/1024 and the final real division is by2omega.
Combine the internal geometric series with
|sin(atan q1)-q1|<=|q1|^3/2 and
|q1-b|<=3|b|q0^2/2. A conservative combined bound is

    |n20_x-b|<=4Ar^3+8A^3r.

Both helicities give physical energy -2omega*n20_z and pressure

    -2[q^2 n20_z+qM n20_x]/(3omega).

These are not the rotating-frame evolution generators. The local
subtractions are rho0=-2omega, rho2=omega*q0^2 and
P0=-2q^2/(3omega), P2=q^2*q0^2/(3omega)-2qM*b/(3omega).
The complete projector energy remainder is bounded by

    16E A^2r^2 + (3/2)E A^4.

The first term compares twenty and one frame; the second follows from
the exact Taylor remainder of (1+q0^2)^(-1/2). Pressure's diagonal
remainder is at most one third this value. M<2m bounds the remaining
two-helicity off-diagonal factor by4m/3.

The S6.170 full half-line Duhamel estimate compares the same exact
state to P20. Its already integrated energy allowance is below1e-1090;
the pressure allowance is at most one sixth that generous bound.
The full-line in/out bound was not substituted for an occupation-only
energy expression: linear coherence is retained.

The physical measure is N q^2 dq/(2pi^2), so the scale factor cancels
when changing comoving to physical radial momentum. Use

    integral q^4/E^7 dq=1/(5m0^2),
    integral q^6/E^11 dq=2/(63m0^4),
    integral q^3/E^6 dq=1/(4m0^2),
    integral q^5/E^10 dq=1/(24m0^4).

All integrals cover [0,infinity). Also B/L<=Delta/tau^2+m,
B<=Delta/tau+m, B/L^3<=Delta/tau^4+m/tau^2.
The resulting uniform rational bounds are implemented in bloch.py,
using pi^2>9. All42 copies are retained, with the same conservative
Delta upper bound even for their constant-mass members.
