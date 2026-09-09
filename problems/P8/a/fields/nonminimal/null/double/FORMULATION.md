# P8-A.22 — a physical double-null scalar inequality

This is a finite, nonoptimal bound for the ACTUAL physical null stress
averaged on a timelike two-plane in four dimensions. It is not a bound
on a single null line, a square of a transversely smeared field, an
effective gravitational stress tensor, or a null singularity theorem.
Original P8 remains open.

Let Phi be one real massless scalar on Minkowski space with signature
(+---), improvement coupling 0<=xi<=1/2, and physical null derivative
D=partial_t+partial_z. Restrict to (t,0,0,z). For any Hadamard target
state, write w=<:Phi^2:> and normal order against the Minkowski vacuum.
Every real f in C_c^infinity(R^2) obeys

    integral dt dz f^2 <:T_DD:>
      >= -hbar/(8pi^2) [(2/3)(1-xi)||f_tt||_2^2
            +(4xi/3)<f_tt,f_tz> +2xi||f_tz||_2^2]
         -2xi integral dt dz w (Df)^2.                    (1)

Here T_DD=(D Phi)^2-xi D^2(Phi^2), including its full improvement.
If w<=Phi_*^2 for a fixed Phi_*^2>=0 on the ENTIRE two-dimensional
sampling region, the last term is bounded below by
-2xi Phi_*^2 ||Df||_2^2. No lower Wick cap or momentum cutoff is assumed.
The range xi<=1/2 is a hypothesis of THIS positive-squares proof.

For xi=1/6, put x_plus=t+z, x_minus=t-z, and choose

    h(s)=sqrt(315/256)(1-s^2)^2 for |s|<=1, zero otherwise,
    h_delta(s)=delta^-1/2 h(s/delta), delta>0,
    f(t,z)=sqrt(2) h_delta_plus(x_plus) h_delta_minus(x_minus).

The sampler is normalized in dt dz. It is a compact H0^2 product,
not C-infinity; (1) extends to it by density from INSIDE the sampling
rectangle. A fixed rational boost gives the explicit consequence

    integral dt dz f^2 <:T_DD:>
      >= -14117 hbar/(1536pi^2 delta_plus^3 delta_minus)
         -4 Phi_*^2/delta_plus^2.                         (2)

Both widths are nonzero. At fixed delta_plus the quantum cost diverges
as delta_minus approaches zero. Thus (2) cannot contradict A.21's
one-sided-state-cap null-line obstruction. Only a one-parameter boost
family is optimized analytically; neither (1) nor (2) is claimed sharp.

For the conformal member ONLY, transport to smooth positive spatially
flat FLRW g=a(eta)^2 eta using proper time dt=a d eta. On the same
coordinate plane use physical volume dvol=a dt dz and the affine
null vector K=a^-2(partial_eta+partial_z). Set

    L=f_tt-3H f_t+(2H^2-2Hdot)f,
    M=f_tz-2H f_z,
    G=f_t+a^-1 f_z-2H f.

In the fixed scalar renormalization prescription of A.20, every
Hadamard target obeys the full physical-stress inequality

    integral dvol f^2 <T_KK>_omega
      >= integral dvol f^2 Tconf_KK
        -hbar/(8pi^2) integral dvol [(5/9)a^-2 L^2
                         +(2/9)a^-3 L M+(1/3)a^-4 M^2]
        -(1/3) integral dvol a^-2 w_physical G^2,           (3)

where w_physical is relative to the actual conformal reference,
and the scalar reference is retained, including its finite coefficient:

    Tconf_KK=hbar/(2880pi^2 a^2) [-4H^2 Hdot
             +beta_S(12Hthird+72Hdot^2+36H Hddot)].         (4)

A one-sided physical Wick upper cap bounds the final cost in (3).
The scale factor, volume and affine normalization in (3) cannot be
replaced by their central values without a NEW error bound. Equation
(2)'s numerical constant is not automatically a proper-width FLRW
constant. Equations (3)-(4) alone are not a single-ray Ricci bound,
an index-form estimate, or a self-consistent semiclassical solution.
