# Proper forests and the scalar-mass transfer

Use the frozen S6.139 total counterterms (the factor 1/e is explicit):

    z_psi=-(Y/2+a C_F)/(Qe),
    eta_m=nu_y=(Y-4a C_F)/(Qe).

For the first fermion scalar two-point function f_D, the complete
proper insertion is

    delta f_D=2(nu_y-z_psi)f_D
              +(eta_m-z_psi)m partial_m f_D at fixed mu.

At zero momentum f_D=2NY m^2 exp(gamma_E e)(mu/m)^(2e)
Gamma(e)[4-2/(e-1)]/Q. Its fixed-mu homogeneity is 2-2e,
so the multiplier is (6-3e)(Y-2a C_F)/(Qe).
This gives both c_s,c_g of anchors.md and counts the two
overlapping vertex subtractions only once each. There is no
product forest of overlapping subgraphs.

There is also the proper whole fermion cycle with four scalar
legs. For a massless exchanged scalar its local counterterm
closes into a scaleless tadpole and vanishes in dimensional
regularization. At scalar mass squared b>0 it does NOT vanish.

Let Gamma_0,D(q^2) be the complete zero-external-momentum four-scalar
kernel with its proper quartic MS subtraction, as in S6.144.
The scalar quadratic primitive has contraction
(1/2) Integral_q Gamma_0,D(q^2)/(q^2+b). Consequently,

    f_1,D(0)-f_0,D(0)
       =-1/2 Integral_0^1 db F_mix,D(b),
    F_mix,D(b)=Integral_q Gamma_0,D(q^2)/(q^2+b)^2.

The six labeled boxes and the one-half contraction give the same
three words as the direct quadratic ledger. The fermion mass,
kinetic and Yukawa MS pole coefficients are independent of b;
their insertions in the first fermion bubble have no internal
boson and cancel in this difference. The quartic MS counterterm
does not cancel and is already included in F_mix,D.

In regulated units Q suppressed,
Tad_D(b)=exp(gamma_E e)mu^(2e)Gamma(e)b^(1-e)/(e-1),
and I2_D(b)=exp(gamma_E e)mu^(2e)Gamma(e)b^(-e).
Thus partial_b Tad_D=-I2_D and Integral_0^1 I2_D=-Tad_D(1).
The term -K I2_D/e in F_mix,D integrates to
-K Tad_D(1)/(2e), the exact required quadratic quartic-counterterm
insertion. Dropping it or integrating a four-dimensional
subtracted kernel would change the finite mass anchor.
