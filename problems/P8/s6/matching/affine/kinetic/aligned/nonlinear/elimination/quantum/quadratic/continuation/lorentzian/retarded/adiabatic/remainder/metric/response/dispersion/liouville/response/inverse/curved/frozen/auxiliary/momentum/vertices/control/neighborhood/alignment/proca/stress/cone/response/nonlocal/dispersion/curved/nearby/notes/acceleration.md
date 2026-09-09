# Actual central acceleration, including the boundary primitive

Use the same dimensionless clock units as S6.83. R is the hat scale;
p and ell are the trace and matter momenta per hat volume. At zero
spatial and vector data the actual Hamiltonian per hat volume is

    Hcal=N[(p-b)^2/(4a)-f0+ell^2/(2U)].

Every lapse derivative below is at fixed canonical phase data.
The already-derived central constraint fixes ell0^2=P2(N0) only
after taking such derivatives. Differentiating along that family
instead would change the problem.

At u=0, p=0, b=0, with even a,U,f0 and odd b. Let a0,U0 be their
central values, b1 the actual first time derivative, f2 the actual
second time derivative, p1 the trace acceleration and h1 the
first derivative of H_hat. Then

    p1=N(f0+ell^2/(2U)),
    h1=N(p1-b1)/(6a),
    ell''=-3ell h1.

The exact primitive has I(u,1)=0 and the fixed integrand I_s.
Independent differentiation of that literal integrand gives

    I_s,u(0,s)=-9(s-1)(s+1)/s^(3/2),
    I_s,uuu(0,s)=189(s-1)(s+1)(s^2+3)/(2s^(7/2)).

Integrating from 1 to 1/N, with the basepoint kept fixed, gives

    I_u=-6(-4N^(3/2)+3N^2+1)/N^(3/2),
    I_uuu=63(16N^(3/2)+9N^4-30N^2+5)/(5N^(3/2)).

The first agrees with S6.83. Both antiderivatives and basepoint
conditions are checked. Since f0=U*F_transformed-I_u/N,
its actual second time derivative contains -I_uuu/N. It is not
obtained by holding the placeholder I_phi constant. The resulting
literal coefficient is

    f2=3[-134400000N^(3/2)+3374997N^6-59209999N^4
         +185190007N^2-26875005]/(2000000N^(5/2)).

Write F=Hcal_N and h0=Hcal_NN at the central fixed phase.
Parity gives N'=ell'=0. Differentiating F=0 twice along the
actual solution gives

    N''=-(F_uu+2F_up p1+F_pp p1^2+F_ell ell'')/h0.

Here F_up=partial_N[-N b1/(2a)] and F_pp=partial_N[N/(2a)].
The F_uu derivative is taken at fixed ell before ell^2=P2(N).
A separate degree-two path expansion of the literal fixed-phase
Hamiltonian reproduces the numerator. Omitting I_uuu changes N''
by I_uuu,N/h0, a retained nonzero off-clock negative control.

The exact answer is

    N''=3N(N^2-1)
      [13921850250011N^6-93228677254991N^4
       +129756291744989N^2-31453170750009]
      /[250000(16874985N^4-14244998N^2-8625003)].

This is not yet the physical Hubble acceleration. For
omega=log(eomega)=-log[(h-1+N^-2)/h]/4 at the center,

    omega_uu=3(1-N^2)/2, omega_N=1/(2N).

Since physical dt=N du and a=eomega R,

    dH_physical/dt=(h1+3(1-N^2)/2+N''/(2N))/N^2.

Its exact rational numerator is

    11390604750009N^8-127400509505002N^6
    +319016194499976N^4-255642003494998N^2+46640697750015,

with denominator
250000N^2(16874985N^4-14244998N^2-8625003).
At N=1 these formulas give N''=0 and physical Hubble
acceleration 4. At N=1+10^-6 they give approximately
-0.0000760488361 and 3.9999569757 respectively; these decimals
are illustrative, not the certificate.

Continuous rational enclosures on the larger symmetric interval
|N-1|<=10^-6 prove P2>1/200, P2<(11/100)^2, h0<-2,
|N''|<1/1000 and physical Hubble acceleration >39999/10000.
The actual family [1,1+10^-6] is a subset. Exact polynomial
coefficient bounds establish every inequality without a sampling
grid or a floating-point zero criterion.
