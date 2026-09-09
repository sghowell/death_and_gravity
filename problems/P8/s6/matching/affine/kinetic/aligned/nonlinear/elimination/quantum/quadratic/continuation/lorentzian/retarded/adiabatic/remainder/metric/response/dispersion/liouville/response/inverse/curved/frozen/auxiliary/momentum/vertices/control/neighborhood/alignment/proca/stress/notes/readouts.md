# New readout, unchanged exact clock state

This is the S6.81 source-free, constant physical mass Proca action,
with m*tau=1000, evaluated on the original complete clock and the
unchanged S6.55 Borel-prepared state. Its exact vector mode operator,
Cauchy data and covariance coincide on that clock with the original
ones. Its functional metric derivatives do not coincide: energy and
lapse contact terms must be recalculated.

For both transverse and longitudinal canonical pairs set
p=v'-d v, with d=c1-lambda/2, lambda=-H*z,
z=(k_com/a)^2/omega^2, omega^2=m^2+(k_com/a)^2.
The new physical energy is

    Q=(|p|^2+omega^2|v|^2)/(2a^3).

There are two transverse polarizations and one longitudinal.
The actual canonical metric variation in S6.81 gives these unit
weights. Pressure, including its fourth-order subtraction and
covariantly matched local terms, is exactly the old pressure.
This is an operator/readout identity, not an assertion based just
on equality of covariances.

The exact mode equations are v'=p+d*v and p'=-omega^2*v-d*p.
For X=(|v|^2,Re(v*conj(p)),|p|^2), X'=M X with

    M=[[2d,2,0],[-omega^2,0,1],[0,-2omega^2,-2d]].

Put r0=(omega^2,0,1) and r_(j+1)=D r_j-3H*r_j+r_j*M.
Here D=partial_u+z'*partial_z+2lambda*omega^2*partial_omega^2
holds the comoving momentum fixed. Then Q^(j)=r_j X/(2a^3).
The code computes all twelve new rows for j=0,...,5. Rational
coefficient boxes on |u|<=1/2, 0<=z<=1 bound their products by
E_j*omega^(j+1), using (omega^-1,2,3omega) for the three
reference product envelopes. Every frequency degree and every
coefficient reconstruction is checked.

The new r1 obeys the ordinary conservation identity
r1+3H*(r0+r_pressure)=0 for each polarization. The new local
energy coefficients, obtained in S6.81 from the constant-mass
limit of the fixed D-dimensional covariant prescription at mu=m,
are

    -5/2,
    -10H^2,
    2(6H^2 H'+2H H''-H'^2).

Their pressure coefficients are unchanged and each adiabatic
order obeys rho'+3H*(rho+p)=0. The new finite ordinary stress is
therefore conserved, including all matched local terms. Do not
reuse the old nonminimal energy coefficients or clock source.
