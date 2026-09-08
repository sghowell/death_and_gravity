# Finite observable transfer and unchanged regulator prescription

Use C_H=10,000,000 and the explicit initial comparison
|B0|<=C_H/nu^6 from the [construction](construction.md).
Both old and new modes are exact solutions of the same physical
equation; their relative Bogoliubov coefficients are constant
in time. Their phases need not be close. The local counterterm
prescription remains exactly S6.53's.

## Real physical observable difference

The old exact mode has coefficient one-norm below two relative
to the frozen reference throughout I. Thus

    |v_old|^2<=4/omega, |p_old|^2<12omega,
    |p_old*v_old|<7, p_old=v_old'-d*v_old.

The initial CCR gives |A0|^2=1+|B0|^2; since |B0|<1,
|A0|<2. The exact squared and mixed product identities of S6.54
then bound the state-change products by 6|B0| times the
corresponding old absolute reference products.

For either physical energy or pressure, the absolute readout
with |A|<=1, |B|<=3/2 is bounded by 9omega/a^3 for the
old mode. Therefore its state difference per polarization is
at most 54 C_H*omega/(a^3 nu^6). Summing all three and using
omega<=Amax*nu gives 162 C_H Amax/nu^5. With the exact integral

    integral d^3k/(2pi)^3 nu^-5=Amax^3/(6pi^2*m^2),

the total value difference is below 20 C_H/m^2. This applies
to either energy or pressure and includes the absolute integral.

For the differentiated readout use S6.54's exact quadratic
time-derivative identity. There is no reference-residual forcing
term in the difference of two exact solutions. The same bounds
give, per polarization,

    |delta Q'| <=816 |B0| omega+84 |B0| omega^2
                 <85 C_H omega^2/nu^6.

Here 816/m+84<85. Summing polarizations and using
integral nu^-4=Amax^3/(8pi*m), with the usual d^3k/(2pi)^3
measure, bounds the integrated derivative difference by
100 C_H/m. The integer bound follows from
255 Amax^5/24<100 and pi>3.

The state change also obeys the same clock-source identity, so

    |delta J_clock| <=100 C_H/m+240 C_H/m^2.

The second term is 3|H| times both value bounds, with 3|H|<=6.
All estimates are continuous on I. They are not an assumed
smallness of arbitrary state data or an uncomputed reference
vacuum energy.

## Formal dimensional extension of the finite state change

To retain the stated regulator prescription, continue the old
S6.53 Cauchy frequency W_4(D) and slope analytically near D=3,
and add the same real, D-independent cutoff correction delta W
and delta W' just constructed. The formal continued state need
not be all-order adiabatic away from D=3; it is an analytic
regulator family, not a Hilbert-space state in noninteger dimension.

The old frequency has modulus >=omega/2 on |D-3|<=1/4.
The added correction is bounded by F/omega^5, so the new one
has modulus >=omega/4 and a consistent analytic square root.
Put sigma=sqrt(1+delta W/W_4), on its branch near one.
Then |sigma|>=1/sqrt(2) and |sigma-1|<=|delta W/W_4|.
The diagonal amplitude shift has the exact form

    (1+sigma^2)/(2sigma)-1=(sigma-1)^2/(2sigma).

The off-diagonal Cauchy formula uses the two analytic mode signs,
without conjugating D. For each transfer column, these identities
give

    |A0-1|+|B0|
      <=[4F+8(G+4F)/m]/nu^6 <= C_H/nu^6.

The rational constant inequality is checked for both polarizations.
The old continued exact pair has absolute mode and physical
derivative envelopes 4/sqrt(omega) and 8sqrt(omega), from
S6.53's reference pair and finite-interval evolution estimate.
With |A|<=5/4, |B|<=3/2, its bilinear readout constant is
52omega. A transfer column differs from the identity by at most
b=C_H/nu^6<1; the product difference is bounded by
2b+b^2<=3b. Summing |D-1|+1<=13/4 yields

    |combined finite state-change readout|
       <=1000 C_H Amax/nu^5.

The dimensionally weighted radial tail is therefore bounded by
k^(-11/4) at infinity and k^(7/4) at zero, with the analytic,
compactly bounded prefactor already controlled in S6.53.
The finite state-change integral is holomorphic near D=3 and
tends to precisely the physical state-change integral. It adds
no new pole or finite counterterm. The first-derivative clock
source limit is justified distributionally against compact time
tests, as in S6.54, and the physical C1 bound identifies its
continuous representative.

## Physical bound

Let L=M*tau and R=m0*tau. Relative to M^2/tau^2, the added
energy or pressure error is 20 C_H/(R^2 L^2). Relative to
M^2/tau^3, the density-derivative error is 100 C_H/(R L^2),
and the clock-source error adds twelve times the value error.
Add these explicit errors to S6.53 and S6.54, without replacing
their already fixed local terms.

At L=10^12,R=1000 the state-change value bound is 2*10^-22,
the derivative bound is 10^-18, and the source bound is
1253/1250000000000000000000. Each of the five new-state
totals (energy, pressure, both density derivatives and clock
source) remains below 10^-14 in its stated reference units.
The certificate retains the exact fractions and decomposition.

This transfers the fixed-background, vector-only Gaussian results
to the explicitly constructed Hadamard state. It does not
establish other field or interacting-loop estimates, quantitative
higher time/functional variations, a quantum-corrected solution,
principal cones or cutoff, vacuum/Regge matching or V/G/B.
