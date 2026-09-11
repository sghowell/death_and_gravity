# Absolute overlapping-time memory remainder

Use the S186 pointwise readout decomposition, not its already
time-smeared estimate:

    u_actual=alpha_initial u_W8+e,
    ||e||<=4000 R nu^(-11/2), R<2e6.

The coefficient alpha_initial is constant in time and belongs
to the original all-order state. Each full stress creation-pair
amplitude obeys

    ||a_ref||F<=Cr sqrt(nu mu), Cr=1e9;
    ||a_actual-a_ref||F<=Ce sqrt(nu mu)(nu^-6+mu^-6), Ce=1e15,

where nu=sqrt(m^2+k^2/A^2), mu=sqrt(m^2+l^2/A^2), A=25/16.
The spatial volume already cancels both physical readout
a^-3/2 factors. Both actual transverse modes, the longitudinal
mode and the constrained temporal readout are retained.

At different times t,s, expanding the complete pair-product
difference gives BOTH mixed terms and the error-error term:

    conjugate(a_D) a_Gamma-conjugate(r_D) r_Gamma
    =conjugate(r_D)e_Gamma+conjugate(e_D)r_Gamma
       +conjugate(e_D)e_Gamma.

Their norm is bounded by

    nu mu [2Cr Ce h+Ce^2 h^2] ||Dhat(t,P)||F||Gammahat(s,P)||F,
    h=nu^-6+mu^-6, P=k+l.

No derivative is applied to rapidly oscillating actual mixing.

## Complete internal integrals

At fixed P, the Euclidean triangle inequality gives
mu<=L nu and nu<=L mu with L=1+|P|/(Am). Hence

    integral nu mu h <=2L J4;
    integral nu mu h^2 <=4L J10,

using h^2<=2(nu^-12+mu^-12), and

    J4=A^3/(8pi m), J10=5A^3/(512pi m^7).

All momentum measures are d^3k/(2pi)^3. The stress creation
pair exchange factor2 and the i/4 commutator give factor1
per ordered polarization pair, hence9 for all pairs.
The complete internal bound is

    36L(Cr Ce J4+Ce^2 J10).

## Retarded time and external momentum

Keep theta(t-s). Its modulus is at most1, so the time triangle
is bounded by the product of the two time L1 norms, and then
by the L2 norms because the slab has length1. This use of an
absolute bound does NOT replace the retarded triangle by a
signed full commutator square.

Since L<2(1+|P|^2), weighted Fourier Cauchy-Schwarz and
Plancherel give

    |R_memory|<=72(Cr Ce J4+Ce^2 J10) M[D]M[Gamma]
               <1e23 M[D]M[Gamma].

The inequalities J4<A^3/(24m) and
J10<5A^3/(1536m^7) use only pi>3. The exact rational
coefficient is checked. The integrable majorant is valid
on overlapping supports and at the time diagonal for this
DIFFERENCE; no reference singular extension has been chosen.
