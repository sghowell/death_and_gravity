# Actual mass-one scalar and the full on-shell proper forest

The frozen first fermion kernel satisfies
f_R,D(-q^2)=-(q^2+1)^2 P_D(q), with

    P_D(q)=Integral_T^infinity rho_D(v)/[(v-1)^2(q^2+v)] dv,
    T=4m^2.

Pairing the whole-fermion-cycle subgraph with its physical scalar
mass/residue counterterm therefore gives
V_OS,D=-(1/2)Integral_q(q^2+1)P_D(q).
Use (q^2+1)/(q^2+v)=1+(1-v)/(q^2+v). The polynomial trace is
scaleless in dimensional regularization. The remaining expression is

    V_OS,D=(1/2)Integral_T^infinity rho_D(v)Tad_D(v)/(v-1) dv.

This includes the nonzero mass-reference tadpole. The slope-reference
part closes into a constant trace, not an extra finite term.
The two proper fermion self-energy counterterms remain to be added.

The subtracted one-loop spectral identity is analytic in e on
0<Re(e)<3/2: at infinity the twice-subtracted kernel is integrable,
and at threshold the density exponent is 3/2-Re(e)>0.
It follows first from the convergent cut calculation and extends
on that strip by analytic continuation. After the polynomial contact
is removed, the vacuum momentum/spectral integrals converge jointly
for 1<Re(e)<3/2. This provides a genuine common integration domain
before continuing the resulting beta representation near e=0.

Set z=T/v and H(e)=exp(2gamma_E e)4^(-e)sqrt(pi)/
[2Gamma(3/2-e)]. In NY/Q^2 units,

    V_OS,D=T^2 Gamma(e-1)H(e)
       Integral_0^1 z^(2e-3)(1-z)^(3/2-e)/(1-z/T) dz.

Split 1/(1-k)=1+k+k^2+h3(k), h3=k^3/(1-k).
Using the exact F0(e) of S6.143, the three beta terms are

    -T^2 F0(e)(e+1/2)(e+3/2)/[(2e-2)(2e-1)],
    T F0(e)(3/2+e)/(1-2e),
    -F0(e).

The first is exactly the raw massless vacuum m^4 A^2 V_s/e^2.
The proof uses Gamma duplication and two explicit half-integer
recurrences; it is also checked by the direct Gaussian masters.
Adding the proper fermion MS counterterm from massless.md gives

    V_scalar,paired=(NY/Q^2){
       -19m^4+T(13/4+pi^2/8)-47/18-pi^2/12+delta_V}.

All three beta terms and the remainder are retained. This is not
the pure-MS scalar subgraph followed by an unrecorded finite mass
adjustment: the named physical inner forest is part of the formula.
