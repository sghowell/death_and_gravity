# One fixed local full-model subtraction scheme

Keep the same renormalized tree parameters M, g and lambda4 throughout
this calculation. The parent's historical key bare_polynomial_quartic
labels its displayed polynomial coefficient; the fixed value now specifies
the renormalized tree parameter in the following explicit scheme. It is
not identified with an unregulated bare coupling.

Let I0 denote the dimensionless regulated zero-external-momentum,
mass-one light bubble, normalized so its four-dimensional radial expression
is integral_0^infinity y/(y+1)^2 dy. Use a common translation-invariant
regulator for all diagrams and counterterms, and take its removal only
after the finite differences. I0 is not assigned an unregulated finite value.

At large loop momentum the two nonlocal internal heavy exchanges vanish;
the vertex tends to C(z)=-lambda4+g/(M-z). Subtract C(z)^2 I0/(32pi^2)
in each of the three channels. The same subtraction follows from local
counterterms in the **full two-field model**, not an arbitrary b2 contact.
Writing F=I0/(32pi^2), choose the total counterterms

delta lambda4=3lambda4^2 F,
delta g=2lambda4 g F,
delta M=g F.

Literal variation of A0=-lambda4+g sum h(z), h(z)=1/(M-z), gives

delta A0=-delta lambda4+delta g sum h(z)-g delta M sum h(z)^2
        =-F sum_z C(z)^2.

Every coefficient and sign is checked symbolically before integration.
These are the allowed local Phi^4, H Phi^2 and H^2 terms. In particular
delta g=2G delta G at this order. The heavy parameters are not defined
by an exact stable-heavy on-shell prescription; the heavy field can decay.

The corresponding constant-field effective quartic counterterm is

delta lambda4 -3delta g/M +3g delta M/M^2
 =3F(lambda4-g/M)^2.

With A=(lambda4-g/M)/2, its contribution to the potential is
A^2 I0 Phi^4/(64pi^2), precisely cancelling the ultraviolet coefficient
of the complete Schur-Hessian logarithm in S6.110. Add the single
momentum-independent **finite** Phi^4 counterterm required to retain
that parent's fixed zero-field quartic subtraction. Its forward second
derivative vanishes, so it cannot tune the b2 bound proved here.

Mass, kinetic and one-point terms are fixed by the S6.112 light on-shell
conditions and the chosen heavy one-point convention. After integrating
H, the one-point term changes only the local light quadratic/vacuum terms.
The explicit potential curvature conversion 1-p from S6.112 is retained;
no extra curvature-one condition is imposed.

The displayed interaction counterterms are TOTAL coefficients in the
canonical renormalized fields. Independent bare coupling relations absorb
the field-normalization counterterms. Because the light propagator has
unit renormalized LSZ residue, adding another wave-function multiple of
the tree amplitude would double count that normalization. Counterterms
inserted inside one-loop integrals are two-loop terms and are not included
or bounded here.

The remaining finite radial integrand is exactly the one in notes/radial.md.
The angular part is already ultraviolet integrable. Both are calculated
in this one scheme, with no later finite b2 subtraction.
