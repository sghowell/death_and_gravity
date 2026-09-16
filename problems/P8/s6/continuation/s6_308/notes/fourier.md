# Independent Fourier reconstruction and order of limits

After the same IR division and the fixed-order dimensional limit,
the D4 impact factor is

exp[i*eta*(2/V-gamma_E)]*(pi*b^2)^(-i*eta).

For a Gaussian Abel parameterh>0, its radial transform is determined by

2pi int_0^infinity b db J0(sqrt(tau)*b)*b^(-2i*eta)*exp(-h*b^2)
=pi*h^(i*eta-1)*Gamma(1-i*eta)
 *1F1(1-i*eta;1;-tau/(4h)).

ExpandJ0 into its entire power series and integrate each Gaussian
moment, or use its angular representation followed by the Gaussian
Mellin transform. The resulting confluent hypergeometric function is
also checked against literal high-precision radial quadrature.

On removingh as a tempered-distribution boundary value, the transform
at every fixednonzerotau is

4pi/tau *Gamma(1-i*eta)/Gamma(i*eta)*(tau/4)^(i*eta).

The Fourier transform of the subtracted constant is supported at
q=0 and is kept separate. Multiply the displayed expression by
2D/i and the impact prefactor. Gamma(1+i*eta)=i*eta*Gamma(i*eta)
then reconstructsV/(kappa*tau) times the exactC_eta in the formulation.

This provides a second derivation with the same finite constants.
In particular the rawln(4pi),gamma_E and2/V factors are not replaced
by an arbitrary inverse-length regulator.

The order of operations matters. At finitee, the infinite absolute
rung sum is not uniformly integrable nearb=0; sufficiently high powers
can meet ultraviolet Mellin poles. The proof does not interchange
that sum and the unregulated integral. It first takes the specified
IR-divided limit at each fixedorder. In D4 with a Gaussian factor,
the resulting logarithmic Taylor series is dominated nearb=0 by
b^(-2|eta|), which is locally integrable for|eta|<1. The Gaussian
controls infinity. This justifies the resummed Fourier definition
for the finite leading class, after which the Abel regulator is
removed distributionally awayfromq0.

A change in unspecified short-distance dynamics or local operators
belongs to the omitted full-theory remainder. The Fourier construction
does not furnish its numerical bound, define the full interacting
quantum measure, or make the exact forward amplitude finite.
