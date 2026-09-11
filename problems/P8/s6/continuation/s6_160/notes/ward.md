# Regulated Ward identity and third map order

At finite formal order, the change of integration variables is
the total-derivative identity

    integral D Psi delta/delta Psi [
      R exp(-S/h + sources/h)
    ] = 0,

with the action variation, Jacobian divergence and transformed
source variation kept together, and with a common regulator.
Iterating the identity through map degree three covers the
two-loop four-Phi support of the preceding note. Local
counterterms, including their poles and positive-epsilon
products, are transformed before finite parts. This is a
perturbative functional statement, not a global integration
cycle or nonperturbative inverse theorem.

For a local differential R, expand Tr log(1+R') at each finite
order. Each closed Jacobian ghost loop has identity propagators
and polynomial momentum numerators; its momentum integral is
scaleless in dimensional regularization, including with
external fields/derivatives. It therefore vanishes before
Laurent products. In finite-dimensional diagnostics the
Jacobian is not scaleless and must not be discarded.

For F=x+r x^3 and Gaussian covariance h/K, the transformed
Gaussian weight through r^3 is

    1-r K x^4/h
    +r^2[K^2 x^8/(2h^2)-K x^6/(2h)]
    +r^3[K^2 x^10/(2h^2)-K^3 x^12/(6h^3)].

Multiply by F^n and 1+3r x^2, then use
mu_(2k)=(2k-1)!!(h/K)^k. Factoring mu_n h^j/K^j and defining
R_k=product_(i=0..k-1)(n+2i+1), the three cancellations are

    (n+3)R_1-R_2=0,
    n(n+5)R_2/2-(n+7/2)R_3+R_4/2=0,
    n(n-1)(n+7)R_3/6-(n^2+6n+3)R_4/2
       +(n+4)R_5/2-R_6/6=0.

These polynomial identities hold for arbitrary even n; odd n
vanishes by parity. They therefore cover each monomial in
arbitrary parent polynomial weights, including counterterm
expansions and the opposite-Yukawa determinant pair. Sample
n=0,2,...,16 values and a coupled two-variable map are checked
independently; they are diagnostics, not the continuum proof.

Omitting the Jacobian gives Z's r coefficient -3h/K. Omitting
the physical-source correction gives the two-point r defect
-6h^2/K^2. Omitting the generated sextic gives Z's r^2 defect
15h^2/(2K^2). If the parent potential is q x^4/24, its generated
octic is q r^2 x^8/4; omitting it changes the connected physical
four-point r^2 coefficient by 1260 q h^5/K^6, a two-loop term.

For the opposite-flavor determinant pair, normalized weight
1-Y F^2/m^2, omission of the Yukawa R term instead uses
1-Y x^2/m^2. To first order in rY the difference in physical
coordinates is 2 rY Phi^4/m^2. After normalizing the moments,
the connected four-point defect is 48 rY h^4/(m^2 K^4),
already a one-loop term. These nonzero controls distinguish
the full prescription from the inadequate truncations.
