# A canonical chart regular at every real momentum

Put mu=sqrt(1+k^2)>=1 and q_reg=mu^2/a^2, while the ORIGINAL
generator still uses q=k^2/a^2=(mu^2-1)/a^2. Exchange
b_reg=-pi_v/(2 a^3 q_reg), P_reg=2 a^3 q_reg v.
With x=(b_reg,chi), p=(P_reg,pi_chi), set
p_tilde=p-a^3 q_reg B_q x, B_q=diag(-2 Theta/Lambda,0), and
Y_reg=(mu a^(3/2) x,a^(-3/2) p_tilde).

The exact map E_reg satisfies E_reg Omega E_reg^T=mu Omega;
E_reg/sqrt(mu) is symplectic. Differentiating it BEFORE changing
q to (mu^2-1)/a^2 gives

    M_reg=(E_reg'+E_reg M_density)E_reg^-1
         =mu J+L0+L_(-1)/mu+L_(-2)/mu^2+L_(-3)/mu^3.

Every coefficient is retained. The leading J is exactly S6.99's
canonical J. Native identities also check the complete k=0 value,
inverse, physical probe force and original observable row. This
Fourier-dependent chart is a proof device, not a new local observable.

On |u|<=h=1/100, set G=-Omega J. Direct exact rational enclosures
of its entries and derivatives prove

    (1/14)I <= G <= 13 I, ||G'||_2 <= 4,
    sum_(j=-3)^0 ||L_j||_2 <= 2.

The rational bounds enclose numerator/denominator polynomials over
the WHOLE interval, rejecting every denominator enclosure through
zero. Matrix norms use both row and column bounds where needed;
G and G' are symmetric. The chart denominators are nonzero on
this interval, not on an assumed global gamma chart.

Because GJ+J^T G=0, the exact energy satisfies
|d log sqrt(Y^* GY)/du|<=392. Consequently

    ||U_reg(t,s,k)||_2 <= sqrt(182) exp(392 |t-s|)

for every real k>=0 on the strip. The original k-chart has the
same G and a smaller lower-order remainder for k>=1, hence the
same bound. The time derivative of the original observable has
row norm at most 16 sqrt(mu) in normalized regular variables;
the spatial-gradient sum has row norm at most sqrt(mu).
For k>=1, the corresponding original-chart constants are 16 and 1.

The whole-interval derivative jets used later employ real/complex
balls with outward enclosures, at a scoped 160-bit precision;
see the [python-flint real-ball documentation](https://python-flint.readthedocs.io/en/latest/arb.html).
No binary diagnostic decimal is used as a certified endpoint.
