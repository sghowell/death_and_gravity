# Full new contractions

## Local six-field terms

Use Fourier derivatives i p_mu, all momenta incoming, external p_i²=1, sum p_i=0, and s+t+u=4. Internal tadpole legs carry k,-k. Transform the literal operators, sum all720 assignments, and multiply by the identical internal factor1/2.

The coincident mass1 contraction in the inherited convention is I=-(Delta+1)/(16pi²). In dimensional regularization, radial powers reduce to I times mass1 powers, since the remaining polynomial integrals are scaleless. Retain

<k_mu k_nu>=eta_mu_nu/d,

<k_mu k_nu k_rho k_sigma>=(eta_mu_nu eta_rho_sigma+eta_mu_rho eta_nu_sigma+eta_mu_sigma eta_nu_rho)/(d(d+2)).

Let sigma2=s²+t²+u² and sigma3=stu. For a raw factor f(d), the finite MS contribution is -[f(4)-2f'(4)]/(16pi²). This follows by multiplying f(4-2epsilon) by -(1/epsilon+1) and subtracting its pole. Premature d4 substitution loses finite terms.

The six finite factors, before multiplying their Lagrangian coefficients and -1/(16pi²), are:

| Operator | Finite factor |
| --- | --- |
| Phi^6 | 360 |
| Phi^4Y | 72 |
| Phi²Y² | 2sigma2+20 |
| Y³ | 15sigma2-60 |
| Phi²(L3-L4) | 14-sigma2/4-3sigma3/2 |
| Y(L3-L4) | -21sigma3/8 |

Independent Wick counting gives raw Phi²Y²=Y4+8(2+4/d), raw Y³=(3+12/d)Y4, where Y4=2sum(channel-2)². The raw Galileon factors are -3sigma3/2+(2-6/d)sigma2-16+80/d and -3(1+2/d)sigma3/2 respectively. A constant quartic coefficient a multiplying L3-L4 has vertex -3a stu/2: its constant term is zero. The old -8gamma arose from a separate Phi^4 potential. Constant-coefficient integration-by-parts identities cannot be applied under a field-dependent coefficient.

Assembling all coefficients gives A_extra,local=-(K0+K2 sigma2+K3 sigma3)/(16pi²), with

K0=-2(300Ck+6000k lambda+1602721)/(25k²),

K2=120lambda/k-2362496/(25k²),  K3=47232/k².

## All three source classes

The heavy source is J2+J4+..., with J2=gPhi²/2 and J4=cPhi²Y, c=-4g/k. Its Gaussian cross term J2 G_H J4 produces three distinct light-contraction classes.

1. **Across the factors.** There are four choices of the single external light attached to J2. The full H Phi²Y vertex is -4c sum_{i<j}p_i.p_j=-2c(P_H²-sum_i p_i²). With three external mass1 legs it becomes 8g(k_H²-q_light²-3)/k. Denominator reduction therefore gives the entire mixed contribution

A_extra,mix=32g²[ A0(1)-A0(n)+(n-4)B0(1;1,n) ]/(16pi² k).

After subtraction, the bracket is Kmix=1-n+n log n+(n-4)B_MS, B_MS=-integral_0^1 log F(x)dx, F=(1-x)²+nx.

An independent Feynman shift q -> q+xp changes the numerator to -2(x+1), since p²=1. Consequently the UV residue is -3 and the finite bracket is exactly integral_0^1 2(x+1)log F dx. Integration by parts using F'=2x+n-2 gives the A0/B0 expression above. This also avoids cancellation-sensitive evaluation of its individually large terms.

2. **Within J4.** Both the Phi² pair and the Y pair contribute. The corrected source is c I H(Y+Phi²), whose full on-shell two-light vertex is c I(4-s). Including both exchange ends gives

A_extra,vertex=8g² sum_channels (4-channel)/(n-channel)/(16pi² k).

This is a full rational heavy exchange correction, not a contact term.

3. **Within J2.** The heavy tadpole gives g I J4/(2n). The original zero-onepoint condition fixes the full source counterterm -j1H with j1=gI/2. Its cross term is -g I J4/(2n), cancelling exactly. This is a new quartic derivative cancellation in addition to the old quadratic one; it cannot be silently discarded.

The old polynomial diagrams plus these classes exhaust the one-loop four-Phi coefficient. The mixed bubble is evaluated at a single external p²=1, below its heavy-light threshold. The local terms have no channel cuts, and the new heavy pole lies outside the physical window. Thus all new pieces are real there and the old complete tree-squared first elastic discontinuity is unchanged. This says nothing about a resummed heavy resonance.
