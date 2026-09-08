# The same reference, now varied along both physical metric sources

Use the S6.67 physical Hamiltonian vertices at D=3, retaining
a_m,NN and b_m,NN in their readout variations. The mode frequency
and moving-normalization variations are

    r_T=delta log omega_T=[1+beta(1-z)/2]n-z zeta,
    r_L=delta log omega_L=[1+(beta-alpha z)/2]n-z zeta,

    rg_T=delta log g_T=(zeta-n)/2,
    rg_L=delta log g_L=[(1+2z)zeta+(alpha z-1)n]/2.

Here z=q/(q+m^2), alpha=4/(9h), beta=28/(81h), h=(1+u^2)^3.
Differentiate all time-dependent coefficients and both source-jet
families. Then delta d=rg', delta lambda=r', and
delta U=(delta d)'+2d delta d.

## Physical outputs

For each baseline mode the physical readout has the form

    Q=(A |p_mode|^2+B omega^2 |v|^2)/(2a^3).

The energy weights are (1,1+beta(1-z)) for T and
(1-alpha z,1+beta) for L. The pressure weights are
(1/3,(2z-1)/3) for T and ((1+2z)/3,-1/3) for L.
Thus |A|<=1 and |B|<=2 on the complete interval/momentum box:
0<=z<=1, 0<alpha<=4/9, 0<beta<=28/81.

The source variations start from the independently checked
S6.67 coordinate-density current weights. For energy, negate
the lapse-current weights and subtract 3 zeta times the baseline
weights. For pressure, divide the scale-current weights by three
and subtract (n+3 zeta) times their baseline weights. These
additional terms vary a_hat^-3 and N^-1; no normalization contact
is silently discarded.

## Eighth-order comparison

Use W8=omega S with S=1+sum_(j=1..4)P_(2j) omega^(-2j).
The P_(2j) are the unchanged original-clock S6.55 coefficients.
If dS_j=delta P_(2j)-2j r P_(2j), differentiate the exact
Riccati residual:

    delta rho=-2S delta S/t+2r(1-S^2)/t-delta U
              -(delta log W)''/2+(log W)'(delta log W)'/2,
    t=omega^-2, delta log W=r+delta S/S.

Its triangular recurrence determines dS_j exactly. The first
two orders independently reproduce S6.67. The varied residual
starts at omega^-8. Its coefficient contains tenth derivatives
of both sources; at u=z=0 their respective fixtures are
-7/20736 and -1/512 in both sectors. The C10 hypothesis is
therefore explicit for this method, not replaced by a C4 claim.

The same exact denominator-clearing algebra as S6.66 gives
delta rho=(delta N S-2N delta S)/S^3. Continuous polynomial
majorants and S>1/2 yield residual constants

    C_delta,T=598084814046630,
    C_delta,L=627487153562466.

The relative varied-reference frequency bounds are below four;
the varied canonical readout rate bounds are below 64. Every
coefficient bound is reconstructed exactly over the continuous
u,z box.

## Integrable reference tail

Vary the full physical readout before subtracting its 0,2,4
coefficients. Its exact denominator is 16S^4, bounded below
by one. All three low numerator coefficients vanish separately
for T/L energy/pressure. The remaining tail is at most T/omega^5
times the joint C10 source norm, with constants

    T_energy,T=1631601633, T_energy,L=2349992276,
    T_pressure,T=591591679, T_pressure,L=737068541.

The physical radial integral of omega^-5 is 1/(6pi^2 m^2).
No ultraviolet cutoff or time/momentum sampling is used.
