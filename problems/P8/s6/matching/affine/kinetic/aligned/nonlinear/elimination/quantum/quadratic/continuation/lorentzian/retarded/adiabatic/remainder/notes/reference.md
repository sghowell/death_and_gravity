# Actual varied eighth-order reference and the integrable readout tail

Let omega^2=m^2+q, z=q/omega^2 and lambda=omega'/omega=-H z.
Write t=omega^-2 and S=1+sum_{j=1}^4 P_j t^j, so W=omega S.
The P_j are exactly the frozen S6.55 recurrence coefficients.
The S6.59 estimates give 1/2<S<3/2, |W'/W|<4 and
|rho|<=C/omega^8, where

    rho=(1-S^2)/t-U-(log W)''/2+((log W)')^2/4.

All time derivatives act on the actual background, on
z'=-2H z(1-z), and on source jets. For r=delta log omega,
the variation of S is a polynomial dS with coefficients
delta P_j-2j r P_j. Linearizing the exact Riccati equation gives

    delta rho=-2S dS/t+2r(1-S^2)/t-delta U
              -(delta log W)''/2+(log W)'(delta log W)'/2,
    delta log W=r+dS/S.

At each order, the new dS coefficient appears with coefficient -2.
Solving that exact recurrence supplies orders six and eight as well
as independently reproducing the frozen second/fourth-order variation.

Let R=D_u S, T=-D_u R/2+lambda R/2 and
Q=(1-S^2)/t+2P_1. Then rho=N/S^2, with

    N=Q S^2+T S+3R^2/4,
    delta rho=(delta N S-2N dS)/S^3.

The coefficients below t^4 vanish exactly in both sectors.
Products of continuous coefficient majorants and S^-3<8 give

    C_delta,T=123763354409866,
    C_delta,L=152418105409587,

in |delta rho|<=C_delta ||n||_C10/omega^8. The coefficients use
at most the tenth source derivative. Its nonzero numerator fixture
at u=0,z=0 is -7/20736 in both sectors. This explains the source
norm used by this residual method; it does not prove that a weaker
norm cannot work by a different method.

The same coefficient bounds give |delta W|<=B omega ||n||_C10,
B<1/2, and |delta c|<7||n||_C10 for
c=d+W'/(2W). No time-dependent cutoff or changed initial state is used.

## Reference readout remainder

For the exact mass-readout weights A,B of S6.65,

    J_ref=omega/4 [A S+B/S+A t(c1+R/(2S))^2/S],
    c1=d+lambda/2.

Differentiate this expression, including A,B,t,S,R and c1.
Subtract the actual varied adiabatic orders 0,2,4.
A common denominator is 16 S^4. All numerator coefficients below
t^3 vanish exactly; the remaining positive majorants give

    |delta J_ref-delta J_ad,0..4|
       <=T_sector ||n||_C10/omega^5,
    T_T=48728198, T_L=85166609.

The physical momentum integral of omega^-5 is 1/(6pi^2 m^2).
Consequently this reference remainder has an explicit integrable
bound, uniform across the whole momentum range, including q=0.
