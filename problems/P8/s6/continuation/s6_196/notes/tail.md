# Actual full stress-vector tail

## Inputs and conventions

Use the unchanged S186 pure all-order state and its complete
ten-feature Proca stress, including the temporal constraint,
all three physical polarizations and all nine creation pairs.
Put A=25/16, m=1000, nu=sqrt(m^2+|k|^2/A^2),
mu=sqrt(m^2+|l|^2/A^2), P=k+l. Fourier measures are
dk/(2pi)^3 and dP/(2pi)^3. For real smooth compact f write

    U_j(P)^2 = sum(r=0..j) integral dt ||partial_t^r fhat(t,P)||F^2.

S186's integrated pair bounds are

    |I_ref| <= R sqrt(nu mu)/(nu+mu)^3 U_3(P), R=1e25;
    |I_err| <= E sqrt(nu mu)(nu^-6+mu^-6) U_0(P), E=1e15.

The pair-exchange factor2, nine polarizations and splitting
|ref+err|^2<=2|ref|^2+2|err|^2 give the common coefficient36.
The actual stress variance is <1e50 N[f]^2. A local c-number
subtraction or finite mean does not affect the centered vector.

## Reference union tail

The discarded region is max(|k|,|l|)>K, not a cutoff on just
one leg. By AM-GM,

    nu mu/(nu+mu)^6 <= 1/(4nu^4), and also <=1/(4mu^4).

Use the union bound and exchange k,l. At fixed P the integral is
<=J4_tail(K)/2, uniformly for all P, where

    Jn_tail(R) = integral(|k|>R) dk/(2pi)^3 nu^-n;
    J4_tail(K) <= A^4/(2pi^2 K).

Consequently the complete reference variance tail is bounded by
36 R^2 A^4/(4pi^2 K) integral U_3^2
< R^2 A^4/K integral U_3^2, using pi>3.
The last numerator is <1e51.

## Remainder and the essential external split

Set L(P)=1+|P|/(Am). The triangle inequality gives mu<=L nu
and nu<=L mu. Therefore

    nu mu (nu^-6+mu^-6)^2
    <= 2L (nu^-10+mu^-10).

If |P|<=K/2, the union tail forces BOTH momenta>K/2.
Its integral is <=4L J10_tail(K/2). Direct radial domination gives

    J10_tail(K/2) <=64 A^10/(7pi^2 K^7)
                  <(64 A^10/63)/K^7 <100/K^7.

Also L<2(1+|P|^2); completing the square verifies this for every P.
Thus the low-external remainder contribution is at most
c_low/K^7 times integral(1+|P|^2)U_0^2, with
c_low=36 E^2*8*(64 A^10/63)<1e35.

For |P|>K/2 it is invalid to force both internal legs large.
Instead retain the entire internal integral4L J10, where
J10=5A^3/(512pi m^7)<m^-7, and use

    L/(1+|P|^2) <=4/K^2+2/(AmK)
                 <=[4/1000+2/(Am)]/K <1/(100K).

The high-external contribution is bounded by
c_high/K times integral(1+|P|^2)U_0^2, where
c_high=36 E^2*4/(100m^7)<1e10.
No external momentum cutoff has been imposed.

## Complete quantitative limit

Plancherel bounds both integral U_3^2 and
integral(1+|P|^2)U_0^2 by the stated N[f]^2.
For K>=1000 replace K^-7<=1000^-6 K^-1. The total numerator

    R^2 A^4 + c_low/1000^6 + c_high <1e52

is checked exactly. This proves tau_f(K)^2<1e52 N[f]^2/K
for nonzero f, with exact zero for f=0. The dominant displayed
tail is K^-1; it is not an assumed high-frequency asymptotic.
