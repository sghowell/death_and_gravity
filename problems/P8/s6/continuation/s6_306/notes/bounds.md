# Explicit whole-amplitude and full-Born interference bound

Let L=ln n, n>=128. Every physical channel is in[-12,0] or[25/4,16].
Then q_a>=n/2, |a|<=16, |a-2|<=14, |a+2|<=18 and|V|<=194.
The master bounds are established in notes/contour.md. The auxiliary
1/4<=E^2<=1 gives|c_E|<3 and|c0|<3.

## Remaining heavy masters

For U use its positive denominator
(n-a)v+r[a v+(1-v)^2]>=nv/2+r(1-v)^2.
On v<=1/2 integrate1/(nv/2+r/4); its contribution is at most
2[ln(n+1)+1]/n. On v>=1/2 the remaining integral is<=2ln2/n.
Thus |U|<2(L+3)/n. This bound includes the integrable r=v=0 corner.
The original mixed bubble denominator1+(n-2)v+v^2 is between1 andn,
so0<=Lmix<=L. The heavy on-shell logarithm obeys|Lh|<L+2.

For the complete K1 expression in S294, its first, constant and squared
terms give900/n+2048/n^2. Its remaining absolute bound is
8L+16(1+32/n+128/n^2)(L+2)+128(L+1)/n.
At n>=128 this is strictly below40(L+2), without assuming cancellation
of the large logarithms. The exact box collapse yields
|ln(1-a/n)|<=|a|/(n-|a|)<=32/n and consequently
|W|<[4(L+9)+512]/n=4(L+137)/n.

## Complete selected nonendpoint and endpoint sums

The first Htree group in F_g is below
6[3(194*128+2*16)+3(3*194*16+2)]/n=615204/n.
The three K1 terms contribute120(L+2). The T/U terms are below
[864L+6624]/n. All six ordered boxes give18624(L+137)/n^2.
The resulting positive polynomial, evaluated at n128 in its decreasing
inverse-n factors, is strictly below6000(L+1). Both coefficients of the
remaining linear polynomial in L are positive exact rationals; the
certificate checks these margins, rather than sampling n or angles.

For the whole endpoint |N_a|<=2+32+256=290, so
|2N_a+a+a^2/2|<=724. The full bounds |f1/a|<c,|f2_g|<2c imply
|A_endpoint|<3(724+36)c/kappa<2500g^2/(16pi^2 kappa).

For the full quartic representative,
|T0|<4657, |T1|<37969, |E1|<433, E0=5/3.
Therefore |F_C|<2*37969+433+2*3*4657+5<105000.
With the stated |C|<=6g^2/n its contribution is below630000/n in the
same g^2/(16pi^2 kappa) units. Adding every component gives

|A_known|<14000(L+1)g^2/(16pi^2 kappa).

All coefficient, contour, root and endpoint margins are checked exactly.
The general bound API explicitly rejects violation of the quartic
premise; it does not assume it for arbitrary unrelated theories.

## Original parameters and all angles

The unchanged n=10^200/512+2 is below2*10^197, hence L<600.
Since pi>3, 14000*601/144<10^5. Thus
|A_known|<10^5g^2/kappa. S297 proves on the entire compact window
A_m>4g^2/n^3 and A_G>0, with A_B=A_m+A_G. It follows that

|2Re(A_known)/A_B|<50000n^3/kappa<4*10^-204.

This uses the full Born denominator, including the pole-dominated
inner angular region. It is not a fixed-angle small-gravity expansion.
At a0 the logarithmic box multiplier and light/heavy functions are
continuous, and the OS difference quotient removes the endpoint1/a.
The known mixed amplitude is therefore finite at either angular
endpoint. The positive Born pole diverges, so the normalized known
interference tends to0. No finite forward cross section, missing loop
square bound or interchange with the fixed-transfer Regge limit follows.
