# Complete finite-reference bulk, endpoint and tail

Put nu=sqrt(m^2+k^2/A^2), mu=sqrt(m^2+l^2/A^2),
A=25/16 and P=k+l. Each detector reference amplitude is
bounded by C sqrt(nu mu), C=1e9. From notes/majorants.md,
the sixth bulk has nine-pair numerator

    Q6=9 C^2 2^6 *97423/r^6.

Keep the retarded time triangle. Its absolute value is
bounded by time L1 products, then by the unit-time L2
norms. The source requires jets through6, and the
detector only its L2 norm.

The j5 endpoint has one g outside L^5. Its numerator is

    Q5=9 C^2 2^6 *7916/r^5.

It has the SAME sixth inverse-frequency power. Its
equal-time product is bounded directly by time L2
Cauchy, requiring source jets only through5.
Both complete internal weights are

    nu mu/(nu+mu)^6 <=1/(4nu^4).

The exact full integral is bounded by J4/4, with
J4=A^3/(8pi m)<A^3/(24m). This is uniform in ALL
external P. Fourier Cauchy and Plancherel therefore
give, without a spatial derivative weight,

    |bulk6|<1e48 ||D||L2 S6[Gamma];
    |endpoint5|<1e43 ||D||L2 S6[Gamma];
    |F|<= (Q6+Q5)J4/4 ||D||L2 S6[Gamma]
         <1e48 ||D||L2 S6[Gamma].

The exact rational sum, not merely the two rounded
displays, proves the last coefficient.

For the removed union max(|k|,|l|)>K, symmetry gives
J4_tail(K)/2<=A^4/(4pi^2 K)<A^4/(36K).
Thus |F-F_K|<1e52 ||D||L2 S6[Gamma]/K uniformly in P.

At P0 the corresponding fifth-bulk weight is1/(32nu^3).
Its three-dimensional radial absolute majorant behaves
as1/(32k) and is not integrable. The sixth has
1/(64nu^4) and is integrable. This motivates the sixth
step; it is not a claim that every exact fifth-order
current diverges after all possible cancellations.
