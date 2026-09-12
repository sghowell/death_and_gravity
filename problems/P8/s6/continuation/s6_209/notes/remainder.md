# Quantitative all-P conversion remainder

Write e=p/K, q=j+d, M=4-q, and let N(Fbar) be the sum of the sup norm on [-1,1] and the first two u-derivative sup norms on |u|<=1/4. First assume K>=2m and p<=K/4. Then K-p>m, so the exact high-band shell formula applies without a lower-band boundary.

On the complex e-disc |e|<=1/2 and real -1<=u<=0, the chosen square root is analytic, |R|<2, and |A_M|<=17 for M1,...,4. For real e<=1/4, Cauchy's remainder after h=M is bounded by

34*2^(M+1) e^(M+1) <=1088 e^(M+1).

Multiplying by Fbar and integrating the hemisphere supplies the first part of the error. The angular coefficient is independent of the cutoff.

For the positive strip put u=e v, 0<=v<=1/2. Write delta=R-1. Then -e^2<=delta<=0 and

|delta-e^2(v-1/2)| <= e^4/3.

Taylor's square-root remainder gives this bound using e<=1/4. The finite binomial expansion through M<=4 consequently gives

|A_M(e v,e)-e^2(1/2-v)| <=4 e^4.

After both strip coefficients are retained, its error is at most e^5[2||Fbar||+||Fbar_second||/384]. For M3 the omitted linear Taylor term contributes e^4||Fbar_prime||/48. For M1 and M2 the unretained strip is bounded directly. In every case the remaining strip is below3N(Fbar)e^(M+1).

Thus 2pi(1088+3)<8192 bounds the complete hemisphere plus strip error. For q4 the exact logarithm obeys log(1/R)<=-log(1-e)<=4e/3 on the lost shell, giving the same conservative constant. For every original slot,

|Delta_q minus all retained shapes|
<=8192 p^(5-q) N(Fbar_jd)/K

before the harmless Fourier factor. No logarithmic UV slot is discarded.

S6.208 gives the full source-jet row bound

S <= B_j epsilon0^(-d) (m+p)^d ||D|| times the full source-jet norm,

with epsilon0=1/100 and the unchanged B_j. The degree-eight angular argument gives N<=6121S. Since p^(5-q)(m+p)^d <=(m+p)^(5-j), the unchanged spatial weight (m+p)^s<=8m^s(1+p^2)^3 controls every small-transfer remainder.

For LARGE TRANSFER p>K/4, do not use the small-e expansion outside its domain. The absolute original lost UV integral is at most its complete one-ball high-band integral. For q<4 this is bounded, after the Fourier angular factor is bounded by one, by

B_j epsilon0^(-d)(m+p)^d K^(4-q)/(4-q).

Use K<4p to replace it by

B_j epsilon0^(-d) 4^(5-q)(m+p)^(5-j)/[(4-q)K].

For q4 retain log(K/m)<=K/m. Multiplying by K and using K^2<16p^2 bounds this term by

16 B_j epsilon0^(-d)(m+p)^(d+2)/(m K).

The maximal degree is six, so the existing X46 norm still suffices.

Each retained shape is bounded by64N(Fbar). The hemisphere polynomial has sup at most17/4 and the strip coefficients are1/8 and1/48. Again using K<4p, its contribution is at most

64*6121 B_j epsilon0^(-d) 4^(5-j-d-h)(m+p)^(5-j)/K.

Sum every original slot and shape, not just those that later cancel. The exact weighted constants in centered.tail_constants are approximately

small transfer: 1.495065730988703744e39,
large-transfer original shell: 4.77080143872e32,
large-transfer retained shapes: 3.3840419963904e35.

Their exact total is 1.495404612268486656e39, strictly below1e40. The full coefficient includes all source jets; time and momentum Cauchy-Schwarz and Plancherel give

|E_K|<1e40 ||D||L2 X46[Gamma]/K

for all P and K>=2m. At P0 both original masks coincide exactly, so the conversion is zero. The complete shape sum has only cubic and linear terms by the actual symmetry proof; the bound was obtained without deleting the individual finite terms first.

This is a quantitative regulator comparison. Its finite cancellation does not identify the remaining one-ball UV/contact current with the fixed covariant local Hessian, and it does not supply a reduced inverse.
