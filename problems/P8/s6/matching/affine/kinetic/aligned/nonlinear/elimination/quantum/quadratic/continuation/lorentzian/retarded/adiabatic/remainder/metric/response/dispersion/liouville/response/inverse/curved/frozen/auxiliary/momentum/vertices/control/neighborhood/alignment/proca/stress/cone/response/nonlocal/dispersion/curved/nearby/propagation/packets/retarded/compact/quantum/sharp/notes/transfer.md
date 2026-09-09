# Uniform real-time transfer and sharpened front error

The unchanged exact packet system is
Y'=(k J+L0+L1/k+L2/k^2+L3/k^3)Y.
Every entry of the four literal L matrices is bounded on the
actual complex disc; complete Laurent row sums are below 10^4.
The moving Hamiltonian, volume and q'= -2Hq contacts were included
before deriving these matrices. No finite-q Legendre inverse is used.

Set W=S^-1 Y. Then W'=(i k Lambda+B0+Q/k)W, where
B0=S^-1 L0 S-S^-1 S' and
Q=S^-1(L1+L2/k+L3/k^2)S.
For all real k>=1 the direct complex bounds yield
norm(B0)<4*10^4*6+4/50<250000,
norm(Q)<3*4*10^4*6=720000<10^6.
In fact norm(B0+Q/k)<960001<10^6.
On real times Lambda is real diagonal. Its phase propagator
has infinity norm one. In that interaction picture, Gronwall
therefore gives exp(T*10^6)=exp(1/10)<=10/9<2.
The inequality follows from the positive exponential series,
termwise bounded by the geometric series for 0<=x<1.
Returning through both basis endpoints gives norm(T_Y)<48.
Since the physical scalar endpoint prefactor is below two,
abs(Ghat)<96/k<100/k for every k>=1.

For 0<=k<=1 retain the original chart, whose generator is
polynomial in q=k^2/R^2. Here q<=4. Summing all coefficient
majorants gives norm(M_X)<=28169/64, including weighted damping.
Thus T norm(M_X)<1/2 and norm(T_X)<2. The original source is
N e^3 e_4 J with N e^3<2; its chi response satisfies abs(Ghat)<4.
The origin is a regular Fourier-multiplier extension, with the
same homogeneous-mode scope boundary as S6.91.

Now improve the high-frequency approximation. The complex B0
bound on |u|<=T and Cauchy radius T/2 around the real inner
interval give norm(B0')<=5*10^12. The actual signed frequency
derivatives are below 1/100, and every gap has modulus >10^-7.
The parent identity diag(B0)=0 and exact matrix
Z_ij=i(B0)_ij/(lambda_i-lambda_j), Z_ii=0 are unchanged.
Row-sum estimates give norm(Z)<3*10^12 and
norm(Z')<=5*10^12/10^-7+
           2*(1/100)*250000/(10^-7)^2<6*10^19.
The second term retains the derivative of the gap.

With C_k=I+Z/k, the exact transformed equation is
V'=i k Lambda V+(1/k) C_k^-1(B0 Z+Q+QZ/k-Z')V.
For k>=10^14, delta=3*10^12/k<=.03 and C_k^-1 has norm <2.
The complete numerator bound is below
2*(250000*3*10^12+(3/2)*10^6+6*10^19)<2*10^20.
Hence eta=T*2*10^20/k<=1/5.
The exact V transfer differs from the phase transfer D by
at most 2 eta, and has norm below two.

For both arbitrary ordered endpoints, use
C_t V C_s^-1-D
 =(C_t-I)V C_s^-1+(V-D)C_s^-1+D(C_s^-1-I).
Its norm is at most 6 delta+4 eta. Both S endpoints give
norm(T_Y-T_Y0)<=24*(6*3*10^12+4*T*2*10^20)/k.
The actual scalar prefactor <2 then yields, uniformly,
abs(Ghat-Ghat0)<10^16/k^2 for k>=10^14.

This is the same two-front leading kernel as S6.91, with a new
error proof valid below the previous threshold. It remains an
absolute error bound at phase zeros. No old frequency-domain
guard is bypassed: this child proves and validates its own domain.
