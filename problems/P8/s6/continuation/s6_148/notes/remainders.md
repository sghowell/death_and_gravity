# Uniform finite remainders, including both endpoints

For 0<=k<=1/16,

    h2(k)=k^2(3-2k)/(1-k)^2 <=4k^2,

because 4k^2-h2=k^2(1-6k+4k^2)/(1-k)^2 and
1-6k+4k^2>=5/8. Thus the finite tadpole remainder obeys

    |delta_T| <= T*4/T^2 [3+2+1]=24/T.

Here |4log2-3|<=3 and the two unit-interval negative-log moments
both equal one; (1-z)^(3/2)<=1. For the unexpanded regulated
remainder, h2(z/T)=O(z^2) cancels z^(-2); on |e|<1/8,
the remaining endpoint powers and every finite e derivative have
integrable majorants. This justifies finite-part extraction.
The simple pole and its finite prefactor product are both kept.

For the bubble, h(z/T)<=3z/T. The first logarithmic group therefore
has absolute integral at most 18/T. Independently put k_M=Mz/T.
Since A(z/T)/(1-k_M)<=(16/15)^3<2 for 1<=M<=T/16,

    |Integral w A(z/T)L(Mz/T) dz|
       <=2(M/T)[log(T/M)+1].

The exact correction is consequently bounded by
18/T+2(M/T)[log(T/M)+1]. Its e-dependent compact subtraction
has integrable endpoint powers for |e|<1/8, as in S6.143,
with the two fixed ratios retained separately.

Use pi<4 and 0<log2<1. Then
13/4+pi^2/8<21/4, 47/9+pi^2/6<8, and
47/18+pi^2/12<4. If n is the least integer with T/M<=2^n,
log(T/M)<n. For rational inputs, the code uses

    |T_insert,MS| <= (12Y_upper/Q_lower^2)[21T/4+8+24/T],
    |F_MS(M)| <= (12Y_upper/Q_lower^2)
                   [4+18/T+2(M/T)(n+1)].

N=6; Q_lower<=144 is a lower bound on physical 16pi^2.
At the actual model n=676 and M/T is about 4.8828125e-204.
No endpoint, subleading term or finite mass-ratio correction
is rounded to zero.
