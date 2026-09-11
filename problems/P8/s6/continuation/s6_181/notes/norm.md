# Summing every band and identifying the same inverse

For j<=0, h=2^j<=1 and the band lies within 0<x<=2.
Since h/x<=2, the threshold estimate gives a valid common
band constant

    A_j<=4 C_small sqrt(2h)<6 C_small 2^(j/2).

The complete geometric sum is
sum_(j<=0)2^(j/2)=2+sqrt(2)<4. Thus

    sum_(j<=0) A_j<24 C_small.

For j>=1 the band has x>=2^(j-1)>=1. The tail estimate and
h/x<=2 give

    A_j<=4 C_tail/[log(1+2^(j-1))]^2.

For j=1 this is less than16 C_tail, since log(2)>1/2.
For j>=2, log(1+2^(j-1))>(j-1)/2, so
A_j<16 C_tail/(j-1)^2. The elementary integral comparison
sum_(n>=1)n^-2<1+integral_1^infinity t^-2 dt=2 gives

    sum_(j>=1) A_j<48 C_tail.

All compact sine kernels therefore sum absolutely in L1:

    sum_j ||K_j||L1
       <35(24*200000000+48*20000000)
       =201600000000<300000000000.

This establishes an L1 limit K_tilde, but identification is
still required: an arbitrary summed kernel would not suffice.

For real Laplace parameter p>0, each compact band admits
ordinary Fubini and the elementary sine transform gives

    Laplace[K_j](p)
      =-2 integral_0^infinity a(x)phi_j(x)(2+x)
                              /[p^2+(2+x)^2] dx.

The absolute L1 sum permits interchange with the time integral.
Because a phi_j is nonnegative, Tonelli permits summation of
the spectral integrals. The resulting integral is finite:
the threshold is O(sqrt(x)), while the tail is bounded by a
constant times 1/[x log(x)^2]. The full partition equals1.
Hence

    Laplace[K_tilde](p)
      =-2 integral_2^infinity rho(Omega^2)Omega
                              /[p^2+Omega^2] dOmega
      =-integral_4^infinity rho(tau)/(p^2+tau) dtau
      =R(p^2)=1/F_1(p^2),

the same source-pinned massive range inverse from S6.86.
That source already supplies an L1 causal kernel with this
Laplace transform and no instantaneous term. Uniqueness in
L1 identifies K_tilde with it. Explicitly, the transform of
the difference is analytic in Re(p)>0 and vanishes on the
positive real axis; the identity theorem and Fourier uniqueness
applied to exp(-ct) times the difference give zero almost
everywhere for each c>0. No pole or contact is removed.

Finally the source density depends on tau and m only through
tau/m^2. Substitution Omega=m omega gives
K_m(t)=m K_1(mt). Changing time variable shows the full
half-line L1 norm is independent of every fixed m>0.
The same bound applies to m=1000 in the current parent.
