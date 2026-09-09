# Complete complex-frequency evolution

The universal Euler-first lapse reduction and weighted canonical source
identities are replayed. They apply after inserting the actual new
a,m,g,r,r_N,h_lapse,ell,alpha,beta,H. Independently, the complete moving
gradient, including B' and the new lapse velocity, is checked against
S6.96's fresh rational clock speed. No frozen old-background bound is
used in this bridge.

For the original phase X=(v,chi,p,P_chi),

    X'=M_X(u,q) X + N e^3 e4 J,     q=k^2/R^2,

where M_X is polynomial in q of degree at most two. The complete packet
map D, with the moving momentum shift and both volume factors retained,
gives for nonzero k

    Y'=(k J0+L0+L1/k+L2/k^2+L3/k^3)Y + R^(3/2) N e^3 e4 J,
    G_hat(t,s,k)=(R_s/R_t)^(3/2) N_s e_s^3 T_Y(t,s,k)[1,3]/k.

The original formula G_hat=N_s e_s^3 U_X(t,s,k)[1,3] defines the same
response without a singular chart. The source comes from the physical
volume probe and the output is chi itself. Indices above are zero-based.

For a complex scalar k with |k|>=1, put Z=S^-1 Y. Its COMPLETE generator
is i*k*Lambda+B(k), Lambda=diag(omega_c,omega_m,-omega_c,-omega_m),

    B(k)=S^-1*(L0+L1/k+L2/k^2+L3/k^3)*S - S^-1*S'.

All frequencies are real on I and 0<omega_c<omega_m. The Hermitian part
of i*k*Lambda is exactly -Im(k)*Lambda. Thus real frequency creates no
norm growth, and the largest Hermitian eigenvalue is at most
|Im(k)|*omega_m+||B(k)||_2. The fresh bounds give

    ||B(k)||_infinity < 5*4*10^4*8+5*2*10^8 < 2*10^9,
    ||B(k)||_2 < 4*10^9,
    ||S||_2<16,    ||S^-1||_2<10.

Euclidean logarithmic-norm Gronwall therefore gives

    ||T_Y(t,s,k)||_2 <=160 exp(S_m(t,s)*|Im(k)|+4*10^9*(t-s)).

The exact scalar endpoint factor is below two, so |G_hat| is bounded
by 320 times the same exponential divided by |k|. No smallness of B,
frequency gap, interaction-picture exponential estimate, or truncated
normal form is needed.

For |k|<=1, the actual |R|>.99 gives |q|<4. The ORIGINAL polynomial
generator then yields |G_hat|<=2*exp(T*28169/64)<4, using
exp(x)<=1/(1-x) for 0<=x<1. This bound applies even when the underlying
complex spatial vector is arbitrarily large. The two scalar-root
domains cover every frequency and agree where they overlap.
