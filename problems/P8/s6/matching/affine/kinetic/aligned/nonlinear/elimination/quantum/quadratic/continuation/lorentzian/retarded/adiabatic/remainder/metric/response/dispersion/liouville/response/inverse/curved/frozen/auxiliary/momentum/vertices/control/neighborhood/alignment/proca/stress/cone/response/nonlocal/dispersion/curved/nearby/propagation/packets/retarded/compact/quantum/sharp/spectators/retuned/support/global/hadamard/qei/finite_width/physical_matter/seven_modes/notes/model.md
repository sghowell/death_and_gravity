# The actual Hamiltonian and the constrained vector

The metric signature is -+++. Here a>0 is a general smooth scale,
k is a comoving Cartesian momentum, and m>0 is constant. Restore
m=1000 and a=(1+u^2)^2 only after deriving the sector equations.
Let J2=[[0,1],[-1,0]] and Omega6=[[0,I],[-I,0]].

The ordinary Proca Hamiltonian following elimination of W0 is

    H_P = [m^2 Pi^2/a + (k.Pi)^2/a^3
           + a W^2 + |k cross W|^2/(m^2 a)]/2.

Every summand is nonnegative. Its momentum and coordinate Hessians are

    A=m^2 I/a+kk^T/a^3,
    B=a I+(q I-kk^T)/(m^2 a),
    M_P=[[0,A],[-B,0]],  A B=(m^2+q/a^2)I.

The product identity is instantaneous algebra, not an assertion that
time-dependent evolution has constant-frequency solutions. The full
evolution uses M_P(u,k). The actual temporal constraint is

    W0=-i k.Pi/a^3.

Differentiating it with the full Hamiltonian gives
W0_dot+3 H W0-a^-2 i k.W=0. Thus the missing fourth configuration
component is reconstructed, and the covariant divergence constraint
holds. This is the actual positive massive one-form theory. Its
three physical polarizations are not four unconstrained oscillators.

For either tensor polarization, normalized by E_ij E_ij=2,

    L_T=a^3 T_dot^2/4-a q T^2/4,
    PT=a^3 T_dot/2,
    H_T=PT^2/a^3+a q T^2/4,
    M_T=[[0,2/a^3],[-a q/2,0]].

The canonical minimal scalar variables are chi=T/sqrt(2),
p_chi=sqrt(2)PT. Dropping these factors changes both the CCR and
the physical S6.102 stress contraction.

The scalar block is the complete S6.98 regular density generator,
with the original retuned functions substituted and its physical
wave number replaced by q/a^2. No large-k truncation or pole
division by the clock-crossing function is used. Assemble

    M7=diag(M_scalar,M_T,M_T,M_P),
    Omega7=diag(Omega_scalar,J2,J2,Omega6).

The native 14 by 14 identity M7 Omega7+Omega7 M7^T=0 holds at
all finite u and Cartesian k. It proves exact preservation of the
joint density CCR under the fundamental solution. The positive
action normalization is still hbar/kappa, not unity by reassignment.
