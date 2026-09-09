# An opposite-sign matrix Fourier estimate

Let A=-ik omega+B, alpha>=0, k>=32, and s=alpha+k.
The diagonal leading matrix has singular values at least
alpha+0.99k. Since ||B||<=2 and 0.99*32>4, a Neumann estimate
gives Q=(A-i alpha I)^(-1), ||Q||<=3/s. Differentiating the
inverse identity, with the whole-interval derivative bounds, gives

    ||Q^(j)|| <= D_j/s, D=(3,18,351,11673), j<=3.

Explicitly D_0=3 and
D_n=3 sum_(r=1)^n binom(n,r) C_r D_(n-r), with
C_r=ceil(omega_r+B_r/32). Here k/s<=1.
This inversion separates opposite signs; it introduces no
clock-minus-matter denominator.

For a compactly supported row a=g l, let a_0=a and
a_(n+1)=-(a_n Q)'. Since
(e^(-i alpha u)F)'=(A-i alpha I)e^(-i alpha u)F,
three exact integrations by parts leave no boundary terms.

If ||l^(j)||<=L_j, define
c_(0,j,r)=binom(j,r)L_(j-r) and
c_(n+1,j,r)=sum_(t=0)^(j+1) binom(j+1,t)
c_(n,t,r)D_(j+1-t), with absent entries zero.
The ordered product rule and submultiplicativity give the bound
s^(-3) sum_r c_(3,0,r)||g^(r)||_1 for the remaining row integral.
A native independent symbolic expansion checks every sampler
derivative coefficient of this recurrence.

For the ORIGINAL O, after the common k^(1/2) factor in derivative
modes, the spatial row is a^(-5/2)(-ell,1). The time row is
-i f omega+k^(-1)(f'+f B), f=a^(-3/2)(-ell,1).
Whole-interval bounds on these COMPLETE rows are

    L_time=(3,2,18,108), L_space=(2,1,13,13).

Thus the two sampler combinations are
T_g=sum_(j=0)^3 T_j ||g^(j)||_1,
S_g=sum_(j=0)^3 S_j ||g^(j)||_1, with

    T=(619407,63666,3078,81),
    S=(405837,41877,2025,54).

The initial configuration factor has norm <2, and exact
approximate evolution has norm <=exp(2h).
The joint positive-frequency integral is exactly

    integral_K^infinity integral_0^infinity
        k^3/(alpha+k)^6 d_alpha d_k = 1/(5K).

It converges at BOTH endpoints of the high-frequency region.
All three lower-order transport derivatives have been retained.
