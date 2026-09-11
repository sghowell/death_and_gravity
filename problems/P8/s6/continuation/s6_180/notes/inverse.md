# Actual Volterra normal form and constructive finite weight

Let z=(eta,w), A=-6delta^2 and gamma=1/(64pi^2 kappa).
Combining the literal current tree, exact local quantum time
channel and matched scale kernel gives

    I4(T+gamma Q)=diag(A(t),gamma F_m(partial_t^2))+V.

I4 is four zero-past ordinary time primitives. Coefficients are
not pulled through derivatives or integrals. In particular

    I4[A eta'''']=A(t)eta(t)+integral K_A(t,s)eta(s) ds,
    K_A=sum_(j=1..4) (-1)^j binom(4,j)
                   (t-s)^(j-1) A^(j)(s)/(j-1)!.

For a lower local term c_j(t) z^(j), j<=3, the exact kernel is

    K_j=sum_(k=0..j) (-1)^k binom(j,k)
                   (t-s)^(3-j+k)c_j^(k)(s)/(3-j+k)!.

All prepared endpoint terms vanish. Independent polynomial
integrals check the full variable-coefficient formulas and show
that removing the commutator changes the operator.

The matched remainder V has an integrable row-norm majorant
C(1+|log r|), 0<r<=1, and each fixed diagonal derivative has
the same class with its own finite constant. Its first row is
stronger: after separating A eta'''', it comes entirely from
lower LOCAL terms, so its kernel and first output derivative
are bounded. The scalar loop is not assigned this stronger class.

The diagonal inverse is

    B0=diag(A^-1 multiplication,gamma^-1 K_m convolution).

The inactive direction of the isolated rank-one loop is supplied
by the actual classical coefficient A, not a pseudoinverse.
Write a0=15625/6144, K=||K_m||_L1, and
beta=C(a0+K/gamma). These constants are finite. No numerical
value for C or K is asserted.

For lambda>=1, direct substitution y=lambda r and
integral_0^1 -log y dy=1 give

    integral_0^1 exp(-lambda r)(1-log r)dr
        <=(2+log lambda)/lambda <=2/sqrt(lambda).

The last inequality follows from log x<=x-1 at
x=sqrt(lambda). Weighted Young's inequality therefore bounds
||B0 V|| by 2beta/sqrt(lambda).
The explicit finite choice lambda=(4beta+1)^2 makes this
less than1/2. No source-dependent state or coefficient is chosen.

The Neumann series for I+B0V converges on the entire fixed
slab in that weighted norm. For smooth prepared adapted forcing,

    ||z||_C0 <=2 exp(lambda) max(a0,K/gamma)
                  ||I4 g_adapted||_C0.

This formula is constructive in the unevaluated majorants; it is
not a numerical or small full inverse bound. It keeps every causal
kernel contribution, including any growing response allowed by
the retained approximation.
