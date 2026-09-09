# Two-time transfer and both exact front amplitudes

S6.90 gives Y=S U, U=C V, C=1+Z/k and

    V'=i k Lambda V+E_k V/k,
    ||E_k||_infinity<10^38,
    ||Z||_infinity<10^21, ||S||,||S^-1||<1000.

Here k>=K0=4*10^31, the interval length is at most T=10^-7,
and Lambda is real diagonal on I. The leading diagonal
transport is exactly zero, not omitted. Set eta=T*10^38/k<=1/4
and epsilon=10^21/k<1/2. If V(t,s) denotes its fundamental
matrix and D(t,s)=exp(i k integral Lambda), D has norm one.
Duhamel in the interaction picture and the exponential series give

    ||V||<=exp(eta)<2, ||V-D||<=2eta,
    ||C||<3/2, ||C^-1||<=2, ||C^-1-1||<=2epsilon.

The physical transfer and its leading term are

    T_Y=S(t) C(t) V(t,s) C(s)^-1 S(s)^-1,
    T_Y^0=S(t) D(t,s) S(s)^-1.

The source endpoint inverse cannot be dropped. The exact
noncommutative decomposition

    C_t V C_s^-1-D
      =(C_t-1)V C_s^-1+(V-D)C_s^-1+D(C_s^-1-1)

bounds its norm by 6epsilon+4eta. Hence

    ||T_Y-T_Y^0|| < 10^39/k.

Both endpoint volume factors are below four and N e^3<16.
The exact observable reconstruction adds 1/k. Therefore the
scalar source multiplier has the uniform absolute bound

    |G_hat(t,s;k)-G_hat^0(t,s;k)| < 10^44/k^2.

It applies to every ordered pair s<t in I, regardless of the
oscillatory phase. There is no claimed interacting EFT domain
extending to these arbitrarily large classical frequencies.

## Mode pair calculation

Write omega_c=N c_clock/(eR), omega_m=N/(eR), with positive
kinetic diagonals kappa_c,kappa_m and positive ell. The exact
action-normalized columns are the S6.90 columns

    (v_j, sign*i*omega_j K v_j)/sqrt(2omega_j*kappa_j),
    v_c=(1,-ell), v_m=(0,1).

fronts.py multiplies the complete target basis, both signed
phase pairs and the source inverse basis. Let

    S_j(t,s)=integral_s^t omega_j(u) du.

The result, including the physical probe prefactor, is

    G_hat^0 = A_c(t,s) sin(k S_c)/k
              + A_m(t,s) sin(k S_m)/k,
    A_c = (R_s/R_t)^(3/2) N_s e_s^3 ell_s ell_t
          /sqrt[(omega_c kappa_c)_s (omega_c kappa_c)_t],
    A_m = (R_s/R_t)^(3/2) N_s e_s^3
          /sqrt[(omega_m kappa_m)_s (omega_m kappa_m)_t].

All roots are the positive actual real-time roots. Native algebra
checks the two-time expressions, the sine rather than cosine,
both endpoint density factors and the same-time inverse kinetic
residue limit. The two amplitudes are smooth and strictly positive.

The actual parent bounds give R,N,e in (0.99,1.01), ell>1/25,
omega_j<4, kappa_c<200 and kappa_m<100. The volume ratio is
above 1/4 and N_s e_s^3>1/2. It follows that

    A_c > 1/10^8, A_m > 1/10^4.

The lower kinetic and frequency bounds also give each amplitude
less than 10^6. These deliberately loose bounds are enough for
non-cancellation; no numerical amplitude optimization is claimed.
