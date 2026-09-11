# Momentum-uniform finite derivative bounds

The domain is a real smooth compactly time-supported
homogeneous tracefree symmetric gamma, zero near the
initial Cauchy surface, with every raw derivative
||gamma^(j)||op<=delta=1/100 for j=0..12 on
I=[-1/2,1/2]. No analyticity or commuting time directions
are assumed. This is a specified finite-jet neighborhood,
not every bounded-amplitude shear. The background is
a=(1+t^2)^2, m=1000.

The ordered-simplex formula for a Frechet derivative
of exp(gamma) has n! terms of simplex volume1/n!.
Its norm is at most exp(delta) times the product of
direction norms. The ordinary time derivative sums
over set partitions, giving exp(delta) B_j(delta),
where B_0=1 and
B_n=delta sum_(j=0..n-1) binom(n-1,j) B_j.
Conjugation by exp(-gamma/2) costs exp(delta).
Thus the relative matrix exponential derivatives are
bounded by e_0=1 and e_j=2 B_j for j>=1,
since exp(2delta)<=1/(1-2delta)<2. The same argument
applies to exp(-gamma).

For p=1,2,3 write (a^-p)^(n)/a^-p=P_(p,n)/(1+t^2)^n.
The exact polynomial recursion is

    P_(p,0)=1,
    P_(p,n+1)=(1+t^2)P_(p,n)'-2(n+2p)t P_(p,n).

The sum of absolute polynomial coefficients weighted by
(1/2)^degree bounds it on I. Denote this c_(p,n).
Loewner bounds for the two positive pieces of K give

    ||K^-1/2 K^(n) K^-1/2||
    <=q_n=max(c_(3,n),
              sum_j binom(n,j)c_(1,n-j)e_j).

The normalized omega-squared derivative is bounded by
the analogous p=2 sum. Direct first derivatives sharpen
these to q_1=5 and s_1=4:
3|H|+2delta<5 and 2|H|+2delta<4, using |H|<=8/5.
The code records all finite rational constants.

Crucially these are not estimates involving cond(K).
Let M_n=B^-1 B^(n), M_0=I. Differentiating B^2=K n
times and balancing by B^-1 on both sides gives

    M_n+M_n^T =
      K^-1/2 K^(n) K^-1/2
      -sum_(j=1..n-1) binom(n,j) M_j M_(n-j)^T.

In an orthonormal B eigenbasis with eigenvalues b_i>0,
the solution is (M_n)_ij=b_j/(b_i+b_j) times this
symmetric right side. Its Frobenius norm is no larger
than that of the right side and hence no larger than
sqrt(3) times its operator norm. This proves the safe
recurrence

    ell_0=1,
    ell_n=2(q_n+sum_(j=1..n-1)binom(n,j)ell_j ell_(n-j)).

No uniform lower bound on a ratio of distinct b_i is
required. This remains uniform at arbitrarily large k.

For N_n=(B^-1)^(n) B, the identity B^-1 B=I gives
N_n=-sum_(j=0..n-1)binom(n,j)N_j M_(n-j).
Leibniz then bounds L^(n) by
sum_j binom(n,j)||N_j|| ell_(n-j+1).

Likewise, differentiating omega^2 and omega^-1 omega=1
gives the positive majorant recurrences recorded in
jets.py. Finally rho^(n) is bounded by half the
Leibniz sum of normalized omega inverse and omega
derivatives. These bound R^(n) by L^(n), and S^(n)
by L^(n)+rho^(n). In particular R_0<=10, S_0<=11.
Twelve prescribed gamma jets provide more than the
eleven needed for the tenth-order reference residual.
