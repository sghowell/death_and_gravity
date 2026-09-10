# Full regulated inserted tadpole

Use e=epsilon, Q=16pi^2, C=2NY/Q, T=4m^2 and mu=m.
The unchanged D-dimensional inner on-shell subtraction gives

    P_D(q)=Integral_T^infinity rho_D(v)/[(v-1)^2(q^2+v)] dv,
    rho_D(v)/C=exp(gamma_E e)mu^(2e)4^e sqrt(pi)
       v^(1-e)(1-T/v)^(3/2-e)/(2Gamma(3/2-e)).

This is a regulated representation of one perturbative insertion,
not an exact normalized spectral measure. Pair the inner mass and
residue references before integrating the remaining loop.

The regulated tadpole is Integral rho_D Tad_D(v)/(v-1)^2,
where Q Tad_D(v)=exp(gamma_E e)mu^(2e)Gamma(e-1)v^(1-e).
With z=T/v and r=1/T,

    T_insert,D/(C/Q)=T Gamma(e-1)H(e)
       Integral_0^1 z^(2e-2)(1-z)^(3/2-e) A(rz) dz,
    H(e)=exp(2gamma_E e)4^(-e)sqrt(pi)/(2Gamma(3/2-e)),
    A(k)=(1-k)^(-2).

The representation initially converges for 1/2<Re(e)<1;
the beta formulas and subtracted remainder then give its meromorphic
continuation near zero. Keeping only A=1 would omit a finite term.
Split A(k)=1+2k+h2(k), h2=k^2(3-2k)/(1-k)^2.

Let F0(e) be the exact leading outer reference of S6.143. The
first two beta integrals are exactly

    T R(e) F0(e)-2F0(e), R(e)=(3/2+e)/(1-2e).

This follows by three explicit Gamma recurrences, with no change
to special-function evaluation. Since
F0=1/(2e^2)-7/(6e)+47/18+pi^2/12+O(e),

    T R F0=T[3/(4e^2)+1/(4e)+13/4+pi^2/8]+O(e),
    -2F0=-1/e^2+7/(3e)-47/9-pi^2/6+O(e).

The remainder still has a nonzero simple pole: its compact integral
is finite at e=0 but Gamma(e-1)H(e)=-1/e+4log2-3+O(e).
Keep the product with its first e derivative before MS subtraction.
The finite remainder is

    delta_T=T Integral_0^1 (1-z)^(3/2)/z^2 h2(z/T)
            [4log2-3-2logz+log(1-z)] dz.

Thus T_insert,MS=(C/Q){T(13/4+pi^2/8)-47/9-pi^2/6+delta_T}.
The complete remainder and endpoint proof are in remainders.md.
The large tadpole is not set to zero because its on-shell affine
momentum remainder vanishes.
