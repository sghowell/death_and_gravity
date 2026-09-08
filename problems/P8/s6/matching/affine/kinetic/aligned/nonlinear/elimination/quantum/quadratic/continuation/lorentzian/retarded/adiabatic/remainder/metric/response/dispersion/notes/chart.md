# Pullback to the regular metric chart

The S6.57/S6.60 chart is

    physical N=1+n,
    physical log-scale= v+omega(N,u),
    omega=-log[(h-1+N^-2)/h]/4, h=(1+u^2)^3.

At N=1 its first and second lapse jets are

    delta=1/(2h), omega_NN=-3/(2h)+1/h^2,
    T=[[1,0],[delta,1]].

For a functional whose background metric currents are (q_N,q_Z),
the second variation in this nonlinear chart is

    T^T K T + diag(q_Z omega_NN,0).

The operator notation T^T K T means that the input is multiplied by
the time-dependent T before K acts. Derivatives of delta are
retained; this is not constant-coefficient matrix multiplication
except for the highest derivative symbol. The displayed extra
term follows by the functional chain rule and is local.

The vector's background current is (-rho_sigma,3p_sigma).
The S6.68 fixed-tadpole action has exactly the opposite gradient.
Consequently the two nonlinear-chart contacts cancel. Omitting the
vector-only contact before including the tadpole would give an
incorrect intermediate Hessian: at h=1 and p_sigma=1 it is -3/2.
No profile, state covariance or subtraction value is reselected.

For the background-cancelled total, physical output normalization
has no remaining background-gradient terms. Its regular force is

    (F_n,F_v)=(3delta delta_p-delta_rho, 3delta_p).

This reproduces the literal S6.57 source convention. The tree
inverse there solves E_tree y+F=0. Thus a formal first feedback
application is G K_total T, with the sign already contained in G;
one should not introduce an extra minus sign.

## Continuous multiplication and one retarded response

Let the regular input sources be smooth, compact in time after
u0=-1/2, and zero on an initial neighborhood. Multiplication by T
preserves that preparation. On I=[-1/2,1/2] define the joint source
norm using both fields and all derivatives zero through ten.

Leibniz' rule gives

    ||T y||_C10 <= C_T ||y||_C10,
    C_T=max_{k<=10}[1+sum_{j<=k} binom(k,j)||delta^(j)||_infinity]
       <=46090764897/8.

Each derivative majorant is a rational positive-denominator bound
on the whole interval; no grid interpolation is used.

Let kappa be S6.68's complete vector-plus-fixed-tadpole physical
C10-to-C0 norm. For unit regular input the stress is bounded by
eta0=C_T kappa. S6.57's original zero-initial forced tree equations
then give phase, lapse, matter and physical-logscale C0 bounds
131 eta0,94 eta0,44 eta0,225 eta0 respectively. They are each
below 10^-25 at L=10^24,m=1000. The source normalization lies within
the frozen small-source gate; linearity extends the displayed
operator coefficients to general amplitudes of the linear problem.

These estimates do not bound the lapse derivative: the old lapse
rate estimate also requires a stress derivative eta1, which S6.68
does not provide. Nor does a C10 input/C0 output estimate license
iteration on the same Banach space. No g^2 residual estimate is
claimed for y0+g G K y0, since it would require control of K G K y0
in the original source domain. The matter invariant retained is
a^3(p_m+3ell v), not p_m separately.

The finite local fourth-derivative matrix in this chart is

    [[-7357/(6561h^2),-562/(243h)],[-562/(243h),-4]].

Its nonzero lapse entry shows why bounded low-derivative local
coefficients alone are not a no-loss inverse argument. It is not
a full nonlocal instability claim.
