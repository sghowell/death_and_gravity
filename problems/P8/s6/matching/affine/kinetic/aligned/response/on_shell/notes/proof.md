# A nontrivial prepared on-shell linear light family

All actions and frozen ancestors are unchanged. This is a solution of
the original linear CD/M1 light equations and its leading second-order
vector response. It is not a constructed full nonlinear parent solution.

## Actual regular physical phase system

Use incoming comoving k_in²=1/16, q=k_in²/a², u in [-1/2,1/2]. The
regular CD/M1 Hamiltonian density before the canonical swap is

    a³*(H0+R0²/J),
    R0=q*(Theta*b+Lambda*v)+w*P/2-3ell*s*Theta/2,
    H0=P²/2+q*ell*b*s+q*s²/2-q*v²-3ell²*s²/4.

The original pairs are (v,p_v=-2a³q*b), (s,p_s=a³P). With
Y=(b,s,P_b,P_s), P_b=2a³q*v, the new Hamiltonian is the substituted
density minus H*b*P_b. This moving-boundary term follows from
(2a³q)'=H*(2a³q) at fixed comoving momentum. The independent check
derives A_new=(T'+T*A_old)*T^-1 from the original phase coordinates.
Dropping T' is detected away from u=0, although it vanishes at the
single center point. The full actual system is Y'=A(u)*Y.

There is no inverse Theta. The frozen J is strictly positive, and all
denominators in this compact chart are nonzero. The symplectic form
is constant and nondegenerate. Lapse n=l_N(u)*Y is the original
n=-R0/J constraint; the matter momentum and shift constraint are
checked directly against the canonical flow.

## Geometric trace readout

The physical curvature perturbation is zeta_metric=v+n/(2h). The
exact S6.41 spatial conformal factor has omega_N=1/(2h), so the
hat-metric curvature perturbation is v. Expanding its ADM trace gives

    delta_K_hat=3v'-3H*n+q*b
               =3(Theta-H)*n-3ell*s/2+q*b.

This defines a second actual observer l_K(u)*Y. No lapse time
derivative is discarded or treated as an independent source profile.

## Nonzero prepared initial data

For an observer l_N, its first four Cauchy rows are recursively
l_(j+1)=l_j'+l_j*A. At u=-1/2, the first three rows have rank three
and the first four have rank four. Exact rational nullspace calculation
gives the unit-infinity-norm initial direction

    Y_star=(1,
      63503417382507211357087625/97994374099531076644486592,
      9248050957415283203125/1531162095305173072570103,
      72813352996311187744140625/3135819971184994452623570944).

It obeys n=n'=n''=0 initially. Its n''' and initial delta_K_hat are
both strictly positive. The actual S6.43 source therefore has

    S2'''(-1/2)=-(2/h)*n'''*delta_K_hat < 0.

Every nonzero multiple of this initial state supplies a genuinely
nonzero prepared source. An independent Cauchy recurrence for Y and
Leibniz derivatives of the lapse confirms the preparation, without
assuming the stored observer-row formulas.

## Continuous solution and jet bounds through the bounce

For each actual rational entry of A and of the observers, including
the needed time derivatives, the denominator has positive constant
and nonnegative even coefficients. Bound its numerator by its absolute
coefficient sum evaluated at |u|=1/2 and its denominator by its constant
term. Summing rows gives continuous infinity-norm bounds A_j for
A,A',A'', and observer bounds L_(N,j), L_(K,j) for j=0..3. In particular

    A_0=4771877971625/20115881984 < 238.

The bounded continuous matrix admits a unique fundamental solution on
the whole interval: its ordered integral series is dominated termwise
by (A_0*length)^j/j!. The same domination gives ||Y(u)||_infinity
<=exp(A_0)*||Y(-1/2)||_infinity. Since e<3, a sufficient exact rational
growth bound is G=3^238. The elementary inequality e<3 follows from
the factorial series and n!>=2^(n-1), strict from n=3 onward.

Define the derivative multipliers

    R_0=1, R_1=A_0, R_2=A_1+A_0²,
    R_3=A_2+3A_0*A_1+A_0³.

They follow by differentiating the actual matrix equation; the norm
bound does not assume that A and A' commute. For either observer set

    B_observer=max_(j<=3) sum_(i=0..j) binomial(j,i)*L_(observer,i)*R_(j-i),
    C=1+B_N+B_K.

All constants are exact rationals recomputed from the actual matrices,
not fitted to numerical trajectories. For any epsilon>0 choose

    rho=epsilon/(G*C), Y(-1/2)=rho*Y_star.

Then rho is strictly positive and all four lapse and trace time-jet
norms on the full interval are <epsilon. This is an explicit, very
conservative small-amplitude existence bound, not an optimal or
observationally motivated amplitude estimate.

## Transfer to the actual prepared vector response

Take the real spatial profile Y(u)*cos(k_in*x), with Y the scaled
solution just constructed. The two incoming
complex Fourier modes each carry half its amplitude. Their summed
time-jet envelopes obey the same epsilon bounds. The quadratic source
has zero and double momentum, with k_out²=4k_in²=1/4, not k_in².
Each nonzero complex source coefficient has one quarter of the
cosine-squared coefficient; the third source derivative remains nonzero.

The actual linear light evolution now satisfies the S6.43 preparation
and source-jet assumptions, with E=(136139/8)*epsilon². Applying its
retarded theorem to zero initial heavy data gives the stated physical
spatial and temporal errors on this actual leading source. The zero
output mode is the already checked algebraic temporal response.

For epsilon=1/1000 and zeta=1/10^9, E=136139/8000000, the normalized
spatial error bound is 408417/4000000000000000 and the temporal error
bound is 136139/16000000000000. Physical vector amplitudes have factor
1/tau. These are bounds, not observed deviations or a cutoff estimate.

This discharges the missing on-shell-linear source realization for
this prepared family. It does not solve the nonlinear light equations,
bound higher perturbative orders or quantum loops, prove nonlinear
secondary constraints or stability, establish an interacting cutoff,
or settle the original V/G/B UV conditions or P8 closure.
