# Quantitative bounds on the actual unchanged sampled Gramian

Let S(t,t*) be the complete original canonical scalar flow and E_sel(t)
the original positive preparation matrix on the fixed bump support.
The S251 initial Gramian is M*=integral f(t)^2 S(t,t*)^T E_sel(t) S(t,t*) dt.
At any target u, rewriting this SAME quadratic form gives
M_u=S(u,t*)^(-T) M* S(u,t*)^(-1).
Further exact canonical coordinate changes are congruences of M_u;
they never minimize again or change its Cauchy state.

Between a sampling time and |u|<=T=10^-60, there is exactly one
chart switch and total duration at most1/2+T.
The absolute root-energy rate12 and the complete conversion cost16000
give both-direction root-energy factor bounded by
2*3^6*16000=23328000<Gamma=10^8.
Here exp(6)<3^6 and exp(12T)<2 follow from e<3 and
12T<1/2 with the positive exponential series.
These estimates use the full finite-P equations and nonzero profiles.

The ORIGINAL selection form is a^3(E_y+y^T K y/2), with its
unchanged mass1 term. It lies between E_y and16 E_y on the support:
a^3<=8 and the added K term is at most the main qG contribution.
Using E_y=P r^T Q r/2 and Q in[10^-4,100], integration and the
both-direction transport bound imply
m I<=M_r/[P integral f^2]<=L I,
m=10^-4/Gamma^2=10^-20,
L=16*100*Gamma^2=1.6*10^19.
The positive integral is a common scalar. The selected covariance
is homogeneous of degree zero in the positive Gramian, so this
normalization cancels exactly.

For any positive real symmetric full4x4 M the complete spectral
formula for the selected pure covariance is
V=(1/2) M^(-1/2) |i M^(1/2) Omega M^(1/2)| M^(-1/2).
It is the same S251 two-mode covariance, including mixing and
degenerate preparation frequencies. Put X=2V>0.
Since the Hermitian modulus squared is
-(M^(1/2) Omega M^(1/2))^2, multiplication gives
X M X=-Omega M Omega.
For m I<=M<=L I and orthogonal Omega,
m X^2<=X M X=-Omega M Omega<=L I.
Hence ||V||<=sqrt(L/m)/2<10^21.
This is a full positive-matrix argument, not a scalar-mode formula.
Symplectic congruence of the minimizer proves that V(M_r) is exactly
the transported original covariance. Rescaling M by a positive scalar
does not change it. Both facts are tested by independent spectral
linear algebra as well as exact noncommuting fixtures.

For each original tensor, the balanced energy metric is diag(a,a^-3),
with eigenvalues in[1/8,2]. Its absolute root rate is at most6.
There is no chart switch. Both-direction path factor
Gamma_T=2*3^3=54 bounds the original sampled transport.
The unchanged mass1 selection is at most twice its energy.
Thus its normalized full2x2 Gramian has eigenvalues in
[(1/8)/54^2,4*54^2], and the same covariance argument gives
V_rT<10^5 I for each polarization.

All four physical modes, including both coupled scalars and both
tensors, are retained. V_x=P V_r restores the weighted CCR factor.
The original H/Proca product states and their full determinants stay
fixed; the theorem does not compute their interacting mean.
