# Exact paired external currents

Use all-outgoing scalar momentap_i and outgoing nullk, in the fixed
physical TT gaugeepsilon_0mu=0,k.epsilon=0,tr epsilon=0.
LetJ_i=p_i.epsilon.p_i/(p_i.k).
The scalar emission factor is exactly-J_i, not merely its leading
smallk approximation: the terms containingk ortr epsilon vanish.

For a partitionL=ij,R=lm setT_L=T(p_i,p_j),
Hnum_L=eta T_L eta-eta tr(eta T_L)/2, and similarly forR.
SetD_L=(p_i+p_j)^2,D_R=(p_l+p_m)^2. Trace-reversal reciprocity gives

N=pair(T_L,Hnum_R)=pair(T_R,Hnum_L).

The full sixteen-component identity

T(p_i+k,p_j;1)=T(p_i,p_j;1)+T(k,p_j;0)

then rewrites the four original external diagrams exactly as

-N[(J_i+J_j)/D_R+(J_l+J_m)/D_L]
-sum_i J_i*pair(Hnum_other,T(k,p_partner;0))/D_other.

The second line is mandatory. All18 exact comparisons
(three rational physical states, two TT polarizations, three channels)
have a nonzero second line, providing explicit omission controls.
The general tensor identity, not these samples, proves the regrouping.

For a future mass-one scalar withE<=2, writeJ=N_p/(omega*d_p),
whereN_p=pvec.epsilon.pvec andd_p=E-pvec.n.
On the convex massive momentum ball,
|N_p|<=4, |N_p-N_r|<=4|pvec-rvec|,
d_p>=1/4 and|d_p-d_r|<=2|pvec-rvec|.
The energy gradient has normbelow1; the TT operator norm is at most1.
The elementary quotient-difference identity therefore gives

|J(p)-J(r)|<=144|pvec-rvec|/omega.

BecauseJ(-p)=-J(p), a mixed incoming/outgoing pair is precisely this
difference. At its Born transfer tau_j, recoil gives
|pvec-rvec|<=sqrt(tau_j)+3omega. Cauchy-Schwarz yields
(sqrt(tau_j)+3omega)^2<=10(tau_j+omega^2), hence

|J_i+J_j|<600sqrt(tau_j+omega^2)/omega.

Foromega<=sqrt(delta)/192<=sqrt(tau_j)/192 the same estimate gives
|J_i+J_j|<300sqrt(tau_j)/omega.
There is no cancellation assumption between different hard channels.
For the timelike channel the old sum_i|J_i|<=64/omega suffices.
