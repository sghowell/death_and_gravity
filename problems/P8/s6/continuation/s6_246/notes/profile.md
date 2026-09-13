# Full profile restoration and nonlinear common clock

Let \(\rho(t),P(t)\) denote the complete fixed actual heavy reference integrals, not freely adjustable functions. Their numerical uniform bounds are inherited from S6.240, and their values are held fixed under every subsequent variation.

## Entire physical profile Hessian

The literal full fixed coefficient has clock jets
\[
 A=-P,\qquad B=-(\rho+P)/2,\qquad \partial_X^2(\kappa\Delta F)|_{X=1}=0.
\]
The complete high-order switch and extra mass-one vacuum term have their required vanishing jets. The packet differentiates the full frozen coefficient to check these facts; it does not replace the coefficient away from the clock.

Thus the physical density's complete quadratic jets are those of
\[
 \mathcal P=a^3N e^{\operatorname{tr}Q/2}[A+B(N^{-2}-1)].
\]
The one-current is
\[
 P_1(D)=a^3[\rho n_D-P\tau_D/2]=-J_H(D).
\]
The entire mixed Hessian is
\[
 P_2(D,G)=a^3\left[-(\rho+P)n_Dn_G
 +\frac{\rho}{2}(n_D\tau_G+n_G\tau_D)
 -\frac P4\tau_D\tau_G\right].
\]
There is no shift term. Its spatial restriction is
\(P_{QQ}=-a^3P\tau_D\tau_G/4\).

The fixed profile is not separately on shell as a metric-only sector with the clock held fixed. Consequently apply Gaussian Ward to
\[
 R_H(D_{\rm syn},G_{\rm syn})
 =R_{245}(D_{\rm syn},G_{\rm syn})-P_{QQ}(D_{\rm syn},G_{\rm syn}),
\]
and only then restore \(P_2(D,G)\). The full actual physical form is
\[
 R_{\rm ADM}=[R_{245}-P_{QQ}](D_{\rm syn},G_{\rm syn})
       +W_s(E)+W_d(E)+P_2(D,G).
\]
This removes/restores exactly one profile contribution. No second finite heat action is added.

## Full nonlinear clock, including general symmetric directions

The unchanged parent has
\(R|_1=1,\ R_X|_1=2\delta,\ R_{XX}|_1=0\), with
\(\delta=1/[2(1+t^2)^3]\). The code checks these jets from the FULL parent function. For an arbitrary symmetric hatted \(Q\),
\[
 Q_{\rm phys}=\widehat Q-\tfrac12\log R(t,N^{-2})I.
\]
The first map is
\[
 T(n,\beta,Q)=(n,\beta,Q+2\delta nI),
\]
and the full mixed second map is
\[
 Q_{\rm phys}^{(2)}=2(4\delta^2-3\delta)n_Dn_GI.
\]
The nonlinear Gaussian contact is therefore
\[
 C_H=3a^3P(4\delta^2-3\delta)n_Dn_G.
\]
The profile one-current gives \(C_{\rm prof}=-C_H\). Retain both before their full-reference cancellation. In particular
\[
 R_{\rm clock}(D,G)=R_{\rm ADM}(TD,TG)
\]
holds only after the complete mean/contact assembly, not after declaring either mean to vanish separately.

## Literal full clock-density check

Differentiating
\[
 a^3N e^{\operatorname{tr}\widehat Q/2}
 R(t,N^{-2})^{-3/4}[A+B(N^{-2}-1)]
\]
gives one-current
\(a^3[(\rho-3\delta P)n_D-P\tau_D/2]\) and mixed Hessian
\[
\begin{aligned}
 a^3\{&
 -[(1-6\delta)\rho+(21\delta^2-9\delta+1)P]n_Dn_G\\
 &+(\rho-3\delta P)(n_D\tau_G+n_G\tau_D)/2
 -P\tau_D\tau_G/4\}.
\end{aligned}
\]
It equals \(P_2(TD,TG)+C_{\rm prof}\) exactly. Setting \(Q_D=2v_DI,\ Q_G=2v_GI\) recovers the frozen S6.241 scalar coefficients, with no new quadratic retuning.

## Original projected mean

For the computational cutoff, the Gaussian pressure is \(P_K\) but the profile remains fixed at \(P\). The retained total clock contact is consequently
\[
 C_K=3a^3(P_K-P)(4\delta^2-3\delta)n_Dn_G,
\]
not zero. An exact control with unequal means makes this term nonzero. Its complete tail bound is included in the cutoff note. No state or finite profile is adjusted as \(K\) changes.
