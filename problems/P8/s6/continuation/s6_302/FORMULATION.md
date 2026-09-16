# Complete known finite minimal-gravity one-loop reference

## Unchanged inputs and stated sector

Keep the complete original QG2-H8A420 action and all fixed profiles.
In mass-one units kappa=10^800, n=10^200/512+2, g=1/8192.
Do not replace the current action by a minimal scalar model: the
minimal Einstein/light-massive-scalar one-loop sector is a specified
component, not the entire current-action amplitude.

Use the S283 raw scalar masters, S284 full-D cut coefficients and Gram
completion, and S288 physical-pole completion. Work atD=4+2epsilon.
The light pole mass and external residue are on shell. Retain the
original raw Gamma_E,4pi and scale conventions. Remove the tagged local
UV pole before dividing the S278 analytic soft factor.

At analytic reference resolution/nu=1 the selected finite representative is

A_finite = [F_known + alpha(s^2+t^2+u^2) + beta*mu^2
             +16pi^2 delta_kappa*T0]/(16pi^2 kappa^2).

T0 is the S288 Newton shape; the tree is -T0/kappa.
The complete explicit F_known is in notes/assembly.md.
alpha,beta,delta_kappa are required, unassigned finite matching
coordinates. Pole subtraction defines their coordinate convention,
not their physical values. The known piece alone is convention-dependent.

## Physical compact domain and quantitative claim

For the bound choose mu=nu=1,25/4<=s<=16,t,u<0,s+t+u=4,
0<delta<=1 and min(-t,-u)>=delta. No exact forward or backward point
belongs to a fixed window. Both physical Coulomb phases and full
crossing are retained; no real-part-only master approximation is used.

The written termwise contour and coefficient estimates give

|F_known| <10^11(1+|ln(delta)|)/delta^2,
|A_known| <10^9(1+|ln(delta)|)/(kappa^2 delta^2).

At the explicit comparison window delta=10^-204, original matter
Born exceeds10^-600 and full Born A_m+A_G is positive and at leastA_m.
Thus |A_known/(A_m+A_G)|<10^-580 and the absolute known interference
2Re[A_known/(A_m+A_G)] is below10^-579.

The window is not a UV cutoff or a choice of physical matching.
It is not an estimate uniform asdelta tends to zero. The analytic
reference resolution/nu=1 is not asserted to be the physical detector
range used by the previous inclusive-rate certificates.

## Closure boundary

Unknown matching pieces, the remaining current-action hard terms,
gravity-Born radiation and higher-loop/operator errors remain separate.
A finite amplitude coefficient is not a quantum S-matrix construction,
a Regge completion, a positive-unitarity theorem, or a cosmological
bounce. Original V/G/B/P8 and all historical qualifications are unchanged.
