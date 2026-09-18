# Excluded-event singular layer and rare-index control

For M_eta=deltaa_eta*ln(x-R_eta), the normalization term in the
density split is bounded by 2100*a*eta*L_eta/kappa using the
R-weighted logarithmic moment and r<=a*eta/x.

On A_eta minus A the omitted tail satisfies T>y_eta=x-R_eta.
Before conditioning T is independent of the retained process, and
E[T]=a*eta; hence its conditional-on-retained probability is
<=min(1,a*eta/y_eta). No independence after the total cut is asserted.

If R_eta<=x/2, the mark is bounded by
1400*(x/2)*[ln(1/x)+ln 2]/kappa. The tail probability is <=2a*eta/x.
This gives <=1400*a*eta*L_eta/kappa.
If R_eta>=x/2, the S330 renewal equation implies
P_eta'(s)/P_eta(x)<=a/s for almost every s<=x.
The factor s in the mark cancels. The high part is therefore bounded by

1400*a/kappa * integral_0^(x/2) (-ln y)*min(1,a*eta/y) dy.

The zero-count atom contributes zero because the Born-subtracted
index mark vanishes there. The positive energy endpoints have no atoms.

For c=a*eta and b=x/2 the integral equals
c*[1-ln c+((ln c)^2-(ln b)^2)/2] if c<=b,
and b*(1-ln b) otherwise. In either case it is at most
c*[1+ln(1/c)+(ln(1/c))^2/2].
Write L=ln(1/eta), z=ln(1/a). Since a*z<=1/e<1/2 and
a*z^2<=4/e^2<1, expansion gives

a*[1+L+z+(L+z)^2/2] <=2+(3/2)*L+L^2/2 <=2*(1+L)^2.

Thus the high part is <=2800*a*eta*L_eta^2/kappa.
Adding low, high and normalization budgets yields
(1400+2800+2100)*a*eta*L_eta^2/kappa.
At a=0 the formula with c=0 is defined by its zero limit, without
evaluating ln 0. This proof is uniform at the original rare index;
it does not leave a spurious ln(1/a) in the final error bound.
