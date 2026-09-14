# Complete physical rows, both normalization factors and finite jets

The S263 evaluated bound on d(a^3 S/P^2)/du is
48640097818.79867...<10^11, where S=(B+B^T)/2 is the full
central reference symmetric boundary. At the bounce S=0, but S'
is not zero. The exact rational bound, not this printed decimal,
and |u|<=10^-60 imply |Sbar_ij|<10^-49 for Sbar=a^3 S/P^2.

For x=(P y_b,P y_m,Pi_b,Pi_m), restore the entire old momentum
with R=[[I,0],[P^2 Sbar,I]] and the exact original central swap.
The complete old canonical vector is T_central^(-1) R diag(P,P,1,1)^(-1)x.
The packet composes this with the original S261 coordinate/lapse
rows that survived S262's physical-binding erratum.

After removing ONLY the common1/sqrt(kappa), the hat log-scale row is
(Sbar00/(2aP),Sbar01/(2aP),1/(2aP^2),0).
For the lapse let p_b_old=Pi_b+P^2(Sbar00 y_b+Sbar01 y_m),
p_m_old=Pi_m+P^2(Sbar01 y_b+Sbar11 y_m). Its complete scaled value is
[-2Theta P^2 y_b/a^2+3Theta ell y_m
 -(2E P^2/a^2+3Tcorr)p_b_old/(2aP^2)
 +ell E p_m_old/a^3]/(2J).
The independent full constraint reconstruction checks every term.
The weighted CCR check is vrow (P Omega) nrow^T=Theta/(2J a^3),
before restoring the explicit1/kappa. These physical rows need not commute.

The actual tiny-slab source re-enters the box
a in[1,2], J in[1,2], E in[-1/2,-1/4], |Theta|<=5T,
0<=ell<=1/10, |Tcorr|<=10^-380 and |Sbar_ij|<=10^-49.
The pivot range is the previous full-profile re-entry result,
not a new bare substitution. Theta's bound retains both terms;
the original profile envelope bounds Tcorr.
All eight row-entry absolute intervals are stored. Their all-entry
sums give safe operator norms10^-112 and10^16.

With the COMPLETE scalar covariance V_x<=10^21 P I and both
1/sqrt(kappa) factors, integrate against the original radial measure.
The exact common integral is
I_B=integral_LOW^HIGH P^3 dP/(2 pi^2)=(HIGH^4-LOW^4)/(8 pi^2).
An outward rational ceiling is(HIGH^4-LOW^4)/72, using pi>3.
Thus Var(v_hat)<=10^-224*10^21*I_B/kappa<10^-740,
and Var(n)<=10^32*10^21*I_B/kappa<10^-490.
Positivity of the symmetrized covariance gives
|Cov_sym(v_hat,n)|<=sqrt(Var(v_hat)Var(n))<10^-615.
The cross term is bounded, not assumed zero.

Each original canonical tensor coordinate has variance at most10^5/P.
Restoring gamma=2h_canonical/sqrt(kappa), two polarizations and the
norm-one TT projector gives total metric Frobenius variance at most
2*10^5*(HIGH^2-LOW^2)/(9kappa)<10^-660.

For each fixed point and time the self-adjoint LINEAR band lapse has
a centered Gaussian spectral measure. The scalar Chernoff inequality
and threshold epsilon=10^-230 give
Pr(|n|>=epsilon)<=2 exp[-epsilon^2/(2 Var(n))]
<=2 exp(-10^29).
The exact lower exponent from the safe variance is5*10^29.
No underflowed numerical exponential is ever interpreted as zero.

The correct U=R^(-3/4) lapse jets at N=1 are
c1=3/(2h), c2=21/(4h^2)-9/(2h), h=(1+u^2)^3.
Their literal source bridges and bounce values3/2,3/4 are checked.
Both absolute values are less than2 on the slab.
For the normalized volume exp(3v_hat)U, its complete second-order
Weyl expectation is(9/2)Cvv+3c1 Cvn+(c2/2)Cnn.
Its absolute value is below
5*10^-740+6*10^-615+10^-490<10^-489.
The linear physical volume variance is bounded by
18*10^-740+8*10^-490<10^-488.
Only these finite polynomials are claimed, not the full nonlinear
volume function on an unbounded Gaussian.
