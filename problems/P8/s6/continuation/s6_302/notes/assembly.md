# Whole finite assembly, tagged poles and symbolic matching

Let channels(a,b,c)=(s,t,u), withc=4mu-a-b. Use all three channel rows
and all six ordered pairs. Define

V_a=a^2-4mu*a+2mu^2,
V_D(a)=V_a+2mu^2*epsilon/(1+epsilon),
T_D=sum_a N_D(a,b,c)/a=T0+epsilon*T1+O(epsilon^2),
T1=2mu^2 sum_a1/a.

Let R_D=R0+epsilon*R1+O(epsilon^2) be S288's crossing-symmetric
physical-pole completion. Its four basis elements are
mu^3 sum1/a,mu^4 sum1/a^2,mu sum bc/a,mu^2 sum(b-c)^2/a^2.
Their coefficients are

R0: (164/15,-2/15,-73/15,1/15),
R1: (-4879/225,-28/225,3128/225,-16/225).

For the S284 whole-D row denote its epsilon0 coefficients by
c00,cmm,b00,bmm and let cmm1=d_epsilon cmm_D atzero.
The full known finite expression is the sum of:

1. All ordered boxes:
   sum_(a!=b)[V_b^2 J_b ell_a/a -4mu^2 V_b J_b/a].
2. All three nonbox rows:
   sum_a[c00*C_a+cmm*(ell J_a-Q_a)/2-cmm1*J_a/2
         +b00*(2+ell_a)+bmm*(ell-L_a)].
3. The S284 compact crossed evanescent rational E_compact.
4. Minus sum_a old_gram_choice(a,b,mu,bubble=ell+2).
5. Four external legs:
   6mu*ell*T0+14mu*T0-6mu*T1.
6. Physical pole completion:
   -ell*R0+R1.
7. Analytic soft-reference division:
   +2*T1*B_soft, whereB_soft=sum_a V_a*J_a/2-mu.

The old finite Gram term in4 is indispensable. E_compact already
includes the evanescent bubble derivative minus the evanescent Gram
choice; it does NOT include the old Gram finite term.
For item5 the complete residue ratio is
r_D=(3+2epsilon)/[(1+epsilon)(1+2epsilon)]
=3-7epsilon+O(epsilon^2).
The four-leg raw term is2mu*r_D*T_D*(-1/epsilon+ell);
its derivative produces both14mu*T0 and-6mu*T1.

## Ultraviolet subtraction before infrared division

Tag the common-dimensional poles by origin before combining them.
The massive triangle and the four-leg infrared contribution give

P_IR=-2*T0*B_soft.

The remaining ultraviolet residue after the physical completion is

P_UV=-U, U=203/40*(s^2+t^2+u^2)-169/3*mu^2.

The equality of zero-momentum Cmm and Bmm Laurent functions does not
identify their infrared and ultraviolet origins. Add+U/epsilon
first. This is a pure local pole subtraction coordinate convention,
not a physical finite matching choice.

The tree amplitude is-T_D/kappa. Expanding division by the unchanged
S278 factor at x=resolution/nu=1 then adds
+2*T_D*B_soft/epsilon inside the overall1/(16pi^2 kappa^2).
Its finite part is+2*T1*B_soft. Dropping the D-tree derivative would
therefore give an incorrect known finite reference.
At another fixed x there is additionally+4*T0*B_soft*ln x;
the quantitative bound here explicitly usesx=1.

Exact symbolic extraction from the full-D coefficients and master
Laurent functions agrees with items1--7. The complete expression is
crossing-symmetric. With ell_a=ell-lambda_a its common-log derivative
is2*T0*B_soft+U, checked independently of finite matching.

## What remains independent

The complete selected finite family contains
alpha*(s^2+t^2+u^2)+beta*mu^2+16pi^2*delta_kappa*T0
in addition to F_known. The Newton sign follows by differentiating
-T0/(kappa+delta_kappa). None of these three coordinates has been
assigned a value or a bound. A change of subtraction coordinates moves
terms between F_known and those constants, leaving the matched family
unchanged. The compact bound is deliberately on the stated known
representative, not on that entire family.
