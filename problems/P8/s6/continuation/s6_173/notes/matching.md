# Constraint and domain data that the parent must match

At the unit clock, h=(1+t^2)^3>=1, the actual coefficients are

    B=F2=-1/2, F2X=-1/(2h),
    A3=1/h, A4=-1/h-7/(4h^2), A5=1/h^2.
    C=4F2X+A3=-1/h,
    Dkin=A3+A4+A5=-3/(4h^2).

The homogeneous trace/lapse velocity block is

    (2B/3)K^2+C K V+Dkin V^2
       =-[K/sqrt(3)+sqrt(3)V/(2h)]^2.

Its two-variable Hessian has rank one. Deleting A4,A5 leaves
determinant -(4h+3)/(3h^2)<0 and rank two. Their vanishing first
background Euler sources do not authorize their removal. This
rank control is not itself a fresh proof of every propagating
mode's stability or a classification of the altered theory.

The original M1 scalar must be varied before imposing its
conserved-momentum solution. For L=a^3 chidot^2/(2N), canonical
momentum p=a^3 chidot/N. Direct substitution chidot=Np/a^3 gives
L_sub=Np^2/(2a^3); varying it as if it were the original action
incorrectly flips both energy and pressure. The fixed-momentum
Routhian L_sub-p*chidot=-Np^2/(2a^3) gives the correct signs.
The calculation uses the original field variation, not the
incorrectly substituted action.

The S6.164 action comparison remains on the same real Schwartz
fields with Fourier support in the Euclidean unit ball and Fourier
L1 norm<=1. All coordinate first jets have sup norm<=1; thus
|X|<=4/kappa. The canonical clock instead has

    partial_t Psi=sqrt(kappa)=1e400, X=1.

For every member of that fixed class, even on a finite interval,

    |partial_t(Psi_clock-Psi)|/sqrt(kappa)
        >=1-1/sqrt(kappa)>1/2.

Scaling the field to reach the clock changes the fixed Fourier-L1
class and its nonlinear estimates. The old tiny common-class
action norm does not supply a clock-neighborhood action, equation
or variation bound. It remains valid in its own declared domain.

Together, the leading source table, exact first jets, constraint
block and domain separation specify necessary data for a common
propagating parent. They do not construct that parent or bound its
additional fields, quantum state, cutoff or omitted operators.
