# S6.278: original finite-gravity M1 cut and conditional massive soft dictionary

Keep the whole current QG2-H8A420 parent, all fixed profiles, all three
vacuum constants, the physical matter frame and kappa=10^800. At its
retained quadratic flat reference the light Phi has tree mass squared1,
and Psi=sqrt(kappa)chi is the original minimally coupled massless M1.
No claim that these tree parameters equal exact quantum poles is made.

## Complete specified first gravitational contribution

For Phi Phi -> Psi Psi the Mandelstam sum is2m^2, whereas the elastic
four-Phi sum is4m^2. The exact finite-kappa production tree is

    A_pair(s,x)=[s-(s-4m^2)x^2]/(4kappa).

The full identical-intermediate-state optical relation gives

    Im A_M1,s(s,t)=[s^2-tu+2m^2s+6m^4]/(960pi kappa^2).

On physical s>=4m^2 this is derived by angular sewing, not by treating
M1 as a massive Phi. The one-loop massless pair branch starts at s=0;
continuation below4m^2 is not a physical external-state positivity claim.
The complete crossing nonlocal contribution is

    A_log=-c sum_z P_z log(-z/nu^2), c=1/(960pi^2 kappa^2),
    P_s=s^2-tu+2m^2s+6m^4, with s,t,u permuted.

A real local crossing polynomial remains unspecified by this cut. The
statement exhausts this first M1-loop nonanalytic sector, not graviton
loops, other species, all orders in short couplings or the full amplitude.
At t0,s=2m^2+v,u=2m^2-v, P_t=v^2+2m^4 is nonzero. Neither the graviton
pole subtraction nor the one-loop universal graviton soft factor removes
this distinct M1 logarithm.

## Massive analytic soft prescription

Define F(z) on the plane cut from4m^2 to infinity by

    F(z)=[(z-2m^2)^2-2m^4]
          integral_0^1 dx/[4m^2-z(1-x^2)].

Its principal-sheet closed form is recorded in soft.py; the apparent
pseudothreshold at zero is removable. Set B=2(F(s)+F(t)+F(u))-m^2.
Under the stated full-amplitude universal-factorization premise, analytic
soft stripping has the exact resolution relation

    A_E1=exp[-B log(E1/E0)/(4pi^2 kappa)] A_E0.

This relates the same prescribed stripped amplitude at two resolutions.
It does not construct the amplitude, make finite-resolution hard and
stripped amplitudes identical, supply a uniform complex-energy bound, or
prove finite-coupling unitarity. The new external detector framework is a
research route with explicit hypotheses, not a numerical bound imported
from its massless shift-symmetric scalar example.

## Exact low-cut convention and retained transfer term

For a squared subtraction cap L>4m^2, define

    A_low=-c sum_z P_z[log(-z/nu^2)-log((L-z)/nu^2)].
    A_log-A_low=-c sum_z P_z log((L-z)/nu^2).

The branch prescription and equality retain all three channels. A_low
has precisely the original logarithmic discontinuity below L and zero
above L. The remainder of this sector is analytic near the massive
crossing center. The original real local polynomial is not changed.

With D=L-2m^2, its coefficient of v^2 at t0 is

    b20_log,sub=-c[2log(D/nu^2)+log(L/nu^2)-12m^2/D-14m^4/D^2].
    d_L b20_log,sub=-c[2/D+1/L+12m^2/D^2+28m^4/D^3].

The1/L term comes from the transfer channel; omitting it is an error.
At m^2=nu^2=1,L=10^196, the nonlocal absolute coefficient is below
2000/(960kappa^2), less than10^-990 of the original4lambda.
This bound does not constrain the independent local coefficient or the
full quantum error, and L is not a certified Wilsonian cutoff.

For illustration only, a q-weighted graviton-pole coefficient integrates
to log(Q/E)/kappa. Its fixed-resolution decoupling and joint detector
scaling limits differ. The weight is not claimed to be a positive UV
functional. Finite counterterms, all other cuts, detector errors, Regge
control, common-parent matching and original P8 closure remain open.
