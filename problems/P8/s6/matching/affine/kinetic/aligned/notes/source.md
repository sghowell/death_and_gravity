# Exact alignment and the nonlinear source that remains

All conventions, the unrestricted all64 stationary connection and the
source-centered mu=11/9 mass are inherited from the read-only S6.41
replay. The new operator is F(W), W=T-B(phi,x)dphi, not F(T).
The source coefficient is B=-15H(1+x)/(8h); B and its phi derivative
vanish on x=-1, and its lapse derivative is -d, d=15H/(4h).
Thus the complete first variation of B*dphi equals that of Tstar:
(-d*n,0,0,0). Both sources vanish on the actual background. This
uses the full source with its coefficient jets, not a fixed trace ansatz.

## Exact nonlinear identity

In the S6.41 spatial chart, h_physical=(4p²)^(-1/2)*h_hat and

    Tstar_normal=(4p²-1)*K_hat+6p*J3*s³,
    s=sqrt(-x), p²=(h-1+s²)/(4h), h=(1+u²)^3,
    J3=-(Q_lower+f_phi)/(4p*x),
    f_phi=h_phi*(1+x)/(2h²), h_phi/h=3H/2.

Substitution, without imposing the background, gives

    Sstar_normal=Tstar_normal-s*B
                =delta*(K_hat-3H*s)+(3/2)*s*Q_lower,
    delta=(s²-1)/h.

No lapse velocity is introduced. W0=T0-B(N,phi) is a triangular
point transformation with unit Jacobian in the timelike unitary chart;
the shared primary kinetic null and the prior ten-velocity rank persist.
This statement does not establish nonlinear secondary constraints.

The coefficient ODE has Q_lower=Q_lower,x=0 on the clock. Setting
s=1-eps*n, K_hat=3H+eps*delta_K and Q_lower=eps²*Q2 gives

    Sstar_normal=eps²*(-2*n*delta_K/h-6H*n²/h+(3/2)*Q2)+O(eps³).

This is generally nonzero. Writing the centered retained mass as
one half (W-Sstar)^T D^-1 (W-Sstar), its quadratic term is the free
clock Proca mass, but its cubic term includes -W^T eta*Sstar_second.
Coefficient and volume variations begin at cubic order as well.
There is no exact nonlinear spectator claim.

## Continuous source bounds

Use the original closed tube -11/10<=x<=-9/10, all real u. Here
h>=1, p²>=9/40 and |h_phi|/h<=3. The unchanged S6.37 ODE is

    Q_x+A*Q=F, Q(u,-1)=0,
    A=1/(2x)+3/(16h*p²),
    F=-3h_phi*(1+x)/(16h³*p²).

At fixed u, |A|<=25/18 and |F|<=(5/2)*|x+1|/h². Variation of
constants on a segment of length at most 1/10, in either direction,
therefore gives

    |Q| <= (5/4)*exp(5/36)*(x+1)²/h²
         <= (45/31)*(x+1)²/h² <= (3/2)*(x+1)²/h².

The exponential estimate uses exp(y)<=1/(1-y) for 0<=y<1,
termwise from their power series. Substitution back into the ODE gives

    |Q_x| <= (335/124)*|x+1|/h² <= 3*|x+1|/h².

Let eta=|x+1|/h<=1/10 and Delta_K=K_hat-3H*s. Since s<21/20,

    |Sstar_normal| <= eta*|Delta_K|+(5/2)*eta².

In unitary gauge Sstar has only a temporal component
Sstar_0=Sstar_normal/s. Spatial derivatives keep phi, h and H fixed.
The physical electric component is E_i=-s*D_i(Sstar_0), so exactly

    -E_i=delta*D_i(Delta_K)
          +[(s+1/s)*Delta_K/h-3s²*Q_x]*D_i(s).

There is no magnetic source. With Euclidean spatial norms in the
physical orthonormal frame, ||D Delta_K||<=epsilon_K and
||D s||/h<=epsilon_s, s+1/s<3 and 9s²<10 imply

    ||E|| <= eta*epsilon_K+(3|Delta_K|+10eta)*epsilon_s.

The normalized normal source has physical factor 1/tau, and E has
factor 1/tau². These explicit bounds do not control an inverse
operator, a forced solution, loops, or an omitted-action remainder.
