# Uniform full local bounds and the actual common-clock pullback

The physical scalar leg graph is

J=sum over f in(n_lapse,zeta,B) and j=0,1,2 of (1+|P|²)|partial_t^j f|,

with norm ||J||L2(dt dP), on compact prepared variations in I=[-1/2,1/2]. Fourier normalization is common to both legs. The exact a=(1+t²)² mixed density in [local.md](local.md) has73 polynomial monomials. Every monomial has exactly one field jet from each leg, time order at most two on each leg, and spatial degree0,2 or4. Split the degree-four factor as |P|² on each leg. This bounds each full monomial by its coefficient times J_D J_G, without a momentum cutoff or inverse spatial Laplacian.

For a rational time coefficient num(t)/den(t), factor it exactly. Each denominator here has strictly positive coefficients at even powers only. On I its denominator is at least den(0), and its numerator is bounded by the absolute coefficient sum weighted by (1/2)^degree. Expand coefficients in n and ell; use ell<462. Summing all actual terms gives the explicit rational C1 and C0 stored in the report, and

L_heat=(C1 n+C0)/(576 kappa0),

using 64 pi²>576. Cauchy-Schwarz on the two full leg graphs proves the integrated bilinear bound. This includes every shift, lapse, lower-curvature and second-metric term, not just a high-frequency symbol.

S240 proves the full state/subtraction reference stress bound
epsilon_s=10^310/(n kappa0)<10^-686. Its physical fixed-profile Hessian has absolute coefficient sum at most(2+6+9)a³ epsilon_s<68epsilon_s. Thus

L_physical=L_heat+68epsilon_s<10^-580.

No bound on the nonlocal state response is being inferred from the reference one-point stress.

## Nonlinear metric-chart contact

The same candidate has the exact physical-to-affine spatial metric map

zeta=v-(1/4)log R(t,N^-2).

Its first lapse jet is zeta=v+delta n_lapse, with delta=1/[2(1+t²)^3]. The second lapse derivative of zeta at N=1 is

4delta²-3delta,

which is not zero. The full R has R_X=2delta and R_XX=0 at X=1; higher R jets cannot be inserted or discarded arbitrarily. Exact delta jet bounds through order two are(1/2,3/2,33/4). They give J_physical<=12J_clock for the otherwise identical graph on fields(n_lapse,v,B).

The matched finite curvature action PLUS its own fixed reference profile has zero first Euler variation. Its second chart contact therefore vanishes after the compact first-variation identity is applied. Its complete bilinear bound becomes144L_heat. This argument is NOT applicable to the fixed STATE-PROFILE term alone.

Let rho_s,p_s denote the fixed normalized S240 state/subtraction reference integrals. Set

DeltaJ_s=(21delta²-3delta)(-p_s)/2
+(1-6delta)(-(rho_s+p_s)/2),
T_s=rho_s-3delta p_s.

The exact whole same-clock state-profile mixed density is

a³[2DeltaJ_s n_D n_G+3T_s(n_D v_G+n_G v_D)-9p_s v_D v_G].

Subtracting its physical Hessian under the first-order zeta map leaves the NONZERO contact

-3a³p_s(4delta²-3delta)n_D n_G.

The independent literal profile expansion verifies the entire second-variation chain rule. At t=0,a=1,p_s=1 and n_D=n_G=1 this contact is3/2, so setting it to zero fails an explicit control.

The bounds |DeltaJ_s|<4epsilon_s, |T_s|<=5epsilon_s/2, |p_s|<=epsilon_s and a³<= (5/4)^6<4 give the full state-profile mixed coefficient bound128epsilon_s. Finally,

L_same_clock=144L_heat+128epsilon_s<10^-580.

This is the complete local heavy finite heat plus full fixed heavy-profile Hessian in the ACTUAL common-clock variables, before quantum constraint elimination. It is not a same-space inverse, a full determinant response bound, or a nonlinear bounce estimate.
