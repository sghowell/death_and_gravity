# Exact actual-state error derivatives

For each fixed momentum, the finite-dimensional
Hamiltonian flow is smooth in the real parameter.
The pure-state graph remains well defined by its CCR.
Let G(r)=2i omega r+[R,r]+S-rSr, rhat be the finite
reference and e=r-rhat. Then

    e'=G(rhat+e)-G(rhat)-F.

The actual state mismatch e0 is not set to zero. Only
its parameter derivatives at the initial surface vanish:
all admitted histories and all their finite initial
jets coincide there, so e1(t0)=e2(t0)=0 exactly.

For a graph variation X,
A_r X=D_rG(r)X=2i omega X+[R,X]-X S r-r S X.
The fast phase and skew rotation are norm isometries.
S6.190 proves ||r||<1/10 and ||rhat||<1/100 in the
common high band. Thus the remaining tangent norm
growth is at most2*11/10<3 and the unit-slab comparison
factor is less than exp(3)<27.

Differentiating the exact error equation places A_r e1
and A_r e2 on their respective left evolutions. All
other terms are forcing; their complete noncommuting
expressions are checked independently. Write E_a for
the error norm and use R_a,S_a for coefficient bounds.

The first forcing norm is at most

    [2|omega1|+2R1
     +S1(||r||+||rhat||)+2S0||rhat1||]E0+||F1||.

For the second derivative, the forcing norm is at most

    [2S0||rhat2||+4S1||rhat1||+2|omega2|+2R2
     +S2(||r||+||rhat||)]E0
    +[2S0(||r1||+||rhat1||)+4|omega1|+4R1
      +4S1||r||]E1+||F2||.

The r1 term is bounded by ||rhat1||+E1. In particular
the quadratic E1^2 term is retained. No fast phase is
differentiated inside a proposed absolutely convergent
momentum integral.

Use |omega1|<=3nu_minus, |omega2|<=6nu_minus,
||rhat_a||<=A_a/nu_minus and nu_minus>=K. Writing the
constant parts relative to the indicated frequency
powers gives exactly the M1,M20,M21 formulas in
variation.py. With the frozen E0 coefficient1e30,

    E1 <=27[1e30 M1+C1/K] nu_minus^-9
       <2e32 nu_minus^-9,

    E2 <=27[1e30 M20/K+(2e32)M21
             +2S0(2e32)^2/K^10+C2/K^2]nu_minus^-8
       <1e35 nu_minus^-8.

The positive rational inequalities are checked before
any display rounding. One inverse-frequency power is
lost per parameter derivative, and the second derivative
still leaves sufficient decay for the infinite
three-dimensional energy-weighted tail.
