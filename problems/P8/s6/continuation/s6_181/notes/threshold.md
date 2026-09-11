# Continuous massive-threshold estimates

All frequencies are first measured in units of the fixed positive
mass. Put Omega=2+x, z=1-4/(2+x)^2 and

    C(z)=3-2z+3z^2, A(z)=16/15+z-3z^2,
    D=A+sqrt(z) C atanh(sqrt(z)), U=sqrt(z) C/2,
    Q=D^2+pi^2 U^2, a(x)=U/Q.

These are literally the S6.86 density used by S6.180; the code
compares the source expressions. No new form factor is chosen.
Since atanh(sqrt(z))>=sqrt(z),

    D>=16/15+4z-5z^2+3z^3
      =16/15+z[3(z-5/6)^2+23/12]>=16/15.

Also C=8/3+3(z-1/3)^2>0. An elementary rational pi bound suffices:
polynomial division and integration give

    integral_0^1 t^4(1-t)^4/(1+t^2) dt=22/7-pi>0.

Thus pi^2<(22/7)^2<10. This is checked symbolically, not inferred
from a floating-point pi value.

For 0<x<=2, 0<z<=3/4, z'=8/(x+2)^3<=1 and
|z''|=24/(x+2)^4<=2. Define the analytic endpoint function

    f(z)=atanh(sqrt(z))/sqrt(z)
        =integral_0^1 (1-z t^2)^-1 dt.

Differentiation under this uniformly bounded integral gives
f<=4, f'<=16 and f''<=128 throughout [0,3/4]. Triangle bounds give

    |C|<=4, |C_z|<=8, |C_zz|=6,
    |A|<=4, |A_z|<=6, |A_zz|=6,
    |zC|<=3, |(zC)_z|<=10, |(zC)_zz|<=21.

With D=A+(zC)f these imply

    D<=16, |D_z|<=94, |D_zz|<=794,
    |D_x|<=94, |D_xx|<=982.

For V=U^2=z C^2/4, direct product differentiation gives
|V_z|<=16, |V_zz|<=65, |V_x|<=16 and |V_xx|<=97.
Consequently Q>1 and

    |Q_x|<=3168, |Q_xx|<=50066,
    Q^-1<=1, |(Q^-1)_x|<=3168,
    |(Q^-1)_xx|<=2*3168^2+50066=20122514.

The nonsmooth threshold factor is isolated, not differentiated as
though it were smooth: a(x)=sqrt(x)b(x), with

    b(x)=c(x) C(z)/(2Q), c(x)=sqrt(x+4)/(x+2).

Here c<=1 since (x+2)^2-(x+4)=x(x+3)>=0.
Using x+2>=2, x+4>=4 and sqrt(x+4)<=sqrt(6)<3,
the literal first and second derivatives give

    |c'|<=1/8+3/4<1,
    |c''|<=1/64+1/8+3/4<1.

The x derivatives of C obey |C_x|<=8, |C_xx|<=22.
Therefore, for d=cC/2, |d|<=2, |d'|<=6, |d''|<=21.
The product d Q^-1 yields

    |b|<=2, |b'|<=6342, |b''|<=40283065.

Since x<=2, the three scaled derivatives of a satisfy

    |a|/sqrt(x)<=2,
    |x a'|/sqrt(x)<=12685,
    |x^2 a''|/sqrt(x)<=161144944.5.

The common outward bound C_small=200000000 is valid on the
entire threshold interval. All constants and product/chain
identities are exact rationals in threshold.py.
