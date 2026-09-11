# Complete nearby-clock source: a spatial Frechet bound

Write u=t+psi on the fixed CD geometry. Use the full
source S_mu=partial_mu u (R-1)Q from S6.178, with
Q=E/X+3R_u/(4R)+C_Z Z,
E=Boxu-3H(u)X,
C_Z=-1/X^2+3R_X/(2RX).
The actual new scalar retuning does not change this
source. Its first scalar-jet derivative vanishes at
the zero perturbation jet.

Take all15 independent psi/value/gradient/symmetric
Hessian coordinates in the complex ball with each
component of modulus<=1/50. Then |u|<=13/25,
|X-1|<=2/50+4/2500<1/20 and
|1+u^2|>=1-(13/25)^2.

To bound R_X, use the larger complex X disc |X-1|<=1/10.
There |(1-X)/X|<=1/9. The full switch obeys
|T-1|<=2(1/9)^1024, |T|<2.
For w=N X^2 exp(-N X^2), Re(X^2)>=4/5;
N|X|^2 exp(-4N/5)<1 at N=1024.
Thus |B|<3. Cauchy on radius1/20 then gives
|B_X|<=60 in the inner X disc.
No full switch or exponential term is truncated.

At |u|<=13/25 and |X-1|<=1/20 these imply
|R-1|<2/5, |R|>1/2, |R_X|<16 and
|R_u|<2, using R_u=-6u(R-1)/(1+u^2).
Also |H(u)|<3 and |C_Z|<52.

The exact covariant CD jet formulas, including every
connection term, give |Boxu|<5 and |Z|<3/100.
Hence |E|<15, |Q|<21. Since every component of du
is at most51/50, every full source component has
modulus below10 throughout the complex jet ball.

About each real psi jet with components<=delta<=1/100,
a polydisc of radius1/100 fits in this complex ball.
Cauchy bounds every second jet derivative of each
source component by2*10/(1/100)^2=200000.
The same bound holds along the real segment to the
zero perturbation jet. Since DS[0]=0, each first
jet derivative at the actual history is bounded by
15*200000*delta=3000000delta.

For a real smooth compact test eta,
partial_i DS_eta has a first-derivative term with15
inputs and a second-derivative term with15^2 inputs.
Every spatial derivative of the background psi jet
is a third coordinate jet bounded by delta.
Each source component is therefore bounded pointwise
by [15*3000000+15^2*200000]delta G_eta.
Four source components and three spatial components
give the safe joint bound

    ||grad_spatial DS_eta||L2
        <=400000000 delta U_eta,

where G_eta^2 sums all coordinate jets through3 and
U_eta(t)=||G_eta||L2(space).
The sharper already source-pinned undifferentiated
bound is ||DS_eta||L2<=2048delta U_eta.
These are full-function bounds. Only psi jets, not
the arbitrary test eta, need to be small.
