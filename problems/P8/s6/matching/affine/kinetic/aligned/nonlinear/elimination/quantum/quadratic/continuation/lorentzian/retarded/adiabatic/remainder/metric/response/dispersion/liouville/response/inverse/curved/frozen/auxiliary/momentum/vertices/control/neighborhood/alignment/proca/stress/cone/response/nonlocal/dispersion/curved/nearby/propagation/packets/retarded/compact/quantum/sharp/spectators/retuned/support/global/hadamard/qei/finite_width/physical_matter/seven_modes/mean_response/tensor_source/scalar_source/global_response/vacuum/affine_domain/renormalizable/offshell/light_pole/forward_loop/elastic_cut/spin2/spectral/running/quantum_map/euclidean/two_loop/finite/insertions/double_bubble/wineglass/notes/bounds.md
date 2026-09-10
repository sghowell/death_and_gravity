# All-radius integrated bounds

First consider the subtracted local core. After
the common outer shift, split its integrand as

    [I_R(Q)-I_R(y)]/(y+Delta)^2
    +I_R(y)[(y+Delta)^-2-(y+1)^-2].

The radial measure is y dy/(16pi^2).
The inner bubble supplies the other loop measure.

For the first term, D(y) from notes/kernel.md gives

    integral_0^infinity y D(y)/(y+1/4)^2 dy
      <= integral_0^infinity
             [8y^(3/2)+4y]/(y+1/4)^3 dy
       = 6pi+8 < 32.

For the second term, |1-Delta|<=3/4 and
|1+Delta|<=11/4 give

    |(y+Delta)^-2-(y+1)^-2|
      <= 33/[16 (y+1/4)^2(y+1)].

Also ln(1+y)<=sqrt(y). With t=sqrt(y), the
difference t-ln(1+t^2) vanishes at zero and has
derivative (t-1)^2/(1+t^2)>=0.
Thus its radial integral is bounded by

    (33/16) integral_0^infinity
                 y^(3/2)/(y+1/4)^3 dy
      =99pi/64 <99/16.

Both elementary integrals are checked using
anchored primitives. For the half power set y=t^2;
a primitive is

    (3/2)atan(2t)-t(20t^2+3)/(4t^2+1)^2.

Its endpoints are zero and 3pi/4. For the integer
power a primitive is
-1/(y+delta)+delta/[2(y+delta)^2], giving
1/(2delta)=2 at delta=1/4.

Since 32+99/16<40, the fully subtracted core obeys

    |J_R(P)-J0| < 40/(16pi^2)^2.

This is an integrated convergent difference, not
an estimate of either divergent term separately.

For the remaining heavy-weighted numerator, keep
the complete propagators. The inherited routing gives

    |T_A(q)| <= 4g/(y+M),  |V_A|<3L,
    |1/(M+Q^2)| <= 2/(y+M).

Using g/M<L/3 then gives

    |C(Q)^2-L^2| < (16/3)Lg/(y+M),
    |L^2T_A+V_A(C(Q)^2-L^2)|
       <20L^2g/(y+M).

Its inner logarithm is bounded by ln(1+y)+3.
Put ell=ln(4M). The inequality

    ln(1+y) <= ln(1+M)+ln(1+y/M)

follows from
(1+M)(1+y/M)-(1+y)=M+y/M>0.
Use ln(1+M)<ell, the exact triangle integral
of S6.123 and the exact logarithmic integral
of S6.122. They imply

    integral_0^infinity
      y[ln(1+y)+3]/[(y+1/4)^2(y+M)] dy
         <2ell(ell+4)/M.

In the logarithmic remainder term, y/(y+M)<=1
reduces the integral directly to the S6.122 one.
No finite upper integration limit or physical
momentum cutoff has been introduced.

Finally |C(P)|<2L. Summing the original total
external-assignment weight three, the local
subtracted and heavy terms give

    sup_disc |A_wineglass|
      <[240+40ell(ell+4)] L^3/(16pi^2)^2.

Cauchy's radius-one estimate bounds the second
forward coefficient by this same whole-amplitude
bound, with no extra factor 1/2.

At the actual parameters ell<600 and pi>3:

    E_wineglass = 14496240 L^3/20736
                  <10^-611,
    E_wineglass/(4lambda) <3 times 10^-12.

The diagnostic value is about 8.38232 times
10^-612. Exact rational comparisons, not that
decimal, prove the inequalities.

The four disjoint raw-graph groups still sum to
an absolute upper bound below 3 times 10^-607
and below 10^-7 of tree. This includes all
192 raw refinements with the stated subtractions,
but not the additional finite-potential insertions
or complete two-loop normalization.
