# Uniform complex discs and the same C10 inverse

Let r0=10^-8. For each real t in I, every point z in the
closed complex disc |z-t|<=r0 obeys |z|<51/100,
7/10<|1+z^2|<13/10. For the exact positive polynomial P,

    sum_l l|P_l|=1363114784,
    |P(z)-P(t)|<=r0 sum_l l|P_l|=13.63114784<215.

The estimate follows by integrating P' along the straight
segment, on which |z|<1. Since P(t)>=1215, |P(z)|>1000
on every disc. There is a uniform analytic extension and

    |1/(2J(z))|<400(13/10)^18/1000<50.

From the actual rational formulas, these same discs give
|H|<3, |delta|<3/2, |ell|<1, |w0|<6, |L|<18 and
|theta|<6. For theta use H-z/(1+z^2)^4. The weighted
phase matrix has row sums below200000; the weighted
physical-force matrix has row sums below30000.
The lapse readout n=c y+d g, with y=(vhat,10p), has
complex row norms |c|<1000 and |d|<=125.
All these bounds use explicit product/triangle inequalities,
not sampling on the discs.

Cauchy's integral formula, applied entrywise and summed
under the integral for each row, bounds every jth
derivative by its matrix/row constant times j! r0^-j.
The same argument gives
|delta^(j)|< (3/2)j!r0^-j.

For a smooth prepared force let G=||g||C10 and
Y_j=||y^(j)||joint. Set x_j=r0^j Y_j/j!.
Differentiating y'=A y+B g exactly j times gives

    x_(j+1) <= r0/(j+1) [
        200000 sum_(i=0)^j x_i
        +30000 G sum_(i=0)^j r0^i/i! ].

No higher derivative of g than j is used. The real bound
gives x0<15000G. Since r0<1 and
r0(200000*15000+30000)<15000, induction proves
x_j<=15000G for j<=10. This does not presume analytic
forces; only the coefficients use complex continuation.

Leibniz's rule in the lapse reconstruction gives

    ||n^(j)|| < (j+1)j!r0^-j(15000000+125)G.

The physical scale v=vhat+delta n satisfies

    ||v^(j)|| < j!r0^-j [
        15000+(3/4)(15000125)(j+1)(j+2)]G.

The triangular factor is the exact sum
(3/2)sum_(k=0)^j(k+1).
At j=10 these bounds are respectively
5.987569896e94G and5.3888673384e95G.
Every j=0,...,10 entry is checked exactly and is below1e96G.
Consequently ||Bcl g||C10<1e96||g||C10.
This intentionally loose bound is a genuine continuous
no-loss CLASSICAL estimate; it is not a quantum C10 norm.
