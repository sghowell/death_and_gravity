# All-N graph enumeration and analytic coefficient bound

Introduce formal markers c for one Phi^4 vertex and b for one H Phi^2
vertex. Give every other allowed vertex weight1. For labeled external
Phi leaves use x and for graviton leaves y. With a root species fixed,
the exponential generating functions obey
 P=x+P(exp(h)-1)+b P H exp(h)+c P^3 exp(h)/6,
 H=H(exp(h)-1)+b P^2 exp(h)/2,
 h=y+[(P^2+H^2)/2+b H P^2/2+c P^4/24]exp(h)
     +exp(h)-1-h.
The factorials come from unordered identical-species branches. The
full canonical vertex is already symmetric; no extra graph symmetry
factor is appended to the labeled root-cut evaluation.

There are no terms independent of x in P,H. Through the only powers
needed for four external Phi legs, write
 P=x P1+x^3 P3+O(x^5), H=x^2 H2+O(x^4),
 h=h0+x^2 h2+O(x^4).
With z=exp(h0), D=2-z, the x-independent equation is
 h0=y+exp(h0)-1-h0.
The next equations give
 P1=1/D, h2=z/(2D^3), H2=b z/(2D^3),
 P3=c z/(6D^4)+(1+b^2)z^2/(2D^5).
Therefore T_N=3!N![x^3 y^N]P is the stated closed EGF.
Differentiating the implicit equation gives
 d/dy=(z/(2-z))*d/dz,
so every finite N count can be computed without truncating higher
vertices. The only coupling sectors are c,b^2,1; no omitted higher
power can enter a connected tree with four external Phi leaves.

Independent dynamic enumeration cuts the fixed Phi-root vertex,
enumerates literal unordered partitions of labeled proper leaf sets,
and sums the allowed current-species assignments. It gives
 N0: c+3b^2+3,
 N1: 5c+21b^2+21,
 N2: 38c+198b^2+198,
 N3: 388c+2364b^2+2364,
 N4: 4972c+34236b^2+34236.
This enumeration uses neither the EGF solution nor its derivative.

## A genuine all-N graph majorant

For |y|<=1/8 and |h|<=1/4, the map
 h -> y+exp(h)-1-h
sends the disc strictly into itself:
 |map|<=1/8+exp(1/4)-1-1/4<5/24<1/4.
Its derivative norm is at most exp(1/4)-1<1/3.
The elementary exponential-series bound exp(1/4)<4/3 suffices.
Uniform contraction constructs a holomorphic fixed point on a
neighborhood of the closed y-disc. Then |z|<4/3 and |2-z|>2/3.

For c=b=1, the absolute EGF is bounded on that circle by
 (4/3)(3/2)^4+6(16/9)(3/2)^5=351/4<88.
Cauchy's coefficient estimate gives
 T_N<88*8^N*N! for every N>=0.
The nonnegative coefficients count topologies, so no sign cancellation
is used in this statement.

This analytic disc belongs to the graph-count generating function,
not to the momentum-space S-matrix or a forward/Regge amplitude.
A graph majorant does not control the internal propagators. Combining
loose factorial upper bounds need not prove summability, and failure
of that crude bound to sum would not prove that the physical series
diverges.
