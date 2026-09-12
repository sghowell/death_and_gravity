# Primal variable-coefficient coordinate inverse

Write Bc=B0+C D in conformal time. Let R0(t;q) be the S224 causal inverse kernel of B0. Its exact formula is

R0=[[sq,-sq/3],[sq-t,-sq/3-2t/3]],
sq=sin(sqrt(q)t)/sqrt(q), with s0=t.

It satisfies R0(0)=0, ||R0'||<=C0=5/2, and ||R0||<=C0 t, uniformly for everyq>=0. The exact Frobenius square of R0' is(20cos^2-14cos+13)/9<=47/9<25/4. These bounds avoid dividing byq, including at the continuous endpoint.

For a causal source f, the equation Bc y=f with zero-past coordinate data is equivalent to

y=R0*(f-Cy').

Set u=y'. Since R0(0)=0, differentiation gives the ordinary Volterra velocity equation

u=R0'*f-R0'*(C u), y=Iu.

Its kernel has norm at mostC0*8=20. A Neumann-Volterra iteration on any finite slab converges by the simplex factorial bound: the nth iterated kernel is at most20^n elapsed^(n-1)/(n-1)! for n>=1. This argument requires neither a small window nor a small norm20T.

The impulse-response velocity kernel V(t,s) therefore satisfies

||V(t,s)||<=C0 exp(20(t-s)).
Y(t,s)=integral_s^t V(v,s)dv,
||Y(t,s)||<=C0(t-s)exp(20(t-s)).

V=partial_t Y andY(s,s)=0. Differentiating the defining Volterra equation gives Bc Y=identity in the causal distribution sense, including the correct leading L^-1 initial derivative. Conversely, any causal homogeneous solution solves the zero-drive velocity equation and hence vanishes. This proves both Bc Y=I and Y Bc=I on the complete zero-past graph, not merely a one-sided formal symbol product.

For each fixedq the smooth-coefficient ordinary differential operator has the same causal inverse on distributions: distributional initial atoms are carried by Y, not deleted. On a compact time slab, higher time derivatives of the kernels follow from Bc=L D^2+Q+C D and the finite smooth coefficient jets. The proven zeroth/first derivative bounds are uniform inq; higher finite derivative bounds grow only polynomially inq. This supports the Sobolev-valued distribution extension used in notes/spaces.md. It is not an exponential high-momentum estimate or a claim of a full quantum-force inverse.
