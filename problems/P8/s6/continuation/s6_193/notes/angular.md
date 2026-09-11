# Full-shear angular contractions in general spatial dimension

In d spatial dimensions the actual constrained Proca Hamiltonian is

    K=a^(2-d)E+kk^t/(a^d m^2),
    V=a^(d-2)m^2Q+a^(d-4)[(k^tQk)Q-Qkk^tQ], Q=E^-1.

Multiplication gives KV=omega^2 I with
omega^2=m^2+a^-2 k^tQk and d physical polarizations. The rank-one
temporal-constraint contribution is retained.

At one evaluation time a CONSTANT unimodular spatial coordinate
change sets E=I. It acts on every time jet by that same constant map;
it is not a time-dependent change of momentum variables. Put

    F=E', G=E'', trF=0, trG=trF^2,
    Pn=nn^t, z=|k|^2/(a^2 omega^2), n=k/|k|, u=2-d.

The local symbols then have
A=F-zFPn+uH I-2HzPn and

    K''K^-1=G-zGPn+2uH(F-zFPn)
      +(uH'+u^2H^2)I+z[-2H'+4(d-1)H^2]Pn.

With f=nFn, f2=nF^2n and g=nGn,

    p=-z(H+f/2),
    p'=z(2H^2-H'+2Hf+f2-g/2)-2z^2(H+f/2)^2.

These full expressions are used in tr(s^2), tr(t^2+s^4).

A trace word containing Pn splits into products of quadratic
projections n^t W_j n, preserving the order of F/G inside each W_j.
For r projections the uniform unit-sphere average is the sum over
all pairings of their2r endpoints, divided by
d(d+2)...(d+2r-2). This follows by taking Cartesian Gaussian moments
and dividing out the independent radial moment. Each contraction
loop is the trace of its oriented matrix-edge words. A reversed
traversal reverses a word. Only cyclic trace rotations, reversal
of symmetric-factor words, and Pn^2=Pn are used.

In particular tr(FGFG) and tr(FFGG) are distinct. The code retains
both in a noncommuting control and enumerates all105 fourth-pair
contractions. Independent expanded Cartesian polynomials integrated
by their exact sphere monomials test these formulas in dimensions3,4,5.

For radial integrals the ratio of the z^n term to the n0 term is

    (d/2)_n/(alpha)_n for a denominator omega^(2alpha).

The second-order case has alpha=1/2; fourth order has alpha=3/2.
This is the meromorphic gamma-integral identity. Do not set d=3
inside these ratios before extracting the finite dimension limit.
