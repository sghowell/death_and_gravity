# Full bare scalar map, not only a fixed-input coefficient

S6.155 gave the finite part of a transformed total vertex
coefficient at fixed MS inputs. To match bare parameters,
one must ALSO re-express its first pole as a function of
the FULL regulated hybrid coordinates q_H,D.

For coefficient X of Phi weight w, write its first total
pole as p1(q)/epsilon. Expanding
X_H+h p1(q_H,D)/epsilon+h^2 C2_H and equating it to
K_D^(-w)[X+h p1(q)/epsilon+h^2 C2_MS] yields

    X2_D = [w(w+1)k_D^2/2-w t_D]X
           +k_D(E_Phi p1-w p1)/epsilon
           +C2_MS-C2_H.

The last difference is pure poles. The required change
of the simple pole is k0(E_Phi p1-w p1), consistent with
the regulated vertex/line covariance. Its finite remnant is

    X2 = [w(w+1)k0^2/2-w t_MS]X
          +k1(E_Phi p1-w p1).

It would be incorrect to reuse the fixed-input term
-w k1 p1 while omitting the parameter derivative of p1.
The complete first MS residues are those checked in S6.155.

For G,L,M, the resulting map is

    G1=-k0 G,   L1=-2k0 L,   M1=0,
    G2=(k0^2-t_MS)G+k1 LG/Q,
    L2=(3k0^2-2t_MS)L+3k1 L^2/Q,
    M2=k1 g/Q.

The fermion Y^2 piece of p1L has the same field weight
as L, so it cancels from this commutator; it is not dropped.
The nonzero M2 comes from re-expressing p1M=g/(2Q), even
though the bare heavy coefficient has Phi weight zero.
M is a renormalized heavy mass squared, not a stable-heavy
on-shell mass.

Using fundamental G before squaring gives

    g2=g[3k0^2-2t_MS+2k1 L/Q].

For the three-channel heavy tree with b2=2g/(M-2)^3,
the full order-two field-direction coefficient relative
to tree is consequently

    3k0^2-2t_MS+(k1/Q)[2L-3g/(M-2)].

This includes the G1 square and the induced heavy-mass
variation. It is a field-direction allowance for the
future four-point assembly, not the whole two-loop
amplitude. The complete scalar map cannot be added a
second time to a differently normalized amplitude.

At the orders relevant to Phi two- and four-point functions,
first fermion proper counterterms remain the assigned MS
ones, with the full Phi-induced Y1_D=-k_D Y retained.
Possible external-fermion second-order field or mass maps
first affect these loop-induced functions at order three.
They are not claimed solved here.
