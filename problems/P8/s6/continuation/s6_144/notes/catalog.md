# Exact mixed operator and all routes

The frozen master expansion contains -Tr[D W2 D F2]/2.
Equivalently vary -Tr[D(W2+h F2)D(W2+h F2)]/4 at h=0.
The code differentiates the actual two-site W2 and an independent
general quartic F potential four times, checking every independent
external-index pattern. Each of three channels has two assignments,
with coefficient one half each. F is the proper potential vertex;
the amplitude is minus W_amplitude times Gamma4_potential.

Fix one internal label A around the fermion box. Permuting B and
the two external labels gives six cyclic words with arc counts
(1,3),(2,2),(3,1), twice each. Both orientations remain. The W
vertex is C(z)+g(h1+h2), C(z)=-L+g/(M-z). Both internal heavy
propagators, including repeated or zero-transfer routes, remain.

The only divergent proper cycle is the four-fermion box, of degree
zero. A cycle through the W vertex uses both light lines and at least
one fermion line and has degree at most -1. Adding an internal-heavy
tree only improves it. The overall local divergence is paired with
the local-parent reference described below. There is no omitted
proper fermion mass, kinetic or Yukawa subgraph in this row.

Choose the light momenta q+K/2 and q-K/2 in each channel. At zero
external-soft momentum their common center q stays fixed. The four
fermion momenta have real bases k and l=k+q. Their shifts are linear
combinations of the four external vectors and K/2, with one-norm
below eighteen on the physical forward disc. A partial-sum routing
starts at A and follows each word; at zero soft momentum only the
A/B internal legs change the real base.

Separate the full proper kernel exactly as Gamma4_MS(q,p)
=Gamma0_MS(q^2)+DeltaGamma(q,p). The local MS box counterterm cancels
in DeltaGamma. Gamma0 has the two external-soft legs zero and internal
momenta q,-q. The split is not a derivative expansion of either light
propagator or W and does not set their physical external momenta to zero.
