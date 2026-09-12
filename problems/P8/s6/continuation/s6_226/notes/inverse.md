# Ordered bounded-channel inverse and its complete causal graph

Let K=Kdiag be the ORIGINAL bounded weighted causal inverse from notes/weight.md. Multiplication by the measurable time-dependent2x2 matrixV(t), with essential-sup norm<=M, has norm at mostM on every H_sigma^r. No commutation withK is assumed.

The chosen weight gives

||KV||,||VK||<=5M/(32+13M)<5/13<1/2.

The Neumann series converges in bounded-operator norm on the weighted space, and every partial product is causal. Its limit is causal because causal support is closed. Define

K_V=(I+KV)^-1 K=K(I+VK)^-1.

The equality follows term by term, (KV)^n K=K(VK)^n; it is not an assertion thatKV=VK. Its norm is at most

[5/(32+13M)]/[1-5M/(32+13M)]
=5/(32+8M).

## Both initial-boundary graph identities

The original forward factors have the explicit paired local/cut distribution definition inherited from S224. They send weighted L2 functions locally to finite-order time distributions with a finite spatial Sobolev loss. Their original left and right causal inverse identities include the complete initial boundary.

Define the graph domain of Fdiag+V as the causal u inH_sigma^r for which (Fdiag+V)u is an ordinary source inH_sigma^r as a distribution, including possible initial atoms. For f inH_sigma^r, let u=K_V f. The equation

u+KVu=Kf

holds inH_sigma^r. Applying the original forward distribution Fdiag gives (Fdiag+V)u=f. Conversely, if u is in the graph domain, apply the original inverseK to its distributional equation to obtain the same bounded equation. Invertibility ofI+KV gives u=K_V f. This proves existence, uniqueness and BOTH inverse identities on the full graph.

One can justify all distributional applications first for smooth compact-momentum inputs, use the locally continuous finite-order forward map, and then pass to the weighted L2 limit. No initial delta is projected away and no forward bounded self-map is assumed.

For a requested finite interval[0,T], extendV by zero outside it and extend a source by zero to the right. Causality makes the response on[0,T] independent of those future extensions. The graph identities are local at the complete initial boundary and throughout the interval. This supplies the same finite-window bound without inventing a final-time boundary condition.

The proof actually needs a bounded causal middle operator on this chosen weighted space; the stated certified class here is the explicitly norm-bounded time-dependent matrix multiplication. An application to another actual remainder must establish that operator bound and its graph compatibility, rather than declare it from differential order or formal matching.
