# Ordered curvature-adapted quantum reference inverse

Let

Fdiag=diag(Ftrace(D_eta^2+q),(8/3)F2(D_eta^2+q)),
Kdiag=Fdiag^-1,
Jdiag=I Kdiag.

These are exactly the fixed flat scalar factors of S224, not an asserted time-dependent mass or state substitution. The shear factor retains its full subthreshold pole and continuum. Their shifted positive static spectral measures give the uniform primitive bound

||Jdiag(t;q)||<=45/2,

for allq>=0 andt>=0. A uniform half-line L1 estimate on Kdiag is neither needed nor claimed.

For the variable-coefficient curvature reference Aref=Bc* Fdiag Bc, the inverse MUST be in the order

Ecurv=Y Kdiag Z, Y=Bc^-1, Z=(Bc*)^-1.

At fixed source times, Z(s,s)=0 allows the middle causal derivative Kdiag=D Jdiag to be moved onto the first time argument ofZ, with no omitted initial term. Therefore Ecurv is the ordinary two-time kernel

Ecurv(t,s)=integral_s^t du integral_s^u dv
          Y(t,u) Jdiag(u-v) partial_v Z(v,s).

Both Y andZ are genuinely two-time kernels. They are not translation-invariant convolutions, and no coefficient is commuted throughJdiag. The exponent in the product is bounded by

20(t-u)+108(v-s)<=108(t-s).

The remaining simplex integral of(t-u) equals(t-s)^3/6. With C0=5/2 and Jmax=45/2,

||Ecurv(t,s)||<=C0^2 Jmax (t-s)^3 exp(108(t-s))/6
             =375(t-s)^3 exp(108(t-s))/16.

The order itself proves both identities as causal distribution compositions:

Aref Ecurv=Bc* Fdiag (Bc Y) Kdiag Z=Bc* Z=I,
Ecurv Aref=Y Kdiag (Z Bc*) Fdiag Bc=Y Bc=I.

The cancellations are of adjacent matching factors; no commuting approximation is used. The elementary factor inverse identities and association hold on compact causal time triangles, and the Sobolev-valued extension is described in notes/spaces.md. The full initial boundary is part of those identities.

The bounds are an existence estimate for this specified reference. They do not show small backreaction, stability or absence of poles in a different full operator. In particular the huge physical normalization is restored, not silently canceled.
