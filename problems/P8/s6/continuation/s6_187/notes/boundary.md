# Exact initial-data compatibility

For each k and TT component write
L_k=d_t^2+3H d_t+k^2/a^2 and choose the actual fundamental
solutions f0(t0)=1,f0'(t0)=0 and f1(t0)=0,f1'(t0)=1.
Here t0=-1/2,T=1/2 and a0=25/16. Their Wronskian obeys

    W=f0 f1'-f0' f1=a0^3/a^3>0.

For compact q, let v=Gadv q solve L_k v=q with zero data
at T. The exact weighted Green identity is

    d_t[a^3(f v'-f' v)]=a^3(f L_k v-v L_k f).

Writing I_j=integral_(t0)^T a^3 f_j q gives

    v(t0)=I1/a0^3,       v'(t0)=-I0/a0^3.

Both moments must vanish. One alone does not remove both
initial data. The Fourier conditions use the actual
expanding-background homogeneous modes, not flat sine
and cosine substituted for them.

Since q is zero near both time endpoints, the two zero
initial data and zero final data imply that v vanishes
near both endpoints by uniqueness. For compact spatial
q, smooth finite-speed propagation on this finite slab
also makes v spatially compact. L and the TT constraints
commute, so v remains TT. Thus v is a legitimate smooth
compact stress test. Conversely q=L psi with compact
smooth TT psi has this property and Gadv q=psi.
This proves the equivalence and the construction in
projector.md proves nonemptiness.

For any distributional solution of the source-only tree
equation, weighted self-adjointness gives

    h[q]:=integral a^3 q_ij h_ij = T[v]/sqrt(kappa).

Every homogeneous tensor contribution pairs to zero on
this detector class. This is an algebraic restricted
observable identity; it does not construct a full quantum
metric state or assume that its initial data vanish.

For a generic compact q the advanced initial data need
not vanish. A sharp zero-extension at t0 can then fail
even continuity, and hence the proven H3-time stress-test
bound is inapplicable. This is a negative applicability
control, NOT a proof of divergence of this Proca/CD model.

Primary-source context: Hu, Roura and Verdaguer,
[gr-qc/0402029v1, Appendix D](https://arxiv.org/pdf/gr-qc/0402029),
analyze a finite-initial-time divergence in a specific
Minkowski massless conformal scalar example and discuss
initial correlations. Their example is not transferred
to Proca on CD. Their conservation concern also motivates
not inserting a finite gravitational coupling switch:
multiplying a conserved stress by an arbitrary time
function adds its derivative to the divergence.
Neither an uncorrelated zero quantum initial state nor
such a switch is used in this construction.
