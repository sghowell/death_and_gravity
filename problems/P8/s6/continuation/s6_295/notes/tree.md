# Complete tensor, Ward identity and improved heavy remainder

Use all-outgoing k_i^2=mu,q^2=0,sum k_i+q=0, eta=+---.
T_m(p,r)=p r+r p-eta(p.r-m). Put d_i=2k_i.q and
T_i=T_mu(k_i,k_i+q). For the three partitions L={1,j},R=complement,
P_L=sum_L k_i,P_R=-P_L-q,D_L=P_L^2-n,D_R=P_R^2-n,
S_L=sum_L T_i/d_i, and T_H=T_n(P_L,P_L+q).

sqrt(kappa) M=C(sum_i T_i/d_i+eta)
-g^2 sum_L[S_L/D_R+S_R/D_L+T_H/(D_LD_R)+eta(1/D_L+1/D_R)].
This is the literal26-graph amplitude, including every metric contact.
T_i.q=d_i k_i and T_H.q=P_LD_R+P_RD_L. Thus each heavy group
and the entire quartic group have zero Ward contraction separately.

For a physical unit-Frobenius helicity tensor, q.epsilon=0 and trace0.
Define J_i=k_i k_i/(k_i.q), Abar=C-g^2/2 sum_L(1/D_L+1/D_R).
Then M_TT=Abar sum_i J_i/sqrt(kappa)+R_TT, where
R=g^2/sqrt(kappa) sum_L[(P_L.q)(J_L-J_R)-2P_LP_L]/(D_LD_R).

The sum of all three numerators is2(k_1 q+q k_1), pure gauge.
Proof: the coefficient of each J_i in sum_L(P_L.q)(J_L-J_R)
is2k_i.q, while sum_L P_LP_L=sum_i k_i k_i-(k_1 q+q k_1).
Consequently each denominator in R_TT can be replaced by
1/(D_LD_R)-1/n^2. This improvement is essential; it is not dropping
a heavy matching term or altering the exact tree.

At original mu1, the six invariants a_ij=(k_i+k_j)^2 sum to8.
Expanding each1/(n-a) about a2 gives the exact positive identity
Abar=g^2/[2(n-2)^2] sum_six(a_ij-2)^2/(n-a_ij).
For the physical four-point Born limit, the timelike invariant is>=4,
giving A0>=4g^2/[n(n-2)^2]>4g^2/n^3 below the heavy pole.
This is tree-level positivity on a real subheavy domain, not a quantum
spectral or forward dispersion result.

Literal D4/5/6 matrices, two rational physical configurations, two
transverse rotations, every group/deletion Ward, the general coefficient
and outer-product proof, and public API checks are retained. Each D4
helicity has norm1; the two-helicity sum equals the full TT projector.
