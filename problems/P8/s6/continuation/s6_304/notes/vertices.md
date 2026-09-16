# Direct canonical action expansion

All momenta are outgoing. Scalar legs obey p_i^2=mu, the emitted k obeys
k^2=0, and sum_i p_i+k=0. H is a lower-index graviton field matrix;
T(p,q)=p q^T+q p^T-eta(p.q+mu) has upper indices.
With kappa scaled out, V3=T:H. Its off-shell Ward contraction is
T eta(p+q)=(p^2-mu)q+(q^2-mu)p.

For two metric fields A,B, write tr(A)=tr(eta A). Their mixed determinant
coefficient at unit metric coupling is
d2=tr(A)tr(B)/4-tr(eta A eta B)/2.
The mixed inverse-metric-density coefficient is
K2=eta d2-[tr(A)B^up+tr(B)A^up]/2+(A eta B+B eta A)^up.
The scalar seagull stripped of i is
V4=-4[p_lower K2 q_lower+mu d2]/kappa.
An independent exact determinant/adjugate expansion of eta+2xA+2yB
checks this coefficient without reusing the seagull implementation.

Modulo the Einstein boundary term, the action density is
sqrt(-g) g^mn (Gamma^r_mn Gamma^s_rs-Gamma^r_ms Gamma^s_nr).
Let C(H,p)^a_mn=eta^ar(p_m H_rn+p_n H_rm-p_r H_mn)/2,
with the factor i from differentiation omitted.
For a nondifferentiated field A, the second connection is
C2(A;B)=-eta A C(B), and the first inverse-density variation is
M1(A)=tr(A)eta/2-A^up.
Define B(M,C,D)=M^mn(C^r_mn D^s_rs-C^r_ms D^s_nr).
The stripped trilinear Einstein vertex is -4/sqrt(kappa) times the sum
over all six permutations(A,B,C) of
B(M1(A),C(B),C(C))+B(eta,C2(A;B),C(C))+B(eta,C(B),C2(A;C)).
The minus sign is the two Fourier derivatives; the factor4 is
(kappa/2)*(2/sqrt(kappa))^3. No nonlinear gauge-fixing term is added.

For a partition ij|lm, let H_ij=P T(pi,pj)/(pi+pj)^2.
Trace reversal simplifies its numerator exactly to
pi_lower pj_lower+pj_lower pi_lower+mu eta.
The external-emission graph is
V3(epsilon;pi,-pi-k) V3(H_lm;pi+k,pj)/(2pi.k).
The two seagulls are -V4(epsilon,H_lm;pi,pj) and its exchanged partner.
The internal term is +Vhhh(H_ij,H_lm,epsilon), with vertex momenta
pi+pj,pl+pm,k. These signs follow the products of i vertices and propagators:
three vertices/two propagators give+i; two vertices/one propagator give-i.

The component implementation sums all21 graphs. Three exact rational recoil
states check shells, conservation and all four gauge directions separately
in every seven-graph channel. Their four-point amplitudes agree with the
unchanged original S297 gravity Born amplitude. These finite calibrations
supplement the general proof in the Ward note.
