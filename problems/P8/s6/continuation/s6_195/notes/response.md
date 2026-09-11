# Actual off-diagonal covariance and both response terms

Use the Fourier covariance C(k,q)=<sym Z(k)Z(q)^dagger>.
The unperturbed actual common state is diagonal in momentum:
Cbar(k,q)=delta(k-q) C_k. Its full physical six-dimensional
covariance includes all three polarizations and both canonical
quadratures. The Fourier formulas are the complexification
of the real-field canonical algebra, with reverse-pair reality
conditions; they do not introduce independent complex fields.

For a prescribed external perturbation with unchanged initial
state, A_k=J M0(k) and

    delta C_kq'=A_k delta C_kq+delta C_kq A_q^t
       +J M_Gamma(k,q) C_q+C_k M_Gamma(k,q) J^t,
    delta C_kq(t0)=0.

The last equality follows from the common zero initial
neighborhood, not from an assumed zero quantum metric state.
Duhamel integration uses U_k(t,s) on the left and U_q(t,s)^t
on the right of the entire source. Replacing the second
momentum or propagator by the first is not permitted.

The current variation has two parts:

    delta J_D(t)=-1/2 integral dk dq
                    tr[M_D(q,k;t) delta C_kq(t)]
                -1/2 integral dk
                    tr[M_DGamma(k,k;t) C_k(t)].

All momentum integrals use the Fourier measures fixed in
notes/hamiltonian.md; the source/detector Fourier pairing
reduces the propagation term to q=k-p at a fixed external p.
The contact is local in time, but contains the complete
spatial Fourier convolution.

On a common finite real-mode regulator, set H_M=Z^t M Z/2.
The exact CCR identity
[H_M,H_N]=(i/2)Z^t(MJN-NJM)Z, symplectic transport and trace
cyclicity reproduce the covariance tangent as

    delta J_D=-<H_DGamma>
               +i integral_(t0)^t <[H_D(t),H_Gamma(s)]> ds.

Both reverse-pair blocks are counted. The sign comes from
perturbing the Hamiltonian by +H_Gamma and reading out -H_D.
The metric contact is not the stress noise anticommutator.

The separate commutator and second-metric-variation kernels in
[Hu and Verdaguer,0802.0658v1,eq.(4.6)](https://arxiv.org/pdf/0802.0658)
provide primary-source context for this separation; their
scalar calculation is not substituted for the Proca algebra.
No renormalized infinite spatial limit is claimed by this note.
