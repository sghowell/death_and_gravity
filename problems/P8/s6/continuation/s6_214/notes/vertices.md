# Complete first and second spatial/lapse vertices

Use the same real-field Fourier pairing, q=k-P, and the unchanged ten-row energy feature map

    F_k Z=(m sqrt(a)A, C_k A/sqrt(a),Pi/sqrt(a),
                  k.Pi/(a^(3/2)m)),

where C_k is the real cross-product matrix. Its omitted explicit Fourier i factors are retained in mixed feature forms when needed. F_k^t F_k=M0(k) is the complete positive unperturbed energy matrix.

For an arbitrary symmetric spatial direction D, define B_D=D-tr(D)I/2. The first feature form is

    S_Q(D)=diag(-B_D,B_D,B_D,-tr(D)/2).

Thus

    M_Q(k,q;D)=diag(-am^2 B_D+C_k^t B_D C_q/a,
                    B_D/a-tr(D)kq^t/(2a^3m^2)).

For independent, potentially noncommuting D,G, put H_B=(B_D B_G+B_G B_D)/2. Expanding exp(B) and exp(-B) to the mixed coefficient gives

    S_QQ(D,G)=diag(H_B,H_B,H_B,tr(D)tr(G)/4).

The mass second vertex is positive, and the scalar constraint second coefficient is nonzero. Both ordered products are required. For tracefree D,G this exactly recovers both frozen S195 vertices.

With linear lapse N=1+n, S_N(n)=n I10, S_NN=0, and the full mixed lapse/spatial form is

    S_DG=S_QQ(D_Q,G_Q)+n_D S_Q(G_Q)+n_G S_Q(D_Q).

There is no lapse-velocity term or hidden A0 oscillator. Shift contacts vanish only in the independent contravariant ADM shift chart described in shift.md; the four-metric chart contact rule in chart.md remains necessary.

The feature formulas extend complex-linearly to individual Fourier amplitudes. Reverse-pair reality is

    M(k,q;D(P))=M(q,k;conj(D(P)))^dagger,

for a real spacetime perturbation with D(-P)=conj(D(P)). A single complex Fourier matrix is not assumed Hermitian. Both the general trace constraint and noncommuting mixed contacts are tested directly from h/sqrt(h),sqrt(h)h^-1 and1/sqrt(h), independently of the B_D definitions.
