# Actual first-order current and canonical normalization

For physical-frame fields E=pi/a^2, B=curl A/a^2,
V=mA/a, the complete spatial Proca stress is

    Tij=-EiEj-BiBj+ViVj
        +deltaij(E^2+B^2+V0^2-V^2)/2.

The temporal constraint is retained in V0 and in the complete
state; its direct trace term drops only after contracting with
a tracefree direction D. The actual first Hamiltonian vertex is

    H_D = [pi.D.pi/a + (curl A).D.(curl A)/a
           -a m^2 A.D.A]/2
        = -a^3 T:D/2.

Thus J_D=-H_D. With source Gamma strictly earlier than D,
the finite-band current from S195 is i<[H_D,H_Gamma]>
with no metric contact. After smearing it is exactly
B_K=i<[T_K[D],T_K[Gamma]]>/4.

The common-projection limit in notes/projection.md proves
existence, uniqueness for this regulator sequence and its
quantitative full-momentum error. On the separated product
supports the Wick commutator is already defined; no extension
of a coincident retarded product is needed. Any unchanged local
curvature counterterm contributes zero here. This identifies
the physical weak first-order response in this sector.

Write c for S186's unrounded complete stress variance coefficient.
It obeys c<1e50. Then |B|<=c N[D]N[Gamma]/2<5e49 NN.
For h=sqrt(kappa)gamma/2, each metric direction contributes
2/sqrt(kappa). The weak proper force is paired with a^3 dt dx;
no extra volume factor is discarded. The resulting bounds are

    |B_can| <=2c NN/kappa <2e-750 NN;
    |B_can-B_can,K| <2e-748 NN/K.

At K1e16 the second display is2e-764 NN. These inequalities
have exact zero when a test is zero. The canonical factors are
those of the tensor sector, not a fully constraint-reduced
mixed scalar/metric norm.
