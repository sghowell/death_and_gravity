# Proper counterterms and the finite common-field dictionary

The total MSbar coefficients in renormalized
fermion fields are

    delta Zpsi=-(Y/2+a Cf) Ibar/Q,
    delta m_total=(Y-4a Cf)m Ibar/Q,
    delta y_total=(Y-4a Cf)y Ibar/Q.

The common scalar field has its frozen
delta ZPhi=-2NY Ibar/Q, with N=6. Equating
bare coefficients gives

    delta m_bare/m=(3Y/2-3a Cf) Ibar/Q,
    delta y_bare/y=[(N+3/2)Y-3a Cf] Ibar/Q.

The one-loop beta coefficient is twice the
latter residue. At N=6,Cf=4/3 it gives
Q beta_y=y(15Y-8a), exactly matching the
independent matrix calculation in S6.128.
This check includes both external fermion
fields and the scalar field, not just the
proper vertex divergence.

At the named scale define finite coefficients

    z=(-Y J1+a Cf/2)/Q,
    eta=(Y J0+2a Cf)/Q,
    upsilon=[Y(J0+2R)-6a Cf]/Q.

The scalar physical field uses the SAME
complete one-loop w=r_scalar-fp of S6.133.
For a formal loop marker h the selected
zero-momentum canonical dictionary is

    m_zero/m_MS=(1+h eta)/(1+h z),
    y_zero/y_MS=(1+h upsilon)/
                 [(1+h z)sqrt(1+h w)].

Its first-order increments are eta-z and
upsilon-z-w/2. The ratios are not an all-order
resummation; the displayed formal coefficients
are what enter order-by-order matching.
In particular m_zero and y_zero are NOT
relabeled as unchanged MS input parameters.

The proper local mass, kinetic and Yukawa
counterterms belong to the fermionic insertion
-Tr[D_F^-1 D_F^(1)] in S6.137. A counterterm
occurrence assigned inside a renormalized
subgraph is removed from the separate ledger;
it cannot be added twice.

This checkpoint does not integrate those
insertions into the two-loop amplitude.
There is no direct tree H Yukawa vertex.
A finite induced H-fermion form factor, if
needed in other graph decompositions, is
not asserted absent or evaluated here. The
twelve inert flavors retain gauge mass/field
corrections but their scalar Yukawa vertex
is exactly zero, rather than a ratio obtained
by dividing by their vanishing Yukawa.
