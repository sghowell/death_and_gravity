# Complete original reciprocal, positive cut and two moment identities

Let z be the unique upper first-sheet zero, R=1/[z D_C'(z)], and E_C(p)=1/[p D_C(p)]. The original cut values are A2(-tau+i0)=D_A(tau)+i pi U(tau), with U>0.

The full reciprocal has a massless pole with residue1/C, the simple conjugate pair z and conjugate(z), and the original cut. Direct algebra gives

rho_E(tau)=-Im E_C(-tau+i0)/pi
 =U/[(C-tau D_A)^2+pi^2 tau^2 U^2]>0.                        (1)

The complete first-sheet Cauchy representation is

E_C(p)=1/(C p)+R/(p-z)+conjugate(R)/(p-conjugate(z))
       +integral_(4m^2)^infinity rho_E(tau)/(p+tau) dtau.     (2)

The cut sign follows from (1), not an analogy to an isolated form factor. At threshold D_C tends to C+688m^2/225>0 and U has square-root onset, so rho_E=O(sqrt(tau-4m^2)). At infinity the original massive-factor expansion gives rho_E=O(1/[tau^2 log^2(tau/m^2)]). Hence both M0=integral rho_E and M1=integral tau rho_E are finite and strictly positive.

The full exterior estimate gives E_C=O(1/[p^2 log(p/m^2)]) uniformly on a sufficiently large closing circle. The large contour therefore contributes no entire subtraction polynomial. The small threshold indentation vanishes. Residues and the discontinuity then give (2); no pole or arc is discarded.

Take p to positive infinity in p times (2). Dominated convergence uses integrable rho_E and gives

1/C+2 Re R+M0=0.                                            (3)

After using (3), the exact identities

1/(p-z)=1/p+z/[p(p-z)],
1/(p+tau)=1/p-tau/[p(p+tau)]

permit multiplying by p^2. Dominated convergence now uses integrable tau rho_E, and p^2 E_C tends to0, giving the second identity

2 Re(R z)-M1=0.                                             (4)

These cancellations are consequences of the full massive reciprocal. A cut-only or massless-plus-cut reciprocal violates (3): its coefficients would all be positive. The conjugate pair is indispensable.

At the real zero p=-r m^2 of the ISOLATED S223 factor A2, the full D_C equals C, so E_C=1/(Cp) is finite. This is the actual tree-plus-loop composition, not a modified contour prescription or manual deletion.

Combining (3) with the exact residue bound in notes/bounds.md,

M0=-1/C-2 Re R<=2|R|-1/C,
1/C+2|R|+M0=2(|R|-Re R)<=4|R|<8/C.                          (5)

These finite-mass sum rules supply the uniform time-domain bounds. Diagnostic numerical integration of the resonant cut and both moments is independent evidence only; the original positive massive density in (1) is retained exactly.
