"""S4-5 review, G1 centre half: kappa-dependence of the t = 0 data (P_0, G_0) of
the scaled 4D system at the regular centre (lincentre.RegularFamily), on the
certified (a*, mu*) centre series.  Prints P_0, the kappa^1 and kappa^2 parts
of G_0, and det(n P_0 - G_0(kappa)) as a polynomial in kappa for n = 1..6.
"""
import sys

from flint import acb, acb_mat, arb

from p4.validated import ec_centre, linscaled
from p4.validated.arbseries import precision
from p4.validated.lincentre import coefficient_series

with precision(256):
    ce = ec_centre.certified_centre_expansion(256, 30)
    L = linscaled.LinSystem() if hasattr(linscaled, "LinSystem") else None
    from p4.validated.linsys import LinSystem
    S4 = linscaled.full_system(LinSystem())
    Pser, Gser = coefficient_series(S4, ce)
    d = S4.d
    P0 = [[sum((Pser[j][r][c][0] for j in range(1)), arb(0)) for c in range(d)] for r in range(d)]
    print("P kappa-degrees:", len(Pser), " G kappa-degrees:", len(Gser))
    for j, M in enumerate(Pser):
        print(f"P_j={j} at t=0:", [[M[r][c][0].str(6) for c in range(d)] for r in range(d)])
    for j, M in enumerate(Gser):
        print(f"G_j={j} at t=0:", [[M[r][c][0].str(6) for c in range(d)] for r in range(d)])

    def mat(Mlist, k):
        kp = [acb(1), acb(k), acb(k) * acb(k)]
        return acb_mat([[sum((kp[j] * Mlist[j][r][c][0] for j in range(len(Mlist))), acb(0))
                         for c in range(d)] for r in range(d)])
    for k in (0, 1, 7, acb(3, 50), acb(0, 1000)):
        P0m, G0m = mat(Pser, k), mat(Gser, k)
        B0 = P0m.inv() * G0m
        # characteristic polynomial via eigenvalues (balls)
        try:
            ev = B0.eig(multiple=True)
        except Exception as e:  # noqa: BLE001
            ev = str(e)
        dets = [(n * P0m - G0m).det().str(6, radius=True) for n in (1, 2, 3, 4, 5, 6)]
        print("kappa", k, "eig(P0^-1 G0):", [e.str(6, radius=True) for e in ev] if not isinstance(ev, str) else ev)
        print("    det(n P0 - G0), n=1..6:", dets)
