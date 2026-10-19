"""The pilot discrepancy (S4 review): the pre-registration's S pilot reported the leading-order margin
min_t q0/|Y'|^2 (eps = 0) for the best softplus design; the certificates bound q^ including the O(eps) terms over
eps in [0, eps_max] and |nu| <= d.  For the stored (compact-family) designs of the top leaves this splits the
difference: (a) leading order, eps = 0, normalised by |Y'|^2, at zeta_hi; (b) the same without normalisation;
(c) the full q^ vertex minimum over the certificate's parameter box (as in indep_S); (d) the certified bound."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import indep_S as S  # noqa: E402

out = []
for side, tag in ((1, "plus"), (-1, "c")):
    cert = json.loads((S.CERT / f"S_PB_{tag}.json").read_text())
    for x in cert["leaves"]:
        if x.get("closing") or Fraction(x["zeta_over_OmegaP"][0]) < Fraction(8, 10):
            continue
        w, ya, phinf, pc, pw = S.parse_design(x["design"])
        Xa, k = float(Fraction(x["X_a"])), float(Fraction(x["k"]))
        zl, zh = (float(Fraction(v)) * S.OmP for v in x["zeta_over_OmegaP"])
        em = float(Fraction(x["eps_max"]))
        t, sv, ra = S.integrate_contour(side, Xa, k, w, ya, -8.0)
        U, Wt, r = S.UW_s(sv, side)
        eta, e1, e2 = S.eta_all(t, w)
        phi, d1 = S.phase_all(t, phinf, pc, pw)
        Yp2 = np.abs(1 + 1j * e1) ** 2
        G = np.exp(-1j * phi) * (1 + 1j * e1) ** 2
        q0_hi = (G * (U / zh ** 2 - 1)).real
        q0_lo = (G * (U / zl ** 2 - 1)).real
        lead = min((q0_hi / Yp2).min(), (q0_lo / Yp2).min())
        lead_raw = min(q0_hi.min(), q0_lo.min())
        full = S.qhat_vertices(t, U, Wt, w, phinf, pc, pw, k, [1 / zh ** 2, 1 / zl ** 2], [0, em * S.d], [0, (em * S.d) ** 2],
                               [0, em * em]).min()
        cert_mu = float(x["eps_boxes"][0]["mu_cells"].strip("[]").split(",")[0])
        rec = {"side": tag, "zeta": x["zeta_over_OmegaP"], "eps_max": em, "leading_order_over_Yp2": float(lead),
               "leading_order": float(lead_raw), "full_vertex_min": float(full), "certified": cert_mu,
               "float_design_margin": x["float_margin"]}
        out.append(rec)
        print(rec)
Path(sys.argv[1] if len(sys.argv) > 1 else "pilot_compare.json").write_text(json.dumps(out, indent=1) + "\n")
