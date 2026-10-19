"""A rigorous count at l = 20 that shares no code or argument with S4 (S4 review).

S3's reviewed counting pipeline (p2.part3, zsol box enclosures + Cauchy discs + argument principle) is run unchanged
at l = 20 on R'_20 = [-a_20, a_20] x [-b, e], with a_20 the CERTIFIED P2-2 threshold (certificates/lg/lg_PB_l20.json),
b >= (1/2 + delta) kappa_-, e = kappa_-/10 rounded down.  part3's working directory is redirected to the output
directory given on the command line (no file under problems/P2/results or certificates is written;
assemble(write=False)).
Winding 0 on dR'_20, with P2-2 (|Re w| >= a_20) and P2-0 (Im w > 0), gives: no zero of W_20 in S_delta --
independently of S4's regimes, contours, W1/W2 and certificates.  l = 20 is S4's hardest l (largest h).
Usage: count_l20.py OUTDIR [workers]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(W / "src"))

from p2 import part3  # noqa: E402


def main():
    outdir = Path(sys.argv[1]).resolve()
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    outdir.mkdir(parents=True, exist_ok=True)
    part3._RES["dir"] = outdir
    part3.run("PB", [20], workers=workers, launch={"note": "S4 review: independent rigorous count at l = 20"})
    res = part3.assemble("PB", 20, write=False)
    summary = {k: v for k, v in res.items() if k != "pieces"}
    (outdir / "count_l20_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    (outdir / "count_l20_full.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
