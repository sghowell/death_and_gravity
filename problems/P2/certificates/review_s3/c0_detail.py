"""S3 review: two details of the C0 deviation (3e7fca05).
(a) The test point dropped from tests/test_part3.py, (7, 40.0, -0.55), with S1's arbw UNWRAPPED at
    256 bits (the 89a2336b test): does it raise?
(b) The three points where unwrapped arbw raises at 256 bits, with the C0 wrapper (more Frobenius
    terms) but still at 256 bits: is the ball too wide for the overlap to mean anything?
Writes c0_detail.json.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from flint import acb, arb  # noqa: E402
from p2 import part3  # noqa: E402
import p2.arbw as arbw_mod  # noqa: E402
from p2.arbutil import precision  # noqa: E402
from p2.params import Point  # noqa: E402

pt = Point("PB", 256)
out = {}
with precision(256):
    w = acb(arb(part3._kq(pt, 40.0)), arb(part3._kq(pt, -0.55)))
    try:
        Wa = arbw_mod.wronskian(pt, 7, w)
        out["a_unwrapped_7_40"] = {"raised": None, "relrad": float(abs(Wa).rad() / abs(Wa).mid())}
    except ArithmeticError as ex:
        out["a_unwrapped_7_40"] = {"raised": str(ex)}
orig = arbw_mod.frobenius


def frob(ptx, which, w, L, t, min_terms=20):
    kap = ptx.kp if which == "+" else ptx.kc
    ra = float((abs(acb(w)) / (2 * kap)).upper())
    return orig(ptx, which, w, L, t, min_terms=max(min_terms, int(2 * ra) + 8))


arbw_mod.frobenius = frob
out["b_wrapped_256"] = []
for l, re_k, im_k in [(10, 56.5747, -0.55), (19, 71.1985, 0.1), (19, 101.7122, -0.55)]:
    with precision(256):
        w = acb(arb(part3._kq(pt, re_k)), arb(part3._kq(pt, im_k)))
        Wa = arbw_mod.wronskian(pt, l, w)
        out["b_wrapped_256"].append({"l": l, "w_kappa": [re_k, im_k], "W": Wa.str(8),
                                     "relrad": float(abs(Wa).rad() / abs(Wa).mid()) if abs(Wa).mid() != 0 else None,
                                     "contains_0": bool(Wa.contains(0))})
arbw_mod.frobenius = orig
print(json.dumps(out, indent=1))
(HERE / "c0_detail.json").write_text(json.dumps(out, indent=1) + "\n")
