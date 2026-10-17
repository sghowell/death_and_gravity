"""S4-5 review: the join of tc_energy's tube region (steps with float(x) - h >= -4.5)
with Lemma D's region y = e^x in (0, YCUT]; exact dyadic step ends vs ln(YCUT)."""
import sys
from flint import arb, ctx
from p4.validated import lintube
from p4.validated.tc_commutator import YCUT
ctx.prec = 256
tube = lintube.Tube.load(sys.argv[1])
inc = [sd for sd in tube.steps if not (float(sd.x) - sd.h < -4.5)]
lowest = min(inc, key=lambda sd: float(sd.x) - sd.h)
xlow = lowest.x - arb(lowest.h)
print("steps total", len(tube.steps), "included", len(inc))
print("lowest included step: x =", lowest.x.str(25, radius=True), " h =", repr(lowest.h))
print("exact lower end x_k - h_k =", xlow.str(25, radius=True), " rad(x_k) =", lowest.x.rad())
ly = arb(YCUT).log()
print("ln(YCUT) =", ly.str(25, radius=True), " YCUT =", repr(YCUT))
print("e^{-4.5} =", arb(-4.5).exp().str(25, radius=True))
print("YCUT >= e^{-4.5}?", arb(YCUT) >= arb(-4.5).exp(), " YCUT < e^{-4.5}?", arb(YCUT) < arb(-4.5).exp())
print("x_low <= ln(YCUT) (no gap)?", xlow <= ly, "  x_low > ln YCUT (gap)?", xlow > ly)
print("gap/overlap size x_low - ln(YCUT):", (xlow - ly).str(10, radius=True))
# next excluded step just below
exc = [sd for sd in tube.steps if float(sd.x) - sd.h < -4.5]
hi_exc = max(exc, key=lambda sd: float(sd.x))
print("highest excluded step: x =", hi_exc.x.str(25), " h =", hi_exc.h, " float(x)-h =", repr(float(hi_exc.x) - hi_exc.h))
