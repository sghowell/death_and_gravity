"""Certified parameters of the frozen points P_A and P_B (FORMULATION.md section 1.3).

Exact definitions (M = 1, Lambda = 3/50):
  r_ext   smallest positive root of p_ext(r) = r^3 - 25 r + 25
          (double-root condition 2 Lambda r^3/3 - r + 1 = 0 at Lambda = 3/50);
  Q_max^2 = r_ext - r_ext^4/50  (= r_ext (3 - r_ext)/2 modulo p_ext; both are evaluated);
  Q^2     = c^2 Q_max^2 with c = 999/1000 (P_B) or 249/250 (P_A);
  Delta(r) = -r^4/50 + r^2 - 2 r + Q^2 = -(1/50)(r - r_n)(r - r_-)(r - r_+)(r - r_c).

Everything below is ball arithmetic (Arb, outward rounding) starting from these exact
rationals: r_ext from a certified interval-Newton enclosure of the integer cubic, the four
roots of Delta from interval-Newton enclosures with ball coefficients (four disjoint
enclosures with a unique root each = all roots of the quartic), and then
  kappa_i = |Delta'(r_i)| / (2 r_i^2)                  (CCDHJ18 eq. 4)
  A_i = r_i^2 / Delta'(r_i),  c_+ = sum_{i != +} A_i ln|r_+ - r_i|,  c_c likewise
  r_P = (3 + sqrt(9 - 8 Q^2))/2,  Omega_P = sqrt(F(r_P))/r_P,
  lambda_P = sqrt(F(r_P)(2 - 3/r_P))/r_P,  gamma0 = lambda_P/2   (D6, HV17 section 2.7)
  beta_X = (1/2 + delta) kappa_-  (strip depth of D2),  eta = kappa_-/10  (D7)
  Hintz (2504.01734) Thm 1.1: varkappa(r_+, r_c), lambda_l^+ = (1 + sqrt(1 + 4 l(l+1)/varkappa))/2.
An independent exact-algebra cross-check: the resultant in s of p_ext(s) and Delta(r; Q^2(s))
is an integer polynomial in r; its certified complex roots (fmpz_poly.complex_roots) must
contain the interval-Newton enclosures of r_n, r_-, r_+, r_c (``defining_polynomials``).

Replay:  .venv/bin/python -m p2.params  [write|verify]   (cwd problems/P2/src; ~2 s)
"""
from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from flint import arb, fmpz_poly

from .arbutil import A, Timer, ball_str, git_commit, mid_str, precision

PREC = 256
LAMBDA = Fraction(3, 50)
POINTS = {
    "PB": {"qfrac": Fraction(999, 1000), "delta": Fraction(1, 20),
           "role": "primary if v1.2 is approved (near-extremal window of CCDHJ18 Fig. 3)"},
    "PA": {"qfrac": Fraction(249, 250), "delta": Fraction(1, 100),
           "role": "second target (CCDHJ18 Table I point); sole target if v1.2 is not approved"},
}
CERT_DIR = Path(__file__).resolve().parents[2] / "certificates"
P_EXT = [25, -25, 0, 1]          # ascending: 25 - 25 r + r^3


# ----------------------------------------------------------------------------- root isolation
def poly_eval(co, x):
    """Horner evaluation, co ascending (arb/acb balls or exact)."""
    acc = 0 * x
    for c in reversed(co):
        acc = acc * x + c
    return acc


def poly_deriv(co):
    return [co[k] * k for k in range(1, len(co))]


def interval_newton(co, x0, rad, iters=60):
    """Certified enclosure of the unique real root of sum co_k r^k in [x0 - rad, x0 + rad].

    co: arb coefficients (balls).  The interval Newton operator N(X) = m - p(m)/p'(X) with
    N(X) strictly inside X proves existence and uniqueness of a root in X for every choice
    of coefficients inside the balls (Moore/Krawczyk).  Returns (enclosure, proof_ok)."""
    d = poly_deriv(co)
    m = arb(x0)
    X = m + arb(0, rad)
    proven = False
    for _ in range(iters):
        dX = poly_eval(d, X)
        if dX.contains(0):
            return X, False
        N = m - poly_eval(co, m) / dX
        if not proven:
            if not (N.lower() > X.lower() and N.upper() < X.upper()):
                return X, False
            proven = True
        # intersect and refine
        lo_ = max(N.lower(), X.lower())
        hi_ = min(N.upper(), X.upper())
        Xn = (lo_ + hi_) / 2 + arb(0, 1) * ((hi_ - lo_) / 2)
        if Xn.rad() >= X.rad() * 0.999:
            X = Xn
            break
        X, m = Xn, arb(Xn.mid())
    return X, proven


def mp_real_roots(co_str, dps=90):
    with mp.workdps(dps):
        rts = mp.polyroots([mp.mpf(c) for c in reversed(co_str)], maxsteps=500, extraprec=400)
        return sorted(mp.re(z) for z in rts if abs(mp.im(z)) < mp.mpf(10) ** (-dps // 2))


def isolate_all_real_roots(co_arb, co_mid, degree):
    """All roots real & simple: certified disjoint enclosures, sorted.  Fails loudly otherwise."""
    approx = mp_real_roots(co_mid)
    if len(approx) != degree:
        raise RuntimeError(f"expected {degree} real roots, float found {len(approx)}")
    encl = []
    for k, x0 in enumerate(approx):
        gaps = [abs(x0 - y) for j, y in enumerate(approx) if j != k]
        ok = False
        for rad in ("1e-50", "1e-30", "1e-15", "1e-8"):
            if mp.mpf(rad) * 4 >= min(gaps):
                break
            X, ok = interval_newton(co_arb, arb(mp.nstr(x0, 80)), arb(rad))
            if ok:
                break
        if not ok:
            raise RuntimeError(f"interval Newton failed near {x0}")
        encl.append(X)
    for a, b in itertools.pairwise(encl):
        if not (a.upper() < b.lower()):
            raise RuntimeError("root enclosures not disjoint")
    return encl


# ----------------------------------------------------------------------------- the point
class Point:
    """Certified Arb balls for a frozen parameter point (M = 1, Lambda = 3/50)."""

    def __init__(self, name: str, prec: int = PREC):
        spec = POINTS[name]
        self.name, self.prec = name, prec
        self.qfrac, self.delta = spec["qfrac"], spec["delta"]
        with precision(prec):
            self._build()

    def _build(self):
        c = A(self.qfrac)
        # r_ext: the three real roots of the integer cubic, smallest positive one
        cub = [A(v) for v in P_EXT]
        cub_roots = isolate_all_real_roots(cub, [str(v) for v in P_EXT], 3)
        pos = [x for x in cub_roots if x > 0]
        self.cubic_roots = cub_roots
        self.r_ext = pos[0]
        self.Qmax2 = self.r_ext - self.r_ext ** 4 / 50
        self.Qmax2_reduced = self.r_ext * (3 - self.r_ext) / 2
        self.Q2 = c * c * self.Qmax2
        self.Qmax = self.Qmax2.sqrt()
        self.Q = c * self.Qmax
        lam = A(LAMBDA)
        self.Lam = lam
        # Delta, ascending coefficients
        self.dco = [self.Q2, arb(-2), arb(1), arb(0), -lam / 3]
        mid = [x.mid().str(85, radius=False) for x in self.dco]
        self.roots = isolate_all_real_roots(self.dco, mid, 4)
        self.rn, self.rm, self.rp, self.rc = self.roots
        self.ddco = poly_deriv(self.dco)
        d1 = [poly_eval(self.ddco, r) for r in self.roots]
        self.dDelta_at = d1
        self.kap = {k: abs(dv) / (2 * r * r) for k, r, dv in zip(("n", "-", "+", "c"), self.roots, d1)}
        self.km, self.kp, self.kc = self.kap["-"], self.kap["+"], self.kap["c"]
        self.kc_dS = (lam / 3).sqrt()
        # partial fractions 1/F = r^2/Delta = sum_i A_i/(r - r_i)
        self.Acoef = [r * r / dv for r, dv in zip(self.roots, d1)]
        self.c_plus = sum((self.Acoef[i] * abs(self.rp - self.roots[i]).log() for i in (0, 1, 3)), arb(0))
        self.c_cos = sum((self.Acoef[i] * abs(self.rc - self.roots[i]).log() for i in (0, 1, 2)), arb(0))
        # photon sphere and trapping constants (D6)
        disc = 9 - 8 * self.Q2
        self.photon_disc = disc
        self.rP = (3 + disc.sqrt()) / 2
        self.FP = self.F(self.rP)
        self.OmegaP = self.FP.sqrt() / self.rP
        self.lamP = (self.FP * (2 - 3 / self.rP)).sqrt() / self.rP
        self.gamma0 = self.lamP / 2
        # strip (D2) and boxes (D7)
        self.beta = (A(Fraction(1, 2)) + A(self.delta)) * self.km
        self.eta = self.km / 10
        # Hintz Thm 1.1 near-horizon limit with r_e = r_+, r_c current
        re_, rc_ = self.rp, self.rc
        self.varkappa = (rc_ ** 2 + 2 * re_ * rc_ - 3 * re_ ** 2) / (rc_ ** 2 + 2 * re_ * rc_ + 3 * re_ ** 2)
        self.lam_plus = [(1 + (1 + 4 * l * (l + 1) / self.varkappa).sqrt()) / 2 for l in range(4)]

    # -- functions of r (ball arithmetic) -------------------------------------------------
    def Delta(self, r):
        return poly_eval(self.dco, r)

    def dDelta(self, r):
        return poly_eval(self.ddco, r)

    def F(self, r):
        return self.Delta(r) / (r * r)

    def V(self, r, l):
        """V_l = F (L/r^2 + F'/r), CCDHJ18 eq. 7."""
        L = l * (l + 1)
        Fr = (self.dDelta(r) * r - 2 * self.Delta(r)) / r ** 3
        return self.F(r) * (L / (r * r) + Fr / r)

    # -- certified statements ---------------------------------------------------------------
    def inequalities(self):
        """Each entry: (statement, certified bool, margin ball = rhs - lhs)."""
        a = self.beta
        out = {}

        def ineq(key, text, lhs, rhs):
            out[key] = {"statement": text, "certified": bool(lhs < rhs), "margin": rhs - lhs}

        ineq("order_rn_0", "r_n < 0", self.rn, arb(0))
        ineq("order_0_rm", "0 < r_-", arb(0), self.rm)
        ineq("order_rm_rp", "r_- < r_+", self.rm, self.rp)
        ineq("order_rp_rc", "r_+ < r_c", self.rp, self.rc)
        ineq("photon_disc", "9 - 8 Q^2 > 0", arb(0), self.photon_disc)
        ineq("rP_exterior_lo", "r_+ < r_P", self.rp, self.rP)
        ineq("rP_exterior_hi", "r_P < r_c", self.rP, self.rc)
        ineq("F_rP_pos", "F(r_P) > 0", arb(0), self.FP)
        ineq("trap_pos", "2 - 3/r_P > 0", arb(0), 2 - 3 / self.rP)
        ineq("strip_vs_kp", "(1/2+delta) kappa_- < kappa_+  (E_+ misses S_delta; Jost margin)", a, self.kp)
        ineq("strip_vs_kc", "(1/2+delta) kappa_- < kappa_c  (E_c misses S_delta; Jost margin)", a, self.kc)
        ineq("strip_vs_gamma0", "(1/2+delta) kappa_- < gamma0  (HV17: alpha < gamma0)", a, self.gamma0)
        ineq("CF17_2kp", "kappa_- < 2 kappa_+  (CF17 Thm 3.1: 2 min{p, kappa_+} > kappa_-)", self.km, 2 * self.kp)
        ineq("CF17_2alpha", "kappa_- < 2 (1/2+delta) kappa_-", self.km, 2 * a)
        ineq("Ros25_rho", "rho = kappa_-/kappa_+ < 2  (Ros25)", self.km / self.kp, arb(2))
        rho = self.km / self.kp
        out["CGNS18_rho"] = {"statement": "rho = kappa_-/kappa_+ < 9/7  (CGNS18; context only)",
                             "certified": bool(rho < A(Fraction(9, 7))),
                             "certified_false": bool(rho > A(Fraction(9, 7))),
                             "margin": A(Fraction(9, 7)) - rho}
        # sum of roots = 0 (no cubic term): consistency of the enclosures
        s = self.rn + self.rm + self.rp + self.rc
        out["vieta_sum"] = {"statement": "r_n + r_- + r_+ + r_c contains 0", "certified": bool(s.contains(0)),
                            "margin": s}
        out["Qmax2_two_forms"] = {"statement": "r_ext - r_ext^4/50 and r_ext(3 - r_ext)/2 overlap",
                                  "certified": bool(self.Qmax2.overlaps(self.Qmax2_reduced)),
                                  "margin": self.Qmax2 - self.Qmax2_reduced}
        return out

    def quantities(self):
        km = self.km
        return {
            "r_ext": self.r_ext, "Qmax2": self.Qmax2, "Qmax": self.Qmax, "Q2": self.Q2, "Q": self.Q,
            "r_n": self.rn, "r_minus": self.rm, "r_plus": self.rp, "r_c": self.rc,
            "Delta_prime_at_roots": dict(zip(("r_n", "r_minus", "r_plus", "r_c"), self.dDelta_at)),
            "kappa_minus": km, "kappa_plus": self.kp, "kappa_c": self.kc, "kappa_n": self.kap["n"],
            "kappa_c_dS": self.kc_dS,
            "kappa_plus_over_kappa_minus": self.kp / km, "kappa_c_over_kappa_minus": self.kc / km,
            "rho": km / self.kp,
            "A_coefficients": dict(zip(("r_n", "r_minus", "r_plus", "r_c"), self.Acoef)),
            "c_plus": self.c_plus, "c_c": self.c_cos,
            "r_P": self.rP, "F_rP": self.FP, "Omega_P": self.OmegaP, "Omega_P_over_kappa_minus": self.OmegaP / km,
            "lambda_P": self.lamP, "gamma0": self.gamma0, "gamma0_over_kappa_minus": self.gamma0 / km,
            "strip_depth_beta": self.beta, "strip_depth_over_gamma0": self.beta / self.gamma0,
            "eta": self.eta,
            "E_plus_first": self.kp, "E_c_first": self.kc,
            "jost_margin_plus_over_kappa_minus": (self.kp - self.beta) / km,
            "varkappa_Hintz": self.varkappa,
            "lambda_plus_l": {str(l): v for l, v in enumerate(self.lam_plus)},
        }


# ----------------------------------------------------------------------------- exact algebra
def defining_polynomial(qfrac: Fraction):
    """Integer polynomial (sympy, factored over Q) whose roots include r_n, r_-, r_+, r_c.

    Res_s( s^3 - 25 s + 25 , Delta(r; Q^2 = (c^2/2) s (3 - s)) ), c = qfrac.  The reduction
    Q_max^2 = s - s^4/50 = s(3 - s)/2 mod p_ext is checked symbolically as well."""
    import sympy as sp
    s, r = sp.symbols("s r")
    pext = s ** 3 - 25 * s + 25
    red = sp.rem(sp.expand(s - s ** 4 / 50 - s * (3 - s) / 2), pext, s)
    assert sp.simplify(red) == 0
    c = sp.Rational(qfrac.numerator, qfrac.denominator)
    Q2 = c ** 2 * s * (3 - s) / 2
    Delta = -sp.Rational(1, 50) * r ** 4 + r ** 2 - 2 * r + Q2
    res = sp.resultant(pext, Delta, s)
    res = sp.Poly(sp.together(res) * sp.denom(sp.together(res)), r)
    _content, prim = res.primitive()
    factors = sp.factor_list(prim.as_expr(), r)[1]
    return [sp.Poly(f, r) for f, _mult in factors]


def check_defining_polynomials(pt: Point):
    """Each root ball of Delta must meet exactly one certified root of exactly one factor."""
    facs = defining_polynomial(pt.qfrac)
    report = []
    with precision(pt.prec):
        for f in facs:
            co = [int(x) for x in reversed(f.all_coeffs())]
            P = fmpz_poly(co)
            rts = [z for z, _m in P.complex_roots()]
            hits = {}
            for name, rb in zip(("r_n", "r_minus", "r_plus", "r_c"), pt.roots):
                n_hit = sum(1 for z in rts if z.imag.contains(0) and z.real.overlaps(rb))
                if n_hit:
                    hits[name] = n_hit
            report.append({"degree": f.degree(), "coefficients_ascending": [str(x) for x in co], "hits": hits})
    names = ("r_n", "r_minus", "r_plus", "r_c")
    ok = all(sum(rep["hits"].get(nm, 0) for rep in report) == 1 for nm in names)
    return ok, report


# ----------------------------------------------------------------------------- certificate I/O
def _ser(x):
    if isinstance(x, arb):
        return {"ball": ball_str(x, 40), "approx": mid_str(x, 20)}
    if isinstance(x, dict):
        return {k: _ser(v) for k, v in x.items()}
    return x


def certificate(prec: int = PREC):
    out = {"claim": "P2 S1 params (FORMULATION section 1.3)", "prec_bits": prec,
           "conventions": "M = 1, Lambda = 3/50; Delta(r) = -r^4/50 + r^2 - 2r + Q^2; kappa_i = |Delta'(r_i)|/(2 r_i^2)",
           "definitions": {
               "r_ext": "smallest positive root of r^3 - 25 r + 25",
               "Qmax2": "r_ext - r_ext^4/50",
               "Q2": "(qfrac)^2 Qmax2",
               "photon_sphere": "r_P = (3 + sqrt(9 - 8 Q^2))/2",
               "gamma0": "sqrt(F(r_P)(2 - 3/r_P))/(2 r_P)",
               "strip_depth": "(1/2 + delta) kappa_-",
               "varkappa_Hintz": "(r_c^2 + 2 r_+ r_c - 3 r_+^2)/(r_c^2 + 2 r_+ r_c + 3 r_+^2)",
               "lambda_plus_l": "(1 + sqrt(1 + 4 l(l+1)/varkappa))/2"},
           "method": "interval Newton (Arb, outward rounding) on the integer cubic and on Delta with "
                     "ball coefficients; four disjoint unique-root enclosures; resultant cross-check",
           "commit": git_commit(), "points": {}}
    for name in ("PB", "PA"):
        with Timer() as tm:
            pt = Point(name, prec)
            with precision(prec):
                ineq = pt.inequalities()
            okdef, rep = check_defining_polynomials(pt)
        with precision(prec):
            entry = _entry(name, pt, ineq, okdef, rep, tm)
        out["points"][name] = entry
    return out


def _entry(name, pt, ineq, okdef, rep, tm):
    return {"qfrac": str(pt.qfrac), "delta": str(pt.delta), "role": POINTS[name]["role"],
            "quantities": _ser(pt.quantities()),
            "inequalities": {k: {"statement": v["statement"], "certified": v["certified"],
                                 **({"certified_false": v["certified_false"]} if "certified_false" in v else {}),
                                 "margin": _ser(v["margin"])} for k, v in ineq.items()},
            "defining_polynomials": {"all_roots_located": okdef, "factors": rep},
            "timing_s": {"wall": round(tm.wall, 3), "cpu": round(tm.cpu, 3)}}


def write(path: Path | None = None):
    path = path or CERT_DIR / "params.json"
    cert = certificate()
    path.write_text(json.dumps(cert, indent=1) + "\n")
    return cert


def verify(path: Path | None = None) -> list[str]:
    """Replay: recompute every ball at the stored precision; each must lie inside the stored
    bounds (identical algorithm => identical balls), and every certified verdict must recur."""
    path = path or CERT_DIR / "params.json"
    stored = json.loads(path.read_text())
    fresh = certificate(stored["prec_bits"])
    problems = []

    def walk(a, b, where):
        if isinstance(a, dict) and "ball" in a and isinstance(a.get("ball"), str):
            sa, sb = a["ball"], b["ball"]
            if sa != sb:
                problems.append(f"{where}: stored {sa} != replay {sb}")
            return
        if isinstance(a, dict):
            for k in a:
                if k in ("timing_s", "commit"):
                    continue
                if k not in b:
                    problems.append(f"{where}.{k}: missing in replay")
                    continue
                walk(a[k], b[k], f"{where}.{k}")
        elif a != b:
            problems.append(f"{where}: stored {a!r} != replay {b!r}")

    walk({k: v for k, v in stored.items() if k not in ("commit",)},
         {k: v for k, v in fresh.items() if k not in ("commit",)}, "params")
    return problems


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "verify"
    if cmd == "write":
        c = write()
        for nm, e in c["points"].items():
            bad = [k for k, v in e["inequalities"].items() if not v["certified"] and k != "CGNS18_rho"]
            print(nm, "uncertified inequalities:", bad, "defining polys ok:", e["defining_polynomials"]["all_roots_located"])
    else:
        probs = verify()
        print("params.json replay:", "OK" if not probs else "MISMATCH")
        for p in probs:
            print("  ", p)
        sys.exit(1 if probs else 0)
