"""Validated rays x = x_P + e^{i theta} s in the tortoise plane and their images r(s).

r(x) is the inverse of r*(r) = sum_i A_i log(r - r_i) (exterior branch, real on (r_+, r_c)),
so along a ray dr/ds = e^{i theta} F(r), F = Delta/r^2, r(0) = r_P (x_P := r*(r_P)).

Each step [s0, s0 + h] (h real, either sign) carries:
  * an a-priori box B with r0 + [0, h] e^{i theta} F(B) contained in B  (Picard: r(s) in B on
    the whole step; B must avoid r = 0, where F has its only finite singularity);
  * Taylor coefficients r_0..r_N at s0 (from the ball r0) and the (N+1)-th coefficient
    evaluated with the box B as initial value, so that, componentwise (Lagrange remainder for
    the real-analytic Re r and Im r),
        r(s0 + t) in sum_{k<=N} r_k t^k + t^{N+1} [r_{N+1}(B)],   t between 0 and h.
The recursion: u = 1/r, F = -(Lambda/3) r^2 + 1 - 2u + Q^2 u^2, r_{k+1} = e^{i theta} F_k/(k+1).

Uncertainty of the initial value (Lohner-type, wrapping-free): for an autonomous scalar ODE the
derivative of the flow with respect to the initial value is dr(t)/dr0 = F(r(t))/F(r0); for w in
the ball r0 with midpoint m0, r(t; w) - r(t; m0) lies in conv{F(r(t; w'))/F(w')}(w - m0), which is
enclosed by the box F(R_t)/F(r0) times the centred ball (R_t encloses all r(t; w')).  The Taylor
model is propagated from the exact midpoint m0.  F is evaluated in factored form
-(Lambda/3)(r - r_n)(r - r_-)(r - r_+)(r - r_c)/r^2 to avoid cancellation near the horizons.

(The sector analyticity used by the certificate is proved by the Re r* identity, not by a fan of
rays; see notes/s2-m3.md, Lemma 1.)
"""
from __future__ import annotations

from flint import acb, arb

N_TAYLOR = 20


def _unit(r):
    """Ball [0, r] for real r (arb with that span)."""
    return (r / 2) + arb(0, 1) * (abs(r) / 2)


def taylor_coeffs(pt, r0, E, N):
    """Coefficients r_0..r_N of the solution through r0 (acb ball or box)."""
    lam3 = pt.Lam / 3
    Q2 = acb(pt.Q2)
    r = [acb(r0)]
    u = [1 / r[0]]
    for k in range(N):
        r2k = sum((r[j] * r[k - j] for j in range(k + 1)), acb(0))
        u2k = sum((u[j] * u[k - j] for j in range(k + 1)), acb(0))
        Fk = -lam3 * r2k - 2 * u[k] + Q2 * u2k + (1 if k == 0 else 0)
        r.append(E * Fk / (k + 1))
        n = k + 1
        u.append(-sum((r[j] * u[n - j] for j in range(1, n + 1)), acb(0)) / r[0])
    return r


def F_of(pt, r):
    """F = Delta/r^2 in factored form (no cancellation near the roots)."""
    prod = acb(-pt.Lam / 3)
    for ri in pt.roots:
        prod = prod * (r - ri)
    return prod / (r * r)


def apriori_box(pt, r0, E, h, tries=12):
    """Box B with r0 + [0,h] E F(B) inside B, or None."""
    f0 = abs(F_of(pt, r0)).upper()
    rad = 2 * abs(arb(h)) * f0 + arb("1e-30")
    for _ in range(tries):
        B = acb(r0.real + arb(0, 1) * rad, r0.imag + arb(0, 1) * rad)
        if B.real.contains(0) and B.imag.contains(0):
            return None
        TB = r0 + _unit(arb(h)) * E * F_of(pt, B)
        if (TB.real.lower() > B.real.lower() and TB.real.upper() < B.real.upper()
                and TB.imag.lower() > B.imag.lower() and TB.imag.upper() < B.imag.upper()):
            return B
        rad = rad * 2
    return None


class Step:
    __slots__ = ("B", "D0", "J", "N", "co", "h", "rem", "s0")

    def __init__(self, s0, h, co, rem, B, N, D0, J):
        self.s0, self.h, self.co, self.rem, self.B, self.N = s0, h, co, rem, B, N
        self.D0, self.J = D0, J          # centred initial ball and flow-derivative box over the step

    def eval_mid(self, t):
        """Enclosure of the midpoint solution r(s0 + t; m0), t an arb between 0 and h."""
        acc = acb(0)
        for c in reversed(self.co):
            acc = acc * t + c
        tn = t ** (self.N + 1)
        rr = self.rem
        return acc + acb(rr.real * tn, rr.imag * tn)

    def eval(self, t):
        """Enclosure of r(s0 + t) for every initial value in the step's initial ball."""
        return self.eval_mid(t) + self.J * self.D0


class Ray:
    """Validated integration of r along x = x_P + E s for s in [0, S] (S may be negative)."""

    def __init__(self, pt, E, S, h, N=N_TAYLOR):
        self.pt, self.E, self.N = pt, acb(E), N
        self.steps = []
        r0 = acb(pt.rP)
        s = arb(0)
        S = arb(S)
        h = arb(h)
        sgn = 1 if S > 0 else -1
        while (S - s) * sgn > 0:
            hh = h if (S - s) * sgn > abs(h) else S - s
            hh = abs(hh) * sgn
            B = apriori_box(pt, r0, self.E, hh)          # valid for every initial value in r0
            if B is None:
                raise ArithmeticError(f"no a-priori box at s = {s.mid()} (h = {hh})")
            m0 = acb(r0.real.mid(), r0.imag.mid())
            D0 = r0 - m0
            co = taylor_coeffs(pt, m0, self.E, N)
            rem = taylor_coeffs(pt, B, self.E, N + 1)[N + 1]
            F0 = F_of(pt, r0)
            Jcrude = F_of(pt, B) / F0                    # flow derivative over the whole step
            st = Step(s, hh, co, rem, B, N, D0, Jcrude)
            Rh = st.eval_mid(hh) + Jcrude * D0
            Jh = F_of(pt, Rh) / F0                       # refined at the step end
            r0 = st.eval_mid(hh) + Jh * D0
            self.steps.append(st)
            s = s + hh
        self.S = S

    def step_index(self, s):
        """Index of the step containing s (float).  Steps are uniform except the last one."""
        h0 = abs(float(self.steps[0].h.mid()))
        k = int(abs(s) / h0 + 1e-12)
        k = min(max(k, 0), len(self.steps) - 1)
        st = self.steps[k]
        a, b = float(st.s0.mid()), float((st.s0 + st.h).mid())
        lo_, hi_ = min(a, b), max(a, b)
        if not (lo_ - 1e-12 <= s <= hi_ + 1e-12):
            for i, st in enumerate(self.steps):
                a, b = float(st.s0.mid()), float((st.s0 + st.h).mid())
                if min(a, b) - 1e-12 <= s <= max(a, b) + 1e-12:
                    return i
            raise IndexError(s)
        return k

    def enclose(self, s_lo, s_hi):
        """Enclosure of r over [s_lo, s_hi] (arbs), which must lie in one step."""
        i = self.step_index(float(((s_lo + s_hi) / 2).mid()))
        st = self.steps[i]
        t_lo, t_hi = s_lo - st.s0, s_hi - st.s0
        T = (t_lo + t_hi) / 2 + arb(0, 1) * (abs(t_hi - t_lo) / 2)
        return st.eval(T)

    def at(self, s):
        i = self.step_index(float(s.mid()))
        st = self.steps[i]
        return st.eval(s - st.s0)


def rstar(pt, r):
    """Re r*(r) = sum_i A_i log|r - r_i| (branch-independent)."""
    return sum((A * abs(r - ri).log() for A, ri in zip(pt.Acoef, pt.roots)), arb(0))
