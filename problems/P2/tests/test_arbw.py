"""Validated Wronskian (arbw.py): enclosures contain the float values; exact symmetry; D4 check.

Fast (~3 s).
"""
import mpmath as mp
import pytest
from flint import acb, arb

from p2.arbutil import precision
from p2.arbw import wronskian
from p2.params import Point
from p2.qnm_float import Geom, Solver


def _acb(z):
    return acb(arb(mp.nstr(mp.re(z), 60)), arb(mp.nstr(mp.im(z), 60)))


@pytest.mark.parametrize("name,l,z", [("PA", 1, ("3.0", "-0.5")), ("PB", 0, ("1.0", "-0.4")),
                                       ("PB", 2, ("12.0", "-0.2"))])
def test_enclosure_contains_float(name, l, z):
    mp.mp.dps = 50
    pt, g = Point(name), Geom.frozen(name)
    w = mp.mpc(*z) * g.km
    with precision(256):
        W = wronskian(pt, l, _acb(w))
        Wf = Solver(g).wronskian(w, l)
        assert W.real.rad() < arb("1e-50") * abs(W).mid()
        assert abs(W - _acb(Wf)) < arb("1e-40") * abs(W).mid()


def test_symmetry_minus_conj():
    """W_l(-conj w) = conj W_l(w) (used by P2-2 for Re w < 0)."""
    pt = Point("PA")
    with precision(256):
        w = acb("0.31", "-0.04")
        W1 = wronskian(pt, 3, w)
        W2 = wronskian(pt, 3, -w.conjugate())
        assert (W2 - W1.conjugate()).real.contains(0) and (W2 - W1.conjugate()).imag.contains(0)
        assert abs(W2 - W1.conjugate()) < arb("1e-50") * abs(W1).mid()


@pytest.mark.parametrize("name", ["PB", "PA"])
def test_D4_derivative_matches_exact_formula(name):
    """C8: W_0'(0) = i (r_+^2 + r_c^2)/(r_+ r_c); central difference of validated values."""
    pt = Point(name)
    with precision(256):
        h = arb("1e-20")
        fd = (wronskian(pt, 0, acb(h)) - wronskian(pt, 0, acb(-h))) / (2 * h)
        exact = acb(0, 1) * (pt.rp ** 2 + pt.rc ** 2) / (pt.rp * pt.rc)
        assert abs(wronskian(pt, 0, acb(0))) < arb("1e-60")
        assert abs(fd - exact) < arb("1e-30")
