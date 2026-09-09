"""Fresh one-real-conformal-scalar reference, not a photon stress substitution."""

from functools import cache

import sympy as sp

H0, H1, H2, H3 = sp.symbols("H Hdot Hddot Hthird", real=True)
BETA_S = sp.Symbol("beta_S", real=True)


def reference_jets(h=H0, hd=H1, hdd=H2, h3=H3, beta_s=BETA_S):
    h, hd, hdd, h3, beta = map(sp.sympify, (h, hd, hdd, h3, beta_s))
    irho = 18 * hd**2 - 108 * h * h * hd - 36 * h * hdd
    itr = -36 * (h3 + 4 * hd * hd + 7 * h * hdd + 12 * h * h * hd)
    rho = 3 * h**4 + beta * irho
    trace = 12 * h * h * (h * h + hd) + beta * itr
    return {
        "rho": sp.expand(rho),
        "trace": sp.expand(trace),
        "pressure": sp.expand((rho - trace) / 3),
        "EED": sp.expand(rho - trace / 2),
        "I_rho": irho,
        "I_trace": itr,
        "I_EED": sp.expand(irho - itr / 2),
    }


@cache
def checks():
    t = sp.Symbol("time", real=True)
    h = sp.Function("H")(t)
    hd, hdd, h3 = [sp.diff(h, t, n) for n in (1, 2, 3)]
    data = reference_jets(h, hd, hdd, h3)
    r = 6 * (hd + 2 * h * h)
    r00 = 3 * (hd + h * h)
    rii = -(hd + 3 * h * h)
    ricci2 = r00 * r00 + 3 * rii * rii
    h3rho = r00 * r00 - sp.Rational(2, 3) * r * r00 - ricci2 / 2 + r * r / 4
    h3trace = -ricci2 + r * r / 3
    boxr = sp.diff(r, t, 2) + 3 * h * sp.diff(r, t)
    # HH has the simultaneous opposite metric/Ricci tensor, so H3_00 flips.
    hh_scalar_rho_pi2 = -2 * sp.Rational(1, 360) * (-h3rho) / 16
    return {
        "fresh_scalar_anomaly_density_contraction": sp.factor(h3rho - 3 * h**4),
        "fresh_scalar_anomaly_trace_contraction": sp.factor(
            h3trace - 12 * h * h * (h * h + hd)
        ),
        "one_real_scalar_HH_to_FK_physical_density": sp.factor(
            hh_scalar_rho_pi2 - h**4 / 960
        ),
        "new_beta_is_coefficient_of_same_covariant_I_variation": sp.factor(
            data["I_rho"] - (2 * r * r00 - r * r / 2 - 6 * h * sp.diff(r, t))
        ),
        "finite_beta_scalar_reference_trace": sp.factor(data["I_trace"] + 6 * boxr),
        "full_physical_scalar_reference_conservation": sp.factor(
            sp.diff(data["rho"], t) + 3 * h * (data["rho"] + data["pressure"])
        ),
        "conformal_reference_EED_keeps_anomaly_trace": sp.factor(
            data["EED"] - data["rho"] + data["trace"] / 2
        ),
        "scalar_Euler_coefficient_is_one_half_not_photon_thirty_one": sp.factor(
            data["trace"].subs(BETA_S, 0)
            - sp.Rational(1, 2) * 24 * h * h * (h * h + hd)
        ),
        "radiation_scalar_finite_beta_vanishes": sp.factor(
            data["I_EED"].subs(h, 1 / (2 * t)).doit()
        ),
        "radiation_actual_scalar_reference_EED": sp.factor(
            data["EED"].subs(h, 1 / (2 * t)).doit() - sp.Rational(9, 16) / t**4
        ),
    }


def loss_polynomial(caps, beta_abs):
    h, hd, hdd, h3 = caps
    return 3 * (h**4 + 2 * h * h * hd) + beta_abs * (
        18 * h3 + 90 * hd * hd + 90 * h * hdd + 108 * h * h * hd
    )
