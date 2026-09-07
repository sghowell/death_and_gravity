"""Rational Maxwell costs with separate past and shrinking future jet caps."""

import sympy as sp
from p8a_maxwell.bounds import reference_loss_polynomial
from p8a_maxwell.domain import nonnegative, rational


def caps(values):
    if not isinstance(values, (tuple, list)) or len(values) != 4:
        raise TypeError("four exact nonnegative proper Hubble jet caps are required")
    return tuple(map(nonnegative, values))


def coefficients(past_caps, future_caps, ratio, *, beta_m):
    """Costs multiply kappa*hbar/(8*pi^2*tau^3).

Past: |H^(j)|<=d_j/tau^(j+1). Future: |H^(j)|<=c_j/(tau-t)^(j+1).
These are independent geometric assumptions, never inferred by this API.
    """
    past, future = caps(past_caps), caps(future_caps)
    ratio = nonnegative(ratio, positive=True)
    beta = rational(beta_m)
    vp = reference_loss_polynomial(*past, abs(beta))
    vf = reference_loss_polynomial(*future, abs(beta))
    qp = sp.Rational(3, 4)*past[0]**2+sp.Rational(3, 2)*past[1]
    qf = sp.Rational(3, 4)*future[0]**2+sp.Rational(3, 2)*future[1]
    # sqrt12<7/2, sqrt(6/5)<11/10, sqrt(13/35)<2/3.
    past_root = sp.Rational(7, 2)+sp.Rational(11, 5)*past[0]*ratio
    past_root += sp.Rational(2, 3)*qp*ratio**2
    past_cost = past_root**2/ratio**3+sp.Rational(13, 12600)*vp*ratio
    # sqrt(13/3)<25/12; cubic endpoint zero removes the relative poles.
    future_root = sp.Rational(7, 2)*(1+2*future[0])+sp.Rational(25, 12)*qf
    future_cost = future_root**2+sp.Rational(13, 1080)*vf
    return {"past_caps": list(past), "future_caps": list(future), "ratio": ratio,
            "beta_M": beta, "past_reference_loss": vp, "future_reference_loss": vf,
            "past_root_after_ratio_factor": past_root, "future_root": future_root,
            "past_cost": past_cost, "future_cost": future_cost,
            "total_cost": past_cost+future_cost,
            "source_weight": sp.Rational(13, 35)*(1+ratio)}


def identities():
    r, tau = sp.symbols("r tau", positive=True)
    d0, d1 = sp.symbols("d0 d1", nonnegative=True)
    v = sp.Symbol("V", nonnegative=True)
    q = sp.Rational(3, 4)*d0*d0+sp.Rational(3, 2)*d1
    dimensional_root = sp.Rational(7, 2)/(r*tau)**sp.Rational(3, 2)
    dimensional_root += sp.Rational(11, 5)*d0/(tau*sp.sqrt(r*tau))
    dimensional_root += sp.Rational(2, 3)*q*sp.sqrt(r*tau)/tau**2
    normalized_root = (sp.Rational(7, 2)+sp.Rational(11, 5)*d0*r
                       +sp.Rational(2, 3)*q*r*r)/r**sp.Rational(3, 2)
    return {"past_exact_duration_scaling": sp.simplify(tau**sp.Rational(3, 2)*dimensional_root-normalized_root),
            "past_anomaly_normalization": sp.Rational(13, 35)*r*v/360-sp.Rational(13, 12600)*r*v,
            "future_anomaly_normalization": sp.Rational(13, 3)*v/360-sp.Rational(13, 1080)*v}


def radical_margins():
    return {"sqrt12": sp.Rational(7, 2)**2-12,
            "sqrt_6_over_5": sp.Rational(11, 10)**2-sp.Rational(6, 5),
            "sqrt_13_over_35": sp.Rational(2, 3)**2-sp.Rational(13, 35),
            "sqrt_13_over_3": sp.Rational(25, 12)**2-sp.Rational(13, 3)}
