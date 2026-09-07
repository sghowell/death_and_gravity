"""Independent canonical action, source weight and asinh endpoint map."""
from functools import cache

import sympy as sp

U = sp.Symbol("proper_time_over_tau", real=True)
DELTA, TAU = sp.symbols("delta tau", positive=True)
RHO = sp.sqrt(7)/2


@cache
def canonical():
    k, spring, q, hidden, y = (sp.Function(name)(U) for name in ("kinetic", "spring", "prescribed_q", "hidden_Q", "canonical_Y"))
    root = sp.sqrt(k)
    pump = sp.diff(root, U, 2)/root
    potential = spring/k-pump
    original = sp.diff(q, U)**2/4+k*sp.diff(hidden, U)**2/4-spring*(hidden-q)**2/4
    replaced = original.subs(hidden, y/root).doit()
    normal = sp.diff(q, U)**2/4+sp.diff(y, U)**2/4-potential*y**2/4+spring*y*q/(2*root)-spring*q**2/4
    boundary = -sp.diff(root, U)*y**2/(4*root)
    original_euler = sp.diff(k*sp.diff(hidden, U), U)+spring*(hidden-q)
    normalized_euler = sp.diff(y, U, 2)+potential*y-spring*q/root
    return {"k": k, "spring": spring, "q": q, "Q": hidden, "Y": y,
            "sqrt_k": root, "pump": pump, "V": potential,
            "original_density": original, "substituted_density": replaced,
            "canonical_density": normal, "boundary": boundary,
            "original_own_Euler": original_euler, "canonical_Euler": normalized_euler,
            "canonical_source_weight": spring/root}


@cache
def asinh_map():
    k, ku = sp.symbols("kinetic_value kinetic_u", real=True)
    # Positivity of k is a separate interval-proved domain premise.
    k = sp.Symbol("kinetic_value", positive=True)
    radial = sp.sqrt(U**2+DELTA/8)
    t = sp.asinh(sp.sqrt(8)*U/sp.sqrt(DELTA))
    psi, psit, psitt = sp.symbols("psi psi_t psi_tt", real=True)
    remainder, source = sp.symbols("potential_remainder canonical_source", real=True)

    def derivative(expression):
        return (sp.diff(expression, U)+sp.diff(expression, psi)*psit/radial
                +sp.diff(expression, psit)*psitt/radial)

    y = sp.sqrt(radial)*psi
    yu, yuu = derivative(y), derivative(derivative(y))
    transformed = radial**sp.Rational(3, 2)*(yuu+(16/(DELTA+8*U**2)+remainder)*y-source)
    reference = RHO**2+3*DELTA/(32*radial**2)
    expected = psitt+(reference+radial**2*remainder)*psi-radial**sp.Rational(3, 2)*source
    endpoint = sp.Matrix([[sp.sqrt(k/radial), 0],
                          [sp.sqrt(k*radial)*(ku/(2*k)-U/(2*radial**2))/RHO,
                           sp.sqrt(k*radial)/RHO]])
    # Q_u=tau*Q_T; both maps are explicit so a physical-time column
    # is not confused with a dimensionless-time column.
    physical_endpoint = endpoint*sp.diag(1, TAU)
    return {"r": radial, "t": t, "rho": RHO, "k": k, "k_u": ku,
            "psi": psi, "psi_t": psit, "psi_tt": psitt,
            "remainder": remainder, "source": source,
            "Y": y, "Y_u": yu, "transformed": transformed, "expected": expected,
            "reference_potential_in_u": reference,
            "source_weight_after_asinh": radial**sp.Rational(3, 2),
            "endpoint_from_Q_Qu": endpoint, "endpoint_from_Q_QT": physical_endpoint}


@cache
def identities():
    c, a = canonical(), asinh_map()
    values = {
        "canonical_action_boundary": c["substituted_density"]-c["canonical_density"]-sp.diff(c["boundary"], U),
        "canonical_equation_and_source": c["original_own_Euler"].subs(c["Q"], c["Y"]/c["sqrt_k"]).doit()/c["sqrt_k"]-c["canonical_Euler"],
        "canonical_action_Euler": -2*(sp.diff(c["canonical_density"], c["Y"])-sp.diff(sp.diff(c["canonical_density"], sp.diff(c["Y"], U)), U))-c["canonical_Euler"],
        "asinh_clock_derivative": sp.diff(a["t"], U)-1/a["r"],
        "asinh_full_potential_and_source": a["transformed"]-a["expected"],
        "reference_frequency": RHO**2-sp.Rational(7, 4),
        "radial_to_sech": a["reference_potential_in_u"].subs(U, sp.sqrt(DELTA/8)*sp.sinh(sp.Symbol("t", real=True)))
                         -RHO**2-sp.Rational(3, 4)/sp.cosh(sp.Symbol("t", real=True))**2,
        "dimensionless_endpoint_determinant": a["endpoint_from_Q_Qu"].det()-a["k"]/RHO,
        "physical_endpoint_determinant": a["endpoint_from_Q_QT"].det()-TAU*a["k"]/RHO,
    }
    return {name: sp.factor(sp.simplify(value)) for name, value in values.items()}


def calibration():
    c, a = canonical(), asinh_map()
    return {"canonical_field": "Y=sqrt(k)*Q", "canonical_source_weight": c["canonical_source_weight"],
            "asinh_source": a["source_weight_after_asinh"]*c["canonical_source_weight"]*c["q"],
            "explicit_u_boundary": c["boundary"],
            "wave_state": "(psi,psi_t/rho); Y=sqrt(r)*psi; r=sqrt(u²+delta/8)",
            "endpoint_Q_Qu": a["endpoint_from_Q_Qu"],
            "endpoint_Q_QT": a["endpoint_from_Q_QT"],
            "physical_metric_not_redefined": True,
            "canonical_source_bound_or_zero_data_response_proved_here": False}
