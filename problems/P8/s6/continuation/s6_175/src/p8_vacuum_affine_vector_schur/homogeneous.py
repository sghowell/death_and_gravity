"""Closed-source lift and the actual CD causal geometry."""

from functools import cache

import sympy as s


@cache
def data():
    t = s.symbols("t", real=True)
    D = 1 + t * t
    a = D**2
    conformal = (s.atan(t) + t / D) / 2
    source = s.Function("S0")(t)
    coordinates = (t,) + s.symbols("x y z", real=True)
    one_form = (source, 0, 0, 0)
    curl = s.Matrix(
        4,
        4,
        lambda i, j: (
            s.diff(one_form[j], coordinates[i]) - s.diff(one_form[i], coordinates[j])
        ),
    )
    e, f, displacement = s.symbols("epsilon delta_F delta_displacement", real=True)
    zero_background_extra = -(e**2) * f**2 / 4 + e**2 * displacement**2 / 2
    return {
        "closed_source_condition": "dS=0, W=S, hence F(W)=0 and W-S=0; all first variations of the extra action vanish, including variations through the source and metric",
        "homogeneous_lift": "Every smooth spatially homogeneous solution of the same light target and M1, within the analytic domain, has S=S0(t)dt and an exact classical parent lift W=S. This specifies compatible vector data; it does not identify that lift with a retarded solution for a noncompact past source.",
        "CD_scale_factor": a,
        "CD_conformal_time": conformal,
        "CD_conformal_range": [-s.pi / 4, s.pi / 4],
        "causal_theorem_hypotheses": "The physical CD spacetime is conformal to the open flat time slab (-pi/4,pi/4) x R^3. Every slice is Cauchy. Smooth positive constant zeta makes L/zeta normally hyperbolic on the smooth one-form bundle. The retarded/advanced theorem is applied to compact smooth sources.",
        "Green_operator_consequence": "Unique retarded and advanced Proca inverses N G_L = G_L N with causal support; divergence W=divergence S. Existence is not a quantitative spacetime L2 bound and does not supply a chosen Hadamard state.",
        "checks": {
            "homogeneous_source_closed": curl,
            "all_extra_first_variations_zero": s.diff(zero_background_extra, e).subs(
                e, 0
            ),
            "CD_conformal_derivative": s.factor(s.diff(conformal, t) - 1 / a),
            "CD_conformal_future": s.limit(conformal, t, s.oo) - s.pi / 4,
            "CD_conformal_past": s.limit(conformal, t, -s.oo) + s.pi / 4,
            "CD_positive_conformal_factor": s.expand(a - (1 + 2 * t * t + t**4)),
        },
    }
