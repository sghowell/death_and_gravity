"""Exact weighted-source, trace, constraint, and actual-Wick map identities.

Primes in this module are dimensionless conformal derivatives. The physical
normalization and the unchanged quantum state are specified in FORMULATION.md.
No numerical quadrature, state-moment replacement, or order reduction is used.
"""

import sympy as sp
from p8a_residual import reconstruction as prior


def local_coefficient(scale):
    """The fixed named prescription, not an adjustable counterterm."""
    return -sp.Rational(19, 60)-sp.log(sp.sympify(scale)/2)/2


def weighted_source(scale, coefficient, time):
    """Return source components in hbar/(pi^2 A^4 eta_star^8) units.

Here scale is the dimensionless scale factor. The common physical prefactor
is hbar/(pi^2 A^4 eta_star^8); equivalently the density is that prefactor
times c/a^4. The supplied c is prescribed; pressure depends on the unknown h.
    """
    h = sp.diff(scale, time)/scale
    density = coefficient/scale**4
    pressure = density/3-sp.diff(coefficient, time)/(3*h*scale**4)
    return {"density": density, "pressure": pressure,
            "trace": sp.diff(coefficient, time)/(h*scale**4)}


def trace_rhs(scale, u, qprime, coefficient_prime, delta, time):
    """Right-hand side of q'' for the full forced dimensionless trace."""
    h = sp.diff(scale, time)/scale
    return (u/(60*delta)-(u**2/4+h**2*(u+h**2)/30)/scale**2
            +8*coefficient_prime/(h*scale**2)-2*h*qprime)


def integrated_local_density(scale, h, u, delta):
    """The nonsingular bulk term in (a^2 q')'."""
    return scale**2*u/(60*delta)-u**2/4-h**2*(u+h**2)/30


def constraint(scale, s, j, ell, coefficient, delta, time):
    """Exact dimensionless energy constraint with fixed ordinary radiation.

rho_rad=3 A^2/(kappa a_phys^4). The quantum state's own radiation
integration constant is the original zero; c is the external source.
    """
    quantum_density = sp.pi**2*scale**4*prior.components(
        scale, s, j, ell, time)["density"]
    return 3*sp.diff(scale, time)**2-3-2880*delta*(quantum_density+coefficient)


def identities():
    t = sp.Symbol("x", real=True)
    delta = sp.Symbol("delta", positive=True)
    a = sp.Function("a", positive=True)(t)
    q, c, j, ell = (sp.Function(name)(t) for name in ("q", "c", "J", "L"))
    h = sp.diff(a, t)/a
    u = -sp.diff(a, t, 2)/a
    s = a**2*q
    src = weighted_source(a, c, t)
    values = prior.components(a, s, j, ell, t)
    p_trace = sp.pi**2*a**4*(values["density"]-3*values["pressure"])
    trace_form = (-a**2*(sp.diff(q, t, 2)+2*h*sp.diff(q, t))/8
                  -u**2/32-h**2*(u+h**2)/240)
    trace_defect = -6*u/a**2-2880*delta*(p_trace/a**4+src["trace"])
    trace_equation = sp.diff(q, t, 2)-trace_rhs(a, u, sp.diff(q, t), sp.diff(c, t), delta, t)
    rules = {sp.diff(j, t): sp.diff(u, t)*s, sp.diff(ell, t): h*u**2}
    energy = constraint(a, s, j, ell, c, delta, t)
    integration = (sp.diff(c/h, t)-c*(1+u/h**2))
    bulk = integrated_local_density(a, h, u, delta)
    born_log, response, raw_log = sp.symbols("BornLog Rmode rawLog", real=True)
    named_s = response-(born_log+u*(sp.log(a/2)+sp.Rational(5, 6)))/2+u/10
    singular_s = response-born_log/2+local_coefficient(a)*u
    residual_source = sp.Symbol("Dc", real=True)
    return {
        "weighted_source_conserved_on_unknown_metric": sp.simplify(
            sp.diff(src["density"], t)+3*h*(src["density"]+src["pressure"])),
        "weighted_source_trace": sp.simplify(src["density"]-3*src["pressure"]-src["trace"]),
        "inverse_hubble_derivative": sp.simplify(sp.diff(1/h, t)-1-u/h**2),
        "source_integration_by_parts": sp.simplify(integration-sp.diff(c, t)/h),
        "actual_stress_trace_in_q": sp.simplify(p_trace-trace_form),
        "full_trace_equation_normalization": sp.simplify(trace_defect-360*delta*trace_equation/a**2),
        "full_constraint_propagation": sp.simplify(
            sp.diff(energy, t).subs(rules)-a**4*h*trace_defect),
        "integrated_auxiliary_equation": sp.simplify(
            sp.diff(a**2*sp.diff(q, t), t)-bulk-8*sp.diff(c, t)/h-a**2*trace_equation),
        "named_prescription_singular_split": sp.simplify(named_s-singular_s),
        "local_coefficient_derivative": sp.simplify(sp.diff(local_coefficient(a), t)+h/2),
        "fixed_radiation_not_adjusted_by_source": sp.diff(energy, c)+2880*delta,
        "density_defect_source_sign": (raw_log-2880*delta*residual_source).subs(
            residual_source, raw_log/(2880*delta)),
    }


def difference_identities():
    """Polynomial check of every term of the shared-history fixed-point RHS."""
    a, ab, h, hb, q, qb, v, vb, u, ub, x, upb, rp, rpb, d, db = sp.symbols(
        "a ab h hb q qb qprime qbprime u ub X ubprime Rprime Rbprime d db", real=True)
    full = a**2*(v+2*h*q)-rp-d*(upb+x)+h*u/2
    baseline = ab**2*(vb+2*hb*qb)-rpb-db*upb+hb*ub/2
    expected = (a**2*(v+2*h*q)-ab**2*(vb+2*hb*qb)-(rp-rpb)
                -d*x-(d-db)*upb+h*(u-ub)/2+(h-hb)*ub/2)
    return {"all_local_and_state_difference_terms": sp.expand(full-baseline-expected),
            "wick_derivative_to_auxiliary": sp.expand(
                a**2*(v+2*h*q)-(a**2*v+2*a*(a*h)*q))}


def negative_controls():
    t = sp.Symbol("x", real=True)
    a = sp.Function("a", positive=True)(t)
    c = sp.Function("c")(t)
    h = sp.diff(a, t)/a
    # Switching density without its compensating pressure violates conservation.
    naive_density, naive_pressure = c/a**4, c/(3*a**4)
    return {
        "naive_radiation_pressure_conservation_defect": sp.simplify(
            sp.diff(naive_density, t)+3*h*(naive_density+naive_pressure)),
        "free_radiation_past_required_c_at_a_one": -sp.Rational(1, 960),
        "trace_only_constraint_integration_constant": sp.Symbol("C_initial", real=True),
    }
