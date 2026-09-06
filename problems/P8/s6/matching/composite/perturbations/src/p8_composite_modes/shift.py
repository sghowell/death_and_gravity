"""Positive 2x2 square root and the finite relative transverse-shift stiffness."""

from functools import cache

import sympy as sp

from . import model as m


@cache
def derive():
    w, L, z = sp.symbols("w L z", real=True)
    # w=b^2(sigma_f-sigma_g)^2/Ng^2; the trace of the positive
    # t-x square root is L=sqrt((c+y)^2-w). Use L algebraically,
    # including its exact positive-branch first derivative at w=0.
    generating = (1+m.y*z)**2*(1+L*z+m.c*m.y*z**2)
    e = sp.Poly(generating, z)
    interaction = -m.M4*sum(m.betas[j]*e.nth(j) for j in range(5))
    determinant = m.alpha**2+m.alpha*m.beta*L+m.beta**2*m.c*m.y
    spatial = m.alpha**2+m.beta**2*m.y**2+2*m.alpha*m.beta*m.y*(m.c+m.y)/L
    matter = m.r**2*(m.s**2*m.n*spatial/(2*determinant)-(m.rho-m.p)*determinant/2)

    def quadratic_coefficient(value):
        return m.cancel(-sp.diff(value, L).subs(L, m.c+m.y)/(2*(m.c+m.y)))

    # Apply the literal matter chain rule before expanding its large rational
    # denominator. D_A(0)=r*s and W(0)=r² give exact cancellations in this
    # factor-preserving route; no field equation or coefficient guess is used.
    da, ww = sp.symbols("effective_D effective_W", positive=True)
    template = m.r**2*(m.s**2*m.n*ww/(2*da)-(m.rho-m.p)*da/2)
    dprime = quadratic_coefficient(determinant)
    wprime = quadratic_coefficient(spatial)
    matter_literal = m.cancel(sp.diff(template, da).subs({da: m.r*m.s, ww: m.r**2})*dprime
                             +sp.diff(template, ww).subs({da: m.r*m.s, ww: m.r**2})*wprime)
    literal = quadratic_coefficient(interaction)+matter_literal
    zv = m.P+m.alpha*m.beta*m.r**2*m.rho+m.alpha*m.beta*m.r*m.s*m.y*m.n/(m.c+m.y)
    xi = 2*m.y**2*zv/(m.c+m.y)
    return {"w": w, "L": L, "generating": generating,
            "interaction": interaction, "matter": matter,
            "D_A": determinant, "W": spatial,
            "D_w_initial": dprime, "W_w_initial": wprime,
            "literal_matter_quadratic": matter_literal,
            "literal_quadratic": literal, "potential_quadratic": quadratic_coefficient(interaction),
            "Z_v": zv, "Xi": xi,
            "sqrt_trace": sp.sqrt((m.c+m.y)**2-w)}


def checks():
    d = derive()
    v = sp.Symbol("v", real=True)  # v=b(sigma_f-sigma_g)/Ng, w=v^2
    block = sp.Matrix([[m.c**2-v**2, -v*m.y], [v*m.y, m.y**2]])
    # This is the coordinate block conjugated by diag(Ng,a): no
    # spectrum, determinant or positive-root trace is changed.
    trace_squared = (m.c+m.y)**2-v**2
    numerator = block+m.c*m.y*sp.eye(2)
    square_residual = numerator*numerator-trace_squared*block
    affine = m.alpha*sp.eye(2)+m.beta*numerator/d["L"]
    effective = sp.diag(1, -1)*(m.alpha**2*sp.eye(2)
                                +2*m.alpha*m.beta*numerator/d["L"]+m.beta**2*block)
    # Reduce only by the already checked positive-root equation. The
    # denominator L is nonzero in the local shift chart.
    determinant_residual = sp.fraction(sp.cancel(affine.det()-d["D_A"]))[0]
    determinant_remainder = sp.rem(determinant_residual, v**2+d["L"]**2-(m.c+m.y)**2, v)
    at_zero = {d["L"]: m.c+m.y}
    result = {"relative_block_determinant": m.cancel(block.det()-m.c**2*m.y**2),
              "positive_trace_at_zero": sp.sqrt((m.c+m.y)**2)-(m.c+m.y),
              "effective_determinant_at_zero": m.cancel(d["D_A"].subs(at_zero)-m.r*m.s),
              "effective_inverse_lapse_at_zero": m.cancel(d["W"].subs(at_zero)/d["D_A"].subs(at_zero)**2-1/m.s**2),
              "literal_affine_effective_determinant": m.cancel(determinant_remainder),
              "literal_effective_inverse_numerator": m.cancel(-effective[1, 1]-d["W"]),
              "literal_D_w_at_zero": m.cancel(d["D_w_initial"]+m.alpha*m.beta/(2*(m.c+m.y))),
              "literal_W_w_at_zero": m.cancel(d["W_w_initial"]-m.alpha*m.beta*m.y/(m.c+m.y)**2),
              "literal_potential_shift": m.cancel(d["potential_quadratic"]-m.P/(2*(m.c+m.y))),
              "literal_matter_plus_potential_shift": m.cancel(d["literal_quadratic"]-d["Z_v"]/(2*(m.c+m.y))),
              "finite_Xi_no_division_by_mu": m.cancel(d["Xi"]-4*m.y**2*d["literal_quadratic"])}
    result.update({f"square_root_matrix_{i}_{j}": m.cancel(square_residual[i, j])
                   for i in range(2) for j in range(2)})
    return result


def positive_factors():
    """Sufficient domain, not a positivity claim for every beta_n or source."""
    pos_p, pos_rho, pos_null = sp.symbols("positive_P positive_rho positive_null", positive=True)
    zv = pos_p+m.alpha*m.beta*m.r**2*pos_rho+m.alpha*m.beta*m.r*m.s*m.y*pos_null/(m.c+m.y)
    return {"positive_Z_v": zv, "positive_Xi": 2*m.y**2*zv/(m.c+m.y)}


def negative_controls():
    d = derive()
    return {"omit_shift_matter": m.cancel(d["literal_quadratic"]-d["potential_quadratic"]),
            "omit_inverse_effective_metric": m.cancel(d["W"].subs(d["L"], m.c+m.y)-m.s**2)}
