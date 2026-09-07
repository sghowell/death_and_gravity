"""Validated physical controllability Gramian on the punctured interval.

No floating-point trajectory is used as evidence. Arb Taylor coefficients
and interval Taylor remainders enclose the exact matrix ODE at K=5/2.
The continuous K extension is the separate energy argument in notes/gramian.md.
"""
from fractions import Fraction
from functools import cache

from flint import arb, arb_mat, arb_series, ctx

ORDER = 18
STEPS = 64
PRECISION = 256


def _integer(value, name, lower, upper):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} requires an exact Python integer")
    if not lower <= value <= upper:
        raise ValueError(f"{name} is outside the audited range")
    return value


def _fraction(value):
    """Convert an exact finite Arb endpoint to a rational, without decimals."""
    mantissa, exponent = value.man_exp()
    power = int(exponent)
    return Fraction(int(mantissa)) * Fraction(2)**power


def _round_up(value, places=24):
    value = _fraction(value.upper())
    scale = 10**places
    return Fraction(-(-value.numerator*scale//value.denominator), scale)


def _round_down(value, places=24):
    value = _fraction(value.lower())
    scale = 10**places
    return Fraction(value.numerator*scale//value.denominator, scale)


def _matrix_pairs(matrix, places=24):
    return [[(_round_down(matrix[i, j], places), _round_up(matrix[i, j], places))
             for j in range(4)] for i in range(4)]


def _generator(center, order):
    """Exact Taylor arithmetic for N in X=(g,ell*g_u,f,ell*f_u), x=(u-a)/ell."""
    x = arb_series([center, 1], prec=order+2)
    u = (x-2)/100
    d = 1+u*u
    # d^4-1 factored without cancellation; denominators exclude zero on J.
    dm = u*u*(4+6*u*u+4*u**4+u**6)
    ug = arb(32)/10000*(1-u*u)/(d**14*dm)
    uf = arb(8)/10000*(1-u*u)/(d**2*dm)
    p = arb(12)/100*u/d
    gg, gf = arb(5)/20000/d**4, arb(5)/20000*d**4
    one = arb_series([1], prec=order+2)
    zero = arb_series([], prec=order+2)
    coefficients = {(0, 1): one, (1, 0): -gg-ug, (1, 1): -p, (1, 2): ug,
                    (2, 3): one, (3, 0): uf, (3, 2): -gf-uf, (3, 3): p}
    return [arb_mat(4, 4, [coefficients.get((i, j), zero)[n]
                          for i in range(4) for j in range(4)])
            for n in range(order+1)]


def _infnorm(matrix):
    return max(sum((matrix[i, j].abs_upper() for j in range(4)), arb(0)).upper()
               for i in range(4))


def _taylor(generator, initial, order):
    """Differentiate G'=NG+GN^T+e2e2^T; all coefficients include 1/n!."""
    zero = arb_mat(4, 4)
    forcing = arb_mat([[0, 0, 0, 0], [0, 1, 0, 0],
                      [0, 0, 0, 0], [0, 0, 0, 0]])
    coefficients = [initial]
    for n in range(order):
        total = forcing if n == 0 else zero
        for j in range(n+1):
            total = (total+generator[j]*coefficients[n-j]
                     +coefficients[n-j]*generator[j].transpose())
        coefficients.append(total/(n+1))
    return coefficients


def _integrate(order, steps):
    h = arb(1)/steps
    value = arb_mat(4, 4)
    max_lambda, max_remainder = arb(0), arb(0)
    for index in range(steps):
        left = arb(index)/steps
        generator = _generator(left, order)
        # This closed ball contains the entire step, not only its midpoint.
        range_generator = _generator(left+h/2+arb(0, h/2), order)
        lam = (_infnorm(range_generator[0])
               +_infnorm(range_generator[0].transpose())).upper()
        if not h*lam < 1:
            raise ValueError("The local matrix-ODE tube condition failed")
        # Integral inequality on the whole step supplies G(xi) for the
        # Lagrange remainder. This is not a guessed trajectory enclosure.
        tube = ((_infnorm(value)+h)/(1-h*lam)).upper()
        range_state = arb_mat(4, 4, [arb(0, tube)]*16)
        remainder = _taylor(range_generator, range_state, order+1)[-1]
        coefficients = _taylor(generator, value, order)
        polynomial = coefficients[-1]
        for n in range(order-1, -1, -1):
            polynomial = polynomial*h+coefficients[n]
        errors = [(remainder[i, j].abs_upper()*h**(order+1)).upper()
                  for i in range(4) for j in range(4)]
        value = polynomial+arb_mat(4, 4, [arb(0, error) for error in errors])
        max_lambda = max(max_lambda, lam)
        max_remainder = max(max_remainder, *errors)
    return value, max_lambda, max_remainder


def _endpoint_transform():
    """Scaled physical X -> scaled outer Z, including all moving derivatives."""
    u = arb_series([arb(-1)/100, 1], prec=2)
    d = 1+u*u
    denominator = 2*d**12+8
    fs = denominator.sqrt()/(2*arb(2).sqrt()*d**3)
    w1, w2 = 2*d**12/denominator, 8/denominator
    fr = arb(2).sqrt()*d**3/denominator.sqrt()
    a, b = fs*w1, fs*w2
    r, mu = arb(1)/100, arb(39).sqrt()/2
    return arb_mat([
        [a[0], 0, b[0], 0],
        [a[1]/100, a[0], b[1]/100, b[0]],
        [-fr[0]/r.sqrt(), 0, fr[0]/r.sqrt(), 0],
        [(r*fr[1]+fr[0]/2)/(mu*r.sqrt()), fr[0]/(mu*r.sqrt()),
         (-r*fr[1]-fr[0]/2)/(mu*r.sqrt()), -fr[0]/(mu*r.sqrt())]])


def _positive_pivots(matrix):
    """Interval LDL on the symmetric exact matrix enclosed by matrix."""
    working = arb_mat(matrix.tolist())
    pivots = []
    for j in range(4):
        pivot = working[j, j]
        if not pivot > 0:
            raise ValueError("A strict positive Gramian pivot was not certified")
        pivots.append(pivot)
        for i in range(j+1, 4):
            for k in range(j+1, 4):
                working[i, k] = working[i, k]-working[i, j]*working[j, k]/pivot
    return pivots


def enclosure(*, order=ORDER, steps=STEPS, precision=PRECISION):
    order = _integer(order, "order", 12, 24)
    steps = _integer(steps, "steps", 64, 256)
    precision = _integer(precision, "precision", 128, 512)
    old_precision, old_cap = ctx.prec, ctx.cap
    try:
        ctx.prec, ctx.cap = precision, order+2
        gramian, lam, remainder = _integrate(order, steps)
        transform = _endpoint_transform()
        # L2(du) normalization: B_x=2ell^2 e2, du=ell dx.
        outer = arb(4)/1000000*transform*gramian*transform.transpose()
        shifted = outer-arb_mat([[arb(36)/10**12 if i == j else 0
                                 for j in range(4)] for i in range(4)])
        pivots = _positive_pivots(shifted)
        if not lam < 26 or not remainder < arb(1)/10**18:
            raise ValueError("The validated step constants exceeded their caps")
        return {
            "order": order, "steps": steps, "precision_bits": precision,
            "midpoint_K": Fraction(5, 2),
            "closed_x_interval": (Fraction(0), Fraction(1)),
            "step_matrix_norm_upper": _round_up(lam),
            "local_entry_remainder_upper": _round_up(remainder, 40),
            "normalized_physical_gramian": _matrix_pairs(gramian),
            "scaled_outer_gramian": _matrix_pairs(outer, 32),
            "endpoint_transform": _matrix_pairs(transform),
            "shifted_LDL_pivot_intervals": [
                (_round_down(pivot, 36), _round_up(pivot, 36)) for pivot in pivots],
            "midpoint_scaled_outer_lower": Fraction(36, 10**12),
            "full_L2_du_scaling": Fraction(4, 10**6),
        }
    finally:
        ctx.prec, ctx.cap = old_precision, old_cap


@cache
def certificate():
    return enclosure()
