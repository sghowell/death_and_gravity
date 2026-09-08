"""Two physical readout tails after the same fourth-order subtraction."""
from functools import cache

import sympy as sp
from p8_vector_mass_remainder.tail import algebra_checks, numerator
from p8_vector_state import wkb

from . import envelopes, source, tangent

__all__ = ["algebra_checks", "reference_tail"]


def reference_tail(sector, component):
    return _reference_tail(source.kind(sector), source.output(component))


@cache
def _reference_tail(sector, component):
    data, s, readout = tangent.reference(sector), source.data(sector), source.readout(sector, component)
    adiabatic = {j: source.adiabatic(sector, component, j) for j in (0, 1, 2)}
    exact = numerator(data["S"], data["dS"], data["R"], data["dR"],
                      readout["A"], readout["delta_A"], readout["B"], readout["delta_B"],
                      s["c1"], s["delta_c1"], s["r"], adiabatic)
    low = {j: source.clean(exact.get(j, 0)) for j in range(3)}
    reference = envelopes.reference(sector)
    upper = reference["polynomial_majorants"]

    def linear(value):
        result = source.linear_bound(source.clean(value))
        if any(value != 0 for value in result["reconstructions"]):
            raise ValueError("A metric source-coefficient reconstruction failed")
        return result["absolute_upper"]

    def base(value):
        result = wkb.box_bound(sp.factor(value))
        if result["reconstruction"] != 0:
            raise ValueError("A reference-coefficient reconstruction failed")
        return result["absolute_upper"]

    majorant = numerator(upper["S"], upper["dS"], upper["R"], upper["dR"],
                         base(readout["A"]), linear(readout["delta_A"]), base(readout["B"]), linear(readout["delta_B"]),
                         base(s["c1"]), linear(s["delta_c1"]), linear(s["r"]),
                         {j: linear(value) for j, value in adiabatic.items()}, absolute=True)
    mass = wkb.MASS_TIME_MIN
    bound = sp.ceiling(sum(value/mass**(2*(j-3)) for j, value in majorant.items() if j >= 3))
    return {"low_tail_numerator_coefficients": low,
            "tail_numerator_majorants": majorant,
            "varied_reference_tail_over_inverse_frequency_fifth_upper": bound,
            "all_majorants_nonnegative": all(value >= 0 for value in majorant.values())}
