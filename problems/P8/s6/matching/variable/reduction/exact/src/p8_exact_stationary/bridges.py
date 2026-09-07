"""Compare separately derived actions and the unchanged frozen profiles."""
from functools import cache

import sympy as sp
from p8_variable_reduction import stationary as ancestor

from . import audit, stationary


@cache
def checks():
    a, e, t = audit.action(), stationary.own_equations(), stationary.tensors()
    mapping = {audit.T: stationary.TIME, audit.M2: stationary.M**2,
               a["b"]: e["b"], a["N"]: e["N"], a["p"]: e["beta1"], a["w"]: e["beta4"],
               a["q"]: t["q"], a["h"]: t["Q_tensor"]}

    def convert(value):
        return value.xreplace(mapping).doit()

    values = {"literal_coordinate_isotropic_action": convert(a["isotropic"])-e["density"],
              "literal_coordinate_lapse_equation": convert(a["C"])-e["C"],
              "literal_coordinate_spatial_equation": convert(a["E"])-e["E"],
              "literal_coordinate_EH_boundary": convert(a["boundary"])-e["EH_boundary"],
              "literal_coordinate_tensor_action": convert(a["quadratic"])-t["density"],
              "literal_coordinate_physical_EH": convert(a["g_EH"])-stationary.M**2*sp.diff(t["q"], stationary.TIME)**2/4,
              "literal_coordinate_tensor_g_Euler": convert(a["tensor_equations"]["g"])-t["g_equation"],
              "literal_coordinate_tensor_f_Euler": convert(a["tensor_equations"]["f"])-t["f_equation"],
              "literal_coordinate_full_EH_to_ADM": convert(a["f_EH"]-a["f_ADM"]-sp.diff(a["boundary"], audit.T)),
              "literal_coordinate_full_potential": convert(a["potential"])-t["literal_potential"].rewrite(sp.exp)}
    p, old = stationary.profile(), ancestor.profile()
    variable_map = {old["v"]: stationary.V, old["c"]: stationary.C}
    for current, prior in (("B1", "b1"), ("B4", "b4"), ("algebraic_r_cubed", "r_cubed")):
        values["unchanged_frozen_"+current] = p[current]-old[prior].subs(variable_map)
    # Clear the cube-root relation algebraically, not through a numerical root.
    j = stationary.joint()
    # Keep the common positive cube root as a single algebraic atom while
    # cancelling a rational identity; expanding it first is unnecessary.
    atom = sp.Symbol("positive_cube_root_atom", positive=True)
    values["whole_algebraic_factor"] = (
        3*atom*stationary.V*j["F"]**2-2*p["B1"]*(1-p["R"]*j["b_cubed"])
        -2*stationary.V*p["Q"]*(3*atom*j["F"]**2/(2*p["Q"])-stationary.ZETA))
    # Each expression is an exact identity on its stated nonzero-denominator
    # domain. Clear those denominators before polynomial expansion; expanding
    # separate high-degree rational terms first wastes a large symbolic GCD.
    values = {name: sp.factor(sp.expand(sp.together(value).as_numer_denom()[0]))
              for name, value in values.items()}
    if any(value != 0 for value in values.values()):
        raise ValueError("An independent coordinate or frozen-profile bridge failed")
    return values
