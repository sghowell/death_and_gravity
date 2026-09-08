"""Rational continuous bounds for the complete source, not an inverse bound."""
from functools import cache

import sympy as sp
from p8_affine_retuned import bounds as exact_parent

RHO = sp.Rational(1, 25)


@cache
def constants():
    n, rho = 1-RHO, RHO
    R = (sp.Rational(5, 2)+sp.Rational(25, 18)*sp.Rational(45, 31))/3
    return {"Q_remainder": R, "Qx_remainder": 3*R,
            "distance_ratio": (2+rho)/n**2,
            "normal_r2k": (3+2*rho)/n**2,
            "normal_geometry": 6*(5+6*rho+2*rho**2)/n**3,
            "normal_Q_geometry": sp.Rational(27, 16)*(16+39*rho+40*rho**2+20*rho**3+4*rho**4)/n**5,
            "normal_Q_remainder": sp.Rational(3, 2)*R/n*sp.Rational(9, 4)**3,
            "coordinate_r2k": 1/n,
            "coordinate_geometry": 6*(3+2*rho)/n**2,
            "coordinate_Q_geometry": sp.Rational(27, 16)*(12+23*rho+16*rho**2+4*rho**3)/n**4,
            "coordinate_Q_remainder": sp.Rational(3, 2)*R*sp.Rational(9, 4)**3,
            "derivative_geometry": 6*(9+11*rho+4*rho**2)/n**3,
            "derivative_Q_geometry": sp.Rational(27, 8)*(18+40*rho+40*rho**2+20*rho**3+4*rho**4)/n**5,
            "derivative_Q_remainder": 3/n**3*(3*R)*sp.Rational(9, 4)**2,
            "electric_error_constant": (230+sp.Rational(83, 2))/n}


@cache
def proof_checks():
    c = constants()
    controls = {"normal_r2k": 4, "normal_geometry": 36, "normal_Q_geometry": 37,
                "normal_Q_remainder": 27, "coordinate_r2k": 2, "coordinate_geometry": 21,
                "coordinate_Q_geometry": 26, "coordinate_Q_remainder": 26,
                "derivative_geometry": 65, "derivative_Q_geometry": 82, "derivative_Q_remainder": 78}
    out = {name+"_continuous_upper": bool(c[name] < upper) for name, upper in controls.items()}
    out.update({"strictly_inside_original_X_tube": bool((1+RHO)**-2 > sp.Rational(9, 10)
                                                         and (1-RHO)**-2 < sp.Rational(11, 10)),
                "original_Q_cubic_remainder": c["Q_remainder"] == sp.Rational(140, 93),
                "original_Qx_quadratic_remainder": c["Qx_remainder"] == sp.Rational(140, 31),
                "distance_ratio_below_nine_fourths": bool(c["distance_ratio"] < sp.Rational(9, 4)),
                "normal_source_error_constant": 4+36+37+27 == 104,
                "coordinate_source_error_constant": 2+21+26+26 == 75,
                "coordinate_r_derivative_constant": 65+82+78 == 225,
                "coordinate_gradient_error_constant": 3+225+2 == 230,
                "quadratic_normal_bound": 2+12+sp.Rational(27, 4) == sp.Rational(83, 4),
                "full_normal_bound_below_25": bool(sp.Rational(83, 4)+104*RHO < 25),
                "electric_error_below_300": bool(c["electric_error_constant"] < 300),
                "full_electric_bound_below_54": bool(sp.Rational(83, 2)+300*RHO < 54),
                "local_Maxwell_density_error_constant": sp.Rational(300, 2)*(54+sp.Rational(83, 2)) == 14325,
                "local_density_is_not_the_nonlocal_induced_action": True})
    return out


def values(epsilon, zeta):
    epsilon = exact_parent.exact(epsilon, "epsilon")
    zeta = exact_parent.exact(zeta, "zeta")
    if epsilon.is_positive is not True or epsilon > RHO or zeta.is_positive is not True:
        raise ValueError("Require exact 0<epsilon<=1/25 and zeta>0")
    return {"normal_source_remainder": 104*epsilon**3,
            "coordinate_source_remainder": 75*epsilon**3,
            "coordinate_spatial_gradient_remainder": 230*epsilon**3,
            "physical_electric_source_remainder": 300*epsilon**3,
            "full_normal_source": 25*epsilon**2, "full_electric_source": 54*epsilon**2,
            "local_first_derivative_Maxwell_density_remainder": 14325*zeta*epsilon**5,
            "nonlocal_inverse_or_nonlinear_solution_bound": False}


def physical_values(epsilon, zeta_physical, mass_squared=1, time_scale=1):
    zeta_physical, mass_squared, time_scale = [exact_parent.exact(value, name) for value, name in
                                               ((zeta_physical, "zeta_physical"),
                                                (mass_squared, "mass_squared"), (time_scale, "time_scale"))]
    if any(value.is_positive is not True for value in (zeta_physical, mass_squared, time_scale)):
        raise ValueError("Physical curl, mass squared and time scale must be positive")
    normalized = values(epsilon, zeta_physical/(mass_squared*time_scale**2))
    return {"normalized_zeta": zeta_physical/(mass_squared*time_scale**2),
            "normal_source_remainder": normalized["normal_source_remainder"]/time_scale,
            "electric_source_remainder": normalized["physical_electric_source_remainder"]/time_scale**2,
            "local_Maxwell_density_remainder": mass_squared/time_scale**2
                *normalized["local_first_derivative_Maxwell_density_remainder"],
            "nonlocal_inverse_or_nonlinear_solution_bound": False}
