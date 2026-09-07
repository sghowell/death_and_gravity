"""Independent source-loading constants and physical scalar phase readout."""
from fractions import Fraction
from functools import cache

import sympy as sp


def constants():
    f = Fraction
    length, a, delta = f(1, 100), f(1, 200), f(1, 10**6)
    potential_upper = 8/length**2+44
    green_derivative = 1-potential_upper*(length/4)**2/2
    green_ratio = 1-potential_upper*(length/4)**2/6
    area, lag = length/16, length/8
    central = 44*(a**2+2*a*f(1, 2000))
    tails = 3*delta/(32*a**2)
    budget = (central+tails)*f(4, 5)
    error = 4*(1/(1-f(1, 400))-1)
    q_gap = (2-2*error)*f(9, 1280)*f(2, 5)
    return {"L": length, "a": a, "delta_max": delta,
            "canonical_potential_upper": potential_upper,
            "Green_derivative_floor": green_derivative,
            "Green_ratio_floor": green_ratio,
            "sharp_loaded_Y_floor_per_eta": green_ratio*lag*area/length**2,
            "sharp_loaded_Y_u_floor_per_eta": green_derivative*area/length**2,
            "central": central, "tails": tails, "budget": budget,
            "transfer_error_upper": error,
            "balanced_pair_separation_factor": 2-2*error,
            "Q_pair_separation_per_eta": q_gap,
            "Q_gap_above_one_over_200": q_gap-f(1, 200),
            "strict_kinetic_upper": f(21, 10)**3/f(19, 10),
            "strict_kinetic_lower": f(19, 10)**3/f(21, 10)}


@cache
def phase_algebra():
    ar, ai, beta, cc, ss, w1, w2 = sp.symbols("ar ai beta cos_phase sin_phase w1 w2", real=True)
    aa = ar+sp.I*ai
    phase = cc+sp.I*ss
    waves = sp.Matrix([[1, 1], [sp.I, -sp.I]])/sp.sqrt(2)
    dressed = sp.Matrix([[sp.conjugate(aa)*phase, sp.I*beta],
                         [-sp.I*beta, aa*sp.conjugate(phase)]])
    real_map = sp.simplify(waves*dressed*waves.H)
    rr, ii = ar*cc+ai*ss, ar*ss-ai*cc
    expected = sp.Matrix([[rr, ii-beta], [-ii-beta, rr]])
    state = sp.Matrix([w1, w2])
    oscillatory = sp.Matrix([[rr, ii], [-ii, rr]])
    scalar_cos, scalar_sin = ar*w1-ai*w2, ai*w1+ar*w2
    flipped = real_map.subs({cc: -cc, ss: -ss}, simultaneous=True)
    identities = {
        "real_time_frequency_map": real_map-expected,
        "constant_mixing_cancels_in_phase_pair": real_map-flipped-2*oscillatory,
        "oscillatory_norm": sp.expand((oscillatory*state).dot(oscillatory*state)
            -(ar**2+ai**2)*(cc**2+ss**2)*(w1**2+w2**2)),
        "scalar_first_component": sp.expand((real_map*state)[0]
            -(cc*scalar_cos+ss*scalar_sin-beta*w2)),
        "scalar_phase_amplitude_squared": sp.expand(scalar_cos**2+scalar_sin**2
            -(ar**2+ai**2)*(w1**2+w2**2)),
    }
    flattened = {}
    for name, value in identities.items():
        if isinstance(value, sp.MatrixBase):
            for row in range(value.rows):
                for col in range(value.cols):
                    flattened[f"{name}_{row}{col}"] = sp.simplify(value[row, col])
        else:
            flattened[name] = sp.simplify(value)
    # Q is the original amplitude, not an arbitrary coordinate projection.
    radial, kinetic = sp.symbols("r k", positive=True)
    kinetic_u, u, rho = sp.symbols("k_u u rho", real=True, nonzero=True)
    actual_map = sp.Matrix([[sp.sqrt(kinetic/radial), 0],
                           [sp.sqrt(kinetic*radial)/rho*(kinetic_u/(2*kinetic)-u/(2*radial**2)),
                            sp.sqrt(kinetic*radial)/rho]])
    first_row = actual_map.inv()[0, :]-sp.Matrix([[sp.sqrt(radial/kinetic), 0]])
    for col in range(2):
        flattened[f"original_Q_readout_{col}"] = sp.simplify(first_row[col])
    return {"identities": flattened, "real_map": real_map,
            "scalar_cos": scalar_cos, "scalar_sin": scalar_sin,
            "state": state, "oscillatory_part": oscillatory}


def checks():
    c, f = constants(), Fraction
    return {
        "Green_positive_flux": c["Green_derivative_floor"] > f(2, 3),
        "Green_positive_ratio": c["Green_ratio_floor"] > f(9, 10),
        "strict_limit_loading_Y": c["sharp_loaded_Y_floor_per_eta"] > f(9, 1280),
        "strict_limit_loading_velocity": c["sharp_loaded_Y_u_floor_per_eta"] > 1/(24*c["L"]),
        "smaller_window_budget": c["budget"] == f(507, 125000) < f(1, 200),
        "same_transfer_error": c["transfer_error_upper"] == f(4, 399),
        "pair_error_retained_twice": c["balanced_pair_separation_factor"] == f(790, 399),
        "parent_kinetic_above_one": c["strict_kinetic_lower"] > 1,
        "parent_kinetic_below_five": c["strict_kinetic_upper"] < 5 < f(25, 4),
        "original_Q_gap": c["Q_pair_separation_per_eta"] == f(237, 42560),
        "original_Q_gap_above_gate": c["Q_gap_above_one_over_200"] == f(121, 212800) > 0,
    }
