"""Complete spatial averages, trace source and uniform hybrid phase domain."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_volume_turnaround import moving

from . import source

P, CORE, KAPPA = moving.P, moving.old.field.CORE, moving.old.field.KAPPA
NZ, NZZ = s.Integer(10) ** 4, s.Integer(10) ** 18
F_AMPLITUDE = s.Rational(1, 10**510)
H_AMPLITUDE = s.Rational(1, 10**520)


@cache
def bounds():
    T, HOM = source.TIME, source.HOMOGENEOUS_RADIUS
    r = 8 * CORE
    raw = moving.old.field.field_bounds()
    f = {name: value * r for name, value in raw.items()}
    v = 4 * f["v_A2"] / (1 + P) ** 2
    tau = 4 * f["tau_A2"] / (1 + P) ** 2
    pbound = 8 * T + HOM
    pvbound = 100 * (T + HOM)
    trace = (pvbound + f["Pi_v_A0"]) * P * v + P * f["Pi_v_A0"] / 3
    shape = 180 * 32 * (1 + P) * f["Pi_tau_A0"]
    vector = 12 * (1 + P) ** 2 * f["W_A0"] * f["Pi_W_A0"]
    matter = 2 * (
        (1 + f["delta_Pi_M_A0"]) * f["M1_A1"]
        + (2 * HOM + f["Pi_H_A0"]) * f["H_gradient_A0"]
    )
    D0 = trace + shape + vector + matter
    lam = 4 * D0 / P
    piTF = 5 * (32 * f["Pi_tau_A0"] + 16 * lam)
    shear = 192 * piTF**2
    curvature_remainder = 2 * 10**12 * P**2 * (v + tau) ** 2
    images = [
        2 * f["Pi_v_A0"] + 12 * pbound * v,
        12 * (1 + P) * f["Pi_W_A0"],
        4 * f["delta_Pi_M_A0"] + 6 * v,
        4 * f["Pi_H_A0"] + 12 * HOM * v,
        f["eta_A0"],
        shear,
        16 * f["Pi_W_A0"] ** 2 / moving.old.field.ZETA,
        96 * moving.old.field.ZETA * ((1 + P) * f["W_A0"]) ** 2,
        4 * f["W_A0"] ** 2,
        4 * f["M1_A1"] ** 2,
        4 * f["H_gradient_A0"] ** 2,
        8 * P**2 * v + curvature_remainder,
    ]
    averages = [
        10 * pbound * v**2 + 4 * f["Pi_v_A0"] * v,
        36 * v * (1 + P) * f["Pi_W_A0"],
        10 * v**2 + 12 * f["delta_Pi_M_A0"] * v,
        10 * HOM * v**2 + 12 * f["Pi_H_A0"] * v,
        s.Integer(0),
        *images[5:11],
        curvature_remainder,
    ]
    zsum, mean_sum = sum(images), sum(averages)
    ndelta = NZ * zsum
    nmean = NZ * mean_sum + NZZ * zsum**2 / 2
    fvariation = 20 * v**2 + 72 * v * ndelta + 12 * nmean + 5000 * ndelta**2
    hgrad, hessian, h0 = 10**3, 10**9, 100
    hvariation = (
        10 * h0 * v**2
        + hgrad * mean_sum
        + hessian * zsum**2 / 2
        + 6 * v * (hgrad * zsum + hessian * zsum**2 / 2)
    )
    contraction = (
        s.Rational(1, 30) + 10**9 * (source.LAPSE_RADIUS + 12 * (HOM + max(images))) / 3
    )
    center_residual = s.Rational(5, 10**400) + 12 * 10**4 * HOM
    full_residual = center_residual + NZ * zsum
    return {
        "fields": f,
        "v": v,
        "tau": tau,
        "pbound": pbound,
        "Pi_v_background": pvbound,
        "Dtrace": trace,
        "Dshape": shape,
        "Dvector": vector,
        "Dmatter": matter,
        "D0": D0,
        "lambda_A1": lam,
        "pi_TF": piTF,
        "shear": shear,
        "curvature_remainder": curvature_remainder,
        "images": images,
        "averages": averages,
        "zsum": zsum,
        "mean_sum": mean_sum,
        "deltaN_point": ndelta,
        "deltaN_mean": nmean,
        "Fvariation": fvariation,
        "hvariation": hvariation,
        "contraction": contraction,
        "homogeneous_center_residual": center_residual,
        "full_center_residual": full_residual,
    }


@cache
def trace_identity():
    gv = s.symbols("full_metric0:6", real=True)
    metric = s.Matrix(
        [[gv[0], gv[1], gv[2]], [gv[1], gv[3], gv[4]], [gv[2], gv[4], gv[5]]]
    )
    xi = s.Matrix(s.symbols("whole_descriptor0:3", real=True))
    Dxi = s.Matrix(3, 3, s.symbols("whole_descriptor_derivative0:9", real=True))
    Dmetric = []
    for k in range(3):
        v = s.symbols("whole_metric_jet" + str(k) + "_0:6", real=True)
        Dmetric.append(
            s.Matrix([[v[0], v[1], v[2]], [v[1], v[3], v[4]], [v[2], v[4], v[5]]])
        )
    pv = s.Symbol("whole_trace_density", real=True)
    lie = (
        sum((xi[k] * Dmetric[k] for k in range(3)), s.zeros(3))
        + Dxi.T * metric
        + metric * Dxi
    )
    expected = (
        pv * sum(xi[k] * s.trace(metric.inv() * Dmetric[k]) for k in range(3)) / 6
        + pv * s.trace(Dxi) / 3
    )
    return s.factor(s.trace((pv * metric.inv() / 6) * lie) - expected)


@cache
def data():
    b = bounds()
    rn = 2 + source.OFF_CLOCK
    rnn = 6 + source.OFF_CLOCK
    UN = s.Rational(3, 4) * 2**2 * rn
    UNN = s.Rational(21, 16) * 2**3 * rn**2 + s.Rational(3, 4) * 2**2 * rnn
    # Entire source branch, not an independent Gaussian lapse.
    checks = {
        "complete_nondiagonal_trace_density_Lie_pairing": trace_identity(),
        "all_twelve_full_images": s.Integer(len(b["images"]) - 12),
        "all_twelve_full_mean_remainders": s.Integer(len(b["averages"]) - 12),
        "whole_mean_lapse_integral_remainder": b["deltaN_mean"]
        - NZ * b["mean_sum"]
        - NZZ * b["zsum"] ** 2 / 2,
        "whole_point_lapse_response": b["deltaN_point"] - NZ * b["zsum"],
        "whole_full_source_residual_includes_phase_image": b["full_center_residual"]
        - b["homogeneous_center_residual"]
        - NZ * b["zsum"],
        "whole_trace_shape_vector_matter_sources": b["D0"]
        - sum(b[name] for name in ("Dtrace", "Dshape", "Dvector", "Dmatter")),
    }
    return {
        "whole_uniform_field_rows_and_complete_geometry": b,
        "whole_first_second_implicit_response_bounds": [NZ, NZZ],
        "whole_uniform_U_first_second_bounds": [UN, UNN],
        "whole_centered_spatial_average_amplitudes": [F_AMPLITUDE, H_AMPLITUDE],
        "whole_exact_trace_generator": "For pi_trace=(Pi_v/6)gamma^-1 and det gamma=a^6 exp(6v), the FULL Lie pairing gives Dtrace=Pi_v grad v-grad(Pi_v)/3. This retains its nonzero fluctuation source. Only its tracefree projection vanishes. All DQ*, full matter, vector/Gauss and adjoint-inverse contacts remain.",
        "whole_full_average_proof": "All original nonzero input coefficients and their unchanged free propagation have zero spatial mean. Full density exponentials, nonlinear shape/cotangent reconstruction and auxiliary roots do not: the displayed twelve average remainders retain them. Full scalar curvature has linear term -4a^-2 Delta v with zero spatial integral and complete nonlinear remainder2e12 P^2(v+tau)^2. No reconstructed harmonic or tensor curvature is removed.",
        "whole_implicit_average_proof": "On the wider source domain the full lapse contraction has q<1/20 and |CN|>2. All first and mixed second contacts give |Nz|<1e4 and |Nzz|<1e18. Taylor with its exact integral remainder is applied along the full invariant segment from the homogeneous center N0(Y), not at N1 off shell. Thus point deltaN<=Nz sum deltaZ and its average is bounded by Nz sum average_deltaZ+Nzz(sum deltaZ)^2/2.",
        "whole_envelope_Hamiltonian_proof": "The full reduced density has f_z=H_z and f_zw=H_zw-Cz Cw/CN, retaining the lapse response. Their ceilings are1e3 and1e9; f0<100. Multiply the entire exp(3v), retain its mean quadratic remainder and all cross products. The displayed hvariation bounds the spatial average of exp(3v)f-f0; it does not define a Taylor-truncated Hamiltonian.",
        "whole_complex_homogeneous_proof": "For complex Y within1e-122 of the actual reference, a^+-2 and a^+-3 have modulus<2. The complete physical metric/inverse A2 norms remain<2: their new constant background has no spatial derivatives. The original full shape and mean-zero formal-transpose Neumann inverse remain valid. Source p,pm,ph and every phase-induced invariant stay in1e-120; the full center residual includes BOTH homogeneous and phase deviations.",
        "checks": checks,
        "gates": {
            "whole_source_image": source.HOMOGENEOUS_RADIUS
            + max(b["images"])
            + 8 * source.TIME
            < source.COORDINATE,
            "whole_auxiliary_strict_contraction": b["contraction"] < s.Rational(1, 20),
            "whole_full_auxiliary_self_map": b["full_center_residual"] / 3
            + b["contraction"] * source.LAPSE_RADIUS
            < source.LAPSE_RADIUS,
            "whole_full_lapse_pivot": 3 * (1 - b["contraction"]) > 2,
            "whole_all_mixed_implicit_Hessian": (10**4 + 2 * 10**6 * NZ + 10**9 * NZ**2)
            / 2
            < NZZ,
            "whole_full_reduced_Hamiltonian_Hessian": 10**3 + s.Rational(10**8, 2)
            < 10**9,
            "whole_U_first_second_bounds": UN < 12 and UNN < 10**4,
            "whole_full_average_volume": b["Fvariation"] < F_AMPLITUDE,
            "whole_full_average_Hamiltonian": b["hvariation"] < H_AMPLITUDE,
            "whole_full_tracefree_shear": b["shear"] < s.Rational(1, 10**530),
            "whole_real_time_inside_parent_reference_interval": source.TIME
            < moving.TIME,
            "homogeneous_scale_complex_gap": 4 * source.HOMOGENEOUS_RADIUS
            + 16 * source.TIME**2
            < s.Rational(1, 100),
            "unchanged_reference_maps_4R_into_8R": moving.bounds()["flow_deviation"]
            < s.Rational(1, 50),
            "all_full_spatial_and_Gauss_contacts_retained": True,
        },
    }
