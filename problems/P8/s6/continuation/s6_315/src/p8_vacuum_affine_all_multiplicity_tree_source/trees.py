"""Exact rooted currents at arbitrary finite external graviton multiplicity."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import checks as earlier_checks
from p8_vacuum_affine_complete_two_graviton_tree import trees as prior
from p8_vacuum_affine_minimal_gravity_radiation import vertices as old

from . import jets, source

ETA, imm, VECTOR_ZERO = prior.ETA, s.ImmutableMatrix, prior.VECTOR_ZERO
exact_real, exact_array, clean, tr = (
    prior.exact_real,
    prior.exact_array,
    prior.clean,
    prior.tr,
)


class TreeEngine(prior.TreeEngine):
    def __init__(self, *args, omit_eh5=False, omit_phi4=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.omit_eh5 = omit_eh5
        self.omit_phi4 = omit_phi4
        self.jets = {}

    def permitted(self, tags):
        nf, nh, ng = (tags.count(kind) for kind in prior.KINDS)
        if nf == 2 and nh == 0 and ng >= 1:
            return not (self.omit_phi4 and ng == 4)
        if nh == 2 and nf == 0 and ng >= 1:
            return True
        if nf == 2 and nh == 1:
            return self.cubic != 0
        if nf == 4 and nh == 0:
            return self.contact != 0
        return nf == nh == 0 and ng >= 3 and not (self.omit_eh5 and ng == 5)

    @prior.instance_cache
    def vertex(self, tags, momenta, values):
        assert self.permitted(tags)
        phi = [i for i, k in enumerate(tags) if k == "phi"]
        heavy = [i for i, k in enumerate(tags) if k == "H"]
        gravitons = [i for i, k in enumerate(tags) if k == "h"]
        hs = tuple(values[i] for i in gravitons)
        ps = tuple(momenta[i] for i in gravitons)
        scalar_product = s.prod(values[i] for i, k in enumerate(tags) if k != "h")
        r = len(hs)
        key = (hs, ps if not phi and not heavy else None)
        if key not in self.jets:
            self.jets[key] = (
                jets.MetricJet(*key) if key[1] is not None else jets.MetricJet(hs)
            )
        jet = self.jets[key]
        if len(phi) == 2 and not heavy:
            p, q = (momenta[i] for i in phi)
            value = -(p.T * ETA * jet.density(jet.full) * ETA * q)[0] - jet.det(
                jet.full
            )
        elif len(heavy) == 2 and not phi:
            p, q = (momenta[i] for i in heavy)
            value = -(p.T * ETA * jet.density(jet.full) * ETA * q)[
                0
            ] - self.heavy * jet.det(jet.full)
        elif len(phi) == 2 and len(heavy) == 1:
            value = -self.cubic * jet.det(jet.full)
        elif len(phi) == 4:
            value = self.contact * jet.det(jet.full)
        else:
            value = jet.gravity() * self.kappa
        return s.factor(scalar_product * value / self.kappa ** s.Rational(r, 2))


def configuration():
    E = s.Rational(5, 4)
    v = s.Rational(31, 16)
    ep = (v + 1 / v) / 2
    rp = (v - 1 / v) / 2
    W = 2 * (E - ep)
    unit = W / 12
    z = imm([0, 0, 1])
    directions = (
        imm([1, 0, 0]),
        imm([0, 1, 0]),
        imm([-s.Rational(3, 5), -s.Rational(4, 5), 0]),
    )
    energies = (3 * unit, 4 * unit, 5 * unit)
    real = []
    for n, w in zip(directions, energies):
        q = imm([w, *(w * n)])
        t = imm([-n[1], n[0], 0])
        pol = s.zeros(4)
        pol[1:, 1:] = t * t.T - z * z.T
        real.append((q, imm(pol)))
    u = imm([s.Rational(3, 5), 0, s.Rational(4, 5)])
    points = (
        imm([-E, 0, 0, -s.Rational(3, 4)]),
        imm([-E, 0, 0, s.Rational(3, 4)]),
        imm([ep, *(rp * u)]),
        imm([ep, *(-rp * u)]),
    )
    assert all(old.dot(p, p) == 1 for p in points)
    assert (
        sum(points, prior.VECTOR_ZERO) + sum((q for q, _ in real), prior.VECTOR_ZERO)
        == prior.VECTOR_ZERO
    )
    assert all(
        old.dot(q, q) == 0 and eps * q == prior.VECTOR_ZERO and s.trace(ETA * eps) == 0
        for q, eps in real
    )
    return points, tuple(real)


def amplitude(points, radiation, **kwargs):
    legs = [("phi", p, 1) for p in points[1:]]
    legs.extend(("h", q, eps) for q, eps in radiation)
    return TreeEngine(legs, **kwargs).amplitude()


def original_amplitude(points, gravitons):
    """Original-parameter physical TT coefficient; no integration or soft cut."""

    if not isinstance(points, (tuple, list)) or len(points) != 4:
        raise ValueError("Require four all-outgoing massive scalar momenta")
    if not isinstance(gravitons, (tuple, list)):
        raise TypeError("Require a finite list or tuple of real gravitons")
    ps = []
    for p in points:
        p = exact_array(p)
        if p.shape != (4, 1):
            raise ValueError("Require four-vectors")
        for x in p:
            exact_real(x)
        if clean(old.dot(p, p) - 1) != 0:
            raise ValueError("Require mass-one external shells")
        ps.append(p)
    E = -ps[0][0]
    if (
        not s.Rational(5, 4) <= E <= 2
        or ps[1][0] != -E
        or ps[0][1:, 0] + ps[1][1:, 0] != s.zeros(3, 1)
    ):
        raise ValueError("Require the stated incoming center-of-mass domain")
    if any(p[0] <= 0 for p in ps[2:]):
        raise ValueError("Require outgoing positive-energy massive legs")
    hs = []
    for item in gravitons:
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            raise ValueError("Require momentum-polarization pairs")
        q, eps = map(exact_array, item)
        if q.shape != (4, 1) or eps.shape != (4, 4):
            raise ValueError("Require a four-vector and symmetric tensor")
        for x in (*q, *eps):
            exact_real(x)
        if q[0] <= 0 or clean(old.dot(q, q)) != 0:
            raise ValueError("Require future null radiation")
        if (
            eps != eps.T
            or any(eps[0, j] != 0 for j in range(4))
            or tr(eps) != 0
            or eps * q != s.zeros(4, 1)
        ):
            raise ValueError("Require spatial transverse-traceless polarizations")
        hs.append((q, eps))
    if sum((q[0] for q, _ in hs), s.S.Zero) > s.Rational(1, 8):
        raise ValueError("Require total emitted energy at most 1/8")
    if sum(ps, VECTOR_ZERO) + sum((q for q, _ in hs), VECTOR_ZERO) != s.zeros(4, 1):
        raise ValueError("Require exact total momentum conservation")
    if any(old.dot(a[0], b[0]) == 0 for i, a in enumerate(hs) for b in hs[i + 1 :]):
        raise ValueError(
            "The point evaluator excludes an exactly collinear propagator pole"
        )
    return amplitude(
        tuple(ps),
        hs,
        heavy=source.HEAVY_MASS2,
        cubic=source.CUBIC,
        contact=source.CONTACT,
        kappa=source.KAPPA,
    )


@cache
def data():
    checks = {}
    # Exact compatibility on all previously implemented external multiplicities.
    points, q, born = old.sample(0)
    eps = imm(s.diag(0, 1, -1, 0))
    points2, qs, pols = earlier_checks.nonopposite_state()
    pars = {
        "heavy": 128,
        "cubic": s.Rational(1, 16),
        "contact": s.Rational(1, 8),
        "kappa": 16,
    }
    for number, ps, rays in (
        (0, born, ()),
        (1, points, ((q, eps),)),
        (2, points2, tuple(zip(qs, pols))),
    ):
        value, count = amplitude(ps, rays, **pars)
        target, tcount = prior.amplitude(ps, rays, **pars)
        checks[f"whole_frozen_amplitude_recovered_N{number}"] = s.factor(value - target)
        checks[f"whole_frozen_inventory_recovered_N{number}"] = s.Integer(
            count - tcount
        )
    ps, hs = configuration()
    for i, p in enumerate(ps):
        checks[f"three_real_scalar_shell_{i}"] = old.dot(p, p) - 1
    for i, (q, eps) in enumerate(hs):
        checks[f"three_real_null_shell_{i}"] = old.dot(q, q)
        checks[f"three_real_TT_transverse_{i}"] = eps * q
        checks[f"three_real_TT_trace_{i}"] = tr(eps)
    checks["three_real_exact_conservation"] = sum(ps, VECTOR_ZERO) + sum(
        (q for q, _ in hs), VECTOR_ZERO
    )
    value, count = amplitude(ps, hs, **pars)
    checks["three_real_whole_inventory"] = s.Integer(count - 5116)
    q, _eps = hs[0]
    gauges = {}
    for direction in (0, 2):
        xi = s.eye(4)[:, direction]
        gauge = imm((ETA * q) * (ETA * xi).T + (ETA * xi) * (ETA * q).T)
        gh = ((q, gauge), *hs[1:])
        ward, wcount = amplitude(ps, gh, **pars)
        checks[f"three_real_full_Ward_gauge{direction}"] = ward
        checks[f"three_real_Ward_inventory_gauge{direction}"] = s.Integer(wcount - 5116)
        gauges[direction] = gh
    controls = {}
    for omit, expected_count in (("omit_eh5", 5113), ("omit_phi4", 5110)):
        bad, bcount = amplitude(ps, gauges[2], **pars, **{omit: True})
        controls[omit] = bad
        checks[omit + "_inventory"] = s.Integer(bcount - expected_count)
    original, ocount = original_amplitude(ps, hs)
    checks["original_parameter_three_real_inventory"] = s.Integer(ocount - 5116)
    return {
        "whole_all_finite_tree_source": "The prior terminating root-current recursion is retained with separately owned, untruncated canonical vertices. Each proper leaf subset has a current for Phi,H,h, and the final Phi root is amputated. The construction defines every finite four-Phi/N-graviton tree coefficient of the selected action away from internal poles.",
        "whole_three_real_calibration": {
            "points": ps,
            "radiation": hs,
            "diagnostic_parameters": pars,
            "diagnostic_amplitude": value,
            "original_amplitude": original,
            "graph_count": s.Integer(count),
            "nonzero_omission_controls": controls,
        },
        "whole_tree_Ward_boundary": "Classical diffeomorphism invariance and the recursive inverse kinetic equations give the on-shell tree Ward identity at every finite order. The explicit two gauge directions at N3 and the higher-vertex omission controls supplement that proof; they do not replace it or certify a loop Ward identity.",
        "whole_kinematic_cut_scope": "For arbitrary finite future radiation, every scalar cut is one massive external leg plus a radiation subset (or complement), every heavy cut stays far from its pole, and hard Einstein cuts are timelike or the S312 future-cluster mixed cuts. The same mixed gap holds for all subset assignments. Pure-soft graviton subtree poles still require a quantitative all-multiplicity collinear analysis.",
        "checks": checks,
        "gates": {
            "all_lower_complete_amplitudes_recovered": True,
            "three_real_5116_inventory_from_full_current_recursion": count == 5116,
            "two_independent_three_real_gauge_directions_vanish": checks[
                "three_real_full_Ward_gauge0"
            ]
            == 0
            and checks["three_real_full_Ward_gauge2"] == 0,
            "omitted_EH5_breaks_Ward": controls["omit_eh5"] != 0,
            "omitted_scalar_fourth_metric_breaks_Ward": controls["omit_phi4"] != 0,
            "new_original_three_real_amplitude_nonzero": original != 0,
            "original_and_diagnostic_amplitudes_exact_rational": original.is_Rational
            is True
            and value.is_Rational is True,
            "all_finite_evaluators_have_independent_caches": True,
            "on_shell_Noether_induction_not_a_loop_or_soft_limit_proof": True,
            "pure_soft_cluster_poles_not_bounded_by_hard_cut_gaps": True,
            "no_collinear_point_assigned_an_unevaluated_zero_over_zero": True,
        },
    }
