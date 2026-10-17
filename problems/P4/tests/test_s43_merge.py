"""Tests of the coverage-checked merge of scripts/s43_parts.py and its S4-4
hooks: layer parts (scripts/s44_layer.py, p4.validated.tc_layer) are joined to
the S4-3 parts by piece id (a valid passing one carries tau > tstar for a
pass-direct piece), and, opt-in (merge --substitute-layer-t14), a strictly
checked passing layer_t14 part carries a tube or sonic piece whose S4-3 part
is missing, stale or invalid (review
docs/review-2026-10-17-p4-layer-engine-adversarial.md, section 2).  Synthetic
plan and parts in a temporary directory; no region code runs (the merge only
imports the modules for their constants)."""
import argparse
import copy
import json
import math
import os
import sys
from fractions import Fraction

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import s43_parts as s43  # noqa: E402

CODE_FILE = "problems/P4/FORMULATION.md"     # stands in for a source file: its current hash
EPS_CELLS = [[0.0, 0.03], [0.03, 0.07], [0.07, 0.12], [0.12, 0.2], [0.2, 0.35], [0.35, 0.6],
             [0.6, 1.0], [1.0, 1.6], [1.6, 2.4], [2.4, 3.4], [3.4, 4.6], [4.6, 6.0]]


def _write(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh)


def _vp(m):
    return s43._vpath(m)


def _enc(x):
    """Tube.save's encoding [man, exp, rad_man, rad_exp] of the exact dyadic
    x (a float or a Fraction), radius 0."""
    q = Fraction(x)
    d = q.denominator
    assert d & (d - 1) == 0
    return [str(q.numerator), str(-(d.bit_length() - 1)), "0", "0"]


class Env:
    """Synthetic plan: centre (-inf, log 0.375], tube [-1, -0.02], gap
    [-0.02, -0.01], sonic [-0.01, 0]; the two tube-kind pieces are steps of
    small tube cache files (exact positions), as the merge's exact cover
    check reads them."""

    def __init__(self, root):
        self.root = root
        self.parts = os.path.join(root, "parts")
        self.cache_dir = os.path.join(root, "cache")
        self.code = {CODE_FILE: s43.current_hash(CODE_FILE)}
        assert self.code[CODE_FILE]
        self.sha = {"full": self.write_cache("full", [(-0.02, 0.98)]),
                    "gap": self.write_cache("gap", [(-0.01, 0.01)])}
        self.pieces = [
            dict(id="centre/c00", region="centre", kind="centre", index=0, lo=0.0, hi=0.375),
            dict(id="tube/full-000", region="tube", kind="tube", cache="full",
                 cache_content_sha256=self.sha["full"], index=0, x=-0.02, h=0.98,
                 x_min=-4.5, x_max=-0.02, degrade=True, tq=True),
            dict(id="gap/gap-000", region="gap", kind="tube", cache="gap",
                 cache_content_sha256=self.sha["gap"], index=0, x=-0.01, h=0.01,
                 x_min=-0.056, x_max=-0.008, degrade=True, tq=False),
            dict(id="sonic/s00", region="sonic", kind="sonic", index=0, lo=0.0, hi=0.01,
                 tq=False),
        ]
        self.write_plan()
        self.checks(True)
        for pid, st in (("centre/c00", "pass-both"), ("tube/full-000", "pass-both"),
                        ("gap/gap-000", "pass-direct"), ("sonic/s00", "pass-direct")):
            self.s43_part(pid, st)

    def write_cache(self, name, steps, rad=None):
        """A tube cache file of the steps [(x, h)] (exact positions; rad
        {k: radius mantissa} for nonzero radii) and its build record;
        returns its content hash."""
        rad = rad or {}
        enc = [_enc(x) for x, h in steps]
        for k, r in rad.items():
            enc[k][2] = str(r)
        path = os.path.join(self.cache_dir, "tube_%s.json" % name)
        _write(path, dict(K=1, info={}, steps=[
            dict(x=e, h=h, eps=["0", "0", "0", "0"], region="tube", co=[])
            for e, (x, h) in zip(enc, steps)]))
        sha = s43.tube_content_hash(path)
        _write(os.path.join(self.parts, "caches", name + ".json"), dict(content_sha256=sha))
        return sha

    def write_plan(self):
        self.byid = {p["id"]: p for p in self.pieces}
        _write(os.path.join(self.parts, "plan.json"),
               dict(schema="s43-plan/1", pieces=self.pieces, caches={}))

    def checks(self, ok):
        _write(os.path.join(self.parts, "checks.json"),
               dict(ok=ok, code=dict(files=self.code),
                    results=dict(E1_U6top=True, E2_layers=[14, True],
                                 E3_sonic_zeros=[True], E5_sonic_vanishing=[True])))

    def s43_part(self, pid, status, **over):
        p = self.byid[pid]
        res = dict(ok=status.startswith("pass"), margin=-1.0,
                   open=("open (tau > tstar not certified on this piece)"
                         if status == "pass-direct" else None))
        d = dict(schema=s43.SCHEMA, id=pid, region=p["region"], spec=s43.spec_of(p),
                 status=status, result=res, extra=dict(centre_balls="certified"),
                 code=dict(files=dict(self.code)), timing=dict(cpu_s=1.0))
        d.update(over)
        _write(os.path.join(self.parts, pid + ".json"), d)

    def drop_s43(self, pid):
        os.remove(os.path.join(self.parts, pid + ".json"))

    def layer_record(self, pid, tau_min=224.0, bound=-0.8):
        """A complete passing layer part as s44_layer.py writes it: tube pieces
        with the combine/piece/mv3 engine, sonic pieces with the default."""
        p = self.byid[pid]
        a, b = s43.x_interval(p)
        tube = p["kind"] == "tube"
        la, lb = (0.0, p["h"]) if tube else (p["lo"], p["hi"])
        xi_la = (-p["x"] + la) if tube else la
        mods = list(s43.LAYER_REQUIRED) + (["lintube", "tc_layer_mv"] if tube else [])
        files = {_vp(m): s43.current_hash(_vp(m)) for m in mods}
        assert all(files.values())
        eng = dict(norm="piece" if tube else "layer", combine=tube,
                   mv="mv3" if tube else None, piece_ladder_evals=None,
                   piece_final_evals=None)
        params = dict(tau_min=tau_min, eps_cells=copy.deepcopy(EPS_CELLS), la=la, lb=lb,
                      s_lo=xi_la * tau_min * (1 - 1e-12),
                      evaluator="sigma_box_mv3" if tube else "sigma_box", combine=tube,
                      tstar=224.0, target=bound / 2)
        return dict(
            schema="s44-layer/1", id=pid, region=p["region"], spec=s43.spec_of(p),
            status="pass", scope=dict(tau=[tau_min, "inf"], eps=[0.0, 6.0], xi=[-b, -a]),
            certified_bound=bound,
            result=dict(ok=True, worst=bound, certified_bound=bound, fails=[], boxes=10,
                        evals=20, params=params, laurent=not tube),
            inputs=(dict(kind="tube", cache=p["cache"],
                         cache_content_sha256=p["cache_content_sha256"])
                    if tube else dict(kind="sonic")),
            code=dict(files=files, fingerprint=s43.fingerprint(files), git_head="x",
                      git_dirty_src=False, hashed="import",
                      write_check=dict(ok=True, changed_sources={}, changed_drivers={})),
            engine=eng, driver_sha256=s43.sha256_file(s43.S44_DRIVER))

    def layer_part(self, pid, dirname="layer", tau_min=224.0, d=None, **over):
        d = d or self.layer_record(pid, tau_min)
        d.update(over)
        _write(os.path.join(self.parts, dirname, pid + ".json"), d)
        return d

    def merge(self, substitute=False, accept=None):
        s43.cmd_merge(argparse.Namespace(accept_code=accept, accept_reason=None,
                                         substitute_layer_t14=substitute))
        for name in ("tc_commutator.json", "tc_commutator_partial.json"):
            f = os.path.join(self.root, name)
            if os.path.exists(f):
                d = json.load(open(f))
                os.remove(f)
                return name, d
        raise AssertionError("merge wrote nothing")


@pytest.fixture
def env(tmp_path, monkeypatch):
    root = str(tmp_path)
    monkeypatch.setattr(s43, "PARTS", os.path.join(root, "parts"))
    monkeypatch.setattr(s43, "PLAN", os.path.join(root, "parts", "plan.json"))
    monkeypatch.setattr(s43, "CERT", os.path.join(root, "tc_commutator.json"))
    monkeypatch.setattr(s43, "CERT_PARTIAL", os.path.join(root, "tc_commutator_partial.json"))
    monkeypatch.setattr(s43, "CACHE_DIR", os.path.join(root, "cache"))
    return Env(root)


def _gaps(d):
    return [tuple(g) for g in d["scope"]["tau_gt_tstar_uncovered_x"]]


def _detail(d, pid):
    return [e for e in d["detail"] if e["id"] == pid][0]


# ------------------------------------------------- the tau > tstar join
def test_merge_without_layer_parts_keeps_s43_semantics(env):
    name, d = env.merge()
    assert name == "tc_commutator.json" and d["ok"] and d["status"] == "complete"
    assert d["scope"]["tau_le_tstar_uncovered_x"] == []
    assert _gaps(d) == [(-0.02, 0.0)]                  # the band: direct phase only
    t = d["scope"]["tau_above_tstar"]
    assert t.startswith("OPEN") and "nothing is certified" not in t
    assert "only tau in [b1, tstar] is certified on x in [(-0.02, 0.0)]" in t
    assert d["open_pieces"] == 2 and not d["strip_exclusion_certified"]
    assert d["scope"]["tau_gt_tstar_by"] == dict(tau_quadratic=2, layer=0, open=2,
                                                 substituted=0)
    assert d["cover"]["exact"] is True and d["cover"]["errors"] == []
    assert d["layer"]["joined_pieces"] == 0 and d["layer"]["dirs"] == []
    assert _detail(d, "gap/gap-000")["tau_gt_tstar"] == "open"
    assert _detail(d, "tube/full-000")["tau_gt_tstar"] == "tau-quadratic"
    assert d["substitution"]["enabled"] is False and d["substitution"]["count"] == 0


def test_layer_parts_cover_the_band(env):
    env.layer_part("gap/gap-000")
    env.layer_part("sonic/s00")
    name, d = env.merge()
    assert name == "tc_commutator.json" and d["ok"]
    assert _gaps(d) == [] and d["scope"]["tau_above_tstar"] == "certified everywhere"
    assert d["strip_exclusion_certified"] and d["open_pieces"] == 0
    assert d["scope"]["tau_gt_tstar_by"] == dict(tau_quadratic=2, layer=2, open=0,
                                                 substituted=0)
    for pid in ("gap/gap-000", "sonic/s00"):
        e = _detail(d, pid)
        assert e["tau_gt_tstar"] == "layer:layer"
        assert e["layer"]["layer"]["state"] == "pass"
        assert e["layer"]["layer"]["certified_bound"] == -0.8
    assert d["layer"]["states"] == {"layer": {"pass": 2}}
    assert d["layer"]["worst_bound"] == {"layer": -0.8}


def test_layer_parts_join_by_id_only(env):
    env.layer_part("sonic/s00")
    _, d = env.merge()
    assert _gaps(d) == [(-0.02, -0.01)]
    assert _detail(d, "gap/gap-000")["tau_gt_tstar"] == "open"
    assert _detail(d, "sonic/s00")["tau_gt_tstar"] == "layer:layer"


def _join_bad():
    """Rejections of the tau > tstar join (gap/gap-000, a tube piece)."""
    def res(**kw):
        def f(d):
            d["result"].update(kw)
        return f

    def prm(**kw):
        def f(d):
            d["result"]["params"].update(kw)
        return f

    def top(**kw):
        def f(d):
            d.update(kw)
        return f
    return {
        "fail": (top(status="fail", certified_bound=None), "fail"),
        "error": (top(status="error", certified_bound=None), "error"),
        "spec": (top(spec=dict(id="gap/gap-000", index=1)), "stale"),
        "code": (top(code=dict(files={CODE_FILE: "0" * 64})), "stale"),
        "no-code": (top(code=dict(files={})), "invalid"),
        "write-check-failed": (lambda d: d["code"].update(write_check=dict(ok=False)),
                               "invalid"),
        "cache-hash": (top(inputs=dict(kind="tube", cache="gap",
                                       cache_content_sha256="ccc")), "stale"),
        "cache-name": (lambda d: d.update(inputs=dict(d["inputs"], cache="full")),
                       "stale"),
        "schema": (top(schema="s44-layer/0"), "invalid"),
        "positive": (lambda d: (d.update(certified_bound=0.1),
                                d["result"].update(certified_bound=0.1, worst=0.1)),
                     "invalid"),
        "not-ok": (res(ok=False), "invalid"),
        "fails-recorded": (res(fails=[dict(box=[0], upper=1.0)]), "invalid"),
        "worst-differs": (res(worst=-0.7), "invalid"),
        "eps-scope": (lambda d: d["scope"].update(eps=[0.0, 5.0]), "invalid"),
        "eps-cells": (prm(eps_cells=EPS_CELLS[:5] + EPS_CELLS[6:]), "invalid"),
        "xi": (lambda d: d["scope"].update(xi=[0.01, 0.021]), "invalid"),
        "tau-bounded": (lambda d: d["scope"].update(tau=[224.0, 1e6]), "invalid"),
        "tau-min-scope": (prm(tau_min=200.0), "invalid"),
        "la": (prm(la=0.001), "invalid"),
        "lb": (prm(lb=0.005), "invalid"),
        "s-lo": (prm(s_lo=10.0), "invalid"),
        "laurent": (res(laurent=True), "invalid"),
        "evaluator": (prm(evaluator="sigma_box_other"), "invalid"),
        "evaluator-engine": (prm(evaluator="sigma_box"), "invalid"),
        "combine": (prm(combine=False), "invalid"),
    }


@pytest.mark.parametrize("case", sorted(_join_bad()))
def test_layer_part_rejected(env, case):
    fn, want = _join_bad()[case]
    env.layer_part("sonic/s00")
    d = env.layer_record("gap/gap-000")
    fn(d)
    env.layer_part("gap/gap-000", d=d)
    p = env.byid["gap/gap-000"]
    stored = json.load(open(os.path.join(env.parts, "layer", "gap/gap-000.json")))
    st, why = s43.layer_checks(p, stored, 224.0, 6.0)
    assert st == want and why
    _, m = env.merge()
    assert _gaps(m) == [(-0.02, -0.01)]
    e = _detail(m, "gap/gap-000")
    assert e["tau_gt_tstar"] == "open" and e["layer"]["layer"]["state"] == want
    assert m["layer"]["not_passing"]["layer"] == {want: ["gap/gap-000"]}


def test_layer_tau_min_above_tstar_does_not_cover(env):
    env.layer_part("gap/gap-000", tau_min=300.0)
    env.layer_part("sonic/s00")
    _, d = env.merge()
    assert _gaps(d) == [(-0.02, -0.01)]
    assert _detail(d, "gap/gap-000")["layer"]["layer"]["state"] == "invalid"
    assert any("tau_min" in r for r in d["layer_join_rejections"]["gap/gap-000 @ layer"])


def test_layer_t14_second_method_and_fallback(env):
    env.layer_part("gap/gap-000")
    env.layer_part("gap/gap-000", dirname="layer_t14", tau_min=14.0)
    env.layer_part("sonic/s00", dirname="layer_t14", tau_min=14.0)
    _, d = env.merge()
    assert _gaps(d) == []
    g, s = _detail(d, "gap/gap-000"), _detail(d, "sonic/s00")
    assert g["tau_gt_tstar"] == "layer:layer"          # parts/layer preferred
    assert s["tau_gt_tstar"] == "layer:layer_t14"      # tau >= 14 covers tau > 224
    assert g["direct_range_second_method"] == "layer_t14"
    assert s["direct_range_second_method"] == "layer_t14"
    assert d["layer"]["direct_range_second_method"] == 2
    assert d["layer"]["dirs"] == ["layer", "layer_t14"]


def test_layer_parts_do_not_replace_s43_parts_without_the_flag(env):
    # without --substitute-layer-t14 a layer part never stands in for a
    # missing or failed S4-3 part, and it does not change what carries
    # tau > tstar on a pass-both piece
    env.drop_s43("gap/gap-000")
    env.layer_part("gap/gap-000", dirname="layer_t14", tau_min=14.0)
    env.layer_part("tube/full-000")
    env.layer_part("sonic/s00")
    name, d = env.merge()
    assert name == "tc_commutator_partial.json" and not d["ok"]
    assert d["missing"] == ["gap/gap-000"]
    assert d["scope"]["tau_le_tstar_uncovered_x"] == [[-0.02, -0.01]]
    assert _gaps(d) == [(-0.02, -0.01)]
    t = d["scope"]["tau_above_tstar"]
    assert "nothing is certified for any tau >= b1 on x in [(-0.02, -0.01)]" in t
    assert "only tau in [b1, tstar]" not in t
    assert _detail(d, "tube/full-000")["tau_gt_tstar"] == "tau-quadratic"
    env.s43_part("gap/gap-000", "fail")
    _, d = env.merge()
    assert d["failed"] == ["gap/gap-000"] and _gaps(d) == [(-0.02, -0.01)]


def test_layer_code_hash_enters_conflicts(env):
    # a layer part computed by a different (accepted) source version shows up
    # in code_conflicts like any other used part
    env.layer_part("sonic/s00")
    old = "1" * 64
    d = env.layer_record("gap/gap-000")
    d["code"]["files"][_vp("tc_layer")] = old
    env.layer_part("gap/gap-000", d=d)
    s43._HASH_CACHE.pop(_vp("tc_layer"), None)
    s43.cmd_merge(argparse.Namespace(accept_code=["%s=%s" % (_vp("tc_layer"), old)],
                                     accept_reason="test", substitute_layer_t14=False))
    m = json.load(open(os.path.join(env.root, "tc_commutator.json")))
    assert m["scope"]["tau_gt_tstar_uncovered_x"] == []
    c = m["code_conflicts"][_vp("tc_layer")]
    assert old in c["hashes"] and c["unaccepted"] == []


def test_s43_part_with_failed_write_check_is_invalid(env):
    env.s43_part("tube/full-000", "pass-both",
                 code=dict(files=dict(env.code), write_check=dict(ok=False)))
    p = env.byid["tube/full-000"]
    assert s43.part_state(p, s43.read_part("tube/full-000")) == "invalid"
    _, d = env.merge()
    assert "tube/full-000" in d["missing"] and not d["ok"]


# --------------------------------------------- the substitution rule (opt-in)
def _sub_setup(env, pid="gap/gap-000"):
    """S4-3 part of pid removed; a passing layer_t14 part for it; the other
    band piece covered above tstar by its layer part."""
    env.drop_s43(pid)
    env.layer_part(pid, dirname="layer_t14", tau_min=14.0)
    other = "sonic/s00" if pid == "gap/gap-000" else "gap/gap-000"
    env.layer_part(other)


def test_substitution_is_off_by_default(env):
    _sub_setup(env)
    name, d = env.merge()
    assert name == "tc_commutator_partial.json"
    assert d["missing"] == ["gap/gap-000"] and d["substitution"]["count"] == 0


@pytest.mark.parametrize("pid", ["gap/gap-000", "sonic/s00"])
def test_substitution_of_a_missing_part(env, pid):
    _sub_setup(env, pid)
    name, d = env.merge(substitute=True)
    assert name == "tc_commutator.json" and d["ok"] and d["strip_exclusion_certified"]
    assert d["missing"] == [] and d["substitution"]["substituted"] == [pid]
    assert d["scope"]["tau_le_tstar_uncovered_x"] == [] and _gaps(d) == []
    e = _detail(d, pid)
    assert e["substituted"] == "layer_t14" and e["s43_state"] == "missing"
    assert e["tau_le_tstar"] == e["tau_gt_tstar"] == "layer:layer_t14"
    assert "dT != 0" in e["nonvanishing"] and "Rbar" in e["nonvanishing"]
    assert d["pieces"] == len(d["detail"]) == 4


@pytest.mark.parametrize("how", ["spec", "code", "write-check"])
def test_substitution_of_a_stale_or_invalid_part(env, how):
    env.layer_part("gap/gap-000", dirname="layer_t14", tau_min=14.0)
    env.layer_part("sonic/s00")
    if how == "spec":
        env.s43_part("gap/gap-000", "pass-direct", spec=dict(id="gap/gap-000", index=7))
        want = "stale"
    elif how == "code":
        env.s43_part("gap/gap-000", "pass-direct", code=dict(files={CODE_FILE: "2" * 64}))
        want = "stale"
    else:
        env.s43_part("gap/gap-000", "pass-direct",
                     code=dict(files=dict(env.code), write_check=dict(ok=False)))
        want = "invalid"
    name, d = env.merge(substitute=True)
    assert name == "tc_commutator.json" and d["strip_exclusion_certified"]
    assert d["substitution"]["substituted"] == ["gap/gap-000"]
    assert _detail(d, "gap/gap-000")["s43_state"] == want


@pytest.mark.parametrize("st43", ["fail", "error"])
def test_failed_s43_part_goes_to_manual_review(env, st43):
    env.s43_part("gap/gap-000", st43)
    env.layer_part("gap/gap-000", dirname="layer_t14", tau_min=14.0)
    env.layer_part("sonic/s00")
    name, d = env.merge(substitute=True)
    assert name == "tc_commutator_partial.json"
    assert d["substitution"]["manual_review"] == ["gap/gap-000"]
    assert d["substitution"]["substituted"] == [] and d["failed"] == ["gap/gap-000"]


def test_centre_is_never_substituted(env):
    env.drop_s43("centre/c00")
    _write(os.path.join(env.parts, "layer_t14", "centre/c00.json"),
           dict(schema="s44-layer/1", id="centre/c00", status="pass"))
    _, d = env.merge(substitute=True)
    assert "centre/c00" in d["missing"] and d["substitution"]["substituted"] == []
    assert "centre/c00" not in d["substitution"]["rejected"]


def _sub_bad():
    """Rejections specific to the substitution rule (strict checks), on the
    tube piece gap/gap-000 unless the key says sonic."""
    def drop(m):
        def f(d, env):
            del d["code"]["files"][_vp(m)]
        return f

    def code(**kw):
        def f(d, env):
            d["code"].update(kw)
        return f

    def prm(**kw):
        def f(d, env):
            d["result"]["params"].update(kw)
        return f

    def top(**kw):
        def f(d, env):
            d.update(kw)
        return f
    out = {
        "schema": top(schema="s44-layer/0"),
        "id": top(id="gap/gap-001"),
        "spec": top(spec=dict(id="gap/gap-000", index=3)),
        "cache-hash": top(inputs=dict(kind="tube", cache="gap", cache_content_sha256="ccc")),
        "cache-name": lambda d, env: d.update(inputs=dict(d["inputs"], cache="deep")),
        "no-code": code(files={}),
        "hash-not-current": lambda d, env: d["code"]["files"].update(
            {_vp("tc_layer"): "3" * 64}),
        "dirty-tree": code(git_dirty_src=True),
        "dirty-unknown": code(git_dirty_src=None),
        # written outside a checkout: before 2026-10-17 git_head() returned
        # (None, False) there (review section 6, new issue 1)
        "git-head-missing": code(git_head=None),
        # a module first loaded after the import-time snapshot, as code_record
        # writes it (review section 6, new issue 2)
        "late-module": code(late_modules=[_vp("tc_reduced_eigs")],
                            write_check=dict(ok=False, changed_sources={},
                                             changed_drivers={},
                                             late_modules=[_vp("tc_reduced_eigs")])),
        "no-write-check": lambda d, env: d["code"].pop("write_check"),
        "write-check-failed": code(write_check=dict(ok=False)),
        "driver-pin": top(driver_sha256="4" * 64),
        "mv-module-without-engine": lambda d, env: d["engine"].update(mv=None),
        "engine-mv-without-module": drop("tc_layer_mv"),
        "status-fail": top(status="fail"),
        "status-error": top(status="error"),
        "not-ok": lambda d, env: d["result"].update(ok=False),
        "fails-recorded": lambda d, env: d["result"].update(fails=[dict(box=[0])]),
        "bound-result": lambda d, env: d["result"].update(certified_bound=-0.7),
        "bound-worst": lambda d, env: d["result"].update(worst=-0.7),
        "bound-nonfinite": lambda d, env: (d.update(certified_bound=float("nan"))),
        "bound-nonnegative": lambda d, env: (d.update(certified_bound=0.0),
                                             d["result"].update(certified_bound=0.0,
                                                                worst=0.0)),
        "tau-min-above-b1": None,
        "tau-min-scope": prm(tau_min=13.0),
        "tau-bounded": lambda d, env: d["scope"].update(tau=[14.0, 1e8]),
        "eps-cells-gap": prm(eps_cells=EPS_CELLS[:3] + EPS_CELLS[4:]),
        "eps-cells-short": prm(eps_cells=EPS_CELLS[:-1]),
        "eps-cells-backward": prm(eps_cells=[[0.0, 3.0], [3.0, 3.0], [3.0, 6.0]]),
        "la": prm(la=1e-4),
        "lb": prm(lb=0.0099),
        "xi-scope": lambda d, env: d["scope"].update(xi=[0.0101, 0.02]),
        "s-lo": prm(s_lo=0.2),
        "laurent": lambda d, env: d["result"].update(laurent=True),
        "evaluator-unknown": prm(evaluator="sigma_box_fast"),
        "evaluator-engine": prm(evaluator="sigma_box_mv"),
        "combine": prm(combine=False),
        "sonic/la": prm(la=1e-9),
        "sonic/lb": prm(lb=0.02),
        "sonic/closing-s-lo": prm(s_lo=1e-9),
        "sonic/mv-evaluator": lambda d, env: (d["result"]["params"].update(
            evaluator="sigma_box_mv3"), d["engine"].update(mv="mv3"),
            d["code"]["files"].update({_vp("tc_layer_mv"): s43.current_hash(
                _vp("tc_layer_mv"))})),
        "sonic/combine": lambda d, env: (d["result"]["params"].update(combine=True),
                                         d["engine"].update(combine=True)),
        "sonic/laurent": lambda d, env: d["result"].update(laurent=False),
    }
    for m in s43.LAYER_REQUIRED + ("lintube",):
        out["module-" + m] = drop(m)
    return out


@pytest.mark.parametrize("case", sorted(_sub_bad()))
def test_substitution_rejected(env, case):
    pid = "sonic/s00" if case.startswith("sonic/") else "gap/gap-000"
    other = "gap/gap-000" if pid == "sonic/s00" else "sonic/s00"
    env.drop_s43(pid)
    env.layer_part(other)
    if case == "tau-min-above-b1":
        d = env.layer_record(pid, tau_min=224.0)
    else:
        d = env.layer_record(pid, tau_min=14.0)
        _sub_bad()[case](d, env)
    env.layer_part(pid, dirname="layer_t14", d=d)
    name, m = env.merge(substitute=True)
    assert name == "tc_commutator_partial.json"
    assert m["substitution"]["substituted"] == [] and pid in m["missing"]
    assert m["substitution"]["rejected"][pid], case


def test_substitution_ignores_accept_code(env):
    # --accept-code makes an old hash acceptable for the S4-3 parts and the
    # join, never for a substituted piece
    env.drop_s43("gap/gap-000")
    env.layer_part("sonic/s00")
    old = "5" * 64
    d = env.layer_record("gap/gap-000", tau_min=14.0)
    d["code"]["files"][_vp("tc_layer")] = old
    env.layer_part("gap/gap-000", dirname="layer_t14", d=d)
    _, m = env.merge(substitute=True, accept=["%s=%s" % (_vp("tc_layer"), old)])
    assert m["substitution"]["substituted"] == []
    assert m["substitution"]["rejected"]["gap/gap-000"] == ["source hashes not current"]


def test_substitution_blocked_by_invalid_checks(env):
    _sub_setup(env)
    env.checks(False)
    _, d = env.merge(substitute=True)
    assert d["substitution"]["substituted"] == []
    assert d["substitution"]["blocked"] == ["checks.json not valid for the current code"]
    assert "gap/gap-000" in d["missing"]


def test_substitution_needs_checks_valid_for_the_current_code(env):
    # checks.json computed by an older (accepted) code version keeps the
    # merge's checks_valid_for_code, but blocks substitution
    _sub_setup(env)
    old = "7" * 64
    _write(os.path.join(env.parts, "checks.json"),
           dict(ok=True, code=dict(files={CODE_FILE: old}),
                results=dict(E1_U6top=True, E2_layers=[14, True],
                             E3_sonic_zeros=[True], E5_sonic_vanishing=[True])))
    s43._HASH_CACHE.pop(CODE_FILE, None)
    _, d = env.merge(substitute=True, accept=["%s=%s" % (CODE_FILE, old)])
    assert d["checks_valid_for_code"] is True
    assert d["substitution"]["blocked"] == ["checks.json not valid for the current code"]
    assert d["substitution"]["substituted"] == []


def test_substitution_blocked_by_changed_constants(env, monkeypatch):
    import p4.validated.tc_layer as tl
    _sub_setup(env)
    monkeypatch.setattr(tl, "EPS_MAX", 5.0)
    _, d = env.merge(substitute=True)
    assert d["substitution"]["blocked"] == ["b1, c or EPS_MAX changed"]
    assert d["substitution"]["substituted"] == []


def test_substituted_code_joins_the_conflict_check(env):
    # two passing parts with different (current vs accepted-old) hashes of one
    # file: the substituted part's hash is current, the join part's old one is
    # accepted; both enter code_union, so the conflict is reported
    _sub_setup(env)
    old = "6" * 64
    d = env.layer_record("sonic/s00")
    d["code"]["files"][_vp("tc_tm")] = old
    env.layer_part("sonic/s00", d=d)
    s43._HASH_CACHE.pop(_vp("tc_tm"), None)
    _, m = env.merge(substitute=True, accept=["%s=%s" % (_vp("tc_tm"), old)])
    assert m["substitution"]["substituted"] == ["gap/gap-000"]
    c = m["code_conflicts"][_vp("tc_tm")]
    assert old in c["hashes"] and s43.current_hash(_vp("tc_tm")) in c["hashes"]


# ------------------------------------- provenance (review section 6, new issues)
def _git_fails(monkeypatch):
    """Every git call exits 128 with empty output, as outside a checkout."""
    import subprocess

    def run(argv, **kw):
        return subprocess.CompletedProcess(argv, 128, stdout="",
                                           stderr="fatal: not a git repository")
    monkeypatch.setattr(s43.subprocess, "run", run)


def test_git_head_in_the_checkout():
    head, dirty = s43.git_head()
    if head is None:
        pytest.skip("not run from a git checkout (e.g. a git archive copy)")
    assert len(head) == 40 and isinstance(dirty, bool)


def test_git_head_when_git_fails(monkeypatch):
    # new issue 1: a failed `git status` used to read as a clean tree
    _git_fails(monkeypatch)
    assert s43.git_head() == (None, None)


def test_git_head_outside_a_checkout(tmp_path, monkeypatch):
    monkeypatch.setattr(s43, "ROOT", str(tmp_path))
    assert s43.git_head() == (None, None)


def test_git_head_of_a_copy_inside_another_checkout(tmp_path, monkeypatch):
    # a `git archive` copy placed (and ignored) inside some other repository:
    # git answers for that repository, not for the copy
    import subprocess
    outer = tmp_path / "outer"
    copy = outer / "copy"
    (copy / "problems" / "P4" / "src").mkdir(parents=True)
    (outer / ".gitignore").write_text("copy/\n")
    g = ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false",
         "-c", "core.hooksPath=/dev/null", "-C", str(outer)]
    subprocess.run(g + ["init", "-q"], check=True)
    subprocess.run(g + ["add", ".gitignore"], check=True)
    subprocess.run(g + ["commit", "-q", "-m", "x"], check=True)
    monkeypatch.setattr(s43, "ROOT", str(copy))
    assert s43.git_head() == (None, None)


def test_part_written_without_git_is_not_substituted(env, monkeypatch):
    # end to end: the code record of a process without git records git_head
    # and git_dirty_src as None, and the strict checks reject the part
    with monkeypatch.context() as mp:
        _git_fails(mp)
        rec = s43.code_record()
    assert rec["git_head"] is None and rec["git_dirty_src"] is None
    env.drop_s43("gap/gap-000")
    env.layer_part("sonic/s00")
    d = env.layer_record("gap/gap-000", tau_min=14.0)
    d["code"].update(git_head=rec["git_head"], git_dirty_src=rec["git_dirty_src"])
    env.layer_part("gap/gap-000", dirname="layer_t14", d=d)
    _, m = env.merge(substitute=True)
    assert m["substitution"]["substituted"] == [] and "gap/gap-000" in m["missing"]
    why = m["substitution"]["rejected"]["gap/gap-000"]
    assert "git_head not recorded (not written from a checkout)" in why
    assert "git_dirty_src not false" in why


def _late_probe(tmp_path, monkeypatch):
    """A p4 module first loaded after the import-time snapshot: a private
    copy of the snapshot, then a module object entered into sys.modules."""
    monkeypatch.setattr(s43, "_IMPORT_HASH", dict(s43._IMPORT_HASH))
    s43.snapshot_code()
    src = tmp_path / "late_probe.py"
    src.write_text("X = 1\n")
    mod = type(sys)("p4._late_probe")
    mod.__file__ = str(src)
    monkeypatch.setitem(sys.modules, "p4._late_probe", mod)
    return s43.rel(str(src))


def test_no_late_module_after_the_snapshot(monkeypatch):
    monkeypatch.setattr(s43, "_IMPORT_HASH", dict(s43._IMPORT_HASH))
    s43.import_p4()
    rec = s43.code_record()
    assert rec["late_modules"] == [] and rec["write_check"]["ok"] is True
    assert s43.write_check_ok(dict(code=rec))


def test_late_module_is_flagged(tmp_path, monkeypatch):
    r = _late_probe(tmp_path, monkeypatch)
    assert s43.late_modules() == [r]
    assert r in s43.loaded_code() and s43.late_modules() == [r]   # no silent snapshot
    rec = s43.code_record()
    assert rec["late_modules"] == [r] and r in rec["files"]
    assert rec["write_check"]["ok"] is False and rec["write_check"]["late_modules"] == [r]
    assert not s43.write_check_ok(dict(code=rec))
    # a record that lists a late module fails even with a passing write check
    assert not s43.write_check_ok(dict(code=dict(rec, write_check=dict(ok=True))))
    s43.snapshot_code()                       # an explicit import-time snapshot
    assert s43.late_modules() == []


def test_late_module_part_is_rejected(tmp_path, env, monkeypatch):
    # a layer part whose code record lists a late module: invalid for the join
    # and rejected by the substitution rule; an S4-3 part with one is invalid
    r = _late_probe(tmp_path, monkeypatch)
    rec = s43.code_record()
    env.drop_s43("gap/gap-000")
    d = env.layer_record("sonic/s00")
    d["code"].update(late_modules=rec["late_modules"], write_check=rec["write_check"])
    env.layer_part("sonic/s00", d=d)
    d = env.layer_record("gap/gap-000", tau_min=14.0)
    d["code"].update(late_modules=rec["late_modules"], write_check=rec["write_check"])
    env.layer_part("gap/gap-000", dirname="layer_t14", d=d)
    _, m = env.merge(substitute=True)
    want = ["modules first loaded after the import-time snapshot: %s" % r]
    assert m["substitution"]["substituted"] == []
    assert m["substitution"]["rejected"]["gap/gap-000"] == want
    assert m["layer_join_rejections"]["sonic/s00 @ layer"] == want
    env.s43_part("tube/full-000", "pass-both",
                 code=dict(files=dict(env.code), late_modules=[r],
                           write_check=dict(ok=True)))
    assert s43.part_state(env.byid["tube/full-000"], s43.read_part("tube/full-000")) == "invalid"


# ----------------------------------- the exact x-cover (S4-5 assembly note B1)
X0 = Fraction(-0.01) + Fraction(1, 2 ** 60)      # an exact step position, not a double


def _two_gap_steps(env, x1=None, h1=0.0051):
    """The gap cache as two steps from X0 (h0 = 0.005; x1 = X0 - h0 unless
    given), the plan's gap pieces gap-000 and gap-001 labelled as the drivers
    label them (float(x), h), both pass-direct."""
    x1 = X0 - Fraction(0.005) if x1 is None else x1
    sha = env.write_cache("gap", [(X0, 0.005), (x1, h1)])
    g0 = dict(env.byid["gap/gap-000"], cache_content_sha256=sha, x=float(X0), h=0.005)
    g1 = dict(g0, id="gap/gap-001", index=1, x=float(x1), h=h1)
    env.pieces = [p for p in env.pieces if p["region"] != "gap"] + [g0, g1]
    env.write_plan()
    for pid in ("gap/gap-000", "gap/gap-001"):
        env.s43_part(pid, "pass-direct")
    return g0, g1


def test_cover_is_exact_where_the_float_labels_do_not_tile(env):
    g0, g1 = _two_gap_steps(env)
    # in exact arithmetic the labels leave a hole between the two steps ...
    assert Fraction(g1["x"]) < Fraction(g0["x"]) - Fraction(g0["h"])
    # ... the exact cache positions join exactly, and the merge sees no gap
    name, d = env.merge()
    assert name == "tc_commutator.json" and d["ok"]
    assert d["cover"]["errors"] == [] and d["cover"]["pieces"] == 5
    assert d["scope"]["tau_le_tstar_uncovered_x"] == []


def test_merge_fails_on_an_inexact_join(env):
    # x1 = X0 - h0 - 2^-70: the second step starts 2^-70 below the first's end
    _two_gap_steps(env, x1=X0 - Fraction(0.005) - Fraction(1, 2 ** 70))
    with pytest.raises(SystemExit, match="exact cover check failed"):
        s43.cmd_merge(argparse.Namespace(accept_code=None, accept_reason=None,
                                         substitute_layer_t14=False))
    d = json.load(open(os.path.join(env.root, "tc_commutator_partial.json")))
    assert not d["ok"] and not d["strip_exclusion_certified"]
    assert any("steps 0 and 1 do not join exactly" in e for e in d["cover"]["errors"])


def _cover_bad():
    """exact_x_cover errors: (plan pieces as changes of the env's, the
    expected message)."""
    def gap(**kw):
        def f(env):
            env.pieces = [dict(p, **kw) if p["id"] == "gap/gap-000" else p for p in env.pieces]
        return f

    def cache(steps, **kw):
        def f(env):
            sha = env.write_cache("gap", steps, **kw)
            gap(cache_content_sha256=sha)(env)
        return f

    def sonic_pair(env):
        s1 = dict(env.byid["sonic/s00"], id="sonic/s01", index=1, lo=0.0100001, hi=0.02)
        env.pieces = env.pieces + [s1]
    return {
        "claims-another-step": (gap(x=-0.0101), "not the cache step it claims"),
        "h-label": (gap(h=0.0099), "not the cache step it claims"),
        "content-hash": (gap(cache_content_sha256="ccc"), "differs from the plan's"),
        "index": (gap(index=3), "outside the 1 steps"),
        "radius": (cache([(-0.01, 0.01)], rad={0: 1}), "nonzero radius"),
        "join": (cache([(-0.01, 0.005), (-0.0150001, 0.005)]), "do not join exactly"),
        "missing-file": (gap(cache="deep"), "missing"),
        "sonic-touch": (sonic_pair, "do not touch exactly"),
    }


@pytest.mark.parametrize("case", sorted(_cover_bad()))
def test_exact_cover_errors(env, case):
    fn, want = _cover_bad()[case]
    fn(env)
    env.write_plan()
    iv, errors = s43.exact_x_cover(s43.load_plan(), [p["id"] for p in env.pieces])
    assert any(want in e for e in errors), (case, errors)


def test_exact_cover_intervals(env):
    iv, errors = s43.exact_x_cover(s43.load_plan(), [p["id"] for p in env.pieces])
    assert errors == []
    q = lambda p: (p.q, p.t)
    assert q(iv["gap/gap-000"][0]) == (Fraction(-0.01) - Fraction(0.01), None)
    assert q(iv["gap/gap-000"][1]) == (Fraction(-0.01), None)
    assert q(iv["sonic/s00"][0]) == (-Fraction(0.01), None)
    assert q(iv["centre/c00"][0]) == (None, 0) and q(iv["centre/c00"][1]) == (None, Fraction(0.375))
    assert s43.exact_gaps(list(iv.values())) == []


def test_exact_gaps_with_centre_pieces():
    X, Q = s43.XPoint, Fraction
    tube = (X(q=Q(-1)), X(q=Q(-1, 50)))
    # log 0.375 = -0.98 > -1: the centre piece overlaps the tube step
    assert s43.exact_gaps([(X(t=Q(0)), X(t=Q(3, 8))), tube], hi=X(q=Q(-1, 50))) == []
    # log 0.25 = -1.386 < -1: a hole (log 0.25, -1), found without tolerance
    g = s43.exact_gaps([(X(t=Q(0)), X(t=Q(1, 4))), tube], hi=X(q=Q(-1, 50)))
    assert [(float(a), float(b)) for a, b in g] == [(math.log(0.25), -1.0)]
    # touching: log 1 = 0 meets [0, 1] exactly
    assert s43.exact_gaps([(X(t=Q(0)), X(t=Q(1))), (X(q=Q(0)), X(q=Q(1)))],
                          hi=X(q=Q(1))) == []
    # a gap of 2^-80 between two dyadic ends is a gap
    g = s43.exact_gaps([(X(q=Q(-2)), X(q=Q(-1))), (X(q=Q(-1) + Q(1, 2 ** 80)), X(q=Q(0)))])
    assert len(g) == 2 and g[1][0].q == Q(-1) and g[1][1].q == Q(-1) + Q(1, 2 ** 80)
