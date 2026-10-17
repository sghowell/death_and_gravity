"""Tests of the coverage-checked merge of scripts/s43_parts.py and its S4-4
hook: layer parts (scripts/s44_layer.py, p4.validated.tc_layer) are joined to
the S4-3 parts by piece id, and a valid passing one carries tau > tstar for a
pass-direct piece.  Synthetic plan and parts in a temporary directory; no
region code runs (the merge only imports the modules for their constants).
"""
import argparse
import copy
import json
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import s43_parts as s43  # noqa: E402

CODE_FILE = "problems/P4/FORMULATION.md"     # stands in for a source file: its current hash


def _write(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh)


class Env:
    def __init__(self, root):
        self.root = root
        self.parts = os.path.join(root, "parts")
        self.code = {CODE_FILE: s43.current_hash(CODE_FILE)}
        assert self.code[CODE_FILE]
        e1 = math.exp(-1.0)
        self.pieces = [
            dict(id="centre/c00", region="centre", kind="centre", index=0, lo=0.0, hi=e1),
            dict(id="tube/full-000", region="tube", kind="tube", cache="full",
                 cache_content_sha256="aaa", index=0, x=-0.02, h=0.98, x_min=-4.5,
                 x_max=-0.05, degrade=True, tq=True),
            dict(id="gap/gap-000", region="gap", kind="tube", cache="gap",
                 cache_content_sha256="bbb", index=0, x=-0.01, h=0.01, x_min=-0.056,
                 x_max=-0.008, degrade=True, tq=False),
            dict(id="sonic/s00", region="sonic", kind="sonic", index=0, lo=0.0, hi=0.01,
                 tq=False),
        ]
        self.byid = {p["id"]: p for p in self.pieces}
        _write(os.path.join(self.parts, "plan.json"),
               dict(schema="s43-plan/1", pieces=self.pieces, caches={}))
        for name, sha in (("full", "aaa"), ("gap", "bbb")):
            _write(os.path.join(self.parts, "caches", name + ".json"),
                   dict(content_sha256=sha))
        _write(os.path.join(self.parts, "checks.json"),
               dict(ok=True, code=dict(files=self.code),
                    results=dict(E1_U6top=True, E2_layers=[14, True],
                                 E3_sonic_zeros=[True], E5_sonic_vanishing=[True])))
        for pid, st in (("centre/c00", "pass-both"), ("tube/full-000", "pass-both"),
                        ("gap/gap-000", "pass-direct"), ("sonic/s00", "pass-direct")):
            self.s43_part(pid, st)

    def s43_part(self, pid, status):
        p = self.byid[pid]
        res = dict(ok=status.startswith("pass"), margin=-1.0,
                   open=("open (tau > tstar not certified on this piece)"
                         if status == "pass-direct" else None))
        _write(os.path.join(self.parts, pid + ".json"),
               dict(schema=s43.SCHEMA, id=pid, region=p["region"], spec=s43.spec_of(p),
                    status=status, result=res, extra=dict(centre_balls="certified"),
                    code=dict(files=self.code), timing=dict(cpu_s=1.0)))

    def layer_part(self, pid, dirname="layer", tau_min=224.0, **over):
        p = self.byid[pid]
        a, b = s43.x_interval(p)
        d = dict(schema="s44-layer/1", id=pid, region=p["region"], spec=s43.spec_of(p),
                 status="pass",
                 scope=dict(tau=[tau_min, "inf"], eps=[0.0, 6.0], xi=[-b, -a]),
                 certified_bound=-0.8,
                 result=dict(ok=True, certified_bound=-0.8, params=dict(tau_min=tau_min)),
                 inputs=(dict(kind="tube", cache=p["cache"],
                              cache_content_sha256=p["cache_content_sha256"])
                         if p["kind"] == "tube" else dict(kind="sonic")),
                 code=dict(files=dict(self.code)),
                 engine=dict(norm="layer", combine=False))
        for k, v in over.items():
            d[k] = v
        _write(os.path.join(self.parts, dirname, pid + ".json"), d)
        return d

    def merge(self):
        s43.cmd_merge(argparse.Namespace(accept_code=None, accept_reason=None))
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
    return Env(root)


def _gaps(d):
    return [tuple(g) for g in d["scope"]["tau_gt_tstar_uncovered_x"]]


def _detail(d, pid):
    return [e for e in d["detail"] if e["id"] == pid][0]


def test_merge_without_layer_parts_keeps_s43_semantics(env):
    name, d = env.merge()
    assert name == "tc_commutator.json" and d["ok"] and d["status"] == "complete"
    assert d["scope"]["tau_le_tstar_uncovered_x"] == []
    assert _gaps(d) == [(-0.02, 0.0)]                  # the band: direct phase only
    assert d["scope"]["tau_above_tstar"].startswith("OPEN")
    assert d["open_pieces"] == 2 and not d["strip_exclusion_certified"]
    assert d["scope"]["tau_gt_tstar_by"] == dict(tau_quadratic=2, layer=0, open=2)
    assert d["layer"]["joined_pieces"] == 0 and d["layer"]["dirs"] == []
    assert _detail(d, "gap/gap-000")["tau_gt_tstar"] == "open"
    assert _detail(d, "tube/full-000")["tau_gt_tstar"] == "tau-quadratic"


def test_layer_parts_cover_the_band(env):
    env.layer_part("gap/gap-000")
    env.layer_part("sonic/s00")
    name, d = env.merge()
    assert name == "tc_commutator.json" and d["ok"]
    assert _gaps(d) == [] and d["scope"]["tau_above_tstar"] == "certified everywhere"
    assert d["strip_exclusion_certified"] and d["open_pieces"] == 0
    assert d["scope"]["tau_gt_tstar_by"] == dict(tau_quadratic=2, layer=2, open=0)
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


def _bad_cases():
    return {
        "fail": (dict(status="fail", certified_bound=None), "fail"),
        "error": (dict(status="error", certified_bound=None), "error"),
        "spec": (dict(spec=dict(id="gap/gap-000", index=1)), "stale"),
        "code": (dict(code=dict(files={CODE_FILE: "0" * 64})), "stale"),
        "no-code": (dict(code=dict(files={})), "invalid"),
        "cache": (dict(inputs=dict(kind="tube", cache="gap", cache_content_sha256="ccc")),
                  "stale"),
        "schema": (dict(schema="s44-layer/0"), "invalid"),
        "positive": (dict(certified_bound=0.1,
                          result=dict(ok=True, certified_bound=0.1,
                                      params=dict(tau_min=224.0))), "invalid"),
        "not-ok": (dict(result=dict(ok=False, certified_bound=-0.8,
                                    params=dict(tau_min=224.0))), "invalid"),
        "eps": (dict(scope=dict(tau=[224.0, "inf"], eps=[0.0, 5.0], xi=[0.01, 0.02])),
                "invalid"),
        "xi": (dict(scope=dict(tau=[224.0, "inf"], eps=[0.0, 6.0], xi=[0.01, 0.021])),
               "invalid"),
        "tau-bounded": (dict(scope=dict(tau=[224.0, 1e6], eps=[0.0, 6.0], xi=[0.01, 0.02])),
                        "invalid"),
    }


@pytest.mark.parametrize("case", sorted(_bad_cases()))
def test_layer_part_rejected(env, case):
    over, want = _bad_cases()[case]
    env.layer_part("sonic/s00")
    env.layer_part("gap/gap-000", **over)
    p = env.byid["gap/gap-000"]
    d = json.load(open(os.path.join(env.parts, "layer", "gap/gap-000.json")))
    assert s43.layer_state(p, d, 224.0, 6.0) == want
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


def test_layer_parts_do_not_replace_s43_parts(env):
    # a layer part never stands in for a missing or failed S4-3 part, and it
    # does not change what carries tau > tstar on a pass-both piece
    os.remove(os.path.join(env.parts, "gap/gap-000.json"))
    env.layer_part("gap/gap-000", dirname="layer_t14", tau_min=14.0)
    env.layer_part("tube/full-000")
    env.layer_part("sonic/s00")
    name, d = env.merge()
    assert name == "tc_commutator_partial.json" and not d["ok"]
    assert d["missing"] == ["gap/gap-000"]
    assert d["scope"]["tau_le_tstar_uncovered_x"] == [[-0.02, -0.01]]
    assert _gaps(d) == [(-0.02, -0.01)]
    assert _detail(d, "tube/full-000")["tau_gt_tstar"] == "tau-quadratic"
    env.s43_part("gap/gap-000", "fail")
    _, d = env.merge()
    assert d["failed"] == ["gap/gap-000"] and _gaps(d) == [(-0.02, -0.01)]


def test_layer_code_hash_enters_conflicts(env):
    # a layer part computed by a different (accepted) source version shows up
    # in code_conflicts like any other used part
    env.layer_part("sonic/s00")
    old = "1" * 64
    env.layer_part("gap/gap-000", code=dict(files={CODE_FILE: old}))
    s43._HASH_CACHE.pop(CODE_FILE, None)
    s43.cmd_merge(argparse.Namespace(accept_code=["%s=%s" % (CODE_FILE, old)],
                                     accept_reason="test"))
    d = json.load(open(os.path.join(env.root, "tc_commutator.json")))
    assert d["scope"]["tau_gt_tstar_uncovered_x"] == []
    c = d["code_conflicts"][CODE_FILE]
    assert old in c["hashes"] and c["unaccepted"] == []
