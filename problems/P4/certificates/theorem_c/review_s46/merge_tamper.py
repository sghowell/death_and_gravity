#!/usr/bin/env python3
"""S4-6 review: end-to-end attacks on `s43_parts.py merge --substitute-centre`.

Runs the real cmd_merge in-process on the committed parts, with three
redirections and nothing else changed: the part reader serves edited copies
of chosen centre parts (or reports them missing), the centre identity
artifact and checks.json are read from edited temporary copies, and the
merge output goes to a temporary directory.  No tracked file is written.

Each scenario records status, strip_exclusion_certified, the centre
substitution (substituted, rejected with reasons, blocked, manual review)
and the expected outcome.  A scenario PASSES when the merge does what the
rule says it must (bad input rejected or blocked, good input accepted).

usage (repository root):
    PYTHONPATH=problems/P4/src python problems/P4/results/theorem_c/review_s46/merge_tamper.py [out.json]
"""
import contextlib
import copy
import io
import json
import os
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "problems", "P4", "scripts"))
import s43_parts as s43  # noqa: E402

MISSING = object()
ORIG_READ = s43.read_layer_part
OVR = {}


def read_layer_part(dirname, pid):
    if (dirname, pid) in OVR:
        v = OVR[(dirname, pid)]
        return None if v is MISSING else copy.deepcopy(v)
    return ORIG_READ(dirname, pid)


s43.read_layer_part = read_layer_part
TMP = tempfile.mkdtemp(prefix="s46_tamper_")
s43.CERT = os.path.join(TMP, "tc_commutator.json")
s43.CERT_PARTIAL = os.path.join(TMP, "tc_commutator_partial.json")
IDENT0 = json.load(open(s43.CENTRE_IDENTITY))
CHECKS0 = json.load(open(s43.CHECKS))
ORIG_IDENT, ORIG_CHECKS = s43.CENTRE_IDENTITY, s43.CHECKS


def part(pid, t=14):
    return json.load(open(os.path.join(s43.PARTS, "centre_t%d" % t, pid + ".json")))


def set_identity(rec):
    path = os.path.join(TMP, "centre_nb_identity.json")
    if os.path.exists(path):
        os.remove(path)
    if rec is not None:
        json.dump(rec, open(path, "w"))
    s43.CENTRE_IDENTITY = path


def set_checks(rec):
    path = os.path.join(TMP, "checks.json")
    json.dump(rec, open(path, "w"))
    s43.CHECKS = path


def run(layer=True, centre=True):
    for p in (s43.CERT, s43.CERT_PARTIAL):
        if os.path.exists(p):
            os.remove(p)
    args = types.SimpleNamespace(accept_code=None, accept_reason=None,
                                 substitute_layer_t14=layer, substitute_centre=centre)
    buf = io.StringIO()
    err = None
    with contextlib.redirect_stdout(buf):
        try:
            s43.cmd_merge(args)
        except SystemExit as e:
            err = str(e)
    out = s43.CERT if os.path.exists(s43.CERT) else s43.CERT_PARTIAL
    r = json.load(open(out))
    cs = r["centre_substitution"]
    return dict(status=r["status"], certified=r["strip_exclusion_certified"],
                substituted_total=len(r["detail"]) and sum(1 for e in r["detail"]
                                                            if e.get("substituted")),
                centre_substituted=cs["count"], centre_rejected=cs["rejected"],
                centre_blocked=cs["blocked"], centre_manual=cs["manual_review"],
                centre_no_part=cs["no_part"], layer_substituted=r["substitution"]["count"],
                missing=len(r["missing"]), failed=len(r["failed"]),
                gaps=r["scope"]["tau_gt_tstar_uncovered_x"],
                identity_valid=(cs["identity"] or {}).get("valid"), merge_error=err,
                stdout=buf.getvalue().strip()[-400:])


SCEN = []


def scenario(name, expect):
    def deco(fn):
        SCEN.append((name, expect, fn))
        return fn
    return deco


def edit(pid, fn, t=14):
    d = part(pid, t)
    fn(d)
    OVR[("centre_t%d" % t, pid)] = d


# ------------------------------------------------------------------ scenarios
@scenario("baseline: both rules, committed parts", "certified")
def s_base():
    pass


@scenario("default merge (no substitution flags)", "partial")
def s_default():
    return dict(layer=False, centre=False)


@scenario("--substitute-centre alone (layer pieces not substituted)", "partial")
def s_centre_only():
    return dict(layer=False, centre=True)


@scenario("--substitute-layer-t14 alone (centre missing)", "partial")
def s_layer_only():
    return dict(layer=True, centre=False)


@scenario("c29: certified_bound changed alone", "reject:centre/c29")
def s_bound():
    edit("centre/c29", lambda d: d.update(certified_bound=-0.5))


@scenario("c29: bound, result.certified_bound, result.worst all set to +0.01", "reject:centre/c29")
def s_bound_pos():
    def f(d):
        d["certified_bound"] = 0.01
        d["result"]["certified_bound"] = d["result"]["worst"] = 0.01
    edit("centre/c29", f)


@scenario("c29: bound fields forged consistently to -1.0 (no recompute in the merge)",
          "certified")
def s_forged():
    def f(d):
        d["certified_bound"] = -1.0
        d["result"]["certified_bound"] = d["result"]["worst"] = -1.0
    edit("centre/c29", f)


@scenario("c12: status fail", "reject:centre/c12")
def s_fail():
    edit("centre/c12", lambda d: d.update(status="fail"))


@scenario("c12: result.ok false", "reject:centre/c12")
def s_resok():
    edit("centre/c12", lambda d: d["result"].update(ok=False))


@scenario("c12: a failing box recorded", "reject:centre/c12")
def s_fails():
    edit("centre/c12", lambda d: d["result"].update(fails=[dict(box=[0] * 8, upper=1.0)]))


@scenario("c07: tc_layer_centre hash wrong (stale code)", "reject:centre/c07")
def s_hash():
    def f(d):
        k = "problems/P4/src/p4/validated/tc_layer_centre.py"
        d["code"]["files"][k] = "0" * 64
    edit("centre/c07", f)


@scenario("c07: tc_layer_centre not recorded", "reject:centre/c07")
def s_nomod():
    edit("centre/c07", lambda d: d["code"]["files"].pop(
        "problems/P4/src/p4/validated/tc_layer_centre.py"))


@scenario("c07: ec_centre not recorded", "reject:centre/c07")
def s_noec():
    edit("centre/c07", lambda d: d["code"]["files"].pop(
        "problems/P4/src/p4/validated/ec_centre.py"))


@scenario("c07: the mean-value module recorded (current hash)", "reject:centre/c07")
def s_mv():
    def f(d):
        k = "problems/P4/src/p4/validated/tc_layer_centre_mv.py"
        d["code"]["files"][k] = s43.sha256_file(os.path.join(ROOT, k))
    edit("centre/c07", f)


@scenario("c07: git_dirty_src true", "reject:centre/c07")
def s_dirty():
    edit("centre/c07", lambda d: d["code"].update(git_dirty_src=True))


@scenario("c07: git_head null", "reject:centre/c07")
def s_head():
    edit("centre/c07", lambda d: d["code"].update(git_head=None))


@scenario("c07: a late module", "reject:centre/c07")
def s_late():
    edit("centre/c07", lambda d: d["code"].update(
        late_modules=["problems/P4/src/p4/validated/tc_layer.py"]))


@scenario("c07: write check failed", "reject:centre/c07")
def s_wc():
    edit("centre/c07", lambda d: d["code"]["write_check"].update(ok=False))


@scenario("c07: driver_sha256 of another driver", "reject:centre/c07")
def s_driver():
    edit("centre/c07", lambda d: d.update(driver_sha256="1" * 64))


@scenario("c20: the tau >= 224 part placed as the tau >= 14 part", "reject:centre/c20")
def s_t224():
    OVR[("centre_t14", "centre/c20")] = part("centre/c20", 224)


@scenario("c20: scope tau changed to [14, inf] on the 224 part (params still 224)",
          "reject:centre/c20")
def s_t224b():
    d = part("centre/c20", 224)
    d["scope"]["tau"] = [14.0, "inf"]
    d["tau_min"] = 14.0
    OVR[("centre_t14", "centre/c20")] = d


@scenario("c20: tau range of the 224 part forged in all three places (s_lo still lo 224)",
          "reject:centre/c20")
def s_t224c():
    d = part("centre/c20", 224)
    d["scope"]["tau"] = [14.0, "inf"]
    d["tau_min"] = 14.0
    d["result"]["params"]["tau_min"] = 14.0
    OVR[("centre_t14", "centre/c20")] = d


@scenario("c20: the 224 part forged throughout (tau_min and s_lo): accepted, the merge "
          "cannot see a consistent forgery", "certified-forgery")
def s_t224d():
    d = part("centre/c20", 224)
    d["scope"]["tau"] = [14.0, "inf"]
    d["tau_min"] = 14.0
    d["result"]["params"]["tau_min"] = 14.0
    d["result"]["params"]["s_lo"] = d["spec"]["lo"] * 14.0 * (1 - 1e-12)
    OVR[("centre_t14", "centre/c20")] = d


@scenario("c20: tau upper end 1e6, not inf", "reject:centre/c20")
def s_tauhi():
    edit("centre/c20", lambda d: d["scope"].update(tau=[14.0, 1e6]))


@scenario("c03: one eps cell removed from the cover", "reject:centre/c03")
def s_eps():
    edit("centre/c03", lambda d: d["result"]["params"].update(
        eps_cells=d["result"]["params"]["eps_cells"][:-1]))


@scenario("c03: eps cells end at 5.0", "reject:centre/c03")
def s_eps2():
    def f(d):
        c = d["result"]["params"]["eps_cells"]
        c[-1] = [c[-1][0], 5.0]
    edit("centre/c03", f)


@scenario("c03: cover y-range [lo, hi/2]", "reject:centre/c03")
def s_y():
    edit("centre/c03", lambda d: d["result"]["params"].update(lb=d["result"]["params"]["lb"] / 2))


@scenario("c00: s_lo = 1e-6 on the closing piece", "reject:centre/c00")
def s_slo():
    edit("centre/c00", lambda d: d["result"]["params"].update(s_lo=1e-6))


@scenario("c05: s_lo above lo tau_min", "reject:centre/c05")
def s_slo2():
    edit("centre/c05", lambda d: d["result"]["params"].update(
        s_lo=d["result"]["params"]["s_lo"] * 1.01))


@scenario("c05: evaluator sigma_box", "reject:centre/c05")
def s_eval():
    def f(d):
        d["result"]["params"]["evaluator"] = "sigma_box"
        d["engine"]["evaluator"] = "sigma_box"
    edit("centre/c05", f)


@scenario("c05: laurent false", "reject:centre/c05")
def s_laur():
    edit("centre/c05", lambda d: d["result"].update(laurent=False))


@scenario("c05: combine on", "reject:centre/c05")
def s_comb():
    edit("centre/c05", lambda d: d["result"]["params"].update(combine=True))


@scenario("c05: data not raised", "reject:centre/c05")
def s_raised():
    edit("centre/c05", lambda d: d["engine"].update(raised=False))


@scenario("c05: identity block raise_by 3", "reject:centre/c05")
def s_idb():
    edit("centre/c05", lambda d: d["identity"].update(raise_by=3))


@scenario("c05: TM domain [lo, hi]", "reject:centre/c05")
def s_tmd():
    edit("centre/c05", lambda d: d["inputs"].update(tm_domain=[d["spec"]["lo"], d["spec"]["hi"]]))


@scenario("c05: centre series from the bare midpoint", "reject:centre/c05")
def s_ser():
    edit("centre/c05", lambda d: d["inputs"].update(
        centre_series="centre.centre_expansion(a3 midpoint)"))


@scenario("c06 served as c07 (id mismatch)", "reject:centre/c07")
def s_id():
    OVR[("centre_t14", "centre/c07")] = part("centre/c06")


@scenario("c08: spec lo changed", "reject:centre/c08")
def s_spec():
    edit("centre/c08", lambda d: d["spec"].update(lo=d["spec"]["lo"] * 1.0000001))


@scenario("c15: part missing", "partial")
def s_missing():
    OVR[("centre_t14", "centre/c15")] = MISSING


@scenario("c00: part missing (no closing piece)", "partial")
def s_missing0():
    OVR[("centre_t14", "centre/c00")] = MISSING


@scenario("identity artifact missing", "blocked")
def s_idmiss():
    set_identity(None)


@scenario("identity artifact ok = false", "blocked")
def s_idok():
    r = copy.deepcopy(IDENT0)
    r["ok"] = False
    set_identity(r)


@scenario("identity artifact: kappa^4 lowest nonzero layer 8", "blocked")
def s_idlay():
    r = copy.deepcopy(IDENT0)
    r["identity"]["U6"]["4"]["lowest_nonzero_layer"] = 8
    set_identity(r)


@scenario("identity artifact: odd terms at kappa^2", "blocked")
def s_idodd():
    r = copy.deepcopy(IDENT0)
    r["identity"]["U6"]["2"]["odd_terms"] = 1
    set_identity(r)


@scenario("identity artifact: kappa list 1..6", "blocked")
def s_idk():
    r = copy.deepcopy(IDENT0)
    r["identity"]["kappa"] = [1, 2, 3, 4, 5, 6]
    set_identity(r)


@scenario("identity artifact: stale tc_energy hash", "blocked")
def s_idhash():
    r = copy.deepcopy(IDENT0)
    r["code"]["files"]["problems/P4/src/p4/validated/tc_energy.py"] = "2" * 64
    set_identity(r)


@scenario("identity artifact: tc_layer_centre not recorded", "blocked")
def s_idmod():
    r = copy.deepcopy(IDENT0)
    r["code"]["files"].pop("problems/P4/src/p4/validated/tc_layer_centre.py")
    set_identity(r)


@scenario("identity artifact: wrong driver", "blocked")
def s_iddrv():
    r = copy.deepcopy(IDENT0)
    r["driver_sha256"] = "3" * 64
    set_identity(r)


@scenario("checks.json not ok", "blocked")
def s_checks():
    r = copy.deepcopy(CHECKS0)
    r["ok"] = False
    set_checks(r)


@scenario("checks.json with a stale code hash", "blocked")
def s_checks2():
    r = copy.deepcopy(CHECKS0)
    k = sorted(r["code"]["files"])[0]
    r["code"]["files"][k] = "4" * 64
    set_checks(r)


def judge(expect, r):
    if expect == "certified":
        return r["status"] == "complete" and r["certified"] and r["centre_substituted"] == 30 \
            and not r["centre_rejected"] and not r["centre_blocked"]
    if expect == "certified-forgery":
        return r["status"] == "complete" and r["certified"]
    if expect == "partial":
        return r["status"] == "partial" and not r["certified"]
    if expect == "blocked":
        return (r["centre_substituted"] == 0 and len(r["centre_blocked"]) > 0
                and not r["certified"] and r["status"] == "partial")
    if expect.startswith("reject:"):
        pid = expect.split(":", 1)[1]
        return (pid in r["centre_rejected"] and r["centre_substituted"] == 29
                and not r["certified"] and r["status"] == "partial")
    raise ValueError(expect)


def main(out_path):
    s43.import_p4()
    rows = []
    for name, expect, fn in SCEN:
        OVR.clear()
        s43.CENTRE_IDENTITY, s43.CHECKS = ORIG_IDENT, ORIG_CHECKS
        kw = fn() or {}
        r = run(**kw)
        ok = judge(expect, r)
        why = {k: v[:3] for k, v in r["centre_rejected"].items()}
        rows.append(dict(scenario=name, expect=expect, passed=ok, status=r["status"],
                         strip_exclusion_certified=r["certified"],
                         centre_substituted=r["centre_substituted"],
                         layer_substituted=r["layer_substituted"],
                         centre_rejected=why, centre_blocked=r["centre_blocked"][:4],
                         missing=r["missing"], failed=r["failed"],
                         gaps_tau_gt_tstar=r["gaps"], identity_valid=r["identity_valid"],
                         merge_error=r["merge_error"]))
        print("%-4s %-75s %s cert=%s centre=%d layer=%d rej=%s blocked=%d missing=%d gaps=%s"
              % ("ok" if ok else "FAIL", name[:75], r["status"], r["certified"],
                 r["centre_substituted"], r["layer_substituted"], list(why)[:2],
                 len(r["centre_blocked"]), r["missing"], r["gaps"][:2]), flush=True)
    rec = dict(schema="review-s46-merge-tamper/1", scenarios=rows,
               passed=sum(1 for r in rows if r["passed"]), total=len(rows),
               ok=all(r["passed"] for r in rows))
    with open(out_path, "w") as fh:
        json.dump(rec, fh, indent=1)
    print("%d/%d scenarios behave as the rule requires -> %s" % (rec["passed"], rec["total"],
                                                                 out_path))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "merge_tamper.json"))
