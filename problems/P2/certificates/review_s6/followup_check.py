#!/usr/bin/env python
"""S6 review follow-up (2026-10-27): mechanical checks of the author's fix commit 43eed599.

  F1  scope: the commit changes only FORMULATION.md and the three S6 notes; nothing in src/, tests/,
      certificates/, CLAIMS.md or docs/ differs from the review commit fa1ad2fa;
  F2  S2: the s6-report §9 rows P2-6, P2-E and the three notes are byte-identical to the review's §8 texts;
  F3  FORMULATION: the file is the review-time file with exactly one v1.1 change ("and he informs the user"
      -> "and informs the user") plus an appended v1.2 entry that equals the review's §8 text, with the date
      2026-10-27 and the interim wording "(P2-6 and P2-E reviewed and approved 2026-10-27; not yet entered
      in the ledger)" (compared after collapsing whitespace, since lines are re-wrapped);
  F4  the pre-registration at fa1ad2fa is a byte prefix of the new one (the erratum is appended only);
  F5  W6: every line number cited in s6-proof §2.2 is one that quotes_check.py found for that passage;
  F6  the v1.4 draft (s6-report §11): numbers against the certificates and the review's own values, and a
      LaTeX sanity check (balanced braces, even count of $);
  F7  writes followup_tex_item_v1.4.txt: the draft with the review's one change in (c), and the item body
      wrapped in {\emergencystretch=3em ... \par} (typesetting only: without it the unbreakable 30-digit
      enclosure leaves a 92 pt overfull line);
  F8  if tectonic is on PATH: compiles, in a temporary directory outside the repository, the problems
      document as it is, with the draft inserted after \item[v1.3], and with the final text inserted, and
      counts overfull boxes (new ones must be none for the final text).

Run (from this directory):  python followup_check.py      -> followup_check_out.json
"""
from __future__ import annotations

import json
import re
import subprocess
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
P2 = HERE.parents[1]
ROOT = HERE.parents[3]
REV, FIX = "fa1ad2fa", "43eed599"


def git(*args) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, cwd=ROOT, check=True).stdout


def show(commit: str, path: str) -> str:
    return git("show", f"{commit}:{path}")


def collapse(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def row(text: str, rid: str) -> str:
    lines = [ln for ln in text.splitlines() if ln.startswith(f"| {rid} | P2 |")]
    if len(lines) != 1:
        raise ValueError(f"{rid}: {len(lines)} rows found")
    return lines[0]


def note(text: str, rid: str) -> str:
    """the backtick-quoted note text after the bullet '- **rid** (...)' (same line or the next one)"""
    m = re.search(r"- \*\*" + re.escape(rid) + r"\*\* \([^)]*\):\s*`( Note [^`]*)`", text)
    if not m:
        raise ValueError(f"note {rid} not found")
    return m.group(1)


def compile_check(draft: str, final: str):
    import shutil
    import tempfile
    if shutil.which("tectonic") is None:
        return {"skipped": "tectonic not on PATH"}
    src = (ROOT / "docs" / "problems" / "open-problems-theoretical-cosmology-2026.tex").read_text()
    lines = src.split("\n")
    i = [k for k, ln in enumerate(lines) if ln.startswith("\\item[v1.3]")][0]
    res = {}
    with tempfile.TemporaryDirectory() as td:
        for name, item in (("as_is", None), ("draft", draft), ("final", final)):
            body = lines if item is None else lines[:i + 1] + [item.replace("⟨Month day⟩", "October 27")] + lines[i + 1:]
            f = Path(td) / f"{name}.tex"
            f.write_text("\n".join(body))
            r = subprocess.run(["tectonic", "--keep-logs", f.name], cwd=td, capture_output=True, text=True)
            log = (Path(td) / f"{name}.log").read_text() if (Path(td) / f"{name}.log").exists() else ""
            over = re.findall(r"Overfull \\hbox \(([0-9.]+)pt too wide\) in paragraph at lines (\d+)--(\d+)", log)
            res[name] = {"returncode": r.returncode, "overfull": [(float(w), int(a), int(b)) for w, a, b in over]}
    base = {(round(w, 3)) for w, a, b in res["as_is"]["overfull"]}
    res["final_new_overfull"] = [o for o in res["final"]["overfull"] if round(o[0], 3) not in base]
    res["draft_new_overfull"] = [o for o in res["draft"]["overfull"] if round(o[0], 3) not in base]
    res["ok"] = all(res[k]["returncode"] == 0 for k in ("as_is", "draft", "final")) and not res["final_new_overfull"]
    return res


def main():
    out = {}
    # F1
    changed = git("diff", "--name-only", REV, FIX).split()
    out["F1_changed_files"] = changed
    out["F1_scope_ok"] = sorted(changed) == sorted([
        "problems/P2/FORMULATION.md", "problems/P2/notes/s6-preregistration.md",
        "problems/P2/notes/s6-proof.md", "problems/P2/notes/s6-report.md"])
    out["F1_tex_and_claims_unchanged_since_main"] = git(
        "diff", "--name-only", "300c22a7", FIX, "--", "CLAIMS.md", "docs/problems") == ""
    # F2
    review = show(REV, "docs/review-2026-10-27-p2-s6-adversarial.md")
    report = (P2 / "notes" / "s6-report.md").read_text()
    f2 = {}
    for rid in ("P2-6", "P2-E"):
        f2[f"row {rid}"] = row(review, rid) == row(report, rid)
    for rid in ("P2-K6", "P2-4", "P2-D1′"):
        f2[f"note {rid}"] = note(review, rid) == note(report, rid)
    out["F2_rows_and_notes_byte_identical"] = f2
    # F3
    old = show(REV, "problems/P2/FORMULATION.md")
    new = (P2 / "FORMULATION.md").read_text()
    a, b = "and he informs the user", "and informs the user"
    out["F3_v11_pronoun_once"] = old.count(a) == 1
    old2 = old.replace(a, b)
    out["F3_prefix_is_old_with_pronoun_fix"] = new.startswith(old2)
    appended = new[len(old2):]
    m = re.search(r"```\n(- \*\*v1\.2 .*?)```", review, flags=re.S)
    v12 = m.group(1).replace("⟨date⟩", "2026-10-27").replace(
        "P2-6 at P_B is assembled.",
        "P2-6 at P_B is assembled (P2-6 and P2-E reviewed and approved 2026-10-27; not yet entered in the ledger).")
    out["F3_appended_is_v12_only"] = collapse(appended) == collapse(v12)
    out["F3_appended_starts_with_blank_line"] = appended.startswith("\n- **v1.2 (2026-10-27,")
    out["F3_header_unchanged"] = new.split("\n## 1.")[0] == old2.split("\n## 1.")[0]
    # F4
    pre_old = show(REV, "problems/P2/notes/s6-preregistration.md").encode()
    pre_new = (P2 / "notes" / "s6-preregistration.md").read_bytes()
    out["F4_prereg_old_is_byte_prefix"] = pre_new.startswith(pre_old)
    out["F4_appended"] = pre_new[len(pre_old):].decode()
    # F5
    q = json.loads((HERE / "quotes_check_out.json").read_text())
    lines = {p["id"]: p["line"] for p in q["passages"]}
    proof = (P2 / "notes" / "s6-proof.md").read_text()
    sec = proof[proof.index("- **Smoothness and the threshold.**"):proof.index("- Near r₊, T = −x + (smooth)")]
    cited = {"V.nullspace_smooth": "ℓ. 1888", "H.states_smooth": "ℓ. 1264", "H.radial_above": "ℓ. 1065",
             "H.radial_below": "1067", "H.beta_j": "ℓ. 850", "H.kappa_j": "857", "H.arbitrary_regularity": "ℓ. 1549",
             "V.thm_half_plane": "ℓ. 1940", "V.commutator_sign": "ℓ. 1627"}
    f5 = {}
    for pid, s in cited.items():
        num = int(re.search(r"\d+", s).group(0))
        f5[pid] = {"cited": s, "in_text": s in sec, "quotes_check_line": lines[pid], "match": lines[pid] == num}
    out["F5_W6_line_numbers"] = f5
    # F6
    tex = re.search(r"```latex\n(.*?)\n```", report[report.index("## 11."):], flags=re.S).group(1)
    out["F6_draft"] = tex
    kz = json.loads((P2 / "certificates" / "part3" / "krawczyk_PB_l0.json").read_text())["enclosure_over_kappa_minus"]
    ind = json.loads((HERE / "indep_assembly_out.json").read_text())["point"]
    km = float(ind["kappa_-"].strip("[]").split("+/-")[0])
    kp_km = float(ind["kappa_+/kappa_-"].strip("[]").split("+/-")[0])
    g_km = float(ind["gamma0/kappa_-"].strip("[]").split("+/-")[0])
    f6 = {
        "beta ball matches P2-K6": "0.909380978994805455524933762622\\pm3.2\\times10^{-31}" in tex
        and kz["im"].startswith("[-0.909380978994805455524933762622 +/- 3.1"),
        "0.92 kappa_- = 0.547 gamma0 ~ 0.55": round(0.92 / g_km, 4) == 0.5466 and "0.92\\,\\kappa_-\\approx0.55\\,\\gamma_0" in tex,
        "-i kappa_+ ~ -0.823 i kappa_-": round(kp_km, 3) == 0.823 and "-0.823\\,i\\kappa_-" in tex,
        "l >= 20 and l <= 19": "\\ell\\ge20" in tex and "\\ell\\le19" in tex,
        "point": "\\Lambda M^2=3/50" in tex and "Q/Q_{\\max}=999/1000" in tex and "249/250" in tex,
        "braces balanced": tex.count("{") == tex.count("}"),
        "dollars even": tex.count("$") % 2 == 0,
        "only citation hintzvasy2017": sorted(set(re.findall(r"\\cite\{([^}]*)\}", tex))) == ["hintzvasy2017"],
        "bibitem exists in the .tex": "\\bibitem{hintzvasy2017}" in (ROOT / "docs" / "problems" /
                                                                   "open-problems-theoretical-cosmology-2026.tex").read_text(),
    }
    out["F6_checks"] = f6
    out["F6_kappa_minus"] = km
    # F7: the review's one change in (c): "certified" large-|Re w| exclusions for l <= 19 are not CERTIFIED-level
    old_c = (r"for $\ell\le19$, certified large-$|\operatorname{Re}\omega|$ exclusions and validated argument-principle "
             r"counts allow for the pole at $-i\kappa_+$;")
    new_c = (r"for $\ell\le19$, large-$|\operatorname{Re}\omega|$ exclusions now use the method of the large-$\ell$ "
             r"exclusion, and validated argument-principle counts allow for the pole at $-i\kappa_+$;")
    out["F7_old_clause_present_once"] = tex.count(old_c) == 1
    final = tex.replace(old_c, new_c)
    head = "\\item[v1.4] "
    assert final.startswith(head)
    final = head + "{\\emergencystretch=3em " + final[len(head):] + "\\par}"
    out["F7_final_braces_balanced"] = final.count("{") == final.count("}") and final.count("$") % 2 == 0
    (HERE / "followup_tex_item_v1.4.txt").write_text(final + "\n")
    out["F7_final"] = final
    out["F8_compile"] = compile_check(tex, final)
    flat = [out["F1_scope_ok"], out["F1_tex_and_claims_unchanged_since_main"], *f2.values(),
            out["F3_v11_pronoun_once"], out["F3_prefix_is_old_with_pronoun_fix"], out["F3_appended_is_v12_only"],
            out["F3_appended_starts_with_blank_line"], out["F3_header_unchanged"], out["F4_prereg_old_is_byte_prefix"],
            *(v["in_text"] and v["match"] for v in f5.values()), *f6.values(), out["F7_old_clause_present_once"],
            out["F7_final_braces_balanced"], out["F8_compile"].get("ok", True)]
    out["pass"] = all(flat)
    (HERE / "followup_check_out.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k not in ("F6_draft", "F7_final", "F4_appended")},
                     indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
