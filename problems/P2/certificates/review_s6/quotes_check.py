#!/usr/bin/env python
"""S6 review: the passages of [HV17] and [Vas13] that the review relies on for Lemma E (the notion of
resonance at the E_+ point -i kappa_+), checked verbatim against the arXiv TeX.

Sources (fetched by the reviewer from arXiv on 2026-10-27):
  HV17   arXiv:1512.08004v2 (2020-05-27; TeX date line "December 25, 2015. Final revision: July 12, 2017."),
         https://arxiv.org/e-print/1512.08004v2, file cauchyhor-r.tex,
         sha256 a0d3302f755222eb493170e1eea60f988a44c96f0e69108c89171e99733cf2be
  Vas13  arXiv:1012.4391v2 (2011-11-09; TeX date line "Revised: November 1, 2011. Original version:
         December 20, 2010"), https://arxiv.org/e-print/1012.4391v2, file kerr-ah-b.tex,
         sha256 b7ead254e1b5aa3654d8ba184a89f5fd4a1fb6a7e2cc79c00356f3030b83bc18

Each passage is searched in the LaTeX source with whitespace collapsed; the line where it starts is
recorded.  Run:  python quotes_check.py [--src DIR]   (DIR holding both .tex files; otherwise the
e-prints are downloaded)  -> quotes_check_out.json
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import tarfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = {
    "HV17": ("1512.08004v2", "cauchyhor-r.tex", "a0d3302f755222eb493170e1eea60f988a44c96f0e69108c89171e99733cf2be"),
    "Vas13": ("1012.4391v2", "kerr-ah-b.tex", "b7ead254e1b5aa3654d8ba184a89f5fd4a1fb6a7e2cc79c00356f3030b83bc18"),
}
PASSAGES = [
    ("HV17", "H.alpha_hypothesis", "the alpha hypothesis is about the Omega_23 forward problem",
     r"The only resonance of the forward problem in $\Omega_{23}$ in $\Im\sigma\geq 0$ is a simple resonance at $\sigma=0$"),
    ("HV17", "H.global_may_differ", "the global problem may have other resonances",
     r"(This does \emph{not} mean that the global problem for $\cP$ does not have other resonances in this half space!)"),
    ("HV17", "H.beta_j", "the radial-point quantity beta_j = 2/|mu'(r_j)|",
     r"\beta_j := -\frac{\tau_0^{-1}\rham_G\tau_0}{\wh\rho^{-1}\rham_G\wh\rho} = \frac{2}{|\mu_*'(r_j)|};"),
    ("HV17", "H.kappa_j", "kappa_j = 1/beta_j is the surface gravity",
     r"is equal to the \emph{surface gravity} of the horizon at $r=r_j$"),
    ("HV17", "H.radial_above", "above threshold: s > 1/2 + beta_j alpha (alpha = -Im sigma for a mode)",
     r"If $s\geq s'$, $s'>1/2+\beta_j\alpha$, and if $u\in\Hb^{-\infty,\alpha}(M)$ then $L_{j,\pm}$"),
    ("HV17", "H.radial_below", "below threshold",
     r"On the other hand, if $s<1/2+\beta_j\alpha$, and if $u\in\Hb^{-\infty,\alpha}(M)$"),
    ("HV17", "H.states_smooth", "resonant states are smooth at r_2, r_3 (above threshold)",
     r"Now the propagation of singularities at radial points implies that $u$ is smooth at $r=r_2$ and $r=r_3$ (where the a priori regularity exceeds the threshold value)"),
    ("HV17", "H.arbitrary_regularity", "arbitrarily high regularity may be used at the event and cosmological horizons",
     r"At the event horizon and the cosmological horizon, this is inconsequential, as we may work in spaces of arbitrary high regularity there"),
    ("HV17", "H.general_clause", "Lemma 2.14(1): the general clause (any sigma not an Omega_23 resonance)",
     r"In fact, this holds more generally for any $\sigma\in\C$ which is not a resonance of the forward problem for the wave equation in a neighborhood $\Omega_{23}$ of the black hole exterior."),
    ("Vas13", "V.fredholm", "the family is Fredholm above threshold",
     r"is Fredholm for each $\sigma$ with $s\geq m>(\difford-1-\beta\im\sigma)/2$"),
    ("Vas13", "V.thm_half_plane", "Thm (thm:classical-absorb): the half-plane C_s",
     r"\Cx_s=\{\sigma\in\Cx:\ \im\sigma>\beta^{-1}(\difford-1-2s)\}."),
    ("Vas13", "V.meromorphy", "analytic Fredholm theory: meromorphic inverse",
     r"Thus, analytic Fredholm theory applies, giving meromorphy of the inverse provided the inverse exists for a particular value of $\sigma$."),
    ("Vas13", "V.nullspace_smooth", "the nullspace above threshold is finite dimensional and smooth",
     r"is also finite dimensional, and its elements are in $\CI(X)$"),
    ("Vas13", "V.indep_of_s", "the inverse does not depend on s above threshold",
     r"Thus, in this sense, $(P_\sigma-\imath Q_\sigma)^{-1}$ is independent of $s$"),
    ("Vas13", "V.commutator_sign", "the radial estimate's sign condition s - (k-1)/2 + beta Im sigma > 0",
     r"is negative if $s-(\difford-1)/2+\beta\im\sigma-\delta>0$"),
    ("Vas13", "V.KdS_beta", "the Kerr-de Sitter beta (at a = 0: 2 r^2/|Delta'| = 1/kappa)",
     r"\beta_\pm=\mp 2\left(\frac{d\mut}{dr}\right)^{-1}(1+\gamma)(r^2+\alpha^2)|_{r=r_\pm}"),
]


def collapse(s: str) -> str:
    return re.sub(r"\s+", " ", s)


def linemap(txt: str) -> list:
    """for each character of collapse(txt), the 1-based line of the original character it comes from"""
    out = []
    line = 1
    prev_space = False
    for ch in txt:
        if ch.isspace():
            if not prev_space:
                out.append(line)
            prev_space = True
        else:
            out.append(line)
            prev_space = False
        if ch == "\n":
            line += 1
    return out


def fetch(arxiv_id: str, name: str) -> bytes:
    req = urllib.request.Request(f"https://arxiv.org/e-print/{arxiv_id}", headers={"User-Agent": "Mozilla/5.0"})
    data = urllib.request.urlopen(req, timeout=60).read()
    with tarfile.open(fileobj=io.BytesIO(data)) as tf:
        return tf.extractfile(name).read()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=None)
    args = ap.parse_args()
    texts, shas = {}, {}
    for key, (aid, name, sha) in SOURCES.items():
        raw = None
        if args.src and (Path(args.src) / name).exists():
            raw = (Path(args.src) / name).read_bytes()
        if raw is None:
            raw = fetch(aid, name)
        shas[key] = {"sha256": hashlib.sha256(raw).hexdigest(), "expected": sha,
                     "match": hashlib.sha256(raw).hexdigest() == sha}
        texts[key] = raw.decode("utf-8", errors="replace")
    out = {"sources": shas, "passages": []}
    for src, pid, what, passage in PASSAGES:
        txt = texts[src]
        flat = collapse(txt)
        pos = flat.find(collapse(passage))
        line = None
        if pos >= 0:
            line = linemap(txt)[pos]
        out["passages"].append({"source": src, "id": pid, "what": what, "found": pos >= 0, "line": line,
                                "passage": passage})
    out["all_found"] = all(p["found"] for p in out["passages"])
    out["sources_match"] = all(v["match"] for v in shas.values())
    (HERE / "quotes_check_out.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"all_found": out["all_found"], "sources_match": out["sources_match"],
                      "lines": {p["id"]: p["line"] for p in out["passages"]}}, indent=1))


if __name__ == "__main__":
    main()
