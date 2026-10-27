#!/usr/bin/env python
"""S5 review: check the passages of [HV17] and [CF17] that P2-5 rests on, against the arXiv TeX.

The review read the sources itself; this script makes the key readings replayable.  It looks for each
passage verbatim (LaTeX source, whitespace collapsed) and records the line where it starts.

Sources (fetched by the review on 2026-10-26; the script fetches them again if --src lacks them):
  HV17  arXiv:1512.08004v2 (posted 2020-05-27; TeX date line "December 25, 2015. Final revision:
        July 12, 2017."), e-print https://arxiv.org/e-print/1512.08004, file cauchyhor-r.tex,
        sha256 a0d3302f755222eb493170e1eea60f988a44c96f0e69108c89171e99733cf2be
  CF17  arXiv:1607.01018v2 (posted 2017-10-23), e-print https://arxiv.org/e-print/1607.01018,
        file RN_deSitterRevised6.tex,
        sha256 dcaab0e03e43622592f6f2d336bf131a36085d627d62efdf252903a37464d453

Run:  python quotes_check.py [--src DIR]     (writes quotes_check_out.json next to this file)
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "quotes_check_out.json"
SOURCES = {
    "HV17": ("1512.08004", "cauchyhor-r.tex", "a0d3302f755222eb493170e1eea60f988a44c96f0e69108c89171e99733cf2be"),
    "CF17": ("1607.01018", "RN_deSitterRevised6.tex", "dcaab0e03e43622592f6f2d336bf131a36085d627d62efdf252903a37464d453"),
}

# (source, id, what the review uses it for, verbatim LaTeX passage)
PASSAGES = [
    ("HV17", "H.thm11_setting", "Thm 1.1: non-degenerate RNdS with non-zero charge",
     r"Let $g$ be a non-degenerate Reissner--Nordstr\"om--de Sitter metric with non-zero charge $Q$"),
    ("HV17", "H.thm11_alpha_exists", "Thm 1.1 is existential in alpha",
     r"Then there exists $\alpha>0$, only depending on the parameters of the spacetime"),
    ("HV17", "H.thm11_regularity", "Thm 1.1: the H^{1/2+alpha/kappa_1-0} conclusion",
     r"lie in the weighted spacetime Sobolev space $e^{-\alpha t_0}H^{1/2+\alpha/\kappa_1-0}$ in $t_0>0$"),
    ("HV17", "H.t0_v_event", "t_0 ~ v near the event horizon (Price-law conversion)",
     r"and to the Eddington--Finkelstein coordinate $v$ near the event horizon"),
    ("HV17", "H.alpha_below_gap", "alpha is any number BELOW the spectral gap (intro)",
     r"Now, $\alpha$ can be any real number less than the \emph{spectral gap} $\alpha_0$"),
    ("HV17", "H.two_ended", "two-ended H_I only by 'a minor extension'",
     r"a minor extension of our arguments yields analogous boundedness, decay and regularity results for the Cauchy problem with a `two-ended' Cauchy surface $H_I$ up to the bifurcation sphere $B$"),
    ("HV17", "H.def21", "Def 2.1: non-degenerate",
     r"is \emph{non-degenerate} if $\mu$ has $3$ simple positive roots $0<r_1<r_2<r_3$"),
    ("HV17", "H.HI", "H_I is one-ended: r_2 - 2 delta <= r <= r_3 + 2 delta",
     r"H_I &= \{ t_*=0,\ r_2-2\delta\leq r\leq r_3+2\delta \}"),
    ("HV17", "H.dod", "the domain of dependence of H_I is Omega n {r > r_1}",
     r"Forward solutions to the wave equation $\Box_g u=f$ in the domain of dependence of $H_I$, i.e.\ in $\Omega\cap\{r>r_1\}$, are not affected by any modifications of the operator $\Box_g$ outside, i.e.\ in $r\leq r_1$."),
    ("HV17", "H.U_support_typo", "the stated support of the absorption region (a typo: it must lie in r < r_1)",
     r"\cU\subset\{r_1-\delta<r<r_2-\delta,\ \tau\leq e^{-(t_{*,0}+2)}\}"),
    ("HV17", "H.prop212", "Prop 2.12: high-energy estimates down to -gamma > -gamma_0",
     r"provided $h^{-1}\Im z\in[-\gamma,c_2]$, where $-\gamma_0<-\gamma<c_2$"),
    ("HV17", "H.lemma214_general", "Lemma 2.14(1), the 'more general' clause",
     r"In fact, this holds more generally for any $\sigma\in\C$ which is not a resonance of the forward problem for the wave equation in a neighborhood $\Omega_{23}$ of the black hole exterior."),
    ("HV17", "H.lemma214_2", "Lemma 2.14(2): R o P(sigma)^-1 at 0, rank 1, constants",
     r"then the only pole of $R\circ\wh\cP(\sigma)^{-1}$ with $\Im\sigma\geq 0$ is at $\sigma=0$, has rank $1$, and the space of resonant states consists of constant functions."),
    ("HV17", "H.lemma214_outline", "the general clause is proved in outline",
     r"The more general statement follows along the same lines (and is in fact much easier to prove, since it does not entail a mode stability statement)"),
    ("HV17", "H.lemma214_order_cited", "the order-1 claim at 0 is taken from the SdS/KdS references",
     r"which is known to be equal to $1$, see the references above"),
    ("HV17", "H.sec25_hyp", "§2.5: the alpha hypothesis is about the forward problem in Omega_23",
     r"The only resonance of the forward problem in $\Omega_{23}$ in $\Im\sigma\geq 0$ is a simple resonance at $\sigma=0$, with resonant states equal to constants, see the references given in the proof of Lemma~\ref{LemmaRNdSResonanceSupp}, and there exists $\alpha>0$ such that $0$ is the only resonance in $\Im\sigma\geq-\alpha$. (This does \emph{not} mean that the global problem for $\cP$ does not have other resonances in this half space!)"),
    ("HV17", "H.sec25_gamma0", "§2.5: alpha < gamma_0 for the high-energy estimates",
     r"we may assume $-\alpha>-\gamma_0$ so that we have high energy estimates in $\Im\sigma\geq-\alpha$."),
    ("HV17", "H.prop217_contour", "proof of Prop 2.17: the expansion runs over global resonances in the OPEN strip",
     r"where the $\sigma_j$ are the resonances of $\cP$ in $-\alpha<\Im\sigma<-\wt\alpha$"),
    ("HV17", "H.prop217_removal", "proof of Prop 2.17: global resonances removed via Lemma 2.14(1)",
     r"all resonant states of $\cP$ which are not resonant states of the forward problem in $\Omega_{23}$ must in fact vanish in $r>r_1$"),
    ("HV17", "H.cor219", "Cor 2.19: |u'| <= C tau^alpha, continuity up to CH+",
     r"there exists a constant $C>0$ such that $|u'(\tau,x)|\leq C\tau^\alpha$. In particular, $u$ is uniformly bounded in $r>r_1$ and extends continuously to $\cCH^+$."),
    ("HV17", "H.thm224", "Thm 2.24: conormal regularity at the Cauchy horizon",
     r"V_1\cdots V_N u'\in\Hb^{1/2+\alpha\beta_1-0,\alpha}"),
    ("HV17", "H.beta1", "beta_1 = 2/|mu'(r_1)|, kappa_j = 1/beta_j",
     r"\beta_j := -\frac{\tau_0^{-1}\rham_G\tau_0}{\wh\rho^{-1}\rham_G\wh\rho} = \frac{2}{|\mu_*'(r_j)|};"),
    ("HV17", "H.gamma0", "§2.7: gamma_0 formula",
     r"\alpha < \gamma_0 := \frac{\mu(r_P)\numin}{4} = \frac{1}{2r_P}\bigl(\mu(r_P)(2-3r_P^{-1}\bhm)\bigr)^{1/2}."),
    ("HV17", "H.finitely_many", "§2.7: finitely many resonances in any strip above -gamma_0",
     r"Since there are only finitely many resonances in any strip $\Im\sigma>-\alpha>-\gamma_0$"),
    ("HV17", "H.photon_unique", "§2.2: (r - r_P)(r^-2 mu)' < 0 for r != r_P in (r_2, r_3)",
     r"and $(r-r_P)(r^{-2}\mu_*)'<0$ for $r\neq r_P$"),
    ("CF17", "C.setting", "RNdS = D(Sigma), Sigma compact S^1 x S^2",
     r"we mean the maximal domain of dependence \mbox{$\cD (\Sigma)=\cM$} of a compact Cauchy hypersurface $\Sigma=\mathbb{S}^1\times\mathbb{S}^2$."),
    ("CF17", "C.v_chart", "the v chart covers M_int u H_A^+ u CH_B^+, with T = d_v",
     r"$$\cM_{int}\cup{\cal H}_A^+\cup{\cal CH}^+_B,$$"),
    ("CF17", "C.T_dv", "T = d_v on H_A^+",
     r"with $-\partial_r$ future directed and transverse  to  $\cH_A^+$ and $T=\partial_v$."),
    ("CF17", "C.u_chart", "the u chart covers M_int u CH_A^+ u H_B^+",
     r"$$\cM_{int}\cup{\cal CH}_A^+\cup{\cal H}_B^+,$$"),
    ("CF17", "C.u_orientation", "u runs against the time orientation; H_A^+ at u = +oo",
     r"Observe also that the orientation of $u$ is contrary to the spacetime orientation and as a consequence $\cH_A^+\supset\{u=+\infty\}$."),
    ("CF17", "C.price", "Assumption 2.2: Price law on H^+ = {r = r_+} in the v chart, 0 <= v_1 <= v_2",
     r"we have, in $\cH^+=\{r=r_+\}$, covered by the Eddington-Finkelstein coordinates~\eqref{EFv},"),
    ("CF17", "C.price_range", "Assumption 2.2: for all 0 <= v_1 <= v_2",
     r"for all $0\leq v_1\leq v_2$."),
    ("CF17", "C.thm31_hyp", "Thm 3.1: setting",
     r"Let $(\Mint, g)$ be a black hole interior region of a subextremal Reissner--Nordstr\"om--de Sitter spacetime. Let $\phi\in C^{\infty}(\Mint\cup \cH^+)$ be a solution of the wave equation~\eqref{waveEq} satisfying Price's Law (Assumption~\ref{pricesAss})."),
    ("CF17", "C.thm31_cond", "Thm 3.1: the condition",
     r"2\min\{p,\kappa_+\}>\kappa_-\;,"),
    ("CF17", "C.thm31_region", "Thm 3.1: the energy estimate is for u >= 0",
     r"for any $r_-\leq r_2\leq r_1\leq r_+$ and any $0\leq u_1\leq u_2$, we have"),
    ("CF17", "C.thm31_ext", "Thm 3.1: C^0 n H^1_loc extension",
     r"Moreover, $\phi$ extends as a function to $C^{0}(\Mint\cup\cH^+ \cup{\cal CH}^+)\cap H_{loc}^1(\Mint\cup\cH^+\cup{\cal CH}^+)$."),
    ("CF17", "C.surfgrav", "kappa_* = |D'(r_*)|/2",
     r"\kappa_*=\frac{1}{2}\left|D'(r_*)\right|\quad,\quad *\in\{-,+,c\}\;."),
    ("CF17", "C.roundgrad", "|grad_S2 phi|^2 = r^-2 sum (Omega_i phi)^2",
     r"|\nabb \phi|^2=\frac{1}{r^2}\sum_{i=1}^{3} \left(\leo_i \phi\right)^2\;,"),
]


def collapse(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def fetch(arxiv_id: str, name: str, dest: Path) -> Path:
    url = f"https://arxiv.org/e-print/{arxiv_id}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (P2 S5 review)"})
    data = urllib.request.urlopen(req, timeout=60).read()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
        member = tf.getmember(name)
        (dest / name).write_bytes(tf.extractfile(member).read())
    return dest / name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=None, help="directory holding cauchyhor-r.tex and RN_deSitterRevised6.tex")
    args = ap.parse_args()
    src = Path(args.src) if args.src else Path(tempfile.mkdtemp(prefix="s5rv_src_"))
    texts, meta = {}, {}
    for key, (aid, name, sha) in SOURCES.items():
        hits = list(src.rglob(name))
        p = hits[0] if hits else fetch(aid, name, src)
        raw = p.read_bytes()
        meta[key] = {"arxiv": aid, "file": name, "sha256": hashlib.sha256(raw).hexdigest(),
                     "sha256_matches_review": hashlib.sha256(raw).hexdigest() == sha}
        texts[key] = raw.decode("utf-8", errors="replace")
    rows = []
    for key, pid, what, passage in PASSAGES:
        text = texts[key]
        flat, pos_map = [], []
        # collapse whitespace while keeping a map to the original offsets (for line numbers)
        prev_space = False
        for i, ch in enumerate(text):
            if ch.isspace():
                if not prev_space:
                    flat.append(" "); pos_map.append(i)
                prev_space = True
            else:
                flat.append(ch); pos_map.append(i)
                prev_space = False
        flat_s = "".join(flat)
        j = flat_s.find(collapse(passage))
        line = text.count("\n", 0, pos_map[j]) + 1 if j >= 0 else None
        rows.append({"source": key, "id": pid, "what": what, "found": j >= 0, "line": line})
        print(f"{key} {pid:26s} {'found at l. ' + str(line) if j >= 0 else 'NOT FOUND'}")
    out = {"tool": "quotes_check.py (S5 review)", "sources": meta, "rows": rows,
           "all_found": all(r["found"] for r in rows)}
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print("all passages found:", out["all_found"], "| sources:", {k: v["sha256_matches_review"] for k, v in meta.items()})
    return 0 if out["all_found"] else 1


if __name__ == "__main__":
    sys.exit(main())
