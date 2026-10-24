"""Orchestrator for Track 2: fault-tree coverage, FMEA ranking, rank stability, readthrough check, figures.

    python -m mvapipe.track2.run [--kb knowledge/bub1b_mva1_kb.json] [--out DIR] [--no-network]

Reads only the curated knowledge file; optionally uses the Track 1 outputs (nmd_check.json,
locus_summary.json) for the variant-specific readthrough check. No coordinates are written to the outputs.
"""
from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import pandas as pd

from . import faulttree as ft
from . import fmea
from . import readthrough as rt
from . import scenarios as scn
from . import checkpoint_sim as sim


def load_kb(path=None) -> dict:
    path = Path(path or os.environ.get("MVA_KB") or Path(__file__).resolve().parents[2] / "knowledge" / "bub1b_mva1_kb.json")
    with open(path) as fh:
        return json.load(fh)


def analyse(kb: dict) -> dict:
    root, cons = kb["root_tree"], kb["consequence_tree"]
    for t in (root, cons):
        problems = ft.validate(t)
        if problems:
            raise ValueError("; ".join(problems))
    cuts = {"root": ft.cut_sets(root), "consequence": ft.cut_sets(cons)}
    rows = []
    for it in kb["interventions"]:
        cov = ft.coverage(cuts[it["tier"]], set(it["addresses"]))
        rows.append(fmea.score(it, cov))
    kept, stopped = fmea.rank(rows)
    for ex in kb["exclusions"]:
        stopped.append({"id": ex["id"], "label": ex["label"], "tier": ex["class"], "E": ex["E"], "S": ex["S"],
                        "stop": fmea.hard_stop({"E": ex["E"], "S": ex["S"], "approval": ex["approval"],
                                                        "opposite_mechanism": ex.get("opposite_mechanism")})
                        or "excluded on mechanistic grounds", "reason": ex["reason"]})
    stab = fmea.rank_stability(kept) if len(kept) > 1 else {}
    return {
        "minimal_cut_sets": {k: sorted(sorted(c) for c in v) for k, v in cuts.items()},
        "fully_mitigating_sets_root": sorted(sorted(s) for s in ft.fully_mitigating_sets(cuts["root"])),
        "ranking": kept, "excluded": stopped, "rank_stability": stab,
        "scenarios": {m: {k: v for k, v in scn.scenario_table(kb, m == "with_modifier").items() if k != "rows"}
                      for m in ("without_modifier", "with_modifier")} if "uncertainty" in kb else {},
    }


def readthrough_check(workdir: Path, network: bool = True) -> dict | None:
    nmd_f, loc_f = workdir / "nmd_check.json", workdir / "locus_summary.json"
    if not (nmd_f.exists() and loc_f.exists()):
        return None
    nmd = json.load(open(nmd_f))["results"]
    strand = json.load(open(loc_f)).get("strand", "+")
    for r in nmd:
        v = r.get("variant", {})
        m = re.match(r"p\.([A-Za-z]{3})\d+(\*|Ter)", str(v.get("hgvs_p", "")))
        if not m or len(str(v.get("ref"))) != 1 or len(str(v.get("alt"))) != 1:
            continue
        plus4 = minus1 = None
        note = ""
        cands = rt.candidate_stop_codons(m.group(1), v["ref"], v["alt"], strand)
        if network and len(cands) == 1:
            try:
                coord = rt.plus4_coordinate(int(v["pos"]), cands[0]["index"], strand)
                plus4 = rt.fetch_base(str(v["chrom"]), coord)
                minus1 = rt.fetch_base(str(v["chrom"]), rt.minus1_coordinate(int(v["pos"]), cands[0]["index"], strand))
            except Exception as e:                      # network or API failure: report, do not guess
                plus4 = minus1 = None
                note = f"+4 base not retrieved ({type(e).__name__})"
        else:
            note = "+4 base not retrieved (network disabled or ambiguous codon)"
        res = rt.assess(m.group(1), v["ref"], v["alt"], strand, bool(r.get("predicted_nmd")), plus4, minus1)
        if plus4 is None:
            res["note"] = note
        return res
    return None


def figure_tree(tree: dict, interventions: list[dict], out: Path, title: str, color="#1f4e8c"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    nodes = tree["nodes"]
    level = {tree["top"]: 0}
    order = [tree["top"]]
    for n in order:
        for c in nodes[n].get("children", []):
            level[c] = level[n] + 1
            order.append(c)
    leaves = [n for n in order if nodes[n]["type"] == "BASIC"]
    xs = {n: i for i, n in enumerate(leaves)}
    def xpos(n):
        if n in xs:
            return xs[n]
        v = [xpos(c) for c in nodes[n]["children"]]
        xs[n] = sum(v) / len(v)
        return xs[n]
    for n in order:
        xpos(n)
    depth = max(level.values())
    fig, ax = plt.subplots(figsize=(7.2, 1.9 + 1.25 * depth))
    for n in order:
        for c in nodes[n].get("children", []):
            ax.plot([xs[n], xs[c]], [-level[n], -level[c]], color="#8a93a6", lw=0.9, zorder=1)
    for n in order:
        kind = nodes[n]["type"]
        w = 0.82
        ax.add_patch(FancyBboxPatch((xs[n] - w / 2, -level[n] - 0.2), w, 0.4, boxstyle="round,pad=0.02",
                                    fc="white" if kind != "BASIC" else "#e6ecf7", ec=color, lw=1, zorder=2))
        lab = nodes[n]["label"]
        ax.text(xs[n], -level[n] + 0.0, "\n".join(_wrap(lab, 22)), ha="center", va="center", fontsize=5.8, zorder=3)
        if kind in ("AND", "OR"):
            ax.text(xs[n], -level[n] - 0.33, kind, ha="center", va="center", fontsize=6.5, color="#b3261e",
                    fontweight="bold", zorder=3)
    for it in interventions:
        for e in it["addresses"]:
            if e in xs and e in leaves:
                ax.text(xs[e], -level[e] - 0.62, it["id"], ha="center", va="center", fontsize=7, color="white",
                        bbox=dict(boxstyle="round,pad=0.25", fc="#b3261e", ec="none"), zorder=3)
    ax.set_xlim(-0.7, len(leaves) - 0.3); ax.set_ylim(-depth - 0.9, 0.5)
    ax.set_title(title, loc="left", fontsize=8.5)
    ax.axis("off"); fig.tight_layout(); fig.savefig(out); plt.close(fig)


def _wrap(s: str, n: int) -> list[str]:
    words, line, out = s.split(), "", []
    for w in words:
        if len(line) + len(w) + 1 > n and line:
            out.append(line); line = w
        else:
            line = (line + " " + w).strip()
    return out + [line]


def figure_benefit_risk(ranking: list[dict], excluded: list[dict], out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.axhspan(fmea.S_STOP - 0.5, 10.5, color="#f6d6d3", alpha=0.6, lw=0)
    ax.axvspan(0.5, fmea.E_MIN - 0.5, color="#e6e6e6", alpha=0.6, lw=0)
    for r in excluded:
        ax.plot(r["E"], r["S"], "x", color="#6b6b6b", ms=6)
        ax.text(r["E"] + 0.06, r["S"] + 0.12, r["id"], fontsize=6.5, color="#6b6b6b")
    for r in ranking:
        ax.scatter(r["E"], r["S"], s=60 + 420 * r["coverage"], color="#1f4e8c", alpha=0.8, edgecolor="#1f2a44")
        ax.text(r["E"] + 0.1, r["S"] + 0.2, f"{r['id']}  (RPN {r['RPN']})", fontsize=7)
    ax.set_xlim(0.5, 5.5); ax.set_ylim(0.5, 10.5)
    ax.set_xlabel("Evidence grade E (1 = speculative, 5 = randomised trial)")
    ax.set_ylabel("Severity S if the hypothesis is wrong")
    ax.text(0.6, 10.1, "hard stop: severity", fontsize=7, color="#b3261e")
    ax.text(0.6, 1.0, "hard stop: evidence", fontsize=7, color="#555", rotation=90)
    ax.set_title("Candidates (circles; area = coverage of failure paths) and exclusions (crosses)", loc="left", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def figure_sim(simres: dict, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = list(simres["interventions"])[::-1]
    pb = [simres["interventions"][n]["p_better"] for n in names]
    pn = [simres["interventions"][n]["p_no_change"] for n in names]
    pw = [simres["interventions"][n]["p_worse"] for n in names]
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    y = range(len(names))
    ax.barh(y, pb, color="#2e7d4f", label="fewer errors")
    ax.barh(y, pn, left=pb, color="#cfd4dc", label="no change")
    ax.barh(y, pw, left=[a + b for a, b in zip(pb, pn)], color="#b3261e", label="more errors")
    ax.set_yticks(list(y)); ax.set_yticklabels(names, fontsize=7)
    ax.set_xlim(0, 1); ax.set_xlabel("share of 5,000 random parameter sets (toy model, not fitted)", fontsize=7)
    ax.tick_params(axis="x", labelsize=7); ax.legend(fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.2), frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def run(kb_path=None, outdir=None, network: bool = True) -> dict:
    kb = load_kb(kb_path)
    res = analyse(kb)
    out = Path(outdir or os.environ.get("MVA_WORKDIR") or "work_track2"); out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(res["ranking"]).to_csv(out / "track2_ranking.csv", index=False)
    pd.DataFrame(res["excluded"]).to_csv(out / "track2_excluded.csv", index=False)
    rtk = readthrough_check(out, network)
    if rtk:
        res["readthrough_check"] = rtk
    simres = sim.run()
    raw = simres.pop("_raw")
    simres["dose_response_median_misseg"] = sim.dose_response(raw["params"])
    res["checkpoint_simulation"] = simres
    figure_sim(simres, out / "fig_t2_checkpoint_sim.pdf")
    figure_tree(kb["root_tree"], kb["interventions"], out / "fig_t2_root_tree.pdf", "Root-cause fault tree (AND of two alleles)")
    figure_tree(kb["consequence_tree"], kb["interventions"], out / "fig_t2_consequence_tree.pdf", "Consequence tree", "#6b4d9b")
    figure_benefit_risk(res["ranking"], [e for e in res["excluded"] if "E" in e and "S" in e], out / "fig_t2_benefit_risk.pdf")
    with open(out / "track2_summary.json", "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kb"); ap.add_argument("--out"); ap.add_argument("--no-network", action="store_true")
    a = ap.parse_args(argv)
    res = run(a.kb, a.out, not a.no_network)
    print(json.dumps({k: res[k] for k in ("minimal_cut_sets", "fully_mitigating_sets_root")}, indent=1))
    print(pd.DataFrame(res["ranking"]).to_string(index=False))
    print("\nExcluded:"); print(pd.DataFrame(res["excluded"])[["id", "label", "stop"]].to_string(index=False))
    print("\nRank stability:", json.dumps(res["rank_stability"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
