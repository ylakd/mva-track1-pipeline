"""Genotype-uncertainty analysis.

The fault tree assumes the two BUB1B alleles are in trans and that the missense allele is hypomorphic. Track 1 could
establish neither (unphased short reads; the missense is a variant of uncertain significance). Let

    r = P(trans) * P(missense hypomorphic)

be the probability that the compound-heterozygous genotype is the disease-causing one. For root-tier interventions
(which act on an allele) the effective coverage is c * r: if the genotype is not causal they have nothing to repair.
Consequence-tier interventions act on harm that is present whatever the cause, so their coverage is unchanged
(they remain gated by a biomarker). An optional modifier (a second tumour-predisposition variant) raises the severity
of interventions tagged with a matching risk.

Outputs: the order of candidates over a grid of r, the threshold at which the order changes, and a minimax view
(worst RPN over the grid).
"""
from __future__ import annotations

from . import faulttree as ft
from . import fmea


def score_under(item: dict, cuts: dict, r: float, modifier: dict | None = None) -> dict:
    """Score one candidate when the compound-heterozygous genotype is causal with probability r.

    Failure-probability proxy f = (10 - 1.2 E - 2 c) / 10 as in the base model. A root-tier hypothesis can only work if
    the genotype is causal, so f_eff = 1 - (1 - f) * r; consequence-tier hypotheses keep f. O = clamp(round(10 f_eff)).
    """
    cov = ft.coverage(cuts[item["tier"]], set(item["addresses"]))
    f = (10 - 1.2 * item["E"] - 2.0 * cov) / 10
    f_eff = 1 - (1 - f) * r if item["tier"] == "root" else f
    it = dict(item)
    if modifier:
        d = max([modifier["severity_delta"].get(t, 0) for t in item.get("tags", [])] + [0])
        it["S"] = min(10, item["S"] + d)
    o = fmea.clamp(10 * f_eff)
    row = {"id": it["id"], "label": it["label"], "tier": it["tier"], "coverage": round(cov, 2), "E": it["E"],
           "S": it["S"], "O": o, "D": it["D"], "RPN": it["S"] * o * it["D"], "stop": fmea.hard_stop(it), "r": r}
    return row


def scenario_table(kb: dict, use_modifier: bool = True) -> dict:
    u = kb.get("uncertainty", {})
    grid = u.get("r_grid", [0.2, 0.5, 0.8, 1.0])
    mod = u.get("modifier") if use_modifier else None
    cuts = {"root": ft.cut_sets(kb["root_tree"]), "consequence": ft.cut_sets(kb["consequence_tree"])}
    rows, orders = [], {}
    for r in grid:
        sc = [score_under(it, cuts, r, mod) for it in kb["interventions"]]
        kept, _ = fmea.rank(sc)
        orders[r] = [x["id"] for x in kept]
        rows += sc
    ids = [it["id"] for it in kb["interventions"]]
    worst = {i: max(x["RPN"] for x in rows if x["id"] == i) for i in ids if any(x["id"] == i and x["stop"] is None for x in rows)}
    best = {i: min(x["RPN"] for x in rows if x["id"] == i) for i in worst}
    changes = [(a, b) for a, b in zip(grid, grid[1:]) if orders[a] != orders[b]]
    return {"r_grid": grid, "orders": {str(k): v for k, v in orders.items()}, "order_changes_between_r": changes,
            "rpn_worst_over_grid": worst, "rpn_best_over_grid": best,
            "minimax_order": sorted(worst, key=worst.get), "modifier_applied": bool(mod), "rows": rows}
